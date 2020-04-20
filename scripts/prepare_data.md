Data Preparation
================
Simon Kucharsky
2020-04-20

## Eye movement data

The eye movement data comes from
(<span class="citeproc-not-found" data-reference-id="Renswoude">**???**</span>),
and is based on the stimulus materials provided in Xu et al. (2014).

First, we clean the data and prepare relevant indexes (i.e., id of
participants need to be recoded from factors to integers for Stan).

NB. All data are in a coordinate system with origin (0, 0) in the
top-left of the picture (i.e., compared to classic cartesian coordinate
system, the y-axis is inverted). This is a standard in eye-tracking and
image processing, but it is mentioned here to avoid confusion.

``` r
library(tidyverse)
library(imager)
library(here)
source(here::here("R", "load_image.R"))
ggplot2::theme_set(ggplot2::theme_classic())
# load data
df <- read.csv(here::here("data", "object_familiarity", "Fixations_all.csv"), sep = ";", dec = ",")

# get relevant variables
df_clean <- dplyr::select(df, PP, Trial, Order, image, mean_x, mean_y, Duration)
names(df_clean) <- c("ppt", "trial", "order", "image", "x", "y", "duration (ms)")

# convert duration from ms to sec
df_clean$duration <- df_clean[['duration (ms)']]/1000

# convert indexes to integers
df_clean$id_ppt <- as.integer(df_clean$ppt)
df_clean$id_img <- as.integer(as.factor(df_clean$image))

df_clean$image  <- as.character(df_clean$image)

# arrange by ppt, img, and order of fixations
df_clean <- dplyr::arrange(df_clean, id_ppt, id_img, order)

# reorder variables
df_clean <- dplyr::select(df_clean, id_ppt, id_img, order, ppt, image, trial, x, y, duration)

# save
readr::write_csv(df_clean, path = here::here("data", "object_familiarity", "fixations.csv"))
rm(df, df_clean)
# load again (check whether it can be loaded correctly)
cols_spec <- readr::cols(
  id_ppt = readr::col_integer(),
  id_img = readr::col_integer(),
  order  = readr::col_integer(),
  ppt    = readr::col_character(),
  image  = readr::col_character(),
  trial  = readr::col_integer(),
  x      = readr::col_double(),
  y      = readr::col_double()
)
df <- readr::read_csv(file = here::here("data", "object_familiarity", "fixations.csv"), col_types = cols_spec)
```

Here, we assign each trial into a sample which is used for model
fitting, or model validation. We want to have some data for each
participant, and some data for each image, in each of the samples.

``` r
n_fixations <- df %>% 
  dplyr::group_by(id_ppt, id_img) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(id_ppt, id_img, fill = list(n = 0)) %>%
  dplyr::mutate(is_na = n == 0, is_not_na = n != 0)

n_fixations %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(id_ppt), y = as.factor(id_img), alpha = n)) + 
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_point(data = n_fixations %>% subset(is_na), ggplot2::aes(x = id_ppt, y = id_img), shape = 4, size = 2, alpha = 1) + 
  ggplot2::coord_equal() + 
  ggplot2::xlab("Participant") + 
  ggplot2::ylab("Image") + 
  ggplot2::ggtitle("Number of fixations") 
```

<img src="prepare_data_files/figure-gfm/n_fixations-1.png" width="100%" style="display: block; margin: auto;" />

First, we see that some participants (5, 16) have less than half of the
trials, and so we will group the data by participants, and split each
subset on two roughly equal in the number of trials.

``` r
set.seed(2020)
df %<>% plyr::ddply(.variables = "id_ppt", .fun = function(d){
  img   <- unique(d$id_img)
  n_img <- length(img)
  n_train <- floor(n_img / 2)
  train <- sample(x = c(rep(TRUE, n_train), rep(FALSE, n_img - n_train)), size = n_img, replace = FALSE)
  
  d$train <- logical(length = nrow(d))
  
  for(i in 1:nrow(d)) d$train[i] <- train[which(img == d$id_img[i])]
  d
})
```

``` r
trains <- df %>%
  dplyr::group_by(id_ppt, id_img) %>%
  dplyr::summarise(train = all(train))

trains %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(id_ppt), y = as.factor(id_img), fill = as.factor(train))) + 
  ggplot2::geom_tile() +
  ggplot2::geom_point(data = n_fixations %>% subset(is_na), ggplot2::aes(x = id_ppt, y = id_img), shape = 4, size = 2, alpha = 1, inherit.aes = FALSE) + 
  ggplot2::coord_equal() + 
  ggplot2::xlab("Participant") + 
  ggplot2::ylab("Image") + 
  ggplot2::ggtitle("Used for estimating parameters") +
  ggplot2::scale_fill_manual(values = c("#E69F00", "#999999"),
                             name   = NULL,
                             breaks = c("TRUE", "FALSE"),
                             labels = c("Yes", "No")
                             )
```

<img src="prepare_data_files/figure-gfm/plot_train-1.png" width="100%" style="display: block; margin: auto;" />

``` r
par(mfrow=c(1, 2))
barplot(df %>% subset(train) %>% .$id_img %>% table(), main = "Fixations in training set", xlab = "Image", ylab = "Count")
barplot(trains %>% subset(train) %>% .$id_img %>% table(), main = "Trials in training set", xlab = "Image", ylab = "Count")
```

<img src="prepare_data_files/figure-gfm/n_per_img-1.png" width="100%" style="display: block; margin: auto;" />

``` r
par(mfrow=c(1, 2))
barplot(df %>% subset(train) %>% .$id_ppt %>% table(), main = "Fixations in training set", xlab = "Participant", ylab = "Count")
barplot(trains %>% subset(train) %>% .$id_ppt %>% table(), main = "Trials in training set", xlab = "Participant", ylab = "Count")
```

<img src="prepare_data_files/figure-gfm/n_per_ppt-1.png" width="100%" style="display: block; margin: auto;" />

## Saliency data

We already preprocessed the stimuli images using the Itti and Koch
algorithm. The processed images are in
[`/data/saliency/`](/data/saliency/). `Python 3.7` script is available
at [`get_saliency.py`](/data/saliency/get_saliency.py) and requires
cloning the Itti and Koch saliency repository from
<https://github.com/tamanobi/saliency-map>, originally created by Mayo
Yamasaki.

We downsample the image saliency by a factor of 40 in each dimension.
First, we apply the gaussian blur with sd = 20 and range of 40, then
take every 40<sup>th</sup> row and column. That means that a picture of
dimensions 800 by 600 pixels will have a downsampled saliency map of
size 20 by 15 aggregated pixels. Each of the aggregated pixels then
subtends an area of 300 original pixels.

``` r
library(imager)
library(OpenImageR)
source(here::here("R", "load_image.R"))

img <- unique(df$image)
img_names <- paste0(img, ".jpg")

# load images that are used in the experiment
saliency <- lapply(img_names, load_image, folder = here::here("data", "saliency"))
# downsample
saliency_downsampled <- lapply(saliency, function(s) imager::as.cimg(OpenImageR::down_sample_image(as.matrix(s), down_pars$factor, TRUE, down_pars$sigma, down_pars$range)))
names(saliency) <- names(saliency_downsampled) <- img

# save downsampled images
for(i in img){
  imager::save.image(im = saliency_downsampled[[i]],
                     file = here::here("data", "saliency", sprintf("%s_downsampled.jpg", i)),
                     quality = 1)
}
```

Here, we normalize the saliency map to that each map sums up to 1 and
reshape for convenience.

``` r
# convert to long format
saliency_normalized <- lapply(saliency_downsampled, as.data.frame)
image_key <- dplyr::select(df, id_img, image) %>% unique()
for(i in img){
  s <- saliency_normalized[[i]]
  # normalize map
  s$value_normalized <- s$value / sum(s$value)
  # keep info about index of pixel
  s$row <- s$x
  s$col <- s$y
  s$idx <- seq_len(nrow(s))
  # rescale to original coordinates
  resc <- function(x, factor = 2) {(factor*(x-1) + 1 + factor*x)/2}
  s$x <- resc(s$x, down_pars$factor)
  s$y <- resc(s$y, down_pars$factor)
  
  s$image <- i
  s$id_img <- image_key$id_img [image_key$image == i]
  # compute the log of the normalized saliency map (discrete)
  s$saliency_log <- log(s$value_normalized)
  
  # compute the log of the density saliency map (standardize by area of pixels)
  s$log_lik_saliency <- s$saliency_log - log(down_pars$area)
  
  saliency_normalized[[i]] <- select(s, id_img, image, row, col, idx, x, y, value, value_normalized, saliency_log, log_lik_saliency)
}

saliency_normalized <- dplyr::bind_rows(saliency_normalized)

readr::write_csv(saliency_normalized, path = here::here("data", "saliency.csv"))
```

## Combination of saliency and eye movement data

For models that contain saliency as one part of the explanatory
variables, we can copy the value of the log-likelihood of the saliency
into the fixation data sets for the model of where. For the model of
when, we consider only aggregated pixels in a certain distance from the
current location. The distance (82.3422716 pixels) is calculated by
converting width of foveal attention (5 degrees of visual angle), given
the distance from the screen (60 cm), width of the monitor (51 cm), and
number of pixels in the horizontal direction (800).

``` r
df$log_lik_saliency <- numeric(length = nrow(df))
saliency_log_list <- list()

pb <- dplyr::progress_estimated(nrow(df))
for(i in seq_len(nrow(df))){
  pb$tick()$print()
  x <- df$x[i]
  y <- df$y[i]
  current_img <- df$id_img[[i]]
  s <- subset(saliency_normalized, id_img == current_img) # get saliency of the correct image
  
  distances <- sqrt((x-s$x)^2 + (y-s$y)^2)
  mean_sq_distances <- distances^2 / 2
  
  which_closest <- which.min(distances) # get index of the aggregated pixel the fixation location is inside of
  df$log_lik_saliency[i] <- s$log_lik_saliency[which_closest] # copy the corresponding log normalized saliency into the dataset
  df$n_neighbors[i] <- sum(distances < radius)
  saliency_log_list[[i]] <- data.frame(
    distances          = distances[distances < radius],
    mean_sq_distances  = mean_sq_distances[distances < radius],
    saliency_log       = s$saliency_log[distances < radius]
    )
}
```

``` r
#prepare elements of saliency_log as structures passable to Stan
max_neighbors     <- max(df$n_neighbors)
mean_sq_distances <- matrix(999, ncol = max_neighbors, nrow = nrow(df))
saliency_log      <- matrix(999, ncol = max_neighbors, nrow = nrow(df))
for(i in seq_len(nrow(df))){
  mean_sq_distances[i, 1:df$n_neighbors[i]] <- saliency_log_list[[i]]$mean_sq_distances
  saliency_log     [i, 1:df$n_neighbors[i]] <- saliency_log_list[[i]]$saliency_log
}
```

## Prepare objects data

We also extracted the information about the objects on the scenes.

``` r
obj <- read.csv(here::here("data", "object_familiarity", "cen700.txt"), sep = "", dec = ".")
obj$image <- as.character(1000 + 1)

objects <- dplyr::tibble(id_img = image_key$id_img[sapply(obj$image, function(i) which(i == image_key$image))],
                         image  = obj$image,
                         x      = obj$x, 
                         y      = obj$y,
                         width  = obj$max_x - obj$min_x,
                         height = obj$max_y - obj$min_y,
                         min_x  = obj$min_x,
                         min_y  = obj$min_y,
                         max_x  = obj$max_x,
                         max_y  = obj$max_y)
```

``` r
readr::write_csv(df,                    path = here::here("data", "fixations.csv"))
readr::write_csv(df %>% subset(train),  path = here::here("data", "fixations_train.csv"))
readr::write_csv(df %>% subset(!train), path = here::here("data", "fixations_validate.csv"))
readr::write_csv(as.data.frame(mean_sq_distances), path = here::here("data", "mean_sq_distances.csv"))
readr::write_csv(as.data.frame(saliency_log),      path = here::here("data", "saliency_log.csv"))
readr::write_csv(objects,               path = here::here("data", "objects.csv"))

save(df, mean_sq_distances, saliency_log, saliency_normalized, image_key, objects, 
     file = here::here("data", "cleaned_data.Rdata"))
```

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Xu2014beyond">

Xu, J., Jiang, M., Wang, S., Kankanhalli, M. S., & Zhao, Q. (2014).
Predicting human gaze beyond pixels. *Journal of Vision*, *14*(1),
28–28. <https://doi.org/10.1167/14.1.28>

</div>

</div>
