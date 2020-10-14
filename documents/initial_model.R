# Script with the initial model simulations: later it will be rewritten in rmarkdown for a nice github document
library(tidyverse)
library(imager)
library(plotrix)
library(gtools)
library(rstan)
library(here)
source(here::here("R", "load_image.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "expose_helpers_stan.R"))

image_nr <- c(1001, 1014, 1049)
image_name  <- paste0(image_nr, ".jpg")
objects <- read.csv(here::here("data", "object_familiarity", "objects.csv"))
saliency <- read.csv(here::here("data", "saliency.csv")) %>% subset(image %in% image_nr)


par(mfcol = c(2, 3))
for(i in seq_along(image_name)) {
  img <- load_image(image_name[i])
  dim_img <- list(min_x = 0, min_y = 0, max_x = imager::width(img), max_y = imager::height(img))
  obj <- subset(objects, image_name == image_nr[i])
  sal <- subset(saliency, image == image_nr[i])
  
  plot(img, axes = FALSE)
  points(obj$x, obj$y, pch = 10, cex = 2)
   for(i in 1:nrow(obj))
     plotrix::draw.ellipse(x = obj$x[i], y = obj$y[i], a = obj$width[i]/2, b = obj$height[i]/2, angle = 0, lwd = 2)
  
  plot(imager::as.cimg(sal$value, x = max(sal$row), y = max(sal$col)), axes = FALSE)
}
par(mfrow = c(1, 1))

