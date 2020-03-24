# example: only saliency
library(tidyverse)
library(imager)
library(OpenImageR)
library(saliency)
library(rstan)
rstan_options(auto_write = TRUE)
library(here)
source(here::here("R", "helpers.R"))

im <- imager::load.image("https://raw.githubusercontent.com/NUS-VIP/predicting-human-gaze-beyond-pixels/master/data/stimuli/1001.jpg")
saliency::view(im)
sal <- saliency::itti_koch(im, c_scale = c(2, 3, 4), d_scale = c(3, 4))
saliency::view(sal)
sal <- saliency::downsample(sal, get_gaussian_kernel())
saliency::view(sal)
saldf <- matrix2df(sal, "top-left")
saldf$s_norm <- saldf$s / sum(saldf$s)

ggplot(saldf, aes(x = x, y = y, fill = s_norm)) + 
  geom_tile() + 
  scale_fill_gradient(low="black", high = "white")


simulate <- function(t_max, sigma_a, alpha){
  t <- 0
  fix <- sample.int(nrow(saldf), 1, prob = saldf$s_norm)
  counter <- 0
  while(t[counter] < t_max){
    
  }
}