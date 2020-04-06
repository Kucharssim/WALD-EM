#!/usr/bin/python
# -*- coding: utf-8 -*-

# saliency_map is taken from https://github.com/tamanobi/saliency-map
# originally created by Mayo Yamasaki: https://github.com/mayoyamasaki/saliency-map
import sys
import numpy as np
import urllib.request as ur
from saliency_map import SaliencyMap
import cv2

# this link points to the repository associated with the article Xu, et al. (2014). Predicting human gaze beyond pixels. Journal of Vision. 
url_base = "https://raw.githubusercontent.com/NUS-VIP/predicting-human-gaze-beyond-pixels/master/data/stimuli/"
n_images = 700

def url_to_image(url):
  resp = ur.urlopen(url)
  image = np.asarray(bytearray(resp.read()), dtype="uint8")
  image = cv2.imdecode(image, cv2.IMREAD_COLOR)
  return image
  
image_index = 0  
while  image_index < n_images:
  image_index = image_index + 1
  image_name = str(1000 + image_index) + ".jpg"
  print("processing image: " + image_name)
  full_url = url_base + image_name
  img = url_to_image(full_url)
  sm  = SaliencyMap(img)
  sm_normalized = 255 * sm.map / sm.map.max()
  cv2.imwrite(image_name, sm_normalized)


