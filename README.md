## Dynamic model of eye movements

This repository provides code and additional materials associated with the article 
> Kucharský, Š., van Renswoude, D., Raijmakers, M.E.J., & Visser, I. (2020). Dynamic model of eye movements (provisional title).

### Getting started

Running the model requires two crucial dependencies:

1. R (version 3.6.3 was used to create the original output). To download R, visit [https://www.r-project.org/](https://www.r-project.org/). Optionally it is nice to use R Studio ([https://rstudio.com/](https://rstudio.com/)).

2. Stan and R package rstan. Visit [https://mc-stan.org/](https://mc-stan.org/) for more information, installation instruction, manual, and tutorials.


To reproduce all code and output, there are additional dependencies, mostly in the form of additional R packages. The list of R packages (and their versions) is available in the [sessionInfo](sessionInfo) file.

Additionally, we calculated saliency maps of images created by Xu, et al. (2014) stored at [https://github.com/NUS-VIP/predicting-human-gaze-beyond-pixels](https://github.com/NUS-VIP/predicting-human-gaze-beyond-pixels) using the Python (version 3.7.7, [https://www.python.org/](https://www.python.org/)) scripts originally created by Mayo Yamasaki [https://github.com/mayoyamasaki/saliency-map](https://github.com/mayoyamasaki/saliency-map) and rewritten by Kohki Yamagiwa for Python 3 [https://github.com/tamanobi/saliency-map](https://github.com/tamanobi/saliency-map). However using these resources is not strictly necessary to produce the main output as the saliency maps were saved in folder [data/saliency](data/saliency/).

### Structure of the repository

1. [`stan/`](stan/) folder contains all necessary `.stan` files used for fitting the model, generating posterior predictives, and cross-validation.
	1. [`helpers/`](stan/helpers) subfolder contains definitions of functions that are used in the stan files. These functions are copied into the `.stan` scripts using the directive `#include "path_to_file.stan"` to be used as user defined functions.
	2. 

2. [`data/`](data/) folder contains data that are used in this project
	1. [`saliency/`](data/saliency/) folder contains all 700 stimuli from the repository by Xu, et al. (2014) ([https://github.com/NUS-VIP/predicting-human-gaze-beyond-pixels](https://github.com/NUS-VIP/predicting-human-gaze-beyond-pixels)) converted to saliency maps and saved as `.jpg` files. The folder contains [`get_saliency.py`](data/saliency/get_saliency.py) script that should (provided necessary dependencies are included, see above) reproduce the output.
	2. 
