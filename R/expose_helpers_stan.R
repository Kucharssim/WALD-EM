# loads stan functions as R functions
library(rstan)

functions_paths <- readLines(here::here("stan", "helpers", "load_functions.stan"))
functions_paths <- gsub("#include stan/helpers/", "", functions_paths)
functions_code  <- lapply(functions_paths, function(p) paste(readLines(here::here("stan", "helpers", p)), collapse = "\n"))
all_functions_code <- paste(functions_code, collapse = "\n")
all_functions_code <- paste0("functions{\n", all_functions_code, "\n}")

all_functions_code_c <- rstan::stanc(model_code = all_functions_code)
rstan::expose_stan_functions(all_functions_code_c)

rm(functions_paths, functions_code, all_functions_code, all_functions_code_c)
