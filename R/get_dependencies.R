library(here)
library(knitr)
library(tidyverse)

get_packages <- function(file){
  lines <- readLines(file, warn = FALSE)
  library_lines <- lines[grep("^library\\(", lines)]
  packages <- gsub("library\\(", "", library_lines)
  packages <- gsub("\\)", "", packages)
  
  return(packages)
}

# list R files
files <- list.files(path = here::here(), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
# list Rmd files
files <- c(files,
           list.files(path = here::here(), pattern = "\\.Rmd$", recursive = TRUE, full.names = TRUE)
           )

# get pacakges
packages <- lapply(files, get_packages)
packages <- unique(unlist(packages))

pkgdata <- tibble::tibble(package = packages, 
                          version = sapply(packages, function(p) paste(packageVersion(p), collapse = ".")))
readr::write_csv(pkgdata, here::here("packages.csv"))

for (p in packages) library(p, character.only = TRUE) 
sink(file = here::here("sessionInfo"))
sessionInfo()
sink()

knitr::write_bib(packages, file = here::here("packages.bib"))
