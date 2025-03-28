pkgs <- c(
    "BiocManager", "biomaRt", "tidyverse", "OmnipathR")

for(x in pkgs){
    usethis::use_package(x, type = "depends")}