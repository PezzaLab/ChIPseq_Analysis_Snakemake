# load libraries
library(tidyverse)
library(cowplot)
library(knitr)
library(rmarkdown)
# library(magick) # This is for image processing, not using for now

wd <- getwd()

(outfile <- 
  str_c(wd, "/", snakemake@output[["html_report"]] ))

# save.image(file = "html_report.RData")
rmarkdown::render("Resources/html_report.Rmd",
                  knit_root_dir = wd,
                  output_file = outfile)

# Check how can I make desitions on which rmakrdown document to use regarding the files present (HS, TSS, single strand or only both strands...)