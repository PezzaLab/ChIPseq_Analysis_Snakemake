library(tidyverse)
save.image("image_modify_annotated_peaks_script.RData")
annotated_table <- 
  read_tsv(snakemake@input[[1]]) %>% 
  mutate("Genomic_feature" = str_extract(Annotation,
                                         "^[:alnum:]+(-|'){0,1} *[:alnum:]*")) %>%
  rename("peak_ID"= 1)  %>% 
  rename("-log(q-value)" = `Peak Score`) %>% 
  mutate("-log(q-value)" = `-log(q-value)` / 10,
         "IGV_coordiates" = str_c(Chr, ":", Start,"-", End))

write_tsv(annotated_table,
          file = snakemake@output[[1]])