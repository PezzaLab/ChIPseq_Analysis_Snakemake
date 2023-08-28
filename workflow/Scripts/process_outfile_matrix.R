# For debugging:
# file_name <- paste0(
#   "process_outfile_matrix_image.",
#   runif(n = 1, min = 0, max = 9999),
#   ".RData")
# save.image(file = file_name)
# 
# print(paste("File name:", file_name))
#  

# Load libraries -----------
library(tidyverse)
# Functions ----------------------------------------------------------
get_strands <- 
  function(path){
    file_name <- 
      path %>% 
      basename
    if (grepl("mit_filt\\.matrix$", file_name)) {
      strand <- as.factor("Both_strands")
    } else {
      strand <- 
        as.factor(str_extract(
          file_name,"(?<=mit_filt\\.)83-163|99-147|inc_16|exc_16(?=\\.matrix)"))
    }
    strand
  }

get_file_name <- function(path) basename(path) %>% str_remove("\\..*") %>% as.factor

get_role <- function(path) {
  role <- path %>%
    str_extract("(?<=B6xCAST_PRDM9_assymetric_hs_)invading|receiving") %>%
    as.factor
  if (is.na(role)) {
    role <- "invading_and_receiving"
  }
  return(role)
  }

average_signal_per_coordinate <- 
  function(matrix_file) {
    # Read file
    suppressMessages(
      matrix <- read_tsv(matrix_file, col_names = FALSE, skip = 3)
    )
    line1 <- read_lines(matrix_file, n_max = 1)
    line2 <- read_lines(matrix_file, n_max = 1, skip = 1)
    # Transform NA to 0
    matrix[is.na(matrix)] <- 0 #is.na recognizes both "NA" and "NaN", while is.nan recognizes only "NaN"
    # Get number of genes
    number_genes <- str_extract(line1, "(?<=genes:).*")
    # Get number of columns
    n_cols <- length(colnames(matrix))
    # Get bin size
    bin_size <- as.integer(str_extract(line2, "(?<=bin size:)\\d*"))
    # Get coordinates
    coordinates <- as.character(seq(-( n_cols*bin_size/2 - bin_size/2),
                                    (n_cols*bin_size/2 - bin_size/2), 
                                    by = bin_size))
    # Do average for each column
    average <- list()
    for (j in 1:n_cols) {
      average[[j]] <- mean(matrix[[j]])
    }
    rm(j)
    average2 <- unlist(average)
    # Generate table with average value
    table <- tibble("Coordinates" = as.double(coordinates),
                    "Average_signal" = average2)
  # Return
  return(table)
}

smooth_fun <-
  function(table) {
    knots = nrow(table) / 30
    a <- smooth.spline(table, nknots = knots)
    b <- tibble("Coordinates" = a$x,
                "Average_signal" = a$y)
    b
  }
    
normalize <- 
  function(final_table) {
    min_max <- 
      final_table %>% 
      group_by(Protein) %>% 
      summarise(max = max(Average_signal), 
                min = min(Average_signal))
    
    norm_table <- list()
    
    for (i in levels(final_table$Protein)) {
      df_i <- filter(final_table, Protein == i)
      Min <- as.double(min_max$min[min_max$Protein == i])
      Max <- as.double(min_max$max[min_max$Protein == i])
      Range <- Max-Min
      df_i_min_corrected_signal <- 
        df_i %>% mutate(Average_signal = Average_signal - Min)
      area <- sum(df_i_min_corrected_signal$Average_signal)
      norm_table[[i]] <- 
        df_i_min_corrected_signal %>% 
        mutate("Mean_Area_Normalized_Coverage" = Average_signal * 100 / area,
               "Mean_Max_Normalized_Coverage" = Average_signal / Range)
    }
    
    rm(i)
    
    final_table_norm <- reduce(norm_table, bind_rows)
    return(final_table_norm)
  }

# Process reads ------------------------------------------------------
strands <- lapply(snakemake@input, get_strands)

file_names <- 
  lapply(snakemake@input, get_file_name)

role <- get_role(snakemake@wildcards[["hs_region"]])

## Get averages per coordinate --------------------------------
averages <- 
  lapply(snakemake@input, average_signal_per_coordinate)

## Smooth --------------------------------
if (tolower(
  snakemake@config[['aggregate_profiles']][['process_outfile']][['smooth']]) == "true") {
  averages <-
    lapply(averages,smooth_fun)
  }

## Add prot, strand, hotspot region, library and role metadata -------------
for (i in seq_along(averages)) {
    averages[[i]][["Protein"]] <- 
      file_names[[i]]
    
    averages[[i]][["Strand"]] <- 
      strands[[i]]
    
    averages[[i]][["Region"]] <- 
      as.factor(snakemake@wildcards[["hs_region"]])
    
    averages[[i]][["Library"]] <- 
      as.factor(snakemake@config$library$name)
    
    averages[[i]][["Role"]] <- 
      role
    }

## Join all tables --------------------------------
table_all <- 
  reduce(averages,bind_rows)

## Normalize --------------------------------
final_table_dsDNA <- 
  normalize(table_all %>% filter(Strand == "Both_strands"))

final_table_ssDNA <- 
  normalize(table_all %>% filter(Strand != "Both_strands"))

final_table <- 
  bind_rows(final_table_dsDNA,
            final_table_ssDNA)

## Save --------------------------------
save(final_table,
     file = snakemake@output[[1]])


# Debugging plot ----
# ggplot(final_table,
#        aes(x = Coordinates,
#            y = Mean_Max_Normalized_Coverage,
#            color = Strand)) +
#   geom_line() +
#   facet_wrap(~Protein)
