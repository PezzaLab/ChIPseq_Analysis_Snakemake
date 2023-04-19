# For debugging:
  # save.image(file = paste0("process_outfile_matrix_image.RData", runif(n=1,min=0, max = 9999)))

# Load libraries -----------
library(tidyverse)
# --------------------------------------------------------------------
# Functions ----------------------------------------------------------
# --------------------------------------------------------------------
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

# --------------------------------------------------------------------
# Process reads ------------------------------------------------------
# --------------------------------------------------------------------
files_grouped_by_strand <- 
  tail(snakemake@input, 3) # snakemake object has all objects individually, and then 3 lists with each of the groups (strands in this case), which makes everything very difficult...

file_names <- 
  lapply(
    files_grouped_by_strand,
    function(strand_group){
      lapply(strand_group, 
             function(file) {
               basename(file) %>% str_remove("\\..*")
             }
      )}
  )

# Get averages per coordinate
averages <- 
  lapply(
    tail(snakemake@input, 3),
    function(input){
      lapply(input, 
             average_signal_per_coordinate
             )}
    )

# Smooth
if (tolower(snakemake@params[["smooth"]]) == "true") {
  averages <-
    lapply(averages,
           function(strand_group){
             lapply(strand_group,
                    smooth_fun)
           })
  }

# Add prot, strand, hotspot region and library info
for (strand in seq_along(averages)) {
  for (file in seq_along(averages[[strand]])) {
    averages[[strand]][[file]][["Protein"]] <- 
      as.factor(file_names[[strand]][[file]])
    
    averages[[strand]][[file]][["Strand"]] <- 
      as.factor(names(averages)[[strand]] %>% 
      str_extract("83-163|99-147|both_strands"))
    
    averages[[strand]][[file]][["Region"]] <- 
      as.factor(snakemake@wildcards[["hs_region"]])
    
    averages[[strand]][[file]][["Library"]] <- 
      as.factor(snakemake@config$library$name)
  }
}

# Join all tables
table_all <- 
  map_df(averages, 
         ~reduce(.,bind_rows)
  )

# Normalize
final_table_dsDNA <- 
  normalize(table_all %>% filter(Strand == "both_strands"))

final_table_ssDNA <- 
  normalize(table_all %>% filter(Strand != "both_strands"))

final_table <- 
  bind_rows(final_table_dsDNA,
            final_table_ssDNA)

# Save
save(final_table,
     file = snakemake@output[[1]])