# For debugging 
# save.image(file = paste0("wragle_X_nonPAR_asymmetric_HS_image.",
#                          runif(n=1, min=0, max = 9999),
#                          ".RData"))
# Load library ----
library(tidyverse)

# Load aggregate profiles ----
file_names <-  
  snakemake@input %>% 
  unlist %>% 
  basename %>% 
  str_remove(".RData") %>% 
  str_remove(paste0(snakemake@wildcards$libary, "_"))

results <- list()

for (i in seq_along(snakemake@input)){
  load(snakemake@input[[i]])
  results[[file_names[i]]] <- final_table
}

top_5000_plus_minus_2000 <- 
  results$top_5000_plus_minus_2000 %>% 
  mutate("Region" = "Top_5000_HSs")

# Join autosomal and XnonPAR ----
xNonPar_autosomal <- 
  bind_rows(results$autosomal_x_non_par_ctrl %>% mutate("Region" = "Autosomal"),
            results$x_non_par %>% mutate("Region" = "X-nonPAR"))

# Asymetric  ----
# Both strands are averaged to get one single, smoother profile.  
# There are 766 Watson sites, and 768 Crick sites, hence a weighted mean is done for each side.

results$asymetric_watson_strong <-
  results$asymetric_watson_strong %>% 
  mutate("Region" = "asymetric",
         "strong_strand" = "watson")

results$asymetric_crick_strong <-
  results$asymetric_crick_strong %>% 
  mutate("Region" = "asymetric",
         "strong_strand" = "crick")

## Merge both watson and crick ----
asymetric <- 
  bind_rows(results$asymetric_watson_strong, results$asymetric_crick_strong)

asymetric <- 
  asymetric %>%
  mutate("Average_signal" = ifelse(strong_strand == "watson",
                                   Average_signal * 766/1534,
                                   Average_signal * 768/1534))
# Strong strand is always on the right side (no need to flip anything) and **83-163 & crick** or **99-147 & watson** are weak.

## Separate weak from strong ----
weak_side_a <-
  filter(asymetric,
         Strand == "83-163",
         strong_strand == "crick")
weak_side_b <-
  filter(asymetric,
         Strand == "99-147",
         strong_strand == "watson")
weak_side <- bind_rows(weak_side_a, weak_side_b) %>%
  mutate(Strength = "weak")

strong_side_a <-
  filter(asymetric,
         Strand == "99-147",
         strong_strand == "crick")

strong_side_b <-
  filter(asymetric,
         Strand == "83-163",
         strong_strand == "watson")

strong_side <- bind_rows(strong_side_a, strong_side_b) %>% 
  mutate(Strength = "strong")
## Merge all weaks/strongs ----
weak_side_merged <- 
  weak_side %>% 
  group_by(Protein, Coordinates, Strength, Library) %>% 
  summarise(
    Average_signal = mean(Average_signal)
  )

strong_side_merged <- 
  strong_side %>% 
  group_by(Protein, Coordinates, Strength, Library) %>% 
  summarise(
    Average_signal = mean(Average_signal)
  )

## Re-join weak and strong ----
asymetric <- 
  bind_rows(strong_side_merged, weak_side_merged)

## Flip strong side ----
# Prepare a set in which strong and weak are on the same side of DSB, for direct 
# comparison. Flip strong side so both are on leeft side of DSB
strong_flipped <- 
  filter(asymetric, Strength == "strong") %>% 
  mutate(Coordinates = Coordinates * -1)

asymetric_flipped <- 
  filter(asymetric, Strength == "weak") %>% 
  bind_rows(strong_flipped)

# Save ----
  # I use "lst()" from tidyverse because it keeps the original name of the 
 # object when creating the list. With list() you loose the objects name.
aggregate_profiles <- 
  lst(asymetric,
       asymetric_flipped,
       xNonPar_autosomal,
       top_5000_plus_minus_2000)

save(aggregate_profiles,
     file = snakemake@output[[1]])