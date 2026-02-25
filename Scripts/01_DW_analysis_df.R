## Bergmann's rule size metric comparison in Caprimulgids ##

# Data wrangling 01 -- Create nested PCA tbl and analysis_df #

# Load data & libraries --------------------------------------------------
library(tidyverse)
library(readr)

capri_df <- read_csv("Derived/Capri_BA_compare02.24.26.csv")

# Format ------------------------------------------------------------------
species_vec <- c("Whip-poor-will", "Nightjar", "Nighthawk")
capri_df2 <- capri_df %>%
  mutate(
    Species = factor(Species, levels = species_vec),
    Age = factor(Age, levels = c("Young", "Adult", "Unk")),
    Sex = factor(Sex)
  )

# Custom function to generate nested tbl object with derived surface-area-to-volume ratio (SAtoV) and the PCA 
prepare_species_data <- function(dat, species_name) {
  d <- dat %>%
    filter(Species == species_name) %>%
    filter(Sex %in% c("M", "F")) %>%
    filter(!is.na(B.Lat)) %>%
    drop_na(Wing.comb, Tail.comb, Mass.comb)
  
  if (species_name == "Nighthawk") {
    # Nighthawks were not reliably aged; keep all ages and model by sex.
    d <- d %>% mutate(Age_grp = "Unk")
  } else {
    d <- d %>%
      filter(Age %in% c("Young", "Adult")) %>%
      mutate(Age_grp = as.character(Age))
  }
  
  # Generate derived surface-area-to-volume ratio 
  d <- d %>%
    mutate(
      SAtoV.comb = (Wing.comb^2) / Mass.comb,
      Sex = droplevels(Sex)
    )
  
  if (species_name == "Nighthawk") {
    d <- d %>% mutate(stratum = as.character(Sex))
  } else {
    d <- d %>% mutate(stratum = paste(Age_grp, Sex, sep = "_"))
  }
  
  d
}

# Generate nested tbl
capri_df3 <- tibble(Species = species_vec) %>%
  mutate(data = map(Species, ~ prepare_species_data(capri_df2, .x)))

# Add pca object
capri_pca <- capri_df3 %>% mutate(pca = map(data, \(df){
  df_pca <- df %>% select(Wing.comb, Tail.comb, Mass.comb) 
  prcomp(df_pca, center = TRUE,  scale. = TRUE)
}))

# Extract pc1
capri_pca1 <- capri_pca %>% mutate(
  data = map2(
    data, pca, 
    ~ .x %>% mutate(PC1 = -.y$x[, 1]))
)

# Unnest tbl for analysis
analysis_df <- capri_pca1 %>%
  mutate(data = map(data, ~ select(.x, -Species))) %>%
  unnest(data) %>% 
  select(-pca)

# Export ------------------------------------------------------------------
# Caprimulgid nested df with PCA objects
saveRDS(capri_pca1, "Rdata/capri_pca1.rds")
# Dataframe that we will use for analysis
write_csv(analysis_df, "Derived/analysis_df.csv")
