## Bergmann's rule size metric comparison in Caprimulgids ##

# Data wrangling 01 -- Create nested PCA tbl and analysis_df #

# Load data & libraries --------------------------------------------------
library(tidyverse)
library(readr)
library(dagitty)
library(lavaan)

capri_df <- read_csv("Derived/Capri_BA_compare03.29.26.csv")

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


# DAG-data consistency ----------------------------------------------------
# See 

# NOTE: Age is placeholder for age / sex 
dag_mod <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"
  Site      [pos="0, -0.45"]
  Latitude  [exposure, pos="-0.25, 0"]
  Age       [pos="0.25,0"]
  Size      [outcome, pos="0,0.45"]

  Site -> Latitude
  Site -> Age
  Latitude -> Size
  Age -> Size
}
')

# Examine the conditional independencies implied by the DAG
impliedConditionalIndependencies(dag_mod)

# Visualize causal assumptions
plot(dag_mod, node.names = c("Age_sex" = "Age / sex"))

# Format so names match the DAG, separate into list
df_dag_tests_l <- analysis_df %>% 
  rename(Site = Site.name,
         Size = Wing.Chord,
         Latitude = B.Lat) %>%
  group_split(Species)

# Select variables, turn factors into numeric 
df_dag_tests_l2 <- map(df_dag_tests_l, \(df){
  df %>% select(Site, Age, Latitude, Size) %>%
    mutate(Site = as.numeric(as.factor(Site)),
           Age = as.numeric(Age))
})

# Use lavaan to generate the polychoric correlation matrix
lav_corr <- map(df_dag_tests_l2, \(df){
  lavCor(df)
})

# Test DAG-data consistency
map2(df_dag_tests_l2, lav_corr, \(df, cor){
  localTests(x = dag_mod , sample.cov = cor, sample.nobs = nrow(df))
})

## Conclusions: This 'works', as in it produces estimates and p-values for the the conditional independencies implied by the DAG. However, as Ankan, Wortel, & Textor (2021, doi:10.1002/cpz1.45) recognize, you can't test with categorical variables that have more than two levels.

# > Reduced (remove site) --------------------------------------------
# What about removing Site and using latitude as a proxy for site? 

# Remove site
dag_mod_red <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"
  Latitude  [exposure, pos="-0.25, 0"]
  Age       [pos="0.25,0"]
  Size      [outcome, pos="0,0.45"]

  Latitude -> Age
  Latitude -> Size
  Age -> Size
}
')
plot(dag_mod_red)

# No conditional independencies 
impliedConditionalIndependencies(dag_mod_red)

# No independencies to test so returns empty dataframes
map2(df_dag_tests_l2, lav_corr, \(df, cor){
  localTests(x = dag_mod_red , sample.cov = cor, sample.nobs = nrow(df))
})

# Export ------------------------------------------------------------------
stop()
# Caprimulgid nested df with PCA objects
saveRDS(capri_pca1, "Rdata/capri_pca1.rds")
# Dataframe that we will use for analysis
write_csv(analysis_df, "Derived/analysis_df.csv")
