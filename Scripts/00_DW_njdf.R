## Bergmann's rule size metric comparison in Caprimulgids ##

# General  -----------------------------------------------------------------
# SessionInfo:
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS 26.2

## Instructions / suggestions for this repository
# All necessary data and code is provided to reproduce all analyses, figures, and tables in the main text. The file structure is provided, and if you open the .Rproj file everything will run as is. Otherwise, all that is necessary is that you set your working directory to the repository folder, like so: 

directory <- getwd()
setwd(directory) #Personalize, depending on where you've stored the repository folder

# Data wrangling 00 -- Create base night jar data frame (njdf) #
## This is a modified version of the data wrangling script 00 used in Skinner et al (2025). That script can be seen on github: 
# https://github.com/LezzGitIt/Bergmann-Round2/blob/f2bf3109aedc63c3dac92bf7d8eeb94f7ef76f63/Scripts/Final_JBI/00_DW_njdf_JBI.R 

# This script takes the combined data and formats, generates the time since sunset (tsss) variable, averages masses & wings from multiple captures by band age, & tidies and filters data to produce the capriBA data frame that will be used in downstream analyses

#Contents
# 1) Load & format data 
# 2) Condense age classes to Young, adult, or unknown
# 3) Calculate the difference between time since sunset & when an individual was captured
# 4) Create data frame for export 

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(stringi)
library(naniar)
library(readxl)
library(chron)
library(lutz)
library(suncalc)
library(zoo) 
library(geosphere)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

# Load and format data ----------------------------------------------------
capri.df <- read_excel("Data/Capri_df_combined_07.25.25.xlsx", sheet = "Data")

#Adjust time & date
capri.df <- capri.df %>% 
  replace_with_na_all(condition = ~.x %in% c(-99,-990, 9999, "<NA>", "-", ".", "na", 'NONABAND')) %>% 
  as.data.frame() %>% 
  filter(!is.na(B.Lat) & !is.na(Band.Number))
capri.df <- capri.df %>% 
  mutate(Species = ifelse(capri.df$Species == "Ceur" | capri.df$Species == "European Nightjar" | capri.df$Species == "European Nigthtjar", "EUNI", capri.df$Species), 
         Band.Number = stri_replace_all_regex(capri.df$Band.Number,
                                              pattern=c('-', ' '),
                                              replacement=c(''),
                                              vectorize=FALSE)) %>%
  mutate(uniqID = paste0(Band.Number, "_", Banding.Date, "_", !is.na(W.Lat))) %>% 
  mutate(Species = case_when(
    Species == "EUNI" ~ "Nightjar",
    Species == "EWPW" ~ "Whip-poor-will",
    Species == "CONI" ~ "Nighthawk",
  )) %>%
  arrange(is.na(W.Lat), Species) 

capri.df$Banding.Time <- str_pad(capri.df$Banding.Time, 4, pad = "0")
capri.df$Banding.Time <- sapply(str_split(parse_date_time(capri.df[,c("Banding.Time")], c("HMS"), truncated = 3), " "), function(x){x[2]})
capri.df$Banding.Time <- chron(times = capri.df$Banding.Time)
capri.df$Year <- str_pad(capri.df$Year, 3, pad = "0")
capri.df$Year <- str_pad(capri.df$Year, 4, pad = "2")

#Format times & dates, data types#
capri.df[,c("Wing.Chord","Mass", "Tail.Length", "B.Lat", "B.Long", "W.Lat", "W.Long", "Mig.dist", "Year")] <- lapply(capri.df[,c("Wing.Chord","Mass", "Tail.Length", "B.Lat", "B.Long", "W.Lat", "W.Long", "Mig.dist", "Year")], as.numeric)
capri.df[c("Banding.Date", "B.dep", "W.arr")] <- lapply(capri.df[c("Banding.Date", "B.dep", "W.arr")], parse_date_time, c("mdy", "ymd"))

## Create month day for understanding timing
capri.df[,c("Band.md","Bdep.md", "Warr.md")] <- lapply(capri.df[,c("Banding.Date", "B.dep", "W.arr")], format, "%m/%d")
capri.df[,c("Band.md","Bdep.md", "Warr.md")] <- lapply(capri.df[,c("Band.md","Bdep.md", "Warr.md")], as.Date, "%m/%d")

#Remove 13 individuals banded during the day (improbable capture time). 
capri.df2 <- capri.df %>% 
  filter(hms(Banding.Time) < hms("08:05:00") | 
           hms(Banding.Time) > hms("18:00:00") | 
           is.na(Banding.Time))

# Age classes condense ----------------------------------------------------
#Caprimulgids reach adult size (at least for feather length) after their second calendar year of life, so we need to condense age classes into Young, Adult, or Unknown

#Don't run this more than once
table(capri.df2$Age) #, capri.df$Species)
# #L & 1 are nestling / pullus, 3 = HY, 4 = AHY, 5 = SY, 6 = ASY, 8 = ASY
capri.df3 <- capri.df2 %>% filter(Age != "1" & Age != "L" & Age != "3") %>% 
  mutate(Age = case_when(Age == "4" ~ "Unk", #This should be adult? 
                         Age == "AHY" ~ "Unk",
                         Age == "ASY?" ~ "Unk",
                         Age == "5" ~ "Young",
                         Age == "SY" ~ "Young",
                         Age == "6" ~ "Adult",
                         Age == "4Y" ~ "Adult",
                         Age == "8" ~ "Adult",
                         Age == "A4Y" ~ "Adult",
                         Age == "A5Y" ~ "Adult",
                         Age == "ASY" ~ "Adult",
                         Age == "ATY" ~ "Adult",
                         Age == "TY" ~ "Adult"),
         Sex = trimws(Sex))

# Outliers ----------------------------------------------------------------
# Visualize
Species <- c("Nightjar", "Whip-poor-will", "Nighthawk")
map(Species, \(sp){
  capri.df3 %>% filter(Species == sp) %>%
    ggplot(aes(x = Tail.Length)) +
    geom_histogram() +
    labs(title = sp)
})

# Remove 5 outliers with unrealistically large tails - I discussed with Elly Knight (provider of CONI data) & she thinks these are indeed errors. 
aj_tbl <- capri.df3 %>% filter(Species == "Nighthawk" & Tail.Length > 180)
capri.df4 <- capri.df3 %>% anti_join(aj_tbl)

# Band.Age -------------------------------------------
# Given that some individuals were captured multiple times we want to incorporate that information to have a better estimate of body size. Sizes differ by age so we want to combine morphological variables of individuals that are the same age. Sex doesn't change overtime so don't need to worry about that. 

# Create Band.Age variable, where individuals of the same age will have the same 'Band.Age', and use mutate to average morphologies (and time since sunset cov) by Band.Age
capri.df5 <- capri.df4 %>% 
  mutate(Band.Age = paste0(Age, "_", Band.Number)) %>%
  group_by(Band.Age) %>% 
  mutate(Wing.comb = mean(Wing.Chord, na.rm = TRUE), 
         Mass.comb = mean(Mass, na.rm = TRUE),
         Tail.comb = mean(Tail.Length, na.rm = TRUE)) 

# CapriBA --------------------------------------------------------
##CapriBA has only a single row for each Band & Age combo
capriBA_compare <- capri.df5 %>% 
  group_by(Band.Age) %>% #BA = band age
  arrange(is.na(W.Lat), Year, .by_group = TRUE) %>% 
  slice_head()

# Export capriBA ---------------------------------------------------------
# Write capriBA.. Note uniqID is truly unique, but this df will still be further filtered
write.csv(capriBA_compare, file = paste0("Derived/Capri_BA_compare", format(Sys.Date(), "%m.%d.%y"), ".csv"), row.names = FALSE)
