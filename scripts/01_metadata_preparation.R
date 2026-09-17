###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 01_metadata_preparation.R
## Purpose : Create clean metadata table for all analyses
###############################################################

rm(list = ls())

source("scripts/00_setup.R")

library(tidyverse)
library(readr)

############################################################
## Locate metadata file
############################################################

meta_file <- list.files(
    METADATA_DIR,
    pattern = "\\.(csv|txt|tsv)$",
    full.names = TRUE
)

if (length(meta_file) == 0) {
    stop("No metadata file found in data/metadata/")
}

meta_file <- meta_file[1]

cat("Reading metadata:\n")
cat(meta_file, "\n")

############################################################
## Import metadata
############################################################

if (grepl("\\.csv$", meta_file)) {

    meta_raw <- read_csv(meta_file)

} else {

    meta_raw <- read_tsv(meta_file)
}

cat("Metadata dimensions:\n")
print(dim(meta_raw))

############################################################
## Inspect column names
############################################################

cat("\nColumn names:\n")
print(colnames(meta_raw))

############################################################
## Rename GEO columns
############################################################
## Adjust these names if GEO uses different labels
############################################################

if (!"geo_accession" %in% names(meta_raw)) {
    stop("Required metadata column not found: geo_accession")
}

meta <- meta_raw
names(meta)[names(meta) == "geo_accession"] <- "SampleID"

############################################################
## Extract animal ID
############################################################

animal_col <- intersect(
    c("animal:ch1", "animal", "animal_id", "subject_id"),
    names(meta)
)[1]

if (is.na(animal_col)) {
    stop("Required metadata column not found for animal ID")
}

meta$Animal <- meta[[animal_col]]

############################################################
## Extract sex
############################################################

sex_col <- intersect(
    c("animal sex:ch1", "sex", "gender"),
    names(meta)
)[1]

if (is.na(sex_col)) {
    stop("Required metadata column not found for sex")
}

meta$Sex <- meta[[sex_col]]

############################################################
## Extract vector group
############################################################

vector_col <- intersect(
    c("treatment:ch1", "vector", "group", "treatment", "vaccine"),
    names(meta)
)[1]

if (is.na(vector_col)) {
    stop("Required metadata column not found for vaccine vector")
}

meta$Vector <- meta[[vector_col]]

############################################################
## Extract day
############################################################

day_col <- intersect(
    c("time point:ch1", "day", "time_point", "time"),
    names(meta)
)[1]

if (is.na(day_col)) {
    stop("Required metadata column not found for time point")
}

meta$Day <- meta[[day_col]]

############################################################
## Clean day variable
############################################################

meta$Day <- gsub(" ", "", meta$Day)

############################################################
## Create phase variable
############################################################
############################################################
## Create numeric day variable
############################################################

meta <- meta %>%
    mutate(
        DayNumeric = as.numeric(
            gsub("D", "", Day)
        )
    )
meta <- meta %>%
    mutate(
        Phase = case_when(

            DayNumeric %in% c(0, 3, 7, 14)
            ~ "Prime",

            DayNumeric %in% c(98, 101, 105, 112)
            ~ "Boost1",

            DayNumeric == 560
            ~ "Memory",

            DayNumeric %in% c(743, 746, 750, 757)
            ~ "Boost2",

            DayNumeric %in% c(1177, 1191, 1219)
            ~ "Durability",

            TRUE
            ~ NA_character_
        )
    )
############################################################
## Convert to factors
############################################################

meta$Phase <- factor(
    meta$Phase,
    levels = c(
        "Prime",
        "Boost1",
        "Memory",
        "Boost2",
        "Durability"
    )
)

meta$Sex <- factor(meta$Sex)

meta$Vector <- factor(meta$Vector)

############################################################
## Keep only required columns
############################################################

metadata_clean <- meta %>%
    select(
        SampleID,
        Animal,
        Sex,
        Vector,
        Day,
        DayNumeric,
        Phase
    )

############################################################
## Save metadata
############################################################

write_csv(
    metadata_clean,
    file.path(
        PROCESSED_DATA,
        "metadata_clean.csv"
    )
)

############################################################
## Summary
############################################################

cat("\n")
cat("=====================================\n")
cat("Metadata preparation complete\n")
cat("=====================================\n")

cat("\nSamples:\n")
print(nrow(metadata_clean))

cat("\nAnimals:\n")
print(length(unique(metadata_clean$Animal)))

cat("\nSex distribution:\n")
print(table(metadata_clean$Sex))

cat("\nVector distribution:\n")
print(table(metadata_clean$Vector))

cat("\nPhase distribution:\n")
print(table(metadata_clean$Phase))

cat("\nOutput file:\n")
print(
    file.path(
        PROCESSED_DATA,
        "metadata_clean.csv"
    )
)