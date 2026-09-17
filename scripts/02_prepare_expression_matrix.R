###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 02_prepare_expression_matrix.R
## Purpose : Prepare expression matrix and rename samples to GSM IDs
###############################################################

rm(list = ls())

source("scripts/00_setup.R")

library(tidyverse)
library(data.table)

############################################################
## Load cleaned metadata
############################################################

meta <- read_csv(
    file.path(
        PROCESSED_DATA,
        "metadata_clean.csv"
    ),
    show_col_types = FALSE
)

############################################################
## Load GEO phenodata
############################################################

geo_meta <- read_csv(
    file.path(
        METADATA_DIR,
        "phenodata.csv"
    ),
    show_col_types = FALSE
)

############################################################
## Load count matrix
############################################################

counts <- fread(
    file.path(
        RAW_DATA,
        "GSE273911_count_matrix.csv"
    )
)

cat("Count matrix dimensions:\n")
print(dim(counts))

############################################################
## Basic checks
############################################################

cat("\nSamples in count matrix:\n")
print(ncol(counts) - 1)

cat("\nSamples in GEO metadata:\n")
print(nrow(geo_meta))

cat("\nSamples in clean metadata:\n")
print(nrow(meta))

if ((ncol(counts) - 1) != nrow(geo_meta)) {
    stop(
        "Number of count-matrix samples does not match GEO metadata."
    )
}

############################################################
## Rename sample columns using GSM IDs
############################################################

new_names <- c(
    "GeneID",
    geo_meta$geo_accession
)

colnames(counts) <- new_names

############################################################
## Verify metadata match
############################################################

count_samples <- colnames(counts)[-1]

cat("\nChecking GSM match with metadata...\n")

perfect_match <- all(
    sort(count_samples) ==
    sort(meta$SampleID)
)

print(perfect_match)

if (!perfect_match) {

    missing_in_counts <- setdiff(
        meta$SampleID,
        count_samples
    )

    missing_in_meta <- setdiff(
        count_samples,
        meta$SampleID
    )

    cat("\nMissing in count matrix:\n")
    print(missing_in_counts)

    cat("\nMissing in metadata:\n")
    print(missing_in_meta)

    stop("Sample IDs do not match.")
}

############################################################
## Check for duplicated sample IDs
############################################################

if (anyDuplicated(count_samples) > 0) {
    stop("Duplicate GSM IDs detected.")
}

############################################################
## Check gene column
############################################################

cat("\nGene column name:\n")
print(colnames(counts)[1])

cat("\nNumber of genes:\n")
print(nrow(counts))

############################################################
## Save expression matrix
############################################################

fwrite(
    counts,
    file.path(
        PROCESSED_DATA,
        "counts_GSM.csv"
    )
)

############################################################
## Summary
############################################################

cat("\n=====================================\n")
cat("Expression matrix preparation complete\n")
cat("=====================================\n")

cat("\nGenes:\n")
print(nrow(counts))

cat("\nSamples:\n")
print(ncol(counts) - 1)

cat("\nOutput file:\n")
print(
    file.path(
        PROCESSED_DATA,
        "counts_GSM.csv"
    )
)

cat("\nMetadata match:\n")
print(perfect_match)