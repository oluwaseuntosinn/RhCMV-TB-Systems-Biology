###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 00_setup.R
## Purpose : Define project-wide paths, options, and package list
## Author  : Oluwaseun T. Oluwatosin
###############################################################

## ============================================================
## Resolve project root
## ============================================================

if (!exists("ROOT_DIR", inherits = FALSE) || !nzchar(ROOT_DIR)) {

    cwd <- normalizePath(getwd(),
                         winslash = "/",
                         mustWork = FALSE)

    if (dir.exists(file.path(cwd, "scripts"))) {
        ROOT_DIR <- cwd
    } else {
        ROOT_DIR <- normalizePath(cwd,
                                  winslash = "/",
                                  mustWork = FALSE)
    }
}

ROOT_DIR <- normalizePath(ROOT_DIR,
                          winslash = "/",
                          mustWork = FALSE)

## ============================================================
## Directory structure
## ============================================================

DATA_DIR        <- file.path(ROOT_DIR, "data")

RAW_DATA        <- file.path(DATA_DIR, "raw")
PROCESSED_DATA  <- file.path(DATA_DIR, "processed")
METADATA_DIR    <- file.path(DATA_DIR, "metadata")

RESULTS_DIR     <- file.path(ROOT_DIR, "results")
FIGURES_DIR     <- file.path(ROOT_DIR, "figures")
TABLES_DIR      <- file.path(ROOT_DIR, "tables")
LOGS_DIR        <- file.path(ROOT_DIR, "logs")

MANUSCRIPT_DIR  <- file.path(ROOT_DIR, "manuscript")
SCRIPTS_DIR     <- file.path(ROOT_DIR, "scripts")

NETWORK_DIR     <- file.path(RESULTS_DIR, "network")
WGCNA_DIR       <- file.path(RESULTS_DIR, "WGCNA")
GSVA_DIR        <- file.path(RESULTS_DIR, "GSVA")
ENRICHMENT_DIR  <- file.path(RESULTS_DIR, "enrichment")
CYTOSCAPE_DIR   <- file.path(RESULTS_DIR, "cytoscape")

dirs <- c(
    DATA_DIR,
    RAW_DATA,
    PROCESSED_DATA,
    METADATA_DIR,
    RESULTS_DIR,
    FIGURES_DIR,
    TABLES_DIR,
    LOGS_DIR,
    MANUSCRIPT_DIR,
    SCRIPTS_DIR,
    NETWORK_DIR,
    WGCNA_DIR,
    GSVA_DIR,
    ENRICHMENT_DIR,
    CYTOSCAPE_DIR
)

for (d in dirs) {
    if (!dir.exists(d)) {
        dir.create(d,
                   recursive = TRUE,
                   showWarnings = FALSE)
    }
}

## ============================================================
## Global options
## ============================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(timeout = 600)
options(repos = c(CRAN = "https://cloud.r-project.org"))

set.seed(12345)

enable_installation <- FALSE

## ============================================================
## Dataset information
## ============================================================

PROJECT_NAME <- "RhCMV-TB Systems Vaccinology"

DATASET_ID   <- "GSE273911"

PROJECT_TITLE <- paste(
    "Preservation and Recall of Vaccine-Induced Immune Networks",
    "Following Sequential RhCMV-TB Immunization"
)

## ============================================================
## Analysis objectives
## ============================================================

PRIMARY_AIM <- paste(
    "Identify preserved co-expression modules associated",
    "with durable RhCMV-TB vaccine responses"
)

SECONDARY_AIMS <- c(
    "Module-trait association analysis",
    "Longitudinal eigengene trajectory analysis",
    "Module preservation across vaccination phases",
    "Hub gene identification",
    "Network rewiring analysis",
    "GSVA pathway validation",
    "Functional enrichment analysis",
    "Cytoscape network visualization"
)

## ============================================================
## Package inventory
## ============================================================

all_packages <- c(

    ## Core
    "tidyverse",
    "data.table",
    "readr",
    "readxl",
    "stringr",
    "glue",
    "janitor",

    ## Visualization
    "ggplot2",
    "ggrepel",
    "patchwork",
    "cowplot",
    "pheatmap",
    "viridis",
    "RColorBrewer",
    "ComplexHeatmap",
    "circlize",

    ## Statistics
    "limma",
    "lme4",
    "lmerTest",
    "broom",

    ## RNA-seq
    "DESeq2",
    "apeglm",

    ## WGCNA
    "WGCNA",
    "dynamicTreeCut",
    "fastcluster",

    ## Enrichment
    "clusterProfiler",
    "enrichplot",
    "DOSE",

    ## Network analysis
    "igraph",
    "ggraph",
    "tidygraph",
    "STRINGdb",

    ## GSVA
    "GSVA",
    "GSEABase",

    ## Annotation
    "AnnotationDbi",
    "org.Hs.eg.db",

    ## Utilities
    "matrixStats",
    "BiocParallel"
)

## ============================================================
## Package installer
## ============================================================

if (enable_installation) {

    cran_packages <- all_packages[
        !all_packages %in%
            rownames(installed.packages())
    ]

    if (length(cran_packages) > 0) {
        install.packages(cran_packages)
    }
}

############################################################
## Load required packages
############################################################

required_packages <- c(
    "tidyverse",
    "data.table",
    "readxl",
    "DESeq2",
    "WGCNA"
)

invisible(
    lapply(required_packages,
           library,
           character.only = TRUE)
)

## ============================================================
## WGCNA settings
## ============================================================

WGCNA::allowWGCNAThreads()

WGCNA_OPTIONS <- list(
    networkType = "signed",
    corType = "bicor",
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    deepSplit = 2
)

## ============================================================
## Expected project outputs
## ============================================================

EXPECTED_FIGURES <- c(
    "Figure1_WGCNA_Dendrogram",
    "Figure2_ModuleTrait_Heatmap",
    "Figure3_Eigengene_Trajectories",
    "Figure4_Module_Preservation",
    "Figure5_HubGene_Network",
    "Figure6_GSVA_Validation"
)

## ============================================================
## Startup message
## ============================================================

cat("\n")
cat("============================================\n")
cat(" RhCMV-TB Systems Vaccinology Project\n")
cat(" Dataset :", DATASET_ID, "\n")
cat("============================================\n")
cat(" Root directory :", ROOT_DIR, "\n")
cat("============================================\n\n")