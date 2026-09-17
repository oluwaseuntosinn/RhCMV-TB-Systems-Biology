###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 08_Hub_Genes.R
## Purpose :
##   Identify hub genes from modules showing
##   significant Vector × Phase interactions
###############################################################

rm(list = ls())

source("scripts/00_setup.R")

suppressPackageStartupMessages({

  library(tidyverse)
  library(WGCNA)

})

options(stringsAsFactors = FALSE)

###############################################################
## Directories
###############################################################

results_dir <- file.path(
  RESULTS_DIR,
  "08_Hub_Genes",
  DATASET_ID
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
## Load Significant Modules
###############################################################

sig_modules <- read.csv(
  file.path(
    RESULTS_DIR,
    "06_Vector_Dynamics",
    DATASET_ID,
    "Significant_VectorPhase_Modules.csv"
  ),
  stringsAsFactors = FALSE
)

###############################################################
## Load WGCNA Objects / Files
###############################################################

## 1. Expression Matrix
if (file.exists(file.path(WGCNA_DIR, "WGCNA_Input.rds"))) {
  datExpr <- readRDS(file.path(WGCNA_DIR, "WGCNA_Input.rds"))
} else {
  stop("Cannot find WGCNA_Input.rds in: ", WGCNA_DIR)
}

## 2. Module Eigengenes
if (file.exists(file.path(WGCNA_DIR, "Module_Eigengenes.rds"))) {
  MEs <- readRDS(file.path(WGCNA_DIR, "Module_Eigengenes.rds"))
} else if (file.exists(file.path(WGCNA_DIR, "Module_Eigengenes.csv"))) {
  MEs <- read.csv(file.path(WGCNA_DIR, "Module_Eigengenes.csv"), row.names = 1, check.names = FALSE)
} else {
  stop("Cannot find Module_Eigengenes in: ", WGCNA_DIR)
}

## 3. Module Assignments
if (file.exists(file.path(WGCNA_DIR, "Gene_Module_Assignment.csv"))) {
  gene_module_df <- read.csv(file.path(WGCNA_DIR, "Gene_Module_Assignment.csv"), stringsAsFactors = FALSE)
  moduleColors <- gene_module_df$Module
  names(moduleColors) <- gene_module_df$Gene
} else if (file.exists(file.path(WGCNA_DIR, "WGCNA_Network.rds"))) {
  net <- readRDS(file.path(WGCNA_DIR, "WGCNA_Network.rds"))
  moduleColors <- net$colors
} else {
  stop("Cannot find Gene_Module_Assignment.csv or WGCNA_Network.rds in: ", WGCNA_DIR)
}

## 4. Optional HGNC Annotation Lookup
annot_file <- file.path(ROOT_DIR, "GSE273911_1.norm_matrix_HGNC.csv")
gene_lookup <- NULL

if (file.exists(annot_file)) {
  cat("\nLoading HGNC annotation from:", annot_file, "\n")
  annot <- read.csv(annot_file, check.names = FALSE, stringsAsFactors = FALSE)
  gene_lookup <- data.frame(
    ENSMMUG = annot[[2]],
    SYMBOL  = annot[[3]],
    stringsAsFactors = FALSE
  )
  gene_lookup <- unique(gene_lookup)
  gene_lookup <- gene_lookup[!is.na(gene_lookup$ENSMMUG) & gene_lookup$ENSMMUG != "", ]
}

###############################################################
## Calculate kME (Module Membership)
###############################################################

kMEtable <- as.data.frame(
  signedKME(
    datExpr,
    MEs
  )
)

###############################################################
## Gene names
###############################################################

gene_names <- colnames(datExpr)
moduleColors <- moduleColors[gene_names]

hub_summary <- list()

###############################################################
## Process each significant module
###############################################################

for(mod in sig_modules$Module){

  cat(
    "\n================================\n",
    mod,
    "\n================================\n"
  )

  colour <- gsub(
    "^ME",
    "",
    mod
  )

  ###########################################################
  ## Genes belonging to module
  ###########################################################

  idx <- moduleColors == colour

  module_genes <- gene_names[idx]

  ###########################################################
  ## kME column
  ###########################################################

  kme_col <- grep(
  colour,
  colnames(kMEtable),
  value = TRUE
)

if(length(kme_col) == 0){

  cat(
    "Skipping ",
    mod,
    ": no matching kME column found\n"
  )

  next
}

kme_col <- kme_col[1]

  ###########################################################
  ## Build hub table
  ###########################################################

  hub_df <- data.frame(
    Gene = module_genes,
    kME = kMEtable[idx, kme_col],
    stringsAsFactors = FALSE
  )

  if (!is.null(gene_lookup)) {
    hub_df <- hub_df %>%
      left_join(gene_lookup, by = c("Gene" = "ENSMMUG")) %>%
      relocate(SYMBOL, .after = Gene)
  }

  hub_df <- hub_df %>%
    arrange(
      desc(abs(kME))
    )

  ###########################################################
  ## Save full hub list
  ###########################################################

  write.csv(
    hub_df,
    file.path(
      results_dir,
      paste0(
        mod,
        "_HubGenes_All.csv"
      )
    ),
    row.names = FALSE
  )

  ###########################################################
  ## Top 50 hubs
  ###########################################################

  top50 <- hub_df %>%
    slice_head(
      n = 50
    )

  write.csv(
    top50,
    file.path(
      results_dir,
      paste0(
        mod,
        "_Top50_HubGenes.csv"
      )
    ),
    row.names = FALSE
  )

  ###########################################################
  ## Summary
  ###########################################################

  top_gene_symbol <- if ("SYMBOL" %in% colnames(top50)) top50$SYMBOL[1] else NA_character_

  hub_summary[[mod]] <- data.frame(
    Module = mod,
    Genes = nrow(hub_df),
    TopHubGene = top50$Gene[1],
    TopHubSymbol = top_gene_symbol,
    TopHub_kME = top50$kME[1],
    stringsAsFactors = FALSE
  )

}

###############################################################
## Save summary
###############################################################

hub_summary_df <- bind_rows(
  hub_summary
)

write.csv(

  hub_summary_df,

  file.path(
    results_dir,
    "HubGene_Summary.csv"
  ),

  row.names = FALSE

)

###############################################################
## Completion
###############################################################

cat(
  "\n====================================\n",
  "Hub gene analysis complete\n",
  "====================================\n"
)