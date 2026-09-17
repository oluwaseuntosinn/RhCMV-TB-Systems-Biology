###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 05_Module_Preservation.R
## Purpose : WGCNA Module Preservation between Wild-Type (WT)
##           and Delta-pp71 (Δpp71) Vaccination Vectors
###############################################################

rm(list = ls())

###############################################################
## Load Project Setup
###############################################################

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(data.table)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

###############################################################
## Directories
###############################################################

GEO <- DATASET_ID

OUTPUT_DIR <- file.path(
  RESULTS_DIR,
  "05_Module_Preservation",
  GEO
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
## SECTION 1: Load WGCNA Network
###############################################################

net <- readRDS(
  file.path(
    WGCNA_DIR,
    "WGCNA_Network.rds"
  )
)

###############################################################
## SECTION 3: Load Expression Matrix and Metadata
###############################################################

## Check for WGCNA_Input.rds or fallback paths
wgcna_input_file <- file.path(WGCNA_DIR, "WGCNA_Input.rds")

if (file.exists(wgcna_input_file)) {
  datExpr_all <- readRDS(wgcna_input_file)
  ## datExpr is samples x genes; transpose to genes x samples
  expr_all <- as.data.frame(t(datExpr_all))
} else if (file.exists(file.path(WGCNA_DIR, "WGCNA_Input", "WGCNA_CountMatrix.rds"))) {
  expr_all <- as.data.frame(readRDS(file.path(WGCNA_DIR, "WGCNA_Input", "WGCNA_CountMatrix.rds")))
} else if (file.exists(file.path(WGCNA_DIR, "vst.rds"))) {
  vsd <- readRDS(file.path(WGCNA_DIR, "vst.rds"))
  expr_all <- as.data.frame(SummarizedExperiment::assay(vsd))
} else {
  stop("Expression input matrix not found in WGCNA directory.")
}

metadata <- read.csv(
  file.path(
    PROCESSED_DATA,
    "metadata_clean.csv"
  ),
  stringsAsFactors = FALSE
)

###############################################################
## SECTION 4: Split Vaccination Groups (WT vs Δpp71)
###############################################################

wt_samples <- metadata$SampleID[
  grepl("68-1 RhCMV/TB-6Ag", metadata$Vector, fixed = TRUE) &
    !grepl("pp71|delta", metadata$Vector, ignore.case = TRUE)
]

delta_samples <- metadata$SampleID[
  grepl("delta|pp71", metadata$Vector, ignore.case = TRUE)
]

## Match available sample columns in expression matrix
wt_samples <- intersect(wt_samples, colnames(expr_all))
delta_samples <- intersect(delta_samples, colnames(expr_all))

cat("\nReference (WT) samples:", length(wt_samples), "\n")
cat("Test (delta_pp71) samples:", length(delta_samples), "\n")

expr.discovery <- expr_all[, wt_samples]
expr.validation <- expr_all[, delta_samples]

###############################################################
## SECTION 5: Common Genes & Module Alignment
###############################################################

common.genes <- names(net$colors)

expr.discovery.common <- expr.discovery[common.genes, ]
expr.validation.common <- expr.validation[common.genes, ]

moduleColors <- net$colors[
  match(
    common.genes,
    names(net$colors)
  )
]

###############################################################
## Remove Zero Variance Genes
###############################################################

gsg.discovery <- goodSamplesGenes(
  t(expr.discovery.common),
  verbose = 3
)

gsg.validation <- goodSamplesGenes(
  t(expr.validation.common),
  verbose = 3
)

goodGenes <- gsg.discovery$goodGenes & gsg.validation$goodGenes

expr.discovery.common <- expr.discovery.common[goodGenes, ]
expr.validation.common <- expr.validation.common[goodGenes, ]
moduleColors <- moduleColors[goodGenes]

cat(
  "\nGenes retained for preservation analysis:",
  sum(goodGenes),
  "\n"
)

###############################################################
## SECTION 6: Module Input Summary
###############################################################

module.summary <- data.frame(
  Module = unique(moduleColors)
)

module.summary$DiscoverySize <- sapply(
  module.summary$Module,
  function(m) sum(net$colors == m)
)

module.summary$PreservationInputSize <- sapply(
  module.summary$Module,
  function(m) sum(moduleColors == m)
)

module.summary <- module.summary[
  order(module.summary$PreservationInputSize, decreasing = TRUE),
]

write.csv(
  module.summary,
  file.path(
    OUTPUT_DIR,
    "Module_Preservation_Input_Summary.csv"
  ),
  row.names = FALSE
)

###############################################################
## Build multiExpr and multiColor
###############################################################

multiExpr <- list()

multiExpr[[1]] <- list(
  data = as.data.frame(t(expr.discovery.common))
)

multiExpr[[2]] <- list(
  data = as.data.frame(t(expr.validation.common))
)

names(multiExpr) <- c(
  "Discovery",
  "Validation"
)

multiColor <- list()
multiColor[[1]] <- moduleColors
names(multiColor) <- "Discovery"

###############################################################
## SECTION 7: Module Preservation Calculation
###############################################################

set.seed(123)

pres <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  networkType = "signed",
  nPermutations = 500,
  randomSeed = 123,
  verbose = 3
)

###############################################################
## SECTION 8: Inspect and Extract Zsummary & Observed Statistics
###############################################################

cat("\nAvailable Preservation Comparison Slots:\n")
z_slot_names <- names(pres$preservation$Z$ref.Discovery)
print(z_slot_names)

## Match the validation slot dynamically
slot_idx <- grep("Validation|delta", z_slot_names, ignore.case = TRUE)
if (length(slot_idx) > 0) {
  slot_target <- slot_idx[1]
} else {
  slot_target <- 2
}

z_stats <- as.data.frame(pres$preservation$Z$ref.Discovery[[slot_target]])
obs_stats <- as.data.frame(pres$preservation$observed$ref.Discovery[[slot_target]])

mod_names <- rownames(z_stats)

pres.stats <- data.frame(
  Module = mod_names,
  moduleSize = if ("moduleSize" %in% colnames(z_stats)) z_stats$moduleSize else (if ("moduleSize" %in% colnames(obs_stats)) obs_stats$moduleSize else NA_real_),
  medianRank.pres = if ("medianRank.pres" %in% colnames(obs_stats)) obs_stats$medianRank.pres else (if ("medianRank.pres" %in% colnames(z_stats)) z_stats$medianRank.pres else NA_real_),
  Zsummary.pres = if ("Zsummary.pres" %in% colnames(z_stats)) z_stats$Zsummary.pres else NA_real_,
  stringsAsFactors = FALSE
)

###############################################################
## SECTIONS 9 & 10: Manuscript Preservation Table & Classes
###############################################################

pres.stats <- pres.stats %>%
  mutate(
    moduleSize = as.numeric(moduleSize),
    medianRank.pres = as.numeric(medianRank.pres),
    Zsummary.pres = as.numeric(Zsummary.pres),
    PreservationClass = case_when(
      Zsummary.pres > 10 ~ "Strong",
      Zsummary.pres >= 2 ~ "Moderate",
      TRUE ~ "Weak"
    )
  ) %>%
  arrange(desc(Zsummary.pres))

###############################################################
## SECTION 11: Add Recall Scores
###############################################################

recall_file <- file.path(
  RESULTS_DIR,
  "04_Module_Eigengene_Trajectories",
  GEO,
  "Module_Recall_Scores_Summary.csv"
)

if (file.exists(recall_file)) {
 recall <- read.csv(
  recall_file,
  stringsAsFactors = FALSE
) %>%
  group_by(Module) %>%
  summarise(
  Mean_RecallScore =
    mean(Mean_RecallScore, na.rm = TRUE),

  RecallClass =
    dplyr::first(RecallClass),

  .groups = "drop"
)
  
  ## Ensure module names align (with or without 'ME' prefix)
  pres.stats <- pres.stats %>%
    mutate(Module_Raw = gsub("^ME", "", Module))
    
  recall <- recall %>%
    mutate(Module_Raw = gsub("^ME", "", Module))
    
  pres.stats <- left_join(
    pres.stats,
    recall,
    by = "Module_Raw"
  ) %>%
    mutate(
      Module = ifelse(!is.na(Module.x), Module.x, Module.y)
    ) %>%
    dplyr::select(-any_of(c("Module.x", "Module.y", "Module_Raw")))
}

write.csv(
  pres.stats,
  file.path(
    OUTPUT_DIR,
    "Module_Preservation_Zsummary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  pres.stats,
  file.path(
    OUTPUT_DIR,
    "Module_Preservation_Zsummary.rds"
  )
)

###############################################################
## SECTION 12: Publication-Quality Preservation Plot
###############################################################

plot_preservation <- function() {

  ## Filter out grey and gold benchmark modules for clearer plot scaling
  plot_df <- pres.stats %>%
    filter(!tolower(Module) %in% c("grey", "megrey", "gold", "megold"))

  y_max <- max(c(12, max(plot_df$Zsummary.pres, na.rm = TRUE) * 1.15), na.rm = TRUE)
  y_min <- min(c(-0.5, min(plot_df$Zsummary.pres, na.rm = TRUE) - 1), na.rm = TRUE)

  par(
    mar = c(5, 6, 4, 2),
    cex.lab = 1.4,
    cex.axis = 1.2,
    font.lab = 2
  )

  plot(
    plot_df$moduleSize,
    plot_df$Zsummary.pres,
    pch = 21,
    bg = gsub("^ME", "", plot_df$Module),
    col = "black",
    lwd = 1.2,
    cex = 2.2,
    ylim = c(y_min, y_max),
    xlab = "Module Size (Genes)",
    ylab = expression(paste("Preservation Statistic (", italic(Z)[summary], ")")),
    main = "Preservation of RhCMV-TB Co-expression Modules in \u0394pp71 Vector",
    cex.main = 1.3,
    font.main = 2
  )

  ###########################################################
  ## Preservation Thresholds
  ###########################################################

  abline(
    h = 2,
    col = "red",
    lty = 2,
    lwd = 2
  )

  abline(
    h = 10,
    col = "darkgreen",
    lty = 2,
    lwd = 2
  )

  ###########################################################
  ## Threshold Labels
  ###########################################################

  text(
    x = max(plot_df$moduleSize, na.rm = TRUE) * 0.98,
    y = 2.4,
    labels = "Moderate (Z = 2)",
    col = "red",
    pos = 2,
    cex = 0.9,
    font = 2
  )

  text(
    x = max(plot_df$moduleSize, na.rm = TRUE) * 0.98,
    y = 10.4,
    labels = "Strong (Z = 10)",
    col = "darkgreen",
    pos = 2,
    cex = 0.9,
    font = 2
  )

  ###########################################################
  ## Module Labels
  ###########################################################

  label_colors <- gsub("^ME", "", plot_df$Module)

  text(
    plot_df$moduleSize,
    plot_df$Zsummary.pres,
    labels = plot_df$Module,
    pos = 3,
    offset = 0.5,
    cex = 0.95,
    font = 2,
    col = ifelse(
      label_colors %in% c("lightcyan", "yellow", "greenyellow", "white"),
      "black",
      label_colors
    )
  )
}

###############################################################
## Export PDF & PNG
###############################################################

pdf(
  file.path(
    OUTPUT_DIR,
    "Figure_Module_Preservation_WT_vs_Delta.pdf"
  ),
  width = 9,
  height = 7
)
plot_preservation()
dev.off()

png(
  file.path(
    OUTPUT_DIR,
    "Figure_Module_Preservation_WT_vs_Delta.png"
  ),
  width = 3000,
  height = 2400,
  res = 300
)
plot_preservation()
dev.off()

cat("\n=============================================\n")
cat(" Module preservation analysis completed.\n")
cat(" Results folder:\n")
cat(OUTPUT_DIR, "\n")
cat("=============================================\n")