###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 07_Module_Enrichment.R
## Purpose :
## Functional enrichment of modules showing
## significant Vector × Phase interactions
###############################################################

rm(list = ls())

source("scripts/00_setup.R")

suppressPackageStartupMessages({

  library(tidyverse)

  library(clusterProfiler)

  library(org.Hs.eg.db)

})

###############################################################
## Directories
###############################################################

results_dir <- file.path(
  RESULTS_DIR,
  "07_Module_Enrichment",
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

sig_modules_file <- file.path(
  RESULTS_DIR,
  "06_Vector_Dynamics",
  DATASET_ID,
  "Significant_VectorPhase_Modules.csv"
)

if (!file.exists(sig_modules_file)) {
  stop("Significant modules file not found.")
}

sig_modules <- read.csv(
  sig_modules_file,
  stringsAsFactors = FALSE
)

###############################################################
## Load Module Assignments
###############################################################

module_assignments_file <- file.path(
  WGCNA_DIR,
  "Gene_Module_Assignment.csv"
)

if (!file.exists(module_assignments_file)) {

  module_assignments_file <- file.path(
    RESULTS_DIR,
    "03_WGCNA",
    DATASET_ID,
    "Gene_Module_Assignment.csv"
  )

}

if (!file.exists(module_assignments_file)) {
  stop("Gene_Module_Assignment.csv not found.")
}

module_assignments <- read.csv(
  module_assignments_file,
  stringsAsFactors = FALSE
)

###############################################################
## Load HGNC Annotation File
###############################################################

annot_file <- "C:/Users/oluwa/Documents/8th TB Forum/GSE273911_1.norm_matrix_HGNC.csv"

if (!file.exists(annot_file)) {
  stop("HGNC annotation file not found.")
}

cat("\nLoading HGNC annotation...\n")

annot <- read.csv(
  annot_file,
  check.names = FALSE
)

###############################################################
## Build ENSMMUG -> HGNC lookup
###############################################################

gene_lookup <- data.frame(
  ENSMMUG = annot[[2]],
  SYMBOL  = annot[[3]],
  stringsAsFactors = FALSE
)

gene_lookup <- unique(gene_lookup)

gene_lookup <- gene_lookup[
  !is.na(gene_lookup$ENSMMUG) &
  !is.na(gene_lookup$SYMBOL) &
  gene_lookup$ENSMMUG != "" &
  gene_lookup$SYMBOL  != "",
]


cat(
  "Mapped genes:",
  nrow(gene_lookup),
  "\n"
)

###############################################################
## Helper
###############################################################

convert_entrez_to_symbol <- function(id_string){

  ids <- unlist(
    strsplit(id_string, "/")
  )

  genes <- tryCatch(
    bitr(
      ids,
      fromType = "ENTREZID",
      toType = "SYMBOL",
      OrgDb = org.Hs.eg.db
    ),
    error = function(e) NULL
  )

  if(is.null(genes) || nrow(genes)==0){
    return(id_string)
  }

  paste(
    unique(genes$SYMBOL),
    collapse = "/"
  )
}

###############################################################
## Enrichment
###############################################################

summary_list <- list()

for(mod in sig_modules$Module){

  cat(
    "\n=====================================\n",
    "Running:",
    mod,
    "\n=====================================\n"
  )

  colour <- gsub("^ME","",mod)

  #############################################################
  ## Module genes
  #############################################################

  module_genes <- module_assignments %>%
    filter(Module == colour) %>%
    pull(Gene) %>%
    unique()

  #############################################################
  ## Convert ENSMMUG -> HGNC
  #############################################################

  symbols <- gene_lookup %>%
    filter(
      ENSMMUG %in% module_genes
    ) %>%
    pull(SYMBOL) %>%
    unique()

  symbols <- symbols[
    symbols != ""
  ]

  cat(
    "HGNC symbols:",
    length(symbols),
    "\n"
  )

  #############################################################
  ## HGNC -> ENTREZ
  #############################################################

  gene_map <- tryCatch(

    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ),

    error = function(e) NULL

  )

  if(is.null(gene_map) || nrow(gene_map)==0){

    cat(
      "No ENTREZ IDs found.\n"
    )

    next
  }

  entrez_genes <- unique(
    gene_map$ENTREZID
  )

  cat(
    "ENTREZ IDs:",
    length(entrez_genes),
    "\n"
  )

  #############################################################
  ## GO BP
  #############################################################

  go_bp <- tryCatch(

    enrichGO(
      gene = entrez_genes,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    ),

    error = function(e) NULL

  )

  #############################################################
  ## GO MF
  #############################################################

  go_mf <- tryCatch(

    enrichGO(
      gene = entrez_genes,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    ),

    error = function(e) NULL

  )

  #############################################################
  ## KEGG
  #############################################################

  kegg <- tryCatch(

    enrichKEGG(
      gene = entrez_genes,
      organism = "hsa",
      pvalueCutoff = 0.05
    ),

    error = function(e) NULL

  )

  #############################################################
  ## Convert to data.frames
  #############################################################

  go_bp_df <- if(!is.null(go_bp))
    as.data.frame(go_bp)
  else
    data.frame()

  go_mf_df <- if(!is.null(go_mf))
    as.data.frame(go_mf)
  else
    data.frame()

  kegg_df <- if(!is.null(kegg))
    as.data.frame(kegg)
  else
    data.frame()

  #############################################################
  ## Save GO BP
  #############################################################

  if(nrow(go_bp_df)>0){

    go_bp_df$geneID <- sapply(
      go_bp_df$geneID,
      convert_entrez_to_symbol
    )

    write.csv(
      go_bp_df,
      file.path(
        results_dir,
        paste0(mod,"_GO_BP.csv")
      ),
      row.names = FALSE
    )

    pdf(
      file.path(
        results_dir,
        paste0(mod,"_GO_BP_Dotplot.pdf")
      ),
      width = 10,
      height = 7
    )

    print(
      dotplot(
        go_bp,
        showCategory = 15
      )
    )

    dev.off()
  }

  #############################################################
  ## Save GO MF
  #############################################################

  if(nrow(go_mf_df)>0){

    go_mf_df$geneID <- sapply(
      go_mf_df$geneID,
      convert_entrez_to_symbol
    )

    write.csv(
      go_mf_df,
      file.path(
        results_dir,
        paste0(mod,"_GO_MF.csv")
      ),
      row.names = FALSE
    )
  }

  #############################################################
  ## Save KEGG
  #############################################################

  if(nrow(kegg_df)>0){

    kegg_df$geneID <- sapply(
      kegg_df$geneID,
      convert_entrez_to_symbol
    )

    write.csv(
      kegg_df,
      file.path(
        results_dir,
        paste0(mod,"_KEGG.csv")
      ),
      row.names = FALSE
    )

    pdf(
      file.path(
        results_dir,
        paste0(mod,"_KEGG_Dotplot.pdf")
      ),
      width = 10,
      height = 7
    )

    print(
      dotplot(
        kegg,
        showCategory = 15
      )
    )

    dev.off()
  }

  #############################################################
  ## Summary
  #############################################################

  summary_list[[mod]] <- data.frame(

    Module = mod,

    ModuleGenes = length(module_genes),

    HGNC_Genes = length(symbols),

    EntrezGenes = length(entrez_genes),

    GO_BP_Terms = nrow(go_bp_df),

    GO_MF_Terms = nrow(go_mf_df),

    KEGG_Pathways = nrow(kegg_df)

  )
}

###############################################################
## Save Summary
###############################################################

summary_table <- bind_rows(
  summary_list
)

write.csv(
  summary_table,
  file.path(
    results_dir,
    "Module_Enrichment_Summary.csv"
  ),
  row.names = FALSE
)

###############################################################
## Completion
###############################################################

cat(
  "\n====================================\n",
  "Module enrichment complete\n",
  "====================================\n"
)