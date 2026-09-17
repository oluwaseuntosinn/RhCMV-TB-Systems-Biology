###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 03_WGCNA_Discovery.R
## Purpose : WGCNA Discovery Analysis
###############################################################

rm(list = ls())

###############################################################
## Load project setup
###############################################################

source("scripts/00_setup.R")

###############################################################
## Packages
###############################################################

options(stringsAsFactors = FALSE)

WGCNA::allowWGCNAThreads()

###############################################################
## Output directory
###############################################################

results_wgcna <- WGCNA_DIR

dir.create(
  results_wgcna,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
## Load metadata
###############################################################

metadata <-
  read_csv(
    file.path(
      PROCESSED_DATA,
      "metadata_clean.csv"
    ),
    show_col_types = FALSE
  ) %>%
  as.data.frame()

###############################################################
## Load count matrix
###############################################################

counts <-
  fread(
    file.path(
      PROCESSED_DATA,
      "counts_GSM.csv"
    )
  )

###############################################################
## Convert to matrix
###############################################################

gene_ids <- counts[[1]]

count_matrix <-
  as.matrix(
    counts[, -1]
  )

rownames(count_matrix) <- gene_ids

###############################################################
## Ensure integer counts
###############################################################

storage.mode(count_matrix) <- "integer"

###############################################################
## Match metadata order
###############################################################

count_samples <- colnames(count_matrix)

metadata <-
  metadata %>%
  filter(
    SampleID %in% count_samples
  )

metadata <-
  metadata[
    match(
      count_samples,
      metadata$SampleID
    ),
  ]

if (nrow(metadata) != length(count_samples) ||
    anyNA(metadata$SampleID) ||
    anyDuplicated(metadata$SampleID) > 0 ||
    !identical(metadata$SampleID, count_samples)) {
  stop("Count-matrix samples do not match metadata_clean.csv.")
}

metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata$SampleID

###############################################################
## Create DESeq2 object
###############################################################

dds <-
  DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ 1
  )

###############################################################
## Filter low-count genes
###############################################################

keep <-
  rowSums(
    counts(dds) >= 10
  ) >= 10

dds <- dds[keep, ]

cat(
  "\nGenes after filtering:",
  nrow(dds),
  "\n"
)

###############################################################
## Variance Stabilizing Transformation
###############################################################

vsd <-
  vst(
    dds,
    blind = TRUE
  )

expr_vst <-
  assay(vsd)

###############################################################
## Save DESeq2 outputs
###############################################################

saveRDS(
  dds,
  file.path(
    results_wgcna,
    "dds_filtered.rds"
  )
)

saveRDS(
  vsd,
  file.path(
    results_wgcna,
    "vst.rds"
  )
)

###############################################################
## Select most variable genes
###############################################################

gene_var <-
  apply(
    expr_vst,
    1,
    var
  )

top_n <- min(
  15000,
  length(gene_var)
)

top_genes <-
  names(
    sort(
      gene_var,
      decreasing = TRUE
    )[1:top_n]
  )

expr_wgcna <-
  expr_vst[
    top_genes,
  ]

cat(
  "\nGenes retained for WGCNA:",
  nrow(expr_wgcna),
  "\n"
)

###############################################################
## Transpose for WGCNA
###############################################################

datExpr <-
  t(
    expr_wgcna
  )

###############################################################
## Sample QC
###############################################################

gsg <-
  goodSamplesGenes(
    datExpr,
    verbose = 3
  )

if(!gsg$allOK){

  datExpr <-
    datExpr[
      gsg$goodSamples,
      gsg$goodGenes,
      drop = FALSE
    ]

}

metadata <-
  metadata[
    match(
      rownames(datExpr),
      metadata$SampleID
    ),
    ,
    drop = FALSE
  ]

rownames(metadata) <- metadata$SampleID

stopifnot(
  identical(
    rownames(datExpr),
    rownames(metadata)
  )
)

###############################################################
## Sample clustering
###############################################################

sampleTree <-
  hclust(
    dist(datExpr),
    method = "average"
  )

pdf(
  file.path(
    results_wgcna,
    "SampleClustering.pdf"
  ),
  width = 10,
  height = 6
)

plot(
  sampleTree,
  main = "Sample Clustering",
  xlab = "",
  sub = ""
)

dev.off()

###############################################################
## Save WGCNA input
###############################################################

saveRDS(
  datExpr,
  file.path(
    results_wgcna,
    "WGCNA_Input.rds"
  )
)

saveRDS(
  metadata,
  file.path(
    results_wgcna,
    "metadata_clean.rds"
  )
)

###############################################################
## Summary
###############################################################

cat("\n=====================================\n")
cat("WGCNA INPUT PREPARATION COMPLETE\n")
cat("=====================================\n")

cat(
  "\nSamples:",
  nrow(datExpr),
  "\n"
)

cat(
  "Genes:",
  ncol(datExpr),
  "\n"
)
###############################################################
## Soft-threshold Analysis
###############################################################

powers <-
c(
1:10,
seq(12,30,2)
)

sft <-
pickSoftThreshold(

datExpr,

powerVector = powers,

networkType = "signed",

verbose = 5

)


###############################################################
## Plot Soft Threshold Results
###############################################################

pdf(

file.path(
results_wgcna,
"SoftThreshold.pdf"
),

width = 11,
height = 5

)

par(mfrow=c(1,2))

plot(

sft$fitIndices[,1],

-sign(sft$fitIndices[,3])*
sft$fitIndices[,2],

type="b",

pch=19,

xlab="Soft Threshold (Power)",

ylab="Scale Free Topology Fit (R²)",

main="Scale Independence"

)

abline(
h=0.85,
col="red",
lty=2
)

plot(

sft$fitIndices[,1],

sft$fitIndices[,5],

type="b",

pch=19,

xlab="Soft Threshold (Power)",

ylab="Mean Connectivity",

main="Mean Connectivity"

)

dev.off()

###############################################################
## Automatically Select Soft Power
###############################################################

fitIndices <- sft$fitIndices

write.csv(

  fitIndices,

  file.path(

    results_wgcna,

    "SoftThreshold_Statistics.csv"

  ),

  row.names = FALSE

)

scaleFreeFit <-
  -sign(fitIndices[,3]) *
  fitIndices[,2]

candidateRows <-
  which(scaleFreeFit >= 0.85)

if(length(candidateRows) == 0){

  softPower <-

    fitIndices[
      which.max(scaleFreeFit),
      1
    ]

  message(
    "No power reached signed R² = 0.85.\n",
    "Using power with highest signed R²: ",
    softPower
  )

}else{

  softPower <-

    fitIndices[
      candidateRows[1],
      1
    ]

  message(
    "Automatically selected soft power: ",
    softPower
  )

}

###############################################################
## Determine Dataset Size
###############################################################

nGenes <- ncol(datExpr)

nSamples <- nrow(datExpr)

###############################################################
## Automatically Choose Module Size
###############################################################

if(nGenes < 5000){

minModuleSize <- 30

}else if(nGenes < 10000){

minModuleSize <- 40

}else{

minModuleSize <- 50

}

###############################################################
## Automatically Choose Block Size
###############################################################

maxBlockSize <-

min(

20000,

nGenes

)

###############################################################
## Report Selected Parameters
###############################################################

cat("\n")

cat("=====================================\n")

cat("WGCNA PARAMETERS\n")

cat("-------------------------------------\n")

cat("Samples          :",nSamples,"\n")

cat("Genes            :",nGenes,"\n")

cat("Soft Power       :",softPower,"\n")

cat("Min Module Size  :",minModuleSize,"\n")

cat("Merge Cut Height :",0.25,"\n")

cat("Max Block Size   :",maxBlockSize,"\n")

cat("=====================================\n")

###############################################################
## Construct WGCNA Network
###############################################################

gc()

net <-

  blockwiseModules(

    datExpr,

    power = softPower,

    networkType = "signed",

    TOMType = "signed",

    deepSplit = 2,

    minModuleSize = minModuleSize,

    reassignThreshold = 0,

    mergeCutHeight = 0.25,

    numericLabels = FALSE,

    pamRespectsDendro = FALSE,

    saveTOMs = FALSE,

    maxBlockSize = maxBlockSize,

    nThreads = 8,

    corType = "pearson",

    verbose = 3

  )

if(!exists("net")){

  stop(
    "blockwiseModules failed. Check memory limits and WGCNA settings."
  )

}
  ###############################################################
## Gene Dendrogram with Module Colours
###############################################################

pdf(
  file.path(
    results_wgcna,
    "Gene_Dendrogram_Modules.pdf"
  ),
  width = 12,
  height = 8
)

for(i in 1:length(net$dendrograms)){

  plotDendroAndColors(
    net$dendrograms[[i]],
    net$colors[net$blockGenes[[i]]],
    "Module",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Gene dendrogram - Block", i)
  )

}

dev.off()


###############################################################
## Save Module Assignment
###############################################################

moduleAssignment <-

data.frame(

Gene = colnames(datExpr),

Module = net$colors

)

write.csv(

moduleAssignment,

file.path(

results_wgcna,

"Gene_Module_Assignment.csv"

),

row.names=FALSE

)

###############################################################
## Save Module Sizes
###############################################################

moduleSizes <-

as.data.frame(

table(net$colors)

)

write.csv(

moduleSizes,

file.path(

results_wgcna,

"Module_Sizes.csv"

),

row.names=FALSE

)

###############################################################
## Save Network Parameters
###############################################################

networkParameters <-

data.frame(

Samples = nSamples,

Genes = nGenes,

SoftPower = softPower,

MinModuleSize = minModuleSize,

MergeCutHeight = 0.25,

MaxBlockSize = maxBlockSize

)

write.csv(

networkParameters,

file.path(

results_wgcna,

"Network_Parameters.csv"

),

row.names=FALSE
)

###############################################################
## End Part 1
###############################################################

###############################################################
## PART 2
## Module Eigengenes and Module-Trait Relationships
###############################################################

###############################################################
## Calculate Module Eigengenes
###############################################################

MEs0 <-
  moduleEigengenes(
    datExpr,
    colors = net$colors
  )$eigengenes

MEs <-
  orderMEs(MEs0)

write.csv(
  MEs,
  file.path(
    results_wgcna,
    "Module_Eigengenes.csv"
  )
)

###############################################################
## Number of Samples
###############################################################

nSamples <- nrow(datExpr)

###############################################################
## Build Trait Matrix
###############################################################

traits_core <-

  metadata %>%

  dplyr::transmute(

    SampleID,

    VectorNonAttenuated = as.numeric(
      !grepl("delta|pp71", Vector, ignore.case = TRUE)
    ),

    VectorAttenuated = as.numeric(
      grepl("delta|pp71", Vector, ignore.case = TRUE)
    ),

    SexMale = as.numeric(
      tolower(as.character(Sex)) == "male"
    ),

    DayNumeric,

    PhasePrime = as.numeric(Phase == "Prime"),

    PhaseBoost1 = as.numeric(Phase == "Boost1"),

    PhaseMemory = as.numeric(Phase == "Memory"),

    PhaseBoost2 = as.numeric(Phase == "Boost2"),

    PhaseDurability = as.numeric(Phase == "Durability")

  ) %>%
  as.data.frame()

###############################################################
## Match Sample Order
###############################################################

rownames(traits_core) <- traits_core$SampleID

traits_core$SampleID <- NULL

traits_core <-
  traits_core[
    match(
      rownames(datExpr),
      rownames(traits_core)
    ),
    ,
    drop = FALSE
  ]

stopifnot(
  identical(
    rownames(datExpr),
    rownames(traits_core)
  )
)
###############################################################
## Save Trait Matrix
###############################################################

write.csv(

  traits_core,

  file.path(

    results_wgcna,

    "Trait_Matrix_Core.csv"

  )

)

saveRDS(

  traits_core,

  file.path(

    results_wgcna,

    "Traits_Core.rds"

  )

)


###############################################################
## Module-Trait Correlations
## (CORE TRAITS ONLY)
###############################################################

moduleTraitCor <-
  cor(
    MEs,
    traits_core,
    use="pairwise.complete.obs"
  )

moduleTraitPvalue <-
  corPvalueStudent(
    moduleTraitCor,
    nSamples
  )

write.csv(
  moduleTraitCor,
  file.path(
    results_wgcna,
    "Module_Trait_Correlations.csv"
  )
)

write.csv(
  moduleTraitPvalue,
  file.path(
    results_wgcna,
    "Module_Trait_Pvalues.csv"
  )
)

###############################################################
## Heatmap Labels
###############################################################

pLabel <- ifelse(
  moduleTraitPvalue < 0.001,
  "<0.001",
  sprintf("%.3f", moduleTraitPvalue)
)

textMatrix <-
paste(
  sprintf("%.2f", moduleTraitCor),
  "\n(",
  pLabel,
  ")",
  sep=""
)

dim(textMatrix) <-
dim(moduleTraitCor)

###############################################################
## Module-Trait Heatmap 
###############################################################

pdf(
  file.path(
    results_wgcna,
    "Module_Trait_Heatmap.pdf"
  ),
  width = 12,
  height = 10
)

par(
  mar = c(8, 10, 3, 3)
)

labeledHeatmap(

  Matrix = moduleTraitCor,

  xLabels = colnames(traits_core),

  yLabels = names(MEs),

  ySymbols = names(MEs),

  colorLabels = FALSE,

  colors = blueWhiteRed(50),

  textMatrix = textMatrix,

  setStdMargins = FALSE,

  cex.text = 0.80,

  cex.lab = 1.00,

  zlim = c(-1, 1),

  main = "Module–trait relationships"

)

dev.off()

cat(
  "\nModule–trait heatmap saved to:\n",
  file.path(results_wgcna, "Module_Trait_Heatmap.pdf"),
  "\n"
)

###############################################################
## Eigengene Clustering
###############################################################

METree <-
hclust(
dist(t(MEs)),
method="average"
)

pdf(

file.path(
results_wgcna,
"Module_Eigengene_Clustering.pdf"
),

width=8,
height=6

)

plot(

METree,

main=
"Module Eigengene Clustering",

xlab="",

sub=""

)

dev.off()

###############################################################
## Eigengene Network
###############################################################

pdf(

file.path(
results_wgcna,
"Eigengene_Network.pdf"
),

width=9,
height=8

)

plotEigengeneNetworks(

MEs,

"",

plotDendrograms=TRUE,

marDendro=c(0,4,2,0),

marHeatmap=c(3,4,2,2),

xLabelsAngle=90

)

dev.off()

###############################################################
## Save Objects
###############################################################

saveRDS(
MEs,
file.path(
results_wgcna,
"Module_Eigengenes.rds"
)
)

saveRDS(
traits_core,
file.path(
results_wgcna,
"Traits_Core.rds"
)
)


saveRDS(
net,
file.path(
results_wgcna,
"WGCNA_Network.rds"
)
)

###############################################################
## END OF PART 2
###############################################################

###############################################################
## PART 3
## Gene Significance, Module Membership,
## Hub Genes and Cytoscape Export
###############################################################

###############################################################
## Module Membership (MM)
###############################################################

moduleColors <- net$colors

MMPvalue <-
as.data.frame(matrix(
NA,
nrow = ncol(datExpr),
ncol = ncol(MEs)
))

geneModuleMembership <-
as.data.frame(

cor(

datExpr,

MEs,

use = "pairwise.complete.obs"

)

)

names(geneModuleMembership) <-
paste0(
"MM_",
substring(names(MEs),3)
)

MMPvalue <-
as.data.frame(

corPvalueStudent(

as.matrix(geneModuleMembership),

nSamples

)

)

names(MMPvalue) <-
paste0(
"p.MM_",
substring(names(MEs),3)
)

###############################################################
## Gene Significance (GS)
###############################################################

geneTraitSignificance <-
as.data.frame(

cor(

datExpr,

traits_core,

use = "pairwise.complete.obs"

)

)

GSPvalue <-
as.data.frame(

corPvalueStudent(

as.matrix(geneTraitSignificance),

nSamples

)

)

names(geneTraitSignificance) <-
paste0(
"GS_",
colnames(traits_core)
)

names(GSPvalue) <-
paste0(
"p.GS_",
colnames(traits_core)
)

###############################################################
## Save MM and GS
###############################################################

write.csv(

geneModuleMembership,

file.path(

results_wgcna,

"ModuleMembership.csv"

)

)

write.csv(

geneTraitSignificance,

file.path(

results_wgcna,

"GeneSignificance.csv"

)

)

###############################################################
## Combine Statistics
###############################################################

geneInfo <-

data.frame(

Gene = colnames(datExpr),

Module = moduleColors

)

geneInfo <-

cbind(

geneInfo,

geneModuleMembership,

MMPvalue,

geneTraitSignificance,

GSPvalue

)

write.csv(

geneInfo,

file.path(

results_wgcna,

"Gene_Information.csv"

),

row.names = FALSE

)

###############################################################
## Intramodular Connectivity
###############################################################

modules <-
  sort(unique(moduleColors))

kWithin <-
  data.frame(
    kWithin = rep(NA_real_, ncol(datExpr)),
    Gene = colnames(datExpr),
    Module = moduleColors,
    stringsAsFactors = FALSE
  )

for (mod in modules) {

  inModule <- moduleColors == mod

  if (sum(inModule) < 2) {
    next
  }

  moduleAdjacency <-
    adjacency(
      datExpr[, inModule, drop = FALSE],
      power = softPower,
      type = "signed"
    )

  kWithin$kWithin[inModule] <-
  colSums(moduleAdjacency, na.rm = TRUE)
}

write.csv(

kWithin,

file.path(

results_wgcna,

"IntramodularConnectivity.csv"

),

row.names = FALSE

)

###############################################################
## Hub Gene Identification
###############################################################

hubGenes <-

  kWithin %>%

  group_by(Module) %>%

  arrange(

    dplyr::desc(kWithin),

    .by_group = TRUE

  ) %>%

  ungroup()

write.csv(

hubGenes,

file.path(

results_wgcna,

"HubGenes_All.csv"

),

row.names = FALSE

)

###############################################################
## Top 20 Hub Genes per Module
###############################################################

topHubList <- list()

for(mod in modules){

tmp <-

geneInfo %>%

filter(

Module == mod

)

conn <-

kWithin %>%

filter(

Module == mod

)

tmp <-

left_join(

tmp,

conn,

by = "Gene"

)

tmp <-

  tmp %>%

  arrange(

    dplyr::desc(kWithin)

  )

topHubList[[mod]] <-

head(

tmp,

20

)

}

topHubGenes <-

bind_rows(

topHubList,

.id = "Module"

)

write.csv(

topHubGenes,

file.path(

results_wgcna,

"Top20_HubGenes_Per_Module.csv"

),

row.names = FALSE

)

###############################################################
## MM vs GS plots
###############################################################

dir.create(

file.path(

results_wgcna,

"MM_GS_Plots"

),

showWarnings = FALSE

)

###############################################################
## Trait to Visualize in MM vs GS Plots
###############################################################

traitOfInterest <- "VectorNonAttenuated"

traitColumn <-

  match(

    paste0("GS_", traitOfInterest),

    colnames(geneTraitSignificance)

  )

if(is.na(traitColumn))
  stop("Trait not found in geneTraitSignificance.")

for(mod in modules){

column <-

match(

mod,

substring(

names(MEs),

3

)

)

if(is.na(column))

next

moduleGenes <-

moduleColors == mod

pdf(

file.path(

results_wgcna,

"MM_GS_Plots",

paste0(

"MM_vs_GS_",

mod,

".pdf"

)

),

width = 7,

height = 7

)

verboseScatterplot(

abs(

geneModuleMembership[moduleGenes,column]

),

abs(

geneTraitSignificance[
  moduleGenes,
  traitColumn
]

),

xlab="Module Membership",

ylab="Gene Significance",

main=paste(

mod,

"module"

)

)

dev.off()

}


###############################################################
## Session Information
###############################################################

capture.output(

sessionInfo(),

file = file.path(

results_wgcna,

"SessionInfo.txt"

)

)

###############################################################
## Finished
###############################################################

cat("\n")
cat("=============================================\n")
cat(" WGCNA analysis completed successfully.\n")
cat(" Results folder:\n")
cat(results_wgcna,"\n")
cat("=============================================\n")


###############################################################
## TOM for Cytoscape
###############################################################

###############################################################
## Export Each Module
###############################################################

dir.create(

file.path(

results_wgcna,

"Cytoscape"

),

showWarnings = FALSE

)

for(mod in modules){

  probes <- colnames(datExpr)

  inModule <- moduleColors == mod

  modProbes <- probes[inModule]

  if(length(modProbes) < 10)
    next

  if(length(modProbes) > 1000)
    next

  moduleAdjacency <-
    adjacency(
      datExpr[, inModule, drop = FALSE],
      power = softPower,
      type = "signed"
    )

  modTOM <-
    TOMsimilarity(
      moduleAdjacency,
      TOMType = "signed"
    )

dimnames(modTOM) <-

list(

modProbes,

modProbes

)

exportNetworkToCytoscape(

modTOM,

edgeFile=

file.path(

results_wgcna,

"Cytoscape",

paste0(

mod,

"_edges.txt"

)

),

nodeFile=

file.path(

results_wgcna,

"Cytoscape",

paste0(

mod,

"_nodes.txt"

)

),

weighted = TRUE,

threshold = 0.10,

nodeNames = modProbes,

nodeAttr = moduleColors[inModule]

)

}

###############################################################
## Save Final Objects
###############################################################

saveRDS(

geneInfo,

file.path(

results_wgcna,

"GeneInfo.rds"

)

)

saveRDS(

kWithin,

file.path(

results_wgcna,

"Connectivity.rds"

)

)

saveRDS(

hubGenes,

file.path(

results_wgcna,

"HubGenes.rds"

)

)

