dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)

files_to_copy <- list(
  "results/WGCNA/Gene_Dendrogram_Modules.png" = "manuscript/figures/Gene_Dendrogram_Modules.png",
  "results/WGCNA/Module_Trait_Heatmap.png"    = "manuscript/figures/Module_Trait_Heatmap.png",
  "results/WGCNA/SoftThreshold.png"           = "manuscript/figures/SoftThreshold.png",
  "results/04_Module_Eigengene_Trajectories/GSE273911/Figures/Module_Activity_Heatmap_Zscore.png" = "manuscript/figures/Module_Activity_Heatmap_Zscore.png",
  "results/04_Module_Eigengene_Trajectories/GSE273911/Figures/Module_Recall_Scores_Barplot.png"   = "manuscript/figures/Module_Recall_Scores_Barplot.png",
  "results/05_Module_Preservation/GSE273911/Figure_Module_Preservation_WT_vs_Delta.png"          = "manuscript/figures/Figure_Module_Preservation_WT_vs_Delta.png",
  "results/06_Vector_Dynamics/GSE273911/Figures/MEturquoise_PhaseTrajectory.png"                  = "manuscript/figures/MEturquoise_PhaseTrajectory.png",
  "results/06_Vector_Dynamics/GSE273911/Figures/MEblue_PhaseTrajectory.png"                       = "manuscript/figures/MEblue_PhaseTrajectory.png",
  "results/06_Vector_Dynamics/GSE273911/Figures/MEbrown_PhaseTrajectory.png"                      = "manuscript/figures/MEbrown_PhaseTrajectory.png",
  "results/07_Module_Enrichment/GSE273911/MEturquoise_GO_BP_Dotplot.png"                         = "manuscript/figures/MEturquoise_GO_BP_Dotplot.png",
  "results/07_Module_Enrichment/GSE273911/MEturquoise_KEGG_Dotplot.png"                          = "manuscript/figures/MEturquoise_KEGG_Dotplot.png",
  "results/07_Module_Enrichment/GSE273911/MEbrown_GO_BP_Dotplot.png"                             = "manuscript/figures/MEbrown_GO_BP_Dotplot.png",
  "results/07_Module_Enrichment/GSE273911/MEbrown_KEGG_Dotplot.png"                             = "manuscript/figures/MEbrown_KEGG_Dotplot.png",
  "manuscript/manuscript_fig4_composite.png"                                                     = "manuscript/figures/manuscript_fig4_composite.png"
)

for (src in names(files_to_copy)) {
  dst <- files_to_copy[[src]]
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat("Copied:", src, "->", dst, "\n")
  } else {
    cat("Warning, not found:", src, "\n")
  }
}
