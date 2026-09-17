library(pdftools)
pdf_files <- c(
  "results/WGCNA/Gene_Dendrogram_Modules.pdf",
  "results/WGCNA/Module_Trait_Heatmap.pdf",
  "results/WGCNA/SoftThreshold.pdf",
  "results/07_Module_Enrichment/GSE273911/MEturquoise_GO_BP_Dotplot.pdf",
  "results/07_Module_Enrichment/GSE273911/MEturquoise_KEGG_Dotplot.pdf",
  "results/07_Module_Enrichment/GSE273911/MEbrown_GO_BP_Dotplot.pdf",
  "results/07_Module_Enrichment/GSE273911/MEbrown_KEGG_Dotplot.pdf"
)
for (f in pdf_files) {
  if (file.exists(f)) {
    png_out <- sub("\\.pdf$", ".png", f)
    pdftools::pdf_convert(f, format = "png", pages = 1, dpi = 300, filenames = png_out)
    cat("Converted:", f, "->", png_out, "\n")
  } else {
    cat("File not found:", f, "\n")
  }
}
