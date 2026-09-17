library(ggplot2)
library(patchwork)
library(magick)

img1_path <- "results/06_Vector_Dynamics/GSE273911/Figures/MEturquoise_PhaseTrajectory.png"
img2_path <- "results/06_Vector_Dynamics/GSE273911/Figures/MEblue_PhaseTrajectory.png"
img3_path <- "results/06_Vector_Dynamics/GSE273911/Figures/MEbrown_PhaseTrajectory.png"

if (file.exists(img1_path) && file.exists(img2_path) && file.exists(img3_path)) {
  img1 <- image_read(img1_path)
  img2 <- image_read(img2_path)
  img3 <- image_read(img3_path)
  
  # Append horizontally or vertically
  combined <- image_append(c(img1, img2, img3), stack = FALSE)
  
  dir.create("manuscript", showWarnings = FALSE)
  image_write(combined, "manuscript/manuscript_fig4_composite.png")
  cat("Successfully created manuscript/manuscript_fig4_composite.png\n")
} else {
  cat("One or more trajectory images not found.\n")
}
