###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 06_Vector_Dynamics_and_MixedModels.R
## Purpose :
##   1. Day-level WT vs Δpp71 eigengene trajectories
##   2. Mixed-effects modelling
##   3. Identify modules with significant
##      Vector × Phase interactions
###############################################################

rm(list = ls())

source("scripts/00_setup.R")

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(data.table)

set.seed(12345)

###############################################################
## Directories
###############################################################

results_dir <- file.path(
  RESULTS_DIR,
  "06_Vector_Dynamics",
  DATASET_ID
)

figure_dir <- file.path(
  results_dir,
  "Figures"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
## Load ME Data
###############################################################

ME_data <- readRDS(
  file.path(
    RESULTS_DIR,
    "04_Module_Eigengene_Trajectories",
    DATASET_ID,
    "ME_data.rds"
  )
)

###############################################################
## Factors
###############################################################

phase_order <- c(
  "Prime",
  "Boost1",
  "Memory",
  "Boost2",
  "Durability"
)

ME_data$Phase <- factor(
  ME_data$Phase,
  levels = phase_order
)

ME_data$Sex <- factor(ME_data$Sex)

###############################################################
## Detect modules
###############################################################

modules <- grep(
  "^ME",
  colnames(ME_data),
  value = TRUE
)

modules <- modules[
  modules != "MEgrey"
]

###############################################################
## Theme
###############################################################

theme_pub <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

###############################################################
## PART A
## Phase-level trajectories
###############################################################

phase_summary_all <- list()

for(mod in modules){

  phase_df <- ME_data %>%
    group_by(
      Vector,
      Phase
    ) %>%
    summarise(
      MeanME =
        mean(.data[[mod]], na.rm = TRUE),

      SE =
        sd(.data[[mod]], na.rm = TRUE) /
        sqrt(n()),

      .groups = "drop"
    )

  phase_summary_all[[mod]] <- phase_df

  p <- ggplot(
    phase_df,
    aes(
      x = Phase,
      y = MeanME,
      colour = Vector,
      group = Vector
    )
  ) +
    geom_line(
      linewidth = 1.2
    ) +
    geom_point(
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = MeanME - SE,
        ymax = MeanME + SE
      ),
      width = 0.15
    ) +
    theme_pub +
    labs(
      title = mod,
      x = "Phase",
      y = "Module Eigengene"
    )

  ggsave(
    file.path(
      figure_dir,
      paste0(mod,
             "_PhaseTrajectory.png")
    ),
    p,
    width = 8,
    height = 5,
    dpi = 600
  )
}

###############################################################
## PART B
## Mixed models
###############################################################

model_results <- list()

for (mod in modules) {

  formula_text <- paste0(
    mod,
    " ~ Vector * Phase + Sex + (1|Animal)"
  )

  model <- lmer(
    as.formula(formula_text),
    data = ME_data
  )

  anova_tab <- anova(model)

  interaction_p <- NA_real_
  vector_p <- NA_real_
  phase_p <- NA_real_
  sex_p <- NA_real_

  if ("Vector:Phase" %in% rownames(anova_tab)) {
    interaction_p <- anova_tab["Vector:Phase", ncol(anova_tab)]
  }

  if ("Vector" %in% rownames(anova_tab)) {
    vector_p <- anova_tab["Vector", ncol(anova_tab)]
  }

  if ("Phase" %in% rownames(anova_tab)) {
    phase_p <- anova_tab["Phase", ncol(anova_tab)]
  }

  if ("Sex" %in% rownames(anova_tab)) {
    sex_p <- anova_tab["Sex", ncol(anova_tab)]
  }

  model_results[[mod]] <- data.frame(
    Module = mod,
    Vector_P = vector_p,
    Phase_P = phase_p,
    Sex_P = sex_p,
    Interaction_P = interaction_p,
    stringsAsFactors = FALSE
  )
}

mixed_results <- bind_rows(model_results)

mixed_results$FDR_Vector <- p.adjust(
  mixed_results$Vector_P,
  method = "BH"
)

mixed_results$FDR_Phase <- p.adjust(
  mixed_results$Phase_P,
  method = "BH"
)

mixed_results$FDR_Interaction <- p.adjust(
  mixed_results$Interaction_P,
  method = "BH"
)

write.csv(
  mixed_results,
  file.path(
    results_dir,
    "Module_MixedModel_Results.csv"
  ),
  row.names = FALSE
)

###############################################################
## PART C
## Significant modules
###############################################################

sig_modules <- mixed_results %>%
  filter(
    FDR_Interaction < 0.05
  ) %>%
  arrange(
    FDR_Interaction
  )

write.csv(
  sig_modules,
  file.path(
    results_dir,
    "Significant_VectorPhase_Modules.csv"
  ),
  row.names = FALSE
)

###############################################################
## Save object
###############################################################

saveRDS(
  mixed_results,
  file.path(
    results_dir,
    "Module_MixedModel_Results.rds"
  )
)

cat(
  "\n===================================\n",
  "Vector dynamics analysis complete\n",
  "===================================\n"
)