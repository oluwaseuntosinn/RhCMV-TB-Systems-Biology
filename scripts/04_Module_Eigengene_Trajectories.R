###############################################################
## Project : RhCMV-TB Systems Vaccinology
## Dataset : GSE273911
## Script  : 04_Module_Eigengene_Trajectories.R
## Purpose : Module Eigengene Longitudinal Trajectories, Kinetics,
##           Vector Statistical Divergence, and Recall Scoring
###############################################################

rm(list = ls())

###############################################################
## Load Project Setup
###############################################################

source("scripts/00_setup.R")

library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)

set.seed(12345)

###############################################################
## Directories
###############################################################

GEO <- DATASET_ID

wgcna_dir <- WGCNA_DIR

results_dir <- file.path(
  RESULTS_DIR,
  "04_Module_Eigengene_Trajectories",
  GEO
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
## Load Eigengenes
###############################################################

if (file.exists(file.path(wgcna_dir, "Module_Eigengenes.rds"))) {
  MEs <- readRDS(file.path(wgcna_dir, "Module_Eigengenes.rds"))
  MEs$SampleID <- rownames(MEs)
} else if (file.exists(file.path(wgcna_dir, "Module_Eigengenes.csv"))) {
  MEs <- read.csv(file.path(wgcna_dir, "Module_Eigengenes.csv"), check.names = FALSE)
  if (colnames(MEs)[1] %in% c("X", "")) {
    colnames(MEs)[1] <- "SampleID"
  }
} else {
  stop("Module_Eigengenes file not found in WGCNA results directory.")
}

###############################################################
## Load Metadata
###############################################################

metadata_file <- file.path(PROCESSED_DATA, "metadata_clean.csv")

if (!file.exists(metadata_file)) {
  stop("Metadata clean file not found at: ", metadata_file)
}

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
) %>%
  as.data.frame()

###############################################################
## Merge Metadata & Eigengenes
###############################################################

ME_data <- metadata %>%
  inner_join(MEs, by = "SampleID")

cat("\nSamples merged:", nrow(ME_data), "\n")

###############################################################
## Ensure Ordered Phases and Day Factors
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

unique_days <- unique(ME_data %>% arrange(DayNumeric) %>% pull(Day))
ME_data$Day <- factor(
  ME_data$Day,
  levels = unique_days
)

###############################################################
## Detect Modules Automatically
###############################################################
modules <- grep(
  "^ME",
  colnames(ME_data),
  value = TRUE
)

modules <- modules[
  modules != "MEgrey"
]

cat("Detected modules (", length(modules), "):\n", sep = "")
cat(paste(modules, collapse = ", "), "\n")

###############################################################
## Publication Theme
###############################################################

theme_pub <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

###############################################################
## Part A — Mean Trajectories (by Phase)
###############################################################

cat("\n=============================================\n")
cat(" Generating Part A: Mean Phase Trajectory Plots\n")
cat("=============================================\n")

for (module in modules) {

  cat("Processing mean phase trajectory:", module, "\n")

  plot_df <- ME_data %>%
    group_by(
      Vector,
      Phase
    ) %>%
    summarise(
      MeanME = mean(.data[[module]], na.rm = TRUE),
      SE = sd(.data[[module]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  p <- ggplot(
    plot_df,
    aes(
      x = Phase,
      y = MeanME,
      group = Vector,
      colour = Vector
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
      title = module,
      y = "Module Eigengene",
      x = NULL
    )

  ggsave(
    file.path(
      figure_dir,
      paste0(module, "_Trajectory.png")
    ),
    p,
    width = 7,
    height = 5,
    dpi = 600
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(module, "_Trajectory.pdf")
    ),
    p,
    width = 7,
    height = 5
  )
}

###############################################################
## Part B — Individual Animal Trajectories (Day-Level)
###############################################################

cat("\n=============================================\n")
cat(" Generating Part B: Individual Animal Trajectories\n")
cat("=============================================\n")

for (module in modules) {

  cat("Processing individual trajectories:", module, "\n")

  p <- ggplot(
    ME_data,
    aes(
      x = DayNumeric,
      y = .data[[module]],
      group = Animal,
      colour = Vector
    )
  ) +
    geom_line(
      alpha = 0.4
    ) +
    geom_smooth(
      aes(group = Vector),
      se = FALSE,
      linewidth = 1.5
    ) +
    theme_pub +
    labs(
      title = paste(
        module,
        "Individual Trajectories"
      ),
      x = "Day",
      y = "Module Eigengene"
    )

  ggsave(
    file.path(
      figure_dir,
      paste0(module, "_Individual_Trajectories.png")
    ),
    p,
    width = 7,
    height = 5,
    dpi = 600
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(module, "_Individual_Trajectories.pdf")
    ),
    p,
    width = 7,
    height = 5
  )
}

###############################################################
## Part C — Export Phase Means & Day Means (Improvement 2)
###############################################################

cat("\n=============================================\n")
cat(" Generating Part C: Export Phase & Day Means\n")
cat("=============================================\n")

## 1. Phase Means
phase_summary <- ME_data %>%
  pivot_longer(
    cols = all_of(modules),
    names_to = "Module",
    values_to = "Eigengene"
  ) %>%
  group_by(
    Module,
    Vector,
    Phase
  ) %>%
  summarise(
    MeanME = mean(Eigengene, na.rm = TRUE),
    SD = sd(Eigengene, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

write.csv(
  phase_summary,
  file.path(
    results_dir,
    "Module_Phase_Means.csv"
  ),
  row.names = FALSE
)

saveRDS(
  phase_summary,
  file.path(
    results_dir,
    "Module_Phase_Means.rds"
  )
)

## 2. Day Means (Captures fine-scale kinetics)
day_summary <- ME_data %>%
  pivot_longer(
    cols = all_of(modules),
    names_to = "Module",
    values_to = "Eigengene"
  ) %>%
  group_by(
    Module,
    Vector,
    Phase,
    Day,
    DayNumeric
  ) %>%
  summarise(
    MeanME = mean(Eigengene, na.rm = TRUE),
    SD = sd(Eigengene, na.rm = TRUE),
    SE = sd(Eigengene, na.rm = TRUE) / sqrt(n()),
    N = n(),
    .groups = "drop"
  ) %>%
  arrange(Module, Vector, DayNumeric)

write.csv(
  day_summary,
  file.path(
    results_dir,
    "Module_Day_Means.csv"
  ),
  row.names = FALSE
)

saveRDS(
  day_summary,
  file.path(
    results_dir,
    "Module_Day_Means.rds"
  )
)

saveRDS(
  ME_data,
  file.path(
    results_dir,
    "ME_data.rds"
  )
)

###############################################################
## Part D — Standardized Heatmap of Module Activity (Improvement 1)
###############################################################

cat("\n=============================================\n")
cat(" Generating Part D: Standardized Module Activity Heatmap\n")
cat("=============================================\n")

heat_df <- phase_summary %>%
  group_by(Module, Vector) %>%
  mutate(
    Zscore = as.vector(scale(MeanME))
  ) %>%
  ungroup()

p_heat_zscore <- ggplot(
  heat_df,
  aes(
    x = Phase,
    y = Module,
    fill = Zscore
  )
) +
  geom_tile(
    colour = "white"
  ) +
  facet_wrap(
    ~Vector
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Z-score (Mean ME)"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Standardized Module Activity Across Vaccination Phases",
    x = "Phase",
    y = "Module"
  )

ggsave(
  file.path(
    figure_dir,
    "Module_Activity_Heatmap_Zscore.png"
  ),
  p_heat_zscore,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    figure_dir,
    "Module_Activity_Heatmap_Zscore.pdf"
  ),
  p_heat_zscore,
  width = 10,
  height = 6
)

###############################################################
## Part E — Vector Difference Statistics per Phase (Improvement 3)
###############################################################

cat("\n=============================================\n")
cat(" Generating Part E: Vector Statistical Testing per Phase\n")
cat("=============================================\n")

vector_stats_list <- list()

for (mod in modules) {
  for (ph in phase_order) {
    sub_df <- ME_data %>%
      filter(Phase == ph) %>%
      select(Animal, Vector, Eigengene = all_of(mod)) %>%
      filter(!is.na(Eigengene))

    n_vec <- length(unique(sub_df$Vector))
    
    if (n_vec >= 2) {
      t_res <- tryCatch(
        t.test(Eigengene ~ Vector, data = sub_df),
        error = function(e) NULL
      )

      w_res <- tryCatch(
        wilcox.test(Eigengene ~ Vector, data = sub_df),
        error = function(e) NULL
      )

      means <- sub_df %>%
        group_by(Vector) %>%
        summarise(m = mean(Eigengene, na.rm = TRUE), .groups = "drop")

      mean_diff <- if (nrow(means) == 2) means$m[1] - means$m[2] else NA_real_

      vector_stats_list[[length(vector_stats_list) + 1]] <- data.frame(
        Module = mod,
        Phase = ph,
        Mean_Vector1 = ifelse(nrow(means) >= 1, means$m[1], NA_real_),
        Vector1_Name = ifelse(nrow(means) >= 1, as.character(means$Vector[1]), NA_character_),
        Mean_Vector2 = ifelse(nrow(means) >= 2, means$m[2], NA_real_),
        Vector2_Name = ifelse(nrow(means) >= 2, as.character(means$Vector[2]), NA_character_),
        MeanDifference = mean_diff,
        Pvalue_TTest = ifelse(!is.null(t_res), t_res$p.value, NA_real_),
        Pvalue_Wilcox = ifelse(!is.null(w_res), w_res$p.value, NA_real_),
        stringsAsFactors = FALSE
      )
    }
  }
}

vector_diff_stats <- bind_rows(vector_stats_list) %>%
  group_by(Phase) %>%
  mutate(
    FDR_TTest = p.adjust(Pvalue_TTest, method = "BH"),
    FDR_Wilcox = p.adjust(Pvalue_Wilcox, method = "BH")
  ) %>%
  ungroup()

write.csv(
  vector_diff_stats,
  file.path(
    results_dir,
    "Vector_Difference_Statistics_Per_Phase.csv"
  ),
  row.names = FALSE
)

saveRDS(
  vector_diff_stats,
  file.path(
    results_dir,
    "Vector_Difference_Statistics_Per_Phase.rds"
  )
)

###############################################################
## Part F — Module Recall Scoring Analysis (Improvement 4)
###############################################################

cat("\n=============================================\n")
cat(" Generating Part F: Module Recall Scoring\n")
cat("=============================================\n")

## Calculate Baseline (D0) per Animal and baseline-adjusted responses
baseline_df <- ME_data %>%
  filter(DayNumeric == 0) %>%
  select(Animal, all_of(modules)) %>%
  pivot_longer(
    cols = all_of(modules),
    names_to = "Module",
    values_to = "BaselineME"
  )

long_me <- ME_data %>%
  select(SampleID, Animal, Vector, Phase, Day, DayNumeric, all_of(modules)) %>%
  pivot_longer(
    cols = all_of(modules),
    names_to = "Module",
    values_to = "Eigengene"
  ) %>%
  left_join(baseline_df, by = c("Animal", "Module")) %>%
  mutate(
    Response_vs_D0 = Eigengene - BaselineME
  )

## Peak response during vaccination phases (Prime, Boost1, Boost2)
peak_responses <- long_me %>%
  filter(Phase %in% c("Prime", "Boost1", "Boost2")) %>%
  group_by(Animal, Vector, Module, Phase) %>%
  summarise(
    PeakResponse =
      Response_vs_D0[which.max(abs(Response_vs_D0))],

    PeakMagnitude =
      max(abs(Response_vs_D0), na.rm = TRUE),

    PeakResponse_Max =
      max(Response_vs_D0, na.rm = TRUE),

    MeanResponse =
      mean(Response_vs_D0, na.rm = TRUE),

    .groups = "drop"
  )

## Animal-level Recall Score
animal_recall <- peak_responses %>%
  pivot_wider(
    id_cols = c(Animal, Vector, Module),
    names_from = Phase,
    values_from = PeakResponse,
    names_prefix = "Peak_"
  ) %>%
  mutate(
    RecallScore =
      (
        abs(Peak_Prime) +
        abs(Peak_Boost1) +
        abs(Peak_Boost2)
      ) / 3
  )


## Vector-level Module Recall Summary
module_recall_summary <- animal_recall %>%
  group_by(Module, Vector) %>%
  summarise(
    Mean_Peak_Prime =
      mean(Peak_Prime, na.rm = TRUE),

    Mean_Peak_Boost1 =
      mean(Peak_Boost1, na.rm = TRUE),

    Mean_Peak_Boost2 =
      mean(Peak_Boost2, na.rm = TRUE),

    Mean_RecallScore =
      mean(RecallScore, na.rm = TRUE),

    SE_RecallScore =
      sd(RecallScore, na.rm = TRUE) /
      sqrt(n()),

    .groups = "drop"
  ) %>%
  mutate(

    RecallClass = case_when(

      Mean_Peak_Prime > 0 &
      Mean_Peak_Boost1 > 0 &
      Mean_Peak_Boost2 > 0
      ~ "RepeatedRecall",

      Mean_Peak_Prime > 0 &
      Mean_Peak_Boost1 > 0
      ~ "PartialRecall",

      TRUE
      ~ "WeakRecall"

    )

  ) %>%
  arrange(desc(Mean_RecallScore))

write.csv(
  animal_recall,
  file.path(
    results_dir,
    "Animal_Level_Recall_Scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_recall_summary,
  file.path(
    results_dir,
    "Module_Recall_Scores_Summary.csv"
  ),
  row.names = FALSE
)


###############################################################
## Top Recall Modules
###############################################################

top_recall_modules <- module_recall_summary %>%
  group_by(Vector) %>%
  slice_max(
    Mean_RecallScore,
    n = 5
  ) %>%
  ungroup()

write.csv(
  top_recall_modules,
  file.path(
    results_dir,
    "Top_Recall_Modules.csv"
  ),
  row.names = FALSE
)

saveRDS(
  module_recall_summary,
  file.path(
    results_dir,
    "Module_Recall_Scores_Summary.rds"
  )
)

## Barplot of Module Recall Scores by Vector
p_recall <- ggplot(
  module_recall_summary,
  aes(
    x = reorder(Module, -Mean_RecallScore),
    y = Mean_RecallScore,
    fill = Vector
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_RecallScore - SE_RecallScore,
      ymax = Mean_RecallScore + SE_RecallScore
    ),
    position = position_dodge(0.8),
    width = 0.25
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Module Recall Scores Across RhCMV-TB Vectors",
    x = "Co-expression Module",
    y = "Recall Score (Mean Peak Response)"
  )

ggsave(
  file.path(
    figure_dir,
    "Module_Recall_Scores_Barplot.png"
  ),
  p_recall,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    figure_dir,
    "Module_Recall_Scores_Barplot.pdf"
  ),
  p_recall,
  width = 10,
  height = 6
)

###############################################################
## Finished
###############################################################

cat("\n=============================================\n")
cat(" Module Eigengene Trajectory analysis completed.\n")
cat(" Results folder:\n")
cat(results_dir, "\n")
cat("=============================================\n")
