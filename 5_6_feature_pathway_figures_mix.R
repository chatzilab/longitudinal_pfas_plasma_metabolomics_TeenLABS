code_dir_candidates <- unique(normalizePath(
  c("1_code", ".", ".."),
  mustWork = FALSE
))
code_dir <- code_dir_candidates[
  file.exists(file.path(
    code_dir_candidates,
    "2_longitudinal_metabolomics",
    "5_6_feature_pathway_figures_mix.R"
  ))
][1]
if (is.na(code_dir)) {
  stop("Cannot locate the 1_code directory from the current working directory.")
}
setwd(code_dir)
rm(code_dir, code_dir_candidates)

here::i_am("2_longitudinal_metabolomics/5_6_feature_pathway_figures_mix.R")

################
# Environment Setup
################

options(scipen = 999)

source(fs::path(here::here(), "2_longitudinal_metabolomics", "!libraries.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!directories.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!functions.R"))

library(kableExtra)
library(ggplot2)
library(nlme)
library(ggh4x)
library(readxl)
library(cowplot)


################
# Shorten Long Path
################

create_shortcut <- function(path, short_cut, drive) {
  sub(
    gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", short_cut)),
    paste0(drive, ":"),
    path
  )
}

if (nchar(dir_data) > 80) {
  system("subst x: /D")
  short_cut <- dirname(dir_home)
  system(paste0("subst x: \"", short_cut, "\""))

  dir_objects_names <- ls(pattern = "^dir_")
  dir_list <- mget(dir_objects_names)
  dir_list_shortened <- lapply(
    dir_list,
    function(y) create_shortcut(y, short_cut, "x")
  )
  list2env(dir_list_shortened, envir = .GlobalEnv)

  rm(dir_objects_names, short_cut, dir_list, dir_list_shortened)
}


################
# Model Parameters
################
mix_names <- c(
  "Sum Sulfonic Acids",
  "Sum Carboxylic Acids"
)

mix_slugs <- c(
  "Sum_Sulfonic_Acids",
  "Sum_Carboxylic_Acids"
)

mix_variables <- c(
  "pfas_sulfonic_acids",
  "pfas_carboxylic_acids"
)

covariates <- c(
  "age_baseline",
  "sex",
  "race_binary",
  "parents_income_new",
  "bmi",
  "site"
)

q_cutoff_sig_feat_pca <- 0.000085


################
# Data Load
################
load(
  fs::path(
    dir_data,
    "tl_long_met_untargeted_pfas_workdata.RData"
  )
)

load(
  fs::path(
    dirname(dir_results),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "model_statistics_untargeted_PFAS_mix_longitudinal_metabo_06012026.RData"
  )
)

platforms <- names(lmm_random_intercept_knot1_covars_PFAS_mix_output)

feat_annot_l <- readRDS(
  fs::path(
    dir_home,
    "4_Projects",
    "TL_POPs AT_plasma metabolome_complete (ZLi)",
    "0_data",
    "confirmed_unique_annotation.rds"
  )
)


################
# Significant Features
################
get_sig_feature_names_by_mix <- function(df, var, threshold) {
  sig_rows <- df %>%
    dplyr::mutate(
      effect_term = rep(c("main", "interaction"), length.out = dplyr::n())
    ) %>%
    dplyr::group_by(exposure, effect_term) %>%
    dplyr::mutate(p_adjust = p.adjust(pvalue, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data[[var]] < threshold)

  feature_list_full <- setNames(
    vector("list", length(mix_variables)),
    mix_variables
  )
  for (mix in mix_variables) {
    feature_list_full[[mix]] <- list(
      main = sig_rows$feature[
        sig_rows$exposure == mix & sig_rows$effect_term == "main"
      ],
      interaction = sig_rows$feature[
        sig_rows$exposure == mix & sig_rows$effect_term == "interaction"
      ]
    )
  }

  feature_list_full
}

sig_feature_annot <- purrr::modify2(
  lapply(
    lmm_random_intercept_knot1_covars_PFAS_mix_output,
    function(df) {
      get_sig_feature_names_by_mix(df, "pvalue", q_cutoff_sig_feat_pca)
    }
  ),
  feat_annot_l,
  ~ lapply(.x, function(exposure_list) {
    exposure_list <- exposure_list[c("main", "interaction")]
    lapply(exposure_list, function(feature_vec) {
      get_annotation(feature_vec, .y)
    })
  })
)

sig_feature_annot <- lapply(sig_feature_annot, function(platform_list) {
  platform_list[vapply(
    platform_list,
    function(exposure_list) {
      any(lengths(exposure_list) > 0)
    },
    logical(1)
  )]
})
sig_feature_annot <- sig_feature_annot[lengths(sig_feature_annot) > 0]
sig_feature_annot_final <- sig_feature_annot[c(
  "c18pos",
  "hilicneg",
  "hilicpos"
)]
sig_feature_annot_final <- sig_feature_annot_final[
  lengths(sig_feature_annot_final) > 0
]


################
# Forest Plot
################
model_summary <- list()
for (i in names(sig_feature_annot_final)) {
  model_df <- met_fts_final[[i]] %>%
    dplyr::select(-key, -visit)
  covariates_df <- tl_final[, covariates]
  key <- tl_final$key
  visit <- tl_final$visit_year

  model_output_by_mix <- list()
  for (k in names(sig_feature_annot_final[[i]])) {
    x <- tl_final[, k]
    x_sd <- sd(tl_final[[k]], na.rm = TRUE)

    for (effect_term in c("main", "interaction")) {
      met_vec <- sig_feature_annot_final[[i]][[k]][[effect_term]]
      met_vec <- met_vec[met_vec %in% colnames(model_df)]

      if (length(met_vec) == 0) {
        next
      }

      y_vec <- model_df[, met_vec, drop = FALSE]
      metabolite_lookup <- setNames(names(met_vec), met_vec)

      model_output <- lapply(
        names(y_vec),
        function(fname) {
          y <- y_vec[[fname]]
          fit_lmm_random_intercept_covars_JT_CI(
            x = x,
            y = y,
            key = key,
            visit = visit,
            covariates_df = covariates_df
          ) %>%
            dplyr::mutate(
              metabolite_name = metabolite_lookup[[fname]],
              effect_term = effect_term,
              exposure_sd = x_sd
            )
        }
      ) %>%
        dplyr::bind_rows()

      model_output_by_mix[[paste(k, effect_term, sep = "_")]] <- model_output
      rm(y_vec, model_output)
    }

    rm(x, x_sd)
  }

  model_summary[[i]] <- dplyr::bind_rows(model_output_by_mix)
  rm(model_df, model_output_by_mix)
}

model_summary <- purrr::modify2(
  model_summary,
  names(model_summary),
  ~ .x %>% dplyr::mutate(mode = .y)
) %>%
  dplyr::bind_rows()

model_summary <- model_summary %>%
  dplyr::mutate(
    show_name = dplyr::case_when(
      metabolite_name == "[C26.0]-Hexacosanoic acid" ~ "Hexacosanoic acid",
      metabolite_name == "Quinolinate; Quinolinic acid" ~ "Quinolinic acid",
      metabolite_name == "Spiroxamine.1; Spiroxamine.2" ~ "Spiroxamine",
      TRUE ~ metabolite_name
    ),
    mixture_name = mix_names[match(exposure, mix_variables)],
    percent_change = 100 * (2^(coef * exposure_sd) - 1),
    percent_change_low = 100 * (2^(conf_low * exposure_sd) - 1),
    percent_change_high = 100 * (2^(conf_high * exposure_sd) - 1)
  )

make_forest_plot <- function(effect_term_name) {
  beta_type_name <- ifelse(
    effect_term_name == "main",
    "main_effect",
    "interaction"
  )
  y_axis_label <- ifelse(
    effect_term_name == "main",
    "Percent change per 1-SD increase",
    "Percent change per 1-SD increase in interaction term"
  )

  ggplot(
    model_summary %>%
      dplyr::filter(
        effect_term == effect_term_name,
        beta_type == beta_type_name
      ),
    aes(x = show_name, y = percent_change)
  ) +
    geom_point(shape = 20, size = 2, color = "black") +
    geom_errorbar(
      aes(ymin = percent_change_low, ymax = percent_change_high),
      width = 0.2,
      alpha = 0.7,
      linewidth = 0.7
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylab(y_axis_label) +
    xlab(NULL) +
    facet_grid(
      mixture_name + mode ~ .,
      scales = "free_y",
      space = "free",
      labeller = labeller(
        mode = c(
          "c18pos" = "C18 (+)",
          "hilicneg" = "HILIC (-)",
          "hilicpos" = "HILIC (+)"
        )
      )
    ) +
    coord_flip() +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    theme_classic() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text = element_text(family = "sans"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, color = "black"),
      axis.line.x = element_line(arrow = arrow(length = unit(0.25, "cm"))),
      legend.position = "none",
      panel.background = element_rect(fill = "grey95"),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(size = 9, angle = 0, color = "black")
    )
}

forest_plot_main <- make_forest_plot("main")
forest_plot_interaction <- make_forest_plot("interaction")

forest_plot_combined_mix <- cowplot::plot_grid(
  forest_plot_main,
  forest_plot_interaction,
  ncol = 1,
  rel_heights = c(2, 3),
  align = "v"
)

ggsave(
  filename = fs::path(
    dir_figures,
    "forest_plot_by_untargeted_PFAS_mix_main_interaction.png"
  ),
  plot = forest_plot_combined_mix,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)


################
# Pathway Plot
################
mummichog_folder <- fs::path(
  dirname(dir_results),
  "1_code",
  "2_longitudinal_metabolomics",
  "temp_data",
  "Mummichog_update_mix"
)

mummichog_subfolders <- list.dirs(
  mummichog_folder,
  recursive = FALSE,
  full.names = TRUE
)

read_mummichog_pathway <- function(subfolder) {
  tables_folder <- list.dirs(subfolder, recursive = TRUE, full.names = TRUE)
  tables_path <- tables_folder[grepl("tables$", tables_folder)]

  if (length(tables_path) == 0) {
    return(NULL)
  }

  xlsx_path <- fs::path(tables_path[1], "mcg_pathwayanalysis_.xlsx")
  tsv_path <- fs::path(tables_path[1], "mcg_pathwayanalysis_.tsv")

  if (file.exists(xlsx_path)) {
    readxl::read_xlsx(xlsx_path)
  } else if (file.exists(tsv_path)) {
    readr::read_tsv(tsv_path, show_col_types = FALSE)
  } else {
    NULL
  }
}

pathway_summary <- lapply(mummichog_subfolders, read_mummichog_pathway)
names(pathway_summary) <- basename(mummichog_subfolders)

pathway_summary <- purrr::imap(
  pathway_summary,
  ~ {
    if (is.null(.x) || nrow(.x) == 0) {
      return(NULL)
    }

    .x %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::any_of(c(
            "overlap_EmpiricalCompounds (id)",
            "overlap_features (id)",
            "overlap_features (name)"
          )),
          as.character
        )
      ) %>%
      dplyr::filter(.data[["p-value"]] < 0.05, overlap_size >= 2) %>%
      dplyr::mutate(run_name = .y)
  }
)
pathway_summary <- pathway_summary[
  !vapply(pathway_summary, is.null, logical(1))
]
pathway_summary <- dplyr::bind_rows(pathway_summary)

pathway_summary <- pathway_summary %>%
  dplyr::mutate(
    effect_term = sub(".*_([^_]+)$", "\\1", run_name),
    mode = sub("^([^_]+)_.*$", "\\1", run_name),
    mix_slug = sub("^[^_]+_(.*)_[^_]+$", "\\1", run_name),
    mixture_name = mix_names[match(mix_slug, mix_slugs)],
    enrichment = as.numeric(overlap_size / pathway_size),
    superclass = dplyr::case_when(
      pathway %in%
        c(
          "Glycosphingolipid metabolism",
          "Butanoate metabolism"
        ) ~ "Others",
      pathway %in%
        c(
          "Arachidonic acid metabolism",
          "Prostaglandin formation from arachidonate",
          "Bile acid biosynthesis",
          "Ascorbate (Vitamin C) and Aldarate Metabolism",
          "C21-steroid hormone biosynthesis and metabolism"
        ) ~ "Lipid metabolism",
      pathway %in%
        c(
          "Histidine metabolism",
          "Aspartate and asparagine metabolism",
          "Tryptophan metabolism"
        ) ~ "Amino acid metabolism",
      TRUE ~ "Others"
    )
  ) %>%
  dplyr::filter(mode %in% c("c18pos", "hilicneg", "hilicpos"))


################
# Pathway Plot Output
################
make_pathway_plot_combined_part <- function(effect_term_name) {
  ggplot(
    data = pathway_summary %>% dplyr::filter(effect_term == effect_term_name),
    aes(
      x = -log10(.data[["p-value"]]),
      y = pathway,
      size = enrichment,
      shape = mode
    )
  ) +
    geom_point(aes(fill = superclass, color = superclass), show.legend = TRUE) +
    scale_size_continuous(
      range = c(2.5, 6),
      name = "Enrichment",
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_fill_manual(
      values = c(
        "Amino acid metabolism" = "#1f77b4",
        "Lipid metabolism" = "#d95f02",
        "Others" = "grey45"
      ),
      name = "Pathway group"
    ) +
    scale_color_manual(
      values = c(
        "Amino acid metabolism" = "#1f77b4",
        "Lipid metabolism" = "#d95f02",
        "Others" = "grey45"
      ),
      name = "Pathway group"
    ) +
    scale_shape_manual(
      values = c(c18pos = 21, hilicneg = 22, hilicpos = 24),
      name = "Mode",
      labels = c(
        c18pos = "C18 (+)",
        hilicneg = "HILIC (-)",
        hilicpos = "HILIC (+)"
      )
    ) +
    ggh4x::facet_nested(
      mixture_name ~ .,
      scales = "free_y",
      space = "free",
      strip = ggh4x::strip_split("left")
    ) +
    theme_classic() +
    xlab(expression(-log[10](`p-value`))) +
    ylab(NULL) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme(
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, color = "black"),
      axis.ticks = element_blank(),
      axis.line.x = element_line(arrow = arrow(length = unit(0.25, "cm"))),
      legend.title = element_text(face = "plain"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.box.spacing = unit(0, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey95"),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text.y.right = element_text(size = 9, angle = 0, hjust = 0),
      strip.text.y.left = element_text(size = 10)
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(size = 6)),
      fill = guide_legend(order = 2, override.aes = list(shape = 21, size = 6)),
      size = guide_legend(order = 3),
      color = "none"
    )
}

pathway_plot_combined_interaction <- make_pathway_plot_combined_part(
  "interaction"
)
pathway_plot_combined_legend <- cowplot::get_legend(
  pathway_plot_combined_interaction
)

pathway_plot_combined_mix <- cowplot::plot_grid(
  pathway_plot_combined_interaction + theme(legend.position = "none"),
  pathway_plot_combined_legend,
  ncol = 1,
  rel_heights = c(1, 0.34)
)

ggsave(
  filename = fs::path(
    dir_figures,
    "pathway_analysis_by_untargeted_PFAS_mix_interaction.png"
  ),
  plot = pathway_plot_combined_mix,
  width = 6,
  height = 2.5,
  units = "in",
  dpi = 600
)


################
# Figure Statistics Output
################
write.csv(
  model_summary,
  fs::path(dir_results, "Figure_1_statistics_mix_main_interaction.csv"),
  row.names = FALSE
)

write.csv(
  pathway_summary,
  fs::path(dir_results, "Figure_2_statistics_mix_main_interaction.csv"),
  row.names = FALSE
)
