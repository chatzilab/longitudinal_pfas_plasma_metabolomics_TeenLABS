setwd("1_code")
here::i_am("2_longitudinal_metabolomics/5_1_MWAS_LMM.R")

source(fs::path(here::here(), "2_longitudinal_metabolomics", "!libraries.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!directories.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!functions.R"))

library(kableExtra)
library(ggplot2)
library(nlme)
library(ggh4x)
library(readxl)
library(cowplot)

create_shortcut <- function(path, short_cut, drive) {
  return(
    sub(
      gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", short_cut)),
      paste0(drive, ":"),
      path
    )
  )
}

if (nchar(dir_data) > 80) {
  system("subst x: /D")
  short_cut <- dir_home |> dirname()
  system(paste0("subst x: \"", short_cut, "\""))

  dir_objects_names <- ls(pattern = "^dir_")

  dir_list <- mget(dir_objects_names)

  dir_list_shortened <- lapply(
    dir_list,
    function(y) create_shortcut(y, short_cut, "x")
  )
  list2env(dir_list_shortened, envir = .GlobalEnv)
}

rm(dir_objects_names, short_cut, dir_list, dir_list_shortened)

####################
# Model Parameters
####################

chemical_names <- c(
  "PFOS",
  "PFHxS",
  "PFHpS",
  "PFOA",
  "PFNA",
  "PFDA",
  "PFHpA",
  "PFUnDa"
)

chemical_variables <- c(
  "pfos_untargeted_plasma",
  "pf_hx_s_untargeted_plasma",
  "pf_hp_s_untargeted_plasma",
  "pfoa_untargeted_plasma",
  "pfna_untargeted_plasma",
  "pfda_untargeted_plasma",
  "pf_hp_a_untargeted_plasma",
  "pf_ud_a_untargeted_plasma"
)

covariates <- c(
  "age_baseline",
  "sex",
  "race_binary",
  "parents_income_new",
  "bmi",
  "site"
)

q_cutoff_sig_feat <- 0.05
q_cutoff_sig_feat_pca <- 0.000085

#####################
# Data Load
#####################

load(
  fs::path(
    dir_data,
    "tl_long_met_untargeted_pfas_workdata.RData"
  )
)

platforms <- names(met_fts_final)

load(
  fs::path(
    dir_results |> dirname(),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "model_statistics_untargeted_PFAS_longitudinal_metabo_04302026.RData"
  )
)

platforms <- names(lmm_random_intercept_knot1_covars_PFAS_output)

feat_annot_l <- readRDS(
  fs::path(
    dir_home,
    "4_Projects",
    "TL_POPs AT_plasma metabolome_complete (ZLi)",
    "0_data",
    "confirmed_unique_annotation.rds"
  )
)

#####################
# Summary
#####################

###############################
# Annotation of Sig. features
###############################

get_sig_feature_names_by_chemical <- function(df, var, threshold) {
  sig_rows <- df %>%
    dplyr::mutate(
      effect_term = rep(c("main", "interaction"), length.out = dplyr::n())
    ) %>%
    dplyr::group_by(exposure, effect_term) %>%
    dplyr::mutate(
      p_adjust = p.adjust(pvalue, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data[[var]] < threshold)

  feature_list_full <- setNames(
    vector("list", length(chemical_variables)),
    chemical_variables
  )
  for (chem in chemical_variables) {
    feature_list_full[[chem]] <- list(
      main = sig_rows$feature[
        sig_rows$exposure == chem & sig_rows$effect_term == "main"
      ],
      interaction = sig_rows$feature[
        sig_rows$exposure == chem & sig_rows$effect_term == "interaction"
      ]
    )
  }

  feature_list_full
}

sig_feature_annot <- purrr::modify2(
  lapply(
    lmm_random_intercept_knot1_covars_PFAS_output,
    function(df) {
      get_sig_feature_names_by_chemical(df, "pvalue", q_cutoff_sig_feat_pca)
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

# After manually removing the external exposures from metabolites, keep HILIC
# columns and C18 positive.
sig_feature_annot_final <- sig_feature_annot[c(
  "c18pos",
  "hilicneg",
  "hilicpos"
)]
sig_feature_annot_final <- sig_feature_annot_final[
  lengths(sig_feature_annot_final) > 0
]

####################################
# Model statistics of sig. features
####################################

model_summary <- list()
for (i in names(sig_feature_annot_final)) {
  model_df <- met_fts_final[[i]] %>%
    dplyr::select(-key, -visit)
  covariates_df <- tl_final[, covariates]
  key <- tl_final$key
  visit <- tl_final$visit_year

  model_output_by_chemical <- list()
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

      model_output_by_chemical[[paste(
        k,
        effect_term,
        sep = "_"
      )]] <- model_output

      rm(y_vec, model_output)
    }

    rm(x)
    rm(x_sd)
  }
  model_summary[[i]] <- dplyr::bind_rows(model_output_by_chemical)

  rm(model_df, model_output_by_chemical)
}

model_summary <- purrr::modify2(
  model_summary,
  names(model_summary),
  ~ .x %>%
    dplyr::mutate(
      mode = .y
    )
) %>%
  dplyr::bind_rows()

# Remove arachidate because it has opposite results across different modes.
model_summary <- model_summary %>%
  dplyr::filter(metabolite_name != "Arachidate")

#####################
# Forest plots
#####################

model_summary <- model_summary %>%
  dplyr::mutate(
    show_name = dplyr::case_when(
      metabolite_name == "[C26.0]-Hexacosanoic acid" ~ "Hexacosanoic acid",
      metabolite_name == "Quinolinate; Quinolinic acid" ~ "Quinolinic acid",
      TRUE ~ metabolite_name
    ),
    pfas_name = chemical_names[match(exposure, chemical_variables)],
    effect_label = dplyr::case_when(
      effect_term == "main" ~ "Main effect selected features",
      effect_term == "interaction" ~ "Interaction selected features",
      TRUE ~ effect_term
    ),
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
    "Percent change in main effect term",
    "Percent change in interaction term"
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
      pfas_name + mode ~ .,
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

# ggsave(
#   filename = fs::path(dir_figures, "forest_plot_by_untargeted_PFAS_main.png"),
#   plot = forest_plot_main,
#   width = 6,
#   height = 5,
#   units = "in",
#   dpi = 600
# )

# ggsave(
#   filename = fs::path(
#     dir_figures,
#     "forest_plot_by_untargeted_PFAS_interaction.png"
#   ),
#   plot = forest_plot_interaction,
#   width = 6,
#   height = 5,
#   units = "in",
#   dpi = 600
# )

forest_plot_combined <- cowplot::plot_grid(
  forest_plot_main,
  forest_plot_interaction,
  ncol = 1,
  rel_heights = c(2, 3),
  align = "v"
)

ggsave(
  filename = fs::path(
    dir_figures,
    "forest_plot_by_untargeted_PFAS_main_interaction.png"
  ),
  plot = forest_plot_combined,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

#####################
# Timepoint plots
#####################

visit_timepoints <- tibble::tibble(
  visit_year = c(0, 0.5, 1, 3),
  timepoint = factor(
    c("Baseline", "0.5y", "1y", "3y"),
    levels = c("Baseline", "0.5y", "1y", "3y")
  )
)

fit_lmm_timepoint_effects <- function(
  x,
  y,
  feature_name,
  key,
  visit,
  covariates_df,
  timepoints
) {
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)

  model <- as.formula(paste(
    "y ~ x*visit",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))

  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  beta <- fixed.effects(fit)
  beta_vcov <- vcov(fit)
  x_name <- grep("^x$", names(beta), value = TRUE)
  interaction_name <- grep("(^x:)|(:x$)", names(beta), value = TRUE)[1]

  dplyr::bind_rows(lapply(seq_len(nrow(timepoints)), function(i) {
    time_value <- timepoints$visit_year[i]
    contrast <- rep(0, length(beta))
    names(contrast) <- names(beta)
    contrast[x_name] <- 1
    contrast[interaction_name] <- time_value

    coef <- sum(contrast * beta)
    se <- sqrt(as.numeric(t(contrast) %*% beta_vcov %*% contrast))
    conf_low <- coef - qnorm(0.975) * se
    conf_high <- coef + qnorm(0.975) * se

    tibble::tibble(
      coef = coef,
      se = se,
      conf_low = conf_low,
      conf_high = conf_high,
      visit_year = time_value,
      timepoint = timepoints$timepoint[i],
      exposure = names(x),
      feature = feature_name
    )
  }))
}

sig_feature_lookup <- purrr::imap_dfr(
  sig_feature_annot_final,
  function(mode_list, mode_name) {
    purrr::imap_dfr(mode_list, function(exposure_list, exposure_name) {
      feature_vec <- unlist(exposure_list, use.names = TRUE)
      if (length(feature_vec) == 0) {
        return(NULL)
      }
      feature_names <- names(feature_vec)
      if (is.null(feature_names)) {
        feature_names <- rep(NA_character_, length(feature_vec))
      }
      feature_names <- sub("^(main|interaction)\\.", "", feature_names)

      tibble::tibble(
        mode = mode_name,
        exposure = exposure_name,
        feature = unname(feature_vec),
        metabolite_name = feature_names
      ) %>%
        dplyr::distinct(feature, .keep_all = TRUE)
    })
  }
) %>%
  dplyr::distinct(mode, exposure, feature, metabolite_name)

feature_timepoint_summary <- list()
for (i in names(sig_feature_annot_final)) {
  model_df <- met_fts_final[[i]] %>%
    dplyr::select(-key, -visit)
  covariates_df <- tl_final[, covariates]
  key <- tl_final$key
  visit <- tl_final$visit_year

  selected_features <- sig_feature_lookup %>%
    dplyr::filter(mode == i, feature %in% colnames(model_df))

  model_output_by_chemical <- list()
  for (k in unique(selected_features$exposure)) {
    x <- tl_final[, k]
    x_sd <- sd(tl_final[[k]], na.rm = TRUE)
    met_lookup <- selected_features %>%
      dplyr::filter(exposure == k) %>%
      dplyr::select(feature, metabolite_name)

    model_output <- lapply(
      met_lookup$feature,
      function(fname) {
        y <- model_df[[fname]]
        names(y) <- fname

        fit_lmm_timepoint_effects(
          x = x,
          y = y,
          feature_name = fname,
          key = key,
          visit = visit,
          covariates_df = covariates_df,
          timepoints = visit_timepoints
        ) %>%
          dplyr::mutate(
            metabolite_name = met_lookup$metabolite_name[
              match(fname, met_lookup$feature)
            ],
            exposure_sd = x_sd
          )
      }
    ) %>%
      dplyr::bind_rows()

    model_output_by_chemical[[k]] <- model_output

    rm(x)
    rm(x_sd)
    rm(model_output)
  }

  feature_timepoint_summary[[i]] <- dplyr::bind_rows(model_output_by_chemical)

  rm(model_df, model_output_by_chemical)
}

feature_timepoint_summary <- purrr::modify2(
  feature_timepoint_summary,
  names(feature_timepoint_summary),
  ~ .x %>%
    dplyr::mutate(
      mode = .y
    )
) %>%
  dplyr::bind_rows()

feature_timepoint_summary <- feature_timepoint_summary %>%
  dplyr::filter(metabolite_name != "Arachidate") %>%
  dplyr::mutate(
    show_name = dplyr::case_when(
      metabolite_name == "[C26.0]-Hexacosanoic acid" ~ "Hexacosanoic acid",
      metabolite_name == "Quinolinate; Quinolinic acid" ~ "Quinolinic acid",
      TRUE ~ metabolite_name
    ),
    pfas_name = chemical_names[match(exposure, chemical_variables)],
    percent_change = 100 * (2^(coef * exposure_sd) - 1),
    percent_change_low = 100 * (2^(conf_low * exposure_sd) - 1),
    percent_change_high = 100 * (2^(conf_high * exposure_sd) - 1)
  )

feature_timepoint_plot <- ggplot(
  feature_timepoint_summary,
  aes(x = timepoint, y = percent_change, group = 1)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_line(linewidth = 0.5, color = "black") +
  geom_pointrange(
    aes(ymin = percent_change_low, ymax = percent_change_high),
    linewidth = 0.45,
    color = "black"
  ) +
  xlab("Timepoint after surgery") +
  ylab("Percent change") +
  facet_grid(
    pfas_name + mode + show_name ~ .,
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
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(family = "sans"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.line.y = element_line(arrow = arrow(length = unit(0.25, "cm"))),
    legend.position = "none",
    panel.background = element_rect(fill = "grey95"),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, angle = 0, color = "black")
  )

feature_timepoint_plot_height <- max(
  6,
  0.45 *
    dplyr::n_distinct(
      paste(
        feature_timepoint_summary$pfas_name,
        feature_timepoint_summary$mode,
        feature_timepoint_summary$show_name
      )
    )
)

ggsave(
  filename = fs::path(
    dir_figures,
    "feature_timepoint_percent_change_by_untargeted_PFAS.png"
  ),
  plot = feature_timepoint_plot,
  width = 6,
  height = feature_timepoint_plot_height + 2,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)

#####################
# Pathway plots
#####################·

mummichog_folder <- fs::path(
  dir_results |> dirname(),
  "1_code",
  "2_longitudinal_metabolomics",
  "temp_data",
  "Mummichog_update"
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
    pfas_name = sub("^[^_]+_(.*)_[^_]+$", "\\1", run_name),
    pfas_binary = dplyr::case_when(
      pfas_name %in% c("PFOS", "PFHxS", "PFHpS") ~ "Sulfonic acids",
      TRUE ~ "Carboxylic acids"
    ),
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
  dplyr::filter(mode %in% c("c18pos", "hilicneg", "hilicpos")) %>%
  dplyr::mutate(
    effect_label = factor(
      dplyr::case_when(
        effect_term == "main" ~ "Main",
        effect_term == "interaction" ~ "Interaction",
        TRUE ~ effect_term
      ),
      levels = c("Main", "Interaction")
    )
  )

make_pathway_plot <- function(effect_term_name) {
  ggplot(
    data = pathway_summary %>%
      dplyr::filter(effect_term == effect_term_name),
    aes(
      x = -log10(.data[["p-value"]]),
      y = pathway,
      size = enrichment,
      shape = mode
    )
  ) +
    geom_point(aes(fill = mode, color = mode), show.legend = TRUE) +
    scale_size_continuous(
      range = c(4, 9),
      name = "Enrichment",
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_shape_manual(
      values = c("c18pos" = 21, "hilicneg" = 22, "hilicpos" = 24),
      name = "Mode",
      labels = c(
        "c18pos" = "C18 (+)",
        "hilicneg" = "HILIC (-)",
        "hilicpos" = "HILIC (+)"
      )
    ) +
    ggh4x::facet_nested(
      pfas_name + superclass ~ .,
      scales = "free_y",
      space = "free",
      strip = ggh4x::strip_split(c("left", "right"))
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
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey95"),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text.y.right = element_text(size = 9, angle = 0, hjust = 0),
      strip.text.y.left = element_text(size = 10)
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 6)),
      fill = "none",
      color = "none"
    )
}

pathway_plot_main <- make_pathway_plot("main")
pathway_plot_interaction <- make_pathway_plot("interaction")

ggsave(
  filename = fs::path(
    dir_figures,
    "pathway_analysis_by_untargeted_PFAS_main.png"
  ),
  plot = pathway_plot_main,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = fs::path(
    dir_figures,
    "pathway_analysis_by_untargeted_PFAS_interaction.png"
  ),
  plot = pathway_plot_interaction,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600
)

make_pathway_plot_combined_part <- function(effect_term_name) {
  ggplot(
    data = pathway_summary %>%
      dplyr::filter(effect_term == effect_term_name),
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
      values = c("c18pos" = 21, "hilicneg" = 22, "hilicpos" = 24),
      name = "Mode",
      labels = c(
        "c18pos" = "C18 (+)",
        "hilicneg" = "HILIC (-)",
        "hilicpos" = "HILIC (+)"
      )
    ) +
    ggh4x::facet_nested(
      pfas_name ~ .,
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
      shape = guide_legend(
        order = 1,
        override.aes = list(size = 6)
      ),
      fill = guide_legend(
        order = 2,
        override.aes = list(shape = 21, size = 6)
      ),
      size = guide_legend(order = 3),
      color = "none"
    )
}

pathway_plot_combined_main <- make_pathway_plot_combined_part("main")
pathway_plot_combined_interaction <- make_pathway_plot_combined_part(
  "interaction"
)
pathway_plot_combined_legend <- cowplot::get_legend(pathway_plot_combined_main)

pathway_plot_combined <- cowplot::plot_grid(
  cowplot::plot_grid(
    pathway_plot_combined_main + theme(legend.position = "none"),
    pathway_plot_combined_interaction + theme(legend.position = "none"),
    ncol = 1,
    # labels = c("Main", "Interaction"),
    # label_size = 11,
    rel_heights = c(9, 6),
    align = "v"
  ),
  pathway_plot_combined_legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

ggsave(
  filename = fs::path(
    dir_figures,
    "pathway_analysis_by_untargeted_PFAS_main_interaction.png"
  ),
  plot = pathway_plot_combined,
  width = 6,
  height = 7,
  units = "in",
  dpi = 600
)

write.csv(
  model_summary,
  fs::path(dir_results, "Figure_1_statistics_main_interaction.csv"),
  row.names = FALSE
)

write.csv(
  pathway_summary,
  fs::path(dir_results, "Figure_2_statistics_main_interaction.csv"),
  row.names = FALSE
)
