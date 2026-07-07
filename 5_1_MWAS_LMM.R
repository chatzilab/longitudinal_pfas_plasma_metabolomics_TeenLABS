setwd("1_code")
here::i_am("2_longitudinal_metabolomics/5_1_MWAS_LMM.R")

source(fs::path(here::here(), "2_longitudinal_metabolomics", "!libraries.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!directories.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!functions.R"))

library(nlme)
library(doParallel)
library(foreach)
library(kableExtra)

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

################
# Model Parameters
################

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

mix_variables <- c(
  "pfas_sulfonic_acids",
  "pfas_carboxylic_acids"
)

################
# Data Load
################

load(
  fs::path(
    dir_data,
    "tl_long_met_untargeted_pfas_workdata.RData"
  )
)

platforms <- names(met_fts_final)

################
# Main Model
################

start_time <- Sys.time()

lmm_random_intercept_knot1_covars_PFAS_output <- {}
for (i in 1:length(platforms)) {
  print(i)
  model_df <- met_fts_final[[i]] %>%
    dplyr::select(-key, -visit)
  covariates_df <- tl_final[, covariates]
  key <- tl_final$key
  visit <- tl_final$visit_year

  output_by_chemical <- {}
  for (k in 1:length(chemical_variables)) {
    print(chemical_variables[k])
    x <- tl_final[, chemical_variables[k]]
    cl <- makeCluster(detectCores() - 2)
    registerDoParallel(cl)
    model_output <- foreach(
      j = 1:dim(model_df)[2],
      .combine = rbind,
      .packages = c("nlme")
    ) %dopar%
      fit_lmm_random_intercept_covars_JT(
        x = x,
        y = model_df[, j],
        key = key,
        visit = visit,
        covariates_df
      )
    stopCluster(cl)

    model_output <- model_output %>%
      as.data.frame %>%
      mutate_at(
        vars(coef, se, statistic, pvalue, pvalue_interaction, L_ratio),
        as.numeric
      )

    model_output$feature <- rep(colnames(model_df), each = 2)

    model_output$p_adjust_interaction <- p.adjust(
      model_output$pvalue_interaction,
      method = "BH"
    )

    output_by_chemical[[k]] <- model_output

    rm(model_output, x)
  }

  lmm_random_intercept_knot1_covars_PFAS_output[[i]] <- output_by_chemical |>
    bind_rows()

  rm(output_by_chemical, model_df)
}
names(lmm_random_intercept_knot1_covars_PFAS_output) <- platforms

save(
  lmm_random_intercept_knot1_covars_PFAS_output,
  file = fs::path(
    dir_results |> dirname(),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "model_statistics_untargeted_PFAS_longitudinal_metabo_04302026.RData"
  )
)

end_time <- Sys.time()
end_time - start_time

# load(
#   fs::path(
#     dir_results |> dirname(),
#     "1_code",
#     "2_longitudinal_metabolomics",
#     "temp_data",
#     "model_statistics_untargeted_PFAS_longitudinal_metabo_04302026.RData"
#   )
# )

#################
# Summary
#################
## Raw p-values \< 0.05

sig_list <- lapply(
  lmm_random_intercept_knot1_covars_PFAS_output,
  function(df) get_sig_by_chemical(df, "pvalue", 0.05)
) |>
  modify2(
    names(lmm_random_intercept_knot1_covars_PFAS_output),
    ~ .x %>%
      dplyr::mutate(mode = .y)
  ) |>
  bind_rows() %>%
  dplyr::mutate(
    exposure = factor(exposure, levels = chemical_variables)
  ) %>%
  dplyr::arrange(mode, exposure)

sig_list$PFAS <- chemical_names[match(sig_list$exposure, chemical_variables)]

sig_list %>%
  dplyr::select(mode, PFAS, pvalue_main, pvalue_interaction) |>
  kbl(
    col.names = c(
      "Mode",
      "PFAS chemicals",
      "# of Sig main (p < 0.05)",
      "# of Sig interaction (p < 0.05)"
    )
  ) |>
  kable_paper("hover")

# PCA-corrected p-values

# p_threshold_pca <- modify(
#   met_fts_final,
#   ~ {
#     pca_mwas <- prcomp(.x[, -c(1:2)])
#     eigenvalues_mwas <- pca_mwas$sdev^2
#     M_eff_mwas <- sum(eigenvalues_mwas > 1)
#     0.05 / M_eff_mwas
#   }
# )
p_threshold_pca <- 8.576329e-05

sig_list <- lapply(
  lmm_random_intercept_knot1_covars_PFAS_output,
  function(df) get_sig_by_chemical(df, "pvalue", 0.000085)
) |>
  modify2(
    names(lmm_random_intercept_knot1_covars_PFAS_output),
    ~ .x %>%
      dplyr::mutate(mode = .y)
  ) |>
  bind_rows() %>%
  dplyr::mutate(
    exposure = factor(exposure, levels = chemical_variables)
  ) %>%
  dplyr::arrange(mode, exposure)

sig_list$PFAS <- chemical_names[match(sig_list$exposure, chemical_variables)]

sig_list %>%
  dplyr::select(mode, PFAS, pvalue_main, pvalue_interaction) |>
  kbl(
    col.names = c(
      "Mode",
      "PFAS chemicals",
      "# of Sig main (nominal p < 0.000085)",
      "# of Sig interation (nominal p < 0.000085)"
    )
  ) |>
  kable_paper("hover")

####################
# Save Sig. Metabolites
####################

sig_feature_list <- lapply(
  lmm_random_intercept_knot1_covars_PFAS_output,
  function(df) {
    get_sig_features_by_chemical(df, "pvalue", 0.000085)
  }
)

saveRDS(
  sig_feature_list,
  file = fs::path(
    dir_results |> dirname(),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "sig_JT_feature_by_chemical_pca_corrected_04302026.rds"
  )
)

####################
# Summarize Sig. Feature Overlap by PFAS
####################

sig_feature_detail_by_pfas <- purrr::imap_dfr(
  lmm_random_intercept_knot1_covars_PFAS_output,
  function(df, mode_name) {
    df %>%
      dplyr::mutate(
        mode = mode_name,
        effect_term = rep(c("main", "interaction"), length.out = dplyr::n())
      ) %>%
      dplyr::group_by(exposure, effect_term) %>%
      dplyr::mutate(
        p_adjust = p.adjust(pvalue, method = "BH")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(pvalue < p_threshold_pca) %>%
      dplyr::transmute(
        mode,
        exposure,
        PFAS = chemical_names[match(exposure, chemical_variables)],
        effect_term,
        feature
      )
  }
) %>%
  dplyr::distinct()

sig_feature_count_by_pfas <- tidyr::crossing(
  mode = names(lmm_random_intercept_knot1_covars_PFAS_output),
  exposure = chemical_variables,
  effect_term = c("main", "interaction")
) %>%
  dplyr::left_join(
    sig_feature_detail_by_pfas %>%
      dplyr::group_by(mode, exposure, effect_term) %>%
      dplyr::summarise(
        n_sig_features = dplyr::n_distinct(feature),
        .groups = "drop"
      ),
    by = c("mode", "exposure", "effect_term")
  ) %>%
  dplyr::mutate(
    n_sig_features = tidyr::replace_na(n_sig_features, 0L),
    PFAS = chemical_names[match(exposure, chemical_variables)]
  ) %>%
  dplyr::mutate(
    exposure = factor(exposure, levels = chemical_variables),
    PFAS = factor(PFAS, levels = chemical_names),
    effect_term = factor(effect_term, levels = c("main", "interaction"))
  ) %>%
  dplyr::arrange(mode, effect_term, exposure)

sig_feature_count_by_pfas %>%
  kbl(
    col.names = c(
      "Mode",
      "Exposure",
      "PFAS",
      "Effect term",
      "# of significant features"
    )
  ) %>%
  kable_paper("hover")

sig_feature_pairwise_overlap <- purrr::imap_dfr(
  lmm_random_intercept_knot1_covars_PFAS_output,
  function(df, mode_name) {
    purrr::map_dfr(
      c("main", "interaction"),
      function(effect_name) {
        chem_pairs <- utils::combn(chemical_variables, 2, simplify = FALSE)

        purrr::map_dfr(
          chem_pairs,
          function(pair) {
            features_1 <- sig_feature_detail_by_pfas %>%
              dplyr::filter(
                mode == mode_name,
                exposure == pair[1],
                effect_term == effect_name
              ) %>%
              dplyr::pull(feature) %>%
              unique()
            features_2 <- sig_feature_detail_by_pfas %>%
              dplyr::filter(
                mode == mode_name,
                exposure == pair[2],
                effect_term == effect_name
              ) %>%
              dplyr::pull(feature) %>%
              unique()
            overlap_features <- intersect(features_1, features_2)
            union_features <- union(features_1, features_2)

            tibble::tibble(
              mode = mode_name,
              effect_term = effect_name,
              exposure_1 = pair[1],
              exposure_2 = pair[2],
              PFAS_1 = chemical_names[match(pair[1], chemical_variables)],
              PFAS_2 = chemical_names[match(pair[2], chemical_variables)],
              n_sig_features_1 = length(features_1),
              n_sig_features_2 = length(features_2),
              n_overlap_features = length(overlap_features),
              n_union_features = length(union_features),
              jaccard_index = dplyr::if_else(
                length(union_features) > 0,
                length(overlap_features) / length(union_features),
                NA_real_
              ),
              overlap_features = paste(overlap_features, collapse = "; ")
            )
          }
        )
      }
    )
  }
) %>%
  dplyr::mutate(
    effect_term = factor(effect_term, levels = c("main", "interaction")),
    PFAS_1 = factor(PFAS_1, levels = chemical_names),
    PFAS_2 = factor(PFAS_2, levels = chemical_names)
  ) %>%
  dplyr::arrange(mode, effect_term, PFAS_1, PFAS_2)

sig_feature_pairwise_overlap %>%
  dplyr::select(
    mode,
    effect_term,
    PFAS_1,
    PFAS_2,
    n_sig_features_1,
    n_sig_features_2,
    n_overlap_features,
    n_union_features,
    jaccard_index
  ) %>%
  kbl(
    col.names = c(
      "Mode",
      "Effect term",
      "PFAS 1",
      "PFAS 2",
      "# sig. features PFAS 1",
      "# sig. features PFAS 2",
      "# overlapping features",
      "# union features",
      "Jaccard index"
    ),
    digits = 3
  ) %>%
  kable_paper("hover")

sig_feature_overlap_matrices <- purrr::imap(
  lmm_random_intercept_knot1_covars_PFAS_output,
  function(df, mode_name) {
    purrr::map(
      c("main", "interaction"),
      function(effect_name) {
        overlap_mat <- matrix(
          0,
          nrow = length(chemical_variables),
          ncol = length(chemical_variables),
          dimnames = list(chemical_names, chemical_names)
        )

        for (chem_1 in chemical_variables) {
          for (chem_2 in chemical_variables) {
            features_1 <- sig_feature_detail_by_pfas %>%
              dplyr::filter(
                mode == mode_name,
                exposure == chem_1,
                effect_term == effect_name
              ) %>%
              dplyr::pull(feature) %>%
              unique()
            features_2 <- sig_feature_detail_by_pfas %>%
              dplyr::filter(
                mode == mode_name,
                exposure == chem_2,
                effect_term == effect_name
              ) %>%
              dplyr::pull(feature) %>%
              unique()
            overlap_mat[
              chemical_names[match(chem_1, chemical_variables)],
              chemical_names[match(chem_2, chemical_variables)]
            ] <- length(intersect(features_1, features_2))
          }
        }

        overlap_mat
      }
    ) %>%
      stats::setNames(c("main", "interaction"))
  }
)

sig_feature_count_by_pfas_all_modes <- tidyr::crossing(
  exposure = chemical_variables,
  effect_term = c("main", "interaction")
) %>%
  dplyr::left_join(
    sig_feature_detail_by_pfas %>%
      dplyr::group_by(exposure, effect_term) %>%
      dplyr::summarise(
        n_sig_features = dplyr::n_distinct(feature),
        modes_with_sig_features = paste(sort(unique(mode)), collapse = "; "),
        .groups = "drop"
      ),
    by = c("exposure", "effect_term")
  ) %>%
  dplyr::mutate(
    mode = "All modes",
    n_sig_features = tidyr::replace_na(n_sig_features, 0L),
    modes_with_sig_features = tidyr::replace_na(modes_with_sig_features, ""),
    PFAS = chemical_names[match(exposure, chemical_variables)],
    exposure = factor(exposure, levels = chemical_variables),
    PFAS = factor(PFAS, levels = chemical_names),
    effect_term = factor(effect_term, levels = c("main", "interaction"))
  ) %>%
  dplyr::select(
    mode,
    exposure,
    PFAS,
    effect_term,
    n_sig_features,
    modes_with_sig_features
  ) %>%
  dplyr::arrange(effect_term, exposure)

sig_feature_count_by_pfas_all_modes %>%
  kbl(
    col.names = c(
      "Mode",
      "Exposure",
      "PFAS",
      "Effect term",
      "# of significant features",
      "Modes with sig. features"
    )
  ) %>%
  kable_paper("hover")

supp_table_feature_count_by_pfas_title <- paste(
  "Supplementary Table X. Number of PFAS-associated metabolomic features",
  "by PFAS and effect term."
)

supp_table_feature_count_by_pfas_note <- paste0(
  "Significant features were selected using the PCA-corrected threshold ",
  "(p < ", format(p_threshold_pca, scientific = TRUE, digits = 3), "). ",
  "Features were collapsed across ionization modes before counting. ",
  "Shared features are features associated with the listed PFAS and at least ",
  "one other PFAS for the same effect term. Unique features are associated ",
  "only with the listed PFAS for the same effect term."
)

supp_table_feature_count_by_pfas_all_modes <- purrr::map_dfr(
  c("main", "interaction"),
  function(effect_name) {
    sig_by_effect <- sig_feature_detail_by_pfas %>%
      dplyr::filter(effect_term == effect_name) %>%
      dplyr::distinct(exposure, PFAS, feature)

    feature_pfas_counts <- sig_by_effect %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(
        n_pfas = dplyr::n_distinct(exposure),
        .groups = "drop"
      )

    tidyr::tibble(
      exposure = chemical_variables,
      PFAS = chemical_names
    ) %>%
      dplyr::left_join(
        sig_by_effect %>%
          dplyr::left_join(feature_pfas_counts, by = "feature") %>%
          dplyr::group_by(exposure, PFAS) %>%
          dplyr::summarise(
            `Total associated features` = dplyr::n_distinct(feature),
            `Associated features also found for other PFAS` =
              dplyr::n_distinct(feature[n_pfas > 1]),
            `Unique associated features` =
              dplyr::n_distinct(feature[n_pfas == 1]),
            .groups = "drop"
          ),
        by = c("exposure", "PFAS")
      ) %>%
      dplyr::mutate(
        `Effect term` = dplyr::case_when(
          effect_name == "main" ~ "Main effect",
          effect_name == "interaction" ~ "PFAS x time interaction",
          TRUE ~ effect_name
        ),
        dplyr::across(
          c(
            `Total associated features`,
            `Associated features also found for other PFAS`,
            `Unique associated features`
          ),
          ~ tidyr::replace_na(.x, 0L)
        )
      ) %>%
      dplyr::select(
        `Effect term`,
        PFAS,
        `Total associated features`,
        `Associated features also found for other PFAS`,
        `Unique associated features`
      )
  }
) %>%
  dplyr::mutate(
    PFAS = factor(PFAS, levels = chemical_names),
    `Effect term` = factor(
      `Effect term`,
      levels = c("Main effect", "PFAS x time interaction")
    )
  ) %>%
  dplyr::arrange(`Effect term`, PFAS)

supp_table_feature_count_by_pfas_all_modes %>%
  kbl(
    caption = paste(
      supp_table_feature_count_by_pfas_title,
      supp_table_feature_count_by_pfas_note
    ),
    booktabs = TRUE
  ) %>%
  kable_paper("hover")

sig_feature_pairwise_overlap_all_modes <- purrr::map_dfr(
  c("main", "interaction"),
  function(effect_name) {
    chem_pairs <- utils::combn(chemical_variables, 2, simplify = FALSE)

    purrr::map_dfr(
      chem_pairs,
      function(pair) {
        features_1 <- sig_feature_detail_by_pfas %>%
          dplyr::filter(
            exposure == pair[1],
            effect_term == effect_name
          ) %>%
          dplyr::pull(feature) %>%
          unique()
        features_2 <- sig_feature_detail_by_pfas %>%
          dplyr::filter(
            exposure == pair[2],
            effect_term == effect_name
          ) %>%
          dplyr::pull(feature) %>%
          unique()
        overlap_features <- intersect(features_1, features_2)
        union_features <- union(features_1, features_2)

        tibble::tibble(
          mode = "All modes",
          effect_term = effect_name,
          exposure_1 = pair[1],
          exposure_2 = pair[2],
          PFAS_1 = chemical_names[match(pair[1], chemical_variables)],
          PFAS_2 = chemical_names[match(pair[2], chemical_variables)],
          n_sig_features_1 = length(features_1),
          n_sig_features_2 = length(features_2),
          n_overlap_features = length(overlap_features),
          n_union_features = length(union_features),
          jaccard_index = dplyr::if_else(
            length(union_features) > 0,
            length(overlap_features) / length(union_features),
            NA_real_
          ),
          overlap_features = paste(overlap_features, collapse = "; ")
        )
      }
    )
  }
) %>%
  dplyr::mutate(
    effect_term = factor(effect_term, levels = c("main", "interaction")),
    PFAS_1 = factor(PFAS_1, levels = chemical_names),
    PFAS_2 = factor(PFAS_2, levels = chemical_names)
  ) %>%
  dplyr::arrange(effect_term, PFAS_1, PFAS_2)

sig_feature_pairwise_overlap_all_modes %>%
  dplyr::select(
    mode,
    effect_term,
    PFAS_1,
    PFAS_2,
    n_sig_features_1,
    n_sig_features_2,
    n_overlap_features,
    n_union_features,
    jaccard_index
  ) %>%
  kbl(
    col.names = c(
      "Mode",
      "Effect term",
      "PFAS 1",
      "PFAS 2",
      "# sig. features PFAS 1",
      "# sig. features PFAS 2",
      "# overlapping features",
      "# union features",
      "Jaccard index"
    ),
    digits = 3
  ) %>%
  kable_paper("hover")

supp_table_pairwise_overlap_title <- paste(
  "Supplementary Table X. Pairwise overlap of PFAS-associated",
  "metabolomic features across ionization modes."
)

supp_table_pairwise_overlap_note <- paste0(
  "Significant features were selected using the PCA-corrected threshold ",
  "(p < ", format(p_threshold_pca, scientific = TRUE, digits = 3), "). ",
  "Features were collapsed across ionization modes before computing pairwise ",
  "overlap. Jaccard index is the number of overlapping features divided by ",
  "the number of unique features associated with either PFAS in the pair."
)

supp_table_sig_feature_pairwise_overlap_all_modes <-
  sig_feature_pairwise_overlap_all_modes %>%
  dplyr::transmute(
    `Effect term` = dplyr::case_when(
      effect_term == "main" ~ "Main effect",
      effect_term == "interaction" ~ "PFAS x time interaction",
      TRUE ~ as.character(effect_term)
    ),
    `PFAS pair` = paste(PFAS_1, PFAS_2, sep = " vs. "),
    `PFAS 1` = as.character(PFAS_1),
    `PFAS 2` = as.character(PFAS_2),
    `N significant features, PFAS 1` = n_sig_features_1,
    `N significant features, PFAS 2` = n_sig_features_2,
    `N overlapping features` = n_overlap_features,
    `N unique features across pair` = n_union_features,
    `Jaccard index` = round(jaccard_index, 3),
    `Overlapping feature IDs` = dplyr::if_else(
      overlap_features == "",
      NA_character_,
      overlap_features
    )
  ) %>%
  dplyr::arrange(`Effect term`, `PFAS 1`, `PFAS 2`)

supp_table_sig_feature_pairwise_overlap_all_modes %>%
  kbl(
    caption = paste(supp_table_pairwise_overlap_title, supp_table_pairwise_overlap_note),
    booktabs = TRUE,
    longtable = TRUE
  ) %>%
  kable_paper("hover")

sig_feature_overlap_matrices_all_modes <- purrr::map(
  c("main", "interaction"),
  function(effect_name) {
    overlap_mat <- matrix(
      0,
      nrow = length(chemical_variables),
      ncol = length(chemical_variables),
      dimnames = list(chemical_names, chemical_names)
    )

    for (chem_1 in chemical_variables) {
      for (chem_2 in chemical_variables) {
        features_1 <- sig_feature_detail_by_pfas %>%
          dplyr::filter(
            exposure == chem_1,
            effect_term == effect_name
          ) %>%
          dplyr::pull(feature) %>%
          unique()
        features_2 <- sig_feature_detail_by_pfas %>%
          dplyr::filter(
            exposure == chem_2,
            effect_term == effect_name
          ) %>%
          dplyr::pull(feature) %>%
          unique()
        overlap_mat[
          chemical_names[match(chem_1, chemical_variables)],
          chemical_names[match(chem_2, chemical_variables)]
        ] <- length(intersect(features_1, features_2))
      }
    }

    overlap_mat
  }
) %>%
  stats::setNames(c("main", "interaction"))

sig_feature_count_by_pfas_main <- sig_feature_count_by_pfas %>%
  dplyr::filter(effect_term == "main") %>%
  dplyr::select(-effect_term)

sig_feature_count_by_pfas_interaction <- sig_feature_count_by_pfas %>%
  dplyr::filter(effect_term == "interaction") %>%
  dplyr::select(-effect_term)

sig_feature_count_by_pfas_all_modes_main <- sig_feature_count_by_pfas_all_modes %>%
  dplyr::filter(effect_term == "main") %>%
  dplyr::select(-effect_term)

sig_feature_count_by_pfas_all_modes_interaction <- sig_feature_count_by_pfas_all_modes %>%
  dplyr::filter(effect_term == "interaction") %>%
  dplyr::select(-effect_term)

sig_feature_pairwise_overlap_main <- sig_feature_pairwise_overlap %>%
  dplyr::filter(effect_term == "main") %>%
  dplyr::select(-effect_term)

sig_feature_pairwise_overlap_interaction <- sig_feature_pairwise_overlap %>%
  dplyr::filter(effect_term == "interaction") %>%
  dplyr::select(-effect_term)

sig_feature_pairwise_overlap_all_modes_main <- sig_feature_pairwise_overlap_all_modes %>%
  dplyr::filter(effect_term == "main") %>%
  dplyr::select(-effect_term)

sig_feature_pairwise_overlap_all_modes_interaction <- sig_feature_pairwise_overlap_all_modes %>%
  dplyr::filter(effect_term == "interaction") %>%
  dplyr::select(-effect_term)

sig_feature_count_by_pfas_all_modes_main %>%
  kbl(
    col.names = c(
      "Mode",
      "Exposure",
      "PFAS",
      "# of significant main-effect features",
      "Modes with sig. features"
    )
  ) %>%
  kable_paper("hover")

sig_feature_count_by_pfas_all_modes_interaction %>%
  kbl(
    col.names = c(
      "Mode",
      "Exposure",
      "PFAS",
      "# of significant interaction features",
      "Modes with sig. features"
    )
  ) %>%
  kable_paper("hover")

sig_feature_pairwise_overlap_all_modes_main %>%
  dplyr::select(
    mode,
    PFAS_1,
    PFAS_2,
    n_sig_features_1,
    n_sig_features_2,
    n_overlap_features,
    n_union_features,
    jaccard_index
  ) %>%
  kbl(
    col.names = c(
      "Mode",
      "PFAS 1",
      "PFAS 2",
      "# sig. features PFAS 1",
      "# sig. features PFAS 2",
      "# overlapping main-effect features",
      "# union main-effect features",
      "Jaccard index"
    ),
    digits = 3
  ) %>%
  kable_paper("hover")

sig_feature_pairwise_overlap_all_modes_interaction %>%
  dplyr::select(
    mode,
    PFAS_1,
    PFAS_2,
    n_sig_features_1,
    n_sig_features_2,
    n_overlap_features,
    n_union_features,
    jaccard_index
  ) %>%
  kbl(
    col.names = c(
      "Mode",
      "PFAS 1",
      "PFAS 2",
      "# sig. features PFAS 1",
      "# sig. features PFAS 2",
      "# overlapping interaction features",
      "# union interaction features",
      "Jaccard index"
    ),
    digits = 3
  ) %>%
  kable_paper("hover")

readr::write_csv(
  supp_table_feature_count_by_pfas_all_modes,
  file = fs::path(
    dir_results |> dirname(),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "supp_table_feature_count_by_pfas_all_modes_pca_corrected_04302026.csv"
  )
)

readr::write_csv(
  supp_table_sig_feature_pairwise_overlap_all_modes,
  file = fs::path(
    dir_results |> dirname(),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "supp_table_sig_feature_pairwise_overlap_all_modes_pca_corrected_04302026.csv"
  )
)
