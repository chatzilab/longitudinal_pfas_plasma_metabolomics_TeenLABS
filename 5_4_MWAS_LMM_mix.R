code_dir_candidates <- unique(normalizePath(c("1_code", ".", ".."), mustWork = FALSE))
code_dir <- code_dir_candidates[
  file.exists(file.path(
    code_dir_candidates,
    "2_longitudinal_metabolomics",
    "5_4_MWAS_LMM_mix.R"
  ))
][1]
if (is.na(code_dir)) {
  stop("Cannot locate the 1_code directory from the current working directory.")
}
setwd(code_dir)
rm(code_dir, code_dir_candidates)

here::i_am("2_longitudinal_metabolomics/5_4_MWAS_LMM_mix.R")

################
# Environment Setup
################

options(scipen = 999)

source(fs::path(here::here(), "2_longitudinal_metabolomics", "!libraries.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!directories.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!functions.R"))

library(nlme)
library(doParallel)
library(foreach)
library(kableExtra)


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

q_cutoff_sig_feat_pca <- 0.000078


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

lmm_random_intercept_knot1_covars_PFAS_mix_output <- list()
for (i in seq_along(platforms)) {
  print(platforms[i])
  model_df <- met_fts_final[[i]] %>%
    dplyr::select(-key, -visit)
  covariates_df <- tl_final[, covariates]
  key <- tl_final$key
  visit <- tl_final$visit_year

  output_by_mix <- list()
  for (k in seq_along(mix_variables)) {
    print(mix_variables[k])
    x <- tl_final[, mix_variables[k]]
    cl <- makeCluster(detectCores() - 2)
    registerDoParallel(cl)
    model_output <- foreach(
      j = seq_len(ncol(model_df)),
      .combine = rbind,
      .packages = c("nlme")
    ) %dopar%
      fit_lmm_random_intercept_covars_JT(
        x = x,
        y = model_df[, j],
        key = key,
        visit = visit,
        covariates_df = covariates_df
      )
    stopCluster(cl)

    model_output <- model_output %>%
      as.data.frame() %>%
      dplyr::mutate_at(
        dplyr::vars(coef, se, statistic, pvalue, pvalue_interaction, L_ratio),
        as.numeric
      )

    model_output$feature <- rep(colnames(model_df), each = 2)
    model_output$p_adjust_interaction <- p.adjust(
      model_output$pvalue_interaction,
      method = "BH"
    )

    output_by_mix[[k]] <- model_output
    rm(model_output, x)
  }

  lmm_random_intercept_knot1_covars_PFAS_mix_output[[i]] <- output_by_mix %>%
    dplyr::bind_rows()

  rm(output_by_mix, model_df)
}
names(lmm_random_intercept_knot1_covars_PFAS_mix_output) <- platforms

save(
  lmm_random_intercept_knot1_covars_PFAS_mix_output,
  file = fs::path(
    dirname(dir_results),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "model_statistics_untargeted_PFAS_mix_longitudinal_metabo_06012026.RData"
  )
)

end_time <- Sys.time()
end_time - start_time


################
# Summary
################
get_sig_features_by_mix_effect <- function(df, var, threshold) {
  sig_rows <- df %>%
    dplyr::mutate(effect_term = rep(c("main", "interaction"), length.out = dplyr::n())) %>%
    dplyr::group_by(exposure, effect_term) %>%
    dplyr::mutate(p_adjust = p.adjust(pvalue, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data[[var]] < threshold)

  feature_list_full <- setNames(vector("list", length(mix_variables)), mix_variables)
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

sig_list <- lapply(
  lmm_random_intercept_knot1_covars_PFAS_mix_output,
  function(df) get_sig_by_chemical(df, "pvalue", q_cutoff_sig_feat_pca)
) %>%
  purrr::modify2(
    names(lmm_random_intercept_knot1_covars_PFAS_mix_output),
    ~ .x %>% dplyr::mutate(mode = .y)
  ) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    exposure = factor(exposure, levels = mix_variables),
    mixture = mix_names[match(exposure, mix_variables)]
  ) %>%
  dplyr::arrange(mode, exposure)

sig_list %>%
  dplyr::select(mode, mixture, pvalue_main, pvalue_interaction) %>%
  kbl(
    col.names = c(
      "Mode",
      "PFAS mixture",
      "# of Sig main (nominal p < 0.000078)",
      "# of Sig interaction (nominal p < 0.000078)"
    )
  ) %>%
  kable_paper("hover")


################
# Significant Feature Output
################
sig_feature_list <- lapply(
  lmm_random_intercept_knot1_covars_PFAS_mix_output,
  function(df) get_sig_features_by_mix_effect(df, "pvalue", q_cutoff_sig_feat_pca)
)

saveRDS(
  sig_feature_list,
  file = fs::path(
    dirname(dir_results),
    "1_code",
    "2_longitudinal_metabolomics",
    "temp_data",
    "sig_feature_by_mix_pca_corrected_main_interaction_06012026.rds"
  )
)

