setwd("1_code")
here::i_am("2_longitudinal_metabolomics/5_2_pathway_analysis.R")

source(fs::path(here::here(), "2_longitudinal_metabolomics", "!libraries.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!directories.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!functions.R"))

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
  "site"
)

q_cutoff_sig_feat <- 0.2
q_cutoff_sig_feat_pca <- 0.000085

################
# Data Load
################

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

effect_terms <- c("main", "interaction")

################
# Pathway Analysis
################

for (i in 1:length(platforms)) {
  for (j in 1:length(chemical_variables)) {
    for (k in 1:length(effect_terms)) {
      mum_wkpath <- fs::path(
        dir_results |> dirname(),
        "1_code",
        "2_longitudinal_metabolomics",
        "temp_data",
        "Mummichog_update",
        paste0(platforms[i], "_", chemical_names[j], "_", effect_terms[k])
      )

      if (file.exists(mum_wkpath)) {
        do.call(
          function(x) unlink(x, recursive = T),
          list(list.files(mum_wkpath, full.names = T))
        )
      } else {
        dir.create(mum_wkpath)
      }

      mum_input <- lmm_random_intercept_knot1_covars_PFAS_output[[i]] %>%
        dplyr::filter(exposure == chemical_variables[j]) %>%
        dplyr::group_by(feature) %>%
        dplyr::slice(c(1:length(effect_terms))[k]) %>%
        dplyr::ungroup() %>%
        dplyr::select(feature, pvalue, statistic)
      mum_input$mz <- as.numeric(sub('_.*', '', mum_input$feature))
      mum_input$rt <- as.numeric(sub('.*_', '', mum_input$feature))
      mum_input <- mum_input %>%
        dplyr::select(mz, rt, pvalue, statistic)
      colnames(mum_input) <- c("mz", "rtime", "p-value", "t-score")

      write.table(
        mum_input,
        fs::path(
          mum_wkpath,
          paste0(
            "sum_table_",
            chemical_variables[j],
            "_",
            platforms[i],
            "_",
            effect_terms[k],
            ".txt"
          )
        ),
        quote = F,
        row.names = F,
        sep = "\t"
      )
    }
  }
}

start_time <- Sys.time()

for (i in 1:length(platforms)) {
  for (j in 1:length(chemical_variables)) {
    for (k in 1:length(effect_terms)) {
      if (
        platforms[i] == "c18pos" &
          chemical_names[j] == "PFOA" &
          effect_terms[k] == "interaction"
      ) {
        next
      }

      mum_wkpath <- fs::path(
        dir_results |> dirname(),
        "1_code",
        "2_longitudinal_metabolomics",
        "temp_data",
        "Mummichog_update",
        paste0(platforms[i], "_", chemical_names[j], "_", effect_terms[k])
      )

      feat_table_name <- platforms[i]
      if (grepl("neg", feat_table_name)) {
        mode <- "negative"
      } else {
        mode <- "positive"
      }

      cmd_code <- paste(
        "python -m mummichog.main -f",
        paste0(
          "sum_table_",
          chemical_variables[j],
          "_",
          platforms[i],
          "_",
          effect_terms[k],
          ".txt"
        ),
        "-p 500 -c ",
        q_cutoff_sig_feat_pca,
        " -m ",
        mode,
        " -z TRUE"
      )
      print(cmd_code)
      setwd(mum_wkpath)
      system(cmd_code)
    }
  }
}

end_time <- Sys.time()
end_time - start_time
