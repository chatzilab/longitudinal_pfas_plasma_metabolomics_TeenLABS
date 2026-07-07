code_dir_candidates <- unique(normalizePath(
  c("1_code", ".", ".."),
  mustWork = FALSE
))
code_dir <- code_dir_candidates[
  file.exists(file.path(
    code_dir_candidates,
    "2_longitudinal_metabolomics",
    "5_5_pathway_analysis_mix.R"
  ))
][1]
if (is.na(code_dir)) {
  stop("Cannot locate the 1_code directory from the current working directory.")
}
setwd(code_dir)
rm(code_dir, code_dir_candidates)

here::i_am("2_longitudinal_metabolomics/5_5_pathway_analysis_mix.R")

################
# Environment Setup
################

options(scipen = 999)

source(fs::path(here::here(), "2_longitudinal_metabolomics", "!libraries.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!directories.R"))
source(fs::path(here::here(), "2_longitudinal_metabolomics", "!functions.R"))


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

q_cutoff_sig_feat_pca <- 0.000085
effect_terms <- c("main", "interaction")


################
# Data Load
################
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


################
# Mummichog Input
################
for (i in seq_along(platforms)) {
  for (j in seq_along(mix_variables)) {
    for (k in seq_along(effect_terms)) {
      mum_wkpath <- fs::path(
        dirname(dir_results),
        "1_code",
        "2_longitudinal_metabolomics",
        "temp_data",
        "Mummichog_update_mix",
        paste0(platforms[i], "_", mix_slugs[j], "_", effect_terms[k])
      )

      if (file.exists(mum_wkpath)) {
        do.call(
          function(x) unlink(x, recursive = TRUE),
          list(list.files(mum_wkpath, full.names = TRUE))
        )
      } else {
        dir.create(mum_wkpath, recursive = TRUE)
      }

      mum_input <- lmm_random_intercept_knot1_covars_PFAS_mix_output[[i]] %>%
        dplyr::filter(exposure == mix_variables[j]) %>%
        dplyr::group_by(feature) %>%
        dplyr::slice(k) %>%
        dplyr::ungroup() %>%
        dplyr::select(feature, pvalue, statistic)
      mum_input$mz <- as.numeric(sub("_.*", "", mum_input$feature))
      mum_input$rt <- as.numeric(sub(".*_", "", mum_input$feature))
      mum_input <- mum_input %>%
        dplyr::select(mz, rt, pvalue, statistic)
      colnames(mum_input) <- c("mz", "rtime", "p-value", "t-score")

      write.table(
        mum_input,
        fs::path(
          mum_wkpath,
          paste0(
            "sum_table_",
            mix_variables[j],
            "_",
            platforms[i],
            "_",
            effect_terms[k],
            ".txt"
          )
        ),
        quote = FALSE,
        row.names = FALSE,
        sep = "\t"
      )
    }
  }
}


################
# Mummichog Run
################
start_time <- Sys.time()

for (i in seq_along(platforms)) {
  for (j in seq_along(mix_variables)) {
    for (k in seq_along(effect_terms)) {
      mum_wkpath <- fs::path(
        dirname(dir_results),
        "1_code",
        "2_longitudinal_metabolomics",
        "temp_data",
        "Mummichog_update_mix",
        paste0(platforms[i], "_", mix_slugs[j], "_", effect_terms[k])
      )

      mode <- ifelse(grepl("neg", platforms[i]), "negative", "positive")
      input_file <- paste0(
        "sum_table_",
        mix_variables[j],
        "_",
        platforms[i],
        "_",
        effect_terms[k],
        ".txt"
      )

      print(fs::path(mum_wkpath, input_file))
      old_wd <- getwd()
      setwd(mum_wkpath)
      system2(
        "python",
        args = c(
          "-m",
          "mummichog.main",
          "-f",
          input_file,
          "-p",
          "500",
          "-c",
          as.character(q_cutoff_sig_feat_pca),
          "-m",
          mode,
          "-z",
          "TRUE"
        )
      )
      setwd(old_wd)
    }
  }
}

end_time <- Sys.time()
end_time - start_time
