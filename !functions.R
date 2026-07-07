fit_lmm_random_intercept_knot1_covars <- function(
  x,
  y,
  key,
  visit,
  covariates_df
) {
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  model <- as.formula(paste(
    "y ~ x*bSpline(visit, knots = 1, degree = 1)",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))
  reduced_model <- as.formula(paste(
    "y ~ x + bSpline(visit, knots = 1, degree = 1)",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))
  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  coef_table <- summary(fit)$tTable
  start_index <- grep("x:bSpline", rownames(coef_table))[1]
  end_index <- start_index + 1
  output <- coef_table[c(2, start_index:end_index), c(1, 2, 4, 5)]

  fit_reduced <- lme(
    reduced_model,
    random = ~ 1 | key,
    data = df,
    na.action = na.omit
  )
  overall_sig <- anova(
    update(fit, . ~ ., method = "ML"),
    update(fit_reduced, . ~ ., method = "ML")
  )
  output <- cbind(output, c(overall_sig$`p-value`[2], rep(NA, time = 2)))
  output <- cbind(output, c(overall_sig$L.Ratio[2], rep(NA, time = 2)))
  colnames(output) <- c(
    "coef",
    "se",
    "statistic",
    "pvalue",
    "pvalue_interaction",
    "L_ratio"
  )

  output <- as.data.frame(output)

  output$exposure <- names(x)

  return(output)
}

# Linear mixed-effect models without spline of time
fit_lmm_random_intercept_covars <- function(x, y, key, visit, covariates_df) {
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
  reduced_model <- as.formula(paste(
    "y ~ x + visit",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))
  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  coef_table <- summary(fit)$tTable
  start_index <- grep("x", rownames(coef_table))[1]
  end_index <- grep("x:", rownames(coef_table))[1]
  output <- coef_table[c(start_index, end_index), c(1, 2, 4, 5)]

  fit_reduced <- tryCatch(
    {
      lme(reduced_model, random = ~ 1 | key, data = df, na.action = na.omit)
    },
    error = function(e) {
      NA
    }
  )

  fit_reduced <- tryCatch(
    {
      update(fit_reduced, . ~ ., method = "ML")
    },
    error = function(e) {
      NA
    }
  )

  if (length(fit_reduced) == 1) {
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
  } else {
    overall_sig <- anova(
      update(fit, . ~ ., method = "ML"),
      update(fit_reduced, . ~ ., method = "ML")
    )
    output <- cbind(output, c(overall_sig$`p-value`[2], rep(NA, time = 1)))
    output <- cbind(output, c(overall_sig$L.Ratio[2], rep(NA, time = 1)))
  }

  colnames(output) <- c(
    "coef",
    "se",
    "statistic",
    "pvalue",
    "pvalue_interaction",
    "L_ratio"
  )

  output <- as.data.frame(output)

  output$exposure <- names(x)

  return(output)
}

# get strings that appear at least twice in a list of elements.
get_intersect_string <- function(string_list, n) {
  word_counts <- table(unlist(string_list))
  common_words <- names(word_counts[word_counts >= n])
  return(common_words)
}

# get annotations of metabolic features
get_annotation <- function(feature_vec, annotation_df) {
  match_index <- match(feature_vec, annotation_df$`colnames(x)`)
  names(feature_vec) <- annotation_df$Metabolite_Name[match_index]
  feature_vec <- feature_vec[
    !is.na(names(feature_vec)) & names(feature_vec) != ""
  ]

  return(feature_vec)
}

# get significant count
get_sig_by_chemical <- function(df, var = "p_adjust", threshold = 0.05) {
  df %>%
    dplyr::mutate(
      term_type = rep(c("main", "interaction"), length.out = dplyr::n())
    ) %>%
    dplyr::group_by(exposure, term_type) %>%
    dplyr::mutate(
      p_adjust = p.adjust(pvalue, method = "BH")
    ) %>%
    dplyr::summarise(
      n_sig = sum(.data[[var]] < threshold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = term_type,
      values_from = n_sig,
      names_glue = paste0(var, "_{term_type}")
    ) %>%
    dplyr::select(
      exposure,
      !!paste0(var, "_main"),
      !!paste0(var, "_interaction")
    )
}

# get significant count by interaction
get_sig_JT_by_chemical <- function(df, var, threshold) {
  # Select only the interaction term rows (every second row)
  interaction_rows <- df[seq(1, nrow(df), by = 2), ]
  # interaction_rows <- interaction_rows %>%
  #   dplyr::group_by(exposure) %>%
  #   dplyr::mutate(
  #     p_adjust = p.adjust(pvalue, method = "BH")) %>%
  #   dplyr::ungroup()
  # Count how many have p_value < threshold
  aggregate(
    as.formula(paste(var, "~ exposure")),
    data = interaction_rows,
    FUN = function(x) sum(x < threshold, na.rm = T)
  )
}

# get significant feature m/z values
get_sig_features_by_chemical <- function(df, var, threshold) {
  sig_rows <- df %>%
    dplyr::mutate(
      term_type = rep(c("main", "interaction"), length.out = dplyr::n()),
      feature_mz = as.numeric(sub("_.*", "", feature))
    ) %>%
    dplyr::group_by(exposure, term_type) %>%
    dplyr::mutate(
      p_adjust = p.adjust(pvalue, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data[[var]] < threshold)

  if (exists("mix_variables") && any(df$exposure %in% mix_variables)) {
    exposure_variables <- mix_variables
  } else if (exists("chemical_variables")) {
    exposure_variables <- chemical_variables
  } else {
    exposure_variables <- unique(df$exposure)
  }

  feature_list_full <- setNames(
    vector("list", length(exposure_variables)),
    exposure_variables
  )
  for (chem in exposure_variables) {
    feature_list_full[[chem]] <- list(
      main = sig_rows$feature_mz[
        sig_rows$exposure == chem & sig_rows$term_type == "main"
      ],
      interaction = sig_rows$feature_mz[
        sig_rows$exposure == chem & sig_rows$term_type == "interaction"
      ]
    )
  }

  return(feature_list_full)
}


# get significant feature names by joint test
get_sig_JT_features_by_chemical <- function(df, var, threshold) {
  # Subset to interaction term rows (every 2nd row)
  interaction_rows <- df[seq(1, nrow(df), by = 2), ]

  # # Filter for significance
  # interaction_rows <- interaction_rows %>%
  #   dplyr::group_by(exposure) %>%
  #   dplyr::mutate(
  #     p_adjust = p.adjust(pvalue, method = "BH")) %>%
  #   dplyr::ungroup()
  sig_rows <- interaction_rows[interaction_rows[[var]] < threshold, ]

  # Group by chemical and collect feature names
  feature_list <- split(sig_rows$feature, sig_rows$exposure)

  # Ensure all chemicals are represented, even with no sigs
  if (any(names(feature_list) %in% mix_variables)) {
    chemical_variables <- mix_variables
  }

  feature_list_full <- setNames(
    vector("list", length(chemical_variables)),
    chemical_variables
  )
  for (chem in chemical_variables) {
    feature_list_full[[chem]] <- feature_list[[chem]] %||% character(0)
  }

  return(feature_list_full)
}

# Linear mixed-effect models without spline of time
# with confidence interval
fit_lmm_random_intercept_covars_CI <- function(
  x,
  y,
  key,
  visit,
  covariates_df
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
  reduced_model <- as.formula(paste(
    "y ~ x + visit",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))
  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  coef_table <- bind_cols(
    summary(fit)$tTable,
    intervals(fit, which = "fixed")$fixed
  )
  coef_table$var_name <- rownames(summary(fit)$tTable)
  start_index <- grep("x", coef_table$var_name)[1]
  end_index <- grep("x:", coef_table$var_name)[1]
  output <- coef_table[c(start_index, end_index), c(1, 2, 4, 5, 6, 8)]

  fit_reduced <- tryCatch(
    {
      lme(reduced_model, random = ~ 1 | key, data = df, na.action = na.omit)
    },
    error = function(e) {
      NA
    }
  )

  fit_reduced <- tryCatch(
    {
      update(fit_reduced, . ~ ., method = "ML")
    },
    error = function(e) {
      NA
    }
  )

  if (length(fit_reduced) == 1) {
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
  } else {
    overall_sig <- anova(
      update(fit, . ~ ., method = "ML"),
      update(fit_reduced, . ~ ., method = "ML")
    )
    output <- cbind(output, c(overall_sig$`p-value`[2], rep(NA, time = 1)))
    output <- cbind(output, c(overall_sig$L.Ratio[2], rep(NA, time = 1)))
  }

  colnames(output) <- c(
    "coef",
    "se",
    "statistic",
    "pvalue",
    "conf_low",
    "conf_high",
    "pvalue_interaction",
    "L_ratio"
  )

  output <- as.data.frame(output)

  output$exposure <- names(x)

  output$feature <- names(y)

  output$beta_type <- c("main_effect", "interaction")

  return(output)
}

# Linear mixed-effect models for metabolite trajectory
fit_lmm_random_intercept_metabo <- function(y, metabo_name, key, visit) {
  df <- as.data.frame(cbind(key, y, visit))
  names(df) <- c("key", "y", "visit")
  df$y <- as.numeric(df$y)
  model <- as.formula(paste("y ~ visit"))
  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  coef_table <- summary(fit)$tTable

  # Extract the second row (visit effect) and coerce it to a 1-row data frame
  values <- coef_table[2, c(1, 2, 4, 5)]
  output <- as.data.frame(t(values)) # transpose to get it as a row
  output$feature <- metabo_name # add feature name

  # Set column names
  colnames(output) <- c("coef", "se", "statistic", "pvalue", "feature")

  return(output)
}

# Linear mixed-effect models without spline of time + joint test
fit_lmm_random_intercept_covars_JT <- function(
  x,
  y,
  key,
  visit,
  covariates_df
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
  reduced_model <- as.formula(paste(
    "y ~ visit",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))
  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  coef_table <- summary(fit)$tTable
  start_index <- grep("x", rownames(coef_table))[1]
  end_index <- grep("x:", rownames(coef_table))[1]
  output <- coef_table[c(start_index, end_index), c(1, 2, 4, 5)]

  fit_reduced <- tryCatch(
    {
      lme(reduced_model, random = ~ 1 | key, data = df, na.action = na.omit)
    },
    error = function(e) {
      NA
    }
  )

  fit_reduced <- tryCatch(
    {
      update(fit_reduced, . ~ ., method = "ML")
    },
    error = function(e) {
      NA
    }
  )

  if (length(fit_reduced) == 1) {
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
  } else {
    overall_sig <- anova(
      update(fit, . ~ ., method = "ML"),
      update(fit_reduced, . ~ ., method = "ML")
    )
    output <- cbind(output, c(overall_sig$`p-value`[2], rep(NA, time = 1)))
    output <- cbind(output, c(overall_sig$L.Ratio[2], rep(NA, time = 1)))
  }

  colnames(output) <- c(
    "coef",
    "se",
    "statistic",
    "pvalue",
    "pvalue_interaction",
    "L_ratio"
  )

  output <- as.data.frame(output)

  output$exposure <- names(x)

  return(output)
}

# Linear mixed-effect models without spline of time + joint test
# with confidence interval
fit_lmm_random_intercept_covars_JT_CI <- function(
  x,
  y,
  key,
  visit,
  covariates_df
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
  reduced_model <- as.formula(paste(
    "y ~ visit",
    "+",
    paste(
      colnames(covariates_df)[!colnames(covariates_df) == "key"],
      collapse = " + "
    )
  ))
  fit <- lme(model, random = ~ 1 | key, data = df, na.action = na.omit)
  coef_table <- bind_cols(
    summary(fit)$tTable,
    intervals(fit, which = "fixed")$fixed
  )
  coef_table$var_name <- rownames(summary(fit)$tTable)
  start_index <- grep("x", coef_table$var_name)[1]
  end_index <- grep("x:", coef_table$var_name)[1]
  output <- coef_table[c(start_index, end_index), c(1, 2, 4, 5, 6, 8)]

  fit_reduced <- tryCatch(
    {
      lme(reduced_model, random = ~ 1 | key, data = df, na.action = na.omit)
    },
    error = function(e) {
      NA
    }
  )

  fit_reduced <- tryCatch(
    {
      update(fit_reduced, . ~ ., method = "ML")
    },
    error = function(e) {
      NA
    }
  )

  if (length(fit_reduced) == 1) {
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
    output <- cbind(output, c(fit_reduced, rep(NA, time = 1)))
  } else {
    overall_sig <- anova(
      update(fit, . ~ ., method = "ML"),
      update(fit_reduced, . ~ ., method = "ML")
    )
    output <- cbind(output, c(overall_sig$`p-value`[2], rep(NA, time = 1)))
    output <- cbind(output, c(overall_sig$L.Ratio[2], rep(NA, time = 1)))
  }

  colnames(output) <- c(
    "coef",
    "se",
    "statistic",
    "pvalue",
    "conf_low",
    "conf_high",
    "pvalue_interaction",
    "L_ratio"
  )

  output <- as.data.frame(output)

  output$exposure <- names(x)

  output$feature <- names(y)

  output$beta_type <- c("main_effect", "interaction")

  return(output)
}
