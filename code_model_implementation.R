########################################################################################
# Masterarbeit: Klassifikation von Chromatografieläufen
########################################################################################

# Vorverarbeitung und Berechnung der Larson-Parameter
source("code_TA/code_data_loading_and_preprocessing_.R")
source("code_TA/code_transition_calculations_for_module_5.R")
source("code_TA/code_transition_calculations_for_module_6.R")

########################################################################################
# 1. Pakete laden
########################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(caret)
library(pROC)
library(xgboost)
library(precrec)
library(tibble)
library(Matrix)
library(broom)
library(lme4)
library(broom.mixed)
library(PRROC)
library(stringr)
library(forcats)

options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 13))

########################################################################################
# 2. Einstellungen
########################################################################################

USE_GROUPED_SPLIT <- FALSE     # TRUE = Split nach ResinID, FALSE = Standard-Split

SEED_MAIN      <- 777
K              <- 5
XGB_TUNE_N     <- 30
XGB_MAX_ROUNDS <- 500
XGB_EARLY_STOP <- 30
LARSON_ALPHA   <- 0.05

set.seed(SEED_MAIN)

split_tag <- ifelse(USE_GROUPED_SPLIT, "grouped_split", "standard_split")

results_dir <- "results"
plots_dir   <- "plots"

thesis_fig_dir   <- "thesis/iu-latex-template-main/figures/chapter5"
thesis_table_dir <- "thesis/iu-latex-template-main/tables/chapter5"

dir.create(results_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)
dir.create(thesis_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(thesis_table_dir, recursive = TRUE, showWarnings = FALSE)

########################################################################################
# 3. Hilfsfunktionen
########################################################################################

make_label <- function(x) {
  factor(x, levels = c("normal", "abnormal"))
}

label01 <- function(y) {
  as.integer(make_label(y) == "abnormal")
}

save_table <- function(df, filename) {
  readr::write_csv(df, file.path(results_dir, filename))
  readr::write_csv(df, file.path(thesis_table_dir, filename))
}

save_plot <- function(plot, filename, width = 8, height = 5, dpi = 300) {
  ggplot2::ggsave(file.path(plots_dir, filename), plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(file.path(thesis_fig_dir, filename), plot, width = width, height = height, dpi = dpi)
}

cm_counts <- function(scores, y, cutoff) {
  y <- make_label(y)
  pred <- ifelse(scores >= cutoff, "abnormal", "normal") |> make_label()
  tab <- table(pred, y)
  
  TP <- ifelse("abnormal" %in% rownames(tab) && "abnormal" %in% colnames(tab), tab["abnormal", "abnormal"], 0)
  TN <- ifelse("normal"   %in% rownames(tab) && "normal"   %in% colnames(tab), tab["normal", "normal"], 0)
  FP <- ifelse("abnormal" %in% rownames(tab) && "normal"   %in% colnames(tab), tab["abnormal", "normal"], 0)
  FN <- ifelse("normal"   %in% rownames(tab) && "abnormal" %in% colnames(tab), tab["normal", "abnormal"], 0)
  
  list(TP = as.numeric(TP), TN = as.numeric(TN), FP = as.numeric(FP), FN = as.numeric(FN))
}

metrics_from_counts <- function(TP, TN, FP, FN) {
  acc  <- (TP + TN) / max(1, TP + TN + FP + FN)
  sens <- TP / max(1, TP + FN)
  spec <- TN / max(1, TN + FP)
  prec <- TP / max(1, TP + FP)
  npv  <- TN / max(1, TN + FN)
  ba   <- (sens + spec) / 2
  
  denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
  mcc <- ifelse(denom > 0, ((TP * TN) - (FP * FN)) / sqrt(denom), 0)
  
  tibble(
    Accuracy = acc,
    Sensitivity = sens,
    Specificity = spec,
    BalancedAccuracy = ba,
    Precision = prec,
    NPV = npv,
    MCC = mcc
  )
}

auc_roc_pr <- function(y, scores) {
  out <- tryCatch({
    y01 <- label01(y)
    mm  <- precrec::evalmod(scores = scores, labels = y01)
    aucs <- precrec::auc(mm)
    
    list(
      ROC_AUC = aucs %>% filter(curvetypes == "ROC") %>% pull(aucs) %>% as.numeric(),
      PR_AUC  = aucs %>% filter(curvetypes == "PRC") %>% pull(aucs) %>% as.numeric()
    )
  }, error = function(e) {
    list(ROC_AUC = NA_real_, PR_AUC = NA_real_)
  })
  out
}

threshold_grid <- function(scores) {
  s <- sort(unique(as.numeric(scores)))
  s <- s[is.finite(s)]
  
  if (length(s) == 0) return(numeric(0))
  if (length(s) == 1) return(s)
  
  mids <- (s[-1] + s[-length(s)]) / 2
  c(s[1] - 1e-12, mids, s[length(s)] + 1e-12)
}

choose_cutoff <- function(scores, y, rule = c("fixed_0_5", "youden")) {
  rule <- match.arg(rule)
  y <- make_label(y)
  
  if (rule == "fixed_0_5") {
    return(list(cutoff = 0.5, rule = "fixed_0_5"))
  }
  
  thr <- threshold_grid(scores)
  if (length(thr) == 0) {
    return(list(cutoff = NA_real_, rule = rule))
  }
  
  eval_df <- bind_rows(lapply(thr, function(t) {
    cc <- cm_counts(scores, y, t)
    metrics_from_counts(cc$TP, cc$TN, cc$FP, cc$FN) %>%
      mutate(threshold = t, YoudenJ = Sensitivity + Specificity - 1)
  }))
  
  best <- eval_df %>%
    arrange(desc(YoudenJ), desc(MCC), desc(BalancedAccuracy)) %>%
    slice(1)
  
  list(cutoff = best$threshold, rule = "youden", summary = best)
}

pick_cutoffs <- function(scores, y) {
  list(
    fixed_0_5 = choose_cutoff(scores, y, "fixed_0_5"),
    youden    = choose_cutoff(scores, y, "youden")
  )
}

evaluate_model <- function(method_name, scores_test, y_test, aucs, cutoff_list) {
  bind_rows(lapply(names(cutoff_list), function(rule_name) {
    cutoff <- as.numeric(cutoff_list[[rule_name]]$cutoff)
    cc <- cm_counts(scores_test, y_test, cutoff)
    mets <- metrics_from_counts(cc$TP, cc$TN, cc$FP, cc$FN)
    
    tibble(
      Methode = method_name,
      Cutoff_Art = rule_name,
      Cutoff_Wert = cutoff
    ) %>%
      bind_cols(mets) %>%
      mutate(
        ROC_AUC = as.numeric(aucs$ROC_AUC),
        PR_AUC = as.numeric(aucs$PR_AUC)
      )
  }))
}

########################################################################################
# 4. Larson-Funktionen
########################################################################################

larson_weights <- c(
  NoIP = 0.6,
  MRoC = 0.7,
  BV = 0.7,
  CE = 0.8,
  A = 0.1,
  NGHETP = 0.8,
  GHETP = 0.1
)

larson_metrics <- c("NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP")

compute_limits_quantile <- function(normal_df, metric_names, alpha = 0.05) {
  out <- list()
  
  for (m in metric_names) {
    x <- normal_df[[m]]
    
    if (m %in% c("MRoC", "BV")) {
      out[[m]] <- list(type = "lower", lower = as.numeric(quantile(x, alpha, na.rm = TRUE)))
    } else if (m == "A") {
      out[[m]] <- list(
        type = "two_sided",
        lower = as.numeric(quantile(x, alpha, na.rm = TRUE)),
        upper = as.numeric(quantile(x, 1 - alpha, na.rm = TRUE))
      )
    } else {
      out[[m]] <- list(type = "upper", upper = as.numeric(quantile(x, 1 - alpha, na.rm = TRUE)))
    }
  }
  
  out
}

compute_oiv_quantile <- function(metrics, limits, weights) {
  pass <- sapply(names(metrics), function(m) {
    val <- metrics[[m]]
    lim <- limits[[m]]
    
    if (!is.finite(val)) return(0)
    
    if (lim$type == "lower") {
      if (val < lim$lower) {
        ifelse(val < lim$lower * 0.9, 0, 0.5)
      } else {
        1
      }
    } else if (lim$type == "upper") {
      if (val > lim$upper) {
        ifelse(val > lim$upper * 1.1, 0, 0.5)
      } else {
        1
      }
    } else {
      if (val >= lim$lower && val <= lim$upper) {
        1
      } else if (val < lim$lower) {
        ifelse(val < lim$lower * 0.9, 0, 0.5)
      } else {
        ifelse(val > lim$upper * 1.1, 0, 0.5)
      }
    }
  })
  
  sum(pass * weights[names(metrics)]) / sum(weights)
}

compute_limits_3sigma <- function(normal_df, metric_names) {
  out <- lapply(metric_names, function(m) {
    x <- normal_df[[m]]
    mu <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    list(KIL = mu - 3 * sd, KIU = mu + 3 * sd)
  })
  names(out) <- metric_names
  out
}

compute_oiv_3sigma <- function(metrics, limits, weights) {
  pass <- sapply(names(metrics), function(m) {
    val <- metrics[[m]]
    KIL <- limits[[m]]$KIL
    KIU <- limits[[m]]$KIU
    
    if (!is.finite(val) || !is.finite(KIL) || !is.finite(KIU)) return(0)
    if (val >= KIL && val <= KIU) return(1)
    if (val < 0.9 * KIL || val > 1.1 * KIU) return(0)
    0.5
  })
  
  sum(pass * weights[names(metrics)]) / sum(weights)
}

get_oof_larson <- function(train_df, folds, metric_names, weights, method = "quantile", alpha = 0.05) {
  oof_scores <- rep(NA_real_, nrow(train_df))
  
  for (f in seq_along(folds)) {
    valid_idx <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(train_df)), valid_idx)
    
    train_fold <- train_df[train_idx, ]
    valid_fold <- train_df[valid_idx, ]
    normal_fold <- train_fold %>% filter(Label == "normal")
    
    if (method == "quantile") {
      limits <- compute_limits_quantile(normal_fold, metric_names, alpha)
      raw_oiv <- apply(valid_fold[, metric_names, drop = FALSE], 1, function(r) {
        compute_oiv_quantile(as.list(r), limits, weights)
      })
    } else {
      limits <- compute_limits_3sigma(normal_fold, metric_names)
      raw_oiv <- apply(valid_fold[, metric_names, drop = FALSE], 1, function(r) {
        compute_oiv_3sigma(as.list(r), limits, weights)
      })
    }
    
    oof_scores[valid_idx] <- 1 - raw_oiv
  }
  
  oof_scores
}

get_test_larson <- function(train_df, test_df, metric_names, weights, method = "quantile", alpha = 0.05) {
  normal_train <- train_df %>% filter(Label == "normal")
  
  if (method == "quantile") {
    limits <- compute_limits_quantile(normal_train, metric_names, alpha)
    raw_oiv <- apply(test_df[, metric_names, drop = FALSE], 1, function(r) {
      compute_oiv_quantile(as.list(r), limits, weights)
    })
  } else {
    limits <- compute_limits_3sigma(normal_train, metric_names)
    raw_oiv <- apply(test_df[, metric_names, drop = FALSE], 1, function(r) {
      compute_oiv_3sigma(as.list(r), limits, weights)
    })
  }
  
  1 - raw_oiv
}

########################################################################################
# 5. Logistische Regression
########################################################################################

make_glm_formula <- function(features, outcome = "Label_num") {
  as.formula(paste(outcome, "~", paste(features, collapse = " + ")))
}

get_logit_class_weights <- function(y) {
  y <- make_label(y)
  weight_abnormal <- sum(y == "normal") / max(1, sum(y == "abnormal"))
  ifelse(y == "abnormal", weight_abnormal, 1)
}

get_oof_logit <- function(train_df, folds, features) {
  oof_scores <- rep(NA_real_, nrow(train_df))
  form <- make_glm_formula(features)
  
  for (f in seq_along(folds)) {
    valid_idx <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(train_df)), valid_idx)
    
    train_fold <- train_df[train_idx, ]
    valid_fold <- train_df[valid_idx, ]
    train_fold$class_weight <- get_logit_class_weights(train_fold$Label)
    
    model <- suppressWarnings(glm(
      formula = form,
      data = train_fold,
      family = binomial(link = "logit"),
      weights = class_weight
    ))
    
    oof_scores[valid_idx] <- predict(model, newdata = valid_fold, type = "response")
  }
  
  oof_scores
}

fit_logit <- function(train_df, features) {
  form <- make_glm_formula(features)
  train_df$class_weight <- get_logit_class_weights(train_df$Label)
  
  suppressWarnings(glm(
    formula = form,
    data = train_df,
    family = binomial(link = "logit"),
    weights = class_weight
  ))
}

run_logit <- function(train_df, test_df, folds, features, method_name) {
  oof_scores <- get_oof_logit(train_df, folds, features)
  cutoffs <- pick_cutoffs(oof_scores, train_df$Label)
  
  final_model <- fit_logit(train_df, features)
  test_scores <- predict(final_model, newdata = test_df, type = "response")
  aucs <- auc_roc_pr(test_df$Label, test_scores)
  
  list(
    method_name = method_name,
    features = features,
    model = final_model,
    oof_scores = oof_scores,
    test_scores = test_scores,
    cutoffs = cutoffs,
    aucs = aucs
  )
}

########################################################################################
# 6. XGBoost
########################################################################################

make_xgb_grid <- function(n = 30) {
  tibble(
    eta = sample(c(0.01, 0.03, 0.05, 0.08, 0.10, 0.15), n, replace = TRUE),
    max_depth = sample(2:6, n, replace = TRUE),
    min_child_weight = sample(c(1, 2, 5, 10), n, replace = TRUE),
    subsample = sample(seq(0.6, 1.0, by = 0.1), n, replace = TRUE),
    colsample_bytree = sample(seq(0.6, 1.0, by = 0.1), n, replace = TRUE),
    gamma = sample(c(0, 0.5, 1), n, replace = TRUE),
    lambda = sample(c(0, 0.5, 1, 2, 5), n, replace = TRUE)
  ) %>% distinct()
}

tune_xgb <- function(X_train, y_train, folds, n_iter = 30) {
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  scale_pos_weight <- sum(y_train == 0) / max(1, sum(y_train == 1))
  grid <- make_xgb_grid(n_iter)
  
  tuning_results <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    
    params <- list(
      objective = "binary:logistic",
      eval_metric = "aucpr",
      scale_pos_weight = scale_pos_weight,
      eta = g$eta,
      max_depth = g$max_depth,
      min_child_weight = g$min_child_weight,
      subsample = g$subsample,
      colsample_bytree = g$colsample_bytree,
      gamma = g$gamma,
      lambda = g$lambda
    )
    
    cv <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = XGB_MAX_ROUNDS,
      folds = folds,
      early_stopping_rounds = XGB_EARLY_STOP,
      maximize = TRUE,
      verbose = 0
    )
    
    best_iter <- cv$best_iteration
    score_col <- names(cv$evaluation_log)[grepl("test_aucpr_mean", names(cv$evaluation_log))][1]
    best_score <- cv$evaluation_log[[score_col]][best_iter]
    
    tuning_results[[i]] <- bind_cols(
      g,
      tibble(iter = i, best_nrounds = best_iter, best_cv_prauc = best_score)
    )
  }
  
  tuning_table <- bind_rows(tuning_results) %>% arrange(desc(best_cv_prauc))
  best <- tuning_table %>% slice(1)
  
  best_params <- list(
    objective = "binary:logistic",
    eval_metric = "aucpr",
    scale_pos_weight = scale_pos_weight,
    eta = best$eta,
    max_depth = best$max_depth,
    min_child_weight = best$min_child_weight,
    subsample = best$subsample,
    colsample_bytree = best$colsample_bytree,
    gamma = best$gamma,
    lambda = best$lambda
  )
  
  list(
    tuning_table = tuning_table,
    best_params = best_params,
    best_nrounds = as.integer(best$best_nrounds)
  )
}

get_oof_xgb <- function(train_df, folds, features, params, nrounds) {
  mm_terms <- terms(~ . - 1, data = train_df[, features, drop = FALSE])
  oof_scores <- rep(NA_real_, nrow(train_df))
  
  for (f in seq_along(folds)) {
    valid_idx <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(train_df)), valid_idx)
    
    train_fold <- train_df[train_idx, ]
    valid_fold <- train_df[valid_idx, ]
    
    X_train <- sparse.model.matrix(mm_terms, data = train_fold[, features, drop = FALSE])
    X_valid <- sparse.model.matrix(mm_terms, data = valid_fold[, features, drop = FALSE])
    
    model <- xgb.train(
      params = params,
      data = xgb.DMatrix(X_train, label = train_fold$Label_num),
      nrounds = nrounds,
      verbose = 0
    )
    
    oof_scores[valid_idx] <- predict(model, xgb.DMatrix(X_valid))
  }
  
  oof_scores
}

fit_xgb <- function(train_df, test_df, features, params, nrounds) {
  mm_terms <- terms(~ . - 1, data = train_df[, features, drop = FALSE])
  X_train <- sparse.model.matrix(mm_terms, data = train_df[, features, drop = FALSE])
  X_test  <- sparse.model.matrix(mm_terms, data = test_df[, features, drop = FALSE])
  
  model <- xgb.train(
    params = params,
    data = xgb.DMatrix(X_train, label = train_df$Label_num),
    nrounds = nrounds,
    verbose = 0
  )
  
  list(
    model = model,
    test_scores = predict(model, xgb.DMatrix(X_test)),
    X_train = X_train
  )
}

run_xgb <- function(train_df, test_df, folds, features, method_name) {
  mm_terms <- terms(~ . - 1, data = train_df[, features, drop = FALSE])
  X_train <- sparse.model.matrix(mm_terms, data = train_df[, features, drop = FALSE])
  y_train <- train_df$Label_num
  
  tuned <- tune_xgb(X_train, y_train, folds, XGB_TUNE_N)
  oof_scores <- get_oof_xgb(train_df, folds, features, tuned$best_params, tuned$best_nrounds)
  cutoffs <- pick_cutoffs(oof_scores, train_df$Label)
  final <- fit_xgb(train_df, test_df, features, tuned$best_params, tuned$best_nrounds)
  aucs <- auc_roc_pr(test_df$Label, final$test_scores)
  
  list(
    method_name = method_name,
    features = features,
    tuning = tuned,
    model = final$model,
    X_train = final$X_train,
    oof_scores = oof_scores,
    test_scores = final$test_scores,
    cutoffs = cutoffs,
    aucs = aucs
  )
}

########################################################################################
# 7. Daten laden und vorbereiten
########################################################################################

df_m5 <- read.csv("data_for_modeling/m5_postwash_larson.csv", sep = ",", header = TRUE)
df_m6 <- read.csv("data_for_modeling/m6_postwash_larson.csv", sep = ",", header = TRUE)

df <- rbind(df_m5, df_m6) %>%
  rename(Module = FacilityID)

df$ResinID <- sub("^[^|]*\\|", "", df$ResinID)
df$Label <- make_label(df$Label)
df$Module <- as.factor(df$Module)
df$Product_Concentration <- as.factor(df$Product_Concentration)
df$ResinID <- as.factor(df$ResinID)
df$Label_num <- label01(df$Label)

cat("\nGesamte Klassenverteilung:\n")
print(table(df$Label))

basic_features <- c("NoIP", "MRoC", "BV", "CE", "A", "NGHETP", "GHETP")

extended_features <- c(
  basic_features,
  "Resin_Cycle_No", "start_pH", "min_flow",
  "Product_Concentration", "Module"
)

needed_vars <- unique(c("Label", "Label_num", "ResinID", basic_features, extended_features))

df_model <- df %>%
  filter(complete.cases(across(all_of(needed_vars))))

cat("\nKlassenverteilung nach Entfernen fehlender Werte:\n")
print(table(df_model$Label))

########################################################################################
# 8. Train-Test-Split und Kreuzvalidierung
########################################################################################

make_grouped_train_test_split <- function(df, group_var = "ResinID", p_train = 0.70, seed = 777) {
  set.seed(seed)
  groups <- unique(df[[group_var]])
  train_groups <- sample(groups, size = floor(p_train * length(groups)), replace = FALSE)
  
  list(
    train_idx = which(df[[group_var]] %in% train_groups),
    test_idx = which(!df[[group_var]] %in% train_groups)
  )
}

make_grouped_folds <- function(df, group_var = "ResinID", k = 5, seed = 777) {
  set.seed(seed)
  groups <- sample(unique(df[[group_var]]))
  fold_ids <- cut(seq_along(groups), breaks = k, labels = FALSE)
  group_folds <- tibble(group = groups, fold = fold_ids)
  
  lapply(seq_len(k), function(i) {
    valid_groups <- group_folds %>% filter(fold == i) %>% pull(group)
    which(df[[group_var]] %in% valid_groups)
  })
}

if (USE_GROUPED_SPLIT) {
  split <- make_grouped_train_test_split(df_model, "ResinID", 0.70, SEED_MAIN)
  train_df <- df_model[split$train_idx, ]
  test_df  <- df_model[split$test_idx, ]
} else {
  idx_train <- caret::createDataPartition(df_model$Label, p = 0.70, list = FALSE)
  train_df <- df_model[idx_train, ]
  test_df  <- df_model[-idx_train, ]
}

train_df$Label <- make_label(train_df$Label)
test_df$Label  <- make_label(test_df$Label)

if (USE_GROUPED_SPLIT) {
  folds <- make_grouped_folds(train_df, "ResinID", K, SEED_MAIN)
} else {
  folds <- caret::createFolds(train_df$Label, k = K, returnTrain = FALSE)
}

split_info <- tibble(
  Split = ifelse(USE_GROUPED_SPLIT, "Grouped split by ResinID", "Standard run-level split"),
  Train_N = nrow(train_df),
  Test_N = nrow(test_df),
  Train_Abnormal = sum(train_df$Label == "abnormal"),
  Test_Abnormal = sum(test_df$Label == "abnormal"),
  Train_Resins = n_distinct(train_df$ResinID),
  Test_Resins = n_distinct(test_df$ResinID)
)

save_table(split_info, paste0("split_info_", split_tag, ".csv"))
print(split_info)

########################################################################################
# 9. Modelle berechnen
########################################################################################

# Larson Quantil
larson_quant_oof <- get_oof_larson(train_df, folds, larson_metrics, larson_weights, method = "quantile", alpha = LARSON_ALPHA)
larson_quant_cutoffs <- pick_cutoffs(larson_quant_oof, train_df$Label)
larson_quant_test <- get_test_larson(train_df, test_df, larson_metrics, larson_weights, method = "quantile", alpha = LARSON_ALPHA)
larson_quant_auc <- auc_roc_pr(test_df$Label, larson_quant_test)

# Larson 3-Sigma
larson_3sigma_oof <- get_oof_larson(train_df, folds, larson_metrics, larson_weights, method = "3sigma")
larson_3sigma_cutoffs <- pick_cutoffs(larson_3sigma_oof, train_df$Label)
larson_3sigma_test <- get_test_larson(train_df, test_df, larson_metrics, larson_weights, method = "3sigma")
larson_3sigma_auc <- auc_roc_pr(test_df$Label, larson_3sigma_test)

# Logistische Regression
logit_basic <- run_logit(train_df, test_df, folds, basic_features, "Logit_basic")
logit_extended <- run_logit(train_df, test_df, folds, extended_features, "Logit_extended")

# XGBoost
xgb_basic <- run_xgb(train_df, test_df, folds, basic_features, "XGBoost_basic")
xgb_extended <- run_xgb(train_df, test_df, folds, extended_features, "XGBoost_extended")

save_table(xgb_basic$tuning$tuning_table, paste0("xgb_tuning_basic_", split_tag, ".csv"))
save_table(xgb_extended$tuning$tuning_table, paste0("xgb_tuning_extended_", split_tag, ".csv"))

########################################################################################
# 10. Ergebnisse auswerten
########################################################################################

results_test <- bind_rows(
  evaluate_model("Larson_OIV_quantile", larson_quant_test, test_df$Label, larson_quant_auc, larson_quant_cutoffs),
  evaluate_model("Larson_OIV_3sigma", larson_3sigma_test, test_df$Label, larson_3sigma_auc, larson_3sigma_cutoffs),
  evaluate_model("Logit_basic", logit_basic$test_scores, test_df$Label, logit_basic$aucs, logit_basic$cutoffs),
  evaluate_model("Logit_extended", logit_extended$test_scores, test_df$Label, logit_extended$aucs, logit_extended$cutoffs),
  evaluate_model("XGBoost_basic", xgb_basic$test_scores, test_df$Label, xgb_basic$aucs, xgb_basic$cutoffs),
  evaluate_model("XGBoost_extended", xgb_extended$test_scores, test_df$Label, xgb_extended$aucs, xgb_extended$cutoffs)
) %>%
  mutate(across(
    c(Cutoff_Wert, Accuracy, Sensitivity, Specificity, BalancedAccuracy, Precision, NPV, MCC, ROC_AUC, PR_AUC),
    ~ round(.x, 3)
  )) %>%
  arrange(Methode, Cutoff_Art)

print(results_test, n = Inf)
save_table(results_test, paste0("results_test_main_models_", split_tag, ".csv"))

cutoff_table <- bind_rows(
  tibble(Methode = "Larson_OIV_quantile", Cutoff_Art = names(larson_quant_cutoffs), Cutoff_Wert = sapply(larson_quant_cutoffs, function(x) x$cutoff)),
  tibble(Methode = "Larson_OIV_3sigma", Cutoff_Art = names(larson_3sigma_cutoffs), Cutoff_Wert = sapply(larson_3sigma_cutoffs, function(x) x$cutoff)),
  tibble(Methode = "Logit_basic", Cutoff_Art = names(logit_basic$cutoffs), Cutoff_Wert = sapply(logit_basic$cutoffs, function(x) x$cutoff)),
  tibble(Methode = "Logit_extended", Cutoff_Art = names(logit_extended$cutoffs), Cutoff_Wert = sapply(logit_extended$cutoffs, function(x) x$cutoff)),
  tibble(Methode = "XGBoost_basic", Cutoff_Art = names(xgb_basic$cutoffs), Cutoff_Wert = sapply(xgb_basic$cutoffs, function(x) x$cutoff)),
  tibble(Methode = "XGBoost_extended", Cutoff_Art = names(xgb_extended$cutoffs), Cutoff_Wert = sapply(xgb_extended$cutoffs, function(x) x$cutoff))
)

save_table(cutoff_table, paste0("chosen_cutoffs_from_oof_", split_tag, ".csv"))

########################################################################################
# 11. Zentrale Grafiken
########################################################################################

pretty_method_name <- function(x) {
  recode(
    x,
    "Larson_OIV_quantile" = "Larson (Quantil)",
    "Larson_OIV_3sigma" = "Larson (3σ)",
    "Logit_basic" = "Logit (Basis)",
    "Logit_extended" = "Logit (Erweitert)",
    "XGBoost_basic" = "XGBoost (Basis)",
    "XGBoost_extended" = "XGBoost (Erweitert)"
  )
}

pretty_cutoff_name <- function(x) {
  recode(x, "fixed_0_5" = "Schwelle 0,5", "youden" = "Youden")
}

results_plot <- results_test %>%
  select(Methode, Cutoff_Art, BalancedAccuracy, MCC) %>%
  pivot_longer(cols = c(BalancedAccuracy, MCC), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Methode = pretty_method_name(Methode),
    Cutoff_Art = pretty_cutoff_name(Cutoff_Art),
    Metric = recode(Metric, "BalancedAccuracy" = "Balanced Accuracy", "MCC" = "Matthews Correlation Coefficient")
  )

p_results <- ggplot(results_plot, aes(x = Methode, y = Value, fill = Cutoff_Art)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(
    title = "Vergleich der Modelle nach Schwellenwahl",
    x = NULL,
    y = "Metrikwert",
    fill = NULL
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(p_results, paste0("model_metric_comparison_", split_tag, ".png"), width = 11, height = 6.5)

# ROC-Kurven
roc_df <- bind_rows(
  tibble(Methode = "Larson (Quantil)", Score = larson_quant_test),
  tibble(Methode = "Larson (3σ)", Score = larson_3sigma_test),
  tibble(Methode = "Logit (Basis)", Score = logit_basic$test_scores),
  tibble(Methode = "Logit (Erweitert)", Score = logit_extended$test_scores),
  tibble(Methode = "XGBoost (Basis)", Score = xgb_basic$test_scores),
  tibble(Methode = "XGBoost (Erweitert)", Score = xgb_extended$test_scores)
)

roc_plot_df <- roc_df %>%
  group_by(Methode) %>%
  group_modify(~ {
    roc_obj <- pROC::roc(test_df$Label, .x$Score, levels = c("normal", "abnormal"), direction = "<", quiet = TRUE)
    tibble(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
  }) %>%
  ungroup()

p_roc <- ggplot(roc_plot_df, aes(x = fpr, y = tpr, color = Methode)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC-Kurven der Modelle", x = "1 - Spezifität", y = "Sensitivität", color = NULL)

save_plot(p_roc, paste0("roc_test_main_models_", split_tag, ".png"), width = 9, height = 5.5)

# Precision-Recall-Kurven
make_pr_df <- function(scores, y, method_name) {
  y01 <- label01(y)
  pr <- PRROC::pr.curve(
    scores.class0 = scores[y01 == 1],
    scores.class1 = scores[y01 == 0],
    curve = TRUE
  )
  
  tibble(Methode = method_name, Recall = pr$curve[, 1], Precision = pr$curve[, 2])
}

pr_df <- bind_rows(
  make_pr_df(larson_quant_test, test_df$Label, "Larson (Quantil)"),
  make_pr_df(larson_3sigma_test, test_df$Label, "Larson (3σ)"),
  make_pr_df(logit_basic$test_scores, test_df$Label, "Logit (Basis)"),
  make_pr_df(logit_extended$test_scores, test_df$Label, "Logit (Erweitert)"),
  make_pr_df(xgb_basic$test_scores, test_df$Label, "XGBoost (Basis)"),
  make_pr_df(xgb_extended$test_scores, test_df$Label, "XGBoost (Erweitert)")
)

baseline_pr <- mean(test_df$Label == "abnormal")

p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision, color = Methode)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = baseline_pr, linetype = "dashed") +
  labs(
    title = "Precision-Recall-Kurven der Modelle",
    subtitle = paste0("PR-Baseline = ", round(baseline_pr, 3)),
    x = "Recall",
    y = "Precision",
    color = NULL
  )

save_plot(p_pr, paste0("pr_test_main_models_", split_tag, ".png"), width = 9, height = 5.5)

########################################################################################
# 12. Odds Ratios und XGBoost-Feature-Importance
########################################################################################

make_or_table <- function(model, data_for_sd, feature_names) {
  tab <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  numeric_features <- feature_names[sapply(data_for_sd[, feature_names, drop = FALSE], is.numeric)]
  sd_values <- sapply(data_for_sd[, numeric_features, drop = FALSE], sd, na.rm = TRUE)
  sd_df <- tibble(term = names(sd_values), SD = as.numeric(sd_values))
  
  tab %>%
    filter(term != "(Intercept)") %>%
    left_join(sd_df, by = "term") %>%
    mutate(
      SD = ifelse(is.na(SD), 1, SD),
      OR_1SD = exp(estimate * SD),
      OR_1SD_lo = exp(conf.low * SD),
      OR_1SD_hi = exp(conf.high * SD)
    ) %>%
    arrange(desc(abs(log(OR_1SD))))
}

logit_basic_or <- make_or_table(logit_basic$model, train_df, basic_features)
logit_extended_or <- make_or_table(logit_extended$model, train_df, extended_features)

save_table(logit_basic_or, paste0("logit_basic_odds_ratios_", split_tag, ".csv"))
save_table(logit_extended_or, paste0("logit_extended_odds_ratios_", split_tag, ".csv"))

xgb_imp_basic <- xgb.importance(model = xgb_basic$model, feature_names = colnames(xgb_basic$X_train)) %>%
  as_tibble() %>%
  arrange(desc(Gain))

xgb_imp_extended <- xgb.importance(model = xgb_extended$model, feature_names = colnames(xgb_extended$X_train)) %>%
  as_tibble() %>%
  arrange(desc(Gain))

save_table(xgb_imp_basic, paste0("xgb_importance_basic_", split_tag, ".csv"))
save_table(xgb_imp_extended, paste0("xgb_importance_extended_", split_tag, ".csv"))

p_xgb_imp <- xgb_imp_extended %>%
  mutate(Feature = reorder(Feature, Gain)) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col() +
  coord_flip() +
  labs(title = "Variablenwichtigkeit von XGBoost (Erweitert)", x = NULL, y = "Gain")

save_plot(p_xgb_imp, paste0("xgb_feature_importance_extended_", split_tag, ".png"), width = 8, height = 6)

########################################################################################
# 13. Gemischtes logistisches Modell
########################################################################################

numeric_features <- extended_features[sapply(df_model[, extended_features, drop = FALSE], is.numeric)]
df_mixed <- df_model

for (v in numeric_features) {
  df_mixed[[paste0(v, "_z")]] <- as.numeric(scale(df_mixed[[v]]))
}

mixed_features <- extended_features
mixed_features[mixed_features %in% numeric_features] <- paste0(numeric_features, "_z")

mixed_formula <- as.formula(
  paste("Label_num ~", paste(mixed_features, collapse = " + "), "+ (1 | ResinID)")
)

mixed_model <- glmer(
  formula = mixed_formula,
  data = df_mixed,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

mixed_or <- broom.mixed::tidy(mixed_model, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    OR_lo = exp(conf.low),
    OR_hi = exp(conf.high)
  )

save_table(mixed_or, "mixed_logit_odds_ratios.csv")

random_effects <- ranef(mixed_model)$ResinID %>%
  rownames_to_column("ResinID") %>%
  rename(RandomIntercept = `(Intercept)`) %>%
  arrange(RandomIntercept)

save_table(random_effects, "mixed_model_random_intercepts_resinid.csv")

p_random <- ggplot(random_effects, aes(x = reorder(ResinID, RandomIntercept), y = RandomIntercept)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  coord_flip() +
  labs(
    title = "Zufällige Interzepte nach Resin-Charge",
    x = "Resin-Charge",
    y = "Zufälliger Interzept"
  )

save_plot(p_random, "mixed_model_random_intercepts_resinid.png", width = 8, height = 6)

var_resin <- as.data.frame(VarCorr(mixed_model)) %>%
  filter(grp == "ResinID", var1 == "(Intercept)") %>%
  pull(vcov)

icc_table <- tibble(
  RandomEffect = "ResinID",
  Variance_RE = var_resin,
  Variance_Logistic = pi^2 / 3,
  ICC = var_resin / (var_resin + pi^2 / 3)
)

save_table(icc_table, "mixed_model_icc.csv")


