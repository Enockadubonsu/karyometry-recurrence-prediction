
library(glmnet)
library(randomForest)
library(MASS)
library(nnet)

set.seed(42)  

bdat<- readRDS("C:/Users/adubo/Desktop/Folders/PhD BIOSTATS/BIOS 648/Project/bladder.rds")

feature_cols <- paste0("f", 1:92)
cat("Dimensions:", dim(bdat), "\n")
cat("Columns:", names(bdat)[93:95], "\n")   
cat("Patients:", length(unique(bdat$id)), "\n")
print(table(bdat$group))                      

patient_ids    <- unique(bdat$id)
patient_labels <- sapply(patient_ids, function(pid)
  bdat$group[bdat$id == pid][1])
names(patient_labels) <- patient_ids

n_pat  <- length(patient_ids)   # 84
y_pat  <- ifelse(patient_labels == "R", 1, 0)  # 1=recurrent, 0=NR
cat("\nPatient-level:\n")
print(table(patient_labels))    # 39 R, 45 NR




standardise <- function(train_mat, test_mat) {
  mu  <- colMeans(train_mat)
  sig <- apply(train_mat, 2, sd)
  sig[sig == 0] <- 1            # guard zero-variance features
  train_s <- sweep(sweep(train_mat, 2, mu, "-"), 2, sig, "/")
  test_s  <- sweep(sweep(test_mat,  2, mu, "-"), 2, sig, "/")
  list(train = train_s, test = test_s, mu = mu, sig = sig)
}

# ── 3. REPRESENTATION A — mean + quantile summaries (460-dim) ───
make_rep_A <- function(bdat, patient_ids, feature_cols) {
  q_levels <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  mat <- t(sapply(patient_ids, function(pid) {
    sub <- as.matrix(bdat[bdat$id == pid, feature_cols])
    c(colMeans(sub),
      apply(sub, 2, sd),
      apply(sub, 2, quantile, probs = q_levels[1]),
      apply(sub, 2, quantile, probs = q_levels[2]),
      apply(sub, 2, quantile, probs = q_levels[3]),
      apply(sub, 2, quantile, probs = q_levels[4]),
      apply(sub, 2, quantile, probs = q_levels[5]))
  }))
  rownames(mat) <- patient_ids
  # 92 * 7 stats = 644 ... actually we want 5 stats per feature = 460
  # mean + sd + Q10 + Q50 + Q90 = 5 * 92 = 460
  mat2 <- t(sapply(patient_ids, function(pid) {
    sub <- as.matrix(bdat[bdat$id == pid, feature_cols])
    c(colMeans(sub),
      apply(sub, 2, sd),
      apply(sub, 2, quantile, probs = 0.10, type = 7),
      apply(sub, 2, quantile, probs = 0.50, type = 7),
      apply(sub, 2, quantile, probs = 0.90, type = 7))
  }))
  rownames(mat2) <- patient_ids
  colnames(mat2) <- c(paste0(feature_cols, "_mean"),
                      paste0(feature_cols, "_sd"),
                      paste0(feature_cols, "_q10"),
                      paste0(feature_cols, "_q50"),
                      paste0(feature_cols, "_q90"))
  mat2
}

cat("\nBuilding Representation A...\n")
Rep_A <- make_rep_A(bdat, patient_ids, feature_cols)
cat("Rep A dimensions:", dim(Rep_A), "\n")   # 84 x 460


# ── 4. Nucleus scorer helper (used in Reps B and C) ─────────────
# Trains RF on training nuclei, returns predicted probs for ALL patients
# (called inside each LOOCV fold)
train_nucleus_scorer <- function(train_nuclei, train_labels_nuc,
                                 all_nuclei_std, feature_cols) {
  rf_scorer <- randomForest(
    x         = train_nuclei[, feature_cols],
    y         = as.factor(train_labels_nuc),
    ntree     = 300,
    importance = FALSE
  )
  # Predict on ALL nuclei (training + test patient)
  predict(rf_scorer, all_nuclei_std[, feature_cols], type = "prob")[, "R"]
}

# ── 5. Softmax weighting ─────────────────────────────────────────
softmax_weights <- function(scores, tau) {
  if (is.infinite(tau)) return(rep(1 / length(scores), length(scores)))
  w <- exp(scores / tau)
  w / sum(w)
}

# ── 6. Build Reps B and C given nucleus scores ───────────────────
make_rep_B <- function(bdat, patient_ids, feature_cols, nucleus_scores, tau) {
  mat <- t(sapply(patient_ids, function(pid) {
    idx    <- which(bdat$id == pid)
    scores <- nucleus_scores[idx]
    w      <- softmax_weights(scores, tau)
    feat   <- as.matrix(bdat[idx, feature_cols])
    colSums(w * feat)
  }))
  rownames(mat) <- patient_ids
  mat
}

make_rep_C <- function(bdat, patient_ids, nucleus_scores) {
  mat <- t(sapply(patient_ids, function(pid) {
    idx    <- which(bdat$id == pid)
    scores <- nucleus_scores[idx]
    c(mean  = mean(scores),
      max   = max(scores),
      q90   = quantile(scores, 0.90, type = 7))
  }))
  rownames(mat) <- patient_ids
  mat
}


# ── 7. Four classifiers ──────────────────────────────────────────

# 7a. Lasso logistic regression
fit_lasso <- function(X_train, y_train, X_test) {
  cv_fit <- cv.glmnet(X_train, y_train, family = "binomial",
                      alpha = 1, nfolds = 5, type.measure = "class")
  prob <- predict(cv_fit, X_test, s = "lambda.min", type = "response")
  as.numeric(prob > 0.5)
}

# 7b. Random forest
fit_rf <- function(X_train, y_train, X_test) {
  rf <- randomForest(x = X_train, y = as.factor(y_train), ntree = 500)
  as.numeric(predict(rf, X_test) == "1")
}

# 7c. PCA + LDA (for Rep A); plain LDA for Reps B and C
fit_pca_lda <- function(X_train, y_train, X_test, use_pca = TRUE) {
  
  # ── Remove near-zero-variance columns (LDA will crash on these) ──
  col_sd <- apply(X_train, 2, sd)
  keep   <- col_sd > 1e-8
  if (sum(keep) < 2) {
    # Degenerate fold: fall back to majority class prediction
    maj <- as.numeric(mean(y_train) >= 0.5)
    return(maj)
  }
  X_train <- X_train[, keep, drop = FALSE]
  X_test  <- X_test[,  keep, drop = FALSE]
  
  # ── Also check within-class variance ─────────────────────────────
  classes <- unique(y_train)
  for (cl in classes) {
    cl_sd  <- apply(X_train[y_train == cl, , drop = FALSE], 2, sd)
    bad    <- cl_sd < 1e-8
    if (any(bad)) {
      X_train <- X_train[, !bad, drop = FALSE]
      X_test  <- X_test[,  !bad, drop = FALSE]
    }
  }
  if (ncol(X_train) < 2) {
    maj <- as.numeric(mean(y_train) >= 0.5)
    return(maj)
  }
  
  # ── PCA path (Rep A) ─────────────────────────────────────────────
  if (use_pca) {
    max_M  <- min(nrow(X_train) - 1, ncol(X_train), 30)
    if (max_M < 1) {
      maj <- as.numeric(mean(y_train) >= 0.5)
      return(maj)
    }
    cv_errs <- sapply(1:max_M, function(M) {
      folds <- cut(seq_len(nrow(X_train)), breaks = 5, labels = FALSE)
      mean(sapply(1:5, function(f) {
        tr  <- X_train[folds != f, , drop = FALSE]
        te  <- X_train[folds == f, , drop = FALSE]
        ytr <- y_train[folds != f]
        yte <- y_train[folds == f]
        # Guard against singular folds within inner CV
        tryCatch({
          pca_tr <- prcomp(tr, center = TRUE, scale. = FALSE)
          Z_tr   <- pca_tr$x[, 1:M, drop = FALSE]
          Z_te   <- predict(pca_tr, te)[, 1:M, drop = FALSE]
          df_tr  <- data.frame(y = as.factor(ytr), Z_tr)
          lda_m  <- lda(y ~ ., data = df_tr)
          pred   <- predict(lda_m, data.frame(Z_te))$class
          mean(as.numeric(pred) - 1 != yte)
        }, error = function(e) 0.5)
      }))
    })
    M_best  <- which.min(cv_errs)
    pca_obj <- prcomp(X_train, center = TRUE, scale. = FALSE)
    Z_train <- pca_obj$x[, 1:M_best, drop = FALSE]
    Z_test  <- predict(pca_obj, X_test)[, 1:M_best, drop = FALSE]
    df_tr   <- data.frame(y = as.factor(y_train), Z_train)
    df_te   <- data.frame(Z_test)
    
  } else {
    # ── Plain LDA path (Reps B and C) ──────────────────────────────
    df_tr <- data.frame(y = as.factor(y_train), X_train)
    df_te <- data.frame(X_test)
  }
  
  # ── Fit LDA with tryCatch in case of remaining singularity ───────
  result <- tryCatch({
    lda_m <- lda(y ~ ., data = df_tr)
    pred  <- predict(lda_m, df_te)$class
    as.numeric(as.character(pred))
  }, error = function(e) {
    # Fallback: majority class
    as.numeric(mean(y_train) >= 0.5)
  })
  
  result
}

# 7d. Single hidden-layer neural network
fit_nnet <- function(X_train, y_train, X_test) {
  p <- ncol(X_train)
  
  # Cap hidden units so total weights stay under a safe limit
  # weights = p*H + H + H*1 + 1 = H*(p+2) + 1
  max_weights <- 3000
  H_max       <- floor((max_weights - 1) / (p + 2))
  H_max       <- max(2, min(H_max, 10))   # at least 2, at most 10
  H_grid      <- unique(c(2, 4, min(6, H_max), H_max))
  H_grid      <- H_grid[H_grid <= H_max]
  
  decay_grid  <- c(0.001, 0.01, 0.1)
  best_err    <- Inf
  best_size   <- H_grid[1]
  best_decay  <- 0.01
  
  folds <- cut(seq_len(nrow(X_train)), breaks = 5, labels = FALSE)
  
  for (H in H_grid) {
    for (d in decay_grid) {
      err <- mean(sapply(1:5, function(f) {
        tr   <- X_train[folds != f, , drop = FALSE]
        te   <- X_train[folds == f, , drop = FALSE]
        ytr  <- y_train[folds != f]
        yte  <- y_train[folds == f]
        nn   <- nnet(tr, ytr, size = H, decay = d,
                     linout   = FALSE,
                     maxit    = 300,
                     MaxNWts  = max_weights,
                     trace    = FALSE)
        pred <- as.numeric(predict(nn, te) > 0.5)
        mean(pred != yte)
      }))
      if (err < best_err) {
        best_err <- err; best_size <- H; best_decay <- d
      }
    }
  }
  
  best_nn  <- NULL
  best_val <- Inf
  for (r in 1:5) {
    nn  <- nnet(X_train, y_train,
                size    = best_size,
                decay   = best_decay,
                linout  = FALSE,
                maxit   = 300,
                MaxNWts = max_weights,
                trace   = FALSE)
    if (nn$value < best_val) { best_val <- nn$value; best_nn <- nn }
  }
  
  as.numeric(predict(best_nn, X_test) > 0.5)
}


# ── 8. LOOCV — THE MAIN LOOP ────────────────────────────────────
cat("\n=== Starting LOOCV (84 folds) ===\n")
cat("This will take several minutes...\n\n")

# Storage: 13 models (M1-M12 + meta) + probabilities for stacking
n_models <- 12
pred_mat  <- matrix(NA, nrow = n_pat, ncol = n_models,
                    dimnames = list(patient_ids,
                                    paste0("M", 1:n_models)))
prob_mat  <- matrix(NA, nrow = n_pat, ncol = n_models,
                    dimnames = list(patient_ids,
                                    paste0("M", 1:n_models)))

tau_grid <- c(0.1, 0.5, 1, 2, 5, Inf)

for (fold in seq_len(n_pat)) {
  
  test_pid   <- patient_ids[fold]
  train_pids <- patient_ids[-fold]
  
  # ── Split nuclei ──────────────────────────────────────────────
  train_idx <- bdat$id %in% train_pids
  test_idx  <- bdat$id == test_pid
  
  train_nuclei_raw <- as.matrix(bdat[train_idx, feature_cols])
  test_nuclei_raw  <- as.matrix(bdat[test_idx,  feature_cols])
  
  # ── Standardise (training stats only) ────────────────────────
  std <- standardise(train_nuclei_raw, test_nuclei_raw)
  train_nuclei <- std$train
  test_nuclei  <- std$test
  
  # Standardise all nuclei for scoring (need test patient scored too)
  all_nuclei_std <- sweep(
    sweep(as.matrix(bdat[, feature_cols]), 2, std$mu, "-"),
    2, std$sig, "/"
  )
  
  # Training nucleus labels (patient group propagated to nuclei)
  train_labels_nuc <- bdat$group[train_idx]
  y_train          <- y_pat[train_pids]
  y_test           <- y_pat[test_pid]
  
  # ── Train nucleus scorer (for Reps B and C) ───────────────────
  nucleus_scores <- train_nucleus_scorer(
    train_nuclei, train_labels_nuc, all_nuclei_std, feature_cols
  )
  
  # ── Rep A: standardised patient-level summaries ───────────────
  # Build for train patients inside fold (using fold-specific scaling)
  make_rep_A_fold <- function(pids, nuclei_std_all) {
    t(sapply(pids, function(pid) {
      idx <- which(bdat$id == pid)
      sub <- nuclei_std_all[idx, , drop = FALSE]
      c(colMeans(sub),
        apply(sub, 2, sd),
        apply(sub, 2, quantile, probs = 0.10, type = 7),
        apply(sub, 2, quantile, probs = 0.50, type = 7),
        apply(sub, 2, quantile, probs = 0.90, type = 7))
    }))
  }
  
  A_train <- make_rep_A_fold(train_pids, all_nuclei_std)
  A_test  <- make_rep_A_fold(test_pid,   all_nuclei_std)
  
  # ── Tune tau for Rep B using 5-fold CV on training patients ───
  best_tau <- tau_grid[1]
  best_tau_err <- Inf
  inner_folds <- cut(seq_along(train_pids), breaks = 5, labels = FALSE)
  
  for (tau in tau_grid) {
    tau_errs <- sapply(1:5, function(f) {
      inner_train <- train_pids[inner_folds != f]
      inner_test  <- train_pids[inner_folds == f]
      # Use pre-computed nucleus scores (approximation within inner CV)
      B_inner_train <- make_rep_B(
        bdat, inner_train, feature_cols, nucleus_scores, tau)
      B_inner_test  <- make_rep_B(
        bdat, inner_test,  feature_cols, nucleus_scores, tau)
      y_inner_train <- y_pat[inner_train]
      y_inner_test  <- y_pat[inner_test]
      # Quick RF to evaluate tau
      rf_tmp <- randomForest(x = B_inner_train,
                             y = as.factor(y_inner_train), ntree = 100)
      pred_tmp <- as.numeric(predict(rf_tmp, B_inner_test) == "1")
      mean(pred_tmp != y_inner_test)
    })
    if (mean(tau_errs) < best_tau_err) {
      best_tau_err <- mean(tau_errs)
      best_tau     <- tau
    }
  }
  
  # ── Build Reps B and C with best tau ─────────────────────────
  all_pids_fold <- c(train_pids, test_pid)
  B_all   <- make_rep_B(bdat, all_pids_fold, feature_cols,
                        nucleus_scores, best_tau)
  B_train <- B_all[train_pids, , drop = FALSE]
  B_test  <- B_all[test_pid,   , drop = FALSE]
  
  C_all   <- make_rep_C(bdat, all_pids_fold, nucleus_scores)
  C_train <- C_all[train_pids, , drop = FALSE]
  C_test  <- C_all[test_pid,   , drop = FALSE]
  
  # ── Run 12 classifiers ──────────────────────────────────────
  # M1-M4: Rep A
  pred_mat[fold, "M1"]  <- fit_lasso(A_train, y_train, A_test)
  pred_mat[fold, "M2"]  <- fit_rf(A_train, y_train, A_test)
  pred_mat[fold, "M3"]  <- fit_pca_lda(A_train, y_train, A_test,
                                       use_pca = TRUE)
  pred_mat[fold, "M4"]  <- fit_nnet(A_train, y_train, A_test)
  
  # M5-M8: Rep B
  pred_mat[fold, "M5"]  <- fit_lasso(B_train, y_train, B_test)
  pred_mat[fold, "M6"]  <- fit_rf(B_train, y_train, B_test)
  pred_mat[fold, "M7"]  <- fit_pca_lda(B_train, y_train, B_test,
                                       use_pca = FALSE)
  pred_mat[fold, "M8"]  <- fit_nnet(B_train, y_train, B_test)
  
  # M9-M12: Rep C
  pred_mat[fold, "M9"]  <- fit_lasso(C_train, y_train, C_test)
  pred_mat[fold, "M10"] <- fit_rf(C_train, y_train, C_test)
  pred_mat[fold, "M11"] <- fit_pca_lda(C_train, y_train, C_test,
                                       use_pca = FALSE)
  pred_mat[fold, "M12"] <- fit_nnet(C_train, y_train, C_test)
  
  # Store predicted probabilities for stacking
  # (re-use model fits — use prob for lasso/RF, pred for LDA/nnet)
  cv_lasso_A <- cv.glmnet(A_train, y_train, family = "binomial",
                          alpha = 1, nfolds = 5)
  prob_mat[fold, "M1"] <- predict(cv_lasso_A, A_test,
                                  s = "lambda.min", type = "response")
  rf_A <- randomForest(x = A_train, y = as.factor(y_train), ntree = 500)
  prob_mat[fold, "M2"] <- as.numeric(
    predict(rf_A, A_test, type = "prob")[, "1"])
  prob_mat[fold, "M3"] <- pred_mat[fold, "M3"]   # LDA: use 0/1 as proxy
  prob_mat[fold, "M4"] <- pred_mat[fold, "M4"]
  
  cv_lasso_B <- cv.glmnet(B_train, y_train, family = "binomial",
                          alpha = 1, nfolds = 5)
  prob_mat[fold, "M5"] <- predict(cv_lasso_B, B_test,
                                  s = "lambda.min", type = "response")
  rf_B <- randomForest(x = B_train, y = as.factor(y_train), ntree = 500)
  prob_mat[fold, "M6"] <- as.numeric(
    predict(rf_B, B_test, type = "prob")[, "1"])
  prob_mat[fold, "M7"] <- pred_mat[fold, "M7"]
  prob_mat[fold, "M8"] <- pred_mat[fold, "M8"]
  
  cv_lasso_C <- cv.glmnet(C_train, y_train, family = "binomial",
                          alpha = 1, nfolds = 5)
  prob_mat[fold, "M9"] <- predict(cv_lasso_C, C_test,
                                  s = "lambda.min", type = "response")
  rf_C <- randomForest(x = C_train, y = as.factor(y_train), ntree = 500)
  prob_mat[fold, "M10"] <- as.numeric(
    predict(rf_C, C_test, type = "prob")[, "1"])
  prob_mat[fold, "M11"] <- pred_mat[fold, "M11"]
  prob_mat[fold, "M12"] <- pred_mat[fold, "M12"]
  
  if (fold %% 10 == 0)
    cat(sprintf("  Fold %d / %d complete\n", fold, n_pat))
}

cat("\nLOOCV complete.\n")

# ── 9. Meta-ensemble (stacking) ──────────────────────────────────
cat("\nFitting meta-ensemble...\n")
meta_preds <- numeric(n_pat)

for (fold in seq_len(n_pat)) {
  phi_train <- prob_mat[-fold, , drop = FALSE]
  phi_test  <- prob_mat[fold,  , drop = FALSE]
  y_meta    <- y_pat[-fold]
  
  meta_cv <- cv.glmnet(phi_train, y_meta, family = "binomial",
                       alpha = 0, nfolds = 5)   # ridge for meta-learner
  meta_p  <- predict(meta_cv, phi_test,
                     s = "lambda.min", type = "response")
  meta_preds[fold] <- as.numeric(meta_p > 0.5)
}


# ── 10. Performance metrics ──────────────────────────────────────
calc_metrics <- function(preds, truth) {
  mce  <- mean(preds != truth)
  ci   <- binom.test(sum(preds != truth), length(truth))$conf.int
  sens <- mean(preds[truth == 1] == 1)   # sensitivity
  spec <- mean(preds[truth == 0] == 0)   # specificity
  list(mce = mce, ci_lo = ci[1], ci_hi = ci[2],
       sens = sens, spec = spec)
}

# M0: majority class (always predict NR = 0)
m0_preds <- rep(0, n_pat)

all_preds <- cbind(
  M0  = m0_preds,
  pred_mat,
  META = meta_preds
)

results <- t(apply(all_preds, 2, function(p) {
  m <- calc_metrics(p, y_pat)
  c(MCE      = round(m$mce, 3),
    CI_lo    = round(m$ci_lo, 3),
    CI_hi    = round(m$ci_hi, 3),
    Sens     = round(m$sens, 3),
    Spec     = round(m$spec, 3))
}))

cat("\n=== RESULTS TABLE ===\n")
print(results)


# ── 11. McNemar's test (pairwise, Bonferroni corrected) ─────────
cat("\n=== McNEMAR'S TEST (selected pairs) ===\n")

mcnemar_pair <- function(pred_a, pred_b, truth) {
  correct_a <- pred_a == truth
  correct_b <- pred_b == truth
  n01 <- sum(!correct_a &  correct_b)   # A wrong, B right
  n10 <- sum( correct_a & !correct_b)   # A right, B wrong
  if (n01 + n10 == 0) return(NA)
  mcnemar.test(matrix(c(sum( correct_a & correct_b),
                        n10, n01,
                        sum(!correct_a & !correct_b)), 2, 2))$p.value
}

# Compare each model vs. M0 (lower bound) and vs. META (best)
n_comp    <- ncol(all_preds) - 1   # exclude M0 from comparisons with itself
p_vals    <- numeric(n_comp)
model_names <- colnames(all_preds)[-1]   # M1 through META

for (i in seq_along(model_names)) {
  p_vals[i] <- mcnemar_pair(
    all_preds[, model_names[i]],
    all_preds[, "M0"],
    y_pat
  )
}
p_adj <- p.adjust(p_vals, method = "bonferroni")
mcnemar_results <- data.frame(
  Model     = model_names,
  p_vs_M0   = round(p_vals, 4),
  p_adj     = round(p_adj, 4),
  sig       = ifelse(p_adj < 0.05, "*", "")
)
print(mcnemar_results)


# ── 12. Leakage demonstration ────────────────────────────────────
cat("\n=== LEAKAGE DEMONSTRATION (M2: RF on Rep A) ===\n")

# Correct: patient-level LOOCV (already computed as M2)
mce_valid <- results["M2", "MCE"]

# Wrong: nucleus-level 10-fold CV (leaky)
all_feat   <- as.matrix(bdat[, feature_cols])
all_labels <- as.factor(ifelse(bdat$group == "R", 1, 0))

nuc_folds  <- sample(rep(1:10, length.out = nrow(bdat)))
nuc_preds  <- character(nrow(bdat))

for (f in 1:10) {
  tr_feat  <- all_feat[nuc_folds != f, ]
  tr_lab   <- all_labels[nuc_folds != f]
  te_feat  <- all_feat[nuc_folds == f, ]
  rf_leak  <- randomForest(x = tr_feat, y = tr_lab, ntree = 300)
  nuc_preds[nuc_folds == f] <- as.character(predict(rf_leak, te_feat))
}
mce_leaky <- mean(nuc_preds != as.character(all_labels))

cat(sprintf("  Valid patient-level LOOCV MCE (M2): %.3f\n", mce_valid))
cat(sprintf("  Leaky nucleus-level 10-fold CV MCE: %.3f\n", mce_leaky))
cat(sprintf("  Leakage bias (overoptimism):        %.3f\n",
            mce_valid - mce_leaky))


# ── 13. Interpretability: feature importance ─────────────────────
cat("\n=== FEATURE IMPORTANCE (RF on Rep A, all 84 patients) ===\n")

Rep_A_all <- make_rep_A(bdat, patient_ids, feature_cols)
std_all   <- standardise(Rep_A_all, Rep_A_all)
Rep_A_std <- std_all$train

rf_interp <- randomForest(
  x          = Rep_A_std,
  y          = as.factor(y_pat),
  ntree      = 500,
  importance = TRUE
)

imp <- importance(rf_interp, type = 1)   # mean decrease in accuracy
imp_df <- data.frame(
  feature    = rownames(imp),
  importance = imp[, 1]
)
imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]
cat("Top 20 features by permutation importance:\n")
print(head(imp_df, 20))


# ── 14. Nucleus-level saliency (Rep B softmax weights) ───────────
# ── 14. Nucleus-level saliency (Rep B softmax weights) ───────────
cat("\n=== NUCLEUS-LEVEL SALIENCY ===\n")

# Use the nucleus scores already stored from the LOOCV loop
# Instead of refitting, aggregate the fold-level scores using
# the last fold's scorer as an approximation for illustration

# Safe standardisation: handle zero-variance columns
feat_mat <- as.matrix(bdat[, feature_cols])
feat_sd  <- apply(feat_mat, 2, sd)
feat_sd[feat_sd == 0] <- 1          # replace 0 SD to avoid NaN
feat_mu  <- colMeans(feat_mat)
all_nuclei_std_full <- sweep(
  sweep(feat_mat, 2, feat_mu, "-"), 2, feat_sd, "/")

# Verify no NAs
cat("NAs in standardised matrix:", sum(is.na(all_nuclei_std_full)), "\n")

rf_scorer_full <- randomForest(
  x     = all_nuclei_std_full,
  y     = as.factor(bdat$group),
  ntree = 300
)
nucleus_scores_full <- predict(rf_scorer_full, all_nuclei_std_full,
                               type = "prob")[, "R"]

tau_final <- 1
max_weights_sal <- sapply(patient_ids, function(pid) {
  idx <- which(bdat$id == pid)
  w   <- softmax_weights(nucleus_scores_full[idx], tau_final)
  max(w)
})

saliency_df <- data.frame(
  patient = patient_ids,
  group   = patient_labels,
  max_w   = max_weights_sal
)

cat("Mean max nucleus weight by group:\n")
print(tapply(saliency_df$max_w, saliency_df$group, mean))


# ── 15. Save results ─────────────────────────────────────────────
save(results, pred_mat, prob_mat, meta_preds,
     imp_df, saliency_df, mcnemar_results,
     file = "karyometry_results.RData")

cat("\n=== Analysis complete. Results saved to karyometry_results.RData ===\n")

getwd()
# ── 16. Basic plots ──────────────────────────────────────────────
# (run interactively)

# --- MCE comparison plot ---
# mce_vals <- results[, "MCE"]
# ci_lo    <- results[, "CI_lo"]
# ci_hi    <- results[, "CI_hi"]
# par(mar = c(5, 5, 3, 2))
# plot(mce_vals, 1:length(mce_vals), xlim = c(0, 0.7),
#      pch = 19, xlab = "MCE (lower is better)",
#      ylab = "", yaxt = "n",
#      main = "LOOCV Misclassification Error with 95% CI")
# axis(2, at = 1:length(mce_vals), labels = names(mce_vals), las = 2)
# segments(ci_lo, 1:length(mce_vals), ci_hi, 1:length(mce_vals))
# abline(v = results["M0", "MCE"], lty = 2, col = "red")

# --- Feature importance plot ---
# barplot(head(imp_df$importance, 20),
#         names.arg = head(imp_df$feature, 20),
#         las = 2, cex.names = 0.7,
#         main = "Top 20 Features: Permutation Importance",
#         ylab = "Mean Decrease in Accuracy",
#         col  = "steelblue")





























# ================================================================
#  BIOS 648 Project — Karyometry Analysis
#  FIGURES SCRIPT
#  Run this AFTER the main analysis script.
#  All objects (results, pred_mat, prob_mat, imp_df, etc.)
#  must already be in your R environment.
#  If starting fresh: load("karyometry_results.RData")
# ================================================================

# ── Output folder ───────────────────────────────────────────────
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# ── Colour palette (consistent across all figures) ──────────────
col_A    <- "#2166AC"   # Rep A  — blue
col_B    <- "#1B7837"   # Rep B  — green
col_C    <- "#D6604D"   # Rep C  — red/orange
col_META <- "#762A83"   # META   — purple
col_M0   <- "#878787"   # M0     — grey
col_ref  <- "#BABABA"   # reference lines

# ================================================================
#  FIGURE 1 — MCE forest plot with 95% CI
#  All 13 models ordered by MCE, colour-coded by representation
# ================================================================
png(file.path(fig_dir, "Fig1_MCE_forestplot.png"),
    width = 2200, height = 1800, res = 220)

# ── Data ──────────────────────────────────────────────────────────
mce   <- results[, "MCE"]
ci_lo <- results[, "CI_lo"]
ci_hi <- results[, "CI_hi"]
mod_names <- rownames(results)

# Order by MCE (best at top)
ord <- order(mce)
mce_o   <- mce[ord]
lo_o    <- ci_lo[ord]
hi_o    <- ci_hi[ord]
nam_o   <- mod_names[ord]

# Colour by representation
rep_col <- ifelse(nam_o == "M0",   col_M0,
                  ifelse(nam_o %in% paste0("M", 1:4),  col_A,
                         ifelse(nam_o %in% paste0("M", 5:8),  col_B,
                                ifelse(nam_o %in% paste0("M", 9:12), col_C,
                                       col_META))))

n <- length(mce_o)
par(mar = c(5, 7, 4, 3))
plot(mce_o, 1:n,
     xlim  = c(0, 0.85),
     ylim  = c(0.5, n + 0.5),
     pch   = 19, cex = 1.3,
     col   = rep_col,
     xlab  = "Misclassification Error (LOOCV)",
     ylab  = "",
     yaxt  = "n",
     main  = "Figure 1. LOOCV Misclassification Error with 95% Clopper-Pearson CI",
     cex.main = 0.95, font.main = 1,
     cex.lab  = 0.9)

# Axis labels
axis(2, at = 1:n, labels = nam_o, las = 2, cex.axis = 0.85)

# CI whiskers
segments(lo_o, 1:n, hi_o, 1:n, col = rep_col, lwd = 1.8)
# End caps
segments(lo_o, 1:n - 0.18, lo_o, 1:n + 0.18, col = rep_col, lwd = 1.5)
segments(hi_o, 1:n - 0.18, hi_o, 1:n + 0.18, col = rep_col, lwd = 1.5)

# Reference line: majority class
abline(v = results["M0", "MCE"], lty = 2, col = col_ref, lwd = 1.5)
text(results["M0", "MCE"] + 0.01, n * 0.15,
     "Majority class\nbaseline", col = col_ref, cex = 0.72, adj = 0)

# Legend
legend("topleft",
       legend = c("M0 — Majority class",
                  "Rep A — Mean summary",
                  "Rep B — Score-weighted",
                  "Rep C — Prob. aggregation",
                  "META — Stacking ensemble"),
       col    = c(col_M0, col_A, col_B, col_C, col_META),
       pch    = 19, lwd = 1.8, pt.cex = 1.1,
       cex    = 0.78, bty = "n")

# MCE values as text
text(mce_o + 0.02, 1:n, labels = sprintf("%.3f", mce_o),
     cex = 0.72, col = rep_col)

dev.off()
cat("Figure 1 saved.\n")


# ================================================================
#  FIGURE 2 — Sensitivity / Specificity comparison
#  Grouped bar chart: all models (excl META), side-by-side
# ================================================================
png(file.path(fig_dir, "Fig2_SensSPec.png"),
    width = 2400, height = 1600, res = 220)

models_show <- c("M0", paste0("M", 1:12), "META")
sens_v <- results[models_show, "Sens"]
spec_v <- results[models_show, "Spec"]

n_m  <- length(models_show)
xpos <- seq_len(n_m)

par(mar = c(6, 5, 4, 3))
bar_mat <- rbind(sens_v, spec_v)
bp <- barplot(bar_mat,
              beside      = TRUE,
              col         = c("#4DAF4A", "#377EB8"),
              names.arg   = models_show,
              las         = 2,
              ylim        = c(0, 1.12),
              ylab        = "Proportion",
              main        = "Figure 2. Sensitivity and Specificity by Model (LOOCV)",
              cex.main    = 0.95, font.main = 1,
              cex.names   = 0.72,
              cex.axis    = 0.82,
              border      = NA)

# Reference line at 0.5
abline(h = 0.5, lty = 2, col = col_ref, lwd = 1.2)

# Legend
legend("topright",
       legend = c("Sensitivity (TPR for R)",
                  "Specificity (TNR for NR)"),
       fill   = c("#4DAF4A", "#377EB8"),
       bty    = "n", cex = 0.82, border = NA)

# Representation group labels
mids_A    <- colMeans(bp[, models_show %in% paste0("M", 1:4)])
mids_B    <- colMeans(bp[, models_show %in% paste0("M", 5:8)])
mids_C    <- colMeans(bp[, models_show %in% paste0("M", 9:12)])
mids_META <- bp[1, models_show == "META"]

mtext("Rep A",  side = 1, line = 4.5,
      at = mean(mids_A),  col = col_A,  cex = 0.75, font = 2)
mtext("Rep B",  side = 1, line = 4.5,
      at = mean(mids_B),  col = col_B,  cex = 0.75, font = 2)
mtext("Rep C",  side = 1, line = 4.5,
      at = mean(mids_C),  col = col_C,  cex = 0.75, font = 2)
mtext("META",   side = 1, line = 4.5,
      at = mean(mids_META), col = col_META, cex = 0.75, font = 2)

dev.off()
cat("Figure 2 saved.\n")


# ================================================================
#  FIGURE 3 — ROC curves
#  Best model from each representation + META
#  Uses prob_mat (predicted probabilities from LOOCV)
# ================================================================

png(file.path(fig_dir, "Fig3_PredProb_Distribution.png"),
    width = 2000, height = 1400, res = 220)

par(mfrow = c(1, 2), mar = c(5, 5, 5, 2), oma = c(0, 0, 3, 0))

# Grouping label — explicit factor with both levels
grp <- factor(ifelse(y_pat == 1, "R (Recurrent)", "NR (Non-recurrent)"),
              levels = c("NR (Non-recurrent)", "R (Recurrent)"))

panel_cols   <- c(adjustcolor(col_A, 0.55), adjustcolor(col_C, 0.55))
panel_border <- c(col_A, col_C)

# Panel 1: M6 — RF on Rep B (best calibrated probabilities)
boxplot(prob_mat[, "M6"] ~ grp,
        col    = panel_cols,
        border = panel_border,
        ylab   = "Predicted Recurrence Probability",
        xlab   = "True Patient Group",
        main   = "M6 — Random Forest on Rep B\n(MCE = 0.452)",
        ylim   = c(0, 1),
        cex.main = 0.88, font.main = 1,
        cex.lab  = 0.85, cex.axis = 0.82,
        outline  = TRUE)
abline(h = 0.5, lty = 2, col = col_ref, lwd = 1.3)
text(1.5, 0.53, "Decision\nthreshold", cex = 0.65, col = col_ref)

# Wilcoxon test p-value
wt6 <- wilcox.test(prob_mat[y_pat == 0, "M6"],
                   prob_mat[y_pat == 1, "M6"])
mtext(sprintf("Wilcoxon p = %.3f", wt6$p.value),
      side = 3, line = 0.2, cex = 0.72, col = "#555555")

# Panel 2: M2 — RF on Rep A (for comparison)
boxplot(prob_mat[, "M2"] ~ grp,
        col    = panel_cols,
        border = panel_border,
        ylab   = "Predicted Recurrence Probability",
        xlab   = "True Patient Group",
        main   = "M2 — Random Forest on Rep A\n(MCE = 0.607)",
        ylim   = c(0, 1),
        cex.main = 0.88, font.main = 1,
        cex.lab  = 0.85, cex.axis = 0.82,
        outline  = TRUE)
abline(h = 0.5, lty = 2, col = col_ref, lwd = 1.3)

wt2 <- wilcox.test(prob_mat[y_pat == 0, "M2"],
                   prob_mat[y_pat == 1, "M2"])
mtext(sprintf("Wilcoxon p = %.3f", wt2$p.value),
      side = 3, line = 0.2, cex = 0.72, col = "#555555")

# Overall title
mtext("Figure 3. Predicted Recurrence Probabilities by True Patient Group\n(LOOCV out-of-fold predictions — Rep B vs. Rep A comparison)",
      side = 3, outer = TRUE, line = 0.5, cex = 0.85, font = 1)

# Shared legend
legend("bottomright",
       legend = c("NR (Non-recurrent)", "R (Recurrent)"),
       fill   = panel_cols,
       border = panel_border,
       bty    = "n", cex = 0.80)

dev.off()
cat("Figure 3 saved.\n")


# Quick check — are probabilities higher for R patients?
tapply(prob_mat[, "M5"], y_pat, mean)
tapply(prob_mat[, "M1"], y_pat, mean)

# ================================================================
#  FIGURE 4 — Feature importance (top 20)
#  Horizontal bar chart, colour-coded by statistic type
# ================================================================
png(file.path(fig_dir, "Fig4_FeatureImportance.png"),
    width = 2200, height = 1800, res = 220)

top20 <- head(imp_df, 20)

# Identify statistic suffix for colouring
stat_type <- sub(".*_", "", top20$feature)   # mean / sd / q10 / q50 / q90
stat_cols  <- c(mean = "#4393C3",
                sd   = "#D6604D",
                q10  = "#74C476",
                q50  = "#FD8D3C",
                q90  = "#9E9AC8")
bar_cols <- stat_cols[stat_type]

par(mar = c(5, 9, 4, 3))
bp4 <- barplot(rev(top20$importance),
               horiz      = TRUE,
               names.arg  = rev(top20$feature),
               las        = 1,
               col        = rev(bar_cols),
               border     = NA,
               xlab       = "Mean Decrease in Accuracy (Permutation Importance)",
               main       = "Figure 4. Top 20 Karyometric Features by Permutation Importance\n(Random Forest on Representation A, all 84 patients)",
               cex.main   = 0.88, font.main = 1,
               cex.names  = 0.78,
               cex.axis   = 0.82)

# Value labels
text(rev(top20$importance) + 0.05, bp4,
     labels = round(rev(top20$importance), 2),
     cex = 0.68, adj = 0)

legend("bottomright",
       legend = c("Mean", "SD", "Q10 (10th pctile)",
                  "Q50 (median)", "Q90 (90th pctile)"),
       fill   = stat_cols,
       border = NA, bty = "n", cex = 0.78)

dev.off()
cat("Figure 4 saved.\n")


# ================================================================
#  FIGURE 5 — Leakage demonstration
#  Two-panel: (a) bar chart comparing MCE; (b) histogram of
#  nucleus-level scores to visualise the leakage mechanism
# ================================================================
png(file.path(fig_dir, "Fig5_Leakage.png"),
    width = 2400, height = 1400, res = 220)

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# ── Panel (a): MCE comparison ────────────────────────────────────
mce_compare <- c("Valid\npatient-level\nLOOCV" = 0.607,
                 "Leaky\nnucleus-level\n10-fold CV" = 0.268)
bar_cols5 <- c("#2166AC", "#D6604D")

bp5 <- barplot(mce_compare,
               col    = bar_cols5,
               border = NA,
               ylim   = c(0, 0.85),
               ylab   = "Misclassification Error",
               main   = "(a) Valid vs. Leaky CV\n(M2: Random Forest on Rep A)",
               cex.main  = 0.88, font.main = 1,
               cex.names = 0.80)

# Value labels on bars
text(bp5, mce_compare + 0.03,
     labels = sprintf("MCE = %.3f", mce_compare),
     cex = 0.82, font = 2,
     col = bar_cols5)

# Bias annotation
arrows(bp5[1], 0.64, bp5[2], 0.30,
       length = 0.12, angle = 25, lwd = 1.8, col = "#666666")
text(mean(bp5), 0.50,
     sprintf("Bias = %.3f\n(%.0f%% overoptimism)",
             0.607 - 0.268, (0.607 - 0.268) / 0.607 * 100),
     cex = 0.78, col = "#333333")

# Reference line: majority class
abline(h = 0.464, lty = 2, col = col_ref, lwd = 1.4)
text(bp5[1] - 0.35, 0.474, "Majority\nbaseline", col = col_ref, cex = 0.68)

# ── Panel (b): nucleus scores distribution by group ──────────────
# Use the nucleus_scores_full computed in Block 14
# If not available, recompute here
if (!exists("nucleus_scores_full")) {
  feat_mat <- as.matrix(bdat[, feature_cols])
  feat_sd  <- apply(feat_mat, 2, sd); feat_sd[feat_sd == 0] <- 1
  feat_mu  <- colMeans(feat_mat)
  all_std  <- sweep(sweep(feat_mat, 2, feat_mu, "-"), 2, feat_sd, "/")
  rf_tmp   <- randomForest(x = all_std, y = as.factor(bdat$group),
                           ntree = 300)
  nucleus_scores_full <- predict(rf_tmp, all_std, type = "prob")[, "R"]
}

scores_R  <- nucleus_scores_full[bdat$group == "R"]
scores_NR <- nucleus_scores_full[bdat$group == "NR"]

# Overlapping histograms
h_NR <- hist(scores_NR, breaks = 30, plot = FALSE)
h_R  <- hist(scores_R,  breaks = 30, plot = FALSE)

plot(h_NR, col = adjustcolor(col_A, 0.55),
     border = "white",
     xlim   = c(0, 1),
     ylim   = c(0, max(c(h_NR$counts, h_R$counts)) * 1.15),
     xlab   = "Nucleus-level Recurrence Risk Score",
     ylab   = "Count",
     main   = "(b) Distribution of Nucleus Scores\nby Patient Group",
     cex.main = 0.88, font.main = 1)
plot(h_R, col = adjustcolor(col_C, 0.55),
     border = "white", add = TRUE)

abline(v = 0.5, lty = 2, col = col_ref, lwd = 1.4)

legend("topright",
       legend = c("NR patients", "R patients"),
       fill   = c(adjustcolor(col_A, 0.55), adjustcolor(col_C, 0.55)),
       border = "white", bty = "n", cex = 0.80)

dev.off()
cat("Figure 5 saved.\n")


# ================================================================
#  FIGURE 6 (BONUS) — Confusion matrices for key models
#  M0, M5 (best individual), META — tiled display
# ================================================================
png(file.path(fig_dir, "Fig7_Heatmap_MCE.png"),
    width = 1800, height = 1400, res = 220)

mce_grid <- matrix(
  results[paste0("M", 1:12), "MCE"],
  nrow = 3, ncol = 4, byrow = TRUE
)
row_labels <- c("Rep A: Mean+Quantile (460-dim)",
                "Rep B: Score-weighted (92-dim)",
                "Rep C: Prob.Aggregation (3-dim)")
col_labels <- c("Lasso LR", "Random Forest", "PCA+LDA", "Neural Net")

par(mar = c(5, 14, 4, 6))

# Colour scale
heat_cols <- colorRampPalette(c("#1B7837","#F7F7F7","#D6604D"))(100)
mce_range <- range(mce_grid)
col_idx <- function(v)
  round((v - mce_range[1]) / diff(mce_range) * 99) + 1

# Draw cells manually for full control
plot(0, type = "n",
     xlim = c(0.5, 4.5), ylim = c(0.5, 3.5),
     axes = FALSE, xlab = "Classifier", ylab = "",
     main = "Figure 7. LOOCV MCE: Representation \u00d7 Classifier",
     cex.main = 0.90, font.main = 1)

for (r in 1:3) for (c in 1:4) {
  val <- mce_grid[r, c]
  rect(c - 0.5, r - 0.5, c + 0.5, r + 0.5,
       col = heat_cols[col_idx(val)], border = "white", lwd = 2)
  text(c, r, sprintf("%.3f", val),
       cex = 1.0, font = 2,
       col = ifelse(val < 0.55, "white", "#222222"))
}

axis(1, at = 1:4, labels = col_labels, cex.axis = 0.82, tick = FALSE)
axis(2, at = 1:3, labels = row_labels, las = 2,
     cex.axis = 0.78, tick = FALSE)

# Colour legend
legend("right", inset = -0.02,
       legend = c(sprintf("High MCE (%.2f)", mce_range[2]),
                  "Medium",
                  sprintf("Low MCE  (%.2f)", mce_range[1])),
       fill   = c(heat_cols[100], heat_cols[50], heat_cols[1]),
       border = NA, bty = "n", cex = 0.75, xpd = TRUE)

dev.off()
cat("Figure 6 saved.\n")


# ================================================================
#  FIGURE 7 (BONUS) — Representation comparison heatmap
#  MCE by Representation × Classifier (3x4 grid)
# ================================================================
png(file.path(fig_dir, "Fig7_Heatmap_MCE.png"),
    width = 1800, height = 1400, res = 220)

mce_grid <- matrix(
  results[paste0("M", 1:12), "MCE"],
  nrow = 3, ncol = 4, byrow = TRUE
)
row_labels <- c("Rep A: Mean+Quantile (460-dim)",
                "Rep B: Score-weighted (92-dim)",
                "Rep C: Prob.Aggregation (3-dim)")
col_labels <- c("Lasso LR", "Random Forest", "PCA+LDA", "Neural Net")

par(mar = c(5, 14, 4, 6))

# Colour scale
heat_cols <- colorRampPalette(c("#1B7837","#F7F7F7","#D6604D"))(100)
mce_range <- range(mce_grid)
col_idx <- function(v)
  round((v - mce_range[1]) / diff(mce_range) * 99) + 1

# Draw cells manually for full control
plot(0, type = "n",
     xlim = c(0.5, 4.5), ylim = c(0.5, 3.5),
     axes = FALSE, xlab = "Classifier", ylab = "",
     main = "Figure 7. LOOCV MCE: Representation \u00d7 Classifier",
     cex.main = 0.90, font.main = 1)

for (r in 1:3) for (c in 1:4) {
  val <- mce_grid[r, c]
  rect(c - 0.5, r - 0.5, c + 0.5, r + 0.5,
       col = heat_cols[col_idx(val)], border = "white", lwd = 2)
  text(c, r, sprintf("%.3f", val),
       cex = 1.0, font = 2,
       col = ifelse(val < 0.55, "white", "#222222"))
}

axis(1, at = 1:4, labels = col_labels, cex.axis = 0.82, tick = FALSE)
axis(2, at = 1:3, labels = row_labels, las = 2,
     cex.axis = 0.78, tick = FALSE)

# Colour legend
legend("bottomright", inset = -0.02,
       legend = c(sprintf("High MCE (%.2f)", mce_range[2]),
                  "Medium",
                  sprintf("Low MCE  (%.2f)", mce_range[1])),
       fill   = c(heat_cols[100], heat_cols[50], heat_cols[1]),
       border = NA, bty = "n", cex = 0.75, xpd = TRUE)

dev.off()
cat("Figure 7 saved.\n")


# ================================================================
#  Summary
# ================================================================
cat("\n========================================\n")
cat("All figures saved to:", normalizePath(fig_dir), "\n")
cat("  Fig1_MCE_forestplot.png   — MCE with 95% CI (all models)\n")
cat("  Fig2_SensSPec.png         — Sensitivity & Specificity\n")
cat("  Fig3_ROC_curves.png       — ROC curves (key models)\n")
cat("  Fig4_FeatureImportance.png— Top 20 features\n")
cat("  Fig5_Leakage.png          — Leakage demonstration\n")
cat("  Fig6_ConfusionMatrices.png— Confusion matrices (M0, M5, META)\n")
cat("  Fig7_Heatmap_MCE.png      — MCE heatmap: Rep x Classifier\n")
cat("========================================\n")

