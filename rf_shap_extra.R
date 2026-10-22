# rf_shap_extra.R - v29 additions to rf_shap_tools.R
#
# 1. Save / load RF + SHAP runs (.RData + a small .meta.rds sidecar)
# 2. Analyses from the Rmd "RF_SHAP_all_regions" script:
#    cross-stage prediction, %IncMSE importance, mark and region totals,
#    stage comparison, correlated feature pairs, CV summary, methods text.
#
# Sourced by app_v29.R after rf_shap_tools.R. Uses only functions from that file
# plus ggplot2 (attached by the app) and base R.

rs_or <- function(a, b) if (is.null(a)) b else a

# ---- 1. Saving and loading -----------------------------------------------------

# A string that identifies a run: stages, conditions, samples, marks, the number of
# bins and segments in each object, and every setting that changes the result.
# If any of these change, the saved run is not reused.
rs_signature <- function(specs, cfg) {
  parts <- vapply(specs, function(sp) {
    paste(sp$label, sp$cond,
          paste(sort(sp$samples), collapse = ","),
          paste(sort(sp$marks), collapse = ","),
          paste0("bins=", length(sp$hmm$bins)),
          paste0("segs=", length(sp$hmm$segments)),
          sep = "|")
  }, character(1))
  paste(c(parts,
          paste0("ntree=", cfg$ntree, ",single=", cfg$ntree_single,
                 ",cv=", cfg$cv_folds, ",nshap=", cfg$n_shap,
                 ",seed=", cfg$seed, ",excl=", isTRUE(cfg$exclude_shared))),
        collapse = " ;; ")
}

rs_save <- function(res, specs, cfg, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base  <- paste0("rf_shap_", paste(make.names(names(res$fits)), collapse = "_vs_"),
                  "_", stamp)
  rdata <- file.path(dir, paste0(base, ".RData"))
  meta  <- list(
    signature = rs_signature(specs, cfg),
    created   = format(Sys.time(), "%Y-%m-%d %H:%M"),
    stages    = names(res$fits),
    n_genes   = sapply(res$fits, function(f) nrow(f$data)),
    marks     = unique(unlist(lapply(res$fits, function(f) f$marks))),
    cfg       = cfg)
  saveRDS(meta, file.path(dir, paste0(base, ".meta.rds")))
  save(res, file = rdata)
  rdata
}

# Lists the saved runs in a folder (newest first). Reads only the small meta files.
rs_list_saved <- function(dir) {
  empty <- data.frame(label = character(0), rdata = character(0),
                      signature = character(0), stringsAsFactors = FALSE)
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(empty)
  mf <- list.files(dir, pattern = "\\.meta\\.rds$", full.names = TRUE)
  if (length(mf) == 0) return(empty)
  metas <- lapply(mf, function(f) tryCatch(readRDS(f), error = function(e) NULL))
  rd    <- file.path(dir, sub("\\.meta\\.rds$", ".RData", basename(mf)))
  ok    <- !vapply(metas, is.null, logical(1)) & file.exists(rd)
  if (!any(ok)) return(empty)
  metas <- metas[ok]; rd <- rd[ok]
  out <- data.frame(
    label = vapply(seq_along(metas), function(i) {
      m <- metas[[i]]
      paste0(m$created, " | ", paste(m$stages, collapse = " vs "), " | ",
             paste(m$n_genes, collapse = "/"), " genes | ", basename(rd[i]))
    }, character(1)),
    rdata     = rd,
    signature = vapply(metas, function(m) rs_or(m$signature, ""), character(1)),
    stringsAsFactors = FALSE)
  out[order(out$label, decreasing = TRUE), , drop = FALSE]
}

# Path of the saved run whose signature matches, or NULL.
rs_find_saved <- function(dir, signature) {
  lst <- rs_list_saved(dir)
  hit <- lst$rdata[lst$signature == signature]
  if (length(hit)) hit[1] else NULL
}

rs_load_saved <- function(path) {
  e <- new.env()
  load(path, envir = e)
  if (!exists("res", envir = e, inherits = FALSE))
    stop("No results object found in ", basename(path))
  e$res
}

# ---- 2. Analyses from the Rmd ----------------------------------------------------

# Cross-stage prediction. Each stage's model is applied to its own held-out test genes
# and to all genes of the other stage (as in the Rmd script).
rs_cross_perf <- function(fits) {
  rows <- list()
  for (m in names(fits)) for (d in names(fits)) {
    fm    <- fits[[m]]
    newd  <- if (m == d) fm$data[fm$test_genes, , drop = FALSE] else fits[[d]]$data
    p     <- predict(fm$rf, newd[, fm$feats, drop = FALSE])
    rows[[length(rows) + 1]] <- data.frame(
      model   = m,
      applied_to = d,
      genes   = if (m == d) "held-out test genes" else "all genes of the other stage",
      n_genes = nrow(newd),
      r       = round(cor(p, newd$rpkm), 3),
      RMSE    = round(sqrt(mean((p - newd$rpkm)^2)), 3),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Classic random forest importance (%IncMSE, unscaled) for each feature.
rs_incmse <- function(fits) {
  do.call(rbind, lapply(names(fits), function(s) {
    imp <- fits[[s]]$rf$importance[, 1]
    data.frame(stage = s, feature = names(imp),
               mark = sub("__.*", "", names(imp)),
               incMSE = unname(imp), stringsAsFactors = FALSE)
  }))
}

# Feature pairs whose training-set Spearman |rho| is at or above the threshold.
rs_pairs_all <- function(fits, thresh = 0.7) {
  out <- lapply(names(fits), function(s) {
    m   <- fits[[s]]$cor
    idx <- which(upper.tri(m), arr.ind = TRUE)
    d   <- data.frame(stage = s,
                      feature_1 = rs_pretty(rownames(m)[idx[, 1]]),
                      feature_2 = rs_pretty(colnames(m)[idx[, 2]]),
                      rho = round(m[idx], 3), stringsAsFactors = FALSE)
    d[abs(d$rho) >= thresh, , drop = FALSE]
  })
  res <- do.call(rbind, out)
  if (is.null(res) || nrow(res) == 0)
    return(data.frame(note = paste0("No feature pairs with |rho| >= ", thresh)))
  res[order(-abs(res$rho)), , drop = FALSE]
}

# Mean and SD over the CV folds, per stage.
rs_cv_summary <- function(cv) {
  g <- split(cv, cv$stage)
  do.call(rbind, lapply(names(g), function(s) {
    d <- g[[s]]
    data.frame(stage = s,
               r_mean = mean(d$r), r_sd = sd(d$r),
               RMSE_mean = mean(d$RMSE), RMSE_sd = sd(d$RMSE),
               R2_mean = mean(d$R2), R2_sd = sd(d$R2),
               folds = nrow(d), stringsAsFactors = FALSE)
  }))
}

# Summed SHAP importance and weighted direction, per stage and mark (over its 8 regions).
rs_mark_totals <- function(fs) {
  d <- split(fs, list(as.character(fs$stage), as.character(fs$mark)), drop = TRUE)
  do.call(rbind, lapply(d, function(x) {
    data.frame(stage = as.character(x$stage[1]),
               mark = as.character(x$mark[1]),
               importance_SHAP = sum(x$mean_abs_SHAP),
               direction_rho = if (all(is.na(x$direction_rho))) NA_real_
                 else weighted.mean(x$direction_rho, x$mean_abs_SHAP, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}

# Summed SHAP importance per stage and region (over its 8 marks).
rs_region_totals <- function(fs) {
  d <- split(fs, list(as.character(fs$stage), as.character(fs$region)), drop = TRUE)
  do.call(rbind, lapply(d, function(x) {
    data.frame(stage = as.character(x$stage[1]),
               region = as.character(x$region[1]),
               importance_SHAP = sum(x$mean_abs_SHAP),
               stringsAsFactors = FALSE)
  }))
}

# ---- Figures -------------------------------------------------------------------

rs_plot_mark_shap <- function(fs) {
  d <- rs_mark_totals(fs)
  ggplot(d, aes(reorder(mark, importance_SHAP), importance_SHAP, fill = direction_rho)) +
    geom_col(colour = "grey30", width = 0.75) +
    coord_flip() +
    facet_wrap(~ stage) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-1, 1), na.value = "grey85", name = "Direction\n(rho)") +
    labs(x = NULL, y = "Importance (sum of mean |SHAP| over 8 regions)",
         title = "Feature importance per mark (SHAP, all regions)") +
    theme_minimal(base_size = 12)
}

rs_plot_incmse <- function(fits) {
  inc <- rs_incmse(fits)
  agg <- aggregate(incMSE ~ stage + mark, data = inc, FUN = sum)
  ggplot(agg, aes(reorder(mark, incMSE), incMSE, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "%IncMSE (summed over 8 regions)",
         title = "Classic random forest importance per mark (no direction)",
         fill = "Stage") +
    theme_minimal(base_size = 12)
}

rs_plot_region_totals <- function(fs) {
  d <- rs_region_totals(fs)
  d$region <- factor(d$region, levels = RS_LABELS)
  ggplot(d, aes(region, importance_SHAP, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    labs(x = NULL, y = "Summed mean |SHAP| (over 8 marks)",
         title = "Importance per region", fill = "Stage") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Mean |SHAP| of each mark @ region in stage A (x) vs stage B (y). Needs two stages.
rs_plot_stage_compare <- function(fs) {
  st <- unique(as.character(fs$stage))
  if (length(st) < 2) return(NULL)
  a <- fs[as.character(fs$stage) == st[1], c("feature", "region", "mark", "mean_abs_SHAP")]
  b <- fs[as.character(fs$stage) == st[2], c("feature", "region", "mean_abs_SHAP")]
  cmp <- merge(a, b, by = c("feature", "region"), suffixes = c("_A", "_B"))
  cmp$region <- as.character(cmp$region)
  rho <- cor(cmp$mean_abs_SHAP_A, cmp$mean_abs_SHAP_B, method = "spearman")
  top <- cmp[order(-pmax(cmp$mean_abs_SHAP_A, cmp$mean_abs_SHAP_B)), , drop = FALSE]
  top <- top[seq_len(min(10, nrow(top))), , drop = FALSE]
  top$label <- rs_pretty(top$feature)
  ggplot(cmp, aes(mean_abs_SHAP_A, mean_abs_SHAP_B, colour = region)) +
    geom_abline(linetype = "dashed", colour = "grey60") +
    geom_point(size = 2.5) +
    geom_text(data = top,
              mapping = aes(x = mean_abs_SHAP_A, y = mean_abs_SHAP_B, label = label),
              size = 2.8, vjust = -0.8, inherit.aes = FALSE) +
    labs(x = paste0("Importance in ", st[1], " (mean |SHAP|)"),
         y = paste0("Importance in ", st[2], " (mean |SHAP|)"),
         title = "Mark @ region importance: stage A vs stage B",
         subtitle = sprintf("Spearman correlation of importances: %.3f", rho)) +
    theme_minimal(base_size = 11)
}

# ---- Methods text (draft, filled from the run) -------------------------------------

rs_methods_text <- function(res) {
  fits <- res$fits
  cfg  <- res$cfg
  stg  <- names(fits)
  mk   <- unique(unlist(lapply(fits, function(f) f$marks)))
  unit <- rs_or(fits[[1]]$unit, "FPKM")
  nfeat <- length(mk) * length(RS_KEYS)
  counts <- paste(vapply(stg, function(s) {
    sprintf("%s: %d genes, %d held-out test genes", s,
            nrow(fits[[s]]$data), length(fits[[s]]$test_genes))
  }, character(1)), collapse = "; ")
  paste0(
    "Histone mark signal was summarised in eight regions per gene (TSS -2 to -1 kb, ",
    "TSS -1 kb to 0, TSS +/-200 bp, three equal gene-body bins, TES 0-1 kb and TES 1-2 kb), ",
    "giving ", nfeat, " features (", length(mk), " marks x ", length(RS_KEYS), " regions). ",
    "For each stage (", paste(stg, collapse = " and "), "), a regression random forest ",
    "(randomForest, ", cfg$ntree, " trees) was trained to predict log2(", unit, " + 1) on 80% ",
    "of genes and evaluated on the held-out 20%. Counts: ", counts, ". ",
    "Genes without any mark but with log2(", unit, " + 1) > 1 were excluded",
    if (isTRUE(cfg$exclude_shared)) ", as were genes whose gene bodies overlap another gene" else "",
    ". Model stability was checked with ", cfg$cv_folds, "-fold cross-validation on the ",
    "80% training set (mean and SD of r, RMSE and R2 across folds). Pairwise Spearman ",
    "correlations between the features were computed on the training set, and strongly ",
    "correlated pairs were taken into account when interpreting importance. Each stage's ",
    "model was also applied to the other stage's genes (cross-stage prediction). ",
    "Feature contributions were computed with exact TreeSHAP (treeshap R package; Lundberg ",
    "et al., 2020) for the test genes. Importance is defined as mean |SHAP|, and direction as ",
    "the Spearman correlation between a feature's value and its SHAP value, following the ",
    "positional SHAP analysis of Ramachandran & Ramakrishnan (2026).")
}

# ---- Output files (called from rs_write_outputs) ----------------------------------

rs_write_extra <- function(res, dir) {
  fits <- res$fits
  wc <- function(x, f) write.csv(x, file.path(dir, f), row.names = FALSE)
  wc(rs_cross_perf(fits), "cross_stage_performance.csv")
  wc(rs_pairs_all(fits, rs_or(res$cfg$cor_threshold, 0.7)), "feature_correlated_pairs.csv")
  wc(rs_cv_summary(res$cv), "cv_summary.csv")
  wc(rs_incmse(fits), "incmse_per_feature.csv")
  wc(rs_mark_totals(res$fs), "mark_importance_summary.csv")
  wc(rs_region_totals(res$fs), "region_importance_summary.csv")
  writeLines(rs_methods_text(res), file.path(dir, "methods_text.txt"))
  ggsave(file.path(dir, "Fig_mark_importance_SHAP.pdf"), rs_plot_mark_shap(res$fs),
         width = 7, height = 5)
  ggsave(file.path(dir, "Fig_incmse_per_mark.pdf"), rs_plot_incmse(fits),
         width = 7, height = 5)
  ggsave(file.path(dir, "Fig_region_totals.pdf"), rs_plot_region_totals(res$fs),
         width = 8, height = 5)
  p <- rs_plot_stage_compare(res$fs)
  if (!is.null(p)) ggsave(file.path(dir, "Fig_stage_comparison.pdf"), p, width = 7, height = 6)
  invisible(dir)
}
