# rf_shap_tools.R - region RF + SHAP engine, sourced by app27.R
#
# Self-contained: calls no helper defined in the app. Uses GenomicRanges,
# IRanges and ggplot2 (already attached by the app), and randomForest,
# caret, treeshap, ggbeeswarm, parallel and zip through :: only, so
# randomForest is never attached and cannot mask ggplot2::margin().

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
})

RS_KEYS   <- c("TSSm2k", "TSSm1k", "TSS200", "Body5", "BodyMid", "Body3", "TES1k", "TES2k")
RS_LABELS <- c("TSS-2kb", "TSS-1kb", "TSS+/-200", "Body 5'", "Body mid", "Body 3'",
               "TES+1kb", "TES+2kb")
RS_MARKS  <- c("H3K23ac", "H3K27me3", "H3K36me3", "H3K4me3",
               "H3K9ac", "H3K9me3", "H4K12ac", "H4K20me1")

rs_pretty <- function(f) {
  f <- as.character(f)
  paste0(sub("__.*", "", f), " @ ", RS_LABELS[match(sub(".*__", "", f), RS_KEYS)])
}

# ---- 1. Region windows in gene orientation (0 = TSS, + into the body) --------
# genes: GRanges with mcols gene_id
rs_make_regions <- function(genes) {
  neg <- as.character(strand(genes)) == "-"
  tss <- ifelse(neg, end(genes), start(genes))
  sgn <- ifelse(neg, -1, 1)
  L   <- width(genes)
  n   <- length(genes)
  bounds <- list(
    TSSm2k  = cbind(rep(-2000, n), rep(-1000, n)),
    TSSm1k  = cbind(rep(-1000, n), rep(0, n)),
    TSS200  = cbind(rep(-200, n),  rep(200, n)),
    Body5   = cbind(rep(0, n), L / 3),
    BodyMid = cbind(L / 3, 2 * L / 3),
    Body3   = cbind(2 * L / 3, L),
    TES1k   = cbind(L, L + 1000),
    TES2k   = cbind(L + 1000, L + 2000))
  out <- lapply(RS_KEYS, function(k) {
    b  <- bounds[[k]]
    a1 <- tss + sgn * b[, 1]
    a2 <- tss + sgn * b[, 2]
    GRanges(seqnames(genes), IRanges(pmin(a1, a2), pmax(a1, a2)),
            key = k, gene = as.character(genes$gene_id))
  })
  unlist(GRangesList(out))
}

# ---- 2. Signal per bin: bins x marks -----------------------------------------
# Accepts counts.rpkm.<mark>.<cond> and replicate-suffixed variants
# (.rep1, .1, ...). Replicates are averaged per bin.
rs_mark_cols <- function(all_cols, mark, cond) {
  p <- paste0("counts.rpkm.", mark, ".", cond)
  all_cols[all_cols == p | startsWith(all_cols, paste0(p, "."))]
}

rs_bin_table <- function(hmm) {
  as.data.frame(mcols(hmm$bins))
}

rs_signal_matrix <- function(bins_df, cond, marks) {
  all_cols <- colnames(bins_df)
  has      <- vapply(marks, function(m) length(rs_mark_cols(all_cols, m, cond)) > 0,
                     logical(1))
  marks <- marks[has]
  if (length(marks) == 0)
    stop("No counts.rpkm columns found for condition '", cond, "'.")
  sig <- vapply(marks, function(m) {
    v <- rowMeans(bins_df[, rs_mark_cols(all_cols, m, cond), drop = FALSE], na.rm = TRUE)
    v[is.na(v)] <- 0
    v
  }, numeric(nrow(bins_df)))
  matrix(sig, nrow = nrow(bins_df), dimnames = list(NULL, marks))
}

# ---- 3. Genes x (mark__region) matrix, overlap-weighted ----------------------
rs_region_matrix <- function(bins, sig, regions, marks) {
  ov <- findOverlaps(regions, bins)
  q  <- queryHits(ov)
  s  <- subjectHits(ov)
  w  <- width(pintersect(regions[q], bins[s]))
  contrib <- sig[s, , drop = FALSE] * (w / width(bins)[s])
  grp  <- paste(regions$gene[q], regions$key[q], sep = "|")
  agg  <- rowsum(contrib, grp)
  keys <- do.call(rbind, strsplit(rownames(agg), "|", fixed = TRUE))
  genes_all <- unique(as.character(regions$gene))
  cn <- paste(rep(marks, times = length(RS_KEYS)),
              rep(RS_KEYS, each = length(marks)), sep = "__")
  M  <- matrix(0, length(genes_all), length(cn), dimnames = list(genes_all, cn))
  gi <- match(keys[, 1], genes_all)
  for (k in RS_KEYS) {
    rows <- keys[, 2] == k
    if (any(rows))
      M[gi[rows], paste(marks, k, sep = "__")] <- agg[rows, , drop = FALSE]
  }
  M
}

# ---- 4. Expression: log2(mean normalised count over samples + 1) -------------
# ct_mat: genes x samples matrix. Matches gene_id first, then gene_name.
rs_gene_expr <- function(genes, ct_mat, samples) {
  sub         <- ct_mat[, samples, drop = FALSE]
  expr_by_row <- rowMeans(sub, na.rm = TRUE)
  ids  <- as.character(mcols(genes)$gene_id)
  nms  <- as.character(mcols(genes)$gene_name)
  idx  <- match(ids, rownames(ct_mat))
  miss <- is.na(idx)
  idx[miss] <- match(nms[miss], rownames(ct_mat))
  val <- expr_by_row[idx]
  val[is.nan(val)] <- NA
  val
}

rs_expr <- function(genes, ct, samples) {
  v <- rs_gene_expr(genes, ct$mat, samples)
  setNames(log2(v + 1), as.character(mcols(genes)$gene_id))
}

# ---- 5. One condition's feature table ----------------------------------------
rs_build_table <- function(hmm, cond, marks, regions, expr_named) {
  bins_df <- rs_bin_table(hmm)
  sig     <- rs_signal_matrix(bins_df, cond, marks)
  rm(bins_df)
  M <- rs_region_matrix(hmm$bins, sig, regions, colnames(sig))
  g <- intersect(rownames(M), names(expr_named))
  if (length(g) == 0)
    stop("No genes shared between the ChromStaR bins and the expression table. ",
         "Check that gene IDs match the count table row names.")
  d <- data.frame(M[g, , drop = FALSE], rpkm = unname(expr_named[g]),
                  check.names = FALSE)
  d[!is.na(d$rpkm), , drop = FALSE]
}

# ---- 6. Cross-validation on the training set ---------------------------------
rs_cv <- function(tr, k, ntree, seed) {
  set.seed(seed)
  folds <- caret::createFolds(tr$rpkm, k = k, list = TRUE)
  do.call(rbind, lapply(seq_along(folds), function(i) {
    va <- tr[folds[[i]], , drop = FALSE]
    t2 <- tr[-folds[[i]], , drop = FALSE]
    set.seed(seed)
    rf <- randomForest::randomForest(rpkm ~ ., data = t2, ntree = ntree)
    p  <- predict(rf, va)
    data.frame(fold = i, r = cor(p, va$rpkm),
               RMSE = sqrt(mean((p - va$rpkm)^2)),
               R2 = 1 - sum((p - va$rpkm)^2) / sum((va$rpkm - mean(va$rpkm))^2))
  }))
}

# ---- 7. Fit one condition: outliers, split, models, CV, SHAP -----------------
rs_fit <- function(d, label, cfg, progress = function(m) NULL) {
  feats <- setdiff(names(d), "rpkm")
  is_outlier <- rowSums(d[, feats, drop = FALSE]) == 0 & d$rpkm > 1
  n_out <- sum(is_outlier)
  d     <- d[!is_outlier, , drop = FALSE]
  if (nrow(d) < 50)
    stop("Only ", nrow(d), " genes left for '", label,
         "'. The model needs at least 50. Check the annotation and count table.")
  mk <- unique(sub("__.*", "", feats))
  
  set.seed(cfg$seed)
  idx <- caret::createDataPartition(d$rpkm, p = 0.8, list = FALSE)
  tr  <- d[idx, , drop = FALSE]
  te  <- d[-idx, , drop = FALSE]
  
  progress(paste("CV:", label))
  cv    <- rs_cv(tr, cfg$cv_folds, cfg$ntree, cfg$seed)
  cor_m <- cor(tr[, feats], method = "spearman")
  
  progress(paste("Combined model:", label))
  set.seed(cfg$seed)
  rf   <- randomForest::randomForest(rpkm ~ ., data = tr, ntree = cfg$ntree,
                                     importance = TRUE)
  p_te <- predict(rf, te)
  perf <- data.frame(stage = label, r_test = cor(p_te, te$rpkm),
                     RMSE_test = sqrt(mean((p_te - te$rpkm)^2)),
                     n_train = nrow(tr), n_test = nrow(te),
                     n_outliers_removed = n_out)
  
  progress(paste("Single regions:", label))
  single <- do.call(rbind, lapply(seq_along(RS_KEYS), function(i) {
    f <- paste0(mk, "__", RS_KEYS[i])
    set.seed(cfg$seed)
    r1 <- randomForest::randomForest(x = tr[, f, drop = FALSE], y = tr$rpkm,
                                     ntree = cfg$ntree_single)
    data.frame(stage = label, region = RS_LABELS[i],
               r = cor(predict(r1, te[, f, drop = FALSE]), te$rpkm))
  }))
  
  progress(paste("SHAP:", label))
  X <- te[, feats, drop = FALSE]
  if (cfg$n_shap > 0 && cfg$n_shap < nrow(X)) {
    set.seed(cfg$seed)
    X <- X[sample(nrow(X), cfg$n_shap), , drop = FALSE]
  }
  uni <- treeshap::randomForest.unify(rf, tr[, feats, drop = FALSE])
  nc  <- if (.Platform$OS.type == "windows") 1 else max(1, min(cfg$n_cores, nrow(X)))
  ch  <- split(seq_len(nrow(X)), cut(seq_len(nrow(X)), nc, labels = FALSE))
  S   <- do.call(rbind, parallel::mclapply(ch, function(ix)
    treeshap::treeshap(uni, X[ix, , drop = FALSE], verbose = FALSE)$shaps,
    mc.cores = nc))
  S    <- as.data.frame(S)
  rownames(S) <- rownames(X)
  pred <- setNames(predict(rf, X), rownames(X))
  y    <- setNames(d[rownames(X), "rpkm"], rownames(X))
  
  list(stage = label, feats = feats, marks = mk, rf = rf, cv = cv, cor = cor_m,
       perf = perf, single = single, S = S, X = X, y = y, pred = pred,
       base = mean(pred - rowSums(S)))
}

# ---- 8. Per-feature summary ---------------------------------------------------
rs_summarise <- function(fit) {
  is_expr <- fit$y >= 1
  do.call(rbind, lapply(fit$feats, function(f) {
    s <- fit$S[[f]]
    x <- fit$X[[f]]
    data.frame(stage = fit$stage, feature = f,
               mark = sub("__.*", "", f),
               region = RS_LABELS[match(sub(".*__", "", f), RS_KEYS)],
               mean_abs_SHAP = mean(abs(s)), mean_SHAP = mean(s),
               direction_rho = if (sd(x) > 0 && sd(s) > 0)
                 cor(x, s, method = "spearman") else NA_real_,
               mean_SHAP_expressed = mean(s[is_expr]),
               mean_SHAP_silent = mean(s[!is_expr]),
               pct_genes_with_mark = 100 * mean(x > 0),
               stringsAsFactors = FALSE)
  }))
}

# ---- 9. Whole run -------------------------------------------------------------
# specs: list of list(label, hmm, cond, marks, ct, samples)
rs_run <- function(specs, genes, cfg, progress = function(m) NULL) {
  regs <- rs_make_regions(genes)
  fits <- list()
  for (sp in specs) {
    progress(paste("Features:", sp$label))
    expr <- rs_expr(genes, sp$ct, sp$samples)
    d    <- rs_build_table(sp$hmm, sp$cond, sp$marks, regs, expr)
    fits[[sp$label]] <- rs_fit(d, sp$label, cfg, progress)
    gc()
  }
  fs <- do.call(rbind, lapply(fits, rs_summarise))
  fs$region <- factor(fs$region, levels = RS_LABELS)
  fs$mark   <- factor(fs$mark, levels = rev(RS_MARKS))
  list(
    fits    = fits,
    perf    = do.call(rbind, lapply(fits, `[[`, "perf")),
    cv      = do.call(rbind, lapply(names(fits), function(s)
      cbind(stage = s, fits[[s]]$cv))),
    single  = do.call(rbind, lapply(fits, `[[`, "single")),
    fs      = fs,
    cor_all = do.call(rbind, lapply(names(fits), function(s) {
      d <- as.data.frame(as.table(fits[[s]]$cor), stringsAsFactors = FALSE)
      names(d) <- c("f1", "f2", "rho")
      cbind(stage = s, d)
    })),
    cfg     = cfg)
}

# ---- Figures ------------------------------------------------------------------
rs_plot_corr <- function(cor_all) {
  lv <- unique(cor_all$f1)
  cor_all$f1 <- factor(cor_all$f1, levels = lv)
  cor_all$f2 <- factor(cor_all$f2, levels = rev(lv))
  ggplot(cor_all, aes(f1, f2, fill = rho)) + geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-1, 1), name = "Spearman rho") +
    scale_x_discrete(labels = rs_pretty) + scale_y_discrete(labels = rs_pretty) +
    facet_wrap(~ stage) +
    labs(x = NULL, y = NULL,
         title = "Correlation between the 64 mark @ region features") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

rs_plot_single <- function(single, perf) {
  all <- rbind(single,
               data.frame(stage = perf$stage, region = "All regions", r = perf$r_test,
                          stringsAsFactors = FALSE))
  all$region <- factor(all$region, levels = c(RS_LABELS, "All regions"))
  ggplot(all, aes(region, r, fill = region == "All regions")) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", r)), vjust = -0.3, size = 3.2) +
    facet_wrap(~ stage) +
    scale_fill_manual(values = c(`FALSE` = "grey65", `TRUE` = "#b2182b"), guide = "none") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Pearson r (test genes)",
         title = "Each region alone vs all regions combined") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

rs_plot_bubble <- function(fs) {
  ggplot(fs, aes(region, mark)) +
    geom_point(aes(size = mean_abs_SHAP, fill = direction_rho),
               shape = 21, colour = "grey30") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-1, 1), na.value = "grey85",
                         name = "Direction\n(rho)") +
    scale_size_area(max_size = 12, name = "Importance\n(mean |SHAP|)") +
    facet_wrap(~ stage) +
    labs(x = NULL, y = NULL, title = "Which mark matters where along the gene") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(colour = "grey92"))
}

rs_plot_heat <- function(fs) {
  hg <- rbind(
    data.frame(stage = fs$stage, mark = fs$mark, region = fs$region,
               group = "Expressed (FPKM >= 1)", mean_SHAP = fs$mean_SHAP_expressed,
               stringsAsFactors = FALSE),
    data.frame(stage = fs$stage, mark = fs$mark, region = fs$region,
               group = "Silent (FPKM < 1)", mean_SHAP = fs$mean_SHAP_silent,
               stringsAsFactors = FALSE))
  lim <- max(abs(hg$mean_SHAP), na.rm = TRUE)
  ggplot(hg, aes(region, mark, fill = mean_SHAP)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", mean_SHAP)), size = 2.6) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-lim, lim), name = "Mean SHAP") +
    facet_grid(group ~ stage) +
    labs(x = NULL, y = NULL, title = "Average push of each mark at each position") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

rs_plot_beeswarm <- function(fit, top = 20) {
  imp <- sort(colMeans(abs(fit$S)), decreasing = TRUE)[seq_len(min(top, ncol(fit$S)))]
  d <- data.frame(
    feature = rep(names(imp), each = nrow(fit$S)),
    shap    = unlist(fit$S[, names(imp), drop = FALSE], use.names = FALSE),
    present = unlist(fit$X[, names(imp), drop = FALSE], use.names = FALSE) > 0,
    stringsAsFactors = FALSE)
  d$feature <- factor(rs_pretty(d$feature), levels = rev(rs_pretty(names(imp))))
  ggplot(d, aes(shap, feature, colour = present)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    ggbeeswarm::geom_quasirandom(orientation = "y", size = 0.4, alpha = 0.5,
                                 width = 0.35) +
    scale_colour_manual(values = c(`FALSE` = "#2c7bb6", `TRUE` = "#d7191c"),
                        labels = c("absent", "present"), name = "Mark") +
    labs(x = "SHAP value (effect on predicted log2 expression)", y = NULL,
         title = paste0("Top ", top, " features - ", fit$stage)) +
    theme_minimal(base_size = 11)
}

rs_plot_waterfall <- function(fit, gene, top = 12) {
  v    <- unlist(fit$S[gene, ])
  o    <- order(abs(v), decreasing = TRUE)
  keep <- v[o[seq_len(top)]]
  vals <- c(keep, sum(v[o[-seq_len(top)]]))
  labs <- c(rs_pretty(names(keep)), "other features")
  end  <- fit$base + cumsum(vals)
  start <- end - vals
  d <- data.frame(label = factor(labs, levels = rev(labs)), start, end,
                  dir = ifelse(vals >= 0, "up", "down"))
  ggplot(d) +
    geom_segment(aes(x = start, xend = end, y = label, yend = label, colour = dir),
                 linewidth = 6) +
    geom_vline(xintercept = fit$base, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = fit$pred[[gene]], colour = "black") +
    scale_colour_manual(values = c(up = "#d7191c", down = "#2c7bb6"), guide = "none") +
    labs(x = "Predicted log2(FPKM+1)", y = NULL,
         title = paste(gene, "-", fit$stage),
         subtitle = sprintf("baseline %.2f -> predicted %.2f (measured %.2f)",
                            fit$base, fit$pred[[gene]], fit$y[[gene]])) +
    theme_minimal(base_size = 11)
}

# ---- Output files -------------------------------------------------------------
rs_write_outputs <- function(res, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  wc <- function(x, f) write.csv(x, file.path(dir, f), row.names = FALSE)
  wc(res$perf,    "performance_combined_model.csv")
  wc(res$cv,      "cv_5fold_results.csv")
  wc(res$single,  "single_region_performance.csv")
  wc(res$fs,      "shap_feature_summary.csv")
  wc(res$cor_all, "feature_correlations.csv")
  saveRDS(res$fits, file.path(dir, "rf_shap_fits.rds"))
  ggsave(file.path(dir, "Fig_feature_correlation.pdf"), rs_plot_corr(res$cor_all),
         width = 11, height = 9)
  ggsave(file.path(dir, "Fig_single_vs_combined.pdf"),
         rs_plot_single(res$single, res$perf), width = 9, height = 5)
  ggsave(file.path(dir, "Fig_mark_by_region_bubble.pdf"), rs_plot_bubble(res$fs),
         width = 11, height = 5.5)
  ggsave(file.path(dir, "Fig_expressed_vs_silent_heatmap.pdf"), rs_plot_heat(res$fs),
         width = 11, height = 7)
  for (s in names(res$fits))
    ggsave(file.path(dir, paste0("Fig_beeswarm_", make.names(s), ".pdf")),
           rs_plot_beeswarm(res$fits[[s]]), width = 9, height = 7)
  invisible(dir)
}

# ---- Condition filter used by the server --------------------------------------
# Keeps only conditions the user has not removed from the Load Data order box.
# Falls back to all of the stage's conditions if nothing overlaps.
rs_kept_conditions <- function(st, kept_order) {
  kept <- intersect(kept_order, st$conditions)
  if (length(kept) == 0) st$conditions else kept
}