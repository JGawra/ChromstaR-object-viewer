# =============================================================================
# ChromstaR Histone Modification Viewer — Shiny App
# Visualizes 8 histone marks across 2 conditions (e.g., Control vs Treatment)
# Supports:
#   - Metagene profiles (around TSS / gene body)
#   - Genomic region browser (zoom in/out)
#   - Side-by-side or overlay condition comparison
# =============================================================================
# Dependencies (run once to install):
#   install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr",
#                      "tidyr", "scales", "shinyWidgets", "DT"))
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("chromstaR", "GenomicRanges", "rtracklayer",
#                          "GenomicFeatures", "Gviz"))
# =============================================================================

library(shiny)
library(shinydashboard)

# Allow uploads up to 500 MB (adjust as needed)
options(shiny.maxRequestSize = 500 * 1024^2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(shinyWidgets)
library(DT)
library(magrittr)   # for %>% pipe in DT formatting
library(openxlsx)   # for Excel export of enrichment data

# Bioconductor
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)

# Optional — used for region browser tracks
# library(Gviz)  # loaded conditionally below

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load a ChromstaR combined multivariate model from .RData
load_chromstar_object <- function(path) {
  env <- new.env()
  load(path, envir = env)
  objs <- ls(env)
  # Return the first ChromstaR-related object found
  for (nm in objs) {
    obj <- get(nm, envir = env)
    if (inherits(obj, c("combinedMultiHMM", "multiHMM", "uniHMM"))) {
      return(obj)
    }
  }
  # Fallback: return the first object
  return(get(objs[1], envir = env))
}

#' Expand ChromstaR bins GRanges to a flat data.frame
#' counts.rpkm and posteriors are list-columns; this unpacks them
#' exactly as as.data.frame(combined.model$bins) does in the Rmd.
expand_bins_df <- function(hmm) {
  bins <- hmm$bins

  # as.data.frame() on the full GRanges expands list-columns correctly
  df <- as.data.frame(bins)  # seqnames, start, end, width, strand + all mcols

  # Rename seqnames -> chr for consistency
  if ("seqnames" %in% colnames(df)) {
    df$chr <- as.character(df$seqnames)
  }
  df
}

#' Extract per-mark, per-condition signal across a genomic region
#' Uses counts.rpkm.{MARK}.{CONDITION}.{replicate} columns from hmm$bins
extract_signal_region <- function(hmm, region_gr, marks, conditions,
                                  genes_gr = NULL, bin_scope = "all") {
  bins    <- hmm$bins
  bins_df <- expand_bins_df(hmm)

  # Find bins overlapping the region
  ov      <- findOverlaps(bins, region_gr)
  if (length(ov) == 0) return(NULL)
  idx     <- queryHits(ov)

  # Optionally restrict by gene overlap: "genic" keeps only bins inside a
  # gene, "intergenic" keeps only bins outside every gene, "all" keeps both
  if (bin_scope %in% c("genic", "intergenic") && !is.null(genes_gr) && length(genes_gr) > 0) {
    gov <- findOverlaps(bins[idx], genes_gr)
    genic_idx_pos <- unique(queryHits(gov))  # positions within idx that are genic

    if (bin_scope == "genic") {
      if (length(genic_idx_pos) == 0) return(NULL)
      idx <- idx[genic_idx_pos]
    } else {  # intergenic
      intergenic_idx_pos <- setdiff(seq_along(idx), genic_idx_pos)
      if (length(intergenic_idx_pos) == 0) return(NULL)
      idx <- idx[intergenic_idx_pos]
    }
  }

  mid_pos <- (bins_df$start[idx] + bins_df$end[idx]) / 2
  chr_vec <- bins_df$chr[idx]

  # Check available rpkm columns once
  rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
  if (length(rpkm_cols_all) == 0) {
    showNotification("No counts.rpkm columns found in bins. Check object structure.", type = "error")
    return(NULL)
  }

  df_list <- list()
  for (cond in conditions) {
    for (mark in marks) {
      pattern <- paste0("counts.rpkm.", mark, ".", cond)
      cols    <- grep(pattern, rpkm_cols_all, value = TRUE, fixed = TRUE)
      if (length(cols) == 0) next

      sig <- rowMeans(bins_df[idx, cols, drop = FALSE], na.rm = TRUE)

      df_list[[paste(cond, mark)]] <- data.frame(
        position  = mid_pos,
        chr       = chr_vec,
        signal    = sig,
        mark      = mark,
        condition = cond,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(df_list) == 0) return(NULL)
  bind_rows(df_list)
}

#' Compute log(observed/expected) enrichment profile around gene boundaries
#' Mirrors the reference figure: x-axis spans upstream-of-TSS -> gene body (%) -> downstream-of-TES
#' "Observed" = mean RPKM in each bin; "Expected" = genome-wide mean RPKM for that mark/condition
compute_enrichment_profile <- function(hmm, genes_gr, marks, conditions,
                                       upstream   = 2000,
                                       downstream = 2000,
                                       n_bins     = 40) {

  bins    <- hmm$bins
  bins_df <- expand_bins_df(hmm)
  rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
  if (length(rpkm_cols_all) == 0) return(NULL)
  bin_mid_all <- (bins_df$start + bins_df$end) / 2

  # Genome-wide "expected" baseline per mark/condition (mean over all bins) —
  # precompute the signal vector for every mark/condition combo once
  expected   <- list()
  signal_mat <- list()  # key "cond mark" -> numeric vector (length = n bins)
  for (cond in conditions) {
    for (mark in marks) {
      pattern <- paste0("counts.rpkm.", mark, ".", cond)
      cols    <- grep(pattern, rpkm_cols_all, value = TRUE, fixed = TRUE)
      if (length(cols) == 0) next
      vals <- rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)
      key  <- paste(cond, mark)
      signal_mat[[key]] <- vals
      expected[[key]]   <- mean(vals, na.rm = TRUE)
    }
  }

  # ---- Vectorised zone construction across ALL genes at once ----------------
  s_g   <- start(genes_gr)
  e_g   <- end(genes_gr)
  str_g <- as.character(strand(genes_gr))
  chr_g <- as.character(seqnames(genes_gr))
  tss   <- ifelse(str_g == "-", e_g, s_g)
  tes   <- ifelse(str_g == "-", s_g, e_g)
  is_minus <- str_g == "-"

  up_start   <- ifelse(is_minus, tss,               pmax(1, tss - upstream))
  up_end     <- ifelse(is_minus, tss + upstream,     tss)
  down_start <- ifelse(is_minus, pmax(1, tes - downstream), tes)
  down_end   <- ifelse(is_minus, tes,                tes + downstream)

  gene_idx <- seq_along(genes_gr)

  build_zone_gr <- function(starts, ends, chrs, idx) {
    ok <- starts >= 1 & ends >= 1 & starts <= ends
    GRanges(
      seqnames = chrs[ok],
      ranges   = IRanges(starts[ok], ends[ok]),
      gene_idx = idx[ok]
    )
  }

  zone_gr <- list(
    upstream   = build_zone_gr(up_start,   up_end,   chr_g, gene_idx),
    body       = build_zone_gr(s_g,        e_g,      chr_g, gene_idx),
    downstream = build_zone_gr(down_start, down_end, chr_g, gene_idx)
  )

  # ---- 3 findOverlaps calls total (not 3 * n_genes) --------------------------
  all_profiles <- list()

  for (zname in names(zone_gr)) {
    zr <- zone_gr[[zname]]
    if (length(zr) == 0) next
    ov <- findOverlaps(bins, zr)
    if (length(ov) == 0) next

    bi <- queryHits(ov)    # bin row index
    zi <- subjectHits(ov)  # row index within zr (maps to a gene via gene_idx)
    g_i <- mcols(zr)$gene_idx[zi]   # gene index for each overlap row

    bin_mid  <- bin_mid_all[bi]
    zone_s   <- start(zr)[zi]
    zone_w   <- width(zr)[zi]
    minus    <- is_minus[g_i]

    if (zname == "upstream") {
      frac  <- (bin_mid - zone_s) / zone_w
      frac  <- ifelse(minus, 1 - frac, frac)
      x_val <- -1 + frac
    } else if (zname == "body") {
      glen  <- (e_g[g_i] - s_g[g_i])
      glen[glen == 0] <- 1
      frac  <- (bin_mid - s_g[g_i]) / glen
      frac  <- ifelse(minus, 1 - frac, frac)
      x_val <- frac
    } else {
      frac  <- (bin_mid - zone_s) / zone_w
      frac  <- ifelse(minus, 1 - frac, frac)
      x_val <- 1 + frac
    }

    slot_idx <- round(x_val * n_bins)

    for (cond in conditions) {
      for (mark in marks) {
        key <- paste(cond, mark)
        sig_vec <- signal_mat[[key]]
        if (is.null(sig_vec)) next
        exp_val <- expected[[key]]
        if (is.na(exp_val) || exp_val <= 0) next

        sig <- sig_vec[bi]
        sig[sig <= 0] <- NA  # avoid log(0)

        prof_key <- paste(cond, mark, zname)
        all_profiles[[prof_key]] <- data.frame(
          slot      = slot_idx,
          log_ratio = log(sig / exp_val),
          mark      = mark,
          condition = cond,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(all_profiles) == 0) return(NULL)

  bind_rows(all_profiles) %>%
    group_by(slot, mark, condition) %>%
    summarise(
      mean_log_ratio = mean(log_ratio, na.rm = TRUE),
      n_genes        = n(),
      .groups        = "drop"
    ) %>%
    filter(is.finite(mean_log_ratio))
}

#' Compute metagene profile around TSS or gene body
#' Mirrors the PNAgene_profile() logic from the Rmd analysis
compute_metagene <- function(hmm, genes_gr, marks, conditions,
                             mode       = "TSS",
                             upstream   = 2000,
                             downstream = 2000,
                             n_bins     = 100) {

  bins    <- hmm$bins
  bins_df <- expand_bins_df(hmm)
  bin_mid_all <- (bins_df$start + bins_df$end) / 2

  rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
  if (length(rpkm_cols_all) == 0) return(NULL)

  # Precompute the signal vector for every mark/condition combo once
  signal_mat <- list()
  for (cond in conditions) {
    for (mark in marks) {
      pattern <- paste0("counts.rpkm.", mark, ".", cond)
      cols    <- grep(pattern, rpkm_cols_all, value = TRUE, fixed = TRUE)
      if (length(cols) == 0) next
      signal_mat[[paste(cond, mark)]] <- rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)
    }
  }

  # Build all reference windows (one per gene) at once
  if (mode == "TSS") {
    ref_points <- promoters(genes_gr, upstream = upstream, downstream = downstream)
  } else if (mode == "TES") {
    # TES-centered window: mirror of promoters() but anchored at the gene's
    # 3' end, respecting strand (TES = end for "+", start for "-")
    str_g <- as.character(strand(genes_gr))
    tes   <- ifelse(str_g == "-", start(genes_gr), end(genes_gr))
    win_start <- ifelse(str_g == "-", tes - downstream, tes - upstream)
    win_end   <- ifelse(str_g == "-", tes + upstream,   tes + downstream)
    win_start <- pmax(1, win_start)
    ref_points <- GRanges(
      seqnames = seqnames(genes_gr),
      ranges   = IRanges(win_start, win_end),
      strand   = strand(genes_gr)
    )
  } else {
    ref_points <- genes_gr
  }
  mcols(ref_points)$gene_idx <- seq_along(ref_points)

  # ---- ONE findOverlaps call for all genes, not one per gene ----------------
  ov <- findOverlaps(bins, ref_points)
  if (length(ov) == 0) return(NULL)

  bi   <- queryHits(ov)
  gi   <- subjectHits(ov)
  g_idx <- mcols(ref_points)$gene_idx[gi]

  region_start <- start(ref_points)[gi]
  region_width <- width(ref_points)[gi]
  bin_mid      <- bin_mid_all[bi]

  rel_pos  <- (bin_mid - region_start) / region_width
  slot_idx <- pmin(pmax(ceiling(rel_pos * n_bins), 1L), n_bins)

  all_profiles <- list()
  for (cond in conditions) {
    for (mark in marks) {
      key <- paste(cond, mark)
      sig_vec <- signal_mat[[key]]
      if (is.null(sig_vec)) next

      sig <- sig_vec[bi]
      if (all(is.na(sig))) next

      all_profiles[[key]] <- data.frame(
        bin_idx   = slot_idx,
        signal    = sig,
        mark      = mark,
        condition = cond,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(all_profiles) == 0) return(NULL)

  bind_rows(all_profiles) %>%
    group_by(bin_idx, mark, condition) %>%
    summarise(
      mean_signal = mean(signal, na.rm = TRUE),
      se_signal   = sd(signal,   na.rm = TRUE) / sqrt(n()),
      .groups     = "drop"
    )
}

#' Compare two gene sets by their per-gene promoter-averaged posterior
#' probabilities, one value per mark/condition combination, then run a
#' Wilcoxon test per mark/condition to see which marks differ most between
#' the two gene sets (e.g. a curated list vs random genes, or up- vs
#' down-regulated genes). Mirrors the user's two-script pipeline:
#' (1) extract posteriors.* columns from hmm$bins, (2) average over each
#' gene's promoter window and compare gene set A vs gene set B.
#' Collapse replicate-level posterior columns (e.g. posteriors.H3K9me3.PA.rep1,
#' posteriors.H3K9me3.PA.rep2) into a single averaged column per mark/condition
#' (posteriors.H3K9me3.PA), the same way counts.rpkm replicates are averaged
#' elsewhere in the app. Returns a list with the collapsed data.frame (only the
#' new averaged columns, same row count) and the new column name vector.
collapse_posterior_replicates <- function(bins_df, post_cols, average_replicates = TRUE) {
  if (!average_replicates) {
    # Keep each replicate as its own column/row in downstream results
    return(list(df = bins_df[, post_cols, drop = FALSE], cols = post_cols))
  }

  # Strip a trailing ".rep<digits>" or ".<digits>" replicate suffix
  base_name <- sub("\\.(rep)?[0-9]+$", "", post_cols)

  if (all(base_name == post_cols)) {
    # No replicate suffix detected — already one column per mark/condition
    return(list(df = bins_df[, post_cols, drop = FALSE], cols = post_cols))
  }

  unique_bases <- unique(base_name)
  out <- lapply(unique_bases, function(b) {
    cols <- post_cols[base_name == b]
    rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)
  })
  out_df <- as.data.frame(out)
  colnames(out_df) <- unique_bases
  list(df = out_df, cols = unique_bases)
}

compute_gene_set_posteriors <- function(hmm, genes_A_gr, genes_B_gr,
                                        upstream = 2000, downstream = 2000,
                                        average_replicates = TRUE,
                                        summary_stat = "mean") {

  stat_fn <- if (identical(summary_stat, "median")) {
    function(x) median(x, na.rm = TRUE)
  } else {
    function(x) mean(x, na.rm = TRUE)
  }

  bins    <- hmm$bins
  bins_df <- expand_bins_df(hmm)

  post_cols_raw <- grep("^posteriors\\.", colnames(bins_df), value = TRUE)
  if (length(post_cols_raw) == 0) {
    return(list(error = "No columns starting with 'posteriors.' were found on hmm$bins. This object may not have posterior probabilities stored."))
  }
  collapsed <- collapse_posterior_replicates(bins_df, post_cols_raw, average_replicates)
  bins_df   <- cbind(bins_df, collapsed$df)
  post_cols <- collapsed$cols

  # Build promoter windows for both gene sets in one go, tagging group + gene
  build_promoters <- function(gr, group_label) {
    if (length(gr) == 0) return(NULL)
    p <- promoters(gr, upstream = upstream, downstream = downstream)
    mcols(p)$gene_label <- mcols(gr)$gene_label
    mcols(p)$group <- group_label
    p
  }

  prom_A <- build_promoters(genes_A_gr, "A")
  prom_B <- build_promoters(genes_B_gr, "B")

  regions <- c(prom_A, prom_B)
  if (is.null(regions) || length(regions) == 0) {
    return(list(error = "No promoter regions could be built for either gene set."))
  }

  # ---- ONE findOverlaps call for both gene sets at once ----------------------
  ov <- findOverlaps(regions, bins)
  if (length(ov) == 0) {
    return(list(error = "No overlaps found between gene promoters and ChromstaR bins. Check that gene and bin chromosome names match."))
  }

  ri <- queryHits(ov)    # region index (promoter row)
  bi <- subjectHits(ov)  # bin index

  overlap_df <- data.frame(
    gene  = mcols(regions)$gene_label[ri],
    group = mcols(regions)$group[ri],
    bins_df[bi, post_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  # Average posterior per gene/group across all overlapping bins
  gene_scores <- overlap_df %>%
    group_by(gene, group) %>%
    summarise(across(all_of(post_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  if (nrow(gene_scores) == 0) {
    return(list(error = "No gene-level posterior scores could be computed."))
  }

  long_scores <- gene_scores %>%
    pivot_longer(cols = all_of(post_cols), names_to = "mark_col", values_to = "posterior")

  # Parse mark_col like "posteriors.H3K4me3.PA" into mark + condition
  long_scores <- long_scores %>%
    mutate(
      mark_col_clean = sub("^posteriors\\.", "", mark_col),
      mark      = sub("\\..*$", "", mark_col_clean),
      condition = sub("^[^.]*\\.", "", mark_col_clean)
    )

  safe_wilcox <- function(x, y) {
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    if (length(x) < 1 || length(y) < 1) return(NA_real_)
    tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  }

  stats <- long_scores %>%
    group_by(mark, condition) %>%
    summarise(
      n_A      = sum(group == "A" & !is.na(posterior)),
      n_B      = sum(group == "B" & !is.na(posterior)),
      mean_A   = mean(posterior[group == "A"], na.rm = TRUE),
      mean_B   = mean(posterior[group == "B"], na.rm = TRUE),
      median_A = median(posterior[group == "A"], na.rm = TRUE),
      median_B = median(posterior[group == "B"], na.rm = TRUE),
      stat_A   = stat_fn(posterior[group == "A"]),
      stat_B   = stat_fn(posterior[group == "B"]),
      delta_A_minus_B = stat_A - stat_B,
      p_value  = safe_wilcox(posterior[group == "A"], posterior[group == "B"]),
      .groups  = "drop"
    ) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR)
  attr(stats, "summary_stat") <- summary_stat

  list(
    error       = NULL,
    long_scores = long_scores,
    gene_scores = gene_scores,
    stats       = stats,
    n_genes_A   = length(unique(long_scores$gene[long_scores$group == "A"])),
    n_genes_B   = length(unique(long_scores$gene[long_scores$group == "B"]))
  )
}

#' Precompute promoter-averaged posterior probabilities for every gene in
#' the annotation (used as the universe to draw random background sets
#' from). This is the expensive part (one findOverlaps + groupby over all
#' genes) and is done exactly once regardless of how many permutations follow.
compute_all_gene_promoter_posteriors <- function(hmm, all_genes_gr, gene_labels, upstream = 2000, downstream = 2000,
                                                 average_replicates = TRUE) {

  bins    <- hmm$bins
  bins_df <- expand_bins_df(hmm)

  post_cols_raw <- grep("^posteriors\\.", colnames(bins_df), value = TRUE)
  if (length(post_cols_raw) == 0) {
    return(list(error = "No columns starting with 'posteriors.' were found on hmm$bins."))
  }
  collapsed <- collapse_posterior_replicates(bins_df, post_cols_raw, average_replicates)
  bins_df   <- cbind(bins_df, collapsed$df)
  post_cols <- collapsed$cols

  if (length(gene_labels) != length(all_genes_gr)) {
    return(list(error = paste0("Internal error: gene_labels length (", length(gene_labels),
                               ") does not match all_genes_gr length (", length(all_genes_gr), ").")))
  }

  prom <- promoters(all_genes_gr, upstream = upstream, downstream = downstream)
  if (length(prom) != length(gene_labels)) {
    return(list(error = paste0("Internal error: promoters() returned ", length(prom),
                               " ranges but gene_labels has ", length(gene_labels),
                               " entries — gene count changed unexpectedly.")))
  }
  mcols(prom)$gene_label <- gene_labels

  ov <- findOverlaps(prom, bins)
  if (length(ov) == 0) {
    return(list(error = "No overlaps found between gene promoters and ChromstaR bins."))
  }

  ri <- queryHits(ov)
  bi <- subjectHits(ov)

  overlap_df <- data.frame(
    gene = mcols(prom)$gene_label[ri],
    bins_df[bi, post_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  gene_scores <- overlap_df %>%
    group_by(gene) %>%
    summarise(across(all_of(post_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  n_dup <- sum(duplicated(gene_scores$gene))
  if (n_dup > 0) {
    return(list(error = paste0("Internal error: ", n_dup, " duplicate gene names found after grouping — gene labels are not unique.")))
  }

  list(error = NULL, gene_scores = gene_scores, post_cols = post_cols)
}

#' Permutation test: compare gene set A's mean promoter posterior (per
#' mark/condition) against the distribution of means from n_perm random
#' draws of the same size from the full gene universe. Returns an empirical
#' p-value (two-sided) per mark/condition, much more robust than comparing
#' against a single arbitrary random set.
compute_gene_set_permutation <- function(all_gene_scores, post_cols, genes_A_labels,
                                          n_perm = 1000, seed = NULL,
                                          summary_stat = "mean") {

  if (!is.null(seed)) set.seed(seed)

  use_median <- identical(summary_stat, "median")

  gene_scores <- all_gene_scores
  n_total <- nrow(gene_scores)
  n_A     <- sum(gene_scores$gene %in% genes_A_labels)

  if (n_A == 0) {
    return(list(error = "None of gene set A's genes were found in the precomputed universe — check gene name matching."))
  }
  if (n_A >= n_total) {
    return(list(error = "Gene set A is the same size as (or larger than) the full gene universe — no room to draw a smaller random comparison set."))
  }

  is_A <- gene_scores$gene %in% genes_A_labels
  mat  <- as.matrix(gene_scores[, post_cols, drop = FALSE])  # n_total x n_marks

  obs_mean_A <- if (use_median) {
    apply(mat[is_A, , drop = FALSE], 2, median, na.rm = TRUE)
  } else {
    colMeans(mat[is_A, , drop = FALSE], na.rm = TRUE)
  }

  # ---- Permutation: draw n_perm random index sets of size n_A ---------------
  # colMeans() is vectorised and fast; median has no base-R vectorised
  # equivalent, so the median path costs more per draw (still fine up to
  # several thousand permutations).
  perm_means <- matrix(NA_real_, nrow = n_perm, ncol = ncol(mat))
  colnames(perm_means) <- post_cols

  for (i in seq_len(n_perm)) {
    idx <- sample.int(n_total, n_A)
    perm_means[i, ] <- if (use_median) {
      apply(mat[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    } else {
      colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
    }
  }

  perm_mean_of_means <- colMeans(perm_means, na.rm = TRUE)
  perm_sd            <- apply(perm_means, 2, sd, na.rm = TRUE)

  # Empirical two-sided p-value: fraction of permutations at least as extreme
  # as the observed deviation from the permutation mean
  emp_p <- vapply(seq_along(post_cols), function(j) {
    obs_dev  <- abs(obs_mean_A[j] - perm_mean_of_means[j])
    perm_dev <- abs(perm_means[, j] - perm_mean_of_means[j])
    (sum(perm_dev >= obs_dev, na.rm = TRUE) + 1) / (n_perm + 1)
  }, numeric(1))

  mark_col_clean <- sub("^posteriors\\.", "", post_cols)
  mark_vec      <- sub("\\..*$", "", mark_col_clean)
  condition_vec <- sub("^[^.]*\\.", "", mark_col_clean)

  stats <- data.frame(
    mark            = mark_vec,
    condition       = condition_vec,
    n_A             = n_A,
    n_random_pool   = n_total,
    n_perm          = n_perm,
    stat_A          = round(obs_mean_A, 4),
    stat_random     = round(perm_mean_of_means, 4),
    sd_random       = round(perm_sd, 4),
    delta_A_minus_random = round(obs_mean_A - perm_mean_of_means, 4),
    z_score         = round((obs_mean_A - perm_mean_of_means) / perm_sd, 3),
    empirical_p     = emp_p,
    stringsAsFactors = FALSE
  )
  stats$FDR <- p.adjust(stats$empirical_p, method = "BH")
  stats <- stats[order(stats$FDR), ]
  attr(stats, "summary_stat") <- summary_stat

  # Long-format for plotting: one random "background" mean draw per
  # permutation (already aggregated) plus the single observed A value
  perm_long <- as.data.frame(perm_means)
  perm_long$perm_id <- seq_len(n_perm)
  perm_long <- perm_long %>%
    pivot_longer(cols = all_of(post_cols), names_to = "mark_col", values_to = "perm_mean") %>%
    mutate(
      mark_col_clean = sub("^posteriors\\.", "", mark_col),
      mark      = sub("\\..*$", "", mark_col_clean),
      condition = sub("^[^.]*\\.", "", mark_col_clean)
    )

  obs_long <- data.frame(
    mark_col  = post_cols,
    obs_mean  = obs_mean_A,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      mark_col_clean = sub("^posteriors\\.", "", mark_col),
      mark      = sub("\\..*$", "", mark_col_clean),
      condition = sub("^[^.]*\\.", "", mark_col_clean)
    )

  # ---- One representative random draw, kept at gene-level resolution -------
  # for an apples-to-apples boxplot against gene set A (same visual style as
  # the manual two-list comparison), separate from the 1000-draw permutation
  # pool used for the empirical p-value above.
  rep_idx <- sample.int(n_total, n_A)
  rep_genes <- gene_scores$gene[rep_idx]

  long_A <- gene_scores[is_A, c("gene", post_cols)] %>%
    pivot_longer(cols = all_of(post_cols), names_to = "mark_col", values_to = "posterior") %>%
    mutate(group = "A")
  long_rand <- gene_scores[rep_idx, c("gene", post_cols)] %>%
    pivot_longer(cols = all_of(post_cols), names_to = "mark_col", values_to = "posterior") %>%
    mutate(group = "Random")

  long_scores <- bind_rows(long_A, long_rand) %>%
    mutate(
      mark_col_clean = sub("^posteriors\\.", "", mark_col),
      mark      = sub("\\..*$", "", mark_col_clean),
      condition = sub("^[^.]*\\.", "", mark_col_clean)
    )

  list(error = NULL, stats = stats, perm_long = perm_long, obs_long = obs_long,
       long_scores = long_scores, n_A = n_A, n_total = n_total)
}
parse_marks_conditions <- function(hmm) {

  marks      <- character(0)
  conditions <- character(0)

  # --- Try hmm$info (most reliable in combinedMultiHMM) ---
  if (!is.null(hmm$info)) {
    info <- hmm$info
    if (!is.null(info$mark))      marks      <- unique(as.character(info$mark))
    if (!is.null(info$condition)) conditions <- unique(as.character(info$condition))
    if (!is.null(info$sample.name) && length(marks) == 0) {
      # Try to split sample names into mark + condition
      parts <- strsplit(as.character(info$sample.name), "[._-]")
      marks      <- unique(sapply(parts, `[`, 1))
      conditions <- unique(sapply(parts, function(x) paste(x[-1], collapse = ".")))
      conditions <- conditions[nchar(conditions) > 0]
    }
  }

  # --- Try hmm$hmms list (named by "mark-condition" or "condition-mark") ---
  if (length(marks) == 0 && !is.null(hmm$hmms)) {
    nms   <- names(hmm$hmms)
    parts <- strsplit(nms, "[._-]")
    # Heuristic: last token = condition if it matches known pattern
    all_last  <- sapply(parts, function(x) x[length(x)])
    all_first <- sapply(parts, `[`, 1)
    conditions <- unique(all_last)
    marks      <- unique(all_first)
  }

  # --- Fallback: parse bin column names like "combination.PA" ---------------
  if (length(conditions) == 0) {
    avail <- colnames(mcols(hmm$bins))
    # Look for columns that end with a condition token after "."
    # e.g. combination.PA, combination.PNA -> conditions PA, PNA
    cond_cols  <- avail[grepl("^combination\\.", avail)]
    if (length(cond_cols) > 0) {
      conditions <- sub("^combination\\.", "", cond_cols)
    } else {
      # Generic: split on "." and take last part
      parts      <- strsplit(avail, ".", fixed = TRUE)
      last_parts <- sapply(parts, function(x) x[length(x)])
      # Keep only tokens that appear more than once (likely conditions)
      tab        <- table(last_parts)
      conditions <- names(tab[tab > 1])
    }
  }

  # --- Marks from posteriors columns or hmm$segments ------------------------
  if (length(marks) == 0 && !is.null(hmm$segments)) {
    seg_cols <- colnames(mcols(hmm$segments))
    marks <- seg_cols[!seg_cols %in% c("combination", "state", "transition.group",
                                        "differential.score", "maxPostInPeak")]
  }

  # Last resort: show all bin column names so user knows what's available
  if (length(marks) == 0) {
    avail  <- colnames(mcols(hmm$bins))
    marks  <- avail  # will be filtered later
  }

  list(marks = marks, conditions = conditions,
       all_cols = colnames(mcols(hmm$bins)))
}

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "ChromstaR Viewer",
    tags$li(
      class = "dropdown",
      style = "padding: 8px 15px;",
      actionButton("show_help", label = "Help", icon = icon("circle-question"),
                  class = "btn-sm btn-info")
    )
  ),

  dashboardSidebar(
    sidebarMenu(
      id = "sidebar_tabs",
      menuItem("Load Data",       tabName = "load",       icon = icon("folder-open")),
      menuItem("Metagene Profile",tabName = "metagene",   icon = icon("chart-area")),
      menuItem("Enrichment Profile", tabName = "enrichment", icon = icon("chart-line")),
      menuItem("Region Browser",  tabName = "browser",    icon = icon("dna")),
      menuItem("Differential Peaks", tabName = "diffpeaks", icon = icon("chart-simple")),
      menuItem("Gene Set Comparison", tabName = "genesetcompare", icon = icon("dna")),
      menuItem("Data Table",      tabName = "table",      icon = icon("table")),
      menuItem("About / Help",    tabName = "about",      icon = icon("circle-info"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #f4f6f9; }
      .box { border-radius: 6px; }
      .shiny-plot-output { background: white; }
      .selectize-control.multi .selectize-input {
        max-height: 210px;
        overflow-y: auto;
        display: flex;
        flex-wrap: wrap;
      }
    "))),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1 — LOAD DATA
      # -----------------------------------------------------------------------
      tabItem(tabName = "load",
        fluidRow(
          box(title = "Load ChromstaR Object (.RData)", width = 6, status = "primary",
            fileInput("chromstar_file", "Upload ChromstaR .RData file",
                      accept = c(".RData", ".rda", ".Rdata")),
            verbatimTextOutput("chromstar_summary")
          ),
          box(title = "Load Gene Annotation (TSV)", width = 6, status = "primary",
            fileInput("gtf_file", 
                      label = "Upload genes.tsv file",
                      accept = c(".tsv", ".txt", ".csv", ".gtf")),
            tags$small(tags$em(
              "Expected columns: chr, start, end, gene_id, gene_name, strand"
            )),
            br(),
            verbatimTextOutput("gtf_summary")
          )
        ),
        fluidRow(
          box(title = "Detected Marks & Conditions", width = 12, status = "info",
            fluidRow(
              column(4,
                h4("Histone Marks"),
                uiOutput("marks_ui")
              ),
              column(4,
                h4("Conditions"),
                uiOutput("conditions_ui")
              ),
              column(4,
                h4("Genome Info"),
                verbatimTextOutput("genome_info")
              )
            )
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2 — METAGENE PROFILE
      # -----------------------------------------------------------------------
      tabItem(tabName = "metagene",
        fluidRow(
          box(title = "Metagene Settings", width = 3, status = "warning",
            checkboxGroupInput("metagene_mode", "Reference (select one or more)",
                         choices = c("TSS" = "TSS", "TES" = "TES", "Gene body" = "gene_body"),
                         selected = "TSS"),
            conditionalPanel(
              condition = "input.metagene_mode.includes('TSS') || input.metagene_mode.includes('TES')",
              numericInput("upstream",   "Upstream (bp)",   value = 2000, min = 0, max = 50000, step = 100),
              numericInput("downstream", "Downstream (bp)", value = 2000, min = 0, max = 50000, step = 100),
              uiOutput("metagene_window_summary")
            ),
            numericInput("n_bins",     "Number of bins",  value = 100,  min = 20,  max = 500,   step = 10),
            hr(),
            h5("Gene scope"),
            radioButtons("gene_scope_meta", NULL,
                         choices = c("Only selected genes" = "selected",
                                     "All genes in GTF"     = "all"),
                         selected = "selected"),
            tags$small(style = "color:#888",
              "This plot is anchored to TSS/gene-body position, so it can only ever include bins that overlap a gene — there's no \"intergenic\" position to plot relative to. To see truly intergenic bins (not part of any gene), use the Region Browser tab with bin scope set to \"All bins, including intergenic.\""),
            hr(),
            h5("Select Genes"),
            fluidRow(
              column(6, actionButton("select_all_genes", "Select all",
                                     class = "btn-sm btn-block")),
              column(6, actionButton("clear_all_genes", "Clear all",
                                     class = "btn-sm btn-block"))
            ),
            br(),
            uiOutput("gene_selector_meta"),
            hr(),
            h5("Marks to display"),
            uiOutput("mark_selector_meta"),
            hr(),
            h5("Display mode"),
            radioButtons("compare_mode_meta", "Condition comparison",
                         choices = c("Side-by-side" = "facet",
                                     "Overlay"       = "overlay"),
                         selected = "facet"),
            checkboxInput("smooth_metagene",
                         "Smooth curve (LOESS) — recommended for many genes",
                         value = TRUE),
            actionButton("run_metagene", "Compute Profile",
                         class = "btn-primary btn-block")
          ),
          box(title = "Metagene Profile", width = 9, status = "primary",
            plotOutput("metagene_plot", height = "550px"),
            downloadButton("dl_metagene", "Download Plot"),
            downloadButton("dl_metagene_xlsx", "Download Data (Excel)")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2.5 — ENRICHMENT PROFILE (log observed/expected, like reference fig)
      # -----------------------------------------------------------------------
      tabItem(tabName = "enrichment",
        fluidRow(
          box(title = "Enrichment Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Shows log(observed/expected) enrichment of each mark around gene boundaries (TSS and TES), one panel per mark."),
            hr(),
            numericInput("enr_upstream",   "Upstream of TSS (bp)",   value = 2000, min = 100, max = 10000, step = 100),
            numericInput("enr_downstream", "Downstream of TES (bp)", value = 2000, min = 100, max = 10000, step = 100),
            numericInput("enr_n_bins",     "Bins per region",        value = 40,   min = 10,  max = 200,  step = 5),
            hr(),
            h5("Gene scope"),
            radioButtons("gene_scope_enr", NULL,
                         choices = c("Only selected genes" = "selected",
                                     "All genes in GTF"     = "all"),
                         selected = "selected"),
            tags$small(style = "color:#888",
              "This plot is anchored to gene boundaries (TSS/TES), so it can only ever include bins that overlap a gene's window — there's no \"intergenic\" position to plot relative to. To see truly intergenic bins, use the Region Browser tab with bin scope set to \"All bins, including intergenic.\""),
            hr(),
            h5("Select Genes"),
            fluidRow(
              column(6, actionButton("select_all_genes_enr", "Select all",
                                     class = "btn-sm btn-block")),
              column(6, actionButton("clear_all_genes_enr", "Clear all",
                                     class = "btn-sm btn-block"))
            ),
            br(),
            uiOutput("gene_selector_enr"),
            hr(),
            h5("Marks to display"),
            uiOutput("mark_selector_enr"),
            hr(),
            checkboxInput("smooth_enrichment", "Smooth curve (LOESS)", value = FALSE),
            actionButton("run_enrichment", "Compute Enrichment",
                         class = "btn-primary btn-block")
          ),
          box(title = "Histone Mark Enrichment", width = 9, status = "primary",
            plotOutput("enrichment_plot", height = "650px"),
            downloadButton("dl_enrichment_plot", "Download Plot"),
            downloadButton("dl_enrichment_xlsx", "Download Data (Excel)")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3 — REGION BROWSER
      # -----------------------------------------------------------------------
      tabItem(tabName = "browser",
        fluidRow(
          box(title = "Region Settings", width = 3, status = "warning",
            radioButtons("chrom_scope", "Chromosome scope",
                         choices = c("Single chromosome (zoom to position)" = "single",
                                     "Selected chromosomes (full length)"    = "multi",
                                     "All chromosomes (full length)"         = "all"),
                         selected = "single"),
            tags$small(style = "color:#888",
              "\"Single\" lets you zoom into a bp range. \"Selected\" or \"All\" show each chosen chromosome's entire length, faceted as separate rows — useful for a genome-wide overview but can be slow/cluttered with hundreds of contigs."),
            br(),
            uiOutput("chrom_selector"),
            conditionalPanel(
              condition = "input.chrom_scope == 'single'",
              numericInput("region_start", "Start (bp)", value = 1000000, min = 1),
              numericInput("region_end",   "End (bp)",   value = 1200000, min = 1)
            ),
            hr(),
            h5("Jump to gene"),
            uiOutput("gene_jump_ui"),
            numericInput("gene_flank", "Flanking (bp)", value = 5000, min = 0),
            actionButton("jump_to_gene", "Go", class = "btn-success"),
            hr(),
            h5("Bin scope"),
            radioButtons("bin_scope_browser", NULL,
                         choices = c("All bins, including intergenic" = "all",
                                     "Only bins inside a gene"         = "genic",
                                     "Only intergenic bins"            = "intergenic"),
                         selected = "all"),
            tags$small(style = "color:#888",
              "\"All bins\" shows every bin in the visible region regardless of gene annotation — this is the only mode that includes truly intergenic bins outside any gene. \"Only bins inside a gene\" or \"Only intergenic bins\" filter using the loaded gene annotation."),
            hr(),
            h5("Marks to display"),
            uiOutput("mark_selector_browser"),
            hr(),
            h5("Display mode"),
            radioButtons("compare_mode_browser", "Condition comparison",
                         choices = c("Side-by-side" = "facet",
                                     "Overlay"       = "overlay"),
                         selected = "facet"),
            actionButton("run_browser", "Load Region",
                         class = "btn-primary btn-block")
          ),
          box(title = "Genome Browser", width = 9, status = "primary",
            plotOutput("browser_plot", height = "600px"),
            downloadButton("dl_browser", "Download Plot")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4 — DATA TABLE
      # -----------------------------------------------------------------------
      tabItem(tabName = "table",
        fluidRow(
          box(title = "Chromatin State Segments", width = 12, status = "primary",
            downloadButton("dl_full_table", "Download full table (CSV)",
                           class = "btn-success"),
            tags$small(style = "color:#888; margin-left:10px;",
                       "Use this for all 574k+ rows — the table's own export buttons only handle small subsets reliably."),
            br(), br(),
            DTOutput("state_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5 — DIFFERENTIAL PEAKS
      # -----------------------------------------------------------------------
      tabItem(tabName = "diffpeaks",
        fluidRow(
          box(title = "Differential Peak Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Reproduces differential chromatin states analysis: filters merged chromatin segments by score, width, and a PA≠PNA state change, then counts how many filtered segments contain each mark in PA but not PNA (and vice versa)."),
            hr(),
            numericInput("diff_score_thresh", "Min differential score",
                        value = 0.9999, min = 0, max = 1, step = 0.0001),
            numericInput("diff_width_thresh", "Min merged region width (bp)",
                        value = 300, min = 0, step = 100),
            hr(),
            actionButton("run_diffpeaks", "Compute Differential Peaks",
                        class = "btn-primary btn-block")
          ),
          box(title = "Pairwise Differential Peaks per Histone Mark", width = 9, status = "primary",
            plotOutput("diffpeaks_plot", height = "550px"),
            downloadButton("dl_diffpeaks_plot", "Download Plot"),
            downloadButton("dl_diffpeaks_xlsx", "Download Data (Excel)")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5.5 — GENE SET COMPARISON (posterior probability comparison)
      # -----------------------------------------------------------------------
      tabItem(tabName = "genesetcompare",
        fluidRow(
          box(title = "Gene Set Comparison Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Compares a gene list's promoter posterior probability per mark against either a second gene list (e.g. up- vs down-regulated genes) or a permutation test against many random gene sets of the same size (a proper random background)."),
            hr(),
            numericInput("gsc_upstream",   "Upstream of TSS (bp)",   value = 2000, min = 0, max = 20000, step = 100),
            numericInput("gsc_downstream", "Downstream of TSS (bp)", value = 2000, min = 0, max = 20000, step = 100),
            hr(),
            radioButtons("gsc_replicate_mode", "Replicates",
                        choices = c("Average across replicates" = "average",
                                    "Show each replicate separately" = "separate"),
                        selected = "average"),
            tags$small(style = "color:#888",
              "\"Average\" gives one result per mark/condition (combines rep1, rep2, etc). \"Separate\" keeps each replicate as its own row/panel, useful for checking replicate consistency."),
            hr(),
            radioButtons("gsc_summary_stat", "Summary statistic",
                        choices = c("Mean" = "mean", "Median" = "median"),
                        selected = "mean"),
            tags$small(style = "color:#888",
              "Used for the group comparison (delta and the displayed stat_A/stat_B columns). Both mean and median are always shown for reference; this controls which one drives delta and the permutation test's central value."),
            hr(),
            radioButtons("gsc_mode", "Comparison type",
                        choices = c("Compare to a random background (permutation test)" = "permutation",
                                    "Compare to a second gene list I provide"             = "manual"),
                        selected = "permutation"),
            hr(),
            h5("Gene set A"),
            tags$textarea(
              id = "gsc_genes_a", rows = 8,
              style = "width:100%; font-family: monospace; font-size: 12px; resize: vertical;",
              placeholder = "Paste gene IDs here, one per line\n(e.g. your curated/stem-cell/up-regulated list)"
            ),
            uiOutput("gsc_match_summary_a"),
            hr(),
            conditionalPanel(
              condition = "input.gsc_mode == 'permutation'",
              numericInput("gsc_n_perm", "Number of random draws", value = 1000, min = 100, max = 10000, step = 100),
              tags$small(style = "color:#888",
                "Each draw picks a random set of genes from the full annotation, the same size as gene set A, and computes its mean promoter posterior per mark. The empirical p-value is the fraction of those random draws as extreme as your real gene set A.")
            ),
            conditionalPanel(
              condition = "input.gsc_mode == 'manual'",
              h5("Gene set B"),
              tags$textarea(
                id = "gsc_genes_b", rows = 8,
                style = "width:100%; font-family: monospace; font-size: 12px; resize: vertical;",
                placeholder = "Paste gene IDs here, one per line\n(e.g. your down-regulated list)"
              ),
              uiOutput("gsc_match_summary_b")
            ),
            hr(),
            actionButton("run_genesetcompare", "Compare Gene Sets",
                        class = "btn-primary btn-block")
          ),
          box(title = "Which histone marks differ most?", width = 9, status = "primary",
            tags$small(style = "color:#888",
              "Ranked by FDR (smallest first) — the marks at the top show the strongest, most significant difference."),
            uiOutput("gsc_summary_stat_label"),
            br(), br(),
            DTOutput("gsc_stats_table"),
            hr(),
            plotOutput("gsc_boxplot", height = "500px"),
            downloadButton("dl_gsc_plot", "Download Plot"),
            downloadButton("dl_gsc_xlsx", "Download Data (Excel)")
          )
        ),
        fluidRow(
          box(title = "Per-gene posterior probabilities (gene set A)", width = 12, status = "primary",
            tags$small(style = "color:#888",
              "One row per gene per mark/condition for gene set A — sort by 'posterior' to see which individual genes drive the signal for a given mark."),
            br(), br(),
            DTOutput("gsc_gene_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6 — ABOUT / HELP
      # -----------------------------------------------------------------------
      tabItem(tabName = "about",
        fluidRow(
          box(title = "About ChromstaR Viewer", width = 12, status = "info",
            tags$div(style = "max-width: 950px;",

              tags$h3("What this tool does"),
              tags$p("ChromstaR Viewer is a point-and-click app for exploring chromatin state and histone modification signal from a ChromstaR ", tags$code("combinedMultiHMM"), " object, alongside a gene annotation. It lets you compare experimental conditions (e.g. PA vs PNA) across multiple histone marks, browse specific genomic regions, and compare gene sets — all without writing any R code. Use the sidebar on the left to move between tabs; each tab is independent, so you can jump around in any order once data is loaded."),

              hr(),
              tags$h3("1. Load Data tab"),
              tags$p(tags$b("Purpose:"), " upload your data. Nothing else works until this is done."),
              tags$ul(
                tags$li(tags$b("\"Upload ChromstaR .RData file\""), " (left box) — click to browse for your ChromstaR ", tags$code("combinedMultiHMM"), " object (a ", tags$code(".RData"), "/", tags$code(".rda"), " file). A text box below it shows a summary once loaded (object class, bin count, etc)."),
                tags$li(tags$b("\"Upload genes.tsv file\""), " (right box) — click to browse for your gene annotation. This should be a simple table with columns ", tags$code("chr, start, end, gene_id, gene_name, strand"), " — convert a GTF to this format first if needed."),
                tags$li(tags$b("Detected Marks & Conditions"), " panel at the bottom fills in automatically once both files are loaded — it lists every histone mark, every condition, and basic genome info found in your object. Use this to sanity-check the upload before moving to other tabs.")
              ),

              hr(),
              tags$h3("2. Metagene Profile tab"),
              tags$p(tags$b("Purpose:"), " plot average signal (RPKM) around a reference point (TSS, TES, or across the gene body) for a chosen set of genes, one panel per histone mark."),
              tags$ul(
                tags$li(tags$b("\"Reference\""), " checkboxes — tick TSS, TES, and/or Gene body. You can tick more than one to compare them side by side in the same plot."),
                tags$li(tags$b("\"Upstream (bp)\" / \"Downstream (bp)\""), " — only shown when TSS or TES is ticked; sets how far before/after the reference point to plot. A live summary line below shows exactly what window this produces."),
                tags$li(tags$b("\"Number of bins\""), " — how many segments to divide the window/gene body into; more bins = finer resolution but a noisier-looking line."),
                tags$li(tags$b("\"Gene scope\""), " — \"Only selected genes\" uses whatever you've pasted below; \"All genes in GTF\" ignores the pasted list and runs on every gene in the annotation."),
                tags$li(tags$b("\"Select all\" / \"Clear all\""), " buttons — fill or empty the gene textbox below with every gene name from the annotation."),
                tags$li(tags$b("Gene textbox"), " — paste gene IDs here, one per line (e.g. copy-pasted from an Excel column). A match count appears underneath showing how many were recognised."),
                tags$li(tags$b("\"Marks to display\""), " checkboxes — tick which histone marks to include in the plot."),
                tags$li(tags$b("\"Condition comparison\""), " — \"Side-by-side\" puts each condition in its own facet panel; \"Overlay\" draws all conditions on the same axes with different colours."),
                tags$li(tags$b("\"Smooth curve (LOESS)\""), " checkbox — draws a smoothed line instead of the raw jagged signal; recommended when plotting many genes at once."),
                tags$li(tags$b("\"Compute Profile\""), " — runs the calculation and draws the plot. Nothing happens until you click this, even if you change settings above."),
                tags$li(tags$b("\"Download Plot\""), " — saves the current plot as a PDF."),
                tags$li(tags$b("\"Download Data (Excel)\""), " — saves the exact numbers behind the plot (one row per bin/mark/condition) as an .xlsx file.")
              ),

              hr(),
              tags$h3("3. Enrichment Profile tab"),
              tags$p(tags$b("Purpose:"), " plot log(observed/expected) enrichment of each mark around gene boundaries (TSS and TES together), one panel per mark — useful for seeing which marks are relatively enriched or depleted across a gene versus the genome-wide average."),
              tags$ul(
                tags$li(tags$b("\"Upstream of TSS (bp)\" / \"Downstream of TES (bp)\""), " — how far before the start and after the end of each gene to include."),
                tags$li(tags$b("\"Bins per region\""), " — resolution of the plot; more bins = finer detail."),
                tags$li(tags$b("\"Gene scope\", \"Select all\"/\"Clear all\", gene textbox, \"Marks to display\""), " — same behaviour as the equivalent controls in the Metagene Profile tab above."),
                tags$li(tags$b("\"Smooth curve (LOESS)\""), " checkbox — same as Metagene tab, smooths the line."),
                tags$li(tags$b("\"Compute Enrichment\""), " — runs the calculation and draws the plot."),
                tags$li(tags$b("\"Download Plot\" / \"Download Data (Excel)\""), " — same as Metagene tab: save the figure as PDF, or the underlying numbers as Excel.")
              ),

              hr(),
              tags$h3("4. Region Browser tab"),
              tags$p(tags$b("Purpose:"), " a genome-browser-style view of raw signal across a chosen chromosome and position range, or across whole chromosomes at once — for zooming into a specific locus."),
              tags$ul(
                tags$li(tags$b("\"Chromosome scope\""), " — \"Single chromosome\" lets you zoom into a specific bp range (most common use); \"Selected chromosomes\" or \"All chromosomes\" show entire chromosome(s) at full length instead, useful for a wide overview but can be slow/cluttered with many small contigs."),
                tags$li(tags$b("Chromosome dropdown"), " — pick which chromosome(s) to view (becomes a multi-select box in \"Selected chromosomes\" mode)."),
                tags$li(tags$b("\"Start (bp)\" / \"End (bp)\""), " — only shown in single-chromosome mode; defines the exact window to zoom into."),
                tags$li(tags$b("\"Jump to gene\""), " dropdown plus the ", tags$b("\"Flanking (bp)\""), " box and the ", tags$b("\"Go\""), " button — instead of typing coordinates by hand, pick a gene name here, set how much flanking sequence you want on each side, and click \"Go\" to jump the Start/End boxes straight to that gene's location."),
                tags$li(tags$b("\"Bin scope\""), " — \"All bins, including intergenic\" shows everything; \"Only bins inside a gene\" or \"Only intergenic bins\" filter the view using the loaded gene annotation."),
                tags$li(tags$b("\"Marks to display\""), " checkboxes — tick which histone marks to show as tracks."),
                tags$li(tags$b("\"Condition comparison\""), " — \"Side-by-side\" or \"Overlay\", same meaning as in the Metagene tab."),
                tags$li(tags$b("\"Load Region\""), " — fetches and draws the signal for the current settings."),
                tags$li(tags$b("\"Download Plot\""), " — saves the current view as a PDF.")
              ),

              hr(),
              tags$h3("5. Differential Peaks tab"),
              tags$p(tags$b("Purpose:"), " counts, per histone mark, how many chromatin segments are confidently present in one condition but not the other (and vice versa) — a bar chart summarising which marks change the most between conditions."),
              tags$ul(
                tags$li(tags$b("\"Min differential score\""), " — only keep chromatin segments with a confidence score at or above this threshold (closer to 1 = stricter, fewer but more confident segments)."),
                tags$li(tags$b("\"Min merged region width (bp)\""), " — discard segments shorter than this, to avoid counting tiny noisy regions."),
                tags$li(tags$b("\"Compute Differential Peaks\""), " — runs the filtering and counting, then draws the bar chart."),
                tags$li(tags$b("\"Download Plot\" / \"Download Data (Excel)\""), " — save the figure as PDF, or the underlying per-mark counts as Excel.")
              ),

              hr(),
              tags$h3("6. Gene Set Comparison tab"),
              tags$p(tags$b("Purpose:"), " compares the average histone-mark \"posterior probability\" (confidence that a mark is present) over the promoter region of one gene list against either a second gene list you provide, or a random background — useful for asking \"is this curated gene set unusual for a given mark?\" or \"do up- and down-regulated genes differ in chromatin state?\""),
              tags$ul(
                tags$li(tags$b("\"Upstream of TSS (bp)\" / \"Downstream of TSS (bp)\""), " — defines the promoter window averaged over for every gene."),
                tags$li(tags$b("\"Replicates\""), " — \"Average across replicates\" merges rep1/rep2/etc into one number per mark/condition; \"Show each replicate separately\" keeps them apart so you can check replicate consistency."),
                tags$li(tags$b("\"Summary statistic\""), " — choose Mean or Median as the main number used to compare the two groups (both are always shown side by side in the results table regardless of this choice)."),
                tags$li(tags$b("\"Comparison type\""), " — \"Compare to a random background\" runs a permutation test against many random gene sets of the same size (statistically rigorous, recommended default); \"Compare to a second gene list I provide\" lets you paste an actual second list instead (e.g. down-regulated genes)."),
                tags$li(tags$b("Gene set A textbox"), " — paste your main gene list here, one ID per line."),
                tags$li(tags$b("\"Number of random draws\""), " — only shown in random-background mode; how many random gene sets to sample for the statistical test. Higher = more precise p-values but slower (1000 is a good default; thousands for a final result)."),
                tags$li(tags$b("Gene set B textbox"), " — only shown in manual mode; paste your second gene list here."),
                tags$li(tags$b("\"Compare Gene Sets\""), " — runs the comparison and fills in the results table, plot, and per-gene table below."),
                tags$li(tags$b("Results table"), " (top right) — one row per mark/condition, sorted by FDR (most significant first); colour-coded green/yellow by significance."),
                tags$li(tags$b("Boxplot"), " — shows the spread of individual gene values for set A vs set B (or vs one representative random draw), one panel per mark."),
                tags$li(tags$b("\"Download Plot\" / \"Download Data (Excel)\""), " — save the boxplot as PDF, or the full results (stats + per-gene values) as a multi-sheet Excel file."),
                tags$li(tags$b("Per-gene table"), " (bottom) — one row per gene in set A, showing its individual posterior value per mark/condition; sortable, filterable by column, and exportable directly as CSV/Excel using the buttons above the table.")
              ),

              hr(),
              tags$h3("7. Data Table tab"),
              tags$p(tags$b("Purpose:"), " the full per-bin chromatin state table, annotated with which gene (if any) and which genomic zone (TSS, gene body thirds, upstream/downstream flanks, or intergenic) each bin falls into."),
              tags$ul(
                tags$li(tags$b("\"Download full table (CSV)\""), " — exports every row (500k+) to CSV. Use this rather than the table's own export, which only reliably handles small subsets."),
                tags$li(tags$b("Table"), " itself — scrollable and sortable by column; use the search boxes to filter by gene name, zone, or other fields.")
              ),

              hr(),
              tags$h3("Key concepts that apply across several tabs"),
              tags$ul(
                tags$li(tags$b("Gene scope"), " (Metagene/Enrichment tabs): \"only selected genes\" vs \"all genes in GTF.\" Both modes only ever include bins overlapping a gene — there's no \"intergenic\" position relative to a TSS or gene body."),
                tags$li(tags$b("Bin/Chromosome scope"), " (Region Browser tab): whether to see every bin (including intergenic), only genic bins, or only intergenic bins; and whether to zoom into one chromosome or view several/all at full length."),
                tags$li(tags$b("RPKM vs log(observed/expected)"), " — Metagene and Region Browser show raw mean RPKM signal. Enrichment Profile shows a log-ratio against the genome-wide average for that mark/condition, better for comparing marks with very different baseline signal levels."),
                tags$li(tags$b("Posterior probability"), " (Gene Set Comparison tab) — a 0-1 confidence score from ChromstaR that a given mark is genuinely present at a bin, distinct from the RPKM signal used elsewhere.")
              ),

              tags$h3("Performance tips"),
              tags$ul(
                tags$li("Computations are vectorised and should stay fast even with all genes selected, but very large gene sets combined with many marks can still take a few seconds."),
                tags$li("The Data Table's gene/zone annotation is computed once per session and cached — the first visit to that tab may take a little longer."),
                tags$li("Selecting \"All chromosomes\" in Region Browser, or running thousands of permutations in Gene Set Comparison, will be noticeably slower — start with smaller values to preview, then scale up for a final result.")
              ),

              hr(),
              tags$h3("Contact"),
              tags$p("Built for chromatin/epigenomics analysis in ", tags$em("Echinococcus multilocularis"), " by Janan Gawra."),
              tags$p(
                tags$a(href = "https://www.linkedin.com/in/janangawra/", target = "_blank",
                      icon("linkedin"), " linkedin.com/in/janangawra")
              )
            )
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # ---- Help button: jump to About/Help tab -----------------------------------
  observeEvent(input$show_help, {
    updateTabItems(session, "sidebar_tabs", "about")
  })

  # ---- Reactive values -------------------------------------------------------
  rv <- reactiveValues(
    hmm            = NULL,
    genes          = NULL,
    marks          = character(0),
    conditions     = character(0),
    chroms         = character(0),
    bin_annotation = NULL   # cached gene context per bin
  )

  # ---- Load ChromstaR object -------------------------------------------------
  observeEvent(input$chromstar_file, {
    req(input$chromstar_file)
    withProgress(message = "Loading ChromstaR object…", {
      tryCatch({
        hmm <- load_chromstar_object(input$chromstar_file$datapath)
        rv$hmm <- hmm

        parsed             <- parse_marks_conditions(hmm)
        rv$marks           <- parsed$marks
        rv$conditions      <- parsed$conditions
        rv$chroms          <- as.character(unique(seqnames(hmm$bins)))
        rv$bin_annotation  <- NULL  # reset cache on new object

        output$chromstar_summary <- renderPrint({
          bins_df   <- expand_bins_df(hmm)
          rpkm_cols <- grep("counts.rpkm", colnames(bins_df), value = TRUE)
          cat("Object class:", class(hmm), "\n")
          cat("Total bins:  ", length(hmm$bins), "\n")
          cat("Chromosomes: ", paste(head(rv$chroms, 5), collapse = ", "), "…\n")
          cat("Conditions:  ", paste(rv$conditions, collapse = ", "), "\n")
          cat("Marks:       ", paste(rv$marks, collapse = ", "), "\n\n")
          if (length(rpkm_cols) > 0) {
            cat("counts.rpkm columns (", length(rpkm_cols), "):\n")
            cat(" ", paste(head(rpkm_cols, 20), collapse = "\n  "), "\n")
          } else {
            cat("WARNING: No counts.rpkm.* columns found after expansion!\n")
            cat("All mcols names: ", paste(head(colnames(bins_df), 30), collapse = ", "), "\n")
          }
          if (!is.null(hmm$info)) {
            cat("\nhmm$info:\n")
            print(hmm$info[, c("mark","condition","replicate"), drop = FALSE])
          }
        })
      }, error = function(e) {
        showNotification(paste("Error loading ChromstaR object:", e$message), type = "error")
      })
    })
  })

  # ---- Load Gene TSV ---------------------------------------------------------
  observeEvent(input$gtf_file, {
    req(input$gtf_file)
    withProgress(message = "Loading gene table…", {
      tryCatch({
        path <- input$gtf_file$datapath

        # Detect separator (tab or comma)
        first_line <- readLines(path, n = 1)
        sep <- if (grepl("\t", first_line)) "\t" else ","

        df <- read.table(path, header = TRUE, sep = sep,
                         stringsAsFactors = FALSE, quote = "")

        # Normalise column names to lowercase
        colnames(df) <- tolower(colnames(df))

        # Required: chr/seqnames, start, end
        # Flexible naming
        if ("seqnames" %in% colnames(df)) df$chr <- df$seqnames
        if (!"chr"   %in% colnames(df)) stop("No 'chr' or 'seqnames' column found.")
        if (!"start" %in% colnames(df)) stop("No 'start' column found.")
        if (!"end"   %in% colnames(df)) stop("No 'end' column found.")

        # gene_name: prefer gene_name, else gene_id, else chr:start-end
        if (!"gene_name" %in% colnames(df)) {
          if ("gene_id" %in% colnames(df)) {
            df$gene_name <- df$gene_id
          } else {
            df$gene_name <- paste0(df$chr, ":", df$start, "-", df$end)
          }
        }
        if (!"gene_id" %in% colnames(df)) df$gene_id <- df$gene_name
        if (!"strand"  %in% colnames(df)) df$strand  <- "*"

        genes <- GRanges(
          seqnames  = df$chr,
          ranges    = IRanges(start = as.integer(df$start),
                              end   = as.integer(df$end)),
          strand    = df$strand,
          gene_id   = df$gene_id,
          gene_name = df$gene_name
        )
        rv$genes <- genes

      }, error = function(e) {
        showNotification(paste("Error loading gene file:", e$message), type = "error")
      })
    })
  })

  # ---- Dynamic UIs -----------------------------------------------------------
  output$marks_ui <- renderUI({
    req(length(rv$marks) > 0)
    tags$ul(lapply(rv$marks, tags$li))
  })

  output$conditions_ui <- renderUI({
    req(length(rv$conditions) > 0)
    tags$ul(lapply(rv$conditions, tags$li))
  })

  output$gtf_summary <- renderPrint({
    req(rv$genes)
    g       <- rv$genes
    nms     <- mcols(g)$gene_name
    has_ref <- sum(grepl("EmuJ_", nms, fixed = TRUE))
    cat("Genes loaded:    ", length(g), "\n")
    cat("With EmuJ name:  ", has_ref, "\n")
    cat("Without ref:     ", length(g) - has_ref, "\n")
    cat("Chromosomes:     ", paste(head(unique(as.character(seqnames(g))), 5),
                                    collapse = ", "), "…\n")
    cat("Example names:   ", paste(head(nms, 5), collapse = ", "), "\n")
  })

  output$genome_info <- renderPrint({
    req(rv$hmm)
    si <- seqinfo(rv$hmm$bins)
    print(si)
  })

  # Mark selectors (each needs its own renderUI call — reusing one renderUI
  # result across two outputs does not work correctly in Shiny)
  output$mark_selector_meta <- renderUI({
    req(length(rv$marks) > 0)
    checkboxGroupInput("marks_meta", NULL,
                       choices  = rv$marks,
                       selected = rv$marks)
  })

  output$mark_selector_browser <- renderUI({
    req(length(rv$marks) > 0)
    checkboxGroupInput("marks_browser", NULL,
                       choices  = rv$marks,
                       selected = rv$marks)
  })

  # Gene selectors
  gene_names_reactive <- reactive({
    req(rv$genes)
    cols <- colnames(mcols(rv$genes))
    # Try common GTF gene name columns in order of preference
    for (col in c("gene_name", "Name", "gene_id", "ID", "name", "transcript_name")) {
      if (col %in% cols) {
        nms <- as.character(mcols(rv$genes)[[col]])
        nms[is.na(nms)] <- paste0("gene_", which(is.na(nms)))
        return(nms)
      }
    }
    # Last resort: row index with chr:start-end label
    paste0(seqnames(rv$genes), ":", start(rv$genes), "-", end(rv$genes))
  })

  # Also store which GTF column was used, shown in Load tab
  output$gtf_summary <- renderPrint({
    req(rv$genes)
    cols     <- colnames(mcols(rv$genes))
    name_col <- NA
    for (col in c("gene_name", "Name", "gene_id", "ID", "name", "transcript_name")) {
      if (col %in% cols) { name_col <- col; break }
    }
    cat("Genes loaded:   ", length(rv$genes), "\n")
    cat("Name column:    ", ifelse(is.na(name_col), "none found — using index", name_col), "\n")
    cat("All columns:    ", paste(cols, collapse = ", "), "\n")
    cat("Chromosomes:    ", paste(head(unique(as.character(seqnames(rv$genes))), 10), collapse = ", "), "\n")
  })

  output$gene_selector_meta <- renderUI({
    req(rv$genes)
    nms      <- gene_names_reactive()
    # Default: first 20 genes selected — fast to compute, user can change
    defaults <- paste(head(nms, 20), collapse = "\n")
    tagList(
      tags$small(style = "color:#888",
        paste0(length(nms), " genes in GTF — paste a list from Excel (one gene per row) or use the buttons above. Only bins overlapping the selected genes' TSS/body window are used; intergenic bins are excluded automatically.")),
      br(), br(),
      tags$textarea(
        id = "selected_genes_text",
        rows = 14,
        style = "width:100%; height:280px; font-family: monospace; font-size: 13px; resize: vertical;",
        placeholder = "Paste gene IDs here, one per line (e.g. from an Excel column)…",
        defaults
      ),
      uiOutput("gene_match_summary")
    )
  })

  # Parse the textarea into a clean vector of gene names, matched against the GTF
  selected_genes_parsed <- reactive({
    req(input$selected_genes_text)
    raw <- strsplit(input$selected_genes_text, "[\r\n,;\t]+")[[1]]
    raw <- trimws(raw)
    raw[nchar(raw) > 0]
  })

  output$gene_match_summary <- renderUI({
    req(rv$genes)
    pasted <- selected_genes_parsed()
    if (length(pasted) == 0) {
      return(tags$small(style = "color:#888", "No genes entered yet."))
    }
    nms      <- gene_names_reactive()
    matched  <- pasted[pasted %in% nms]
    unmatched <- setdiff(pasted, nms)

    tagList(
      br(),
      tags$small(style = "color:#28a745",
        paste0("✓ ", length(matched), " of ", length(pasted), " gene names matched.")),
      if (length(unmatched) > 0) {
        tags$div(
          tags$small(style = "color:#dc3545",
            paste0("✗ ", length(unmatched), " not found, e.g.: ",
                   paste(head(unmatched, 5), collapse = ", "),
                   if (length(unmatched) > 5) "…" else ""))
        )
      }
    )
  })

  observeEvent(input$select_all_genes, {
    req(rv$genes)
    nms <- gene_names_reactive()
    updateTextAreaInput(session, "selected_genes_text",
                        value = paste(nms, collapse = "\n"))
  })

  observeEvent(input$clear_all_genes, {
    updateTextAreaInput(session, "selected_genes_text", value = "")
  })

  # ---- Enrichment tab: mark + gene selectors (mirrors metagene tab) ---------
  output$mark_selector_enr <- renderUI({
    req(length(rv$marks) > 0)
    checkboxGroupInput("marks_enr", NULL,
                       choices  = rv$marks,
                       selected = rv$marks)
  })

  output$gene_selector_enr <- renderUI({
    req(rv$genes)
    nms      <- gene_names_reactive()
    defaults <- paste(head(nms, 20), collapse = "\n")
    tagList(
      tags$small(style = "color:#888",
        paste0(length(nms), " genes in GTF — paste a list from Excel (one gene per row) or use the buttons above. Only bins overlapping the selected genes (plus flanks) are used; intergenic bins are excluded automatically.")),
      br(), br(),
      tags$textarea(
        id = "selected_genes_text_enr",
        rows = 14,
        style = "width:100%; height:280px; font-family: monospace; font-size: 13px; resize: vertical;",
        placeholder = "Paste gene IDs here, one per line (e.g. from an Excel column)…",
        defaults
      ),
      uiOutput("gene_match_summary_enr")
    )
  })

  selected_genes_parsed_enr <- reactive({
    req(input$selected_genes_text_enr)
    raw <- strsplit(input$selected_genes_text_enr, "[\r\n,;\t]+")[[1]]
    raw <- trimws(raw)
    raw[nchar(raw) > 0]
  })

  output$gene_match_summary_enr <- renderUI({
    req(rv$genes)
    pasted <- selected_genes_parsed_enr()
    if (length(pasted) == 0) {
      return(tags$small(style = "color:#888", "No genes entered yet."))
    }
    nms       <- gene_names_reactive()
    matched   <- pasted[pasted %in% nms]
    unmatched <- setdiff(pasted, nms)

    tagList(
      br(),
      tags$small(style = "color:#28a745",
        paste0("✓ ", length(matched), " of ", length(pasted), " gene names matched.")),
      if (length(unmatched) > 0) {
        tags$div(
          tags$small(style = "color:#dc3545",
            paste0("✗ ", length(unmatched), " not found, e.g.: ",
                   paste(head(unmatched, 5), collapse = ", "),
                   if (length(unmatched) > 5) "…" else ""))
        )
      }
    )
  })

  observeEvent(input$select_all_genes_enr, {
    req(rv$genes)
    nms <- gene_names_reactive()
    updateTextAreaInput(session, "selected_genes_text_enr",
                        value = paste(nms, collapse = "\n"))
  })

  observeEvent(input$clear_all_genes_enr, {
    updateTextAreaInput(session, "selected_genes_text_enr", value = "")
  })

  # ---- Enrichment computation & plot -----------------------------------------
  enrichment_data <- eventReactive(input$run_enrichment, {
    req(rv$hmm, rv$genes, input$marks_enr)

    if (identical(input$gene_scope_enr, "all")) {
      genes_sel <- rv$genes
    } else {
      nms <- gene_names_reactive()
      sel <- selected_genes_parsed_enr()
      sel <- sel[sel %in% nms]

      if (length(sel) == 0) {
        showNotification(
          "Please enter at least one valid gene name, or switch scope to \"All genes in GTF\".",
          type = "warning", duration = 5
        )
        return(NULL)
      }
      genes_sel <- rv$genes[nms %in% sel]
    }

    if (length(genes_sel) > 500) {
      showNotification(
        paste0(length(genes_sel), " genes — running vectorised computation, should still be fast."),
        type = "message", duration = 5
      )
    }

    withProgress(message = paste0("Computing enrichment (", length(genes_sel), " genes)…"), {
      compute_enrichment_profile(
        hmm        = rv$hmm,
        genes_gr   = genes_sel,
        marks      = input$marks_enr,
        conditions = rv$conditions,
        upstream   = input$enr_upstream,
        downstream = input$enr_downstream,
        n_bins     = input$enr_n_bins
      )
    })
  })

  build_enrichment_plot <- reactive({
    df <- enrichment_data()
    req(df)

    n_bins <- input$enr_n_bins
    df$x <- df$slot / n_bins

    # Clip to the intended range; smoothing/binning artifacts can occasionally
    # produce points slightly outside [-1, 2]
    df <- df[df$x >= -1 & df$x <= 2, ]

    x_breaks <- c(-1, -0.5, 0, 0.25, 0.5, 0.75, 1, 1.5, 2)
    x_labels <- c(
      paste0("-", input$enr_upstream, "bp"), paste0("-", round(input$enr_upstream/2), "bp"),
      "TSS\n(0%)", "25%", "50%", "75%", "TES\n(100%)",
      paste0("+", round(input$enr_downstream/2), "bp"), paste0("+", input$enr_downstream, "bp")
    )

    conds <- rv$conditions
    line_types <- setNames(rep(c("solid", "dotted", "dashed", "dotdash"),
                               length.out = length(conds)), conds)
    line_cols  <- setNames(c("#E63946", "#2196F3", "#FF9800", "#4CAF50")[seq_along(conds)], conds)

    # Shade the gene-body zone (x in [0,1]) so it's visually distinct from
    # the upstream/downstream flanks even before looking at axis labels
    body_shade <- data.frame(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf)

    p <- ggplot(df, aes(x = x, y = mean_log_ratio,
                        colour = condition, linetype = condition)) +
      geom_rect(data = body_shade,
               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
               inherit.aes = FALSE, fill = "grey85", alpha = 0.35)

    if (isTRUE(input$smooth_enrichment)) {
      p <- p + geom_smooth(method = "loess", span = 0.2, se = FALSE, linewidth = 0.9)
    } else {
      p <- p + geom_line(linewidth = 0.8)
    }

    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_vline(xintercept = 0, linetype = "solid",  colour = "grey30", linewidth = 0.6) +
      geom_vline(xintercept = 1, linetype = "solid",  colour = "grey30", linewidth = 0.6) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels,
                         expand = expansion(mult = 0.02)) +
      scale_colour_manual(values = line_cols) +
      scale_linetype_manual(values = line_types) +
      facet_wrap(~ mark, scales = "free_y", ncol = 4) +
      labs(
        x = NULL,
        y = "log(observed/expected)",
        colour = "Condition", linetype = "Condition",
        title    = "Histone Mark Enrichment — Per Mark",
        subtitle = paste0(paste(conds, collapse = " vs "),
                          "  |  grey band = gene body, vertical lines = TSS and TES")
      ) +
      theme_bw(base_size = 13) +
      theme(
        strip.background = element_rect(fill = "#34495e"),
        strip.text       = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "top",
        axis.text.x      = element_text(size = 8)
      )
    p
  })

  output$enrichment_plot <- renderPlot({ build_enrichment_plot() })

  output$dl_enrichment_plot <- downloadHandler(
    filename = function() paste0("enrichment_profile_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggsave(file, plot = build_enrichment_plot(),
             width = 14, height = 8, device = "pdf")
    }
  )

  output$dl_enrichment_xlsx <- downloadHandler(
    filename = function() paste0("enrichment_data_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      df <- enrichment_data()
      req(df)
      df_out <- df
      df_out$x_position <- df_out$slot / input$enr_n_bins
      df_out <- df_out[, c("mark", "condition", "slot", "x_position",
                           "mean_log_ratio", "n_genes")]
      write.xlsx(df_out, file, overwrite = TRUE)
    }
  )

  output$gene_jump_ui <- renderUI({
    req(rv$genes, input$region_chrom)
    chrom_for_genes <- input$region_chrom[1]  # jump-to-gene only makes sense for one chrom at a time
    nms     <- gene_names_reactive()
    on_chr  <- as.character(seqnames(rv$genes)) == chrom_for_genes
    choices <- nms[on_chr]
    tagList(
      tags$small(style = "color:#888",
        paste0(length(choices), " genes on ", chrom_for_genes)),
      selectizeInput("jump_gene", NULL,
                     choices  = choices,
                     multiple = FALSE,
                     options  = list(maxOptions = 11200))
    )
  })

  output$chrom_selector <- renderUI({
    req(length(rv$chroms) > 0)
    if (identical(input$chrom_scope, "all")) {
      tags$small(style = "color:#888",
        paste0("All ", length(rv$chroms), " chromosomes/contigs will be shown — this can be slow with hundreds of sequences."))
    } else if (identical(input$chrom_scope, "multi")) {
      selectizeInput("region_chrom", "Chromosomes (select one or more)",
                     choices  = rv$chroms,
                     selected = rv$chroms[1],
                     multiple = TRUE,
                     options  = list(maxOptions = length(rv$chroms) + 10))
    } else {
      selectInput("region_chrom", "Chromosome",
                  choices  = rv$chroms,
                  selected = rv$chroms[1])
    }
  })

  # ---- Jump to gene ----------------------------------------------------------
  observeEvent(input$jump_to_gene, {
    req(rv$genes, input$jump_gene, input$gene_flank)
    nms   <- gene_names_reactive()
    idx   <- which(nms == input$jump_gene)[1]
    if (is.na(idx)) return()
    g     <- rv$genes[idx]
    flank <- input$gene_flank
    updateNumericInput(session, "region_start", value = max(1, start(g) - flank))
    updateNumericInput(session, "region_end",   value = end(g) + flank)
  })

  # ---- METAGENE PLOT ---------------------------------------------------------
  output$metagene_window_summary <- renderUI({
    req(input$metagene_mode, input$upstream, input$downstream)
    ref_labels <- intersect(c("TSS", "TES"), input$metagene_mode)
    if (length(ref_labels) == 0) return(NULL)

    lines <- lapply(ref_labels, function(ref_label) {
      tags$div(
        tags$b(paste0(ref_label, " window: ")),
        paste0(
          ref_label, " −", format(input$upstream, big.mark = ","), "bp",
          "  to  ",
          ref_label, " +", format(input$downstream, big.mark = ","), "bp",
          "  (", format(input$upstream + input$downstream, big.mark = ","), "bp total)"
        )
      )
    })
    tags$small(style = "color:#555; display:block; margin-top:6px;", lines)
  })

  metagene_data <- eventReactive(input$run_metagene, {
    req(rv$hmm, rv$genes, input$marks_meta, input$metagene_mode)

    if (identical(input$gene_scope_meta, "all")) {
      genes_sel <- rv$genes
    } else {
      nms <- gene_names_reactive()
      sel <- selected_genes_parsed()
      sel <- sel[sel %in% nms]   # keep only genes that actually matched

      if (length(sel) == 0) {
        showNotification(
          "Please enter at least one valid gene name, or switch scope to \"All genes in GTF\".",
          type = "warning", duration = 5
        )
        return(NULL)
      }
      genes_sel <- rv$genes[nms %in% sel]
    }

    # Warn if too many genes (still runs, just slow)
    if (length(genes_sel) > 500) {
      showNotification(
        paste0(length(genes_sel), " genes — running vectorised computation, should still be fast."),
        type = "message", duration = 5
      )
    }

    withProgress(message = paste0("Computing metagene (", length(genes_sel), " genes)…"), {
      modes_sel <- input$metagene_mode
      n_modes   <- length(modes_sel)
      results   <- list()

      for (i in seq_along(modes_sel)) {
        m <- modes_sel[i]
        incProgress(1 / n_modes, detail = paste0("Reference: ", m))
        res <- compute_metagene(
          hmm        = rv$hmm,
          genes_gr   = genes_sel,
          marks      = input$marks_meta,
          conditions = rv$conditions,
          mode       = m,
          upstream   = input$upstream,
          downstream = input$downstream,
          n_bins     = input$n_bins
        )
        if (!is.null(res) && nrow(res) > 0) {
          res$ref_type <- m
          results[[m]] <- res
        }
      }

      if (length(results) == 0) return(NULL)
      bind_rows(results)
    })
  })

  build_metagene_plot <- reactive({
    df <- metagene_data()
    req(df)

    n_bins  <- input$n_bins
    modes_sel <- unique(df$ref_type)
    multi_ref <- length(modes_sel) > 1

    x_breaks <- c(1, round(n_bins / 2), n_bins)

    label_for_mode <- function(m) {
      if (m == "TSS") {
        c(paste0("-", input$upstream / 1000, "kb"), "TSS", paste0("+", input$downstream / 1000, "kb"))
      } else if (m == "TES") {
        c(paste0("-", input$upstream / 1000, "kb"), "TES", paste0("+", input$downstream / 1000, "kb"))
      } else {
        c("Start", "50%", "End")
      }
    }

    p <- ggplot(df, aes(x = bin_idx, y = mean_signal,
                        colour = condition, fill = condition)) +
      geom_ribbon(aes(ymin = mean_signal - se_signal,
                      ymax = mean_signal + se_signal),
                  alpha = 0.2, colour = NA)

    if (isTRUE(input$smooth_metagene)) {
      p <- p + geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 0.9)
    } else {
      p <- p + geom_line(linewidth = 0.7)
    }

    if (multi_ref) {
      # When multiple reference modes are shown together, x-axis meaning differs
      # per column (TSS/TES = bp distance, gene_body = %), so we can't use one
      # global scale_x_continuous with fixed labels — instead facet by ref_type
      # and rely on the facet strip text ("TSS"/"TES"/"gene_body") for context,
      # using a single consistent breaks/labels set that's "close enough"
      # (Start/Mid/End) across all reference types for visual alignment.
      p <- p + scale_x_continuous(breaks = x_breaks, labels = c("Start", "Mid", "End"))
      facet_formula <- if (input$compare_mode_meta == "facet") {
        mark ~ ref_type + condition
      } else {
        mark ~ ref_type
      }
      p <- p + facet_grid(facet_formula, scales = "free_y")
      title_txt <- paste("Metagene profile —", paste(modes_sel, collapse = " / "))
    } else {
      x_labels <- label_for_mode(modes_sel[1])
      p <- p + scale_x_continuous(breaks = x_breaks, labels = x_labels)
      if (input$compare_mode_meta == "facet") {
        p <- p + facet_grid(mark ~ condition, scales = "free_y")
      } else {
        p <- p + facet_wrap(~ mark, scales = "free_y", ncol = 2)
      }
      title_txt <- paste("Metagene profile —", modes_sel[1])
    }

    p <- p +
      scale_colour_manual(values = c("#E63946", "#2196F3",
                                     "#FF9800", "#4CAF50")[seq_along(rv$conditions)]) +
      scale_fill_manual(values   = c("#E63946", "#2196F3",
                                     "#FF9800", "#4CAF50")[seq_along(rv$conditions)]) +
      labs(x = NULL, y = "Mean RPKM (averaged across replicates)",
           colour = "Condition", fill = "Condition",
           title = title_txt,
           caption = paste0(length(unique(df$bin_idx)), " bins | ",
                            "Ribbon = ±1 SE across selected genes")) +
      theme_bw(base_size = 13) +
      theme(
        strip.background = element_rect(fill = "#34495e"),
        strip.text       = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "top"
      )

    p
  })

  output$metagene_plot <- renderPlot({ build_metagene_plot() })

  output$dl_metagene <- downloadHandler(
    filename = function() paste0("metagene_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggsave(file, plot = build_metagene_plot(),
             width = 14, height = 10, device = "pdf")
    }
  )

  output$dl_metagene_xlsx <- downloadHandler(
    filename = function() paste0("metagene_data_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      df <- metagene_data()
      req(df)
      write.xlsx(df, file, overwrite = TRUE)
    }
  )

  # ---- REGION BROWSER PLOT ---------------------------------------------------
  browser_data <- eventReactive(input$run_browser, {
    req(rv$hmm)

    marks_use <- input$marks_browser
    if (is.null(marks_use) || length(marks_use) == 0) {
      marks_use <- rv$marks  # fallback: use all marks if checkboxes haven't registered
      if (length(marks_use) == 0) {
        showNotification("No histone marks available — load a ChromstaR object first.",
                         type = "error", duration = 6)
        return(NULL)
      }
    }

    seqinfo_obj <- seqinfo(rv$hmm$bins)

    if (identical(input$chrom_scope, "all")) {
      chroms_use <- rv$chroms
    } else if (identical(input$chrom_scope, "multi")) {
      req(input$region_chrom)
      chroms_use <- input$region_chrom
    } else {
      req(input$region_chrom, input$region_start, input$region_end)
      chroms_use <- input$region_chrom
    }

    if (identical(input$chrom_scope, "single")) {
      region_gr <- GRanges(seqnames = chroms_use,
                           ranges   = IRanges(start = input$region_start,
                                              end   = input$region_end))
    } else {
      # Full length of each selected chromosome
      lens <- seqlengths(seqinfo_obj)[chroms_use]
      lens[is.na(lens)] <- 1e9  # fallback if seqlengths missing
      region_gr <- GRanges(seqnames = chroms_use,
                           ranges   = IRanges(start = 1, end = lens))
    }

    scope <- input$bin_scope_browser
    if (scope %in% c("genic", "intergenic") && (is.null(rv$genes) || length(rv$genes) == 0)) {
      showNotification(
        "No gene annotation loaded — showing all bins instead.",
        type = "warning", duration = 5
      )
      scope <- "all"
    }

    if (length(chroms_use) > 20) {
      showNotification(
        paste0(length(chroms_use), " chromosomes selected — this may be slow and the plot may be cluttered."),
        type = "warning", duration = 6
      )
    }

    withProgress(message = paste0("Extracting signal (", length(chroms_use), " chromosome(s))…"), {
      extract_signal_region(
        hmm        = rv$hmm,
        region_gr  = region_gr,
        marks      = marks_use,
        conditions = rv$conditions,
        genes_gr   = rv$genes,
        bin_scope  = scope
      )
    })
  })

  build_browser_plot <- reactive({
    df <- browser_data()
    req(df)

    n_chroms_shown <- length(unique(df$chr))
    multi_chrom    <- n_chroms_shown > 1

    # Gene annotation track (optional overlay) — only meaningful in single-chrom mode,
    # since position scales differ across chromosomes
    gene_annot <- NULL
    if (!multi_chrom && !is.null(rv$genes) && identical(input$chrom_scope, "single")) {
      region_gr <- GRanges(seqnames = input$region_chrom,
                           ranges   = IRanges(start = input$region_start,
                                              end   = input$region_end))
      ov  <- findOverlaps(rv$genes, region_gr)
      if (length(ov) > 0) {
        g_sub <- rv$genes[queryHits(ov)]
        nms   <- gene_names_reactive()
        gene_annot <- data.frame(
          start = start(g_sub),
          end   = end(g_sub),
          y     = 0,
          label = nms[queryHits(ov)],
          strand = as.character(strand(g_sub)),
          stringsAsFactors = FALSE
        )
      }
    }

    p <- ggplot(df, aes(x = position, y = signal,
                        colour = condition, fill = condition)) +
      geom_area(alpha = 0.35, position = "identity") +
      geom_line(linewidth = 0.7) +
      scale_x_continuous(labels = label_comma()) +
      scale_colour_manual(values = c("#E63946", "#2196F3",
                                     "#FF9800", "#4CAF50")[seq_along(rv$conditions)]) +
      scale_fill_manual(values   = c("#E63946", "#2196F3",
                                     "#FF9800", "#4CAF50")[seq_along(rv$conditions)]) +
      theme_bw(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "#2c3e50"),
        strip.text       = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "top"
      )

    if (multi_chrom) {
      p <- p +
        labs(x = "Position (bp, per chromosome)",
             y = "RPKM (averaged across replicates)",
             colour = "Condition", fill = "Condition",
             title  = paste0(n_chroms_shown, " chromosomes — ", input$bin_scope_browser, " bins"))
      # Facet by chromosome AND mark; condition stays as colour
      if (input$compare_mode_browser == "facet") {
        p <- p + facet_grid(mark ~ chr, scales = "free")
      } else {
        p <- p + facet_wrap(chr ~ mark, scales = "free", ncol = length(unique(df$mark)))
      }
    } else {
      p <- p +
        labs(x = paste("Position on", unique(df$chr)[1]),
             y = "RPKM (averaged across replicates)",
             colour = "Condition", fill = "Condition",
             title  = if (identical(input$chrom_scope, "single")) {
               paste0(input$region_chrom, ":",
                     format(input$region_start, big.mark = ","), "–",
                     format(input$region_end,   big.mark = ","))
             } else {
               paste0(unique(df$chr)[1], " (full length)")
             })
      if (input$compare_mode_browser == "facet") {
        p <- p + facet_grid(mark ~ condition, scales = "free_y")
      } else {
        p <- p + facet_wrap(~ mark, scales = "free_y", ncol = 2)
      }

      # Add gene annotation track only for single-chromosome zoomed view
      if (!is.null(gene_annot)) {
        p <- p +
          geom_segment(data = gene_annot,
                       aes(x = start, xend = end, y = -Inf, yend = -Inf),
                       colour = "#27ae60", linewidth = 3,
                       inherit.aes = FALSE) +
          geom_text(data = gene_annot,
                    aes(x = (start + end) / 2, y = -Inf, label = label),
                    vjust = -0.3, size = 3, colour = "#27ae60",
                    inherit.aes = FALSE)
      }
    }

    p
  })

  output$browser_plot <- renderPlot({ build_browser_plot() })

  output$dl_browser <- downloadHandler(
    filename = function() paste0("region_browser_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggsave(file, plot = build_browser_plot(),
             width = 16, height = 12, device = "pdf")
    }
  )

  # ---- DATA TABLE ------------------------------------------------------------
  # Gene annotation is computed once and cached in rv$bin_annotation
  observe({
    req(rv$hmm, rv$genes)
    if (!is.null(rv$bin_annotation)) return()  # already computed

    withProgress(message = "Annotating bins (one-time, ~5 sec)…", {

      genes   <- rv$genes
      g_names <- mcols(genes)$gene_name
      s_g     <- start(genes)
      e_g     <- end(genes)
      str_g   <- as.character(strand(genes))
      glen    <- e_g - s_g

      # TSS / TES vectorised
      tss <- ifelse(str_g == "-", e_g, s_g)
      tes <- ifelse(str_g == "-", s_g, e_g)

      chr_g <- as.character(seqnames(genes))

      # Helper: build a GRanges for all genes at once, clip negatives
      make_zone_gr <- function(starts, ends, names_vec, chr_vec) {
        ok <- starts >= 1 & ends >= 1 & starts <= ends
        GRanges(
          seqnames  = chr_vec[ok],
          ranges    = IRanges(pmax(1L, starts[ok]), pmax(1L, ends[ok])),
          gene_name = names_vec[ok],
          zone_idx  = which(ok)    # index back to original gene row
        )
      }

      # Build all 8 zone GRanges in one shot (vectorised over all genes)
      zone_list <- list(
        TSS_pm200          = make_zone_gr(tss - 200,        tss + 200,        g_names, chr_g),
        gene_body_5prime   = make_zone_gr(s_g,              s_g + floor(glen/3),      g_names, chr_g),
        gene_body_middle   = make_zone_gr(s_g + floor(glen/3)+1, s_g + floor(2*glen/3), g_names, chr_g),
        gene_body_3prime   = make_zone_gr(s_g + floor(2*glen/3)+1, e_g,               g_names, chr_g),
        upstream_1000_0    = make_zone_gr(tss - 1000,       tss - 1,          g_names, chr_g),
        upstream_2000_1000 = make_zone_gr(tss - 2000,       tss - 1001,       g_names, chr_g),
        TES_0_1000         = make_zone_gr(tes + 1,          tes + 1000,       g_names, chr_g),
        TES_1000_2000      = make_zone_gr(tes + 1001,       tes + 2000,       g_names, chr_g)
      )

      # Priority (lower = higher priority)
      zone_priority <- c(
        TSS_pm200 = 1L, gene_body_5prime = 2L, gene_body_middle = 3L,
        gene_body_3prime = 4L, upstream_1000_0 = 5L, upstream_2000_1000 = 6L,
        TES_0_1000 = 7L, TES_1000_2000 = 8L
      )

      bins    <- rv$hmm$bins
      n_bins  <- length(bins)
      best_pri   <- rep(9L,              n_bins)
      best_gene  <- rep(NA_character_,   n_bins)
      best_zone  <- rep("Intergenic",    n_bins)

      # 8 findOverlaps calls total — fast
      for (zname in names(zone_list)) {
        zr  <- zone_list[[zname]]
        if (length(zr) == 0) next
        ov  <- findOverlaps(bins, zr, ignore.strand = TRUE)
        bi  <- queryHits(ov)
        pri <- zone_priority[zname]
        upd <- bi[best_pri[bi] > pri]
        if (length(upd) == 0) next
        zi  <- subjectHits(ov)[best_pri[bi] > pri]
        best_pri[upd]  <- pri
        best_gene[upd] <- mcols(zr)$gene_name[zi]
        best_zone[upd] <- zname
      }

      zone_labels <- c(
        TSS_pm200          = "TSS ±200bp",
        gene_body_5prime   = "Gene body 5′ (1/3)",
        gene_body_middle   = "Gene body middle (2/3)",
        gene_body_3prime   = "Gene body 3′ (3/3)",
        upstream_1000_0    = "Upstream 1000–0bp",
        upstream_2000_1000 = "Upstream 2000–1000bp",
        TES_0_1000         = "TES 0–1000bp",
        TES_1000_2000      = "TES 1000–2000bp",
        Intergenic         = "Intergenic"
      )

      rv$bin_annotation <- data.frame(
        gene_name   = best_gene,
        genomic_zone = zone_labels[best_zone],
        stringsAsFactors = FALSE
      )
    })
  })

  full_table_df <- reactive({
    req(rv$hmm, rv$bin_annotation)

    bins_df <- expand_bins_df(rv$hmm)
    state_cols <- intersect(
      c("chr", "start", "end", "transition.group", "state",
        "combination.PA", "combination.PNA",
        "differential.score", "maxPostInPeak"),
      colnames(bins_df)
    )

    df <- cbind(
      bins_df[, state_cols, drop = FALSE],
      rv$bin_annotation
    )
    df$width  <- df$end - df$start + 1
    df$bin_id <- paste0(df$chr, ":", df$start, "-", df$end)

    # Add mean RPKM per mark/condition (averaged across replicates) so the
    # exported table includes the actual quantitative signal, not just
    # categorical state/combination columns
    rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
    for (cond in rv$conditions) {
      for (mark in rv$marks) {
        pattern <- paste0("counts.rpkm.", mark, ".", cond)
        cols    <- grep(pattern, rpkm_cols_all, value = TRUE, fixed = TRUE)
        if (length(cols) == 0) next
        colname <- paste0("rpkm_", mark, "_", cond)
        df[[colname]] <- round(rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE), 4)
      }
    }

    # Move gene/identifier cols to front
    front <- c("bin_id", "chr", "start", "end", "width", "gene_name", "genomic_zone")
    front <- intersect(front, colnames(df))
    df[, c(front, setdiff(colnames(df), front)), drop = FALSE]
  })

  output$state_table <- renderDT({
    df <- full_table_df()
    req(df)

    datatable(
      df,
      filter     = "top",
      rownames   = FALSE,
      options    = list(
        pageLength = 20,
        scrollX    = TRUE
      )
    ) %>%
      formatStyle("genomic_zone",
        backgroundColor = styleEqual(
          c("TSS ±200bp", "Gene body 5′ (1/3)", "Gene body middle (2/3)",
            "Gene body 3′ (3/3)", "Upstream 1000–0bp", "Upstream 2000–1000bp",
            "TES 0–1000bp", "TES 1000–2000bp", "Intergenic"),
          c("#fff3cd", "#d4edda", "#c3e6cb", "#b1dfbb",
            "#d1ecf1", "#bee5eb", "#f8d7da", "#f5c6cb", "#e2e3e5")
        )
      )
  })

  output$dl_full_table <- downloadHandler(
    filename = function() paste0("chromatin_state_table_", Sys.Date(), ".csv"),
    content  = function(file) {
      write.csv(full_table_df(), file, row.names = FALSE)
    }
  )

  # ---- DIFFERENTIAL PEAKS -----------------------------------------------------
  #' For each mark, find bins where it's present in one condition's combination
  #' but absent in the other, merge adjacent such bins into regions, and filter
  #' by differential score and merged region width.
  #' Replicates the exact "Galaxy-equivalent" pipeline from the user's Rmd:
  #' uses hmm$segments (pre-merged chromatin segments, NOT raw 200bp bins),
  #' filters by differential.score, combination inequality, and width, then
  #' counts mark presence via simple grepl() on the string combination columns.
  # ---- DIFFERENTIAL PEAKS (GENERIC VERSION) -----------------------------------
  
  compute_differential_peaks <- function(hmm, marks, score_thresh, width_thresh) {
    
    cat("\n--- compute_differential_peaks() called ---\n")
    cat("Object class:", paste(class(hmm), collapse = ", "), "\n")
    cat("Has $segments:", !is.null(hmm$segments), "\n")
    
    if (is.null(hmm$segments)) {
      showNotification(
        "hmm$segments not found — differential peaks requires merged segments.",
        type = "error", duration = 10
      )
      return(NULL)
    }
    
    segs_df <- as.data.frame(hmm$segments)
    
    cat("Segments:", nrow(segs_df), "rows\n")
    cat("Columns:", paste(colnames(segs_df), collapse = ", "), "\n")
    
    # --------------------------------------------------------------------------
    # Detect condition columns automatically
    # --------------------------------------------------------------------------
    combination_cols <- grep("^combination\\.", colnames(segs_df), value = TRUE)
    
    if (length(combination_cols) != 2) {
      showNotification(
        paste0(
          "Expected exactly 2 conditions, found ",
          length(combination_cols),
          ": ",
          paste(combination_cols, collapse = ", ")
        ),
        type = "error", duration = 12
      )
      return(NULL)
    }
    
    conditions <- sub("^combination\\.", "", combination_cols)
    cond1 <- conditions[1]
    cond2 <- conditions[2]
    
    cat("Conditions detected:", paste(conditions, collapse = ", "), "\n")
    
    # --------------------------------------------------------------------------
    # Required columns check
    # --------------------------------------------------------------------------
    required_cols <- c("differential.score", "width")
    missing_cols <- setdiff(required_cols, colnames(segs_df))
    
    if (length(missing_cols) > 0) {
      showNotification(
        paste0("Missing columns: ", paste(missing_cols, collapse = ", ")),
        type = "error", duration = 12
      )
      return(NULL)
    }
    
    # --------------------------------------------------------------------------
    # Extract condition values
    # --------------------------------------------------------------------------
    combo1 <- as.character(segs_df[[combination_cols[1]]])
    combo2 <- as.character(segs_df[[combination_cols[2]]])
    
    # --------------------------------------------------------------------------
    # Filter segments
    # --------------------------------------------------------------------------
    filtered <- segs_df[
      segs_df$differential.score >= score_thresh &
        combo1 != combo2 &
        segs_df$width >= width_thresh,
    ]
    
    if (nrow(filtered) == 0) {
      
      score_range <- paste0(
        round(min(segs_df$differential.score, na.rm = TRUE), 4),
        " to ",
        round(max(segs_df$differential.score, na.rm = TRUE), 4)
      )
      
      showNotification(
        paste0(
          "No segments passed filters (score≥", score_thresh,
          ", width≥", width_thresh, "). Range: ", score_range
        ),
        type = "warning", duration = 12
      )
      
      return(NULL)
    }
    
    f1 <- as.character(filtered[[combination_cols[1]]])
    f2 <- as.character(filtered[[combination_cols[2]]])
    
    # --------------------------------------------------------------------------
    # Count differential regions per mark
    # --------------------------------------------------------------------------
    results <- lapply(marks, function(m) {
      
      cond1_not_cond2 <- sum(
        grepl(m, f1, fixed = TRUE) &
          !grepl(m, f2, fixed = TRUE)
      )
      
      cond2_not_cond1 <- sum(
        !grepl(m, f1, fixed = TRUE) &
          grepl(m, f2, fixed = TRUE)
      )
      
      data.frame(
        mark = rep(m, 2),
        direction = c(
          paste0(cond1, "-not-", cond2),
          paste0(cond2, "-not-", cond1)
        ),
        n_regions = c(cond1_not_cond2, cond2_not_cond1),
        stringsAsFactors = FALSE
      )
    })
    
    out <- dplyr::bind_rows(results)
    
    attr(out, "total_filtered") <- nrow(filtered)
    attr(out, "conditions") <- conditions
    
    out
  }
  
  # ------------------------------------------------------------------------------
  # Reactive: compute differential peaks
  # ------------------------------------------------------------------------------
  
  diffpeaks_data <- eventReactive(input$run_diffpeaks, {
    
    req(rv$hmm, rv$marks)
    
    tryCatch({
      
      withProgress(message = "Computing differential peaks…", {
        
        compute_differential_peaks(
          hmm = rv$hmm,
          marks = rv$marks,
          score_thresh = input$diff_score_thresh,
          width_thresh = input$diff_width_thresh
        )
        
      })
      
    }, error = function(e) {
      
      showNotification(
        paste0("Error: ", conditionMessage(e)),
        type = "error", duration = 15
      )
      
      NULL
    })
  })
  
  # ------------------------------------------------------------------------------
  # Plot builder
  # ------------------------------------------------------------------------------
  
  build_diffpeaks_plot <- reactive({
    
    df <- diffpeaks_data()
    req(df)
    
    total <- attr(df, "total_filtered")
    if (is.null(total)) total <- sum(df$n_regions)
    
    conditions <- attr(df, "conditions")
    cond1 <- conditions[1]
    cond2 <- conditions[2]
    
    # order marks
    mark_order <- df %>%
      dplyr::group_by(mark) %>%
      dplyr::summarise(total = sum(n_regions), .groups = "drop") %>%
      dplyr::arrange(desc(total)) %>%
      dplyr::pull(mark)
    
    df$mark <- factor(df$mark, levels = rev(mark_order))
    
    # dynamic colors (works for any condition names)
    dirs <- unique(df$direction)
    
    fill_cols <- setNames(
      c("#d6604d", "#4393c3")[seq_along(dirs)],
      dirs
    )
    
    ggplot(df, aes(x = mark, y = n_regions, fill = direction)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = fill_cols) +
      labs(
        x = "Histone Mark",
        y = "Number of Regions",
        fill = NULL,
        title = "Pairwise Differential Peaks per Histone Mark",
        subtitle = paste0(
          cond1, " vs ", cond2,
          " | score≥", input$diff_score_thresh,
          ", width≥", input$diff_width_thresh,
          "bp | Total=", format(total, big.mark = ",")
        )
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position = "top",
        panel.grid.minor = element_blank()
      )
  })
  
  # ------------------------------------------------------------------------------
  # Outputs
  # ------------------------------------------------------------------------------
  
  output$diffpeaks_plot <- renderPlot({
    build_diffpeaks_plot()
  })
  
  output$dl_diffpeaks_plot <- downloadHandler(
    filename = function() {
      paste0("differential_peaks_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(file, build_diffpeaks_plot(), width = 10, height = 7, device = "pdf")
    }
  )
  
  output$dl_diffpeaks_xlsx <- downloadHandler(
    filename = function() {
      paste0("differential_peaks_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- diffpeaks_data()
      req(df)
      openxlsx::write.xlsx(df, file, overwrite = TRUE)
    }
  )

  # ---- GENE SET COMPARISON ---------------------------------------------------
  parse_gene_text <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(character(0))
    raw <- strsplit(txt, "[\r\n,;\t]+")[[1]]
    raw <- trimws(raw)
    raw[nchar(raw) > 0]
  }

  gsc_genes_a_parsed <- reactive(parse_gene_text(input$gsc_genes_a))
  gsc_genes_b_parsed <- reactive(parse_gene_text(input$gsc_genes_b))

  make_match_summary <- function(pasted, nms) {
    if (length(pasted) == 0) {
      return(tags$small(style = "color:#888", "No genes entered yet."))
    }
    matched   <- pasted[pasted %in% nms]
    unmatched <- setdiff(pasted, nms)
    tagList(
      tags$small(style = "color:#28a745",
        paste0("✓ ", length(matched), " of ", length(pasted), " matched.")),
      if (length(unmatched) > 0) {
        tags$div(tags$small(style = "color:#dc3545",
          paste0("✗ ", length(unmatched), " not found, e.g.: ",
                paste(head(unmatched, 5), collapse = ", "),
                if (length(unmatched) > 5) "…" else "")))
      }
    )
  }

  output$gsc_match_summary_a <- renderUI({
    req(rv$genes)
    make_match_summary(gsc_genes_a_parsed(), gene_names_reactive())
  })

  output$gsc_match_summary_b <- renderUI({
    req(rv$genes)
    make_match_summary(gsc_genes_b_parsed(), gene_names_reactive())
  })

  gsc_result <- eventReactive(input$run_genesetcompare, {
    req(rv$hmm, rv$genes)

    nms <- gene_names_reactive()
    sel_a <- gsc_genes_a_parsed(); sel_a <- sel_a[sel_a %in% nms]

    if (length(sel_a) == 0) {
      showNotification(
        "Gene set A needs at least one valid, matched gene.",
        type = "warning", duration = 6
      )
      return(NULL)
    }

    if (identical(input$gsc_mode, "permutation")) {

      res <- tryCatch({
        withProgress(message = "Building gene universe…", {
          nms_all <- gene_names_reactive()

          universe <- compute_all_gene_promoter_posteriors(
            hmm                 = rv$hmm,
            all_genes_gr        = rv$genes,
            gene_labels         = nms_all,
            upstream            = input$gsc_upstream,
            downstream          = input$gsc_downstream,
            average_replicates  = identical(input$gsc_replicate_mode, "average")
          )
          if (!is.null(universe$error)) return(list(error = universe$error))

          incProgress(0.3, detail = paste0("Running ", input$gsc_n_perm, " random draws…"))
          perm_res <- compute_gene_set_permutation(
            all_gene_scores = universe$gene_scores,
            post_cols       = universe$post_cols,
            genes_A_labels  = sel_a,
            n_perm          = input$gsc_n_perm,
            summary_stat    = input$gsc_summary_stat
          )
          if (!is.null(perm_res$error)) return(list(error = perm_res$error))
          perm_res$mode <- "permutation"
          perm_res
        })
      }, error = function(e) list(error = paste0("Error: ", conditionMessage(e))))

      if (!is.null(res$error)) {
        showNotification(res$error, type = "error", duration = 12)
        return(NULL)
      }
      return(res)

    } else {
      sel_b <- gsc_genes_b_parsed(); sel_b <- sel_b[sel_b %in% nms]
      if (length(sel_b) == 0) {
        showNotification(
          "Gene set B needs at least one valid, matched gene (or switch to the random background mode).",
          type = "warning", duration = 6
        )
        return(NULL)
      }

      genes_A_gr <- rv$genes[nms %in% sel_a]
      mcols(genes_A_gr)$gene_label <- nms[nms %in% sel_a]
      genes_B_gr <- rv$genes[nms %in% sel_b]
      mcols(genes_B_gr)$gene_label <- nms[nms %in% sel_b]

      res <- tryCatch({
        withProgress(message = "Comparing gene sets…", {
          compute_gene_set_posteriors(
            hmm                 = rv$hmm,
            genes_A_gr          = genes_A_gr,
            genes_B_gr          = genes_B_gr,
            upstream            = input$gsc_upstream,
            downstream          = input$gsc_downstream,
            average_replicates  = identical(input$gsc_replicate_mode, "average"),
            summary_stat        = input$gsc_summary_stat
          )
        })
      }, error = function(e) list(error = paste0("Error: ", conditionMessage(e))))

      if (!is.null(res$error)) {
        showNotification(res$error, type = "error", duration = 12)
        return(NULL)
      }
      res$mode <- "manual"
      res
    }
  })

  output$gsc_stats_table <- renderDT({
    res <- gsc_result()
    req(res, res$stats)
    df <- res$stats

    if (identical(res$mode, "permutation")) {
      df$p_value <- signif(df$empirical_p, 3)
      df$empirical_p <- NULL
    } else {
      df$p_value <- signif(df$p_value, 3)
    }
    df$FDR <- signif(df$FDR, 3)

    datatable(df, options = list(pageLength = 10, order = list(list(which(colnames(df) == "FDR") - 1, "asc"))),
             rownames = FALSE) %>%
      formatStyle("FDR", backgroundColor = styleInterval(c(0.01, 0.05),
                                                         c("#c6efce", "#ffeb9c", "white")))
  })

  output$gsc_summary_stat_label <- renderUI({
    res <- gsc_result()
    req(res)
    stat_used <- attr(res$stats, "summary_stat")
    if (is.null(stat_used)) return(NULL)
    tags$small(style = "color:#555",
      paste0("Comparison statistic: ", tools::toTitleCase(stat_used),
            " (delta and stat_A/stat_B columns use this; mean and median are both always shown where applicable)."))
  })

  build_gsc_plot <- reactive({
    res <- gsc_result()
    req(res, res$long_scores)

    df <- res$long_scores
    df$facet_label <- paste0(df$mark, " (", df$condition, ")")

    if (identical(res$mode, "permutation")) {
      group_colors <- c("A" = "#4393c3", "Random" = "#999999")
      df$group <- factor(df$group, levels = c("A", "Random"))
      subtitle_txt <- paste0(res$n_A, " genes in set A  |  ", res$n_A,
                            " genes in one representative random draw  |  ",
                            "empirical p-value uses ", nrow(res$perm_long) / length(unique(res$perm_long$mark)),
                            " random draws total (see stats table)")
    } else {
      group_colors <- c("A" = "#4393c3", "B" = "#d6604d")
      subtitle_txt <- paste0("Set A: ", res$n_genes_A, " genes  |  Set B: ", res$n_genes_B, " genes")
    }

    ggplot(df, aes(x = group, y = posterior, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
      scale_fill_manual(values = group_colors) +
      facet_wrap(~ facet_label, scales = "free_y") +
      labs(
        x = "Gene set", y = "Mean posterior probability over promoter",
        fill = "Gene set",
        title    = "Posterior probability by gene set and mark",
        subtitle = subtitle_txt
      ) +
      theme_bw(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "#34495e"),
        strip.text       = element_text(colour = "white", face = "bold"),
        legend.position  = "top"
      )
  })

  output$gsc_boxplot <- renderPlot({ build_gsc_plot() })

  output$dl_gsc_plot <- downloadHandler(
    filename = function() paste0("gene_set_comparison_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggsave(file, plot = build_gsc_plot(), width = 12, height = 8, device = "pdf")
    }
  )

  output$dl_gsc_xlsx <- downloadHandler(
    filename = function() paste0("gene_set_comparison_data_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      res <- gsc_result()
      req(res)
      if (identical(res$mode, "permutation")) {
        openxlsx::write.xlsx(
          list(stats = res$stats, per_gene_boxplot = res$long_scores, permutation_draws = res$perm_long),
          file, overwrite = TRUE
        )
      } else {
        openxlsx::write.xlsx(
          list(stats = res$stats, per_gene = res$long_scores),
          file, overwrite = TRUE
        )
      }
    }
  )

  output$gsc_gene_table <- renderDT({
    res <- gsc_result()
    req(res, res$long_scores)

    df <- res$long_scores[res$long_scores$group == "A",
                          c("gene", "mark", "condition", "posterior")]
    df$posterior <- round(df$posterior, 4)
    colnames(df)[colnames(df) == "posterior"] <- "mean_posterior_set_A"
    datatable(df, extensions = "Buttons",
             options = list(pageLength = 15, dom = "Bfrtip", buttons = c("csv", "excel"),
                            order = list(list(which(colnames(df) == "mean_posterior_set_A") - 1, "desc"))),
             rownames = FALSE, filter = "top")
  })
}
# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)
