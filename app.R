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
# BUILD STAMP + MASKING GUARD
# =============================================================================
# randomForest (and a few other packages) export their own margin(). Attached
# after ggplot2 it wins the search path, so a bare margin() call resolves to
# randomForest::margin(x, observed, ...) and fails with:
#     argument "observed" is missing, with no default
# Every margin() call in this file is already written as ggplot2::margin().
# This binding is belt-and-braces in case a bare one is ever reintroduced.
margin <- ggplot2::margin

APP_BUILD <- "2026-09-21f / per-condition colour pickers"
message("---------------------------------------------------------------")
message("ChromstaR Viewer build: ", APP_BUILD)
message("margin() resolves to: ", environmentName(environment(margin)),
        "   (must be 'ggplot2')")
message("---------------------------------------------------------------")

# =============================================================================
# PLOT CUSTOMIZATION HELPER FUNCTIONS (NEW - ADD TEXT & STYLING CONTROLS)
# =============================================================================

#' Create a customizable theme based on user-selected text sizes and styling
build_custom_theme <- function(
    base_size = 13,
    title_size = 14,
    axis_title_size = 12,
    axis_text_size = 10,
    legend_title_size = 11,
    legend_text_size = 10,
    strip_text_size = 11,
    line_width = 0.8) {
  
  list(
    theme_bw(base_size = base_size),
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 8)),
      plot.subtitle = element_text(size = axis_title_size - 1, hjust = 0.5, margin = ggplot2::margin(b = 8)),
      axis.title.x = element_text(size = axis_title_size, face = "bold", margin = ggplot2::margin(t = 10)),
      axis.title.y = element_text(size = axis_title_size, face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text.x = element_text(size = axis_text_size, angle = 45, hjust = 1),
      axis.text.y = element_text(size = axis_text_size),
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      strip.text = element_text(size = strip_text_size, colour = "white", face = "bold"),
      strip.background = element_rect(fill = "#34495e"),
      panel.grid.minor = element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  )
}

#' Apply one tab's "Plot Appearance" slider settings to a finished ggplot.
#' Text sizes and gridlines are layered on with theme(); line width, point size
#' and transparency are geom-level parameters, so they are written directly into
#' the plot's existing layers. `prefix` selects which tab's sliders to read,
#' e.g. "mg_" for the metagene tab.
apply_plot_customization <- function(p, input, prefix = "") {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)

  gv <- function(suffix, default) {
    v <- input[[paste0(prefix, suffix)]]
    if (is.null(v) || length(v) != 1) return(default)
    if (is.numeric(v) && !is.finite(v)) return(default)
    v
  }
  title_size      <- gv("plot_title_size",      14)
  axis_title_size <- gv("plot_axis_title_size", 12)
  axis_text_size  <- gv("plot_axis_text_size",  10)
  legend_size     <- gv("plot_legend_size",     10)
  line_width      <- gv("plot_line_width",     0.8)
  point_size      <- gv("plot_point_size",       2)
  alpha_val       <- gv("plot_alpha",          0.2)
  show_grid       <- gv("plot_show_grid",     TRUE)

  # --- geom-level settings: these cannot be reached through theme() ---
  lw_arg <- if (utils::packageVersion("ggplot2") >= "3.4.0") "linewidth" else "size"
  for (i in seq_along(p$layers)) {
    g <- class(p$layers[[i]]$geom)[1]
    if (g %in% c("GeomLine", "GeomPath", "GeomStep", "GeomSmooth", "GeomSegment")) {
      p$layers[[i]]$aes_params[[lw_arg]] <- line_width
    }
    if (g %in% c("GeomPoint", "GeomJitter")) {
      p$layers[[i]]$aes_params$size <- point_size
    }
    # only shaded bands - deliberately NOT tiles/rects, so heatmaps stay readable
    if (g %in% c("GeomRibbon", "GeomArea")) {
      p$layers[[i]]$aes_params$alpha <- alpha_val
    }
  }

  # --- text sizes: preserve any x-axis label rotation the plot already set ---
  ang <- p$theme$axis.text.x$angle
  hj  <- p$theme$axis.text.x$hjust
  p <- p + ggplot2::theme(
    plot.title   = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
    axis.title.x = ggplot2::element_text(size = axis_title_size, face = "bold"),
    axis.title.y = ggplot2::element_text(size = axis_title_size, face = "bold"),
    axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = ang, hjust = hj),
    axis.text.y  = ggplot2::element_text(size = axis_text_size),
    legend.text  = ggplot2::element_text(size = legend_size),
    legend.title = ggplot2::element_text(size = legend_size + 1, face = "bold")
  )

  if (!isTRUE(show_grid)) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank())
  }
  p
}

#' Run `expr` and, if it fails, print a full call stack to the R console and
#' surface the innermost failing calls in the UI message. Used to pin down
#' errors that come from a function masked by another attached package.
with_stack_trace <- function(expr, label = "") {
  trace_rows <- NULL
  tryCatch(
    withCallingHandlers(
      expr,
      error = function(e) {
        skip <- c(".handleSimpleError", "h", "stop", "signalCondition", "try",
                  "tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch",
                  "withCallingHandlers", "with_stack_trace",
                  "..stacktraceon..", "..stacktraceoff..")
        rows <- character(0)
        for (i in seq_len(sys.nframe())) {
          cl <- tryCatch(sys.call(i), error = function(...) NULL)
          if (is.null(cl) || !is.call(cl)) next
          hd <- tryCatch(paste(deparse(cl[[1]]), collapse = ""),
                         error = function(...) "?")
          if (grepl("^function", hd)) next        # anonymous handler frames
          if (hd %in% skip) next
          fn  <- tryCatch(sys.function(i), error = function(...) NULL)
          env <- if (is.function(fn))
                   tryCatch(environmentName(environment(fn)),
                            error = function(...) "") else ""
          txt <- substr(paste(deparse(cl), collapse = " "), 1, 110)
          rows <- c(rows, paste0(txt, if (nzchar(env)) paste0("   <", env, ">") else ""))
        }
        trace_rows <<- rows
      }
    ),
    error = function(e) {
      message("\n=== ChromstaR Viewer traceback (", label, ") ===")
      message(conditionMessage(e))
      if (!is.null(trace_rows))
        message(paste(sprintf("%2d. %s", seq_along(trace_rows), trace_rows),
                      collapse = "\n"))
      message("=== end traceback ===\n")
      hint <- if (!is.null(trace_rows))
                paste(sprintf("%d) %s", seq_along(utils::tail(trace_rows, 8)),
                              utils::tail(trace_rows, 8)), collapse = "    ")
              else ""
      stop(paste0(conditionMessage(e),
                  "   ||  FAILING CALLS (innermost last, <package> in angle brackets):  ",
                  hint, "  ||  (full traceback also printed to the R console)"),
           call. = FALSE)
    }
  )
}

#' Per-mark domain counts from a chromstaR $frequencies table (genome-wide,
#' unfiltered). Mirrors the "Genome-wide Domain-Level Analysis" section of the
#' reference Rmd: for each mark, sum the `domains` column over the combination
#' pairs where that mark is present in one condition and absent in the other.
#' `cond_b` is the "second" condition, so gained = present in cond_b only.
compute_mark_ranking_freq <- function(freq_df, marks, cond_a, cond_b) {
  if (is.null(freq_df)) return(NULL)
  freq_df <- as.data.frame(freq_df)
  col_a <- paste0("combination.", cond_a)
  col_b <- paste0("combination.", cond_b)
  if (!all(c(col_a, col_b, "domains") %in% colnames(freq_df))) return(NULL)

  ca <- as.character(freq_df[[col_a]])
  cb <- as.character(freq_df[[col_b]])
  dm <- suppressWarnings(as.numeric(freq_df$domains))
  total <- sum(dm, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NULL)

  out <- do.call(rbind, lapply(marks, function(m) {
    in_a <- grepl(m, ca, fixed = TRUE)
    in_b <- grepl(m, cb, fixed = TRUE)
    gained <- sum(dm[!in_a &  in_b], na.rm = TRUE)   # present in cond_b only
    lost   <- sum(dm[ in_a & !in_b], na.rm = TRUE)   # present in cond_a only
    data.frame(mark          = m,
               gained        = gained,
               lost          = lost,
               constant      = sum(dm[in_a & in_b], na.rm = TRUE),
               total_changed = gained + lost,
               stringsAsFactors = FALSE)
  }))
  out$total_domains <- total
  out$pct_changed   <- round(out$total_changed / total * 100, 2)
  out$pct_gained    <- round(out$gained        / total * 100, 2)
  out$pct_lost      <- round(out$lost          / total * 100, 2)
  out$net           <- out$gained - out$lost
  out$cond_a        <- cond_a
  out$cond_b        <- cond_b
  out$direction     <- ifelse(out$net > 0, paste0("Gained in ", cond_b),
                                           paste0("Lost in ",   cond_b))
  out <- out[order(-out$total_changed), ]
  rownames(out) <- NULL
  out
}

#' Shared horizontal ranked-bar chart for the two "marks ranked by change"
#' figures. `rank_df` needs: mark, total_changed, pct_changed, direction — and
#' optionally `panel`, which facets when it has more than one level. Marks are
#' ordered by total change summed across panels so the ordering is comparable.
mark_ranking_ggplot <- function(rank_df, title, subtitle,
                                ylab  = "Total Domains Changed",
                                style = "stacked",
                                cond_colors = NULL) {
  if (is.null(rank_df) || nrow(rank_df) == 0) return(NULL)
  if (!"panel" %in% colnames(rank_df)) rank_df$panel <- ""

  # Marks ordered by total change, summed across panels so the order is shared.
  ord <- stats::aggregate(total_changed ~ mark, data = rank_df, FUN = sum)
  ord <- ord[order(ord$total_changed), ]
  lev <- as.character(ord$mark)

  # One row per mark PER DIRECTION, so the bar is never a single colour
  # standing for a total that mixes gains and losses.
  long <- rbind(
    data.frame(mark  = as.character(rank_df$mark),
               panel = rank_df$panel,
               part  = paste0("Gained in ", rank_df$cond_b),
               cond  = as.character(rank_df$cond_b),   # mark present in cond_b
               n     = rank_df$gained,
               pct   = rank_df$pct_gained,
               stringsAsFactors = FALSE),
    data.frame(mark  = as.character(rank_df$mark),
               panel = rank_df$panel,
               part  = paste0("Lost in ", rank_df$cond_b),
               cond  = as.character(rank_df$cond_a),   # mark present in cond_a
               n     = rank_df$lost,
               pct   = rank_df$pct_lost,
               stringsAsFactors = FALSE))
  long$mark <- factor(long$mark, levels = lev)

  # Each half of the bar is "the mark is present in THIS condition", so it is
  # coloured with that condition's colour — the same one every other figure
  # in the app uses. Falls back to the red/blue default when none is set.
  key  <- unique(long[, c("part", "cond")])
  cols <- stats::setNames(vapply(seq_len(nrow(key)), function(i) {
    cn <- key$cond[i]
    if (!is.null(cond_colors) && !is.na(cn) && cn %in% names(cond_colors)) {
      unname(cond_colors[[cn]])
    } else if (grepl("^Gained", key$part[i])) "#d6604d" else "#4393c3"
  }, character(1)), key$part)
  lv <- key$part
  long$part <- factor(long$part, levels = lv)

  tot <- rank_df
  tot$mark <- factor(as.character(tot$mark), levels = lev)

  if (identical(style, "diverging")) {
    # Losses left of zero, gains right of zero: direction is unmistakable and
    # no number on the chart stands for a mixture of the two.
    long$n_signed <- ifelse(grepl("^Lost", as.character(long$part)), -long$n, long$n)
    p <- ggplot(long, aes(x = mark, y = n_signed, fill = part)) +
      geom_col() +
      geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
      geom_text(aes(label = ifelse(n == 0, "",
                                   paste0(formatC(pct, format = "f", digits = 2), "%")),
                    hjust = ifelse(n_signed < 0, 1.12, -0.12)),
                size = 3.1) +
      scale_y_continuous(labels = function(v) scales::comma(abs(v)),
                         expand = expansion(mult = c(0.18, 0.18)))
    ylab_use <- paste0(ylab, "   (left = lost, right = gained)")
  } else {
    # Stacked: total bar length still ranks the marks, but the gained and lost
    # parts are drawn separately. The label is the combined percentage.
    p <- ggplot(long, aes(x = mark, y = n, fill = part)) +
      geom_col() +
      geom_text(data = tot,
                aes(x = mark, y = total_changed,
                    label = paste0(formatC(pct_changed, format = "f", digits = 2), "%")),
                inherit.aes = FALSE, hjust = -0.12, size = 3.4) +
      scale_y_continuous(labels = scales::comma,
                         expand = expansion(mult = c(0, 0.18)))
    ylab_use <- paste0(ylab, "   (bar = gained + lost)")
  }

  p <- p +
    coord_flip(clip = "off") +
    scale_fill_manual(values = cols) +
    labs(title = title, subtitle = subtitle,
         x = "Histone Mark", y = ylab_use, fill = "") +
    theme_bw(base_size = 13) +
    theme(legend.position = "top")

  if (length(unique(rank_df$panel)) > 1) p <- p + facet_wrap(~ panel, scales = "free_x")
  p
}

#' Stable, syntactically safe Shiny input id for a condition name.
cond_colour_id <- function(cond) paste0("cond_col_", gsub("[^A-Za-z0-9]", "_", cond))

#' A colour picker that needs no extra package: a standard Shiny text input
#' switched to the browser's native <input type="color">. Shiny still reads it
#' as a character value like "#1f77b4", so it works like any other input.
condition_colour_input <- function(inputId, label, value) {
  ti <- shiny::textInput(inputId, label, value = value)
  for (i in seq_along(ti$children)) {
    ch <- ti$children[[i]]
    if (!is.null(ch) && !is.null(ch$name) && identical(ch$name, "input")) {
      ti$children[[i]]$attribs$type  <- "color"
      ti$children[[i]]$attribs$style <- "width:66px; height:36px; padding:2px; cursor:pointer;"
    }
  }
  ti
}

#' Get condition colors
get_condition_colors <- function(conditions) {
  color_palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf"
  )
  stats::setNames(color_palette[1:length(conditions)], conditions)
}

#' Get condition line types
get_condition_linetypes <- function(conditions) {
  ltypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  stats::setNames(ltypes[1:length(conditions)], conditions)
}



# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load a ChromstaR combinedMultiHMM from .RData or .rds and immediately
#' slim it to only the four slots the app actually queries, freeing the
#' rest of RAM before returning. This keeps memory usage on shinyapps.io
#' well within the 1 GB free-tier limit.
load_chromstar_object <- function(path) {

  ext <- tolower(tools::file_ext(path))

  if (ext == "rds") {
    obj <- readRDS(path)
  } else {
    env  <- new.env()
    load(path, envir = env)
    objs <- ls(env)
    obj  <- NULL
    for (nm in objs) {
      candidate <- get(nm, envir = env)
      if (inherits(candidate, c("combinedMultiHMM", "multiHMM", "uniHMM"))) {
        obj <- candidate
        break
      }
    }
    if (is.null(obj)) obj <- get(ls(env)[1], envir = env)
    rm(env)
    gc()
  }

  # Slots the app needs:
  #   bins     -> RPKM + posteriors + combinations (all analysis tabs)
  #   segments -> merged chromatin segments (Differential Peaks tab)
  #   info     -> mark / condition metadata (Load Data detection)
  #   hmms     -> fallback for mark/condition detection on some object types
  #   frequencies -> per-combination domain counts, used by the genome-wide
  #                  "marks ranked by change" figure on the Differential Peaks tab
  needed <- c("bins", "segments", "info", "hmms", "frequencies")
  slim   <- list()
  for (slot in needed) {
    val <- tryCatch(obj[[slot]], error = function(e) NULL)
    if (!is.null(val)) slim[[slot]] <- val
  }
  class(slim) <- class(obj)

  rm(obj)
  gc()

  slim
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
#' `ratio_mode` decides HOW log(observed/expected) is formed:
#'   "means"  — (default, and what the Galaxy chromstaR profile does) average the
#'              raw signal of every bin in a slot, zero bins included, then take
#'              one log ratio of that mean against the genome-wide mean:
#'                  log( mean(signal in slot) / mean(signal genome-wide) )
#'   "legacy" — take log(signal/expected) for each bin separately and average the
#'              logs. That is a geometric mean, so bins with zero signal have to
#'              be discarded (log(0) is undefined). Discarding them biases every
#'              slot upward, and the bias is largest exactly where a mark is
#'              genuinely ABSENT (that is where most bins are zero) — so real
#'              depletion is flattened away. Kept only to reproduce old figures.
compute_enrichment_profile <- function(hmm, genes_gr, marks, conditions,
                                       upstream   = 2000,
                                       downstream = 2000,
                                       n_bins     = 40,
                                       ratio_mode = "means") {

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
      expected[[key]]   <- if (identical(ratio_mode, "legacy")) {
        mean(vals[vals > 0], na.rm = TRUE)
      } else {
        mean(vals, na.rm = TRUE)   # zero bins belong in the baseline
      }
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
  gene_slots   <- list()

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

    # A 200-250bp bin can overlap a zone while its midpoint sits outside it.
    # Those rows produced fractions slightly below 0 or above 1, i.e. slots
    # beyond the intended [-1, 2] axis (x = 2.025, 2.05 ...) that the plot then
    # clipped but the Excel export still carried, each averaged over a handful
    # of bins. Drop them at source instead.
    keep  <- !is.na(frac) & frac >= 0 & frac <= 1
    bi    <- bi[keep]
    g_i   <- g_i[keep]
    x_val <- x_val[keep]
    if (length(bi) == 0) next

    slot_idx <- round(x_val * n_bins)

    # Which genes land in which slot depends only on geometry, not on the mark
    # or condition, so record it once here and join it on at the end.
    gene_slots[[zname]] <- data.frame(slot = slot_idx, gene_idx = g_i)

    for (cond in conditions) {
      for (mark in marks) {
        key <- paste(cond, mark)
        sig_vec <- signal_mat[[key]]
        if (is.null(sig_vec)) next
        exp_val <- expected[[key]]
        if (is.na(exp_val) || exp_val <= 0) next

        sig <- sig_vec[bi]

        if (identical(ratio_mode, "legacy")) {
          sig[sig <= 0] <- NA          # a geometric mean cannot take log(0)
          val <- log(sig / exp_val)
        } else {
          val <- sig                   # keep raw signal; the log is taken after
        }                              # averaging, in the summarise() below

        prof_key <- paste(cond, mark, zname)
        all_profiles[[prof_key]] <- data.frame(
          slot      = slot_idx,
          value     = val,
          expected  = exp_val,
          mark      = mark,
          condition = cond,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(all_profiles) == 0) return(NULL)

  out <- bind_rows(all_profiles) %>%
    group_by(slot, mark, condition) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      expected   = dplyr::first(expected),
      n_bins_avg = sum(!is.na(value)),
      n_genes    = n(),
      .groups    = "drop"
    )

  # "means": one log ratio per slot, taken AFTER averaging the raw signal, so
  # empty bins pull the mean down the way they should and real depletion
  # survives. "legacy": mean_value is already a mean of per-bin logs.
  out$mean_log_ratio <- if (identical(ratio_mode, "legacy")) {
    out$mean_value
  } else {
    log(out$mean_value / out$expected)
  }

  # Attach the genes behind every position: how many distinct genes contribute,
  # and their names. The list is capped so one cell cannot blow past Excel's
  # 32,767-character limit when the whole annotation is in scope.
  if (length(gene_slots) > 0) {
    g_all <- mcols(genes_gr)$gene_name
    if (is.null(g_all)) g_all <- as.character(seq_along(genes_gr))

    gmap <- bind_rows(gene_slots)
    gmap$gene <- g_all[gmap$gene_idx]
    gmap <- gmap[!is.na(gmap$gene), c("slot", "gene")]
    gmap <- gmap[!duplicated(gmap), ]

    # Total distinct genes contributing anywhere in the profile. Reported next
    # to the per-position counts so a figure caption states a real sample size
    # instead of just the single best-covered position.
    n_genes_total <- dplyr::n_distinct(gmap$gene)

    gmap <- gmap %>%
      group_by(slot) %>%
      summarise(
        n_genes_distinct = dplyr::n(),
        genes = {
          g <- sort(gene)
          if (length(g) > 100) {
            paste0(paste(g[1:100], collapse = ";"), ";(+", length(g) - 100, " more)")
          } else {
            paste(g, collapse = ";")
          }
        },
        .groups = "drop"
      )

    out <- dplyr::left_join(out, gmap, by = "slot")
    out$n_genes_total <- n_genes_total
  }

  out %>% filter(is.finite(mean_log_ratio))
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

  rel_pos <- (bin_mid - region_start) / region_width

  # STRAND. The window itself is placed strand-aware above (promoters() for TSS,
  # the mirrored construction for TES), but rel_pos is measured left-to-right in
  # genomic coordinates. For a gene on the minus strand the 5' end is on the
  # RIGHT, so the fraction has to be flipped — otherwise minus-strand genes are
  # averaged in backwards and the profile comes out mirror-symmetric (a promoter
  # mark then appears at BOTH ends of the gene body instead of only at the TSS).
  # The Enrichment Profile tab already does this; the metagene did not.
  is_minus <- as.character(strand(ref_points))[gi] == "-"
  rel_pos  <- ifelse(is_minus, 1 - rel_pos, rel_pos)

  # Bins are 200–250bp wide, so a bin can overlap the window while its midpoint
  # falls outside it. Clamping those into the first/last slot piles flanking
  # signal onto the two end points; drop them instead.
  inside <- rel_pos >= 0 & rel_pos <= 1
  bi     <- bi[inside]
  gi     <- gi[inside]
  rel_pos <- rel_pos[inside]
  if (length(bi) == 0) return(NULL)

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


#' Per-gene contributions to a metagene profile.
#'
#' Same windows and the same strand-aware positioning as compute_metagene(),
#' but WITHOUT averaging across genes: returns one row per gene per position.
#' This is what lets you ask whether an average profile reflects most genes or
#' is carried by a handful of very strong ones.
#'
#' Returns a data.frame: gene, bin_idx, signal (mean RPKM over the bins of that
#' gene falling in that position slot).
compute_gene_contributions <- function(hmm, genes_gr, mark, condition,
                                       mode       = "TSS",
                                       upstream   = 2000,
                                       downstream = 2000,
                                       n_bins     = 100) {

  bins        <- hmm$bins
  bins_df     <- expand_bins_df(hmm)
  bin_mid_all <- (bins_df$start + bins_df$end) / 2

  rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
  pattern <- paste0("counts.rpkm.", mark, ".", condition)
  cols    <- grep(pattern, rpkm_cols_all, value = TRUE, fixed = TRUE)
  if (length(cols) == 0) return(NULL)
  sig_vec <- rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)

  if (mode == "TSS") {
    ref_points <- promoters(genes_gr, upstream = upstream, downstream = downstream)
  } else if (mode == "TES") {
    str_g <- as.character(strand(genes_gr))
    tes   <- ifelse(str_g == "-", start(genes_gr), end(genes_gr))
    win_start <- ifelse(str_g == "-", tes - downstream, tes - upstream)
    win_end   <- ifelse(str_g == "-", tes + upstream,   tes + downstream)
    win_start <- pmax(1, win_start)
    ref_points <- GRanges(seqnames = seqnames(genes_gr),
                          ranges   = IRanges(win_start, win_end),
                          strand   = strand(genes_gr))
  } else {
    ref_points <- genes_gr
  }

  g_names <- mcols(genes_gr)$gene_name
  if (is.null(g_names)) g_names <- as.character(seq_along(genes_gr))
  mcols(ref_points) <- NULL
  mcols(ref_points)$gene_name <- g_names

  ov <- findOverlaps(bins, ref_points)
  if (length(ov) == 0) return(NULL)

  bi <- queryHits(ov)
  gi <- subjectHits(ov)

  rel_pos  <- (bin_mid_all[bi] - start(ref_points)[gi]) / width(ref_points)[gi]
  is_minus <- as.character(strand(ref_points))[gi] == "-"
  rel_pos  <- ifelse(is_minus, 1 - rel_pos, rel_pos)

  inside  <- rel_pos >= 0 & rel_pos <= 1
  bi      <- bi[inside]; gi <- gi[inside]; rel_pos <- rel_pos[inside]
  if (length(bi) == 0) return(NULL)

  slot_idx <- pmin(pmax(ceiling(rel_pos * n_bins), 1L), n_bins)

  data.frame(
    gene    = mcols(ref_points)$gene_name[gi],
    bin_idx = slot_idx,
    signal  = sig_vec[bi],
    stringsAsFactors = FALSE
  ) %>%
    group_by(gene, bin_idx) %>%
    summarise(signal = mean(signal, na.rm = TRUE), .groups = "drop")
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
    title = tags$span(
      tags$div("ChromstaR Viewer"),
      tags$div(style = "font-size:11px; font-weight:normal; color:#dce6ec;", "by Janan Gawra")
    ),
    tags$li(
      class = "dropdown",
      style = "padding: 4px 15px;",
      tags$img(src = "ihpe_logo.jpg", height = "42px")
    ),
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
      menuItem("Gene Contributions", tabName = "genecontrib", icon = icon("layer-group")),
      menuItem("Enrichment Profile", tabName = "enrichment", icon = icon("chart-line")),
      menuItem("Region Browser",  tabName = "browser",    icon = icon("dna")),
      menuItem("Differential Peaks", tabName = "diffpeaks", icon = icon("chart-simple")),
      menuItem("Gene Set Comparison", tabName = "genesetcompare", icon = icon("dna")),
      menuItem("Data Table",      tabName = "table",      icon = icon("table")),
      menuItem("About / Help",    tabName = "about",      icon = icon("circle-info"))
    ),
    tags$div(
      style = "position:absolute; bottom:0; width:100%; padding:10px 15px; color:#b8c7ce; font-size:11px; text-align:center;",
      HTML("&copy; 2026 Janan Gawra &mdash; IHPE")
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
      .skin-blue .main-header .logo {
        height: 58px;
        line-height: 1.1;
        padding-top: 8px;
        white-space: normal;
        overflow: visible;
      }
      .skin-blue .main-header .navbar {
        min-height: 58px;
      }
      .main-header { max-height: 58px; }
      .main-header .navbar > .dropdown { display: flex; align-items: center; }
      .main-sidebar, .left-side { top: 58px; }
      .content-wrapper, .main-footer { margin-top: 58px; }

      /* Keep the header fixed/visible while scrolling */
      .main-header { position: fixed; width: 100%; top: 0; z-index: 1030; }
      .main-header .navbar { position: relative; }
      body { padding-top: 0; }
    "))),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1 — LOAD DATA
      # -----------------------------------------------------------------------
      tabItem(tabName = "load",
        fluidRow(
          column(12,
            tags$div(style = "text-align:center; padding: 10px 0 20px 0;",
              tags$img(src = "ihpe_logo.jpg", height = "70px"),
              tags$p(style = "color:#888; font-size:12px; margin-top:8px;",
                HTML("&copy; 2026 Janan Gawra &mdash; IHPE (Institut des Sciences de l'Evolution, UMR 5244 CNRS-UPVD)"))
            )
          )
        ),
        fluidRow(
          box(title = "Load ChromstaR Object — Stage A (.RData / .rds)", width = 6, status = "primary",
            textInput("stage_a_label", "Label for this dataset (used in plots)", value = "Stage A"),
            fileInput("chromstar_file", "Upload ChromstaR .RData or .rds file",
                      accept = c(".RData", ".rda", ".Rdata", ".rds")),
            verbatimTextOutput("chromstar_summary")
          ),
          box(title = "Load ChromstaR Object — Stage B (optional)", width = 6, status = "primary",
            tags$small(style = "color:#888",
              "Optional second object — e.g. a different developmental stage or physiological condition, loaded as an independent ChromstaR run. Once loaded, every analysis tab gains a “stages to compare” option so Stage A and Stage B can be viewed side by side."),
            br(), br(),
            textInput("stage_b_label", "Label for this dataset (used in plots)", value = "Stage B"),
            fileInput("chromstar_file_b", "Upload ChromstaR .RData or .rds file",
                      accept = c(".RData", ".rda", ".Rdata", ".rds")),
            verbatimTextOutput("chromstar_summary_b")
          )
        ),
        fluidRow(
          box(title = "Load Gene Annotation (TSV)", width = 12, status = "primary",
            tags$small(style = "color:#888",
              "Shared across both stages — both ChromstaR objects are assumed to use the same genome/gene annotation."),
            br(), br(),
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
              column(3,
                h4(textOutput("stage_a_label_header")),
                h5("Histone Marks"),
                uiOutput("marks_ui"),
                h5("Conditions"),
                uiOutput("conditions_ui")
              ),
              column(3,
                h4(textOutput("stage_b_label_header")),
                h5("Histone Marks"),
                uiOutput("marks_ui_b"),
                h5("Conditions"),
                uiOutput("conditions_ui_b")
              ),
              column(3,
                h4("Genome Info — Stage A"),
                verbatimTextOutput("genome_info")
              ),
              column(3,
                h4("Genome Info — Stage B"),
                verbatimTextOutput("genome_info_b")
              )
            )
          )
        ),
        fluidRow(
          box(title = "Condition / Life-cycle Order", width = 12, status = "warning",
            tags$small(style = "color:#888",
              HTML("Sets which conditions are used and the order they appear in <b>every</b> plot legend, facet strip, and results table across the whole app. Drag the boxes below to arrange them into a biologically meaningful order &mdash; e.g. <b>TwoDO &rarr; Meta &rarr; PNA &rarr; PA</b> to follow the life cycle &mdash; instead of the default alphabetical order. Colours follow this order too, so the first condition always gets the first colour. <b>Removing a condition (the &times; on its box) hides it everywhere in the app</b>; \"Reset to detected order\" brings them all back.")),
            br(), br(),
            uiOutput("condition_order_ui"),
            actionButton("reset_condition_order", "Reset to detected order",
                        class = "btn-sm btn-default"),
            hr(),
            h5("Condition colours"),
            tags$small(style = "color:#888",
              "Pick a colour for each condition. That colour is then used for the condition everywhere in the app — metagene, region browser, differential peaks, gene set comparison and the ranked figures — so one condition looks the same in every figure you export."),
            br(), br(),
            uiOutput("condition_colour_ui"),
            actionButton("reset_condition_colours", "Reset colours",
                        class = "btn-sm btn-default"),
            br(), br(),
            uiOutput("condition_order_preview")
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
                         choices = c("All genes in GTF"    = "all",
                                     "Only selected genes" = "selected"),
                         selected = "all"),
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
            h5("Stages to compare"),
            uiOutput("stage_selector_meta"),
            hr(),
            h5("Display mode"),
            radioButtons("compare_mode_meta", "Condition comparison",
                         choices = c("Side-by-side" = "facet",
                                     "Overlay"       = "overlay"),
                         selected = "facet"),
            checkboxInput("smooth_metagene",
                         "Smooth curve (LOESS) — recommended for many genes",
                         value = TRUE),
            hr(),
            # ========== NEW: PLOT CUSTOMIZATION PANEL (METAGENE) ==========
            box(
              title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
              width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
              h5("Text Sizes"),
              sliderInput("mg_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
              sliderInput("mg_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
              sliderInput("mg_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
              sliderInput("mg_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
              hr(),
              h5("Lines & Visual Elements"),
              sliderInput("mg_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
              sliderInput("mg_plot_point_size", "Point/marker size", min = 1, max = 5, value = 2, step = 0.5),
              sliderInput("mg_plot_alpha", "Transparency of fills", min = 0.1, max = 1, value = 0.2, step = 0.1),
              hr(),
              h5("Other Elements"),
              checkboxInput("mg_plot_show_grid", "Show gridlines", value = TRUE),
              tags$small(style = "color:#888; display:block; margin-top:10px;",
                "Adjust and click 'Compute Profile' to apply."
              )
            ),
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
      # TAB 2.5 — GENE CONTRIBUTIONS
      # -----------------------------------------------------------------------
      tabItem(tabName = "genecontrib",
        fluidRow(
          box(title = "Gene Contribution Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Breaks a metagene profile back down into the individual genes behind it: which genes carry the signal, how concentrated it is, and what the average line hides."),
            hr(),
            uiOutput("gc_stage_ui"),
            uiOutput("gc_mark_ui"),
            uiOutput("gc_condition_ui"),
            radioButtons("gc_mode", "Reference",
                         choices = c("TSS" = "TSS", "TES" = "TES",
                                     "Gene body" = "gene_body"),
                         selected = "TSS"),
            conditionalPanel(
              condition = "input.gc_mode != 'gene_body'",
              numericInput("gc_upstream",   "Upstream (bp)",   value = 2000, min = 0, step = 100),
              numericInput("gc_downstream", "Downstream (bp)", value = 2000, min = 0, step = 100)
            ),
            numericInput("gc_n_bins", "Number of bins", value = 40, min = 10, max = 200, step = 5),
            hr(),
            h5("Gene scope"),
            radioButtons("gc_gene_scope", NULL,
                         choices = c("All genes in annotation" = "all",
                                     "Only selected genes"     = "selected"),
                         selected = "all"),
            conditionalPanel(
              condition = "input.gc_gene_scope == 'selected'",
              uiOutput("gene_selector_gc")
            ),
            hr(),
            numericInput("gc_top_n", "Genes to show in heatmap (top N by signal)",
                         value = 200, min = 10, max = 5000, step = 50),
            checkboxInput("gc_log_scale", "Log-scale the heatmap colours", value = TRUE),
        
        # ========== NEW: PLOT CUSTOMIZATION PANEL ==========
        box(
          title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
          width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
          h5("Text Sizes"),
          sliderInput("gc_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
          sliderInput("gc_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
          sliderInput("gc_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
          sliderInput("gc_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
          hr(),
          h5("Lines & Visual Elements"),
          sliderInput("gc_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
          sliderInput("gc_plot_point_size", "Point/marker size", min = 1, max = 5, value = 2, step = 0.5),
          sliderInput("gc_plot_alpha", "Transparency", min = 0.1, max = 1, value = 0.2, step = 0.1),
          hr(),
          h5("Other Elements"),
          checkboxInput("gc_plot_show_grid", "Show gridlines", value = TRUE),
          tags$small(style = "color:#888; display:block; margin-top:10px;",
            "Adjust and click 'Compute' to apply."
          )
          ),
            hr(),
            actionButton("run_genecontrib", "Compute Contributions",
                         class = "btn-primary btn-block")
          ),
          column(9,
            box(title = "Per-gene signal heatmap", width = 12, status = "primary",
              tags$small(style = "color:#888",
                "One row per gene, sorted by total signal in the window. The metagene line on the previous tab is essentially the column-wise average of this matrix — this is what that average is hiding."),
              plotOutput("gc_heatmap", height = "600px"),
              downloadButton("dl_gc_heatmap", "Download Heatmap (PDF)")
            ),
            box(title = "How concentrated is the signal?", width = 12, status = "primary",
              tags$small(style = "color:#888",
                "Genes ranked from strongest to weakest. A curve that shoots up immediately means a few genes carry the profile; a diagonal means every gene contributes equally."),
              plotOutput("gc_cumulative", height = "380px"),
              uiOutput("gc_concentration_text"),
              downloadButton("dl_gc_cumulative", "Download Curve (PDF)")
            ),
            box(title = "Gene ranking", width = 12, status = "primary",
              tags$small(style = "color:#888",
                "Every gene with its mean and peak signal in the window, its share of the total, and its running cumulative share."),
              DTOutput("gc_table"),
              br(),
              downloadButton("dl_gc_xlsx", "Download Gene Table (Excel)")
            )
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
              "Shows log(observed/expected) enrichment around gene boundaries (TSS and TES). Panels can be one per mark (conditions overlaid) or one per condition with all marks overlaid, as in the Galaxy chromstaR output."),
            hr(),
            numericInput("enr_upstream",   "Upstream of TSS (bp)",   value = 2000, min = 100, max = 10000, step = 100),
            numericInput("enr_downstream", "Downstream of TES (bp)", value = 2000, min = 100, max = 10000, step = 100),
            numericInput("enr_n_bins",     "Bins per region",        value = 40,   min = 10,  max = 200,  step = 5),
            hr(),
            h5("Gene scope"),
            radioButtons("gene_scope_enr", NULL,
                         choices = c("All genes in GTF"    = "all",
                                     "Only selected genes" = "selected"),
                         selected = "all"),
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
            h5("Stages to compare"),
            uiOutput("stage_selector_enr"),
            hr(),
            h5("Panel layout"),
            radioButtons("enr_layout", NULL,
                         choices = c("One panel per mark (colour = condition)"    = "by_mark",
                                     "One panel per condition (colour = mark)"    = "by_condition"),
                         selected = "by_mark"),
            hr(),
            h5("Ratio calculation"),
            radioButtons("enr_expected", NULL,
                         choices = c("Ratio of means (recommended)"        = "means",
                                     "Mean of per-bin log ratios (legacy)" = "legacy"),
                         selected = "means"),
            tags$small(style = "color:#888",
              "\"Ratio of means\" averages the raw signal at each position — zero bins included — then takes one log ratio against the genome-wide mean, as the Galaxy chromstaR profile does. \"Legacy\" logs each bin first and averages the logs, which forces zero bins to be discarded; that inflates exactly the positions where a mark is absent, so genuine depletion is flattened out."),
            hr(),
            checkboxInput("smooth_enrichment", "Smooth curve (LOESS)", value = FALSE),
        
        # ========== NEW: PLOT CUSTOMIZATION PANEL (ENRICHMENT) ==========
        box(
          title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
          width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
          h5("Text Sizes"),
          sliderInput("enr_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
          sliderInput("enr_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
          sliderInput("enr_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
          sliderInput("enr_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
          hr(),
          h5("Lines & Visual Elements"),
          sliderInput("enr_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
          sliderInput("enr_plot_point_size", "Point size", min = 1, max = 5, value = 2, step = 0.5),
          sliderInput("enr_plot_alpha", "Transparency", min = 0.1, max = 1, value = 0.2, step = 0.1),
          hr(),
          h5("Other Elements"),
          checkboxInput("enr_plot_show_grid", "Show gridlines", value = TRUE),
          tags$small(style = "color:#888; display:block; margin-top:10px;",
            "Adjust and click 'Compute Enrichment' to apply."
          )
          ),
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
            h5("Stages to compare"),
            uiOutput("stage_selector_browser"),
            hr(),
            h5("Display mode"),
            radioButtons("compare_mode_browser", "Condition comparison",
                         choices = c("Side-by-side" = "facet",
                                     "Overlay"       = "overlay"),
                         selected = "facet"),
        
        # ========== NEW: PLOT CUSTOMIZATION PANEL ==========
        box(
          title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
          width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
          h5("Text Sizes"),
          sliderInput("br_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
          sliderInput("br_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
          sliderInput("br_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
          sliderInput("br_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
          hr(),
          h5("Lines & Visual Elements"),
          sliderInput("br_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
          sliderInput("br_plot_point_size", "Point/marker size", min = 1, max = 5, value = 2, step = 0.5),
          sliderInput("br_plot_alpha", "Transparency", min = 0.1, max = 1, value = 0.2, step = 0.1),
          hr(),
          h5("Other Elements"),
          checkboxInput("br_plot_show_grid", "Show gridlines", value = TRUE),
          tags$small(style = "color:#888; display:block; margin-top:10px;",
            "Adjust and click 'Compute' to apply."
          )
          ),
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
            uiOutput("stage_selector_table"),
            tags$small(style = "color:#888",
              "When both stages are selected, rows from each are stacked with a \"stage\" column added — this table is not joined bin-for-bin across stages, since the two ChromstaR objects may not share identical binning."),
            br(), br(),
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
              "Filters merged chromatin segments by differential score and width, then — for each pair of conditions you choose — counts how many of those segments carry each mark in one condition but not the other. Works with any number of conditions in the object."),
            hr(),
            numericInput("diff_score_thresh", "Min differential score",
                        value = 0.9999, min = 0, max = 1, step = 0.0001),
            numericInput("diff_width_thresh", "Min merged region width (bp)",
                        value = 300, min = 0, step = 100),
            hr(),
            h5("Gene scope"),
            radioButtons("diff_gene_scope", NULL,
                         choices = c("Whole genome (all segments)" = "all",
                                     "Only selected genes"         = "selected"),
                         selected = "all"),
            conditionalPanel(
              condition = "input.diff_gene_scope == 'selected'",
              fluidRow(
                column(6, actionButton("select_all_genes_diff", "Select all",
                                       class = "btn-sm btn-block")),
                column(6, actionButton("clear_all_genes_diff", "Clear all",
                                       class = "btn-sm btn-block"))
              ),
              br(),
              uiOutput("gene_selector_diff")
            ),
            hr(),
            h5("Conditions to compare"),
            uiOutput("diff_condition_ui"),
            tags$small(style = "color:#888",
              "Pick 2–5 conditions. Drag to reorder; the × removes one."),
            br(), br(),
            radioButtons("diff_mode", "Comparison mode",
              choices = c("All pairwise combinations" = "pairwise",
                          "One reference vs the others" = "reference"),
              selected = "pairwise"),
            conditionalPanel(
              condition = "input.diff_mode == 'reference'",
              uiOutput("diff_ref_ui")
            ),
            hr(),
            h5("Stages to compare"),
            uiOutput("stage_selector_diffpeaks"),
            tags$small(style = "color:#888",
              "Each stage is filtered and counted independently and shown in its own row of panels."),
            hr(),
            h5("Ranked figure style"),
            radioButtons("dp_rank_style", NULL,
                         choices = c("Stacked (gained + lost)"        = "stacked",
                                     "Diverging (lost | gained)"      = "diverging"),
                         selected = "stacked"),
            tags$small(style = "color:#888",
              "Both styles draw gains and losses separately, so no single bar stands for a mixture of the two. Stacked keeps the ranking by total change; diverging puts losses left of zero and gains right."),
            hr(),
            # ========== NEW: PLOT CUSTOMIZATION PANEL (DIFFERENTIAL PEAKS) ==========
            box(
              title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
              width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
              h5("Text Sizes"),
              sliderInput("dp_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
              sliderInput("dp_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
              sliderInput("dp_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
              sliderInput("dp_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
              hr(),
              h5("Lines & Visual Elements"),
              sliderInput("dp_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
              sliderInput("dp_plot_point_size", "Point/marker size", min = 1, max = 5, value = 2, step = 0.5),
              sliderInput("dp_plot_alpha", "Transparency", min = 0.1, max = 1, value = 0.2, step = 0.1),
              hr(),
              h5("Other Elements"),
              checkboxInput("dp_plot_show_grid", "Show gridlines", value = TRUE),
              tags$small(style = "color:#888; display:block; margin-top:10px;",
                "Adjust and click 'Compute Differential Peaks' to apply."
              )
            ),
            actionButton("run_diffpeaks", "Compute Differential Peaks",
                        class = "btn-primary btn-block")
          ),
          box(title = "Differential Peaks per Histone Mark", width = 9, status = "primary",
            uiOutput("diffpeaks_plot_container"),
            downloadButton("dl_diffpeaks_plot", "Download Plot"),
            downloadButton("dl_diffpeaks_xlsx", "Download Data (Excel)")
          )
        ),
        fluidRow(
          box(title = "Per-mark summary — differential regions", width = 12,
              status = "primary", collapsible = TRUE,
            tags$small(style = "color:#888",
              "Gained / lost counts per mark across the differential regions, with each direction as its own percentage. The \"Direction\" label only says which side is larger — read the Gained and Lost columns for the actual numbers."),
            br(), br(),
            uiOutput("diffpeaks_overlap_note"),
            br(),
            DTOutput("diffpeaks_mark_table"),
            br(),
            downloadButton("dl_diffpeaks_mark_table", "Download Table (Excel)")
          )
        ),
        fluidRow(
          box(title = "Marks ranked by change — differential regions", width = 12,
              status = "primary", collapsible = TRUE,
            tags$small(style = "color:#888",
              "The same differential regions as the chart above, ranked by how many regions each mark changes in. Each bar is split into gained and lost, so the colour never stands for a total that mixes the two. The label is the combined percentage of differential regions — read the table above for the split."),
            plotOutput("diffpeaks_rank_filtered", height = "430px"),
            downloadButton("dl_diffpeaks_rank_filtered", "Download Plot (PDF)")
          )
        ),
        fluidRow(
          box(title = "Marks ranked by change — genome-wide (all domains)", width = 12,
              status = "primary", collapsible = TRUE,
            tags$small(style = "color:#888",
              "Unfiltered, domain-level view built from the object's combination-frequency table: every domain in the genome, with no score or width filter. Bars are split into gained and lost as above. Percentages are of all domains, so they are not comparable with the filtered figure."),
            uiOutput("diffpeaks_rank_genome_note"),
            plotOutput("diffpeaks_rank_genome", height = "430px"),
            downloadButton("dl_diffpeaks_rank_genome", "Download Plot (PDF)")
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
            h5("Stages to compare"),
            uiOutput("stage_selector_gsc"),
            tags$small(style = "color:#888",
              "Each stage is run independently through the same comparison (and, in permutation mode, its own random background), then combined into one results table/plot tagged by stage."),
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
        
        # ========== NEW: PLOT CUSTOMIZATION PANEL ==========
        box(
          title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
          width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
          h5("Text Sizes"),
          sliderInput("gsc_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
          sliderInput("gsc_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
          sliderInput("gsc_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
          sliderInput("gsc_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
          hr(),
          h5("Lines & Visual Elements"),
          sliderInput("gsc_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
          sliderInput("gsc_plot_point_size", "Point/marker size", min = 1, max = 5, value = 2, step = 0.5),
          sliderInput("gsc_plot_alpha", "Transparency", min = 0.1, max = 1, value = 0.2, step = 0.1),
          hr(),
          h5("Other Elements"),
          checkboxInput("gsc_plot_show_grid", "Show gridlines", value = TRUE),
          tags$small(style = "color:#888; display:block; margin-top:10px;",
            "Adjust and click 'Compute' to apply."
          )
          ),
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
              tags$p(tags$b("Purpose:"), " upload your data. Nothing else works until Stage A is loaded."),
              tags$ul(
                tags$li(tags$b("\"Load ChromstaR Object — Stage A\""), " — click to browse for your ChromstaR ", tags$code("combinedMultiHMM"), " object. Accepts the original ", tags$code(".RData"), "/", tags$code(".rda"), " format or a pre-slimmed ", tags$code(".rds"), " file (recommended for faster loading). Give it a short label (e.g. \"Protoscolex\", \"Metacestode\") — this label appears throughout the app's plots and tables. A text box below shows a summary once loaded."),
                tags$li(tags$b("\"Load ChromstaR Object — Stage B\""), " (optional) — load a second, independent ChromstaR object here to compare across developmental stages/conditions over time (e.g. protoscolex vs. metacestode, or dormant vs. activated). It does not need to share bin coordinates with Stage A — every comparison tab runs its own computation on each stage separately and combines the results, tagged by stage. Once loaded, every analysis tab gains a \"Stages to compare\" control."),
                tags$li(tags$b("\"Upload genes.tsv file\""), " — click to browse for your gene annotation, shared by both stages. This should be a simple table with columns ", tags$code("chr, start, end, gene_id, gene_name, strand"), " — convert a GTF to this format first if needed."),
                tags$li(tags$b("Detected Marks & Conditions"), " panel fills in automatically once each object is loaded — it lists every histone mark, every condition, and basic genome info found per stage. Use this to sanity-check the upload before moving to other tabs."),
                tags$li(tags$b("\"Condition / Life-cycle Order\""), " box at the bottom — by default conditions are ordered alphabetically, which rarely matches the biology. Drag the condition boxes into the order you want (e.g. ", tags$code("TwoDO → Meta → PNA → PA"), " to follow the life cycle) and that order is applied everywhere in the app at once: legend order, facet strip order, the order of ", tags$code("rpkm_*"), " columns in the exported table, and which colour each condition gets. A coloured preview strip below the box shows the current order and its colours. \"Reset to detected order\" puts it back to how the objects were read in.")
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
              tags$h3("2.5 Gene Contributions tab"),
              tags$p(tags$b("Purpose:"), " breaks a metagene profile back down into the individual genes behind it. An average curve cannot tell you whether every gene looks like that or whether ten genes carry the whole signal — this tab answers that."),
              tags$ul(
                tags$li(tags$b("\"Histone mark\" / \"Condition\" / \"Reference\""), " — one mark and one condition at a time, anchored at the TSS, TES or across the gene body, with the same windows and the same strand handling as the Metagene tab."),
                tags$li(tags$b("Per-gene signal heatmap"), " — one row per gene, columns are positions, sorted with the strongest gene at the top. Colours can be log-scaled and are capped at the 99th percentile so a single extreme gene does not wash out the rest."),
                tags$li(tags$b("Concentration curve"), " — genes ranked strongest first against their cumulative share of the total signal. A curve hugging the top-left means a few genes carry the profile; the dashed diagonal is what perfectly even contribution would look like."),
                tags$li(tags$b("Gini coefficient"), " — 0 means every gene contributes equally, 1 means a single gene carries everything. Above roughly 0.6, the mean profile describes a minority of genes and should be reported alongside a heatmap or a median line."),
                tags$li(tags$b("Gene ranking table"), " — every gene with its mean and peak signal, the position of its peak, its share of the total and the running cumulative share. Sortable, filterable and exportable."),
                tags$li(tags$b("\"Download Gene Table (Excel)\""), " — two sheets: the ranking, and the full gene-by-position matrix behind the heatmap.")
              ),

              hr(),
              tags$h3("3. Enrichment Profile tab"),
              tags$p(tags$b("Purpose:"), " plot log(observed/expected) enrichment of each mark around gene boundaries (TSS and TES together) — useful for seeing which marks are relatively enriched or depleted across a gene versus the genome-wide average."),
              tags$ul(
                tags$li(tags$b("\"Panel layout\""), " — \"One panel per mark\" overlays the conditions inside each mark's panel (best for asking how one mark changes across the life cycle); \"One panel per condition\" overlays all marks inside each condition's panel, the layout used by the Galaxy chromstaR output (best for asking which marks dominate at a given stage)."),
                tags$li(tags$b("\"Expected (baseline)\""), " — the denominator of log(observed/expected). The observed side has to drop bins with zero signal, because log(0) is undefined, so \"Mean over covered bins\" uses the same set of bins on both sides. \"Mean over all genomic bins\" divides by a mean that includes every empty bin, which shifts each curve upward by log(mean covered / mean all); because coverage sparsity differs per mark and per condition, that shift differs per curve and conditions can no longer be compared by their vertical position.")
              ),
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
              tags$p(tags$b("Purpose:"), " counts, per histone mark, how many chromatin segments are confidently present in one condition but not another (and vice versa) — a bar chart summarising which marks change the most between conditions. Any number of conditions is supported: you choose which ones to compare and how they are paired up."),
              tags$ul(
                tags$li(tags$b("\"Min differential score\""), " — only keep chromatin segments with a confidence score at or above this threshold (closer to 1 = stricter, fewer but more confident segments)."),
                tags$li(tags$b("\"Min merged region width (bp)\""), " — discard segments shorter than this, to avoid counting tiny noisy regions."),
                tags$li(tags$b("\"Gene scope\""), " — \"Whole genome\" counts every chromatin segment (the default, and what the Galaxy differential tool does); \"Only selected genes\" keeps only segments overlapping the genes you paste in, so you can ask whether a specific gene set is remodelled between stages."),
                tags$li(tags$b("\"Conditions to compare\""), " — pick between 2 and 5 conditions from the loaded object. Only these are used on this tab."),
                tags$li(tags$b("\"Comparison mode\""), " — \"All pairwise combinations\" makes one panel for every possible pair of the chosen conditions (5 conditions = 10 panels); \"One reference vs the others\" compares every chosen condition against a single reference you pick (5 conditions = 4 panels), which is usually what you want for a life-cycle baseline."),
                tags$li("Within each panel, the two bars per mark are coloured by the condition the mark is present in — same colours as everywhere else in the app."),
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
                tags$li(tags$b("Condition / life-cycle order"), " — set once on the Load Data tab, applied globally. Conditions are drawn, coloured, faceted and exported in the order you drag them into, so a plot spanning both stages reads left-to-right along the life cycle rather than alphabetically. Any condition you don't explicitly place is appended at the end rather than dropped."),
                tags$li(tags$b("Stages to compare"), " — appears on every analysis tab once a Stage B object is loaded. Each selected stage is run through the exact same computation independently (its own marks, its own conditions, its own bins), then the results are combined and labelled by stage — plots gain an extra facet row/column, and tables gain a \"stage\" column. Nothing requires the two ChromstaR objects to share identical binning."),
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
    bin_annotation = NULL,   # cached gene context per bin (Stage A)

    hmm_b            = NULL,
    marks_b          = character(0),
    conditions_b     = character(0),
    chroms_b         = character(0),
    bin_annotation_b = NULL  # cached gene context per bin (Stage B)
  )

  # ---- Load ChromstaR object — Stage A ---------------------------------------
  observeEvent(input$chromstar_file, {
    req(input$chromstar_file)
    withProgress(message = "Loading ChromstaR object (Stage A)…", {
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
        showNotification(paste("Error loading ChromstaR object (Stage A):", e$message), type = "error")
      })
    })
  })

  # ---- Load ChromstaR object — Stage B (optional) ----------------------------
  observeEvent(input$chromstar_file_b, {
    req(input$chromstar_file_b)
    withProgress(message = "Loading ChromstaR object (Stage B)…", {
      tryCatch({
        hmm_b <- load_chromstar_object(input$chromstar_file_b$datapath)
        rv$hmm_b <- hmm_b

        parsed              <- parse_marks_conditions(hmm_b)
        rv$marks_b          <- parsed$marks
        rv$conditions_b     <- parsed$conditions
        rv$chroms_b         <- as.character(unique(seqnames(hmm_b$bins)))
        rv$bin_annotation_b <- NULL  # reset cache on new object

        output$chromstar_summary_b <- renderPrint({
          bins_df   <- expand_bins_df(hmm_b)
          rpkm_cols <- grep("counts.rpkm", colnames(bins_df), value = TRUE)
          cat("Object class:", class(hmm_b), "\n")
          cat("Total bins:  ", length(hmm_b$bins), "\n")
          cat("Chromosomes: ", paste(head(rv$chroms_b, 5), collapse = ", "), "…\n")
          cat("Conditions:  ", paste(rv$conditions_b, collapse = ", "), "\n")
          cat("Marks:       ", paste(rv$marks_b, collapse = ", "), "\n\n")
          if (length(rpkm_cols) > 0) {
            cat("counts.rpkm columns (", length(rpkm_cols), "):\n")
            cat(" ", paste(head(rpkm_cols, 20), collapse = "\n  "), "\n")
          } else {
            cat("WARNING: No counts.rpkm.* columns found after expansion!\n")
            cat("All mcols names: ", paste(head(colnames(bins_df), 30), collapse = ", "), "\n")
          }
          if (!is.null(hmm_b$info)) {
            cat("\nhmm_b$info:\n")
            print(hmm_b$info[, c("mark","condition","replicate"), drop = FALSE])
          }
        })
      }, error = function(e) {
        showNotification(paste("Error loading ChromstaR object (Stage B):", e$message), type = "error")
      })
    })
  })

  # ---- Stage labels (deduplicated so two stages never collide) --------------
  `%||%` <- function(a, b) if (is.null(a)) b else a
  stage_a_label <- reactive({
    lbl <- trimws(input$stage_a_label %||% "")
    if (!nzchar(lbl)) lbl <- "Stage A"
    lbl
  })
  stage_b_label <- reactive({
    lbl <- trimws(input$stage_b_label %||% "")
    if (!nzchar(lbl)) lbl <- "Stage B"
    if (!is.null(rv$hmm_b) && identical(lbl, stage_a_label())) lbl <- paste0(lbl, " (B)")
    lbl
  })

  # ---- Unified view of whichever stage(s) are currently loaded ---------------
  # Returns a named list keyed by the stage's display label, so downstream code
  # (mark/stage selectors, per-stage compute loops, plot facets) can treat
  # Stage A and Stage B uniformly and stay correct whether one or both are loaded.
  active_stages <- reactive({
    stages <- list()
    if (!is.null(rv$hmm)) {
      stages[[stage_a_label()]] <- list(
        hmm = rv$hmm, marks = rv$marks, conditions = rv$conditions,
        chroms = rv$chroms, key = "A"
      )
    }
    if (!is.null(rv$hmm_b)) {
      stages[[stage_b_label()]] <- list(
        hmm = rv$hmm_b, marks = rv$marks_b, conditions = rv$conditions_b,
        chroms = rv$chroms_b, key = "B"
      )
    }
    stages
  })

  # Reusable "Stages to compare" checkbox UI — only rendered with a visible
  # choice once Stage B is actually loaded; with only Stage A loaded it's
  # skipped entirely and every compute path silently just uses Stage A.
  stage_selector_ui <- function(input_id) {
    stages <- names(active_stages())
    if (length(stages) < 2) return(NULL)
    checkboxGroupInput(input_id, NULL, choices = stages, selected = stages, inline = TRUE)
  }

  # Given a (possibly NULL/empty) selectize/checkbox input value, resolve it to
  # the actual stage names to use — falls back to "all currently loaded stages"
  # so tabs work fine before the user has touched the stage selector.
  resolve_selected_stages <- function(input_val) {
    stages <- names(active_stages())
    if (is.null(input_val) || length(input_val) == 0) return(stages)
    sel <- intersect(input_val, stages)
    if (length(sel) == 0) stages else sel
  }

  # Run `compute_fn(stage_info)` once per selected stage, tag each non-NULL
  # data.frame result with a "stage" column, and row-bind the results. This is
  # the core mechanism that turns every single-object compute_* helper into a
  # multi-stage one without having to touch those helpers themselves.
  compute_across_stages <- function(selected_stage_names, compute_fn) {
    stages  <- active_stages()
    results <- list()
    for (nm in selected_stage_names) {
      st <- stages[[nm]]
      if (is.null(st)) next
      res <- compute_fn(st)
      if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
        res$stage <- nm
        results[[nm]] <- res
      }
    }
    if (length(results) == 0) return(NULL)
    out <- bind_rows(results)
    # Stages keep load order (A then B) rather than alphabetical
    out$stage <- factor(out$stage, levels = names(stages))
    apply_condition_order(out)
  }

  # ---- Condition (life-cycle) ordering --------------------------------------
  # Every condition found across both loaded objects, in raw detected order.
  detected_conditions <- reactive({
    unique(unlist(lapply(active_stages(), function(s) s$conditions)))
  })

  # The user's chosen conditions, in their chosen order. This drives BOTH the
  # display order and which conditions are used at all: a condition removed
  # from the box is dropped from every plot, table and analysis in the app.
  # Clearing the box entirely falls back to "all detected conditions" so the
  # app is never left with nothing to plot.
  ordered_conditions <- reactive({
    all_conds <- detected_conditions()
    chosen    <- intersect(input$condition_order %||% character(0), all_conds)
    if (length(chosen) == 0) return(all_conds)
    chosen
  })

  # The conditions of one stage that survive the user's Condition/Life-cycle
  # Order box, in that order. Falls back to the stage's own conditions if the
  # kept set doesn't overlap this stage at all (e.g. Stage B has different
  # condition names), so a stage is never silently reduced to nothing.
  stage_conditions_kept <- function(st) {
    kept <- intersect(ordered_conditions(), st$conditions)
    if (length(kept) == 0) st$conditions else kept
  }

  # Turn a data.frame's `condition` column into a factor with the user's chosen
  # level order, so ggplot legends, facet strips and DT tables all follow the
  # life cycle instead of sorting alphabetically.
  apply_condition_order <- function(df) {
    if (is.null(df) || !is.data.frame(df) || !"condition" %in% colnames(df)) return(df)
    keep     <- ordered_conditions()
    cond_chr <- as.character(df$condition)
    # Drop conditions the user removed from the order box. Guarded so that a
    # data.frame whose `condition` column holds something unrelated (nothing
    # matches) is passed through untouched rather than emptied.
    if (any(cond_chr %in% keep)) {
      df       <- df[cond_chr %in% keep, , drop = FALSE]
      cond_chr <- as.character(df$condition)
    }
    lv <- c(intersect(keep, unique(cond_chr)),
            setdiff(unique(cond_chr), keep))
    df$condition <- factor(cond_chr, levels = lv)
    df
  }

  output$condition_order_ui <- renderUI({
    conds <- detected_conditions()
    if (length(conds) == 0) {
      return(tags$em(style = "color:#888", "Load a ChromstaR object to see its conditions here."))
    }
    # `selected` is isolated: without it this renderUI would re-run on every
    # edit of the box (via ordered_conditions()) and rebuild the widget with
    # every condition selected again, so removals appeared to bounce back.
    selectizeInput(
      "condition_order", NULL,
      choices  = conds,
      selected = isolate(ordered_conditions()),
      multiple = TRUE,
      width    = "100%",
      options  = list(plugins = list("drag_drop", "remove_button"))
    )
  })

  observeEvent(input$reset_condition_order, {
    updateSelectizeInput(session, "condition_order",
                        choices = detected_conditions(),
                        selected = detected_conditions())
  })

  output$condition_order_preview <- renderUI({
    conds <- ordered_conditions()
    if (length(conds) == 0) return(NULL)
    pal <- condition_palette()
    chips <- lapply(seq_along(conds), function(i) {
      tags$span(
        style = paste0(
          "display:inline-block; margin:2px 4px; padding:3px 10px; border-radius:12px;",
          "background:", pal[[conds[i]]], "; color:white; font-weight:bold; font-size:12px;"
        ),
        paste0(i, ". ", conds[i])
      )
    })
    tagList(
      tags$small(style = "color:#555", "Current order (and the colour each condition will get):"),
      br(),
      tags$div(chips)
    )
  })

  # A colour palette keyed by *condition name*, following the user's chosen
  # life-cycle order, so a given condition (e.g. "PA") always gets the same
  # colour whether it's plotted from Stage A alone or alongside Stage B.
  base_palette <- c("#E63946", "#2196F3", "#FF9800", "#4CAF50",
                    "#9C27B0", "#00BCD4", "#795548", "#607D8B")
  # Default colour for a condition, purely from its position in the order.
  default_condition_palette <- function(conds) {
    stats::setNames(base_palette[((seq_along(conds) - 1) %% length(base_palette)) + 1], conds)
  }

  # The palette every figure in the app uses. Starts from the position-based
  # default and lets the per-condition colour pickers override it, so one
  # condition keeps the same colour across every tab.
  condition_palette <- function() {
    conds <- ordered_conditions()
    pal   <- default_condition_palette(conds)
    for (cnd in conds) {
      v <- input[[cond_colour_id(cnd)]]
      if (!is.null(v) && length(v) == 1 && grepl("^#[0-9A-Fa-f]{6}$", v)) pal[[cnd]] <- v
    }
    pal
  }

  # One picker per condition. isolate() on the current value means re-rendering
  # (e.g. after reordering conditions) keeps whatever the user already chose.
  output$condition_colour_ui <- renderUI({
    conds <- ordered_conditions()
    req(length(conds) > 0)
    defaults <- default_condition_palette(conds)
    tagList(lapply(conds, function(cnd) {
      id  <- cond_colour_id(cnd)
      cur <- isolate(input[[id]])
      val <- if (!is.null(cur) && grepl("^#[0-9A-Fa-f]{6}$", cur)) cur else unname(defaults[[cnd]])
      div(style = "display:inline-block; margin-right:18px; vertical-align:top;",
          condition_colour_input(id, cnd, val))
    }))
  })

  observeEvent(input$reset_condition_colours, {
    conds    <- ordered_conditions()
    defaults <- default_condition_palette(conds)
    for (cnd in conds) {
      updateTextInput(session, cond_colour_id(cnd), value = unname(defaults[[cnd]]))
    }
  })
  # A colour palette keyed by *mark name*, used when marks are the thing being
  # compared inside a panel (Galaxy-style "one panel per condition" layout).
  # Okabe-Ito qualitative palette — distinguishable in colour-blind vision and
  # in greyscale print.
  mark_base_palette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                         "#D55E00", "#56B4E9", "#F0E442", "#000000")
  mark_palette <- function(marks) {
    marks <- unique(marks)
    setNames(mark_base_palette[((seq_along(marks) - 1) %% length(mark_base_palette)) + 1],
             marks)
  }

  # Union of marks across every loaded stage — used to build "Marks to
  # display" checkboxes so a mark present in only one stage still shows up.
  all_marks_reactive <- reactive({
    unique(unlist(lapply(active_stages(), function(s) s$marks)))
  })
  all_chroms_reactive <- reactive({
    unique(unlist(lapply(active_stages(), function(s) s$chroms)))
  })

  output$stage_a_label_header <- renderText(stage_a_label())
  output$stage_b_label_header <- renderText(if (!is.null(rv$hmm_b)) stage_b_label() else "Stage B (not loaded)")

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
        # Drop the cached bin annotation so it is recomputed against THIS gene
        # file — otherwise a stale annotation from a previously loaded (or
        # mismatched) annotation stays in the Data Table for the whole session.
        rv$bin_annotation   <- NULL
        rv$bin_annotation_b <- NULL

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

  output$marks_ui_b <- renderUI({
    if (is.null(rv$hmm_b)) return(tags$p(tags$em(style = "color:#888", "Not loaded.")))
    req(length(rv$marks_b) > 0)
    tags$ul(lapply(rv$marks_b, tags$li))
  })

  output$conditions_ui_b <- renderUI({
    if (is.null(rv$hmm_b)) return(tags$p(tags$em(style = "color:#888", "Not loaded.")))
    req(length(rv$conditions_b) > 0)
    tags$ul(lapply(rv$conditions_b, tags$li))
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

  output$genome_info_b <- renderPrint({
    if (is.null(rv$hmm_b)) {
      cat("Not loaded.\n")
      return(invisible(NULL))
    }
    si <- seqinfo(rv$hmm_b$bins)
    print(si)
  })

  # Mark selectors (each needs its own renderUI call — reusing one renderUI
  # result across two outputs does not work correctly in Shiny)
  output$mark_selector_meta <- renderUI({
    marks <- all_marks_reactive()
    req(length(marks) > 0)
    checkboxGroupInput("marks_meta", NULL,
                       choices  = marks,
                       selected = marks)
  })

  output$mark_selector_browser <- renderUI({
    marks <- all_marks_reactive()
    req(length(marks) > 0)
    checkboxGroupInput("marks_browser", NULL,
                       choices  = marks,
                       selected = marks)
  })

  output$stage_selector_meta     <- renderUI(stage_selector_ui("stages_meta"))
  output$stage_selector_enr      <- renderUI(stage_selector_ui("stages_enr"))
  output$stage_selector_browser  <- renderUI(stage_selector_ui("stages_browser"))
  output$stage_selector_diffpeaks <- renderUI(stage_selector_ui("stages_diffpeaks"))
  output$stage_selector_gsc      <- renderUI(stage_selector_ui("stages_gsc"))
  output$stage_selector_table    <- renderUI(stage_selector_ui("stages_table"))

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
    defaults <- ""   # deliberately empty — a pre-filled list looks like a
                     # placeholder but is real input, and silently reduced
                     # every profile to those few genes
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
    marks <- all_marks_reactive()
    req(length(marks) > 0)
    checkboxGroupInput("marks_enr", NULL,
                       choices  = marks,
                       selected = marks)
  })

  output$gene_selector_enr <- renderUI({
    req(rv$genes)
    nms      <- gene_names_reactive()
    defaults <- ""   # deliberately empty — a pre-filled list looks like a
                     # placeholder but is real input, and silently reduced
                     # every profile to those few genes
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

    stage_names <- resolve_selected_stages(input$stages_enr)

    withProgress(message = paste0("Computing enrichment (", length(genes_sel), " genes)…"), {
      compute_across_stages(stage_names, function(st) {
        compute_enrichment_profile(
          hmm        = st$hmm,
          genes_gr   = genes_sel,
          marks      = input$marks_enr,
          conditions = stage_conditions_kept(st),
          upstream      = input$enr_upstream,
          downstream    = input$enr_downstream,
          n_bins        = input$enr_n_bins,
          ratio_mode    = input$enr_expected %||% "means"
        )
      })
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

    # Keep the user's life-cycle order for the legend and line styles
    conds       <- levels(droplevels(as.factor(df$condition)))
    multi_stage <- length(unique(df$stage)) > 1
    by_cond     <- identical(input$enr_layout, "by_condition")

    # Whichever variable is compared *inside* a panel gets the colour and the
    # line type; the other one becomes the facet.
    if (by_cond) {
      grp_levels <- unique(as.character(df$mark))
      df$grp     <- factor(as.character(df$mark), levels = grp_levels)
      line_cols  <- mark_palette(grp_levels)
      legend_lab <- "Mark"
    } else {
      grp_levels <- conds
      df$grp     <- factor(as.character(df$condition), levels = grp_levels)
      line_cols  <- condition_palette()
      legend_lab <- "Condition"
    }
    line_types <- setNames(rep(c("solid", "dotted", "dashed", "dotdash"),
                              length.out = length(grp_levels)), grp_levels)

    # Shade the gene-body zone (x in [0,1]) so it's visually distinct from
    # the upstream/downstream flanks even before looking at axis labels
    body_shade <- data.frame(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf)

    p <- ggplot(df, aes(x = x, y = mean_log_ratio,
                        colour = grp, linetype = grp)) +
      geom_rect(data = body_shade,
               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
               inherit.aes = FALSE, fill = "grey85", alpha = 0.35)

    if (isTRUE(input$smooth_enrichment)) {
      p <- p + geom_smooth(method = "loess", span = 0.2, se = FALSE, linewidth = input$enr_plot_line_width %||% 0.9)
    } else {
      p <- p + geom_line(linewidth = input$enr_plot_line_width %||% 0.8)
    }

    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_vline(xintercept = 0, linetype = "solid",  colour = "grey30", linewidth = 0.6) +
      geom_vline(xintercept = 1, linetype = "solid",  colour = "grey30", linewidth = 0.6) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels,
                         expand = expansion(mult = 0.02)) +
      scale_colour_manual(values = line_cols) +
      scale_linetype_manual(values = line_types)

    if (by_cond) {
      # Galaxy-style: one panel per condition, every mark overlaid. A shared
      # y scale here, so panels can be compared to each other directly.
      if (multi_stage) {
        p <- p + facet_grid(stage ~ condition)
      } else {
        p <- p + facet_wrap(~ condition, ncol = 2)
      }
    } else {
      if (multi_stage) {
        p <- p + facet_grid(stage ~ mark, scales = "free_y")
      } else {
        p <- p + facet_wrap(~ mark, scales = "free_y", ncol = 4)
      }
    }

    p <- p +
      labs(
        x = NULL,
        y = "log(observed/expected)",
        colour = legend_lab, linetype = legend_lab,
        title    = if (by_cond) "Histone Mark Enrichment — Per Condition"
                   else         "Histone Mark Enrichment — Per Mark",
        subtitle = paste0(
          if (multi_stage) paste(unique(df$stage), collapse = " vs ")
          else             paste(conds, collapse = " vs "),
          if ("n_genes_distinct" %in% colnames(df)) {
            v   <- df$n_genes_distinct[is.finite(df$n_genes_distinct)]
            tot <- if ("n_genes_total" %in% colnames(df))
                     suppressWarnings(max(df$n_genes_total, na.rm = TRUE)) else NA_real_
            paste0(
              if (is.finite(tot))
                paste0("  |  ", format(tot, big.mark = ","), " genes total") else "",
              if (length(v))
                paste0("  |  ", format(min(v), big.mark = ","), "\u2013",
                       format(max(v), big.mark = ","), " genes per position") else "")
          } else "",
          "  |  grey band = gene body, vertical lines = TSS and TES",
          if (identical(input$enr_expected %||% "means", "legacy"))
            "  |  legacy per-bin log ratios" else "  |  ratio of means"
        )
      ) +
      build_custom_theme(
        title_size = input$enr_plot_title_size %||% 14,
        axis_title_size = input$enr_plot_axis_title_size %||% 12,
        axis_text_size = input$enr_plot_axis_text_size %||% 10,
        legend_text_size = input$enr_plot_legend_size %||% 10,
        line_width = input$enr_plot_line_width %||% 0.8
      ) +
      theme(legend.position = "top")
    
    if (!isTRUE(input$enr_plot_show_grid)) {
      p <- p + theme(panel.grid.major = element_blank())
    }
    p
  })

  output$enrichment_plot <- renderPlot({ with_stack_trace(build_enrichment_plot(), "Enrichment Profile") })

  output$dl_enrichment_plot <- downloadHandler(
    filename = function() paste0("enrichment_profile_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = build_enrichment_plot(),
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
      keep_cols <- intersect(
        c("stage", "mark", "condition", "slot", "x_position",
          "mean_log_ratio", "mean_value", "expected", "n_bins_avg", "n_genes",
          "n_genes_distinct", "n_genes_total", "genes"),
        colnames(df_out))
      df_out <- df_out[, keep_cols, drop = FALSE]
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
    chroms <- all_chroms_reactive()
    req(length(chroms) > 0)
    if (identical(input$chrom_scope, "all")) {
      tags$small(style = "color:#888",
        paste0("All ", length(chroms), " chromosomes/contigs (union across loaded stages) will be shown — this can be slow with hundreds of sequences."))
    } else if (identical(input$chrom_scope, "multi")) {
      selectizeInput("region_chrom", "Chromosomes (select one or more)",
                     choices  = chroms,
                     selected = chroms[1],
                     multiple = TRUE,
                     options  = list(maxOptions = length(chroms) + 10))
    } else {
      selectInput("region_chrom", "Chromosome",
                  choices  = chroms,
                  selected = chroms[1])
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

    stage_names <- resolve_selected_stages(input$stages_meta)
    modes_sel   <- input$metagene_mode
    n_steps     <- max(1, length(modes_sel) * length(stage_names))

    withProgress(message = paste0("Computing metagene (", length(genes_sel), " genes)…"), {
      compute_across_stages(stage_names, function(st) {
        results <- list()
        for (m in modes_sel) {
          incProgress(1 / n_steps, detail = paste0("Reference: ", m))
          res <- compute_metagene(
            hmm        = st$hmm,
            genes_gr   = genes_sel,
            marks      = input$marks_meta,
            conditions = stage_conditions_kept(st),
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
        res_all <- bind_rows(results)
        res_all$n_genes_used <- length(genes_sel)
        res_all
      })
    })
  })

  build_metagene_plot <- reactive({
    df <- metagene_data()
    req(df)

    n_bins      <- input$n_bins

    # Reference panels follow the biology (TSS -> gene body -> TES) rather than
    # the alphabetical order a bare character column would give.
    ref_order   <- c("TSS", "gene_body", "TES")
    present     <- unique(as.character(df$ref_type))
    modes_sel   <- c(intersect(ref_order, present), setdiff(present, ref_order))
    df$ref_type <- factor(as.character(df$ref_type), levels = modes_sel,
                          labels = ifelse(modes_sel == "gene_body",
                                          "Gene body", modes_sel))
    multi_ref   <- length(modes_sel) > 1
    multi_stage <- length(unique(df$stage)) > 1

    x_breaks <- c(1, round(n_bins / 2), n_bins)

    label_for_mode <- function(m) {
      if (m == "TSS") {
        c(paste0("-", input$upstream / 1000, "kb"), "TSS", paste0("+", input$downstream / 1000, "kb"))
      } else if (m == "TES") {
        c(paste0("-", input$upstream / 1000, "kb"), "TES", paste0("+", input$downstream / 1000, "kb"))
      } else {
        c("TSS", "50% of gene", "TES")
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

    # Build the column-facet dimensions dynamically: ref_type (when more than
    # one reference is selected), stage (when Stage B is included in the
    # comparison), and condition (only in "Side-by-side" mode — "Overlay"
    # keeps condition as colour instead of a facet column).
    col_terms <- character(0)
    if (multi_ref)                          col_terms <- c(col_terms, "ref_type")
    if (multi_stage)                        col_terms <- c(col_terms, "stage")
    if (input$compare_mode_meta == "facet") col_terms <- c(col_terms, "condition")

    if (multi_ref) {
      # When multiple reference modes are shown together, x-axis meaning differs
      # per column (TSS/TES = bp distance, gene_body = %), so we can't use one
      # global scale_x_continuous with fixed labels — instead facet by ref_type
      # and rely on the facet strip text ("TSS"/"TES"/"gene_body") for context,
      # using a single consistent breaks/labels set that's "close enough"
      # (Start/Mid/End) across all reference types for visual alignment.
      p <- p + scale_x_continuous(breaks = x_breaks, labels = c("Start", "Mid", "End"))
      title_txt <- paste("Metagene profile —", paste(modes_sel, collapse = " / "))
    } else {
      x_labels <- label_for_mode(modes_sel[1])
      p <- p + scale_x_continuous(breaks = x_breaks, labels = x_labels)
      title_txt <- paste("Metagene profile —", modes_sel[1])
    }

    if (length(col_terms) == 0) {
      p <- p + facet_wrap(~ mark, scales = "free_y", ncol = 2)
    } else {
      facet_formula <- as.formula(paste("mark ~", paste(col_terms, collapse = " + ")))
      p <- p + facet_grid(facet_formula, scales = "free_y")
    }

    if (multi_stage) {
      title_txt <- paste0(title_txt, "  |  ", paste(unique(df$stage), collapse = " vs "))
    }
    if ("n_genes_used" %in% colnames(df)) {
      title_txt <- paste0(title_txt, "  |  ",
                          format(max(df$n_genes_used, na.rm = TRUE), big.mark = ","),
                          " genes")
    }

    pal <- condition_palette()
    p <- p +
      scale_colour_manual(values = pal) +
      scale_fill_manual(values   = pal) +
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

  output$metagene_plot <- renderPlot({ apply_plot_customization(build_metagene_plot(), input, "mg_") })

  output$dl_metagene <- downloadHandler(
    filename = function() paste0("metagene_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = apply_plot_customization(build_metagene_plot(), input, "mg_"),
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


  # ---- GENE CONTRIBUTIONS ----------------------------------------------------
  # Breaks a metagene profile back down into its individual genes.

  output$gc_stage_ui <- renderUI({
    stages <- names(active_stages())
    if (length(stages) < 2) return(NULL)
    selectInput("gc_stage", "Stage", choices = stages, selected = stages[1])
  })

  output$gc_mark_ui <- renderUI({
    marks <- all_marks_reactive()
    req(length(marks) > 0)
    selectInput("gc_mark", "Histone mark", choices = marks, selected = marks[1])
  })

  output$gc_condition_ui <- renderUI({
    conds <- ordered_conditions()
    req(length(conds) > 0)
    selectInput("gc_condition", "Condition", choices = conds, selected = conds[1])
  })

  output$gene_selector_gc <- renderUI({
    req(rv$genes)
    nms <- gene_names_reactive()
    tagList(
      tags$small(style = "color:#888",
        paste0(length(nms), " genes in the annotation — paste a list, one per row.")),
      br(), br(),
      tags$textarea(
        id = "selected_genes_text_gc", rows = 10,
        style = "width:100%; height:180px; font-family: monospace; font-size: 13px; resize: vertical;",
        placeholder = "Paste gene IDs here, one per line…",
        ""
      )
    )
  })

  gc_genes_selected <- reactive({
    req(rv$genes)
    if (!identical(input$gc_gene_scope, "selected")) return(rv$genes)
    raw <- strsplit(input$selected_genes_text_gc %||% "", "[\r\n,;\t]+")[[1]]
    raw <- trimws(raw); raw <- raw[nchar(raw) > 0]
    nms <- gene_names_reactive()
    sel <- raw[raw %in% nms]
    if (length(sel) == 0) return(NULL)
    rv$genes[nms %in% sel]
  })

  gc_data <- eventReactive(input$run_genecontrib, {
    req(rv$hmm, rv$genes, input$gc_mark, input$gc_condition)

    genes_sel <- gc_genes_selected()
    if (is.null(genes_sel) || length(genes_sel) == 0) {
      showNotification("No valid genes selected.", type = "error", duration = 8)
      return(NULL)
    }

    stages <- active_stages()
    st_nm  <- input$gc_stage %||% names(stages)[1]
    st     <- stages[[st_nm]]
    req(st)

    if (!(input$gc_condition %in% st$conditions)) {
      showNotification(
        paste0("Condition \"", input$gc_condition, "\" is not present in ", st_nm, "."),
        type = "error", duration = 10
      )
      return(NULL)
    }

    withProgress(message = paste0("Splitting profile across ", length(genes_sel), " genes…"), {
      df <- compute_gene_contributions(
        hmm        = st$hmm,
        genes_gr   = genes_sel,
        mark       = input$gc_mark,
        condition  = input$gc_condition,
        mode       = input$gc_mode,
        upstream   = input$gc_upstream,
        downstream = input$gc_downstream,
        n_bins     = input$gc_n_bins
      )
      if (is.null(df) || nrow(df) == 0) {
        showNotification("No bins overlapped the selected genes' windows.",
                         type = "warning", duration = 10)
        return(NULL)
      }
      df$stage <- st_nm

      # Coordinates for every gene in the run, so the ranking table can say
      # WHERE each contributing gene is, not just its name.
      gi_df <- data.frame(
        gene   = mcols(genes_sel)$gene_name,
        chr    = as.character(seqnames(genes_sel)),
        start  = start(genes_sel),
        end    = end(genes_sel),
        strand = as.character(strand(genes_sel)),
        stringsAsFactors = FALSE
      )
      gi_df$TSS         <- ifelse(gi_df$strand == "-", gi_df$end,   gi_df$start)
      gi_df$TES         <- ifelse(gi_df$strand == "-", gi_df$start, gi_df$end)
      gi_df$gene_length <- gi_df$end - gi_df$start + 1
      gi_df$locus       <- paste0(gi_df$chr, ":", gi_df$start, "-", gi_df$end)

      list(profile = df, coords = gi_df)
    })
  })

  # One row per gene: mean and peak signal in the window, share of the total.
  gc_gene_summary <- reactive({
    res <- gc_data(); req(res)
    df  <- res$profile
    out <- df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        mean_signal = mean(signal, na.rm = TRUE),
        peak_signal = max(signal,  na.rm = TRUE),
        peak_bin    = bin_idx[which.max(signal)],
        n_bins      = dplyr::n(),
        total       = sum(signal, na.rm = TRUE),
        .groups     = "drop"
      ) %>%
      dplyr::arrange(desc(total))
    grand <- sum(out$total, na.rm = TRUE)
    out$pct_of_total     <- if (grand > 0) 100 * out$total / grand else NA_real_
    out$cum_pct_of_total <- cumsum(out$pct_of_total)
    out$rank             <- seq_len(nrow(out))

    # Attach where each gene actually is
    out <- dplyr::left_join(out, res$coords, by = "gene")
    out[, c("rank", "gene", "chr", "start", "end", "strand", "TSS", "TES",
            "gene_length", "locus", "mean_signal", "peak_signal", "peak_bin",
            "n_bins", "total", "pct_of_total", "cum_pct_of_total")]
  })

  build_gc_heatmap <- reactive({
    res <- gc_data(); req(res)
    df  <- res$profile
    sm  <- gc_gene_summary()

    top_n <- min(nrow(sm), max(10, input$gc_top_n %||% 200))
    keep  <- sm$gene[seq_len(top_n)]

    d <- df[df$gene %in% keep, , drop = FALSE]
    d$gene <- factor(d$gene, levels = rev(keep))   # strongest at the top

    fill_val <- if (isTRUE(input$gc_log_scale)) log1p(d$signal) else d$signal
    d$fill_val <- fill_val
    cap <- stats::quantile(d$fill_val, 0.99, na.rm = TRUE)
    d$fill_val <- pmin(d$fill_val, cap)

    x_lab <- if (identical(input$gc_mode, "gene_body")) "TSS → TES (% of gene)" else
             paste0("position relative to ", input$gc_mode)

    ggplot(d, aes(x = bin_idx, y = gene, fill = fill_val)) +
      geom_raster() +
      scale_fill_viridis_c(
        option = "magma", direction = -1,
        name = if (isTRUE(input$gc_log_scale)) "log(1+RPKM)" else "RPKM"
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(
        x = x_lab, y = NULL,
        title = paste0("Per-gene signal — ", input$gc_mark, " / ", input$gc_condition),
        subtitle = paste0("top ", top_n, " of ", nrow(sm),
                          " genes by total signal, strongest at the top")
      ) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.y      = if (top_n <= 60) element_text(size = 6) else element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "right"
      )
  })

  build_gc_cumulative <- reactive({
    sm <- gc_gene_summary(); req(sm)
    d  <- data.frame(
      pct_genes  = 100 * sm$rank / nrow(sm),
      pct_signal = sm$cum_pct_of_total
    )
    ggplot(d, aes(x = pct_genes, y = pct_signal)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(linewidth = 1, colour = "#E63946") +
      scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      labs(
        x = "% of genes (ranked strongest first)",
        y = "% of total signal",
        title = "Signal concentration across genes",
        subtitle = "dashed line = every gene contributes equally"
      ) +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })

  output$gc_heatmap    <- renderPlot({ apply_plot_customization(build_gc_heatmap(), input, "gc_") })
  output$gc_cumulative <- renderPlot({ apply_plot_customization(build_gc_cumulative(), input, "gc_") })

  output$gc_concentration_text <- renderUI({
    sm <- gc_gene_summary(); req(sm)
    n  <- nrow(sm)
    at <- function(pct) {
      k <- max(1L, ceiling(n * pct / 100))
      round(sm$cum_pct_of_total[k], 1)
    }
    # Gini coefficient of the per-gene totals: 0 = perfectly even, 1 = one gene
    x    <- sort(pmax(sm$total, 0))
    gini <- if (sum(x) > 0) {
      (2 * sum(seq_along(x) * x) / (length(x) * sum(x))) - (length(x) + 1) / length(x)
    } else NA_real_

    tagList(
      br(),
      tags$b("Top 1% of genes carry "), tags$b(paste0(at(1), "%")), " of the signal; ",
      tags$b("top 10% carry "), tags$b(paste0(at(10), "%")), "; ",
      tags$b("top 50% carry "), tags$b(paste0(at(50), "%")), ".",
      br(),
      tags$small(style = "color:#888",
        paste0("Gini coefficient = ", round(gini, 3),
               " (0 = every gene contributes equally, 1 = a single gene carries everything). ",
               "Above ~0.6 the average profile is telling you about a minority of genes, ",
               "so report the median gene or a heatmap alongside the mean line."))
    )
  })

  output$gc_table <- renderDT({
    sm <- gc_gene_summary(); req(sm)
    out <- sm[, c("rank", "gene", "locus", "strand", "gene_length",
                  "mean_signal", "peak_signal", "peak_bin",
                  "pct_of_total", "cum_pct_of_total")]
    out$mean_signal      <- round(out$mean_signal, 3)
    out$peak_signal      <- round(out$peak_signal, 3)
    out$pct_of_total     <- round(out$pct_of_total, 4)
    out$cum_pct_of_total <- round(out$cum_pct_of_total, 2)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_gc_heatmap <- downloadHandler(
    filename = function() paste0("gene_contributions_heatmap_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, build_gc_heatmap(), width = 10, height = 12, device = "pdf", limitsize = FALSE)
    }
  )

  output$dl_gc_cumulative <- downloadHandler(
    filename = function() paste0("gene_contributions_curve_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, build_gc_cumulative(), width = 8, height = 6, device = "pdf")
    }
  )

  output$dl_gc_xlsx <- downloadHandler(
    filename = function() paste0("gene_contributions_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      sm <- gc_gene_summary(); req(sm)
      df <- gc_data()$profile
      wide <- tidyr::pivot_wider(df[, c("gene", "bin_idx", "signal")],
                                 names_from = bin_idx, values_from = signal,
                                 names_prefix = "bin_")
      write.xlsx(list(`gene ranking` = as.data.frame(sm),
                      `per-gene matrix` = as.data.frame(wide)),
                 file, overwrite = TRUE)
    }
  )

  # ---- REGION BROWSER PLOT ---------------------------------------------------
  browser_data <- eventReactive(input$run_browser, {
    req(rv$hmm)

    marks_use <- input$marks_browser
    if (is.null(marks_use) || length(marks_use) == 0) {
      marks_use <- all_marks_reactive()  # fallback: use all marks if checkboxes haven't registered
      if (length(marks_use) == 0) {
        showNotification("No histone marks available — load a ChromstaR object first.",
                         type = "error", duration = 6)
        return(NULL)
      }
    }

    if (identical(input$chrom_scope, "all")) {
      chroms_use <- all_chroms_reactive()
    } else if (identical(input$chrom_scope, "multi")) {
      req(input$region_chrom)
      chroms_use <- input$region_chrom
    } else {
      req(input$region_chrom, input$region_start, input$region_end)
      chroms_use <- input$region_chrom
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

    stage_names <- resolve_selected_stages(input$stages_browser)

    withProgress(message = paste0("Extracting signal (", length(chroms_use), " chromosome(s))…"), {
      compute_across_stages(stage_names, function(st) {
        if (identical(input$chrom_scope, "single")) {
          region_gr <- GRanges(seqnames = chroms_use,
                               ranges   = IRanges(start = input$region_start,
                                                  end   = input$region_end))
        } else {
          # Full length of each selected chromosome, per this stage's own seqinfo
          # (the two stages' ChromstaR objects need not share identical seqlengths)
          lens <- seqlengths(seqinfo(st$hmm$bins))[chroms_use]
          lens[is.na(lens)] <- 1e9  # fallback if seqlengths missing / chrom absent here
          region_gr <- GRanges(seqnames = chroms_use,
                               ranges   = IRanges(start = 1, end = lens))
        }

        extract_signal_region(
          hmm        = st$hmm,
          region_gr  = region_gr,
          marks      = marks_use,
          conditions = stage_conditions_kept(st),
          genes_gr   = rv$genes,
          bin_scope  = scope
        )
      })
    })
  })

  build_browser_plot <- reactive({
    df <- browser_data()
    req(df)

    n_chroms_shown <- length(unique(df$chr))
    multi_chrom    <- n_chroms_shown > 1
    multi_stage    <- length(unique(df$stage)) > 1

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

    pal <- condition_palette()
    p <- ggplot(df, aes(x = position, y = signal,
                        colour = condition, fill = condition)) +
      geom_area(alpha = 0.35, position = "identity") +
      geom_line(linewidth = 0.7) +
      scale_x_continuous(labels = label_comma()) +
      scale_colour_manual(values = pal) +
      scale_fill_manual(values   = pal) +
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
             title  = paste0(n_chroms_shown, " chromosomes — ", input$bin_scope_browser, " bins",
                             if (multi_stage) paste0("  |  ", paste(unique(df$stage), collapse = " vs ")) else ""))
      # Facet by chromosome AND mark (and stage, if two stages are shown);
      # condition stays as colour
      row_terms <- c("mark", if (multi_stage) "stage")
      if (input$compare_mode_browser == "facet") {
        facet_formula <- as.formula(paste(paste(row_terms, collapse = " + "), "~ chr"))
        p <- p + facet_grid(facet_formula, scales = "free")
      } else {
        wrap_formula <- as.formula(paste("chr", if (multi_stage) "+ stage" else "", "~ mark"))
        p <- p + facet_wrap(wrap_formula, scales = "free", ncol = length(unique(df$mark)))
      }
    } else {
      p <- p +
        labs(x = paste("Position on", unique(df$chr)[1]),
             y = "RPKM (averaged across replicates)",
             colour = "Condition", fill = "Condition",
             title  = if (identical(input$chrom_scope, "single")) {
               paste0(input$region_chrom, ":",
                     format(input$region_start, big.mark = ","), "–",
                     format(input$region_end,   big.mark = ","),
                     if (multi_stage) paste0("  |  ", paste(unique(df$stage), collapse = " vs ")) else "")
             } else {
               paste0(unique(df$chr)[1], " (full length)")
             })
      col_terms <- character(0)
      if (multi_stage)                            col_terms <- c(col_terms, "stage")
      if (input$compare_mode_browser == "facet")   col_terms <- c(col_terms, "condition")
      if (length(col_terms) == 0) {
        p <- p + facet_wrap(~ mark, scales = "free_y", ncol = 2)
      } else {
        facet_formula <- as.formula(paste("mark ~", paste(col_terms, collapse = " + ")))
        p <- p + facet_grid(facet_formula, scales = "free_y")
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

  output$browser_plot <- renderPlot({ apply_plot_customization(build_browser_plot(), input, "br_") })

  output$dl_browser <- downloadHandler(
    filename = function() paste0("region_browser_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = apply_plot_customization(build_browser_plot(), input, "br_"),
             width = 16, height = 12, device = "pdf")
    }
  )

  # ---- DATA TABLE ------------------------------------------------------------
  # Bin -> gene/zone annotation, extracted as a reusable function so it can be
  # cached once per stage instead of just once globally.
  annotate_bins <- function(bins, genes) {
    g_names <- mcols(genes)$gene_name
    s_g     <- start(genes)
    e_g     <- end(genes)
    str_g   <- as.character(strand(genes))
    glen    <- e_g - s_g

    tss <- ifelse(str_g == "-", e_g, s_g)
    tes <- ifelse(str_g == "-", s_g, e_g)
    chr_g <- as.character(seqnames(genes))

    make_zone_gr <- function(starts, ends, names_vec, chr_vec) {
      ok <- starts >= 1 & ends >= 1 & starts <= ends
      GRanges(
        seqnames  = chr_vec[ok],
        ranges    = IRanges(pmax(1L, starts[ok]), pmax(1L, ends[ok])),
        gene_name = names_vec[ok],
        zone_idx  = which(ok)
      )
    }

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

    zone_priority <- c(
      TSS_pm200 = 1L, gene_body_5prime = 2L, gene_body_middle = 3L,
      gene_body_3prime = 4L, upstream_1000_0 = 5L, upstream_2000_1000 = 6L,
      TES_0_1000 = 7L, TES_1000_2000 = 8L
    )

    n_bins    <- length(bins)
    best_pri  <- rep(9L,              n_bins)
    best_gene <- rep(NA_character_,   n_bins)
    best_zone <- rep("Intergenic",    n_bins)

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

    if (all(best_zone == "Intergenic")) {
      # Nothing overlapped at all — almost always a chromosome-naming mismatch
      # between the ChromstaR object and the gene table (e.g. "1" vs "chr1"),
      # or coordinates from a different genome assembly.
      showNotification(
        paste0(
          "Gene annotation matched 0 bins — every bin was labelled Intergenic. ",
          "Check that the chromosome names agree. Object: ",
          paste(utils::head(unique(as.character(seqnames(bins))), 3), collapse = ", "),
          " | gene file: ",
          paste(utils::head(unique(chr_g), 3), collapse = ", ")
        ),
        type = "warning", duration = 20
      )
    }

    data.frame(
      gene_name    = best_gene,
      genomic_zone = zone_labels[best_zone],
      stringsAsFactors = FALSE
    )
  }

  # Gene annotation is computed once per stage and cached
  observe({
    req(rv$hmm, rv$genes)
    if (!is.null(rv$bin_annotation)) return()
    withProgress(message = "Annotating bins for Stage A (one-time, ~5 sec)…", {
      rv$bin_annotation <- annotate_bins(rv$hmm$bins, rv$genes)
    })
  })

  observe({
    req(rv$hmm_b, rv$genes)
    if (!is.null(rv$bin_annotation_b)) return()
    withProgress(message = "Annotating bins for Stage B (one-time, ~5 sec)…", {
      rv$bin_annotation_b <- annotate_bins(rv$hmm_b$bins, rv$genes)
    })
  })

  # Build the per-bin state table for one stage
  build_stage_table <- function(hmm, bin_annotation, marks, conditions) {
    bins_df <- expand_bins_df(hmm)
    state_cols <- intersect(
      c("chr", "start", "end", "transition.group", "state",
        grep("^combination\\.", colnames(bins_df), value = TRUE),
        "differential.score", "maxPostInPeak"),
      colnames(bins_df)
    )

    df <- cbind(
      bins_df[, state_cols, drop = FALSE],
      bin_annotation
    )
    df$width  <- df$end - df$start + 1
    df$bin_id <- paste0(df$chr, ":", df$start, "-", df$end)

    rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
    # Emit the rpkm_<mark>_<condition> columns in the user's life-cycle order,
    # and only for the conditions they kept in the Condition/Life-cycle box.
    kept <- intersect(ordered_conditions(), conditions)
    if (length(kept) > 0) conditions <- kept
    for (cond in conditions) {
      for (mark in marks) {
        pattern <- paste0("counts.rpkm.", mark, ".", cond)
        cols    <- grep(pattern, rpkm_cols_all, value = TRUE, fixed = TRUE)
        if (length(cols) == 0) next
        colname <- paste0("rpkm_", mark, "_", cond)
        df[[colname]] <- round(rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE), 4)
      }
    }
    df
  }

  full_table_df <- reactive({
    req(rv$hmm, rv$bin_annotation)

    stage_names <- resolve_selected_stages(input$stages_table)
    stages      <- active_stages()
    annots      <- list(); annots[[stage_a_label()]] <- rv$bin_annotation
    if (!is.null(rv$hmm_b)) annots[[stage_b_label()]] <- rv$bin_annotation_b

    tables <- list()
    for (nm in stage_names) {
      st <- stages[[nm]]
      if (is.null(st) || is.null(annots[[nm]])) next
      tdf <- build_stage_table(st$hmm, annots[[nm]], st$marks, stage_conditions_kept(st))
      if (!is.null(rv$hmm_b)) tdf$stage <- nm  # only add the column when a 2nd object is loaded at all
      tables[[nm]] <- tdf
    }
    if (length(tables) == 0) return(NULL)
    df <- bind_rows(tables)

    front <- c("stage", "bin_id", "chr", "start", "end", "width", "gene_name", "genomic_zone")
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
  #' Generic, N-condition differential peak counting.
  #'
  #' Works on hmm$segments (the pre-merged chromatin segments, NOT the raw 200bp
  #' bins). Segments are first filtered globally by differential.score and width;
  #' then, for each requested PAIR of conditions, only the segments whose state
  #' actually differs between those two conditions are kept, and for every mark
  #' we count the segments carrying that mark in condition A but not B, and in B
  #' but not A.
  #'
  #' The object may hold any number of conditions (2, 5, ...); `conds_sel` says
  #' which ones to use and `mode` says how to pair them up:
  #'   "pairwise"  — every combination of the selected conditions
  #'   "reference" — every selected condition against `ref`
  #'
  #' Returns a long data.frame: comparison, mark, present_in, absent_in,
  #' direction, n_regions, n_changed.

  # A ChromstaR combination string looks like "[H3K4me3+H3K9ac]". Matching with a
  # plain grepl() would let one mark name match inside another (e.g. H3K27me3 vs
  # a hypothetical H3K27me); anchoring on the "[", "+" and "]" separators makes
  # the test exact while staying vectorised (and therefore fast on 100k+ rows).
  mark_in_combination <- function(x, mark) {
    # escape every non-word character so a mark name is matched literally
    m <- gsub("(\\W)", "\\\\\\1", mark, perl = TRUE)
    grepl(paste0("(^|\\[|\\+)", m, "($|\\]|\\+)"), x, perl = TRUE)
  }

  compute_differential_peaks <- function(hmm, marks, score_thresh, width_thresh,
                                         conds_sel = NULL, mode = "pairwise",
                                         ref = NULL, stage_label = NULL,
                                         genes_gr = NULL) {

    tag <- if (is.null(stage_label)) "" else paste0(" [", stage_label, "]")

    cat("\n--- compute_differential_peaks()", tag, "---\n")
    cat("Object class:", paste(class(hmm), collapse = ", "), "\n")
    cat("Has $segments:", !is.null(hmm$segments), "\n")

    if (is.null(hmm$segments)) {
      showNotification(
        paste0("hmm$segments not found", tag,
               " — differential peaks requires merged segments."),
        type = "error", duration = 10
      )
      return(NULL)
    }

    segs <- hmm$segments

    # Optional gene scope: keep only segments overlapping the chosen genes.
    if (!is.null(genes_gr) && length(genes_gr) > 0) {
      keep_idx <- unique(queryHits(findOverlaps(segs, genes_gr, ignore.strand = TRUE)))
      if (length(keep_idx) == 0) {
        showNotification(
          paste0("None of the selected genes overlap any chromatin segment", tag, "."),
          type = "warning", duration = 12
        )
        return(NULL)
      }
      segs <- segs[keep_idx]
      cat("Restricted to", length(genes_gr), "genes ->", length(segs), "segments\n")
    }

    segs_df <- as.data.frame(segs)

    cat("Segments:", nrow(segs_df), "rows\n")

    # --------------------------------------------------------------------------
    # Detect condition columns automatically — any number of them
    # --------------------------------------------------------------------------
    combination_cols <- grep("^combination\\.", colnames(segs_df), value = TRUE)

    if (length(combination_cols) < 2) {
      showNotification(
        paste0("Need at least 2 conditions in the object", tag, ", found ",
               length(combination_cols), "."),
        type = "error", duration = 12
      )
      return(NULL)
    }

    all_conds <- sub("^combination\\.", "", combination_cols)
    names(combination_cols) <- all_conds

    cat("Conditions in object:", paste(all_conds, collapse = ", "), "\n")

    conds <- if (is.null(conds_sel) || length(conds_sel) == 0) {
      all_conds
    } else {
      intersect(conds_sel, all_conds)
    }

    if (length(conds) < 2) {
      showNotification(
        paste0("Select at least 2 conditions present in this object", tag,
               ". Available: ", paste(all_conds, collapse = ", ")),
        type = "error", duration = 12
      )
      return(NULL)
    }

    cat("Conditions used:", paste(conds, collapse = ", "), "| mode:", mode, "\n")

    # --------------------------------------------------------------------------
    # Required columns check
    # --------------------------------------------------------------------------
    required_cols <- c("differential.score", "width")
    missing_cols  <- setdiff(required_cols, colnames(segs_df))

    if (length(missing_cols) > 0) {
      showNotification(
        paste0("Missing columns", tag, ": ", paste(missing_cols, collapse = ", ")),
        type = "error", duration = 12
      )
      return(NULL)
    }

    # --------------------------------------------------------------------------
    # Which pairs of conditions to compare
    # --------------------------------------------------------------------------
    pairs <- if (identical(mode, "reference")) {
      r <- if (!is.null(ref) && ref %in% conds) ref else conds[1]
      lapply(setdiff(conds, r), function(o) c(r, o))
    } else {
      utils::combn(conds, 2, simplify = FALSE)
    }

    # --------------------------------------------------------------------------
    # Global filter (score + width), applied once for all pairs
    # --------------------------------------------------------------------------
    keep <- segs_df$differential.score >= score_thresh & segs_df$width >= width_thresh
    keep[is.na(keep)] <- FALSE

    if (!any(keep)) {
      score_range <- paste0(
        round(min(segs_df$differential.score, na.rm = TRUE), 4), " to ",
        round(max(segs_df$differential.score, na.rm = TRUE), 4)
      )
      showNotification(
        paste0("No segments passed the filters", tag, " (score>=", score_thresh,
               ", width>=", width_thresh, "). Score range in this object: ", score_range),
        type = "warning", duration = 12
      )
      return(NULL)
    }

    kept <- segs_df[keep, , drop = FALSE]
    cat("Segments passing score/width:", nrow(kept), "\n")

    # --------------------------------------------------------------------------
    # Per pair: keep segments whose state differs between those two conditions,
    # then count mark presence in each direction
    # --------------------------------------------------------------------------
    results <- list()

    # Segments that actually change state in at least one compared pair — i.e.
    # exactly the set the bar chart is built from. With two conditions this is
    # the single pair's count; with more, a segment counts once even if it
    # changes in several pairs.
    changed_union <- integer(0)

    for (pr in pairs) {
      a <- pr[1]
      b <- pr[2]

      combo_a <- as.character(kept[[combination_cols[[a]]]])
      combo_b <- as.character(kept[[combination_cols[[b]]]])

      changed <- which(combo_a != combo_b)
      changed_union <- union(changed_union, changed)
      cat("  ", a, "vs", b, "— segments changing state:", length(changed), "\n")

      if (length(changed) == 0) next

      fa <- combo_a[changed]
      fb <- combo_b[changed]

      for (m in marks) {
        in_a <- mark_in_combination(fa, m)
        in_b <- mark_in_combination(fb, m)

        results[[length(results) + 1L]] <- data.frame(
          comparison = paste0(a, " vs ", b),
          mark       = m,
          present_in = c(a, b),
          absent_in  = c(b, a),
          direction  = c(paste0(a, "-not-", b), paste0(b, "-not-", a)),
          n_regions  = c(sum(in_a & !in_b), sum(!in_a & in_b)),
          n_changed  = length(changed),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(results) == 0) {
      showNotification(
        paste0("No segments changed chromatin state between the selected conditions",
               tag, " at these thresholds. Try lowering the score or width filter."),
        type = "warning", duration = 12
      )
      return(NULL)
    }

    out <- dplyr::bind_rows(results)

    cat("Differential regions used in the bar chart:", length(changed_union), "\n")

    attr(out, "total_filtered")  <- nrow(kept)
    attr(out, "n_differential")  <- length(changed_union)
    attr(out, "conditions")      <- conds

    out
  }

  # ------------------------------------------------------------------------------
  # Condition / reference pickers for this tab
  # ------------------------------------------------------------------------------

  # Conditions available on this tab = those present in the selected stage(s),
  # restricted to (and ordered by) the user's life-cycle order from the Load tab.
  diff_available_conditions <- reactive({
    stages      <- active_stages()
    stage_names <- resolve_selected_stages(input$stages_diffpeaks)
    conds <- unique(unlist(lapply(stage_names, function(nm) stages[[nm]]$conditions)))
    conds <- conds[!is.na(conds)]
    ord   <- ordered_conditions()
    c(intersect(ord, conds), setdiff(conds, ord))
  })

  output$gene_selector_diff <- renderUI({
    req(rv$genes)
    nms <- gene_names_reactive()
    tagList(
      tags$small(style = "color:#888",
        paste0(length(nms), " genes in the annotation — paste a list (one per row). ",
               "Only chromatin segments overlapping these genes are counted.")),
      br(), br(),
      tags$textarea(
        id = "selected_genes_text_diff",
        rows = 10,
        style = "width:100%; height:200px; font-family: monospace; font-size: 13px; resize: vertical;",
        placeholder = "Paste gene IDs here, one per line…",
        ""
      ),
      uiOutput("gene_match_summary_diff")
    )
  })

  selected_genes_parsed_diff <- reactive({
    req(input$selected_genes_text_diff)
    raw <- strsplit(input$selected_genes_text_diff, "[\r\n,;\t]+")[[1]]
    raw <- trimws(raw)
    raw[nchar(raw) > 0]
  })

  output$gene_match_summary_diff <- renderUI({
    req(rv$genes)
    pasted <- selected_genes_parsed_diff()
    if (length(pasted) == 0) {
      return(tags$small(style = "color:#888", "No genes entered yet."))
    }
    nms       <- gene_names_reactive()
    matched   <- pasted[pasted %in% nms]
    unmatched <- setdiff(pasted, nms)
    tagList(
      br(),
      tags$small(style = "color:#28a745",
        paste0("\u2713 ", length(matched), " of ", length(pasted), " gene names matched.")),
      if (length(unmatched) > 0) {
        tags$div(tags$small(style = "color:#dc3545",
          paste0("\u2717 ", length(unmatched), " not found, e.g.: ",
                 paste(utils::head(unmatched, 5), collapse = ", "),
                 if (length(unmatched) > 5) "\u2026" else "")))
      }
    )
  })

  observeEvent(input$select_all_genes_diff, {
    req(rv$genes)
    updateTextAreaInput(session, "selected_genes_text_diff",
                        value = paste(gene_names_reactive(), collapse = "\n"))
  })

  observeEvent(input$clear_all_genes_diff, {
    updateTextAreaInput(session, "selected_genes_text_diff", value = "")
  })

  # Genes the Differential Peaks tab should be restricted to — NULL means the
  # whole genome (every segment), which is the default.
  diff_genes_selected <- reactive({
    if (!identical(input$diff_gene_scope, "selected")) return(NULL)
    req(rv$genes)
    nms <- gene_names_reactive()
    sel <- selected_genes_parsed_diff()
    sel <- sel[sel %in% nms]
    if (length(sel) == 0) return(NULL)
    rv$genes[nms %in% sel]
  })

  output$diff_condition_ui <- renderUI({
    conds <- diff_available_conditions()
    if (length(conds) == 0) {
      return(tags$em(style = "color:#888",
                     "Load a ChromstaR object to choose conditions."))
    }
    prev <- intersect(isolate(input$diff_conditions) %||% character(0), conds)
    sel  <- if (length(prev) >= 2) prev else utils::head(conds, min(5L, length(conds)))
    selectizeInput(
      "diff_conditions", NULL,
      choices  = conds,
      selected = sel,
      multiple = TRUE,
      width    = "100%",
      options  = list(maxItems = 5, plugins = list("drag_drop", "remove_button"))
    )
  })

  output$diff_ref_ui <- renderUI({
    sel <- intersect(input$diff_conditions %||% character(0), diff_available_conditions())
    if (length(sel) < 2) return(NULL)
    selectInput("diff_ref", "Reference condition",
                choices = sel, selected = sel[1], width = "100%")
  })

  # ------------------------------------------------------------------------------
  # Reactive: compute differential peaks
  # ------------------------------------------------------------------------------

  diffpeaks_data <- eventReactive(input$run_diffpeaks, {

    req(rv$hmm, rv$marks)

    stage_names <- resolve_selected_stages(input$stages_diffpeaks)
    conds_sel   <- intersect(input$diff_conditions %||% character(0),
                             diff_available_conditions())
    mode        <- input$diff_mode %||% "pairwise"
    ref         <- input$diff_ref
    genes_gr    <- diff_genes_selected()

    if (identical(input$diff_gene_scope, "selected") &&
        (is.null(genes_gr) || length(genes_gr) == 0)) {
      showNotification(
        "Gene scope is \"Only selected genes\" but no valid gene names were entered.",
        type = "error", duration = 8
      )
      return(NULL)
    }

    if (length(conds_sel) < 2) {
      showNotification("Pick at least 2 conditions to compare.",
                       type = "error", duration = 8)
      return(NULL)
    }

    tryCatch({

      withProgress(message = "Computing differential peaks…", {

        compute_across_stages(stage_names, function(st) {
          df <- compute_differential_peaks(
            hmm          = st$hmm,
            marks        = st$marks,
            score_thresh = input$diff_score_thresh,
            width_thresh = input$diff_width_thresh,
            conds_sel    = conds_sel,
            mode         = mode,
            ref          = ref,
            genes_gr     = genes_gr
          )
          if (is.null(df)) return(NULL)
          # attr()s don't survive bind_rows() across stages, so carry the
          # per-stage bookkeeping (total filtered segments, which conditions
          # this stage actually contributed) as ordinary columns instead.
          df$total_filtered   <- attr(df, "total_filtered")
          df$n_differential   <- attr(df, "n_differential")
          df$stage_conditions <- paste(attr(df, "conditions"), collapse = ", ")
          df
        })

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

  # Shared prep so the plot and the sizing logic agree on facet counts.
  diffpeaks_prepared <- reactive({
    df <- diffpeaks_data()
    req(df)

    ord <- ordered_conditions()

    # marks ordered by total number of differential regions (biggest at top)
    mark_order <- df %>%
      dplyr::group_by(mark) %>%
      dplyr::summarise(total = sum(n_regions), .groups = "drop") %>%
      dplyr::arrange(desc(total)) %>%
      dplyr::pull(mark)

    df$mark <- factor(df$mark, levels = rev(mark_order))

    # Bars are coloured by the condition the mark is present in, so the colours
    # match the rest of the app no matter how many pairs are on screen.
    pres <- unique(as.character(df$present_in))
    lv   <- c(intersect(ord, pres), setdiff(pres, ord))
    df$present_in <- factor(as.character(df$present_in), levels = lv)

    # Facets follow the life-cycle order of the first, then second condition.
    comps <- unique(as.character(df$comparison))
    r1 <- match(sub(" vs .*$", "", comps), ord)
    r2 <- match(sub("^.* vs ", "", comps), ord)
    r1[is.na(r1)] <- length(ord) + 1L
    r2[is.na(r2)] <- length(ord) + 1L
    comps <- comps[order(r1, r2, comps)]
    df$comparison <- factor(as.character(df$comparison), levels = comps)

    list(
      df          = df,
      n_comps     = length(comps),
      multi_stage = "stage" %in% colnames(df) && length(unique(df$stage)) > 1,
      n_stages    = if ("stage" %in% colnames(df)) length(unique(df$stage)) else 1L,
      levels      = lv
    )
  })

  build_diffpeaks_plot <- reactive({

    prep <- diffpeaks_prepared()
    df   <- prep$df

    fill_cols <- condition_palette()[prep$levels]
    names(fill_cols) <- prep$levels
    fill_cols[is.na(fill_cols)] <- "#9E9E9E"

    info <- df %>% dplyr::distinct(stage, stage_conditions, total_filtered, n_differential)
    mode_txt <- if (identical(input$diff_mode, "reference")) {
      paste0("reference = ", input$diff_ref)
    } else {
      "all pairwise"
    }
    scope_txt <- if (identical(input$diff_gene_scope, "selected")) {
      paste0(" | ", length(diff_genes_selected()), " selected genes")
    } else {
      " | whole genome"
    }

    subtitle_txt <- paste0(
      paste0(info$stage, " (", info$stage_conditions, ", ",
             format(info$n_differential, big.mark = ","),
             " differential regions)", collapse = "   |   "),
      "\nscore>=", input$diff_score_thresh,
      ", width>=", input$diff_width_thresh, "bp | ", mode_txt, scope_txt
    )

    p <- ggplot(df, aes(x = mark, y = n_regions, fill = present_in)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = fill_cols, drop = FALSE) +
      labs(
        x = "Histone Mark",
        y = "Number of Regions",
        fill = "Mark present in",
        title = "Differential Peaks per Histone Mark",
        subtitle = subtitle_txt
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")
      )

    if (prep$multi_stage) {
      p <- p + facet_grid(stage ~ comparison)
    } else {
      p <- p + facet_wrap(~ comparison, ncol = min(3L, max(1L, prep$n_comps)))
    }

    p
  })

  # With up to 10 pairwise panels (5 conditions) a fixed 550px is unreadable,
  # so the canvas grows with the number of facet rows.
  diffpeaks_plot_height <- reactive({
    prep <- diffpeaks_prepared()
    if (prep$multi_stage) {
      max(450L, as.integer(300 * prep$n_stages) + 120L)
    } else {
      rows <- ceiling(prep$n_comps / min(3L, max(1L, prep$n_comps)))
      max(550L, as.integer(300 * rows) + 120L)
    }
  })

  # ------------------------------------------------------------------------------
  # Outputs
  # ------------------------------------------------------------------------------

  output$diffpeaks_plot_container <- renderUI({
    h <- tryCatch(diffpeaks_plot_height(), error = function(e) 550L)
    plotOutput("diffpeaks_plot", height = paste0(h, "px"))
  })

  output$diffpeaks_plot <- renderPlot({
    apply_plot_customization(build_diffpeaks_plot(), input, "dp_")
  })

  output$dl_diffpeaks_plot <- downloadHandler(
    filename = function() {
      paste0("differential_peaks_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      h <- tryCatch(diffpeaks_plot_height() / 70, error = function(e) 7)
      ggplot2::ggsave(file, apply_plot_customization(build_diffpeaks_plot(), input, "dp_"),
             width = 12, height = max(7, min(30, h)), device = "pdf", limitsize = FALSE)
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

  # ---- Per-mark summary + ranked figures (Differential Peaks tab) ------------
  # Built from the same differential regions the bar chart uses, so every
  # percentage here has n_differential as its denominator.
  diffpeaks_mark_summary <- reactive({
    df <- diffpeaks_data()
    req(df)
    rows <- list()
    for (st in unique(as.character(df$stage))) {
      sdf   <- df[as.character(df$stage) == st, , drop = FALSE]
      ndiff <- suppressWarnings(max(sdf$n_differential, na.rm = TRUE))
      if (!is.finite(ndiff) || ndiff <= 0) next
      for (cmp in unique(as.character(sdf$comparison))) {
        cdf <- sdf[as.character(sdf$comparison) == cmp, , drop = FALSE]
        a_c <- sub(" vs .*$", "", cmp)
        b_c <- sub("^.* vs ", "", cmp)
        for (m in unique(as.character(cdf$mark))) {
          mdf    <- cdf[as.character(cdf$mark) == m, , drop = FALSE]
          gained <- sum(mdf$n_regions[as.character(mdf$present_in) == b_c], na.rm = TRUE)
          lost   <- sum(mdf$n_regions[as.character(mdf$present_in) == a_c], na.rm = TRUE)
          rows[[length(rows) + 1L]] <- data.frame(
            stage          = st,
            comparison     = cmp,
            mark           = m,
            gained         = gained,
            lost           = lost,
            total_changed  = gained + lost,
            n_differential = ndiff,
            cond_a         = a_c,
            cond_b         = b_c,
            stringsAsFactors = FALSE)
        }
      }
    }
    if (length(rows) == 0) return(NULL)
    out <- dplyr::bind_rows(rows)
    out$pct_changed <- round(out$total_changed / out$n_differential * 100, 2)
    out$pct_gained  <- round(out$gained        / out$n_differential * 100, 2)
    out$pct_lost    <- round(out$lost          / out$n_differential * 100, 2)
    out$direction   <- ifelse(out$gained > out$lost,
                              paste0("Gained in ", out$cond_b),
                              paste0("Lost in ",   out$cond_b))
    out$panel       <- paste0(out$stage, " — ", out$comparison)
    out[order(out$stage, out$comparison, -out$total_changed), ]
  })

  # Spells out why the per-mark percentages sum to more than 100%.
  output$diffpeaks_overlap_note <- renderUI({
    out <- diffpeaks_mark_summary()
    req(out)
    parts <- vapply(split(out, out$panel, drop = TRUE), function(g) {
      sprintf("%s: %.2f marks change per differential region",
              g$panel[1], sum(g$total_changed) / g$n_differential[1])
    }, character(1))
    tags$div(style = "color:#555; font-size:12px; line-height:1.5;",
      tags$b("These percentages overlap and will not sum to 100%. "),
      "Most differential regions have several marks changing at once, and a region is counted once for every mark that changes in it. ",
      tags$br(), paste(parts, collapse = "   |   "))
  })

  output$diffpeaks_mark_table <- renderDT({
    out <- diffpeaks_mark_summary()
    req(out)
    disp <- out[, c("stage", "comparison", "mark", "gained", "lost",
                    "total_changed", "pct_gained", "pct_lost", "pct_changed",
                    "n_differential", "direction")]
    names(disp) <- c("Stage", "Comparison", "Mark", "Gained", "Lost",
                     "Total changed", "% gained", "% lost", "% of differential",
                     "Differential regions", "Direction")
    DT::datatable(disp, rownames = FALSE,
                  options = list(pageLength = 16, dom = "tip", scrollX = TRUE))
  })

  output$dl_diffpeaks_mark_table <- downloadHandler(
    filename = function() paste0("differential_peaks_per_mark_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      out <- diffpeaks_mark_summary()
      req(out)
      openxlsx::write.xlsx(out[, setdiff(colnames(out), "panel")], file, overwrite = TRUE)
    }
  )

  # ---- Ranked figure: filtered differential regions --------------------------
  build_diffpeaks_rank_filtered <- reactive({
    out <- diffpeaks_mark_summary()
    req(out)
    nd <- format(max(out$n_differential, na.rm = TRUE), big.mark = ",")
    mark_ranking_ggplot(
      out,
      title    = "Histone Marks Ranked by Magnitude of Change",
      subtitle = paste0("Differential regions only (score/width filtered) | n = ", nd),
      ylab        = "Differential Regions Changed",
      style       = input$dp_rank_style %||% "stacked",
      cond_colors = condition_palette())
  })

  output$diffpeaks_rank_filtered <- renderPlot({
    apply_plot_customization(build_diffpeaks_rank_filtered(), input, "dp_")
  })

  output$dl_diffpeaks_rank_filtered <- downloadHandler(
    filename = function() paste0("marks_ranked_differential_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file,
        apply_plot_customization(build_diffpeaks_rank_filtered(), input, "dp_"),
        width = 10, height = 6, device = "pdf")
    }
  )

  # ---- Ranked figure: genome-wide, from the $frequencies table ---------------
  diffpeaks_pairs_for <- function(conds) {
    if (identical(input$diff_mode, "reference")) {
      r <- if (!is.null(input$diff_ref) && input$diff_ref %in% conds) input$diff_ref else conds[1]
      lapply(setdiff(conds, r), function(o) c(r, o))
    } else {
      utils::combn(conds, 2, simplify = FALSE)
    }
  }

  diffpeaks_rank_genome_data <- reactive({
    stages <- active_stages()
    sel    <- resolve_selected_stages(input$stages_diffpeaks)
    req(length(sel) > 0)
    rows <- list(); missing_freq <- character(0)
    for (nm in sel) {
      st <- stages[[nm]]
      if (is.null(st)) next
      fq <- tryCatch(st$hmm$frequencies, error = function(e) NULL)
      if (is.null(fq)) { missing_freq <- c(missing_freq, nm); next }
      conds <- stage_conditions_kept(st)
      if (length(conds) < 2) next
      for (pr in diffpeaks_pairs_for(conds)) {
        rk <- compute_mark_ranking_freq(fq, st$marks, pr[1], pr[2])
        if (is.null(rk)) next
        rk$stage      <- nm
        rk$comparison <- paste0(pr[1], " vs ", pr[2])
        rk$panel      <- paste0(nm, " — ", rk$comparison)
        rows[[length(rows) + 1L]] <- rk
      }
    }
    list(data    = if (length(rows)) dplyr::bind_rows(rows) else NULL,
         missing = missing_freq)
  })

  output$diffpeaks_rank_genome_note <- renderUI({
    r <- diffpeaks_rank_genome_data()
    if (length(r$missing) == 0) return(NULL)
    tags$div(style = "color:#b8860b; font-size:12px; margin-bottom:8px;",
      tags$b("No combination-frequency table found for: "),
      paste(r$missing, collapse = ", "), ". ",
      "This genome-wide view reads the object's $frequencies slot. Re-load the chromstaR object with this version of the app — earlier versions discarded that slot when trimming the object for memory.")
  })

  build_diffpeaks_rank_genome <- reactive({
    r <- diffpeaks_rank_genome_data()
    d <- r$data
    req(!is.null(d), nrow(d) > 0)
    tot <- format(max(d$total_domains, na.rm = TRUE), big.mark = ",")
    mark_ranking_ggplot(
      d,
      title    = "Histone Marks Ranked by Magnitude of Change",
      subtitle = paste0("Genome-wide, unfiltered | Total domains = ", tot),
      ylab        = "Total Domains Changed",
      style       = input$dp_rank_style %||% "stacked",
      cond_colors = condition_palette())
  })

  output$diffpeaks_rank_genome <- renderPlot({
    apply_plot_customization(build_diffpeaks_rank_genome(), input, "dp_")
  })

  output$dl_diffpeaks_rank_genome <- downloadHandler(
    filename = function() paste0("marks_ranked_genomewide_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file,
        apply_plot_customization(build_diffpeaks_rank_genome(), input, "dp_"),
        width = 10, height = 6, device = "pdf")
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

    mode        <- input$gsc_mode
    stage_names <- resolve_selected_stages(input$stages_gsc)
    stages      <- active_stages()

    sel_b <- character(0)
    if (identical(mode, "manual")) {
      sel_b <- gsc_genes_b_parsed(); sel_b <- sel_b[sel_b %in% nms]
      if (length(sel_b) == 0) {
        showNotification(
          "Gene set B needs at least one valid, matched gene (or switch to the random background mode).",
          type = "warning", duration = 6
        )
        return(NULL)
      }
    }

    per_stage <- list()
    err <- NULL

    for (nm in stage_names) {
      st <- stages[[nm]]
      if (is.null(st)) next

      if (identical(mode, "permutation")) {
        res <- tryCatch({
          withProgress(message = paste0("Building gene universe (", nm, ")…"), {
            nms_all  <- gene_names_reactive()
            universe <- compute_all_gene_promoter_posteriors(
              hmm                 = st$hmm,
              all_genes_gr        = rv$genes,
              gene_labels         = nms_all,
              upstream            = input$gsc_upstream,
              downstream          = input$gsc_downstream,
              average_replicates  = identical(input$gsc_replicate_mode, "average")
            )
            if (!is.null(universe$error)) return(list(error = paste0(nm, ": ", universe$error)))

            incProgress(0.3, detail = paste0(nm, ": running ", input$gsc_n_perm, " random draws…"))
            perm_res <- compute_gene_set_permutation(
              all_gene_scores = universe$gene_scores,
              post_cols       = universe$post_cols,
              genes_A_labels  = sel_a,
              n_perm          = input$gsc_n_perm,
              summary_stat    = input$gsc_summary_stat
            )
            if (!is.null(perm_res$error)) return(list(error = paste0(nm, ": ", perm_res$error)))
            perm_res
          })
        }, error = function(e) list(error = paste0("Error (", nm, "): ", conditionMessage(e))))
      } else {
        genes_A_gr <- rv$genes[nms %in% sel_a]
        mcols(genes_A_gr)$gene_label <- nms[nms %in% sel_a]
        genes_B_gr <- rv$genes[nms %in% sel_b]
        mcols(genes_B_gr)$gene_label <- nms[nms %in% sel_b]

        res <- tryCatch({
          withProgress(message = paste0("Comparing gene sets (", nm, ")…"), {
            compute_gene_set_posteriors(
              hmm                 = st$hmm,
              genes_A_gr          = genes_A_gr,
              genes_B_gr          = genes_B_gr,
              upstream            = input$gsc_upstream,
              downstream          = input$gsc_downstream,
              average_replicates  = identical(input$gsc_replicate_mode, "average"),
              summary_stat        = input$gsc_summary_stat
            )
          })
        }, error = function(e) list(error = paste0("Error (", nm, "): ", conditionMessage(e))))
      }

      if (!is.null(res$error)) { err <- res$error; break }
      per_stage[[nm]] <- res
    }

    if (!is.null(err)) {
      showNotification(err, type = "error", duration = 12)
      return(NULL)
    }
    if (length(per_stage) == 0) return(NULL)

    # Combine a given field across stages, tagging each stage's rows with a
    # "stage" column before binding (per-object attr()s like summary_stat
    # don't survive bind_rows, so those are looked up from `input` instead).
    combine_field <- function(field) {
      parts <- lapply(names(per_stage), function(nm) {
        d <- per_stage[[nm]][[field]]
        if (is.null(d)) return(NULL)
        d$stage <- nm
        d
      })
      parts <- Filter(Negate(is.null), parts)
      if (length(parts) == 0) return(NULL)
      out <- bind_rows(parts)
      out$stage <- factor(out$stage, levels = names(per_stage))
      apply_condition_order(out)
    }

    n_info <- if (identical(mode, "permutation")) {
      data.frame(
        stage   = names(per_stage),
        n_A     = vapply(per_stage, function(x) x$n_A, numeric(1)),
        n_total = vapply(per_stage, function(x) x$n_total, numeric(1)),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        stage     = names(per_stage),
        n_genes_A = vapply(per_stage, function(x) x$n_genes_A, numeric(1)),
        n_genes_B = vapply(per_stage, function(x) x$n_genes_B, numeric(1)),
        stringsAsFactors = FALSE
      )
    }

    list(
      error       = NULL,
      mode        = mode,
      stats       = combine_field("stats"),
      long_scores = combine_field("long_scores"),
      perm_long   = combine_field("perm_long"),
      n_info      = n_info
    )
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

    # Put stage first if present, for readability when two stages are shown
    if ("stage" %in% colnames(df)) {
      df <- df[, c("stage", setdiff(colnames(df), "stage")), drop = FALSE]
    }

    datatable(df, options = list(pageLength = 10, order = list(list(which(colnames(df) == "FDR") - 1, "asc"))),
             rownames = FALSE) %>%
      formatStyle("FDR", backgroundColor = styleInterval(c(0.01, 0.05),
                                                         c("#c6efce", "#ffeb9c", "white")))
  })

  output$gsc_summary_stat_label <- renderUI({
    res <- gsc_result()
    req(res)
    stat_used <- input$gsc_summary_stat
    n_txt <- if (identical(res$mode, "permutation")) {
      paste(apply(res$n_info, 1, function(r)
        paste0(r["stage"], ": n=", r["n_A"], " (universe ", r["n_total"], ")")), collapse = "  |  ")
    } else {
      paste(apply(res$n_info, 1, function(r)
        paste0(r["stage"], ": A=", r["n_genes_A"], ", B=", r["n_genes_B"])), collapse = "  |  ")
    }
    tags$small(style = "color:#555",
      paste0("Comparison statistic: ", tools::toTitleCase(stat_used),
            " (delta and stat_A/stat_B columns use this; mean and median are both always shown where applicable).  ",
            n_txt))
  })

  build_gsc_plot <- reactive({
    res <- gsc_result()
    req(res, res$long_scores)

    df <- res$long_scores
    # Order facets by mark, then by the user's life-cycle condition order
    df$facet_label <- paste0(df$mark, " (", as.character(df$condition), ")")
    lvl_grid <- expand.grid(
      cond = intersect(ordered_conditions(), unique(as.character(df$condition))),
      mark = sort(unique(as.character(df$mark))),
      stringsAsFactors = FALSE
    )
    facet_levels <- paste0(lvl_grid$mark, " (", lvl_grid$cond, ")")
    df$facet_label <- factor(df$facet_label,
                             levels = intersect(facet_levels, unique(df$facet_label)))
    multi_stage <- "stage" %in% colnames(df) && length(unique(df$stage)) > 1

    if (identical(res$mode, "permutation")) {
      group_colors <- c("A" = "#4393c3", "Random" = "#999999")
      df$group <- factor(df$group, levels = c("A", "Random"))
      n_txt <- paste(apply(res$n_info, 1, function(r) paste0(r["stage"], ": n=", r["n_A"])), collapse = "  |  ")
      subtitle_txt <- paste0(n_txt, " genes in set A (vs. one representative random draw of the same size) — ",
                            "empirical p-value uses ", input$gsc_n_perm, " random draws total per stage (see stats table)")
    } else {
      group_colors <- c("A" = "#4393c3", "B" = "#d6604d")
      subtitle_txt <- paste(apply(res$n_info, 1, function(r)
        paste0(r["stage"], ": A=", r["n_genes_A"], " genes, B=", r["n_genes_B"], " genes")), collapse = "  |  ")
    }

    p <- ggplot(df, aes(x = group, y = posterior, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
      scale_fill_manual(values = group_colors) +
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

    if (multi_stage) {
      p <- p + facet_grid(stage ~ facet_label, scales = "free_y")
    } else {
      p <- p + facet_wrap(~ facet_label, scales = "free_y")
    }

    p
  })

  output$gsc_boxplot <- renderPlot({ apply_plot_customization(build_gsc_plot(), input, "gsc_") })

  output$dl_gsc_plot <- downloadHandler(
    filename = function() paste0("gene_set_comparison_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = apply_plot_customization(build_gsc_plot(), input, "gsc_"), width = 12, height = 8, device = "pdf")
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

    cols <- c("gene", "mark", "condition", "posterior")
    if ("stage" %in% colnames(res$long_scores)) cols <- c("stage", cols)

    df <- res$long_scores[res$long_scores$group == "A", cols]
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
