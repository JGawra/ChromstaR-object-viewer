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
options(shiny.maxRequestSize = 2000 * 1024^2)
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

APP_BUILD <- "2026-10-11a / v27: RF + SHAP region model tabs"
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
        if (inherits(e, "shiny.silent.error")) return()   # req() waiting for input
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
      # req() and validate() signal silent errors: let Shiny handle them as
      # usual (blank plot until inputs are ready) instead of reporting a failure.
      if (inherits(e, "shiny.silent.error")) stop(e)
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

# =============================================================================
# EXPRESSION-STRATIFIED PROFILES  (normalised count table + expression deciles)
# =============================================================================

#' Load a normalised count table (genes x samples).
#'
#' Expected layout: first column = gene IDs (its header may be blank), then one
#' numeric column per sample. Comma-, tab- or semicolon-separated.
#' Accepts DESeq2 normalised counts (counts(dds, normalized = TRUE)), TPM or CPM.
#' Raw counts are also read, but flagged (looks_raw = TRUE).
load_count_table <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  n_tab   <- lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE)))
  n_comma <- lengths(regmatches(first_line, gregexpr(",",  first_line, fixed = TRUE)))
  n_semi  <- lengths(regmatches(first_line, gregexpr(";",  first_line, fixed = TRUE)))
  sep <- if (n_tab > 0 && n_tab >= max(n_comma, n_semi)) "\t"
         else if (n_semi > n_comma) ";"
         else ","

  df <- read.table(path, header = TRUE, sep = sep, row.names = 1,
                   quote = "\"", comment.char = "", check.names = FALSE,
                   stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))

  if (ncol(df) < 1) stop("No sample columns found.")
  is_num <- vapply(df, function(x) is.numeric(x) || all(is.na(x)), logical(1))
  if (!all(is_num)) {
    stop("Non-numeric columns: ",
         paste(head(names(df)[!is_num], 5), collapse = ", "),
         ". Expected: first column = gene IDs, every other column = numeric counts.")
  }

  mat <- as.matrix(df)
  storage.mode(mat) <- "double"
  rownames(mat) <- rownames(df)

  dup <- duplicated(rownames(mat))
  if (any(dup)) mat <- mat[!dup, , drop = FALSE]

  vals <- as.vector(mat)
  vals <- vals[is.finite(vals)]

  list(
    mat          = mat,
    n_genes      = nrow(mat),
    samples      = colnames(mat),
    n_dup_removed = sum(dup),
    n_na         = sum(is.na(mat)),
    range        = range(vals),
    median_value = median(vals),
    looks_raw    = length(vals) > 0 && all(vals == round(vals)),
    sep          = sep
  )
}

#' Mean normalised expression per annotated gene, from the loaded count table.
#' Matches gene_id first, then gene_name. Returns NA for genes not in the table.
gene_expression_values <- function(genes_gr, count_mat, samples) {
  sub          <- count_mat[, samples, drop = FALSE]
  expr_by_row  <- rowMeans(sub, na.rm = TRUE)
  ids          <- as.character(mcols(genes_gr)$gene_id)
  nms          <- as.character(mcols(genes_gr)$gene_name)
  idx          <- match(ids, rownames(count_mat))
  miss         <- is.na(idx)
  idx[miss]    <- match(nms[miss], rownames(count_mat))
  val          <- expr_by_row[idx]
  val[is.nan(val)] <- NA
  val
}

#' Ten colours for expression deciles: blue (decile 1, lowest) to red (10, highest).
expression_decile_colours <- function() {
  colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#E0F3F8",
                     "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(10)
}

#' Gene-body profile with flanks, stratified by a per-gene group.
#'
#' Same idea as deepTools computeMatrix "scale-regions": the upstream flank is
#' split into n_flank slots, the gene body is scaled into n_body slots, and the
#' downstream flank is split into n_flank slots. Signal comes from the ChromStaR
#' bins (counts.rpkm columns), which are averaged per gene and slot first, then
#' averaged across genes within each group, so long genes do not dominate.
#'
#' genes_gr must carry mcols()[[group_col]] (integer decile 1..10).
compute_scaled_profile <- function(hmm, genes_gr, marks, conditions,
                                   upstream = 1500, downstream = 1500,
                                   n_flank = 10, n_body = 50,
                                   group_col = "expr_group") {
  if (length(genes_gr) == 0) return(NULL)
  bins        <- hmm$bins
  bins_df     <- expand_bins_df(hmm)
  bin_mid_all <- (bins_df$start + bins_df$end) / 2

  rpkm_cols_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
  if (length(rpkm_cols_all) == 0) return(NULL)

  signal_list <- list()
  for (cond in conditions) {
    for (mark in marks) {
      cols <- grep(paste0("counts.rpkm.", mark, ".", cond), rpkm_cols_all,
                   value = TRUE, fixed = TRUE)
      if (length(cols) == 0) next
      signal_list[[length(signal_list) + 1]] <- list(
        cond = cond, mark = mark,
        vec  = rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)
      )
    }
  }
  if (length(signal_list) == 0) return(NULL)

  # Strand-aware TSS (5' end). Unknown strand is treated as "+".
  str_g <- as.character(strand(genes_gr))
  str_g[!str_g %in% c("+", "-")] <- "+"
  tss   <- ifelse(str_g == "-", end(genes_gr), start(genes_gr))
  len   <- width(genes_gr)
  grp   <- mcols(genes_gr)[[group_col]]

  # Search window: wide enough to contain both flanks for every gene
  pad <- max(upstream, downstream)
  win <- GRanges(seqnames = seqnames(genes_gr),
                 ranges   = IRanges(pmax(1L, start(genes_gr) - pad),
                                    end(genes_gr) + pad))
  ov <- findOverlaps(bins, win)
  if (length(ov) == 0) return(NULL)

  bi  <- queryHits(ov)
  gi  <- subjectHits(ov)
  # Distance of the bin midpoint from the TSS, along transcription.
  # Negative = upstream, 0..len = gene body, > len = downstream.
  sgn <- ifelse(str_g[gi] == "-", -1, 1)
  d   <- (bin_mid_all[bi] - tss[gi]) * sgn
  Lg  <- len[gi]

  slot  <- rep(NA_integer_, length(d))
  m_up  <- d >= -upstream & d < 0
  m_bd  <- d >= 0 & d <= Lg
  m_dn  <- d > Lg & d <= Lg + downstream
  slot[m_up] <- pmin(floor((d[m_up] + upstream) / upstream * n_flank) + 1L, n_flank)
  slot[m_bd] <- n_flank + pmin(floor(d[m_bd] / Lg[m_bd] * n_body) + 1L, n_body)
  slot[m_dn] <- n_flank + n_body +
    pmin(floor((d[m_dn] - Lg[m_dn]) / downstream * n_flank) + 1L, n_flank)

  keep <- !is.na(slot) & !is.na(grp[gi])
  if (!any(keep)) return(NULL)
  bi   <- bi[keep]
  gi   <- gi[keep]
  slot <- slot[keep]
  grp_k <- grp[gi]

  pieces <- lapply(signal_list, function(s) {
    data.frame(expr_group = grp_k, gene = gi, slot = slot,
               signal = s$vec[bi], mark = s$mark, condition = s$cond,
               stringsAsFactors = FALSE)
  })
  long <- bind_rows(pieces)

  # Step 1: one value per gene per slot (a gene counts once per slot)
  per_gene <- long %>%
    group_by(mark, condition, expr_group, gene, slot) %>%
    summarise(signal = mean(signal, na.rm = TRUE), .groups = "drop")

  # Step 2: mean and SE across the genes of each decile
  per_gene %>%
    group_by(mark, condition, expr_group, slot) %>%
    summarise(mean_signal = mean(signal, na.rm = TRUE),
              se_signal   = sd(signal, na.rm = TRUE) / sqrt(sum(!is.na(signal))),
              n_genes     = n_distinct(gene),
              .groups     = "drop")
}

# =============================================================================
# MARK COMBINATIONS (UpSet-style figure)
# =============================================================================

#' Per-bin TRUE/FALSE matrix (bins x marks): does this bin carry each mark in
#' its chromStaR combination state for `cond`?
state_mark_matrix <- function(hmm, cond, marks) {
  combo <- get_combination_vector(hmm, cond)
  if (is.null(combo)) return(NULL)
  lv    <- if (is.factor(combo)) levels(combo) else sort(unique(as.character(combo)))
  codes <- if (is.factor(combo)) as.integer(combo) else match(as.character(combo), lv)
  sets  <- combo_mark_sets(lv)
  state_m <- vapply(marks, function(m) vapply(sets, function(s) m %in% s, logical(1)),
                    logical(length(lv)))
  state_m <- matrix(state_m, nrow = length(lv), dimnames = list(NULL, marks))
  out <- state_m[codes, , drop = FALSE]
  out[is.na(out)] <- FALSE
  out
}

#' Per-gene TRUE/FALSE matrix (genes x marks): a gene is "marked" by a mark when
#' any bin in its window carries that mark.
#' window = "promoter": TSS +/- up/down bp.   window = "body": whole gene body.
gene_mark_presence <- function(hmm, cond, genes_gr, marks,
                               window = "promoter", up = 1000, down = 1000) {
  Mb <- state_mark_matrix(hmm, cond, marks)
  if (is.null(Mb)) return(NULL)
  win <- if (identical(window, "body")) genes_gr
         else suppressWarnings(promoters(genes_gr, upstream = up, downstream = down))
  ov  <- suppressWarnings(findOverlaps(win, hmm$bins, ignore.strand = TRUE))
  gi  <- queryHits(ov)
  bi  <- subjectHits(ov)
  P <- matrix(FALSE, nrow = length(genes_gr), ncol = length(marks),
              dimnames = list(NULL, marks))
  for (k in seq_along(marks)) {
    hit <- Mb[bi, k]
    P[, k] <- tabulate(gi[hit], nbins = length(genes_gr)) > 0
  }
  P
}

#' UpSet-style figure: exact mark combinations (each gene counted once, in the
#' combination of exactly the marks it carries), with
#'   top:    expression (log10) of the genes in each combination (optional)
#'   middle: number of genes per combination
#'   bottom: dot matrix showing which marks make up each combination
#'   right:  set size = total genes carrying each mark (any combination)
#' P: genes x marks logical matrix. expr_log: numeric per gene (or NULL).
upset_figure <- function(P, expr_log = NULL, n_show = 15, min_n = 1,
                         base_size = 11, title = NULL) {
  marks <- colnames(P)
  m     <- length(marks)
  key   <- apply(P, 1, function(r) paste(marks[r], collapse = "+"))
  key[key == ""] <- NA                      # genes with no mark are not drawn

  tab <- sort(table(key), decreasing = TRUE)
  tab <- tab[tab >= min_n]
  tab <- utils::head(tab, n_show)
  if (length(tab) == 0) return(NULL)

  combos <- names(tab)
  n_c    <- length(combos)
  x_idx  <- seq_len(n_c)
  ypos   <- function(idx) m - idx + 1       # first mark drawn at the top

  bar_df <- data.frame(x = x_idx, n = as.integer(tab), stringsAsFactors = FALSE)

  mem <- do.call(rbind, lapply(combos, function(cb) marks %in% strsplit(cb, "+", fixed = TRUE)[[1]]))
  mat_df <- data.frame(
    x       = rep(x_idx, times = m),
    mark    = rep(marks, each = n_c),
    present = as.vector(mem),
    stringsAsFactors = FALSE
  )
  mat_df$ypos <- ypos(match(mat_df$mark, marks))

  # Lines joining the filled dots of each combination
  seg_df <- do.call(rbind, lapply(seq_len(n_c), function(i) {
    idx <- which(mem[i, ])
    if (length(idx) < 2) return(NULL)
    data.frame(x = i, ymin = ypos(max(idx)), ymax = ypos(min(idx)))
  }))

  set_df <- data.frame(mark = marks, size = as.integer(colSums(P)),
                       ypos = ypos(seq_len(m)), stringsAsFactors = FALSE)

  th <- theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = base_size))

  # ---- expression boxes (top) ----
  p_expr <- NULL
  if (!is.null(expr_log)) {
    ex_df <- do.call(rbind, lapply(seq_len(n_c), function(i) {
      v <- expr_log[!is.na(key) & key == combos[i] & !is.na(expr_log)]
      if (length(v) == 0) return(NULL)
      data.frame(x = i, value = v)
    }))
    if (!is.null(ex_df)) {
      p_expr <- ggplot(ex_df, aes(x = x, y = value, group = x)) +
        geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey92", colour = "grey20") +
        geom_jitter(width = 0.15, height = 0, size = 0.3, alpha = 0.25, colour = "grey30") +
        scale_x_continuous(limits = c(0.5, n_c + 0.5), breaks = NULL, expand = c(0, 0)) +
        labs(x = NULL, y = "log10(count + 1)", title = NULL) + th
    }
  }

  # ---- intersection bars (middle) ----
  p_bar <- ggplot(bar_df, aes(x = x, y = n)) +
    geom_col(fill = "grey25", width = 0.7) +
    geom_text(aes(label = n), vjust = -0.4, size = base_size / 3.5) +
    scale_x_continuous(limits = c(0.5, n_c + 0.5), breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Intersection size (genes)") + th

  # ---- dot matrix (bottom) ----
  stripe_df <- data.frame(ypos = ypos(seq_len(m)))
  stripe_df <- stripe_df[seq_len(m) %% 2 == 1, , drop = FALSE]
  p_mat <- ggplot() +
    geom_rect(data = stripe_df, aes(xmin = 0.5, xmax = n_c + 0.5,
                                    ymin = ypos - 0.5, ymax = ypos + 0.5),
              fill = "grey93", inherit.aes = FALSE) +
    geom_segment(data = seg_df, aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 1.1, colour = "black") +
    geom_point(data = mat_df, aes(x = x, y = ypos, colour = present), size = base_size / 3.2) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), guide = "none") +
    scale_x_continuous(limits = c(0.5, n_c + 0.5), breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(breaks = ypos(seq_len(m)), labels = marks, expand = c(0, 0.5)) +
    labs(x = NULL, y = NULL) + th +
    theme(panel.border = element_blank(), axis.ticks = element_blank())

  # ---- set sizes (right of the dot matrix) ----
  p_set <- ggplot(set_df, aes(x = size, y = ypos)) +
    geom_col(orientation = "y", fill = "grey35", width = 0.6) +
    geom_text(aes(label = format(size, big.mark = ",")), hjust = -0.15, size = base_size / 3.5) +
    scale_y_continuous(limits = c(0.5, m + 0.5), breaks = NULL, expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3)),
                       labels = function(x) format(x, big.mark = ",")) +
    labs(x = "Set size (genes)", y = NULL) + th +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  if (!is.null(p_expr)) {
    out <- patchwork::wrap_plots(p_expr, patchwork::plot_spacer(),
                                 p_bar,  patchwork::plot_spacer(),
                                 p_mat,  p_set,
                                 ncol = 2, widths = c(4, 1.2), heights = c(3, 2, 1.8))
  } else {
    out <- patchwork::wrap_plots(p_bar, patchwork::plot_spacer(),
                                 p_mat, p_set,
                                 ncol = 2, widths = c(4, 1.2), heights = c(2, 1.8))
  }
  if (!is.null(title)) out <- out + patchwork::plot_annotation(title = title)
  out
}


#' Load the sample metadata (colData). Accepts xlsx, csv or tsv. Uses the column
#' named SampleName (or the first column) for sample IDs and Condition (or the
#' second column) for groups. Whitespace around labels is trimmed.
load_coldata <- function(path, file_name = path) {
  ext <- tolower(tools::file_ext(file_name))
  df <- if (ext %in% c("xlsx", "xlsm")) {
    openxlsx::read.xlsx(path, sheet = 1)
  } else {
    first <- readLines(path, n = 1, warn = FALSE)
    sep   <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
    read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
               stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(df) <- trimws(names(df))
  nm <- tolower(names(df))

  s_col <- which(nm %in% c("samplename", "sample", "sample_name", "sampleid", "id"))
  c_col <- which(nm %in% c("condition", "group"))
  if (length(s_col) == 0) s_col <- 1L
  if (length(c_col) == 0) c_col <- if (ncol(df) >= 2) 2L else NA_integer_
  if (is.na(c_col) || c_col == s_col[1]) stop("No Condition column found.")

  out <- data.frame(
    sample    = trimws(as.character(df[[s_col[1]]])),
    condition = trimws(as.character(df[[c_col[1]]])),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$sample) & nzchar(out$sample) &
             !is.na(out$condition) & nzchar(out$condition), , drop = FALSE]
  if (nrow(out) == 0) stop("No rows with sample names and conditions.")
  out[!duplicated(out$sample), , drop = FALSE]
}

#' Samples used to rank genes: the samples of the chosen colData groups that are
#' present in the count table. With no colData, every sample in the table is used.
decile_samples_for <- function(ct, coldata, groups) {
  if (is.null(coldata)) return(ct$samples)
  cd <- coldata[coldata$condition %in% groups, , drop = FALSE]
  intersect(cd$sample, ct$samples)
}

#' Plot the expression-stratified profile: one line per decile, coloured from
#' blue (lowest expression) to red (highest).
#'   display = "separate": one panel per condition (columns), rows = marks.
#'   display = "overlay":  conditions share one panel per mark; each condition
#'                         is drawn with its own line type (solid, dashed, ...).
#' prof needs: slot, mean_signal, expr_group, mark, condition, stage.
plot_decile_profile <- function(prof, n_flank, n_body, up_bp, down_bp,
                                decile_cols, caption = NULL, display = "separate") {
  prof <- prof[!is.na(prof$mean_signal), , drop = FALSE]
  prof$expr_decile <- factor(prof$expr_group, levels = 1:10)
  prof$condition   <- factor(prof$condition, levels = unique(prof$condition))
  multi_stage <- length(unique(prof$stage)) > 1
  prof$row_lab <- if (multi_stage) paste(prof$stage, prof$mark, sep = " | ") else prof$mark
  prof$row_lab <- factor(prof$row_lab, levels = unique(prof$row_lab))

  overlay <- identical(display, "overlay")
  n_cond  <- length(levels(prof$condition))
  prof$grp_id <- interaction(prof$expr_decile, prof$condition, prof$row_lab, drop = TRUE)

  n_total <- n_flank + n_body + n_flank
  tss_x   <- n_flank + 0.5
  tes_x   <- n_flank + n_body + 0.5

  mapping <- if (overlay) {
    aes(x = slot, y = mean_signal, colour = expr_decile,
        linetype = condition, group = grp_id)
  } else {
    aes(x = slot, y = mean_signal, colour = expr_decile, group = grp_id)
  }
  facet_f <- if (overlay) as.formula("row_lab ~ .") else as.formula("row_lab ~ condition")

  p <- ggplot(prof, mapping) +
    geom_vline(xintercept = c(tss_x, tes_x), linetype = "dotted", colour = "grey55") +
    geom_line(linewidth = 0.7) +
    facet_grid(facet_f, scales = "free_y") +
    scale_colour_manual(values = decile_cols, drop = FALSE,
                        name = "Expression decile\n(10 = highest)") +
    guides(colour = guide_legend(reverse = TRUE)) +
    scale_x_continuous(
      limits = c(0.5, n_total + 0.5),
      breaks = c(0.5, tss_x, tes_x, n_total + 0.5),
      labels = c(paste0("-", up_bp / 1000, "kb"), "TSS", "TES",
                 paste0("+", down_bp / 1000, "kb")),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Mean signal (ChromStaR RPKM)",
         title = "Signal across gene body, stratified by gene expression",
         caption = caption) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#34495e"),
      strip.text       = element_text(colour = "white", face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      panel.spacing.x  = grid::unit(1.2, "lines")
    )

  if (overlay) {
    lt <- rep(c("solid", "dashed", "dotdash", "longdash", "twodash"), length.out = n_cond)
    p <- p + scale_linetype_manual(values = lt, name = "Condition")
  }
  p
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

# =============================================================================
# BIVALENT / MARK-COMBINATION HELPERS
# =============================================================================
# chromstaR stores, per bin and per condition, the set of marks called present
# as a string in mcols(bins)$combination.<condition>, e.g. "[H3K27me3+H3K4me3]"
# or "[]" for no mark. Everything below that works on "states" reads those
# strings; the "sum" mode reads the counts.rpkm columns instead.

#' Split chromstaR combination strings into sorted character vectors of marks.
combo_mark_sets <- function(x) {
  x <- gsub("^\\s*\\[|\\]\\s*$", "", as.character(x))
  lapply(strsplit(x, "+", fixed = TRUE), function(v) {
    v <- trimws(v)
    sort(unique(v[nzchar(v)]))
  })
}

#' Display label for a set of marks, in chromstaR's bracket style.
combo_label <- function(marks) paste0("[", paste(sort(marks), collapse = "+"), "]")

#' Pull one condition's combination column off the bins as a plain vector.
get_combination_vector <- function(hmm, cond) {
  mc  <- S4Vectors::mcols(hmm$bins)
  col <- paste0("combination.", cond)
  if (!col %in% colnames(mc)) return(NULL)
  v <- mc[[col]]
  if (methods::is(v, "Rle")) v <- as.vector(v)
  v
}

#' Per-bin TRUE/FALSE for every line the Bivalent tab draws.
#'
#' rule = "contains": the combination line is every bin where ALL chosen marks
#'                    are present (other marks allowed); a single-mark line is
#'                    every bin where that mark is present.
#' rule = "exact":    the combination line is every bin whose state is exactly
#'                    the chosen marks and nothing else; a single-mark line is
#'                    every bin whose state is exactly that one mark.
#' "[]" is always the no-mark state.
#'
#' Works on the factor levels (a few hundred states at most), then maps back to
#' the bins, so it stays fast on 500k+ bins.
bivalent_state_flags <- function(combo, marks_sel, rule = "contains",
                                 include_singles = TRUE, include_none = TRUE) {
  f    <- as.factor(combo)
  sets <- combo_mark_sets(levels(f))
  code <- as.integer(f)

  lvl_flag <- function(test) {
    lf  <- vapply(sets, test, logical(1))
    out <- lf[code]
    out[is.na(out)] <- FALSE
    out
  }

  res <- list()
  res[[combo_label(marks_sel)]] <- lvl_flag(function(s) {
    if (identical(rule, "exact")) setequal(s, marks_sel) else all(marks_sel %in% s)
  })
  if (include_singles) {
    for (m in marks_sel) {
      res[[combo_label(m)]] <- lvl_flag(function(s) {
        if (identical(rule, "exact")) setequal(s, m) else m %in% s
      })
    }
  }
  if (include_none) res[["[]"]] <- lvl_flag(function(s) length(s) == 0)
  res
}

#' Every chromstaR state that contains at least `min_marks` marks, with how many
#' bins carry it, across the given (hmm, conditions) pairs. Used to let the user
#' pick "a bivalent state that actually occurs in my data".
found_combination_states <- function(stage_list, min_marks = 2) {
  counts <- list()
  for (st in stage_list) {
    for (cond in st$conditions) {
      v <- get_combination_vector(st$hmm, cond)
      if (is.null(v)) next
      tb <- table(as.character(v))
      for (nm in names(tb)) {
        prev <- if (is.null(counts[[nm]])) 0 else counts[[nm]]
        counts[[nm]] <- prev + as.numeric(tb[[nm]])
      }
    }
  }
  if (length(counts) == 0) return(NULL)
  lab  <- names(counts)
  sets <- combo_mark_sets(lab)
  df <- data.frame(state  = vapply(sets, combo_label, character(1)),
                   n_marks = lengths(sets),
                   n_bins  = unlist(counts, use.names = FALSE),
                   stringsAsFactors = FALSE)
  df <- stats::aggregate(n_bins ~ state + n_marks, data = df, FUN = sum)
  df <- df[df$n_marks >= min_marks, , drop = FALSE]
  df[order(-df$n_bins), , drop = FALSE]
}

#' Which bins fall at which x-position of the TSS / gene body / TES profile.
#' Same geometry as compute_enrichment_profile() (2 kb-style flanks, gene body
#' scaled 0-1, strand-aware, bins whose midpoint falls outside a zone dropped),
#' returned as a plain table so any per-bin quantity can be averaged on it.
enrichment_slot_hits <- function(bins, genes_gr, upstream = 2000,
                                 downstream = 2000, n_bins = 40) {
  bin_mid_all <- (start(bins) + end(bins)) / 2

  s_g   <- start(genes_gr)
  e_g   <- end(genes_gr)
  str_g <- as.character(strand(genes_gr))
  chr_g <- as.character(seqnames(genes_gr))
  is_minus <- str_g == "-"
  tss   <- ifelse(is_minus, e_g, s_g)
  tes   <- ifelse(is_minus, s_g, e_g)

  up_start   <- ifelse(is_minus, tss,                        pmax(1, tss - upstream))
  up_end     <- ifelse(is_minus, tss + upstream,             tss)
  down_start <- ifelse(is_minus, pmax(1, tes - downstream),  tes)
  down_end   <- ifelse(is_minus, tes,                        tes + downstream)
  gene_idx   <- seq_along(genes_gr)

  zone <- function(starts, ends) {
    ok <- starts >= 1 & ends >= 1 & starts <= ends
    GRanges(seqnames = chr_g[ok], ranges = IRanges(starts[ok], ends[ok]),
            gene_idx = gene_idx[ok])
  }
  zones <- list(upstream   = zone(up_start, up_end),
                body       = zone(s_g, e_g),
                downstream = zone(down_start, down_end))

  out <- list()
  for (zname in names(zones)) {
    zr <- zones[[zname]]
    if (length(zr) == 0) next
    ov <- suppressWarnings(findOverlaps(bins, zr))
    if (length(ov) == 0) next
    bi  <- queryHits(ov)
    zi  <- subjectHits(ov)
    g_i <- mcols(zr)$gene_idx[zi]
    mid <- bin_mid_all[bi]
    minus <- is_minus[g_i]

    if (zname == "body") {
      glen <- e_g[g_i] - s_g[g_i]; glen[glen == 0] <- 1
      frac <- (mid - s_g[g_i]) / glen
    } else {
      frac <- (mid - start(zr)[zi]) / width(zr)[zi]
    }
    frac  <- ifelse(minus, 1 - frac, frac)
    x_val <- switch(zname, upstream = -1 + frac, body = frac, downstream = 1 + frac)

    keep <- !is.na(frac) & frac >= 0 & frac <= 1
    if (!any(keep)) next
    out[[zname]] <- data.frame(bin = bi[keep],
                               slot = round(x_val[keep] * n_bins),
                               gene_idx = g_i[keep])
  }
  if (length(out) == 0) return(NULL)
  dplyr::bind_rows(out)
}

#' Average a per-bin numeric vector over each slot and compare it with its
#' genome-wide mean: log( mean(v in slot) / mean(v genome-wide) ).
#' For a TRUE/FALSE state vector this is log(observed / expected fraction of
#' bins in that state) — the quantity in chromstaR's combination enrichment plot.
slot_log_ratio <- function(v, hits) {
  expected <- mean(v, na.rm = TRUE)
  if (!is.finite(expected) || expected <= 0) return(NULL)
  vv  <- v[hits$bin]
  obs <- tapply(vv, hits$slot, mean, na.rm = TRUE)
  nb  <- tapply(!is.na(vv), hits$slot, sum)
  ng  <- tapply(hits$gene_idx, hits$slot, function(g) length(unique(g)))
  data.frame(slot       = as.integer(names(obs)),
             observed   = as.numeric(obs),
             expected   = expected,
             log_ratio  = log(as.numeric(obs) / expected),
             n_bins_avg = as.integer(nb),
             n_genes    = as.integer(ng),
             stringsAsFactors = FALSE)
}

#' Profile for the Bivalent tab.
#'
#' mode = "cooccur": lines are the chromstaR STATES — the chosen combination,
#'                   each chosen mark on its own, and "[]" — each as
#'                   log(observed/expected fraction of bins in that state).
#' mode = "sum":     lines are RPKM signal — the summed signal of the chosen
#'                   marks and each mark's own signal — each as
#'                   log(mean in slot / genome-wide mean).
#'   sum_norm = "normalised": each mark is divided by its own genome-wide mean
#'                   before summing, so every mark carries equal weight.
#'   sum_norm = "raw": plain RPKM sum; the mark with the most reads dominates.
compute_bivalent_profile <- function(hmm, genes_gr, marks_sel, conditions,
                                     mode = "cooccur", rule = "contains",
                                     sum_norm = "normalised", show_none = TRUE,
                                     upstream = 2000, downstream = 2000,
                                     n_bins = 40) {
  hits <- enrichment_slot_hits(hmm$bins, genes_gr, upstream, downstream, n_bins)
  if (is.null(hits) || nrow(hits) == 0) return(NULL)

  comb_name <- combo_label(marks_sel)
  sum_name  <- paste0("Sum: ", paste(marks_sel, collapse = " + "))
  rows <- list()

  if (identical(mode, "sum")) {
    bins_df  <- expand_bins_df(hmm)
    rpkm_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
    for (cond in conditions) {
      sigs <- list()
      for (m in marks_sel) {
        cols <- grep(paste0("counts.rpkm.", m, ".", cond), rpkm_all,
                     value = TRUE, fixed = TRUE)
        if (length(cols) == 0) next
        sigs[[m]] <- rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)
      }
      if (length(sigs) < length(marks_sel)) next   # a chosen mark is missing here
      parts <- if (identical(sum_norm, "raw")) sigs else
        lapply(sigs, function(s) { mu <- mean(s, na.rm = TRUE); if (mu > 0) s / mu else s })
      # single marks go through the same normalisation as the sum, so on the
      # "observed" axis all lines share one unit (log ratios are unaffected)
      lines <- c(stats::setNames(list(Reduce(`+`, parts)), sum_name), parts)
      for (ln in names(lines)) {
        r <- slot_log_ratio(lines[[ln]], hits)
        if (is.null(r)) next
        r$line <- ln; r$line_kind <- if (ln == sum_name) "combination" else "single"
        r$condition <- cond
        rows[[length(rows) + 1L]] <- r
      }
    }
  } else {
    for (cond in conditions) {
      combo <- get_combination_vector(hmm, cond)
      if (is.null(combo)) next
      flags <- bivalent_state_flags(combo, marks_sel, rule, include_none = show_none)
      for (ln in names(flags)) {
        r <- slot_log_ratio(as.numeric(flags[[ln]]), hits)
        if (is.null(r)) next          # state never occurs in this condition
        r$line <- ln
        r$line_kind <- if (ln == comb_name) "combination" else if (ln == "[]") "none" else "single"
        r$condition <- cond
        rows[[length(rows) + 1L]] <- r
      }
    }
  }

  if (length(rows) == 0) return(NULL)
  out <- dplyr::bind_rows(rows)
  out$n_genes_total <- length(unique(hits$gene_idx))
  first <- if (identical(mode, "sum")) sum_name else comb_name
  lv <- c(first, setdiff(unique(out$line[out$line_kind == "single"]), first),
          if (any(out$line_kind == "none")) "[]")
  out$line <- factor(out$line, levels = unique(lv))
  # log(0) for a slot where the state never occurs: keep the row (its observed
  # value of 0 is real and is drawn on the "% of bins" axis), blank the log.
  out$log_ratio[!is.finite(out$log_ratio)] <- NA_real_
  out
}

#' Which genes carry the chosen combination at their promoter, per condition.
#'
#' gene_rule = "any":  a gene is called if ANY bin in its promoter window is in
#'                     the combination state.
#' gene_rule = "half": a gene is called if AT LEAST HALF of its promoter bins are.
#'
#' Returns one row per gene with, per condition, the fraction of promoter bins
#' in the state and the TRUE/FALSE call, plus gained / lost / kept / absent
#' between consecutive conditions in the order given. With sum_norm supplied,
#' also the mean summed signal over the promoter (same normalisation as the plot).
compute_bivalent_genes <- function(hmm, genes_gr, marks_sel, conditions,
                                   rule = "contains", gene_rule = "any",
                                   prom_up = 1000, prom_down = 1000,
                                   sum_norm = NULL) {
  bins <- hmm$bins
  prom <- suppressWarnings(promoters(genes_gr, upstream = prom_up, downstream = prom_down))
  ov   <- suppressWarnings(findOverlaps(prom, bins, ignore.strand = TRUE))
  gi   <- queryHits(ov)
  bi   <- subjectHits(ov)

  g_names <- mcols(genes_gr)$gene_name
  if (is.null(g_names)) g_names <- as.character(seq_along(genes_gr))
  n_g <- length(genes_gr)

  out <- data.frame(
    gene   = g_names,
    chr    = as.character(seqnames(genes_gr)),
    start  = start(genes_gr),
    end    = end(genes_gr),
    strand = as.character(strand(genes_gr)),
    stringsAsFactors = FALSE)
  out$locus <- paste0(out$chr, ":", out$start, "-", out$end)
  out$n_promoter_bins <- tabulate(gi, nbins = n_g)

  per_gene_mean <- function(v) {
    m <- rep(NA_real_, n_g)
    if (length(gi) == 0) return(m)
    agg <- tapply(v[bi], gi, mean, na.rm = TRUE)
    m[as.integer(names(agg))] <- as.numeric(agg)
    m
  }

  bins_df <- NULL; rpkm_all <- NULL
  if (!is.null(sum_norm)) {
    bins_df  <- expand_bins_df(hmm)
    rpkm_all <- grep("counts.rpkm", colnames(bins_df), value = TRUE, fixed = TRUE)
  }

  used <- character(0)
  for (cond in conditions) {
    combo <- get_combination_vector(hmm, cond)
    if (is.null(combo)) next
    f <- bivalent_state_flags(combo, marks_sel, rule,
                              include_singles = FALSE, include_none = FALSE)[[1]]
    frac <- per_gene_mean(as.numeric(f))
    call <- if (identical(gene_rule, "half")) !is.na(frac) & frac >= 0.5
            else                               !is.na(frac) & frac > 0
    out[[paste0("frac_", cond)]]      <- round(frac, 3)
    out[[paste0("bivalent_", cond)]]  <- call

    if (!is.null(sum_norm)) {
      sigs <- list()
      for (m in marks_sel) {
        cols <- grep(paste0("counts.rpkm.", m, ".", cond), rpkm_all, value = TRUE, fixed = TRUE)
        if (length(cols) == 0) next
        s <- rowMeans(bins_df[, cols, drop = FALSE], na.rm = TRUE)
        if (!identical(sum_norm, "raw")) { mu <- mean(s, na.rm = TRUE); if (mu > 0) s <- s / mu }
        sigs[[m]] <- s
      }
      if (length(sigs) == length(marks_sel)) {
        out[[paste0("sum_signal_", cond)]] <- round(per_gene_mean(Reduce(`+`, sigs)), 4)
      }
    }
    used <- c(used, cond)
  }
  if (length(used) == 0) return(NULL)

  calls <- as.matrix(out[, paste0("bivalent_", used), drop = FALSE])
  out$n_conditions_bivalent <- rowSums(calls)
  out$bivalent_in <- apply(calls, 1, function(r) paste(used[r], collapse = ", "))

  if (length(used) >= 2) {
    for (k in seq_len(length(used) - 1)) {
      a <- calls[, k]; b <- calls[, k + 1]
      out[[paste0(used[k], " -> ", used[k + 1])]] <-
        ifelse(a & !b, "lost", ifelse(!a & b, "gained", ifelse(a & b, "kept", "absent")))
    }
  }
  attr(out, "conditions") <- used
  out
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
#'
#' pool_labels: optional background list (e.g. expression-matched genes). When
#' given, every random draw comes only from these genes instead of the whole
#' annotation. exclude_A_from_pool drops set A's own genes from that pool.
compute_gene_set_permutation <- function(all_gene_scores, post_cols, genes_A_labels,
                                          n_perm = 1000, seed = NULL,
                                          summary_stat = "mean",
                                          pool_labels = NULL,
                                          exclude_A_from_pool = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  use_median <- identical(summary_stat, "median")

  gene_scores <- all_gene_scores
  n_total <- nrow(gene_scores)
  n_A     <- sum(gene_scores$gene %in% genes_A_labels)

  if (n_A == 0) {
    return(list(error = "None of gene set A's genes were found in the precomputed universe — check gene name matching."))
  }
  is_A <- gene_scores$gene %in% genes_A_labels

  # ---- Which genes random draws may come from --------------------------------
  pool_overlap_A <- 0L
  if (is.null(pool_labels)) {
    if (n_A >= n_total) {
      return(list(error = "Gene set A is the same size as (or larger than) the full gene universe — no room to draw a smaller random comparison set."))
    }
    pool_idx <- seq_len(n_total)
  } else {
    in_pool <- gene_scores$gene %in% pool_labels
    pool_overlap_A <- sum(in_pool & is_A)
    if (exclude_A_from_pool) in_pool <- in_pool & !is_A
    pool_idx <- which(in_pool)
    if (length(pool_idx) == 0) {
      return(list(error = "None of the background list (B) genes have promoter scores (after removing set A's genes) — check gene name matching."))
    }
    if (length(pool_idx) <= n_A) {
      return(list(error = paste0("The background list (B) has ", length(pool_idx),
                                 " usable genes but gene set A has ", n_A,
                                 ". Each random draw needs as many genes as A, so B must be larger than A (ideally several times larger).")))
    }
  }
  n_pool <- length(pool_idx)
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
    idx <- pool_idx[sample.int(n_pool, n_A)]
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
    n_random_pool   = n_pool,
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
  rep_idx <- pool_idx[sample.int(n_pool, n_A)]
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
       long_scores = long_scores, n_A = n_A, n_total = n_pool,
       n_pool_overlap_A = pool_overlap_A)
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
source("rf_shap_tools.R")
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
      menuItem("Expression Deciles", tabName = "deciles", icon = icon("th-large")),
      menuItem("Enrichment Profile", tabName = "enrichment", icon = icon("chart-line")),
      menuItem("Bivalent",        tabName = "bivalent",   icon = icon("object-group")),
      menuItem("Region Browser",  tabName = "browser",    icon = icon("dna")),
      menuItem("Differential Peaks", tabName = "diffpeaks", icon = icon("chart-simple")),
      menuItem("Gene Set Comparison", tabName = "genesetcompare", icon = icon("dna")),
      menuItem("Mark Combinations", tabName = "upset", icon = icon("chart-bar")),
      menuItem("RF Region Model",    tabName = "rs_model", icon = icon("tree")),
      menuItem("SHAP Mark x Region", tabName = "rs_shap",  icon = icon("chart-bar")),
      menuItem("SHAP Gene Explorer", tabName = "rs_genes", icon = icon("dna")),
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
      tabItem(tabName = "rs_model",
              fluidRow(
                box(title = "RF + SHAP input", width = 4, status = "warning", solidHeader = TRUE,
                    tags$small(style = "color:#888",
                               "Uses the objects already loaded on the Load Data tab (ChromStaR objects, count tables and colData). Pick two conditions from any stage."),
                    br(), br(),
                    selectInput("rs_pick1", "Condition 1", choices = NULL),
                    selectInput("rs_pick2", "Condition 2", choices = NULL),
                    numericInput("rs_ntree", "Trees (combined model)", 500, 100, 2000, 50),
                    numericInput("rs_ntree_single", "Trees (single region)", 300, 100, 2000, 50),
                    numericInput("rs_folds", "CV folds", 5, 2, 10),
                    numericInput("rs_ncores", "Cores for SHAP", 4, 1, 16),
                    numericInput("rs_n_shap", "Test genes to explain (0 = all)", 0, 0),
                    actionButton("rs_run", "Run RF + SHAP", class = "btn-primary btn-block"),
                    br(),
                    downloadButton("rs_zip", "Download all outputs (zip)")),
                box(title = "Model performance", width = 8, status = "info", solidHeader = TRUE,
                    tableOutput("rs_perf_tab"), tableOutput("rs_cv_tab"),
                    plotOutput("rs_single_plot", height = "380px"))),
              fluidRow(box(title = "Feature correlation (training set)", width = 12,
                           plotOutput("rs_corr_plot", height = "750px")))
      ),
      
      tabItem(tabName = "rs_shap",
              fluidRow(box(title = "Which mark matters where", width = 12, status = "primary",
                           solidHeader = TRUE, plotOutput("rs_bubble", height = "420px"))),
              fluidRow(box(title = "Expressed vs silent genes", width = 12,
                           plotOutput("rs_heat", height = "620px"))),
              fluidRow(box(title = "Top features (beeswarm)", width = 12,
                           uiOutput("rs_beeswarm_ui")))
      ),
      
      tabItem(tabName = "rs_genes",
              fluidRow(box(title = "Genes", width = 4, status = "primary", solidHeader = TRUE,
                           selectInput("rs_genes_stage", "Model", choices = NULL),
                           textInput("rs_genes_txt", "Gene IDs (first matched one is plotted)", ""))),
              fluidRow(box(title = "Waterfall", width = 12,
                           plotOutput("rs_waterfall", height = "480px")))
      ),
      
      tabItem(tabName = "load",
        fluidRow(
          column(12,
            tags$div(style = "text-align:center; padding: 10px 0 20px 0;",
              tags$img(src = "ihpe_logo.jpg", height = "70px"),
              tags$p(style = "color:#888; font-size:12px; margin-top:8px;",
                HTML("&copy; 2026 Janan Gawra &mdash; IHPE (Interactions Hôtes-Pathogènes-Environnements, UMR 5244 CNRS-IFREMER-UM-UPVD)"))
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
          box(title = "Sample Metadata (colData) — all samples", width = 12, status = "primary",
            tags$small(style = "color:#888",
              HTML("Sample sheet with columns <b>SampleName</b>, <b>Condition</b> and optional <b>Replicate</b> (xlsx, csv or tsv). Used by the Expression Deciles tab to choose which groups of samples rank the genes.")),
            br(), br(),
            fileInput("coldata_file", "Upload colData (.xlsx / .csv / .tsv)",
                      accept = c(".xlsx", ".xlsm", ".csv", ".tsv", ".txt")),
            verbatimTextOutput("coldata_summary")
          )
        ),
        fluidRow(
          box(title = "Normalised Count Table — Stage A (optional)", width = 6, status = "primary",
            tags$small(style = "color:#888",
              HTML("First column = gene IDs (header can be blank), then one numeric column per sample. Use DESeq2 normalised counts (<code>counts(dds, normalized = TRUE)</code>).")),
            br(), br(),
            fileInput("counts_file", "Upload count table (.csv / .tsv / .txt)",
                      accept = c(".csv", ".tsv", ".txt")),
            verbatimTextOutput("counts_summary")
          ),
          box(title = "Normalised Count Table — Stage B (optional)", width = 6, status = "primary",
            tags$small(style = "color:#888",
              "Same format as Stage A. Used when Stage B is selected in the Expression Deciles tab."),
            br(), br(),
            fileInput("counts_file_b", "Upload count table (.csv / .tsv / .txt)",
                      accept = c(".csv", ".tsv", ".txt")),
            verbatimTextOutput("counts_summary_b")
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
      # -----------------------------------------------------------------------
      # TAB 2.2 — EXPRESSION DECILES (gene body profile per expression decile)
      # -----------------------------------------------------------------------
      # TAB 2.2 — EXPRESSION DECILES (gene body profile per expression decile)
      # -----------------------------------------------------------------------
      tabItem(tabName = "deciles",
        fluidRow(
          box(title = "Expression Decile Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Genes are ranked by mean normalised expression and split into 10 equal groups (deciles). Decile 1 = lowest, 10 = highest."),
            hr(),
            h5("Stages"),
            uiOutput("decile_stage_selector"),
            h5("Rank genes using these groups"),
            tags$small(style = "color:#888",
              "Groups come from the colData sheet. Only samples present in the count table are used."),
            uiOutput("decile_group_selector"),
            numericInput("decile_min_expr", "Exclude genes with mean count below",
                         value = 1, min = 0, step = 0.5),
            hr(),
            h5("Marks to display"),
            uiOutput("decile_mark_selector"),
            h5("Conditions to show"),
            uiOutput("decile_condition_selector"),
            radioButtons("decile_display", "Display",
                         choices = c("Separate panels (one per condition)" = "separate",
                                     "Overlay (conditions in one panel: solid / dashed)" = "overlay"),
                         selected = "separate"),
            hr(),
            h5("Window (scaled gene body)"),
            numericInput("decile_upstream",   "Upstream flank (bp)",   value = 1500, min = 0, step = 100),
            numericInput("decile_downstream", "Downstream flank (bp)", value = 1500, min = 0, step = 100),
            numericInput("decile_n_body",  "Gene body bins", value = 20, min = 10, max = 200, step = 5),
            numericInput("decile_n_flank", "Bins per flank", value = 5, min = 2, max = 50, step = 1),
            actionButton("run_decile", "Compute Profile", class = "btn-primary btn-block"),
            hr(),
            box(
              title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance & Size"),
              width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
              h5("Figure size on screen"),
              sliderInput("dec_plot_height", "Plot height (px)", min = 300, max = 1400, value = 720, step = 20),
              sliderInput("dec_plot_width", "Plot width (% of box)", min = 40, max = 100, value = 100, step = 5),
              h5("Text sizes"),
              sliderInput("dec_plot_title_size", "Title size", min = 10, max = 24, value = 14, step = 1),
              sliderInput("dec_plot_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
              sliderInput("dec_plot_axis_text_size", "Axis labels size", min = 6, max = 18, value = 10, step = 1),
              sliderInput("dec_plot_legend_size", "Legend text size", min = 8, max = 16, value = 10, step = 1),
              h5("Lines"),
              sliderInput("dec_plot_line_width", "Line width", min = 0.3, max = 3, value = 0.7, step = 0.1),
              checkboxInput("dec_plot_show_grid", "Show gridlines", value = TRUE),
              h5("Download size"),
              numericInput("dec_dl_width",  "Width (inches)",  value = 12, min = 4, max = 30, step = 0.5),
              numericInput("dec_dl_height", "Height (inches)", value = 8, min = 3, max = 30, step = 0.5),
              tags$small(style = "color:#888; display:block;",
                "Appearance changes apply immediately. Click 'Compute Profile' only after changing the data settings.")
            )
          ),
          box(title = "Signal Across Gene Body, by Expression Decile", width = 9, status = "primary",
            uiOutput("decile_plot_ui"),
            downloadButton("dl_decile_plot", "Download Plot (PNG)"),
            verbatimTextOutput("decile_summary")
          )
        )
      ),

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
      # TAB 2.7 — BIVALENT / MARK COMBINATIONS
      # -----------------------------------------------------------------------
      tabItem(tabName = "bivalent",
        fluidRow(
          box(title = "Bivalent Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Treat two or more marks as one combination (e.g. H3K4me3 + H3K27me3, the classic bivalent promoter) and compare it with each mark on its own, around TSS, gene body and TES."),
            hr(),
            h5("1. Which combination?"),
            radioButtons("bv_define", NULL,
                         choices = c("Pick the marks myself"                   = "pick",
                                     "Choose a combination found in my data"   = "found"),
                         selected = "pick"),
            conditionalPanel(
              condition = "input.bv_define == 'pick'",
              uiOutput("bv_mark_ui"),
              tags$small(style = "color:#888", "Tick 2 or more marks.")
            ),
            conditionalPanel(
              condition = "input.bv_define == 'found'",
              uiOutput("bv_found_ui"),
              tags$small(style = "color:#888",
                "Every chromstaR state with 2 or more marks that occurs in the loaded object(s), most common first, with the number of bins carrying it.")
            ),
            hr(),
            h5("2. How to combine them?"),
            radioButtons("bv_mode", NULL,
                         choices = c("Co-occurrence: regions where the marks are called together" = "cooccur",
                                     "Sum: add the marks' signal together"                        = "sum"),
                         selected = "cooccur"),
            conditionalPanel(
              condition = "input.bv_mode == 'sum'",
              radioButtons("bv_sum_norm", "Before summing",
                           choices = c("Normalise each mark to its genome-wide mean (equal weight)" = "normalised",
                                       "Raw RPKM (stronger mark dominates)"                          = "raw"),
                           selected = "normalised"),
              tags$small(style = "color:#888",
                "A sum does NOT show that the marks sit on the same region: a high sum can come from one mark alone. Use co-occurrence to find true bivalent regions.")
            ),
            conditionalPanel(
              condition = "input.bv_mode == 'cooccur'",
              checkboxInput("bv_show_none", "Also draw the no-mark state \"[]\"", value = TRUE)
            ),
            radioButtons("bv_rule", "What counts as \"the combination\"?",
                         choices = c("Contains all chosen marks (others allowed)" = "contains",
                                     "Exactly the chosen marks, nothing else"      = "exact"),
                         selected = "contains"),
            tags$small(style = "color:#888",
              "Applies to co-occurrence lines and to the bivalent gene calls. With \"exactly\", each single-mark line is also that mark alone."),
            hr(),
            h5("3. Calling bivalent genes"),
            radioButtons("bv_gene_rule", "A gene is bivalent when…",
                         choices = c("Any bin in its promoter is in the state"        = "any",
                                     "At least half of its promoter bins are"         = "half"),
                         selected = "any"),
            numericInput("bv_prom_up",   "Promoter: upstream of TSS (bp)",   value = 1000, min = 0, max = 10000, step = 100),
            numericInput("bv_prom_down", "Promoter: downstream of TSS (bp)", value = 1000, min = 0, max = 10000, step = 100),
            hr(),
            h5("4. Profile window"),
            numericInput("bv_upstream",   "Upstream of TSS (bp)",   value = 2000, min = 100, max = 10000, step = 100),
            numericInput("bv_downstream", "Downstream of TES (bp)", value = 2000, min = 100, max = 10000, step = 100),
            numericInput("bv_n_bins",     "Bins per region",        value = 40,   min = 10,  max = 200,  step = 5),
            hr(),
            h5("Gene scope"),
            radioButtons("bv_gene_scope", NULL,
                         choices = c("All genes in annotation" = "all",
                                     "Only selected genes"     = "selected"),
                         selected = "all"),
            conditionalPanel(
              condition = "input.bv_gene_scope == 'selected'",
              fluidRow(
                column(6, actionButton("select_all_genes_bv", "Select all", class = "btn-sm btn-block")),
                column(6, actionButton("clear_all_genes_bv",  "Clear all",  class = "btn-sm btn-block"))
              ),
              br(),
              uiOutput("gene_selector_bv")
            ),
            hr(),
            h5("Stages to compare"),
            uiOutput("stage_selector_bv"),
            hr(),
            h5("Panel layout"),
            radioButtons("bv_layout", NULL,
                         choices = c("One panel per condition (colour = combination / mark)" = "by_condition",
                                     "One panel per combination / mark (colour = condition)" = "by_line"),
                         selected = "by_condition"),
            checkboxInput("bv_smooth", "Smooth curve (LOESS)", value = FALSE),
            radioButtons("bv_yaxis", "Y axis",
                         choices = c("log(observed/expected)"                        = "log_ratio",
                                     "Observed: % of bins in the state (or mean signal)" = "observed"),
                         selected = "log_ratio"),
            tags$small(style = "color:#888",
              "log(obs/exp) compares each position with how common the state is genome-wide in that condition, so a state that becomes rarer everywhere can still look enriched. The observed axis shows the plain level, where a loss between conditions is seen directly."),

            # ========== PLOT CUSTOMIZATION PANEL (BIVALENT) ==========
            box(
              title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
              width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
              h5("Text Sizes"),
              sliderInput("bv_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
              sliderInput("bv_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
              sliderInput("bv_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
              sliderInput("bv_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
              hr(),
              h5("Lines & Visual Elements"),
              sliderInput("bv_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
              sliderInput("bv_plot_point_size", "Point/marker size", min = 1, max = 5, value = 2, step = 0.5),
              sliderInput("bv_plot_alpha", "Transparency", min = 0.1, max = 1, value = 0.2, step = 0.1),
              hr(),
              h5("Other Elements"),
              checkboxInput("bv_plot_show_grid", "Show gridlines", value = TRUE),
              tags$small(style = "color:#888; display:block; margin-top:10px;",
                "Changes apply immediately to the plot and the PDF."
              )
            ),
            actionButton("run_bivalent", "Compute Bivalent Analysis",
                         class = "btn-primary btn-block")
          ),
          column(9,
            box(title = "Combination vs single marks", width = 12, status = "primary",
              plotOutput("bivalent_plot", height = "600px"),
              downloadButton("dl_bivalent_plot", "Download Plot (PDF)"),
              downloadButton("dl_bivalent_xlsx", "Download Data (Excel)")
            ),
            box(title = "Which genes carry the combination?", width = 12, status = "primary",
              tags$small(style = "color:#888",
                "Genes whose promoter carries the chosen combination in each condition, using chromstaR's state calls and the rule you set in step 3. The arrow columns compare consecutive conditions in your life-cycle order: gained, lost, kept or absent."),
              br(), br(),
              uiOutput("bivalent_summary"),
              br(),
              radioButtons("bv_table_filter", "Show",
                           choices = c("Genes with the combination in at least one condition" = "any",
                                       "Only genes that gain or lose it"                     = "changing",
                                       "All genes"                                            = "all"),
                           selected = "any", inline = TRUE),
              DTOutput("bivalent_gene_table")
            )
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
      # -----------------------------------------------------------------------
      # TAB — MARK COMBINATIONS (UpSet-style figure)
      # -----------------------------------------------------------------------
      tabItem(tabName = "upset",
        fluidRow(
          box(title = "Mark Combinations Settings", width = 3, status = "warning",
            tags$small(style = "color:#888",
              "Counts genes by the exact set of marks they carry (UpSet style). Each gene is counted once, in the combination of exactly its marks."),
            hr(),
            radioButtons("upset_stage", "Stage", choices = c("Stage A" = "A", "Stage B" = "B"),
                         selected = "A"),
            uiOutput("upset_cond_selector"),
            h5("Marks"),
            uiOutput("upset_mark_selector"),
            radioButtons("upset_window", "Gene window",
                         choices = c("Promoter (TSS +/- flank)" = "promoter",
                                     "Whole gene body" = "body"),
                         selected = "promoter"),
            numericInput("upset_up",   "Upstream flank (bp)",   value = 1000, min = 0, step = 100),
            numericInput("upset_down", "Downstream flank (bp)", value = 1000, min = 0, step = 100),
            hr(),
            checkboxInput("upset_show_expr", "Show expression boxes (top)", value = TRUE),
            uiOutput("upset_expr_groups"),
            hr(),
            numericInput("upset_n_show", "Combinations to show", value = 15, min = 1, max = 60, step = 1),
            numericInput("upset_min_n", "Minimum genes per combination", value = 1, min = 1, step = 1),
            actionButton("run_upset", "Compute Combinations", class = "btn-primary btn-block"),
            hr(),
            box(
              title = HTML("<i class='fa fa-sliders-h'></i> Figure Size & Text"),
              width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
              sliderInput("upset_plot_height", "Plot height (px)", min = 400, max = 1600, value = 900, step = 20),
              sliderInput("upset_plot_width", "Plot width (% of box)", min = 40, max = 100, value = 100, step = 5),
              sliderInput("upset_base_size", "Text size", min = 7, max = 20, value = 11, step = 1),
              h5("Download size"),
              numericInput("upset_dl_width",  "Width (inches)",  value = 10, min = 4, max = 30, step = 0.5),
              numericInput("upset_dl_height", "Height (inches)", value = 9,  min = 3, max = 30, step = 0.5)
            )
          ),
          box(title = "Mark Combinations (UpSet)", width = 9, status = "primary",
            uiOutput("upset_plot_ui"),
            downloadButton("dl_upset_plot", "Download Plot (PNG)"),
            verbatimTextOutput("upset_summary")
          )
        )
      ),

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
            h5("Genomic Zones (optional filter)"),
            uiOutput("diff_zone_filter_ui"),
            tags$small(style = "color:#888",
              "Leave empty to show all zones. Select specific zones to filter the results."),
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
          box(title = "Differential regions by genomic zone", width = 12,
              status = "primary", collapsible = TRUE,
            tags$small(style = "color:#888",
              "Breakdown of the differential regions by genomic zone (TSS, gene body, upstream, TES, intergenic). Shows how many regions in each zone changed between conditions."),
            br(), br(),
            DTOutput("diffpeaks_zone_table"),
            br()
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
                        choices = c("Permutation: random genes from the whole annotation"               = "permutation",
                                    "Permutation: random genes from my own background list (list B)"     = "pool",
                                    "Direct test: gene set A vs gene set B (Wilcoxon, no permutation)"   = "manual"),
                        selected = "permutation"),
            tags$small(style = "color:#888",
              "Use the second option when B is a matched background, e.g. all genes with expression similar to A: each random draw then comes only from B, so expression level is controlled for."),
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
              condition = "input.gsc_mode == 'manual' || input.gsc_mode == 'pool'",
              h5(textOutput("gsc_b_title", inline = TRUE)),
              tags$textarea(
                id = "gsc_genes_b", rows = 8,
                style = "width:100%; font-family: monospace; font-size: 12px; resize: vertical;",
                placeholder = "Paste gene IDs here, one per line\n(e.g. your down-regulated list, or all genes\nwith expression similar to set A)"
              ),
              uiOutput("gsc_match_summary_b")
            ),
            conditionalPanel(
              condition = "input.gsc_mode == 'pool'",
              checkboxInput("gsc_pool_exclude_a", "Remove gene set A's genes from the background list before drawing", value = TRUE),
              tags$small(style = "color:#888",
                "Recommended: if A's own genes stay in B, some random draws partly re-sample A and the test becomes conservative.")
            ),
            conditionalPanel(
              condition = "input.gsc_mode == 'permutation' || input.gsc_mode == 'pool'",
              numericInput("gsc_n_perm", "Number of random draws", value = 1000, min = 100, max = 10000, step = 100),
              uiOutput("gsc_perm_help")
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
                             tags$p("ChromstaR Viewer is a point-and-click app for exploring chromatin state and histone modification signal from a ChromstaR ", tags$code("combinedMultiHMM"), " object, alongside a gene annotation and, optionally, expression data. It lets you compare experimental conditions (for example PNA vs PA) across histone marks, browse genomic regions, compare gene sets, and model how the marks at different positions along a gene relate to expression — all without writing R code. Use the sidebar to move between tabs. Each tab is independent, so you can work in any order once data is loaded."),
                             tags$p(tags$b("Tip:"), " no data of your own? A synthetic test set is provided with the app. Generate it with ", tags$code("Rscript make_test_data.R"), " and load the files from the ", tags$code("test_data/"), " folder."),
                             
                             hr(),
                             tags$h3("1. Load Data tab"),
                             tags$p(tags$b("Purpose:"), " upload your data. Nothing else works until Stage A is loaded."),
                             tags$ul(
                               tags$li(tags$b("\"Load ChromstaR Object — Stage A\""), " — browse for your ChromstaR ", tags$code("combinedMultiHMM"), " object. Accepts ", tags$code(".RData"), " / ", tags$code(".rda"), " or a pre-slimmed ", tags$code(".rds"), " file (recommended for speed). Give it a short label (for example \"Protoscolex\") — the label appears in every plot and table. A summary appears below once it has loaded."),
                               tags$li(tags$b("\"Load ChromstaR Object — Stage B\""), " (optional) — a second, independent object, for example a different developmental stage or physiological state. It does not need to share bin coordinates with Stage A. Once loaded, every analysis tab gains a \"Stages to compare\" control."),
                               tags$li(tags$b("\"Upload genes.tsv file\""), " — the gene annotation, shared by both stages. It must be a delimited table with a header line, not a raw GTF. Required columns: ", tags$code("chr"), ", ", tags$code("start"), ", ", tags$code("end"), ". Strongly recommended: ", tags$code("strand"), " (without it minus-strand genes are read backwards and TSS-anchored profiles are wrong), plus ", tags$code("gene_name"), " and ", tags$code("gene_id"), ". Chromosome names must match the ChromstaR object."),
                               tags$li(tags$b("\"Sample Metadata (colData)\""), " — a sample sheet with columns ", tags$code("SampleName"), " and ", tags$code("Condition"), " (xlsx, csv or tsv). Used by the expression-based tabs (Expression Deciles, Mark Combinations and RF Region Model) to choose which samples define expression. Condition names must match those in the ChromstaR object."),
                               tags$li(tags$b("\"Normalised Count Table\""), " — one count table per stage, genes × samples, first column gene IDs. DESeq2 normalised counts are preferred. The app warns if the values look like raw integer counts, which are not comparable between samples."),
                               tags$li(tags$b("\"Detected Marks & Conditions\""), " — lists every mark, condition and chromosome found in each object. Check this after every upload; a missing mark or condition is almost always a naming mismatch."),
                               tags$li(tags$b("\"Condition / Life-cycle Order\""), " — conditions are detected alphabetically, which rarely matches the biology. Drag them into order (for example ", tags$code("TwoDO → Meta → PNA → PA"), "). That order controls legend order, facet order, exported column order and colour assignment across the whole app. Removing a condition with its × hides it from every analysis. \"Reset to detected order\" restores everything."),
                               tags$li(tags$b("\"Condition colours\""), " — one colour per condition, used for that condition in every figure, so a condition looks the same everywhere. \"Reset colours\" restores the defaults.")
                             ),
                             
                             hr(),
                             tags$h3("2. Metagene Profile tab"),
                             tags$p(tags$b("Purpose:"), " average signal (RPKM) around a reference point — TSS, TES or across the gene body — for a chosen set of genes, one panel per mark."),
                             tags$ul(
                               tags$li(tags$b("\"Reference\""), " — tick TSS, TES and/or Gene body. Several can be ticked; each gets its own panel."),
                               tags$li(tags$b("\"Upstream (bp)\" / \"Downstream (bp)\""), " — shown when TSS or TES is ticked. A summary line shows the exact window produced."),
                               tags$li(tags$b("\"Number of bins\""), " — how many segments the window is divided into. More bins give finer detail and a noisier line."),
                               tags$li(tags$b("\"Gene scope\""), " — \"Only selected genes\" uses the list you paste; \"All genes in GTF\" uses every gene in the annotation."),
                               tags$li(tags$b("\"Select all\" / \"Clear all\""), " — fill or empty the gene box with every gene in the annotation."),
                               tags$li(tags$b("Gene box"), " — paste gene IDs, one per line. A count below shows how many were recognised."),
                               tags$li(tags$b("\"Marks to display\""), " — which marks to include."),
                               tags$li(tags$b("\"Condition comparison\""), " — \"Side-by-side\" gives each condition its own panel; \"Overlay\" draws them on shared axes."),
                               tags$li(tags$b("\"Smooth curve (LOESS)\""), " — a smoothed line instead of the raw one. Recommended with many genes."),
                               tags$li(tags$b("\"Compute Profile\""), " — draws the plot. Nothing changes until you press it, including when you move the appearance sliders."),
                               tags$li(tags$b("\"Download Plot\" / \"Download Data (Excel)\""), " — the figure as PDF, or the per-bin numbers behind it as Excel.")
                             ),
                             
                             hr(),
                             tags$h3("2.5 Gene Contributions tab"),
                             tags$p(tags$b("Purpose:"), " break a metagene profile back into the genes behind it. An average cannot tell you whether most genes look like the average or whether a handful carry the whole signal — this tab answers that."),
                             tags$ul(
                               tags$li(tags$b("Stage, mark, condition and reference"), " — one of each at a time, with the same windows and strand handling as the Metagene tab."),
                               tags$li(tags$b("Per-gene heatmap"), " — one row per gene, strongest at the top. Colours can be log-scaled and are capped at the 99th percentile so one extreme gene does not wash out the rest."),
                               tags$li(tags$b("Concentration curve"), " — genes ranked strongest first against their cumulative share of the signal. A curve hugging the top-left means a few genes carry the profile; the dashed diagonal is perfectly even contribution."),
                               tags$li(tags$b("Gini coefficient"), " — 0 means every gene contributes equally; 1 means one gene carries everything. Above about 0.6 the mean describes a minority of genes, so report it alongside the heatmap or a median line."),
                               tags$li(tags$b("Gene ranking table"), " — mean and peak signal, peak position, share of total and running cumulative share for every gene. Sortable, filterable, exportable."),
                               tags$li(tags$b("\"Download Gene Table (Excel)\""), " — two sheets: the ranking, and the gene-by-position matrix behind the heatmap.")
                             ),
                             
                             hr(),
                             tags$h3("3. Enrichment Profile tab"),
                             tags$p(tags$b("Purpose:"), " log(observed/expected) enrichment of each mark around gene boundaries — shows which marks are relatively enriched or depleted across a gene compared with the genome-wide average."),
                             tags$ul(
                               tags$li(tags$b("\"Panel layout\""), " — \"One panel per mark\" overlays conditions within each mark (best for how one mark changes across the life cycle). \"One panel per condition\" overlays marks within each condition (matches the Galaxy chromstaR output; best for which marks dominate at a stage)."),
                               tags$li(tags$b("\"Ratio calculation\""), " — \"Ratio of means\" (recommended) averages the raw signal, empty bins included, then takes one log ratio against the genome-wide mean. \"Legacy\" averages per-bin logs, which drops empty bins and flattens real depletion; it exists only to reproduce older figures."),
                               tags$li(tags$b("\"Upstream of TSS\" / \"Downstream of TES\" (bp)"), " and ", tags$b("\"Bins per region\""), " — window size and resolution."),
                               tags$li(tags$b("Gene scope, gene box, marks, stages"), " — as in the Metagene tab."),
                               tags$li(tags$b("\"Smooth curve (LOESS)\""), " — smooths each curve."),
                               tags$li(tags$b("\"Compute Enrichment\""), " — draws the plot."),
                               tags$li(tags$b("Subtitle"), " — reports the total number of genes contributing, and the range per position. Quote the total as your sample size."),
                               tags$li(tags$b("\"Download Plot\" / \"Download Data (Excel)\""), " — the figure as PDF, or per-position log ratios, gene counts and contributing gene names as Excel.")
                             ),
                             
                             hr(),
                             tags$h3("3.5 Bivalent tab"),
                             tags$p(tags$b("Purpose:"), " treat two or more marks as one combination (for example H3K4me3 with H3K27me3 at bivalent promoters), compare it with each mark alone around the TSS, gene body and TES, and list the genes that carry it."),
                             tags$ul(
                               tags$li(tags$b("\"Which combination?\""), " — tick the marks yourself, or pick a state that occurs in your data, ranked by how many bins carry it."),
                               tags$li(tags$b("Co-occurrence"), " — uses chromstaR's own state calls (", tags$code("combination.<condition>"), "). Each line is log(observed/expected) for the fraction of bins in that state. This is the mode that shows real co-occurrence, because the marks are called on the same bins."),
                               tags$li(tags$b("Sum"), " — adds the RPKM of the chosen marks. \"Normalise\" gives each mark equal weight; \"Raw\" lets the strongest mark dominate. A high sum can come from one mark alone, so a sum is not evidence of co-occurrence."),
                               tags$li(tags$b("\"Contains all\" vs \"Exactly\""), " — \"Contains all\" counts any bin with the chosen marks, whatever else is present. \"Exactly\" requires the chosen marks and nothing else; single-mark lines then mean that mark alone."),
                               tags$li(tags$b("Gene calls"), " — a gene is bivalent in a condition when any bin, or at least half the bins (your choice), of its promoter window is in the state. The table shows the fraction per condition and whether the combination is gained, lost, kept or absent between consecutive conditions in your life-cycle order."),
                               tags$li(tags$b("Y axis"), " — log(observed/expected), or the observed level. The observed axis shows losses directly; log(obs/exp) compares each position with how common the state is genome-wide in that condition."),
                               tags$li(tags$b("\"Download Data (Excel)\""), " — three sheets: settings, profile numbers and the full gene table.")
                             ),
                             
                             hr(),
                             tags$h3("4. Region Browser tab"),
                             tags$p(tags$b("Purpose:"), " raw signal across a chosen chromosome and position range, or across whole chromosomes — for checking a specific locus."),
                             tags$ul(
                               tags$li(tags$b("\"Chromosome scope\""), " — \"Single chromosome\" zooms to a bp range (the usual choice). \"Selected chromosomes\" or \"All chromosomes\" show full lengths; useful for an overview, but slow and cluttered with many contigs."),
                               tags$li(tags$b("\"Start (bp)\" / \"End (bp)\""), " — the exact window in single-chromosome mode."),
                               tags$li(tags$b("\"Jump to gene\""), " with ", tags$b("\"Flanking (bp)\""), " and ", tags$b("\"Go\""), " — pick a gene, set the flanking distance, and the start and end boxes are set to that gene's location."),
                               tags$li(tags$b("\"Bin scope\""), " — all bins including intergenic, only bins inside a gene, or only intergenic bins. This is the only place in the app that shows truly intergenic positions."),
                               tags$li(tags$b("\"Marks to display\""), " and ", tags$b("\"Condition comparison\""), " — which marks to draw, and side-by-side or overlaid conditions."),
                               tags$li(tags$b("\"Load Region\""), " — draws the view. Nothing appears until you press it."),
                               tags$li(tags$b("\"Download Plot\""), " — the current view as PDF.")
                             ),
                             
                             hr(),
                             tags$h3("5. Differential Peaks tab"),
                             tags$p(tags$b("Purpose:"), " count, per mark, how many chromatin segments are present in one condition but not another, and in which direction. Any number of conditions is supported."),
                             tags$ul(
                               tags$li(tags$b("\"Min differential score\""), " — keep only segments at or above this confidence (closer to 1 is stricter). Galaxy default 0.9999."),
                               tags$li(tags$b("\"Min merged region width (bp)\""), " — discard shorter segments. Galaxy default 300."),
                               tags$li(tags$b("\"Gene scope\""), " — whole genome (the default, matching the Galaxy tool), or only segments overlapping the genes you list."),
                               tags$li(tags$b("\"Conditions to compare\""), " — 2 to 5 conditions, in your life-cycle order."),
                               tags$li(tags$b("\"Comparison mode\""), " — all pairwise combinations, or one reference against each of the others (usually the right choice for a life-cycle baseline)."),
                               tags$li(tags$b("\"Ranked figure style\""), " — stacked (gained and lost parts) or diverging (lost to the left of zero, gained to the right)."),
                               tags$li(tags$b("\"Genomic zones\""), " — optionally restrict the counts to chosen zones (TSS, gene body thirds, flanks, TES, intergenic). Needs the gene annotation."),
                               tags$li(tags$b("\"Compute Differential Peaks\""), " — runs the filter and count, then draws the figures."),
                               tags$li(tags$b("Outputs"), " — the per-mark bar chart, a per-mark summary table with gained, lost, total and percentages, a zone breakdown, and two ranked figures: one for the filtered differential regions and one genome-wide from the object's frequency table."),
                               tags$li(tags$b("Reading the numbers"), " — \"total changed\" is gained plus lost, not a directional count, so quote the Gained or Lost column when you mean a direction. The percentages overlap and will not sum to 100 %, because one region can change for several marks at once."),
                               tags$li(tags$b("\"Download\" buttons"), " — plots as PDF, the per-mark counts and summary tables as Excel.")
                             ),
                             
                             hr(),
                             tags$h3("6. Gene Set Comparison tab"),
                             tags$p(tags$b("Purpose:"), " compare the promoter posterior probability (confidence that a mark is present) of one gene list against a background or a second list. Useful for asking whether a curated set is unusual for a mark, or whether up- and down-regulated genes differ in chromatin state."),
                             tags$ul(
                               tags$li(tags$b("\"Upstream / Downstream of TSS\""), " — the promoter window averaged for each gene."),
                               tags$li(tags$b("\"Replicates\""), " — average across replicates, or keep each one to check consistency."),
                               tags$li(tags$b("\"Summary statistic\""), " — mean or median drives the comparison; both are always shown in the table."),
                               tags$li(tags$b("\"Comparison type\""), " — (1) permutation against random genes from the whole annotation; (2) permutation with random draws only from your own background list B, which controls for a matched background such as expression-similar genes (set A's genes are removed from B by default, and B should be several times larger than A); (3) a direct Wilcoxon test of A against B."),
                               tags$li(tags$b("Gene set A and B boxes"), " — one ID per line. A match summary reports how many were found. Always check it before reading results."),
                               tags$li(tags$b("\"Number of random draws\""), " — permutation modes only. 1000 is a good default; use thousands for a final result."),
                               tags$li(tags$b("\"Compare Gene Sets\""), " — fills in the statistics table (sorted by FDR), the box plot and the per-gene table."),
                               tags$li(tags$b("\"Download\" buttons"), " — the box plot as PDF, or the full results as a multi-sheet Excel file.")
                             ),
                             
                             hr(),
                             tags$h3("7. RF Region Model and SHAP tabs"),
                             tags$p(tags$b("Purpose:"), " ask how well the histone marks at eight positions along a gene predict its expression, and which mark at which position pushes expression up or down. Each gene is described by 64 features (eight marks × eight regions: TSS −2 to −1 kb, TSS −1 kb to 0, TSS ±200 bp, three gene-body thirds, and TES 0–1 and 1–2 kb). A random forest predicts log₂(FPKM+1); TreeSHAP attributes each prediction to individual features."),
                             tags$ul(
                               tags$li(tags$b("RF Region Model — inputs"), " — pick two stage-and-condition pairs, set the forest size, CV folds and cores, and choose how many held-out genes to explain (0 = all; 2000 is a quick first run). Nothing runs until you press \"Run RF + SHAP\"."),
                               tags$li(tags$b("RF Region Model — outputs"), " — held-out performance (correlation and RMSE), a cross-validation summary showing whether performance depends on the split, a comparison of each single region with all regions combined, and the feature correlation matrix. Read the correlation matrix before interpreting SHAP: correlated features share importance."),
                               tags$li(tags$b("SHAP Mark x Region"), " — a bubble map (size = importance, colour = direction), a heatmap of the average push for expressed (FPKM ≥ 1) and silent genes, and a beeswarm of the top features."),
                               tags$li(tags$b("SHAP Gene Explorer"), " — choose a model and a gene; the waterfall shows how each of the largest feature contributions moves the prediction for that gene."),
                               tags$li(tags$b("\"Download all outputs (zip)\""), " — tables as CSV, the fitted models as RDS, and the figures as PDF."),
                               tags$li(tags$b("Caution"), " — this is an association and attribution study, not a causal one. A high SHAP value means the model relies on that feature, not that the mark regulates the gene.")
                             ),
                             
                             hr(),
                             tags$h3("8. Expression Deciles and Mark Combinations tabs"),
                             tags$ul(
                               tags$li(tags$b("Expression Deciles"), " — genes are ranked by mean normalised expression and split into ten groups. The gene-body profile of each mark is drawn per decile, showing how signal changes with expression. Uses the count table and colData."),
                               tags$li(tags$b("Mark Combinations"), " — an UpSet-style figure counting genes by the exact set of marks they carry at their promoter or across the gene body, optionally with each combination's expression level.")
                             ),
                             
                             hr(),
                             tags$h3("9. Data Table tab"),
                             tags$p(tags$b("Purpose:"), " the per-bin chromatin state table, annotated with the gene (if any) and genomic zone each bin falls into."),
                             tags$ul(
                               tags$li(tags$b("\"Download full table (CSV)\""), " — every row. Use this rather than the table's own export, which only handles small subsets reliably."),
                               tags$li(tags$b("Table"), " — scrollable and sortable; use the column boxes to filter by gene, zone or other fields.")
                             ),
                             
                             hr(),
                             tags$h3("Key concepts across several tabs"),
                             tags$ul(
                               tags$li(tags$b("Condition / life-cycle order"), " — set once on Load Data and applied everywhere. Conditions are drawn, coloured, faceted and exported in your order. A condition you do not place is appended rather than dropped."),
                               tags$li(tags$b("Stages to compare"), " — once Stage B is loaded, each analysis runs on each selected stage independently, then the results are combined and tagged by stage. The two objects need not share binning."),
                               tags$li(tags$b("Gene scope"), " — every profile includes only bins overlapping a gene, because there is no intergenic position relative to a TSS or gene body. Only the Region Browser shows intergenic bins."),
                               tags$li(tags$b("Strand"), " — profiles are read from the 5′ end of each gene. This needs the strand column in the gene table."),
                               tags$li(tags$b("RPKM vs log(observed/expected)"), " — Metagene and Region Browser show mean RPKM. Enrichment and Bivalent show log ratios against the genome-wide average, which compares marks with very different baseline levels."),
                               tags$li(tags$b("Compute buttons"), " — analyses are not recomputed when a setting changes, including appearance sliders. Press the tab's Compute button after any change."),
                               tags$li(tags$b("Posterior probability"), " — a 0–1 chromstaR confidence that a mark is present at a bin, separate from the RPKM signal.")
                             ),
                             
                             tags$h3("Performance tips"),
                             tags$ul(
                               tags$li("Most analyses stay fast with all genes selected, but very large gene sets with many marks can take a few seconds to a minute."),
                               tags$li("The Data Table's gene and zone annotation is computed once per stage and cached, so the first visit can take a little longer."),
                               tags$li("\"All chromosomes\" in the Region Browser and thousands of permutations in Gene Set Comparison are slow — preview with small values, then scale up."),
                               tags$li("The RF + SHAP tabs are the slowest. Use a small \"Test genes to explain\" value for a first run.")
                             ),
                             
                             hr(),
                             tags$h3("Contact"),
                             tags$p("Built for chromatin and epigenomics analysis in ", tags$em("Echinococcus multilocularis"), " by Janan Gawra."),
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
    counts         = NULL,   # normalised count table (genes x samples)
    counts_b       = NULL,   # Stage B normalised count table
    coldata        = NULL,   # sample metadata: sample, condition

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

  # ---- BIVALENT / MARK COMBINATIONS -------------------------------------------
  output$stage_selector_bv <- renderUI(stage_selector_ui("stages_bv"))

  output$bv_mark_ui <- renderUI({
    marks <- all_marks_reactive()
    req(length(marks) > 0)
    classic <- intersect(c("H3K4me3", "H3K27me3"), marks)
    sel <- if (length(classic) == 2) classic else utils::head(marks, 2)
    checkboxGroupInput("bv_marks", NULL, choices = marks, selected = sel)
  })

  bv_found_states <- reactive({
    stages <- active_stages()
    req(length(stages) > 0)
    sel <- resolve_selected_stages(input$stages_bv)
    stl <- lapply(stages[sel], function(st) {
      st$conditions <- stage_conditions_kept(st); st
    })
    found_combination_states(stl, min_marks = 2)
  })

  output$bv_found_ui <- renderUI({
    fs <- bv_found_states()
    if (is.null(fs) || nrow(fs) == 0) {
      return(tags$em(style = "color:#888",
        "No state with 2 or more marks found (or no combination.* columns on the bins)."))
    }
    ch <- stats::setNames(fs$state,
                          paste0(fs$state, "  —  ", format(fs$n_bins, big.mark = ","), " bins"))
    selectInput("bv_found_state", NULL, choices = ch, selected = fs$state[1], width = "100%")
  })

  # The marks that define "the combination", from whichever picker is active.
  bv_marks_used <- reactive({
    if (identical(input$bv_define, "found")) {
      s <- input$bv_found_state
      if (is.null(s) || !nzchar(s)) return(character(0))
      combo_mark_sets(s)[[1]]
    } else {
      input$bv_marks %||% character(0)
    }
  })

  output$gene_selector_bv <- renderUI({
    req(rv$genes)
    tagList(
      tags$small(style = "color:#888",
        paste0(length(gene_names_reactive()), " genes in the annotation — paste a list, one per row.")),
      br(), br(),
      tags$textarea(id = "selected_genes_text_bv", rows = 10,
        style = "width:100%; height:180px; font-family: monospace; font-size: 13px; resize: vertical;",
        placeholder = "Paste gene IDs here, one per line…", ""),
      uiOutput("gene_match_summary_bv")
    )
  })

  selected_genes_parsed_bv <- reactive({
    txt <- input$selected_genes_text_bv %||% ""
    raw <- trimws(strsplit(txt, "[\r\n,;\t]+")[[1]])
    raw[nchar(raw) > 0]
  })

  output$gene_match_summary_bv <- renderUI({
    req(rv$genes)
    pasted <- selected_genes_parsed_bv()
    if (length(pasted) == 0) return(tags$small(style = "color:#888", "No genes entered yet."))
    nms <- gene_names_reactive()
    n_ok <- sum(pasted %in% nms)
    bad  <- setdiff(pasted, nms)
    tagList(br(),
      tags$small(style = "color:#28a745", paste0("✓ ", n_ok, " of ", length(pasted), " gene names matched.")),
      if (length(bad) > 0) tags$div(tags$small(style = "color:#dc3545",
        paste0("✗ ", length(bad), " not found, e.g.: ", paste(utils::head(bad, 5), collapse = ", ")))))
  })

  observeEvent(input$select_all_genes_bv, {
    req(rv$genes)
    updateTextAreaInput(session, "selected_genes_text_bv",
                        value = paste(gene_names_reactive(), collapse = "\n"))
  })
  observeEvent(input$clear_all_genes_bv, {
    updateTextAreaInput(session, "selected_genes_text_bv", value = "")
  })

  bivalent_data <- eventReactive(input$run_bivalent, {
    req(rv$hmm, rv$genes)

    marks_sel <- bv_marks_used()
    if (length(marks_sel) < 2) {
      showNotification("Choose at least 2 marks to form a combination.",
                       type = "warning", duration = 8)
      return(NULL)
    }

    if (identical(input$bv_gene_scope, "selected")) {
      nms <- gene_names_reactive()
      sel <- selected_genes_parsed_bv(); sel <- sel[sel %in% nms]
      if (length(sel) == 0) {
        showNotification("Enter at least one valid gene name, or switch scope to \"All genes\".",
                         type = "warning", duration = 8)
        return(NULL)
      }
      genes_sel <- rv$genes[nms %in% sel]
    } else {
      genes_sel <- rv$genes
    }

    mode     <- input$bv_mode      %||% "cooccur"
    rule     <- input$bv_rule      %||% "contains"
    sum_norm <- input$bv_sum_norm  %||% "normalised"
    stage_names <- resolve_selected_stages(input$stages_bv)
    stages      <- active_stages()

    # Warn about a chosen mark missing from a stage rather than failing silently.
    for (nm in stage_names) {
      miss <- setdiff(marks_sel, stages[[nm]]$marks)
      if (length(miss) > 0 && length(stages[[nm]]$marks) > 0) {
        showNotification(paste0(nm, ": mark(s) not in this object: ", paste(miss, collapse = ", ")),
                         type = "warning", duration = 10)
      }
    }

    withProgress(message = paste0("Bivalent analysis (", length(genes_sel), " genes)…"), {
      incProgress(0.2, detail = "profile")
      prof <- compute_across_stages(stage_names, function(st) {
        compute_bivalent_profile(
          hmm        = st$hmm,
          genes_gr   = genes_sel,
          marks_sel  = marks_sel,
          conditions = stage_conditions_kept(st),
          mode       = mode,
          rule       = rule,
          sum_norm   = sum_norm,
          show_none  = isTRUE(input$bv_show_none),
          upstream   = input$bv_upstream,
          downstream = input$bv_downstream,
          n_bins     = input$bv_n_bins)
      })

      incProgress(0.5, detail = "gene calls")
      gene_tabs <- list(); gene_conds <- list()
      for (nm in stage_names) {
        st <- stages[[nm]]
        if (is.null(st)) next
        g <- compute_bivalent_genes(
          hmm        = st$hmm,
          genes_gr   = genes_sel,
          marks_sel  = marks_sel,
          conditions = stage_conditions_kept(st),
          rule       = rule,
          gene_rule  = input$bv_gene_rule %||% "any",
          prom_up    = input$bv_prom_up,
          prom_down  = input$bv_prom_down,
          sum_norm   = if (identical(mode, "sum")) sum_norm else NULL)
        if (is.null(g)) next
        gene_conds[[nm]] <- attr(g, "conditions")
        if (length(stage_names) > 1) g <- cbind(stage = nm, g, stringsAsFactors = FALSE)
        gene_tabs[[nm]] <- g
      }
    })

    if (is.null(prof) && length(gene_tabs) == 0) {
      showNotification(
        if (identical(mode, "sum"))
          "No signal found for the chosen marks — check the counts.rpkm columns for these marks."
        else
          "No combination.* columns were found on the bins, or the chosen state never occurs.",
        type = "error", duration = 12)
      return(NULL)
    }

    list(profile    = prof,
         genes      = if (length(gene_tabs)) dplyr::bind_rows(gene_tabs) else NULL,
         gene_conds = gene_conds,
         marks      = marks_sel,
         mode       = mode,
         rule       = rule,
         sum_norm   = sum_norm,
         n_bins     = input$bv_n_bins,
         upstream   = input$bv_upstream,
         downstream = input$bv_downstream,
         gene_rule  = input$bv_gene_rule %||% "any",
         prom_up    = input$bv_prom_up,
         prom_down  = input$bv_prom_down)
  })

  build_bivalent_plot <- reactive({
    res <- bivalent_data()
    req(res, res$profile)
    df <- res$profile
    n_bins <- res$n_bins
    df$x <- df$slot / n_bins
    df <- df[df$x >= -1 & df$x <= 2, , drop = FALSE]

    y_obs <- identical(input$bv_yaxis, "observed")
    if (y_obs) {
      df$y <- if (identical(res$mode, "sum")) df$observed else 100 * df$observed
      y_lab <- if (identical(res$mode, "sum")) {
        if (identical(res$sum_norm, "raw")) "Mean RPKM" else "Signal / genome-wide mean of each mark"
      } else "% of bins in the state"
    } else {
      df$y  <- df$log_ratio
      y_lab <- "log(observed/expected)"
    }
    df <- df[is.finite(df$y), , drop = FALSE]
    req(nrow(df) > 0)

    x_breaks <- c(-1, -0.5, 0, 0.25, 0.5, 0.75, 1, 1.5, 2)
    x_labels <- c(paste0("-", res$upstream, "bp"), paste0("-", round(res$upstream / 2), "bp"),
                  "TSS\n(0%)", "25%", "50%", "75%", "TES\n(100%)",
                  paste0("+", round(res$downstream / 2), "bp"), paste0("+", res$downstream, "bp"))

    multi_stage <- "stage" %in% colnames(df) && length(unique(df$stage)) > 1
    by_line     <- identical(input$bv_layout, "by_line")
    # Level order: combination first, then single marks, then "[]" (bind_rows
    # across two stages turns the factor into character, so rebuild it here).
    ln_chr  <- as.character(df$line)
    ord     <- order(match(df$line_kind, c("combination", "single", "none")))
    line_lv <- unique(ln_chr[ord])
    df$line <- factor(ln_chr, levels = line_lv)

    # Colours for the lines: the combination stands out, single marks use the
    # mark palette, the no-mark state is dark navy (as in the chromstaR figure).
    kinds   <- stats::setNames(df$line_kind[match(line_lv, ln_chr)], line_lv)
    singles <- line_lv[kinds == "single"]
    single_cols <- unname(mark_palette(gsub("^\\[|\\]$", "", singles)))
    line_cols <- c(stats::setNames(rep("#E7298A", sum(kinds == "combination")), line_lv[kinds == "combination"]),
                   stats::setNames(single_cols, singles),
                   stats::setNames(rep("#1B2A6B", sum(kinds == "none")), line_lv[kinds == "none"]))
    line_types <- stats::setNames(ifelse(kinds == "none", "dashed", "solid"), line_lv)

    body_shade <- data.frame(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf)

    if (by_line) {
      p <- ggplot(df, aes(x = x, y = y, colour = condition, group = condition))
    } else {
      p <- ggplot(df, aes(x = x, y = y, colour = line, linetype = line, group = line))
    }
    p <- p + geom_rect(data = body_shade, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       inherit.aes = FALSE, fill = "grey85", alpha = 0.35)
    p <- if (isTRUE(input$bv_smooth)) {
      p + geom_smooth(method = "loess", span = 0.2, se = FALSE, formula = y ~ x)
    } else {
      p + geom_line()
    }
    p <- p +
      geom_hline(yintercept = if (y_obs) numeric(0) else 0, linetype = "dashed", colour = "grey40") +
      geom_vline(xintercept = c(0, 1), colour = "grey30", linewidth = 0.6) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = 0.02))

    if (by_line) {
      p <- p + scale_colour_manual(values = condition_palette(), name = "Condition")
      p <- if (multi_stage) p + facet_grid(stage ~ line) else p + facet_wrap(~ line, ncol = 2)
    } else {
      leg <- if (identical(res$mode, "sum")) "Signal" else "Combination"
      p <- p + scale_colour_manual(values = line_cols, name = leg) +
               scale_linetype_manual(values = line_types, name = leg)
      p <- if (multi_stage) p + facet_grid(stage ~ condition) else p + facet_wrap(~ condition, ncol = 2)
    }

    mode_txt <- if (identical(res$mode, "sum")) {
      paste0("summed RPKM (", if (identical(res$sum_norm, "raw")) "raw" else "each mark normalised to its genome-wide mean", ")")
    } else {
      paste0("fraction of bins in each state, ",
             if (identical(res$rule, "exact")) "exact states" else "state contains the marks")
    }
    tot <- suppressWarnings(max(df$n_genes_total, na.rm = TRUE))

    p +
      labs(x = NULL, y = y_lab,
           title = paste0("Combination ", combo_label(res$marks), " vs single marks"),
           subtitle = paste0(mode_txt,
                             if (is.finite(tot)) paste0("  |  ", format(tot, big.mark = ","), " genes") else "",
                             "  |  grey band = gene body")) +
      theme_bw(base_size = 13) +
      theme(strip.background = element_rect(fill = "#34495e"),
            strip.text       = element_text(colour = "white", face = "bold"),
            panel.grid.minor = element_blank(),
            panel.spacing.x  = grid::unit(1.6, "lines"),
            legend.position  = "top")
  })

  output$bivalent_plot <- renderPlot({
    apply_plot_customization(build_bivalent_plot(), input, "bv_")
  })

  output$dl_bivalent_plot <- downloadHandler(
    filename = function() paste0("bivalent_profile_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(file, apply_plot_customization(build_bivalent_plot(), input, "bv_"),
                      width = 14, height = 8, device = "pdf")
    }
  )

  output$bivalent_summary <- renderUI({
    res <- bivalent_data()
    req(res, res$genes)
    g <- res$genes
    blocks <- lapply(names(res$gene_conds), function(nm) {
      conds <- res$gene_conds[[nm]]
      sub   <- if ("stage" %in% colnames(g)) g[g$stage == nm, , drop = FALSE] else g
      n_cov <- sum(sub$n_promoter_bins > 0)
      per   <- vapply(conds, function(cn) {
        paste0(cn, ": ", format(sum(sub[[paste0("bivalent_", cn)]]), big.mark = ","))
      }, character(1))
      trans <- if (length(conds) >= 2) vapply(seq_len(length(conds) - 1), function(k) {
        col <- paste0(conds[k], " -> ", conds[k + 1]); v <- sub[[col]]
        sprintf("%s → %s: %s lost, %s gained, %s kept", conds[k], conds[k + 1],
                format(sum(v == "lost"), big.mark = ","), format(sum(v == "gained"), big.mark = ","),
                format(sum(v == "kept"), big.mark = ","))
      }, character(1)) else character(0)
      tags$div(style = "margin-bottom:6px;",
        tags$b(if (length(res$gene_conds) > 1) paste0(nm, " — ") else "",
               "genes with ", combo_label(res$marks), " at the promoter: "),
        paste(per, collapse = "   |   "),
        tags$span(style = "color:#888", paste0("   (of ", format(n_cov, big.mark = ","),
                                               " genes with promoter bins)")),
        if (length(trans)) tags$div(style = "color:#555;", paste(trans, collapse = "   |   ")))
    })
    rule_txt <- paste0("Rule: ",
      if (identical(res$rule, "exact")) "exactly these marks" else "contains all these marks",
      "; gene called if ",
      if (identical(res$gene_rule, "half")) "at least half" else "any",
      " of its promoter bins (TSS −", res$prom_up, " / +", res$prom_down, " bp) are in the state.")
    tagList(blocks, tags$small(style = "color:#888", rule_txt))
  })

  bivalent_gene_table_filtered <- reactive({
    res <- bivalent_data()
    req(res, res$genes)
    g <- res$genes
    arrow_cols <- grep(" -> ", colnames(g), value = TRUE, fixed = TRUE)
    f <- input$bv_table_filter %||% "any"
    if (identical(f, "any")) {
      g <- g[g$n_conditions_bivalent > 0, , drop = FALSE]
    } else if (identical(f, "changing")) {
      if (length(arrow_cols) == 0) {
        g <- g[0, , drop = FALSE]
      } else {
        ch <- Reduce(`|`, lapply(arrow_cols, function(cl) g[[cl]] %in% c("gained", "lost")))
        g  <- g[ch, , drop = FALSE]
      }
    }
    g[order(-g$n_conditions_bivalent, g$gene), , drop = FALSE]
  })

  output$bivalent_gene_table <- renderDT({
    g <- bivalent_gene_table_filtered()
    drop_cols <- c("chr", "start", "end")
    disp <- g[, setdiff(colnames(g), drop_cols), drop = FALSE]
    datatable(disp, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_bivalent_xlsx <- downloadHandler(
    filename = function() paste0("bivalent_analysis_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      res <- bivalent_data()
      req(res)
      prof <- res$profile
      if (!is.null(prof)) {
        prof$x_position <- prof$slot / res$n_bins
        prof$line <- as.character(prof$line)
      }
      settings <- data.frame(
        setting = c("combination", "mode", "state rule", "sum normalisation",
                    "gene call rule", "promoter upstream (bp)", "promoter downstream (bp)",
                    "profile upstream (bp)", "profile downstream (bp)", "bins per region"),
        value   = c(combo_label(res$marks), res$mode, res$rule,
                    if (identical(res$mode, "sum")) res$sum_norm else "n/a",
                    res$gene_rule, res$prom_up, res$prom_down,
                    res$upstream, res$downstream, res$n_bins),
        stringsAsFactors = FALSE)
      sheets <- list(settings = settings)
      if (!is.null(prof))      sheets$profile <- as.data.frame(prof)
      if (!is.null(res$genes)) sheets$genes   <- as.data.frame(res$genes)
      openxlsx::write.xlsx(sheets, file, overwrite = TRUE)
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

  # ---- Shared text summary for a loaded count table --------------------------
  count_table_text <- function(ct) {
    cat("Genes (rows):    ", format(ct$n_genes, big.mark = ","), "\n")
    cat("Samples (cols):  ", length(ct$samples), "\n")
    cat("Value range:     ", signif(ct$range[1], 3), "to", signif(ct$range[2], 3), "\n")
    cat("Median value:    ", signif(ct$median_value, 3), "\n")
    if (ct$n_dup_removed > 0) cat("Duplicate gene IDs removed: ", ct$n_dup_removed, "\n")
    if (ct$n_na > 0)          cat("Missing values (NA):        ", ct$n_na, "\n")
    if (!is.null(rv$coldata)) {
      cat("Samples found in colData:   ", sum(ct$samples %in% rv$coldata$sample), "of", length(ct$samples), "\n")
    }
    if (isTRUE(ct$looks_raw)) {
      cat("\nNote: every value is a whole number, so this looks like RAW counts.\n",
          "Deciles still run, but raw counts are not comparable between samples.\n",
          "DESeq2 normalised counts are preferred.\n", sep = "")
    }
  }

  # ---- Load colData (sample metadata) ----------------------------------------
  observeEvent(input$coldata_file, {
    req(input$coldata_file)
    withProgress(message = "Loading colData…", {
      tryCatch({
        rv$coldata <- load_coldata(input$coldata_file$datapath, input$coldata_file$name)
      }, error = function(e) {
        showNotification(paste("Error loading colData:", e$message), type = "error", duration = 10)
      })
    })
  })

  output$coldata_summary <- renderPrint({
    req(rv$coldata)
    cd <- rv$coldata
    cat("Samples: ", nrow(cd), "   Groups: ", length(unique(cd$condition)), "\n\n", sep = "")
    print(as.data.frame(table(Condition = cd$condition)), row.names = FALSE)
  })

  # ---- Load normalised count tables (Stage A and Stage B) --------------------
  observeEvent(input$counts_file, {
    req(input$counts_file)
    withProgress(message = "Loading count table (Stage A)…", {
      tryCatch({
        rv$counts <- load_count_table(input$counts_file$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading count table (Stage A):", e$message),
                         type = "error", duration = 10)
      })
    })
  })

  observeEvent(input$counts_file_b, {
    req(input$counts_file_b)
    withProgress(message = "Loading count table (Stage B)…", {
      tryCatch({
        rv$counts_b <- load_count_table(input$counts_file_b$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading count table (Stage B):", e$message),
                         type = "error", duration = 10)
      })
    })
  })

  output$counts_summary <- renderPrint({
    req(rv$counts)
    count_table_text(rv$counts)
  })

  output$counts_summary_b <- renderPrint({
    req(rv$counts_b)
    count_table_text(rv$counts_b)
  })

  # ---- Expression Deciles: dynamic selectors ---------------------------------
  output$decile_stage_selector <- renderUI({
    req(rv$hmm)
    ch <- "A"
    names(ch) <- stage_a_label()
    if (!is.null(rv$hmm_b)) {
      chb <- "B"
      names(chb) <- stage_b_label()
      ch <- c(ch, chb)
    }
    checkboxGroupInput("decile_stages", NULL, choices = ch, selected = "A")
  })

  output$decile_group_selector <- renderUI({
    cd <- rv$coldata
    if (is.null(cd)) {
      return(tags$small(style = "color:#888",
        "No colData loaded: all samples in the count table are used."))
    }
    present <- unique(c(rv$counts$samples, rv$counts_b$samples))
    groups  <- unique(cd$condition[cd$sample %in% present])
    req(length(groups) > 0)
    checkboxGroupInput("decile_groups", NULL, choices = groups, selected = groups)
  })

  output$decile_mark_selector <- renderUI({
    marks <- unique(c(rv$marks, rv$marks_b))
    req(length(marks) > 0)
    checkboxGroupInput("decile_marks", NULL, choices = marks, selected = marks)
  })

  output$decile_condition_selector <- renderUI({
    conds <- unique(c(if (length(rv$conditions) > 0) ordered_conditions(), rv$conditions_b))
    req(length(conds) > 0)
    checkboxGroupInput("decile_conditions", NULL, choices = conds, selected = conds)
  })

  # ---- Expression Deciles: compute on button press ---------------------------
  decile_result <- eventReactive(input$run_decile, {
    shiny::validate(
      need(rv$hmm,   "Load a ChromstaR object first (Load Data tab)."),
      need(rv$genes, "Load the gene annotation first (Load Data tab)."),
      need(length(input$decile_stages) > 0, "Select at least one stage."),
      need(length(input$decile_marks) > 0,  "Select at least one mark."),
      need(length(input$decile_conditions) > 0, "Select at least one condition.")
    )
    if ("A" %in% input$decile_stages) {
      shiny::validate(need(rv$counts, "Load the Stage A count table (Load Data tab)."))
    }
    if ("B" %in% input$decile_stages) {
      shiny::validate(need(rv$hmm_b, "Load the Stage B object (Load Data tab)."),
                      need(rv$counts_b, "Load the Stage B count table (Load Data tab)."))
    }

    withProgress(message = "Computing expression deciles…", value = 0.05, {
      genes     <- rv$genes
      prof_list <- list()
      dec_list  <- list()
      notes     <- character(0)
      n_st      <- length(input$decile_stages)

      for (st in input$decile_stages) {
        is_a   <- identical(st, "A")
        lab    <- if (is_a) stage_a_label() else stage_b_label()
        hmm    <- if (is_a) rv$hmm else rv$hmm_b
        ct     <- if (is_a) rv$counts else rv$counts_b
        st_conds <- intersect(input$decile_conditions, if (is_a) rv$conditions else rv$conditions_b)
        st_marks <- intersect(input$decile_marks,      if (is_a) rv$marks      else rv$marks_b)
        if (length(st_conds) == 0 || length(st_marks) == 0) {
          notes <- c(notes, paste0(lab, ": none of the selected marks/conditions exist in this object."))
          next
        }

        samp <- decile_samples_for(ct, rv$coldata, input$decile_groups)
        if (length(samp) == 0) {
          notes <- c(notes, paste0(lab, ": no samples in the chosen groups."))
          next
        }

        expr    <- gene_expression_values(genes, ct$mat, samp)
        keep    <- !is.na(expr) & expr >= input$decile_min_expr
        if (sum(keep) < 10) {
          notes <- c(notes, paste0(lab, ": fewer than 10 genes pass the expression filter."))
          next
        }

        genes_k <- genes[keep]
        e_k     <- expr[keep]
        grp     <- dplyr::ntile(e_k, 10L)
        mcols(genes_k)$expr_group <- grp

        prof <- compute_scaled_profile(
          hmm, genes_k,
          marks      = st_marks,
          conditions = st_conds,
          upstream   = input$decile_upstream,
          downstream = input$decile_downstream,
          n_flank    = input$decile_n_flank,
          n_body     = input$decile_n_body
        )
        if (is.null(prof)) {
          notes <- c(notes, paste0(lab, ": no ChromStaR bins overlap these genes."))
          next
        }
        prof$stage <- lab
        prof_list[[st]] <- prof

        f <- factor(grp, levels = 1:10)
        dec_list[[st]] <- data.frame(
          stage     = lab,
          decile    = 1:10,
          n_genes   = as.integer(table(f)),
          min_mean  = as.numeric(tapply(e_k, f, min)),
          max_mean  = as.numeric(tapply(e_k, f, max)),
          n_samples = length(samp),
          n_matched = sum(!is.na(expr)),
          n_kept    = sum(keep),
          stringsAsFactors = FALSE
        )
        incProgress(0.85 / n_st, detail = lab)
      }

      shiny::validate(need(length(prof_list) > 0,
        paste(c("No profile could be computed.", notes), collapse = "\n")))

      list(
        profile    = bind_rows(prof_list),
        deciles    = bind_rows(dec_list),
        notes      = notes,
        n_annot    = length(genes),
        min_expr   = input$decile_min_expr,
        n_flank    = input$decile_n_flank,
        n_body     = input$decile_n_body,
        upstream   = input$decile_upstream,
        downstream = input$decile_downstream
      )
    })
  })

  decile_plot_obj <- function(res) {
    shiny::validate(need(!is.null(res$profile) && nrow(res$profile) > 0,
      "No ChromStaR bins overlap these genes. Check that the gene annotation matches the ChromStaR object."))
    d <- res$deciles
    cap <- paste(vapply(unique(d$stage), function(s) {
      paste0(s, " — genes per decile (1 to 10): ",
             paste(d$n_genes[d$stage == s], collapse = " | "))
    }, character(1)), collapse = "\n")
    plot_decile_profile(res$profile, res$n_flank, res$n_body,
                        res$upstream, res$downstream,
                        expression_decile_colours(), caption = cap,
                        display = if (is.null(input$decile_display)) "separate" else input$decile_display)
  }

  output$decile_plot_ui <- renderUI({
    h <- if (is.null(input$dec_plot_height)) 720 else input$dec_plot_height
    w <- if (is.null(input$dec_plot_width))  100 else input$dec_plot_width
    plotOutput("decile_plot", height = paste0(h, "px"), width = paste0(w, "%"))
  })

  output$decile_plot <- renderPlot({
    res <- decile_result()
    apply_plot_customization(decile_plot_obj(res), input, "dec_")
  })

  output$decile_summary <- renderPrint({
    res <- decile_result()
    d <- res$deciles
    for (s in unique(d$stage)) {
      ds <- d[d$stage == s, ]
      cat(s, ": in count table ", format(ds$n_matched[1], big.mark = ","),
          " | kept (mean >= ", res$min_expr, ") ", format(ds$n_kept[1], big.mark = ","),
          " | samples used ", ds$n_samples[1], "\n", sep = "")
    }
    cat("\n")
    print(d[, c("stage", "decile", "n_genes", "min_mean", "max_mean")], row.names = FALSE)
    if (length(res$notes) > 0) cat("\n", paste(res$notes, collapse = "\n"), "\n", sep = "")
  })

  output$dl_decile_plot <- downloadHandler(
    filename = function() paste0("expression_decile_profile_", Sys.Date(), ".png"),
    content  = function(file) {
      res <- decile_result()
      p   <- apply_plot_customization(decile_plot_obj(res), input, "dec_")
      w   <- if (is.null(input$dec_dl_width))  12 else input$dec_dl_width
      h   <- if (is.null(input$dec_dl_height)) 8  else input$dec_dl_height
      ggsave(file, plot = p, width = w, height = h, dpi = 300, bg = "white")
    }
  )

  # ---- Mark Combinations (UpSet): selectors ----------------------------------
  upset_stage_hmm   <- reactive(if (identical(input$upset_stage, "B")) rv$hmm_b else rv$hmm)
  upset_stage_marks <- reactive(if (identical(input$upset_stage, "B")) rv$marks_b else rv$marks)
  upset_stage_conds <- reactive(if (identical(input$upset_stage, "B")) rv$conditions_b else rv$conditions)
  upset_stage_ct    <- reactive(if (identical(input$upset_stage, "B")) rv$counts_b else rv$counts)

  output$upset_cond_selector <- renderUI({
    conds <- upset_stage_conds()
    req(length(conds) > 0)
    selectInput("upset_cond", "Condition (cell type)", choices = conds, selected = conds[1])
  })

  output$upset_mark_selector <- renderUI({
    marks <- upset_stage_marks()
    req(length(marks) > 0)
    checkboxGroupInput("upset_marks", NULL, choices = marks, selected = marks)
  })

  output$upset_expr_groups <- renderUI({
    cd <- rv$coldata
    if (is.null(cd)) {
      return(tags$small(style = "color:#888",
        "Expression boxes use all samples of the count table (no colData loaded)."))
    }
    ct <- upset_stage_ct()
    samp <- if (is.null(ct)) character(0) else ct$samples
    groups <- unique(cd$condition[cd$sample %in% samp])
    req(length(groups) > 0)
    checkboxGroupInput("upset_groups", "Expression from these groups",
                       choices = groups, selected = groups)
  })

  # ---- Mark Combinations: compute on button press ----------------------------
  upset_result <- eventReactive(input$run_upset, {
    shiny::validate(
      need(upset_stage_hmm(), "Load the ChromStaR object for this stage (Load Data tab)."),
      need(rv$genes, "Load the gene annotation first (Load Data tab)."),
      need(length(input$upset_cond) == 1, "Select a condition."),
      need(length(input$upset_marks) >= 1, "Select at least one mark.")
    )
    withProgress(message = "Computing mark combinations…", {
      hmm   <- upset_stage_hmm()
      marks <- input$upset_marks
      P <- gene_mark_presence(hmm, input$upset_cond, rv$genes, marks,
                              window = input$upset_window,
                              up = input$upset_up, down = input$upset_down)
      shiny::validate(need(!is.null(P), "No combination states found for this condition."))

      expr_log <- NULL
      ct <- upset_stage_ct()
      if (isTRUE(input$upset_show_expr) && !is.null(ct)) {
        samp <- decile_samples_for(ct, rv$coldata, input$upset_groups)
        if (length(samp) > 0) {
          ev <- gene_expression_values(rv$genes, ct$mat, samp)
          expr_log <- log10(ev + 1)
        }
      }
      list(P = P, expr_log = expr_log, cond = input$upset_cond,
           stage = if (identical(input$upset_stage, "B")) stage_b_label() else stage_a_label(),
           n_genes = nrow(P), n_with_mark = sum(rowSums(P) > 0),
           expr_n = if (is.null(expr_log)) 0L else sum(!is.na(expr_log)))
    })
  })

  upset_plot_obj <- function(res) {
    shiny::validate(need(!is.null(res), "Compute the combinations first."))
    upset_figure(res$P, res$expr_log,
                 n_show = input$upset_n_show, min_n = input$upset_min_n,
                 base_size = input$upset_base_size,
                 title = paste0(res$cond, " (", res$stage, "): genes marked by histone acylation peaks"))
  }

  output$upset_plot_ui <- renderUI({
    h <- if (is.null(input$upset_plot_height)) 900 else input$upset_plot_height
    w <- if (is.null(input$upset_plot_width))  100 else input$upset_plot_width
    plotOutput("upset_plot", height = paste0(h, "px"), width = paste0(w, "%"))
  })

  output$upset_plot <- renderPlot({
    upset_plot_obj(upset_result())
  })

  output$upset_summary <- renderPrint({
    res <- upset_result()
    cat("Genes in annotation:       ", format(res$n_genes, big.mark = ","), "\n")
    cat("Genes with >= 1 mark:      ", format(res$n_with_mark, big.mark = ","), "\n")
    cat("Genes with expression:     ", format(res$expr_n, big.mark = ","), "\n")
    cat("\nSet size per mark (genes carrying the mark, any combination):\n")
    print(data.frame(mark = colnames(res$P), set_size = colSums(res$P)), row.names = FALSE)
  })

  output$dl_upset_plot <- downloadHandler(
    filename = function() paste0("mark_combinations_upset_", Sys.Date(), ".png"),
    content  = function(file) {
      p <- upset_plot_obj(upset_result())
      w <- if (is.null(input$upset_dl_width))  10 else input$upset_dl_width
      h <- if (is.null(input$upset_dl_height)) 9  else input$upset_dl_height
      ggsave(file, plot = p, width = w, height = h, dpi = 300, bg = "white")
    }
  )

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

    make_zone_gr <- function(starts, ends, names_vec, chr_vec, strand_vec = NULL) {
      ok <- starts >= 1 & ends >= 1 & starts <= ends
      gr <- GRanges(
        seqnames  = chr_vec[ok],
        ranges    = IRanges(pmax(1L, starts[ok]), pmax(1L, ends[ok])),
        gene_name = names_vec[ok],
        zone_idx  = which(ok)
      )
      if (!is.null(strand_vec)) {
        mcols(gr)$strand <- strand_vec[ok]
      }
      gr
    }

    zone_list <- list(
      TSS_pm200          = make_zone_gr(tss - 200,        tss + 200,        g_names, chr_g, str_g),
      gene_body_5prime   = make_zone_gr(s_g,              s_g + floor(glen/3),      g_names, chr_g, str_g),
      gene_body_middle   = make_zone_gr(s_g + floor(glen/3)+1, s_g + floor(2*glen/3), g_names, chr_g, str_g),
      gene_body_3prime   = make_zone_gr(s_g + floor(2*glen/3)+1, e_g,               g_names, chr_g, str_g),
      upstream_1000_0    = make_zone_gr(tss - 1000,       tss - 1,          g_names, chr_g, str_g),
      upstream_2000_1000 = make_zone_gr(tss - 2000,       tss - 1001,       g_names, chr_g, str_g),
      TES_0_1000         = make_zone_gr(tes + 1,          tes + 1000,       g_names, chr_g, str_g),
      TES_1000_2000      = make_zone_gr(tes + 1001,       tes + 2000,       g_names, chr_g, str_g)
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
    best_strand <- rep(NA_character_, n_bins)  # Track strand for label swapping

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
      best_strand[upd] <- mcols(zr)$strand[zi]  # Capture strand from zones
    }

    # Swap 5prime/3prime labels for minus-strand genes (strand-aware fix)
    for (i in which(!is.na(best_strand))) {
      if (best_strand[i] == "-") {
        if (best_zone[i] == "gene_body_5prime") {
          best_zone[i] <- "gene_body_3prime"
        } else if (best_zone[i] == "gene_body_3prime") {
          best_zone[i] <- "gene_body_5prime"
        }
      }
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

  # ---- Segment annotation (for Differential Peaks breakdown by genomic zone) ----
  # Similar to annotate_bins() but works on merged segments instead of individual bins.
  # Returns a data.frame with segment-level zone annotations.
  annotate_segments <- function(segs, genes) {
    req(segs, genes)

    # Use the existing annotate_bins logic but apply it to segments
    # Convert segments to GRanges if needed
    if (!inherits(segs, "GRanges")) {
      segs_gr <- GRanges(
        seqnames = seqs$seqnames,
        ranges   = IRanges(start = segs$start, end = segs$end)
      )
    } else {
      segs_gr <- segs
    }

    g_names <- mcols(genes)$gene_name
    s_g     <- start(genes)
    e_g     <- end(genes)
    str_g   <- as.character(strand(genes))
    glen    <- e_g - s_g

    tss <- ifelse(str_g == "-", e_g, s_g)
    tes <- ifelse(str_g == "-", s_g, e_g)
    chr_g <- as.character(seqnames(genes))

    make_zone_gr <- function(starts, ends, names_vec, chr_vec, strand_vec = NULL) {
      ok <- starts >= 1 & ends >= 1 & starts <= ends
      gr <- GRanges(
        seqnames  = chr_vec[ok],
        ranges    = IRanges(pmax(1L, starts[ok]), pmax(1L, ends[ok])),
        gene_name = names_vec[ok],
        zone_idx  = which(ok)
      )
      if (!is.null(strand_vec)) {
        mcols(gr)$strand <- strand_vec[ok]
      }
      gr
    }

    zone_list <- list(
      TSS_pm200          = make_zone_gr(tss - 200,        tss + 200,        g_names, chr_g, str_g),
      gene_body_5prime   = make_zone_gr(s_g,              s_g + floor(glen/3),      g_names, chr_g, str_g),
      gene_body_middle   = make_zone_gr(s_g + floor(glen/3)+1, s_g + floor(2*glen/3), g_names, chr_g, str_g),
      gene_body_3prime   = make_zone_gr(s_g + floor(2*glen/3)+1, e_g,               g_names, chr_g, str_g),
      upstream_1000_0    = make_zone_gr(tss - 1000,       tss - 1,          g_names, chr_g, str_g),
      upstream_2000_1000 = make_zone_gr(tss - 2000,       tss - 1001,       g_names, chr_g, str_g),
      TES_0_1000         = make_zone_gr(tes + 1,          tes + 1000,       g_names, chr_g, str_g),
      TES_1000_2000      = make_zone_gr(tes + 1001,       tes + 2000,       g_names, chr_g, str_g)
    )

    zone_priority <- c(
      TSS_pm200 = 1L, gene_body_5prime = 2L, gene_body_middle = 3L,
      gene_body_3prime = 4L, upstream_1000_0 = 5L, upstream_2000_1000 = 6L,
      TES_0_1000 = 7L, TES_1000_2000 = 8L
    )

    n_segs    <- length(segs_gr)
    best_pri  <- rep(9L,              n_segs)
    best_gene <- rep(NA_character_,   n_segs)
    best_zone <- rep("Intergenic",    n_segs)
    best_strand <- rep(NA_character_, n_segs)

    for (zname in names(zone_list)) {
      zr  <- zone_list[[zname]]
      if (length(zr) == 0) next
      ov  <- findOverlaps(segs_gr, zr, ignore.strand = TRUE)
      bi  <- queryHits(ov)
      pri <- zone_priority[zname]
      upd <- bi[best_pri[bi] > pri]
      if (length(upd) == 0) next
      zi  <- subjectHits(ov)[best_pri[bi] > pri]
      best_pri[upd]  <- pri
      best_gene[upd] <- mcols(zr)$gene_name[zi]
      best_zone[upd] <- zname
      best_strand[upd] <- mcols(zr)$strand[zi]
    }

    # Swap 5prime/3prime labels for minus-strand genes
    for (i in which(!is.na(best_strand))) {
      if (best_strand[i] == "-") {
        if (best_zone[i] == "gene_body_5prime") {
          best_zone[i] <- "gene_body_3prime"
        } else if (best_zone[i] == "gene_body_3prime") {
          best_zone[i] <- "gene_body_5prime"
        }
      }
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

    data.frame(
      gene_name    = best_gene,
      genomic_zone = zone_labels[best_zone],
      zone_type    = best_zone,
      stringsAsFactors = FALSE
    )
  }

  compute_differential_peaks <- function(hmm, marks, score_thresh, width_thresh,
                                         conds_sel = NULL, mode = "pairwise",
                                         ref = NULL, stage_label = NULL,
                                         genes_gr = NULL, genes_for_annotation = NULL) {

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

    # Annotate kept segments with genomic zones (if genes available)
    # Use genes_for_annotation if provided (for whole-genome zone annotation),
    # otherwise fall back to genes_gr (for selected-gene mode)
    segment_zones <- NULL
    genes_to_annotate <- genes_for_annotation %||% genes_gr
    if (!is.null(genes_to_annotate)) {
      kept_segs_gr <- GRanges(
        seqnames = kept$seqnames,
        ranges   = IRanges(start = kept$start, end = kept$end)
      )
      segment_zones <- annotate_segments(kept_segs_gr, genes_to_annotate)
    }

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
      changed_zones <- if (!is.null(segment_zones)) segment_zones[changed, ] else NULL

      for (m in marks) {
        in_a <- mark_in_combination(fa, m)
        in_b <- mark_in_combination(fb, m)

        # If we have zone info, break down by zone; otherwise just report totals
        if (!is.null(changed_zones) && nrow(changed_zones) > 0) {
          # Group results by zone
          unique_zones <- unique(changed_zones$zone_type)
          for (z in unique_zones) {
            z_idx <- changed_zones$zone_type == z
            z_label <- unique(changed_zones$genomic_zone[z_idx])[1]

            results[[length(results) + 1L]] <- data.frame(
              comparison = paste0(a, " vs ", b),
              mark       = m,
              genomic_zone = z_label,
              zone_type    = z,
              present_in = c(a, b),
              absent_in  = c(b, a),
              direction  = c(paste0(a, "-not-", b), paste0(b, "-not-", a)),
              n_regions  = c(sum(in_a[z_idx] & !in_b[z_idx]), sum(!in_a[z_idx] & in_b[z_idx])),
              n_changed  = sum(z_idx),
              stringsAsFactors = FALSE
            )
          }
        } else {
          # No zone info, just report totals
          results[[length(results) + 1L]] <- data.frame(
            comparison = paste0(a, " vs ", b),
            mark       = m,
            genomic_zone = "Genome-wide",
            zone_type    = "genome_wide",
            present_in = c(a, b),
            absent_in  = c(b, a),
            direction  = c(paste0(a, "-not-", b), paste0(b, "-not-", a)),
            n_regions  = c(sum(in_a & !in_b), sum(!in_a & in_b)),
            n_changed  = length(changed),
            stringsAsFactors = FALSE
          )
        }
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

  # Genes for zone annotation — always use all loaded genes for annotation,
  # regardless of the gene scope setting (diff_genes_selected filters for the
  # "selected genes" case, but we still want zone annotation in all cases)
  diff_genes_for_annotation <- reactive({
    rv$genes  # Return all genes if loaded, NULL otherwise
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

  output$diff_zone_filter_ui <- renderUI({
    zone_choices <- c(
      "TSS ±200bp" = "TSS_pm200",
      "Gene body 5′" = "gene_body_5prime",
      "Gene body middle" = "gene_body_middle",
      "Gene body 3′" = "gene_body_3prime",
      "Upstream 1k–0bp" = "upstream_1000_0",
      "Upstream 2k–1kbp" = "upstream_2000_1000",
      "TES 0–1kbp" = "TES_0_1000",
      "TES 1k–2kbp" = "TES_1000_2000",
      "Intergenic" = "Intergenic"
    )
    checkboxGroupInput(
      "diff_zone_filter",
      NULL,
      choices = zone_choices,
      selected = character(0),
      inline = FALSE
    )
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
    genes_gr    <- diff_genes_selected()  # For filtering (if "Only selected genes")
    genes_annot <- diff_genes_for_annotation()  # For zone annotation (always, if loaded)
    zone_filter <- input$diff_zone_filter %||% character(0)

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

        result_data <- compute_across_stages(stage_names, function(st) {
          df <- compute_differential_peaks(
            hmm                  = st$hmm,
            marks                = st$marks,
            score_thresh         = input$diff_score_thresh,
            width_thresh         = input$diff_width_thresh,
            conds_sel            = conds_sel,
            mode                 = mode,
            ref                  = ref,
            genes_gr             = genes_gr,
            genes_for_annotation = genes_annot
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

        # Apply zone filter if specified
        if (length(zone_filter) > 0 && !is.null(result_data) && nrow(result_data) > 0) {
          # Only apply zone filter if zone_type column exists and contains actual zone data
          if ("zone_type" %in% colnames(result_data)) {
            # Check if we have any real zone data (not just "genome_wide")
            has_zone_data <- any(!is.na(result_data$zone_type) & result_data$zone_type != "genome_wide")

            if (has_zone_data) {
              # Filter rows where zone_type is in the selected zones
              result_data <- result_data[result_data$zone_type %in% zone_filter, ]

              if (nrow(result_data) == 0) {
                showNotification(
                  paste0("No differential peaks found in the selected genomic zones. ",
                         "Make sure a gene annotation file is loaded on the Load tab."),
                  type = "warning", duration = 12
                )
                return(NULL)
              }

              # n_differential was counted over ALL zones. Recount it for the
              # selected zones only, so the subtitle and every percentage use the
              # regions actually on screen. Each segment sits in exactly one zone,
              # so per comparison the zone counts (n_changed) add up exactly. With
              # several pairs a stage keeps the largest per-pair count.
              result_data$n_differential_all_zones <- result_data$n_differential
              zc <- unique(result_data[, c("stage", "comparison", "zone_type", "n_changed")])
              per_cmp <- stats::aggregate(n_changed ~ stage + comparison, data = zc, FUN = sum)
              per_stg <- stats::aggregate(n_changed ~ stage, data = per_cmp, FUN = max)
              result_data$n_differential <- per_stg$n_changed[match(result_data$stage, per_stg$stage)]
            } else {
              # No zone data available (genes not loaded)
              showNotification(
                paste0("Zone filter cannot be applied: no gene annotation loaded. ",
                       "Upload a gene file on the Load tab to enable zone-based filtering."),
                type = "message", duration = 10
              )
              # Continue without zone filter
            }
          }
        }

        result_data

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

    # With a gene annotation loaded the data has one row per mark x direction x
    # genomic ZONE. geom_col(position_dodge) does not add rows that share a bar
    # position — it draws them on top of each other, so each bar used to show
    # only its largest zone. Collapse the zones first so each bar is the total.
    keep_cols <- intersect(c("stage", "stage_conditions", "total_filtered",
                             "n_differential", "comparison", "mark",
                             "present_in", "absent_in", "direction"),
                           colnames(df))
    df <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(keep_cols))) %>%
      dplyr::summarise(n_regions = sum(n_regions, na.rm = TRUE), .groups = "drop") %>%
      as.data.frame()

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

    # Zone filter info
    # Zone filter info. The zones actually applied are read from the data, not
    # the checkboxes, so ticking a box without pressing Compute changes nothing.
    zones_in_data <- unique(as.character(diffpeaks_data()$zone_type))
    zones_in_data <- zones_in_data[!is.na(zones_in_data)]
    filtered <- "n_differential_all_zones" %in% colnames(diffpeaks_data())
    zone_txt <- if (filtered) {
      paste0(" | zones: ", paste(zones_in_data, collapse = ", "))
    } else {
      " | all zones"
    }
    region_word <- if (filtered) " differential regions in these zones" else " differential regions"
    bar_tot <- tapply(df$n_regions, as.character(df$stage), sum, na.rm = TRUE)
    n_bar_total <- unname(bar_tot[as.character(info$stage)])

    subtitle_txt <- paste0(
      paste0(info$stage, " (", info$stage_conditions, ": ",
             format(info$n_differential, big.mark = ","), region_word,
             "; bars sum to ", format(n_bar_total, big.mark = ","),
             " mark changes)", collapse = "   |   "),
      "\nscore>=", input$diff_score_thresh,
      ", width>=", input$diff_width_thresh, "bp | ", mode_txt, scope_txt, zone_txt,
      "\nA region where 3 marks change is counted once per mark, so the bars add up to more than the number of regions."
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

  # Zone-aware summary table that breaks down differential peaks by genomic zone
  diffpeaks_zone_summary <- reactive({
    df <- diffpeaks_data()
    req(df)

    # If no zone information, return NULL
    if (!"genomic_zone" %in% colnames(df)) return(NULL)

    rows <- list()
    for (st in unique(as.character(df$stage))) {
      sdf <- df[as.character(df$stage) == st, , drop = FALSE]
      for (cmp in unique(as.character(sdf$comparison))) {
        cdf <- sdf[as.character(sdf$comparison) == cmp, , drop = FALSE]
        for (z in unique(as.character(cdf$genomic_zone))) {
          zdf <- cdf[as.character(cdf$genomic_zone) == z, , drop = FALSE]
          # n_changed is the same value repeated on every mark x direction row
          # of this zone, so take it once (summing it multiplied it by 2 x marks).
          n_zones <- suppressWarnings(max(zdf$n_changed, na.rm = TRUE))
          if (!is.finite(n_zones)) next
          if (n_zones == 0) next
          rows[[length(rows) + 1L]] <- data.frame(
            stage = st,
            comparison = cmp,
            genomic_zone = z,
            n_regions = n_zones,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(rows) == 0) return(NULL)
    dplyr::bind_rows(rows)
  })

  output$diffpeaks_zone_table <- renderDT({
    zone_summary <- diffpeaks_zone_summary()
    req(zone_summary)
    names(zone_summary) <- c("Stage", "Comparison", "Genomic Zone", "Number of Regions")
    DT::datatable(zone_summary, rownames = FALSE,
                  options = list(pageLength = 12, dom = "tip", scrollX = TRUE))
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

  output$gsc_b_title <- renderText({
    if (identical(input$gsc_mode, "pool")) "Background list B (genes to draw from)" else "Gene set B"
  })

  output$gsc_perm_help <- renderUI({
    src <- if (identical(input$gsc_mode, "pool")) "only from your background list B" else "from the full annotation"
    tags$small(style = "color:#888", paste0(
      "Each draw picks a random set of genes ", src, ", the same size as gene set A, and computes its ",
      "promoter posterior per mark. The empirical p-value is the fraction of those random draws as extreme as your real gene set A."))
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

    # "pool" is the permutation test with draws restricted to list B; below it
    # runs through the permutation path, and the result remembers which
    # background was used for the labels.
    background  <- if (identical(input$gsc_mode, "pool")) "pool" else "all"
    mode        <- if (identical(input$gsc_mode, "manual")) "manual" else "permutation"
    stage_names <- resolve_selected_stages(input$stages_gsc)
    stages      <- active_stages()

    sel_b <- character(0)
    if (identical(mode, "manual") || identical(background, "pool")) {
      sel_b <- gsc_genes_b_parsed(); sel_b <- sel_b[sel_b %in% nms]
      if (length(sel_b) == 0) {
        showNotification(
          "Gene set B needs at least one valid, matched gene (or switch to the whole-annotation random background).",
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
              summary_stat    = input$gsc_summary_stat,
              pool_labels     = if (identical(background, "pool")) sel_b else NULL,
              exclude_A_from_pool = isTRUE(input$gsc_pool_exclude_a)
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
        n_pool_overlap_A = vapply(per_stage, function(x) as.numeric(x$n_pool_overlap_A %||% 0), numeric(1)),
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

    # Tiny background pools give coarse, overlapping draws: say so once.
    if (identical(background, "pool")) {
      small <- n_info$stage[n_info$n_total < 3 * n_info$n_A]
      if (length(small) > 0) {
        showNotification(paste0(
          "Background list B is less than 3\u00d7 the size of set A (", paste(small, collapse = ", "),
          "): the random draws overlap a lot, so the p-values are less informative. A larger matched list is better."),
          type = "warning", duration = 12)
      }
    }

    list(
      error       = NULL,
      mode        = mode,
      background  = background,
      exclude_A   = isTRUE(input$gsc_pool_exclude_a),
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
      pool_word <- if (identical(res$background, "pool")) "background list B" else "universe"
      paste0(
        if (identical(res$background, "pool"))
          paste0("Random draws taken only from gene set B",
                 if (isTRUE(res$exclude_A)) " (set A's genes removed from B)" else "", ".  ") else "",
        paste(apply(res$n_info, 1, function(r)
          paste0(r["stage"], ": n=", r["n_A"], " (", pool_word, " ", r["n_total"], " genes",
                 if (identical(res$background, "pool") && as.numeric(r["n_pool_overlap_A"]) > 0)
                   paste0("; ", as.numeric(r["n_pool_overlap_A"]), " genes were in both A and B") else "",
                 ")")), collapse = "  |  "))
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
      rand_lab <- if (identical(res$background, "pool")) "Random from B" else "Random"
      df$group <- ifelse(df$group == "Random", rand_lab, df$group)
      group_colors <- stats::setNames(c("#4393c3", "#999999"), c("A", rand_lab))
      df$group <- factor(df$group, levels = c("A", rand_lab))
      n_txt <- paste(apply(res$n_info, 1, function(r) paste0(r["stage"], ": n=", r["n_A"])), collapse = "  |  ")
      subtitle_txt <- paste0(n_txt, " genes in set A (vs. one representative random draw of the same size",
                            if (identical(res$background, "pool")) ", taken from list B" else "", ") — ",
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
        settings <- data.frame(
          setting = c("comparison", "random draws taken from", "set A removed from background",
                      "number of draws", "summary statistic"),
          value   = c("permutation test",
                      if (identical(res$background, "pool")) "background list B" else "whole annotation",
                      if (identical(res$background, "pool")) as.character(res$exclude_A) else "n/a",
                      input$gsc_n_perm, input$gsc_summary_stat),
          stringsAsFactors = FALSE)
        openxlsx::write.xlsx(
          list(settings = settings, stats = res$stats, per_gene_boxplot = res$long_scores,
               permutation_draws = res$perm_long, sizes = res$n_info),
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

  # =========================================================================
  # RF + SHAP REGION MODEL (uses objects already loaded in this app)
  # =========================================================================
  rs_choices <- reactive({
    stages <- active_stages()
    ch <- character(0)
    for (nm in names(stages)) {
      st <- stages[[nm]]
      for (cn in rs_kept_conditions(st, ordered_conditions()))
        ch[paste0(nm, ": ", cn)] <- paste0(st$key, "|", cn)
    }
    ch
  })
  
  observe({
    ch <- rs_choices()
    req(length(ch) >= 2)
    s1 <- isolate(input$rs_pick1); s2 <- isolate(input$rs_pick2)
    if (is.null(s1) || !s1 %in% ch) s1 <- ch[1]
    if (is.null(s2) || !s2 %in% ch || s2 == s1) s2 <- ch[min(2, length(ch))]
    updateSelectInput(session, "rs_pick1", choices = ch, selected = s1)
    updateSelectInput(session, "rs_pick2", choices = ch, selected = s2)
  })
  
  rs_spec <- function(val, label) {
    key  <- sub("\\|.*$", "", val)
    cond <- sub("^[^|]*\\|", "", val)
    st   <- Find(function(s) identical(s$key, key), active_stages())
    ct   <- if (identical(key, "B")) rv$counts_b else rv$counts
    if (is.null(ct)) stop("Load the normalised count table for ", label, " first (Load Data tab).")
    if (is.null(rv$coldata)) stop("Load the colData sample sheet first (Load Data tab).")
    samp <- intersect(rv$coldata$sample[rv$coldata$condition == cond], ct$samples)
    if (length(samp) == 0)
      stop("No samples of '", cond, "' in the count table for ", label, ".")
    list(label = label, hmm = st$hmm, cond = cond,
         marks = intersect(RS_MARKS, st$marks), ct = ct, samples = samp)
  }
  
  rs_results <- eventReactive(input$rs_run, {
    shiny::validate(
      need(rv$genes,   "Load the gene annotation first (Load Data tab)."),
      need(rv$coldata, "Load the colData sample sheet first (Load Data tab)."),
      need(isTRUE(input$rs_pick1 != input$rs_pick2), "Pick two different conditions.")
    )
    ch   <- rs_choices()
    vals <- c(input$rs_pick1, input$rs_pick2)
    cfg  <- list(ntree = input$rs_ntree, ntree_single = input$rs_ntree_single,
                 cv_folds = input$rs_folds, n_cores = input$rs_ncores,
                 n_shap = input$rs_n_shap %||% 0, seed = 123)
    tryCatch({
      specs <- lapply(vals, function(v) rs_spec(v, names(ch)[ch == v][1]))
      withProgress(message = "RF + SHAP", value = 0.05, {
        rs_run(specs, rv$genes, cfg,
               progress = function(m) setProgress(detail = m))
      })
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 12)
      NULL
    })
  })
  
  observe({
    req(rs_results())
    updateSelectInput(session, "rs_genes_stage", choices = names(rs_results()$fits))
  })
  
  output$rs_perf_tab <- renderTable({
    req(rs_results()); rs_results()$perf
  }, digits = 3)
  
  output$rs_cv_tab <- renderTable({
    req(rs_results())
    rs_results()$cv %>% group_by(stage) %>%
      summarise(across(c(r, RMSE, R2), list(mean = mean, sd = sd)), .groups = "drop")
  }, digits = 3)
  
  output$rs_single_plot <- renderPlot({
    req(rs_results()); rs_plot_single(rs_results()$single, rs_results()$perf)
  })
  output$rs_corr_plot <- renderPlot({
    req(rs_results()); rs_plot_corr(rs_results()$cor_all)
  })
  output$rs_bubble <- renderPlot({
    req(rs_results()); rs_plot_bubble(rs_results()$fs)
  })
  output$rs_heat <- renderPlot({
    req(rs_results()); rs_plot_heat(rs_results()$fs)
  })
  
  output$rs_beeswarm_ui <- renderUI({
    req(rs_results())
    lapply(names(rs_results()$fits), function(s)
      plotOutput(paste0("rs_bee_", make.names(s)), height = "520px"))
  })
  
  observe({
    req(rs_results())
    for (s in names(rs_results()$fits)) local({
      fit <- rs_results()$fits[[s]]
      output[[paste0("rs_bee_", make.names(s))]] <- renderPlot({ rs_plot_beeswarm(fit) })
    })
  })
  
  output$rs_waterfall <- renderPlot({
    req(rs_results(), input$rs_genes_stage, nzchar(input$rs_genes_txt))
    fit <- rs_results()$fits[[input$rs_genes_stage]]
    g   <- trimws(strsplit(input$rs_genes_txt, "[,; \n\t]+")[[1]])
    g   <- g[g %in% rownames(fit$S)]
    shiny::validate(need(length(g) > 0, "Gene not in the test set of this model."))
    rs_plot_waterfall(fit, g[1])
  })
  
  output$rs_zip <- downloadHandler(
    filename = function() paste0("RF_SHAP_outputs_", Sys.Date(), ".zip"),
    content = function(file) {
      req(rs_results())
      d <- file.path(tempdir(), "rf_shap_outputs")
      unlink(d, recursive = TRUE)
      rs_write_outputs(rs_results(), d)
      zip::zipr(file, list.files(d, full.names = TRUE))
    }
  ) 
  
  }
# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)
