# volcano_tools.R - v30 add-on: volcano plots of a chromatin-defined gene set
# using Wald results (log2FC, adjusted p) per comparison and stage.
# Sourced by app_v29.R after rf_shap_extra.R. Needs GenomicRanges (attached by
# the app), openxlsx and DT. Self-contained: does not call any app helper except
# apply_plot_customization(), which has a local fallback below.

VOLC_UP_COL   <- "#3B82C4"
VOLC_DOWN_COL <- "#E8590C"
VOLC_NS_COL   <- "#9E9E9E"
VOLC_BG_COL   <- "#D5D5D5"

vol_or <- function(a, b) if (is.null(a)) b else a

# ---- Plot styling -------------------------------------------------------------

#' Plot styling from the Plot Appearance sliders (text sizes, points, gridlines).
volc_style <- function(p, input, prefix = "volc_") {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)
  gv <- function(s, d) {
    v <- input[[paste0(prefix, s)]]
    if (is.null(v) || length(v) != 1) d else v
  }
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint"))
      p$layers[[i]]$aes_params$size <- gv("plot_point_size", 1.8)
  }
  p <- p + ggplot2::theme(
    plot.title   = ggplot2::element_text(size = gv("plot_title_size", 14), face = "bold", hjust = 0.5),
    axis.title   = ggplot2::element_text(size = gv("plot_axis_title_size", 12), face = "bold"),
    axis.text    = ggplot2::element_text(size = gv("plot_axis_text_size", 10)),
    legend.text  = ggplot2::element_text(size = gv("plot_legend_size", 10)),
    legend.title = ggplot2::element_text(size = gv("plot_legend_size", 10) + 1, face = "bold"))
  if (!isTRUE(gv("plot_show_grid", TRUE)))
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  p
}

# Use the app's apply_plot_customization() if it exists, otherwise the local one.
if (!exists("apply_plot_customization", mode = "function")) {
  apply_plot_customization <- function(p, input, prefix = "") volc_style(p, input, prefix)
}

# ---- Chromatin calls (self-contained) -------------------------------------------

#' Split chromstaR combination strings ("[H3K4me3+H3K27me3]") into mark vectors.
volc_combo_sets <- function(x) {
  x <- gsub("^\\s*\\[|\\]\\s*$", "", as.character(x))
  lapply(strsplit(x, "+", fixed = TRUE), function(v) {
    v <- trimws(v)
    sort(unique(v[nzchar(v)]))
  })
}

#' One condition's combination column from the bins, as a plain vector.
volc_get_combo <- function(hmm, cond) {
  mc  <- S4Vectors::mcols(hmm$bins)
  col <- paste0("combination.", cond)
  if (!col %in% colnames(mc)) return(NULL)
  v <- mc[[col]]
  if (methods::is(v, "Rle")) v <- as.vector(v)
  v
}

#' Per-bin TRUE/FALSE: does the bin's state contain the chosen marks
#' (rule "contains") or equal them exactly (rule "exact")?
volc_state_flags <- function(combo, marks_sel, rule = "contains") {
  f    <- as.factor(combo)
  sets <- volc_combo_sets(levels(f))
  code <- as.integer(f)
  lf <- vapply(sets, function(s) {
    if (identical(rule, "exact")) setequal(s, marks_sel) else all(marks_sel %in% s)
  }, logical(1))
  out <- lf[code]
  out[is.na(out)] <- FALSE
  out
}

#' Per-gene calls for the chosen marks in each condition. Same logic as the
#' Bivalent tab: fraction of promoter bins in the state, then a call by the
#' "any" or "half" rule. Row i is annotation row i of `genes`.
volc_gene_calls <- function(hmm, genes, conds, marks, rule, gene_rule, prom_up, prom_down) {
  bins <- hmm$bins
  prom <- suppressWarnings(promoters(genes, upstream = prom_up, downstream = prom_down))
  ov   <- suppressWarnings(findOverlaps(prom, bins, ignore.strand = TRUE))
  gi   <- queryHits(ov)
  bi   <- subjectHits(ov)
  n_g  <- length(genes)
  
  out <- data.frame(
    gene    = as.character(mcols(genes)$gene_name),
    gene_id = as.character(mcols(genes)$gene_id),
    chr     = as.character(seqnames(genes)),
    start   = start(genes),
    end     = end(genes),
    strand  = as.character(strand(genes)),
    stringsAsFactors = FALSE)
  out$locus <- paste0(out$chr, ":", out$start, "-", out$end)
  out$n_promoter_bins <- tabulate(gi, nbins = n_g)
  
  used <- character(0)
  for (cc in conds) {
    combo <- volc_get_combo(hmm, cc)
    if (is.null(combo)) next
    per_bin <- volc_state_flags(combo, marks, rule)
    frac <- rep(NA_real_, n_g)
    if (length(gi) > 0) {
      agg <- tapply(as.numeric(per_bin[bi]), gi, mean)
      frac[as.integer(names(agg))] <- as.numeric(agg)
    }
    thr  <- if (identical(gene_rule, "half")) 0.5 else 0
    call <- !is.na(frac) & frac > thr
    if (identical(gene_rule, "half")) call <- !is.na(frac) & frac >= 0.5
    out[[paste0("frac_", cc)]]     <- round(frac, 3)
    out[[paste0("bivalent_", cc)]] <- call
    used <- c(used, cc)
  }
  if (length(used) == 0) return(NULL)
  out
}

volc_set_flag <- function(g, conds, selection, cx, cy) {
  cl <- function(cc) g[[paste0("bivalent_", cc)]]
  n  <- nrow(g)
  f <- switch(selection,
              present = cl(cx),
              only    = Reduce(`&`, lapply(conds, function(cc) {
                if (identical(cc, cx)) cl(cc) else !cl(cc)
              })),
              gained  = !cl(cx) & cl(cy),
              lost    = cl(cx) & !cl(cy),
              rep(FALSE, n))
  f <- as.logical(f)
  f[is.na(f)] <- FALSE
  f
}

# ---- Reading Wald results -------------------------------------------------------

#' "PA_vs_PNA_wald.csv" -> "PA_vs_PNA"
volc_clean_label <- function(file_name) {
  lbl <- sub("\\.[^.]*$", "", file_name)
  lbl <- sub("[ _.-]*(wald|results?|deseq2?|shrunk|lfc)$", "", lbl, ignore.case = TRUE)
  trimws(lbl)
}

#' "PA_vs_PNA" -> numerator "PA", denominator "PNA"
volc_parse_contrast <- function(label) {
  m <- regmatches(label, regexec("^(.+?)[ _.-]+vs[ _.-]+(.+)$", label, ignore.case = TRUE))[[1]]
  if (length(m) == 3) list(num = m[2], den = m[3])
  else list(num = NA_character_, den = NA_character_)
}

volc_read_wald <- function(path, file_name) {
  ext <- tolower(tools::file_ext(file_name))
  if (ext %in% c("xlsx", "xlsm")) {
    df <- openxlsx::read.xlsx(path, sheet = 1)
  } else {
    first <- readLines(path, n = 1, warn = FALSE)
    n_tab <- lengths(regmatches(first, gregexpr("\t", first, fixed = TRUE)))
    n_com <- lengths(regmatches(first, gregexpr(",",  first, fixed = TRUE)))
    sep   <- if (n_tab > n_com) "\t" else ","
    df <- read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                     check.names = FALSE, stringsAsFactors = FALSE,
                     na.strings = c("NA", "NaN", ""))
  }
  # write.csv() of a DESeq2 table leaves the first header blank, so read.table()
  # turns the gene IDs into row names. Move them back into a real column.
  if (!identical(rownames(df), as.character(seq_len(nrow(df)))))
    df <- data.frame(gene = rownames(df), df, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  
  low <- tolower(trimws(names(df)))
  pick <- function(opts) {
    i <- match(opts, low); i <- i[!is.na(i)]
    if (length(i)) i[1] else NA_integer_
  }
  gcol <- pick(c("gene_name", "gene", "symbol", "gene_id", "x", ""))
  if (is.na(gcol)) gcol <- 1L
  lcol <- pick(c("log2foldchange", "log2fc", "lfc"))
  pcol <- pick(c("pvalue", "p.value", "pval", "p"))
  acol <- pick(c("padj", "fdr", "qvalue", "p.adj", "padj.bh"))
  if (is.na(lcol))
    stop("no log2FoldChange column (columns: ",
         paste(utils::head(names(df), 8), collapse = ", "), ")")
  if (is.na(acol) && is.na(pcol)) stop("no padj or pvalue column")
  
  out <- data.frame(
    gene   = trimws(as.character(df[[gcol]])),
    lfc    = suppressWarnings(as.numeric(df[[lcol]])),
    pvalue = if (is.na(pcol)) NA_real_ else suppressWarnings(as.numeric(df[[pcol]])),
    padj   = if (is.na(acol)) NA_real_ else suppressWarnings(as.numeric(df[[acol]])),
    stringsAsFactors = FALSE)
  out$padj_used <- ifelse(is.na(out$padj), out$pvalue, out$padj)
  out <- out[!is.na(out$gene) & nzchar(out$gene), , drop = FALSE]
  out <- out[!duplicated(out$gene), , drop = FALSE]
  
  lbl <- volc_clean_label(file_name)
  ct  <- volc_parse_contrast(lbl)
  list(label = lbl, num = ct$num, den = ct$den, file = file_name, data = out)
}

# ---- Matching Wald rows to annotation ---------------------------------------------

#' Adds ann_idx (annotation row, or NA) and in_set (TRUE if that gene is in the set).
volc_annotate <- function(w, set_flag, annot_name, annot_id) {
  m1  <- match(w$data$gene, annot_name)
  m2  <- match(w$data$gene, annot_id)
  idx <- ifelse(!is.na(m1), m1, m2)
  ok  <- !is.na(idx)
  w$data$ann_idx <- idx
  w$data$in_set  <- ok & set_flag[ifelse(ok, idx, 1L)]
  w
}

#' Direction label. Positive log2FC = higher in the numerator, so it is
#' "Down in <denominator>" (orange, right of zero).
volc_category <- function(d, p_cut, lfc_cut, num, den) {
  sig     <- d$padj_used < p_cut & abs(d$lfc) >= lfc_cut
  has_den <- !is.na(den)
  lab_dn  <- if (has_den) paste0("Down in ", den) else "Down (log2FC > 0)"
  lab_up  <- if (has_den) paste0("Up in ", den)   else "Up (log2FC < 0)"
  ifelse(!sig, "NS", ifelse(d$lfc > 0, lab_dn, lab_up))
}

# ---- Figures ----------------------------------------------------------------------

volc_plot <- function(df, lfc_cut, p_cut, title, subtitle, xlab) {
  labs_ <- unique(as.character(df$category))
  lv <- c(sort(labs_[startsWith(labs_, "Down")]),
          sort(labs_[startsWith(labs_, "Up")]),
          intersect(c("NS", "Not in set"), labs_))
  cols <- vapply(lv, function(l) {
    if (startsWith(l, "Down")) VOLC_DOWN_COL
    else if (startsWith(l, "Up")) VOLC_UP_COL
    else if (l == "NS") VOLC_NS_COL
    else VOLC_BG_COL
  }, character(1))
  df$category <- factor(df$category, levels = lv)
  # draw background first so set genes sit on top
  ord <- match(as.character(df$category),
               c("Not in set", "NS", setdiff(lv, c("NS", "Not in set"))))
  df  <- df[order(ord), , drop = FALSE]
  
  p <- ggplot(df, aes(x = lfc, y = neglog10p, colour = category)) +
    geom_point(size = 1.8, alpha = 0.7) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    labs(x = xlab, y = "-log10(adjusted p-value)", title = title, subtitle = subtitle) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#34495e"),
          strip.text = element_text(colour = "white", face = "bold"))
  if (length(unique(df$facet)) > 1)
    p <- p + facet_wrap(~ facet, ncol = min(3, length(unique(df$facet))))
  p
}

volc_concord <- function(a, b, p_cut, lfc_cut, title) {
  da <- a$data[a$data$in_set, c("gene", "lfc", "padj_used"), drop = FALSE]
  db <- b$data[b$data$in_set, c("gene", "lfc", "padj_used"), drop = FALSE]
  m  <- merge(da, db, by = "gene", suffixes = c("_a", "_b"))
  m  <- m[is.finite(m$lfc_a) & is.finite(m$lfc_b), , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  sa <- m$padj_used_a < p_cut & abs(m$lfc_a) >= lfc_cut
  sb <- m$padj_used_b < p_cut & abs(m$lfc_b) >= lfc_cut
  sa[is.na(sa)] <- FALSE
  sb[is.na(sb)] <- FALSE
  m$category <- ifelse(sa & sb, "Significant in both",
                       ifelse(sa | sb, "Significant in one", "NS"))
  rho <- suppressWarnings(stats::cor(m$lfc_a, m$lfc_b, method = "spearman"))
  ggplot(m, aes(x = lfc_a, y = lfc_b, colour = category)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
    geom_point(size = 1.8, alpha = 0.7) +
    scale_colour_manual(values = c("Significant in both" = "#7B2D8E",
                                   "Significant in one"  = "#F4A261",
                                   "NS"                  = "#BDBDBD"),
                        name = NULL) +
    labs(x = paste0(a$display, "  log2FC"),
         y = paste0(b$display, "  log2FC"),
         title = title,
         subtitle = sprintf("Spearman rho = %.3f, %d genes in the set with data in both",
                            rho, nrow(m))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

volc_set_table <- function(g, set_flag, files) {
  idx <- which(set_flag)
  out <- data.frame(gene = g$gene[idx], gene_id = g$gene_id[idx],
                    locus = g$locus[idx], n_promoter_bins = g$n_promoter_bins[idx],
                    stringsAsFactors = FALSE)
  for (cc in grep("^(frac|bivalent)_", names(g), value = TRUE))
    out[[cc]] <- g[[cc]][idx]
  for (w in files) {
    k <- match(idx, w$data$ann_idx)
    out[[paste0(w$display, " log2FC")]] <- round(w$data$lfc[k], 3)
    out[[paste0(w$display, " padj")]]   <- signif(w$data$padj_used[k], 3)
  }
  out
}

# ---- UI ---------------------------------------------------------------------------

volcano_tab_ui <- function() {
  fluidRow(
    box(title = "Volcano Settings", width = 3, status = "warning",
        tags$small(style = "color:#888",
                   "Pick genes that carry a chromatin feature, then overlay their Wald results (log2 fold change and adjusted p-value) from one or more comparisons."),
        hr(),
        h5("1. Gene set: chromatin feature"),
        uiOutput("volc_stage_ui"),
        uiOutput("volc_mark_ui"),
        radioButtons("volc_selection", "Genes to highlight",
                     choices = c("Carry the combination in one condition" = "present",
                                 "Carry it ONLY in that condition"        = "only",
                                 "Gained between two conditions"          = "gained",
                                 "Lost between two conditions"            = "lost"),
                     selected = "only"),
        uiOutput("volc_cond_ui"),
        conditionalPanel(
          condition = "input.volc_selection == 'gained' || input.volc_selection == 'lost'",
          uiOutput("volc_cond2_ui")),
        radioButtons("volc_rule", "State rule",
                     choices = c("Contains all chosen marks" = "contains",
                                 "Exactly the chosen marks"  = "exact"),
                     selected = "contains"),
        radioButtons("volc_gene_rule", "Gene is called when",
                     choices = c("Any promoter bin is in the state" = "any",
                                 "At least half of promoter bins are" = "half"),
                     selected = "any"),
        numericInput("volc_prom_up",   "Promoter upstream of TSS (bp)",   value = 1000, min = 0, step = 100),
        numericInput("volc_prom_down", "Promoter downstream of TSS (bp)", value = 1000, min = 0, step = 100),
        hr(),
        h5("2. Wald results"),
        tags$small(style = "color:#888",
                   "One file per comparison (csv, tsv, txt or xlsx; DESeq2 results() columns). Name the file like PA_vs_PNA.csv: the name sets the comparison and its direction (positive log2FC = higher in the first condition)."),
        br(), br(),
        fileInput("volc_wald_A", "Stage A comparisons", multiple = TRUE,
                  accept = c(".csv", ".tsv", ".txt", ".xlsx")),
        fileInput("volc_wald_B", "Stage B comparisons", multiple = TRUE,
                  accept = c(".csv", ".tsv", ".txt", ".xlsx")),
        tableOutput("volc_wald_summary"),
        actionButton("volc_wald_clear", "Clear all Wald results", class = "btn-sm btn-default"),
        hr(),
        h5("3. Figure"),
        checkboxGroupInput("volc_files", "Comparisons to plot", choices = NULL),
        radioButtons("volc_bg", "Other genes",
                     choices = c("Grey background" = "grey", "Hide them" = "hide"),
                     selected = "grey"),
        radioButtons("volc_plot_type", "Figure type",
                     choices = c("Volcano (one panel per comparison)" = "volcano",
                                 "Stage A vs Stage B log2FC (gene set)" = "concord"),
                     selected = "volcano"),
        conditionalPanel(
          condition = "input.volc_plot_type == 'concord'",
          selectInput("volc_cc_a", "Stage A comparison", choices = NULL),
          selectInput("volc_cc_b", "Stage B comparison", choices = NULL)),
        numericInput("volc_p_cut",   "Adjusted p-value cut-off", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("volc_lfc_cut", "|log2 fold change| cut-off", value = 1, min = 0, step = 0.1),
        box(title = HTML("<i class='fa fa-sliders-h'></i> Plot Appearance"),
            width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
            sliderInput("volc_plot_title_size", "Title size", min = 10, max = 20, value = 14, step = 1),
            sliderInput("volc_plot_axis_title_size", "Axis title size", min = 8, max = 18, value = 12, step = 1),
            sliderInput("volc_plot_axis_text_size", "Axis labels size", min = 6, max = 16, value = 10, step = 1),
            sliderInput("volc_plot_legend_size", "Legend text size", min = 8, max = 14, value = 10, step = 1),
            sliderInput("volc_plot_line_width", "Line width", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
            sliderInput("volc_plot_point_size", "Point size", min = 0.5, max = 5, value = 1.8, step = 0.1),
            sliderInput("volc_plot_alpha", "Transparency (fills)", min = 0.1, max = 1, value = 0.2, step = 0.1),
            checkboxInput("volc_plot_show_grid", "Show gridlines", value = TRUE)),
        hr(),
        actionButton("run_volcano", "Compute Volcano Plots", class = "btn-primary btn-block")
    ),
    box(title = "Volcano plot", width = 9, status = "primary",
        plotOutput("volc_plot", height = "560px"),
        tableOutput("volc_counts"),
        downloadButton("dl_volc_plot", "Download Plot (PDF)"),
        downloadButton("dl_volc_xlsx", "Download Genes & Results (Excel)")
    ),
    box(title = "Gene set with Wald results", width = 12, status = "primary",
        DT::DTOutput("volc_table"))
  )
}

# ---- Server -----------------------------------------------------------------------
# ctx: list of values/functions from the app server:
#   stages      = active_stages              (reactive, named list of stages)
#   genes       = reactive(rv$genes)
#   stage_conds = stage_conditions_kept      (function(st) -> conditions)
#   label_a     = stage_a_label              (reactive)
#   label_b     = stage_b_label              (reactive)

volcano_server <- function(input, output, session, ctx) {
  
  wald <- reactiveVal(list())
  
  add_wald <- function(stage_label, files) {
    if (is.null(files) || nrow(files) == 0) return(invisible(NULL))
    cur <- wald()
    ok  <- character(0)
    for (i in seq_len(nrow(files))) {
      res <- tryCatch(volc_read_wald(files$datapath[i], files$name[i]),
                      error = function(e) {
                        showNotification(paste0(files$name[i], ": ", conditionMessage(e)),
                                         type = "error", duration = 12)
                        NULL
                      })
      if (is.null(res)) next
      res$stage   <- stage_label
      res$display <- paste0(stage_label, ": ", res$label)
      cur[[res$display]] <- res
      ok <- c(ok, res$display)
    }
    wald(cur)
    if (length(ok))
      showNotification(paste("Loaded:", paste(ok, collapse = ", ")),
                       type = "message", duration = 8)
  }
  observeEvent(input$volc_wald_A, add_wald(ctx$label_a(), input$volc_wald_A))
  observeEvent(input$volc_wald_B, add_wald(ctx$label_b(), input$volc_wald_B))
  observeEvent(input$volc_wald_clear, wald(list()))
  
  output$volc_wald_summary <- renderTable({
    w <- wald()
    if (length(w) == 0)
      return(data.frame(Note = "No Wald results loaded yet.", stringsAsFactors = FALSE))
    data.frame(
      Comparison = names(w),
      Genes      = vapply(w, function(x) as.numeric(nrow(x$data)), numeric(1)),
      Contrast   = vapply(w, function(x)
        if (is.na(x$num)) "(not parsed)" else paste0(x$num, " / ", x$den),
        character(1)),
      stringsAsFactors = FALSE)
  })
  
  observe({
    w  <- wald()
    ch <- vol_or(names(w), character(0))
    updateCheckboxGroupInput(session, "volc_files", choices = ch, selected = ch)
    a <- ch[vapply(w, function(x) identical(x$stage, ctx$label_a()), logical(1))]
    b <- ch[vapply(w, function(x) identical(x$stage, ctx$label_b()), logical(1))]
    updateSelectInput(session, "volc_cc_a", choices = a, selected = if (length(a)) a[1] else NULL)
    updateSelectInput(session, "volc_cc_b", choices = b, selected = if (length(b)) b[1] else NULL)
  })
  
  output$volc_stage_ui <- renderUI({
    nm <- names(ctx$stages())
    req(length(nm) > 0)
    selectInput("volc_stage", "Chromatin stage that defines the gene set",
                choices = nm, selected = nm[1])
  })
  
  output$volc_mark_ui <- renderUI({
    stages <- ctx$stages()
    req(input$volc_stage)
    st <- stages[[input$volc_stage]]
    req(st)
    mk  <- st$marks
    sel <- intersect(c("H3K4me3", "H3K27me3"), mk)
    if (length(sel) < 1) sel <- utils::head(mk, 2)
    checkboxGroupInput("volc_marks", "Marks (combination)", choices = mk, selected = sel)
  })
  
  output$volc_cond_ui <- renderUI({
    stages <- ctx$stages()
    req(input$volc_stage)
    st <- stages[[input$volc_stage]]
    req(st)
    cn <- ctx$stage_conds(st)
    selectInput("volc_cond", "Condition", choices = cn, selected = cn[1])
  })
  
  output$volc_cond2_ui <- renderUI({
    stages <- ctx$stages()
    req(input$volc_stage)
    st <- stages[[input$volc_stage]]
    req(st)
    cn <- ctx$stage_conds(st)
    selectInput("volc_cond2", "Second condition", choices = cn,
                selected = cn[min(2, length(cn))])
  })
  
  vres <- eventReactive(input$run_volcano, {
    stages <- ctx$stages()
    validate(need(length(stages) > 0, "Load a ChromStaR object first (Load Data tab)."))
    validate(need(!is.null(input$volc_stage) && input$volc_stage %in% names(stages),
                  "Pick a stage."))
    st    <- stages[[input$volc_stage]]
    genes <- ctx$genes()
    validate(need(!is.null(genes), "Load the gene annotation first (Load Data tab)."))
    marks <- input$volc_marks
    validate(need(length(marks) >= 1, "Tick at least one mark."))
    conds <- ctx$stage_conds(st)
    sel   <- vol_or(input$volc_selection, "only")
    cx    <- input$volc_cond
    cy    <- input$volc_cond2
    validate(need(!is.null(cx) && cx %in% conds, "Pick a condition."))
    if (sel %in% c("gained", "lost"))
      validate(need(!is.null(cy) && cy %in% conds && cy != cx,
                    "Pick a second, different condition."))
    
    files_all <- wald()
    sel_files <- intersect(vol_or(input$volc_files, character(0)), names(files_all))
    validate(need(length(sel_files) >= 1,
                  "Upload Wald results (section 2) and tick at least one comparison."))
    
    p_cut   <- vol_or(input$volc_p_cut, 0.05)
    lfc_cut <- vol_or(input$volc_lfc_cut, 1)
    bg      <- vol_or(input$volc_bg, "grey")
    ptype   <- vol_or(input$volc_plot_type, "volcano")
    
    g <- volc_gene_calls(st$hmm, genes, conds, marks,
                         rule      = vol_or(input$volc_rule, "contains"),
                         gene_rule = vol_or(input$volc_gene_rule, "any"),
                         prom_up   = input$volc_prom_up,
                         prom_down = input$volc_prom_down)
    validate(need(!is.null(g), "No combination states were found for this stage."))
    set_flag <- volc_set_flag(g, conds, sel, cx, cy)
    validate(need(any(set_flag),
                  "No genes match this gene set. Try 'any promoter bin' or another condition."))
    
    annot_name <- as.character(mcols(genes)$gene_name)
    annot_id   <- as.character(mcols(genes)$gene_id)
    files <- lapply(files_all[sel_files], volc_annotate, set_flag = set_flag,
                    annot_name = annot_name, annot_id = annot_id)
    
    warn <- character(0)
    cat_rows <- list(); count_rows <- list(); facet_lab <- character(0)
    for (nm in names(files)) {
      w <- files[[nm]]
      if (sum(!is.na(w$data$ann_idx)) == 0) warn <- c(warn, w$display)
      d <- w$data[!is.na(w$data$lfc) & is.finite(w$data$lfc) &
                    !is.na(w$data$padj_used), , drop = FALSE]
      d$category <- volc_category(d, p_cut, lfc_cut, w$num, w$den)
      n_sd <- sum(d$in_set)
      n_dn <- sum(d$in_set & startsWith(d$category, "Down"))
      n_up <- sum(d$in_set & startsWith(d$category, "Up"))
      d$category[!d$in_set] <- "Not in set"
      count_rows[[nm]] <- data.frame(Comparison = w$display,
                                     `Set genes with data` = n_sd,
                                     Down = n_dn, Up = n_up,
                                     check.names = FALSE, stringsAsFactors = FALSE)
      facet_lab[nm] <- sprintf("%s\n%d genes with data | %d down | %d up",
                               w$display, n_sd, n_dn, n_up)
      if (identical(bg, "hide")) d <- d[d$in_set, , drop = FALSE]
      if (nrow(d) == 0) next
      cat_rows[[nm]] <- data.frame(
        file = nm, facet = facet_lab[nm], gene = d$gene, lfc = d$lfc,
        padj = d$padj_used,
        neglog10p = -log10(pmax(d$padj_used, 1e-300)),
        in_set = d$in_set, category = d$category,
        stringsAsFactors = FALSE)
    }
    if (length(warn))
      showNotification(paste("No gene names or IDs matched the annotation in:",
                             paste(warn, collapse = ", ")),
                       type = "warning", duration = 12)
    validate(need(length(cat_rows) > 0,
                  "None of the selected comparisons has genes with a valid log2FC and p-value."))
    
    df <- do.call(rbind, cat_rows)
    rownames(df) <- NULL
    df$facet <- factor(df$facet, levels = unique(facet_lab[names(cat_rows)]))
    
    nums <- unique(vapply(files, function(w) vol_or(w$num, NA_character_), character(1)))
    dens <- unique(vapply(files, function(w) vol_or(w$den, NA_character_), character(1)))
    common <- length(nums) == 1 && length(dens) == 1 && !is.na(nums) && !is.na(dens)
    xlab <- if (common) paste0("log2 Fold Change (", nums, " / ", dens, ")") else "log2 Fold Change"
    
    mk_txt <- paste(sort(marks), collapse = "+")
    if (length(marks) >= 2) mk_txt <- paste0("Bivalent ", mk_txt)
    sel_txt <- switch(sel,
                      present = paste0(cx, " genes"),
                      only    = paste0(cx, " only genes"),
                      gained  = paste0("gained ", cx, " -> ", cy),
                      lost    = paste0("lost ", cx, " -> ", cy),
                      "genes")
    rule_txt <- paste0(if (identical(input$volc_rule, "exact")) "exact states" else "state contains marks",
                       "; promoter: ",
                       if (identical(input$volc_gene_rule, "half")) "half of bins" else "any bin")
    title    <- paste0(mk_txt, " \u2014 ", sel_txt)
    subtitle <- paste0(input$volc_stage, "  |  ", sum(set_flag), " genes in the set  |  ",
                       rule_txt, "  |  adj. p < ", p_cut, ", |log2FC| >= ", lfc_cut)
    
    p_vol <- volc_plot(df, lfc_cut, p_cut, title, subtitle, xlab)
    
    p_cc <- NULL
    if (identical(ptype, "concord")) {
      a <- wald()[[vol_or(input$volc_cc_a, "")]]
      b <- wald()[[vol_or(input$volc_cc_b, "")]]
      validate(need(!is.null(a) && !is.null(b), "Pick a Stage A and a Stage B comparison."))
      a <- volc_annotate(a, set_flag, annot_name, annot_id)
      b <- volc_annotate(b, set_flag, annot_name, annot_id)
      p_cc <- volc_concord(a, b, p_cut, lfc_cut, paste0(title, ": stage A vs stage B"))
      validate(need(!is.null(p_cc), "No gene in the set has data in both comparisons."))
    }
    
    tbl <- volc_set_table(g, set_flag, files)
    
    settings <- data.frame(
      setting = c("stage", "marks", "gene set", "condition", "second condition",
                  "state rule", "promoter gene call", "promoter window (bp)",
                  "adjusted p cut-off", "|log2FC| cut-off", "other genes", "comparisons"),
      value = c(input$volc_stage, paste(sort(marks), collapse = "+"), sel_txt, cx,
                if (sel %in% c("gained", "lost")) cy else "n/a",
                vol_or(input$volc_rule, "contains"), vol_or(input$volc_gene_rule, "any"),
                paste0("-", input$volc_prom_up, " / +", input$volc_prom_down),
                p_cut, lfc_cut,
                if (identical(bg, "hide")) "hidden" else "grey background",
                paste(names(files), collapse = "; ")),
      stringsAsFactors = FALSE)
    
    list(plot = p_vol, cc = p_cc, plot_type = ptype, table = tbl,
         points = df[df$in_set, , drop = FALSE],
         counts = do.call(rbind, count_rows),
         n_facets = length(cat_rows), settings = settings)
  })
  
  output$volc_plot <- renderPlot({
    r <- vres()
    p <- if (identical(r$plot_type, "concord")) r$cc else r$plot
    apply_plot_customization(p, input, "volc_")
  })
  
  output$volc_counts <- renderTable({ vres()$counts })
  
  output$volc_table <- DT::renderDT({
    DT::datatable(vres()$table, rownames = FALSE, filter = "top",
                  options = list(pageLength = 15, scrollX = TRUE))
  })
  
  output$dl_volc_plot <- downloadHandler(
    filename = function() paste0("volcano_", Sys.Date(), ".pdf"),
    content = function(file) {
      r <- vres()
      p <- if (identical(r$plot_type, "concord")) r$cc else r$plot
      h <- if (identical(r$plot_type, "concord")) 7 else 4.5 * ceiling(r$n_facets / 3) + 1.5
      ggplot2::ggsave(file, apply_plot_customization(p, input, "volc_"),
                      width = 12, height = h, device = "pdf", limitsize = FALSE)
    }
  )
  
  output$dl_volc_xlsx <- downloadHandler(
    filename = function() paste0("volcano_genes_", Sys.Date(), ".xlsx"),
    content = function(file) {
      r <- vres()
      openxlsx::write.xlsx(list(settings = r$settings,
                                gene_set_results = r$table,
                                volcano_points = r$points,
                                comparisons = r$counts),
                           file, overwrite = TRUE)
    }
  )
}