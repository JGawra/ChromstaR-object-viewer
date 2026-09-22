# ChromstaR Viewer

An interactive Shiny application for exploring chromatin state analyses of ChIP-seq data.

ChromstaR Viewer opens a ChromstaR combined model — the object produced by the [chromstaR](https://bioconductor.org/packages/chromstaR/) Bioconductor package or by the Galaxy chromstaR workflow — and turns it into interactive figures and tables without writing any R code. Everything is computed live from the object you upload; nothing is precomputed or cached between sessions.

Developed for *Echinococcus multilocularis* epigenomics (8 histone modifications), but it works with any chromstaR combined model.

---

## Contents

- [What it does](#what-it-does)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input file 1 — the ChromstaR object](#input-file-1--the-chromstar-object)
- [Input file 2 — the gene annotation table](#input-file-2--the-gene-annotation-table)
- [The Load Data tab](#the-load-data-tab)
- [Controls that appear on every tab](#controls-that-appear-on-every-tab)
- [Metagene Profile](#metagene-profile)
- [Gene Contributions](#gene-contributions)
- [Enrichment Profile](#enrichment-profile)
- [Region Browser](#region-browser)
- [Differential Peaks](#differential-peaks)
- [Gene Set Comparison](#gene-set-comparison)
- [Data Table](#data-table)
- [Exports](#exports)
- [Troubleshooting](#troubleshooting)

---

## What it does

The app is built around one idea: a chromatin state is a combination of histone marks present at a genomic position, and the interesting biology is in how those combinations differ between conditions. Each tab answers one question about that difference.

Two datasets can be loaded at once, called **Stage A** and **Stage B**, so two independent comparisons can sit side by side in a single figure — for example a protoscolex comparison (PNA vs PA) next to a metacestode comparison (TwoDO vs Meta). Every tab that produces a figure facets by stage automatically when both are loaded.

| Tab | What it answers |
|---|---|
| **Load Data** | Which object and annotation am I working with, and in what condition order and colours? |
| **Metagene Profile** | Where does each mark sit relative to the TSS, the TES, or across the gene body? |
| **Gene Contributions** | Which individual genes produce the average metagene line, and how concentrated is the signal? |
| **Enrichment Profile** | How enriched or depleted is each mark around gene boundaries, relative to the genome-wide average? |
| **Region Browser** | What does the raw signal look like at one specific locus or chromosome? |
| **Differential Peaks** | Which marks change between conditions, in which direction, and in how many regions? |
| **Gene Set Comparison** | Does my gene list differ from a random background, or from a second list? |
| **Data Table** | What are the underlying per-bin values, and can I export them? |
| **About / Help** | In-app reference for each tab. |

---

## Installation

From CRAN:

```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr",
                   "tidyr", "scales", "shinyWidgets", "DT",
                   "magrittr", "openxlsx", "zoo"))
```

From Bioconductor:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer", "GenomicFeatures"))
```

---

## Quick start

```r
shiny::runApp("app_14_FIXED.R")
```

On startup a banner is printed to the R console. Check it whenever behaviour does not match this document — it is the quickest way to confirm which file is actually running:

```
---------------------------------------------------------------
ChromstaR Viewer build: 2026-09-21f / per-condition colour pickers
margin() resolves to: ggplot2   (must be 'ggplot2')
---------------------------------------------------------------
```

A first session that will surface most problems early:

1. Upload the ChromstaR object as **Stage A** and read the summary panel. Confirm the marks, conditions and chromosome names are what you expect.
2. Upload the gene table and confirm the gene count is what you expect.
3. Set the condition order to follow your biology, and set the condition colours.
4. Go to **Metagene Profile**, keep the defaults, press **Compute Profile**. If a sensible TSS profile appears, both input files are correct and the rest of the app will work.
5. Go to **Differential Peaks**, press **Compute Differential Peaks**, and read the per-mark summary table before interpreting any figure.

---

## Input file 1 — the ChromstaR object

The main input. Without it the app can do nothing.

### Accepted formats

`.RData`, `.rda`, `.Rdata`, `.rds`. For an `.RData` file the app searches the saved environment for an object of class `combinedMultiHMM`, `multiHMM` or `uniHMM` and uses the first one found; if none matches it falls back to the first object in the file. The object does not need a particular name, but a file containing several unrelated objects may load the wrong one.

Upload limit is **500 MB** per file.

### What the object must contain

To keep memory use manageable the app immediately discards everything it does not need and keeps five slots. If a slot is missing, the tabs that depend on it stop working but the rest of the app still runs.

| Slot | What it holds | Which tabs need it |
|---|---|---|
| `bins` | Genomic bins with RPKM signal, HMM posteriors and per-condition state combination | Metagene, Gene Contributions, Enrichment, Region Browser, Data Table, Gene Set Comparison |
| `segments` | Merged chromatin domains with `differential.score`, `width` and combination per condition | Differential Peaks |
| `frequencies` | Domain counts for each unique pair of state combinations | Genome-wide ranked figure on the Differential Peaks tab |
| `info` | Sample metadata: mark, condition, replicate | Automatic detection of marks and conditions on load |
| `hmms` | Fallback metadata where `info` is absent | Detection fallback only |

> [!IMPORTANT]
> Earlier builds discarded `frequencies` when trimming the object. If you are continuing a session started with an older build, re-upload the object so that slot is kept — otherwise the genome-wide ranked figure shows an explanatory note instead of a plot.

### Column naming the app depends on

Signal and state columns are located by name, so the object must follow the chromstaR convention. In `bins`:

```
counts.rpkm.<MARK>.<CONDITION>.<replicate>

e.g.  counts.rpkm.H3K4me3.Meta.rep1
```

In `segments` and `frequencies`:

```
combination.<CONDITION>

e.g.  combination.Meta   containing values like  [H3K4me3+H3K9ac]
```

Marks and conditions are detected from these names, so a mark appearing under an inconsistent spelling is treated as a separate mark.

### Stage A and Stage B

Stage A is required, Stage B optional. Each has a free-text label used in figure legends and facet strips — give them names that read well in a figure, since the defaults are generic. The two objects need not share conditions; every tab keeps them separate and labels them by stage.

---

## Input file 2 — the gene annotation table

Supplies gene coordinates. Optional for the Region Browser and the Data Table; **required** for every gene-anchored analysis (metagene, gene contributions, enrichment, gene set comparison).

### It must be a delimited table with a header — not a GTF

> [!WARNING]
> This is the most common source of confusion. The upload control is labelled "Upload genes.tsv file" and its file filter lists `.gtf` among accepted extensions, but the file is read with a plain table reader that expects a header row. **A real GTF file — nine unnamed columns with attributes packed into the ninth — will not load correctly.** Convert your annotation to a simple table first.

The separator is detected from the first line: tab if one is present, otherwise comma. Column names are lower-cased before matching, so capitalisation in the header does not matter.

### Required and optional columns

| Column | Required | Behaviour if absent |
|---|---|---|
| `chr` (or `seqnames`) | **Yes** | Loading fails. Values must match the chromosome names in the ChromstaR object. |
| `start` | **Yes** | Loading fails. |
| `end` | **Yes** | Loading fails. |
| `strand` | No | Every gene is treated as unstranded — see warning below. |
| `gene_name` | No | Falls back to `gene_id`; if that is absent too, a label of the form `chr:start-end` is generated. |
| `gene_id` | No | Falls back to `gene_name`. |

> [!WARNING]
> **Include the strand column.** Without it every gene loads as unstranded, so the app cannot tell the TSS from the TES. Minus-strand genes are then read backwards, and every profile anchored on the TSS or TES — metagene, gene contributions, enrichment — is smeared by mixing correctly and incorrectly oriented genes. The plots are still produced; they are simply wrong, with no warning.

### A minimal valid file

```
chr	start	end	strand	gene_name	gene_id
chr1	12400	15800	+	EmuJ_000123400	EmuJ_000123400
chr1	23100	26050	-	EmuJ_000123500	EmuJ_000123500
chr2	8700	11200	+	EmuJ_000201100	EmuJ_000201100
```

### Preparing the file from a GTF

```r
library(rtracklayer)
gtf   <- import("sample_with_utr_final.gtf")
genes <- subset(gtf, type == "gene")
genes <- genes[width(genes) >= 200]   # drop artefactual tiny features

out <- data.frame(
  chr       = as.character(seqnames(genes)),
  start     = start(genes),
  end       = end(genes),
  strand    = as.character(strand(genes)),
  gene_name = mcols(genes)$gene_name,
  gene_id   = mcols(genes)$gene_id)

write.table(out, "genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

After uploading, the panel beside the control reports how many genes were read. Check that against what you expect — a mismatch is usually a header or separator problem, and is far easier to catch now than to explain later in a figure caption.

### Chromosome names must match

Gene coordinates are intersected with the ChromstaR bins by chromosome name. If the annotation calls a sequence `chr1` and the ChIP-seq object calls it `1`, or if the annotation includes unplaced scaffolds that were never binned, those genes silently contribute nothing. Nothing errors — the gene counts in figure captions simply come out lower than the number of genes you loaded.

---

## The Load Data tab

The only tab that changes global settings. Everything set here applies everywhere at once.

### Uploads and summaries

Each upload is followed by a summary panel reporting the object class, the number of bins and segments, the marks and conditions found, and the chromosomes present. Read these before moving on. A missing mark or condition is almost always a naming mismatch in the object's columns rather than a problem with the app.

### Condition / life-cycle order

Conditions are detected alphabetically, which rarely matches the biology. Drag them into a meaningful sequence — for instance `TwoDO → Meta → PNA → PA`. That order then controls legend order, facet order, the order of exported columns, and which colour each condition receives, everywhere in the app.

Removing a condition with its `×` hides it from every analysis — a quick way to focus a figure on two of four conditions without re-uploading. **Reset to detected order** brings them all back.

### Condition colours

One colour swatch per condition. The colour you pick is used for that condition in **every** figure the app produces — metagene, region browser, differential peaks, gene set comparison and the ranked figures — so a condition looks identical across every panel of a multi-figure layout. Set your journal's or lab's scheme once and export every figure without editing them afterwards.

**Reset colours** restores the defaults, which are assigned by position in the condition order.

In the ranked figures each half of a bar means "the mark is present in this condition", so the half labelled *Gained in Meta* takes Meta's colour and the half labelled *Lost in Meta* takes the colour of the condition being compared against — because that is where the mark is present.

---

## Controls that appear on every tab

### The Compute button — and why nothing updates without it

Every analysis tab has a button such as **Compute Profile**, **Compute Enrichment**, **Load Region** or **Compute Differential Peaks**. These analyses intersect whole-genome bin tables with gene annotations and can take from seconds to a minute or two, so the app deliberately does not recompute when a setting changes.

> [!NOTE]
> Changing a setting and then looking at the plot shows you the **old** result. This includes the appearance sliders. After any change, press the tab's Compute button again.

### The Plot Appearance panel

Each plotting tab has a collapsible panel with eight controls. Each tab has its own independent set, so changing the title size on the metagene tab does not disturb the enrichment figure.

| Control | Range | What it does |
|---|---|---|
| Title size | 10–20 | Point size of the plot title |
| Axis title size | 8–18 | Point size of the x and y axis labels |
| Axis labels size | 6–16 | Point size of tick labels; any rotation already applied is preserved |
| Legend text size | 8–14 | Point size of legend entries; the legend title scales with it |
| Line width | 0.3–2.5 | Thickness of lines in line-based plots |
| Point/marker size | 1–5 | Size of plotted points |
| Transparency | 0.1–1.0 | Alpha of shaded bands. Deliberately not applied to heatmap tiles, which it would wash out |
| Show gridlines | on/off | Removes major and minor gridlines — useful for print figures |

These exist because a figure that reads well on screen rarely reads well at column width in a manuscript. All PDF exports honour the current settings, so set sizes before exporting rather than scaling afterwards.

### Gene scope and stage selectors

Most tabs offer a choice between all genes in the annotation and a list you select; choosing the second reveals a gene selector, and on some tabs **Select all** / **Clear all** buttons. Where two datasets are loaded, a stage selector controls whether the figure shows Stage A, Stage B or both.

---

## Metagene Profile

Averages signal for each mark across all selected genes, producing the classic profile around a reference point.

| Control | Purpose |
|---|---|
| Reference (one or more) | TSS, TES or gene body. Selecting several produces a panel each |
| Upstream / Downstream (bp) | Window either side of the reference point; a summary line reports the result |
| Number of bins | How finely the window is divided. More bins give more detail and more noise |
| Gene scope | All genes, or only genes you select |
| Marks to display | Which marks to draw |
| Condition comparison | Side-by-side puts each condition in its own panel; Overlay draws them on shared axes |
| Smooth curve (LOESS) | Smooths the profile. Recommended when many genes are included |

**Output:** a line plot of mean signal against position.
**Buttons:** *Download Plot* (PDF at current appearance settings) · *Download Data (Excel)* (the per-bin values behind it).

---

## Gene Contributions

A metagene line is an average, and averages hide their composition. This tab breaks one profile back down into the genes behind it — the honest check on whether a profile reflects a general trend or a handful of loci.

Choose a single stage, mark and condition, then the same reference and window controls as the metagene tab. **Genes to show in heatmap (top N by signal)** keeps the heatmap legible; **Log-scale the heatmap colours** compresses the range when a few genes dominate.

**Outputs:**

- A per-gene heatmap, one row per gene, sorted by total signal in the window. The column-wise average of this matrix *is* the metagene line from the previous tab.
- A cumulative contribution curve. A steep rise means a few genes carry the profile; a diagonal means every gene contributes roughly equally. A short text summary states how concentrated the signal is.
- A gene ranking table with mean and peak signal per gene, each gene's share of the total, and a running cumulative share.

**Buttons:** *Download Heatmap (PDF)* · *Download Curve (PDF)* · *Download Gene Table (Excel)* — use the last one when you need the loci driving a profile for downstream analysis.

---

## Enrichment Profile

Plots `log(observed/expected)` enrichment across a composite gene: upstream flank, gene body scaled 0–100 %, downstream flank. Expected is the genome-wide mean signal for that mark and condition, so zero means "as expected by chance" and the sign gives enrichment or depletion.

| Control | Purpose |
|---|---|
| Upstream of TSS / Downstream of TES (bp) | Size of the two flanks |
| Bins per region | Resolution of the composite profile |
| Gene scope | All genes, or a selected list |
| Marks to display / Stages to compare | Which marks and datasets to include |
| Panel layout | One panel per mark (conditions overlaid), or one panel per condition (marks overlaid) — the second matches the Galaxy chromstaR output |
| Ratio calculation | Ratio of means (recommended) or legacy per-bin log ratio — see below |
| Smooth curve (LOESS) | Smooths each curve |

### Why the ratio calculation matters

**Ratio of means** averages the raw signal at each position, zero bins included, then takes a single log ratio against the genome-wide mean. This is what the Galaxy chromstaR profile does and is the correct default.

**Legacy** takes the log of each bin separately and averages the logs. Because `log(0)` is undefined, empty bins must be discarded, which biases every position upward — most severely exactly where a mark is genuinely absent, since that is where most bins are zero. Real depletion is flattened away. The option exists only to reproduce older figures.

### Reading the gene count in the caption

The subtitle reports the total number of distinct genes contributing anywhere in the profile, and the range of genes contributing per position. These differ because a gene is counted at a position only if one of its bins falls there; short genes populate only some positions. **Quote the total as your sample size**; the range tells a reviewer how evenly the profile is supported.

**Buttons:** *Download Plot* · *Download Data (Excel)* (per-position log ratios, gene counts and contributing gene names).

---

## Region Browser

Shows raw signal at a locus rather than an average — the right tool for checking whether a genome-wide pattern is real at a place you know.

| Control | Purpose |
|---|---|
| Chromosome scope | A single chromosome zoomed to a range, selected chromosomes at full length, or all chromosomes |
| Start / End (bp) | The visible window when a single chromosome is selected |
| Jump to gene + **Go** | Centres the view on a named gene plus a flanking distance, instead of typing coordinates |
| Flanking (bp) | How much context to show either side when jumping |
| Bin scope | All bins including intergenic, only bins inside a gene, or only intergenic bins. **The only place in the app that can show truly intergenic positions** |
| Marks to display / Stages to compare | Which marks and datasets to draw |
| Condition comparison | Side-by-side panels or overlaid conditions |
| **Load Region** | Runs the query. Nothing is drawn until this is pressed |

**Buttons:** *Download Plot*.

---

## Differential Peaks

The quantitative core: which marks change between conditions, in which direction, and in how many regions. Reproduces the Galaxy chromstaR differential analysis and extends it.

| Control | Purpose |
|---|---|
| Min differential score | chromstaR's confidence threshold. Galaxy default `0.9999` |
| Min merged region width (bp) | Discards very short domains. Galaxy default `300` |
| Gene scope | Whole genome, or only segments overlapping selected genes |
| Conditions to compare | Which conditions enter the analysis, in your chosen order |
| Comparison mode | All pairwise combinations, or one reference condition against each of the others |
| Ranked figure style | Stacked (gained + lost) or Diverging (lost left of zero, gained right) |
| **Compute Differential Peaks** | Runs the analysis |

### How the filtering works, and the two numbers in the console

Filtering happens in two steps, and the console reports both:

1. Segments are kept if their differential score **and** width pass the thresholds.
2. Of those, only segments whose chromatin state actually **differs** between the two conditions are retained — these are the *differential regions* everything on the tab is built from.

The two counts are usually almost identical but need not be. chromstaR's differential score is a sum across marks of the difference in posterior probability between conditions, so a segment can accumulate enough total difference to pass the threshold without any single mark crossing the boundary that defines its state. Such a segment passes the score filter yet has an identical state in both conditions, and is correctly excluded. **The figure caption reports the number of differential regions — the number the bars are actually built from.**

### Outputs

- **Differential Peaks per Histone Mark** — grouped bars, one pair per mark, giving the number of regions where the mark is present in one condition and absent in the other, in each direction.
- **Per-mark summary table** — Gained, Lost and Total changed for each mark, each with its own percentage, plus the direction label.
- **Marks ranked by change — differential regions** — marks ranked by how many differential regions they change in, each bar split into gained and lost parts, coloured by the condition the mark is present in.
- **Marks ranked by change — genome-wide** — the same ranking across every domain in the genome with no score or width filter, from the object's frequency table. Percentages here are of all domains and are **not** comparable with the filtered figure. This panel needs the `frequencies` slot.

### Two things to be careful about when quoting these numbers

> [!CAUTION]
> **The total is not a directional count.** "Total changed" for a mark is gained **plus** lost; the Direction label only says which is larger. If H4K20me1 shows 30,284 total changed regions labelled "Gained in Meta", it does **not** mean 30,284 regions gained the mark — the split might be 24,229 gained and 6,055 lost. Always quote the Gained or Lost column, not the total, when you mean a direction.

> [!CAUTION]
> **The percentages overlap and will not sum to 100 %.** Most differential regions have several marks changing at once, and each region is counted once for every mark that changes in it. The note above the summary table reports the mean number of marks changing per region for your data.

**Buttons:** *Download Plot* · *Download Data (Excel)* · *Download Table (Excel)* · *Download Plot (PDF)* for each ranked figure.

---

## Gene Set Comparison

Tests whether a gene list of your own carries different promoter chromatin from a background — the tab to use when you arrive with a list from a differential expression analysis.

| Control | Purpose |
|---|---|
| Upstream / Downstream of TSS (bp) | Defines the promoter window scored for each gene |
| Stages to compare | Each stage is run independently and tagged in the results |
| Replicates | Average across replicates, or keep each separate — useful for checking replicate consistency |
| Summary statistic | Mean or median. Drives the delta and the permutation test's central value; both are always shown |
| Comparison type | A permutation test against random gene sets, or a direct comparison with a second list |
| Gene set A / Gene set B | Paste gene identifiers, one per line. A match summary reports how many were found |
| Number of random draws | Permutation mode only. Each draw takes a random gene set the same size as set A |
| **Compare Gene Sets** | Runs the analysis |

**Outputs:** a statistics table ranked by FDR (strongest, most significant differences at the top); a box plot of the distributions; and a per-gene table of promoter posterior probabilities for gene set A, showing which genes drive a significant result. In permutation mode the empirical p-value is the fraction of random draws as extreme as your real gene set.

> [!NOTE]
> Always read the match summary before interpreting anything. If only a fraction of your identifiers matched the annotation, the test ran on that fraction — identifier-format mismatches are the usual cause.

**Buttons:** *Download Plot* · *Download Data (Excel)*.

---

## Data Table

The per-bin values behind every other tab, as a searchable, sortable table. Use it to check a specific position or to confirm a column exists under the name you expect.

> [!NOTE]
> Use **Download full table (CSV)** for a full export. The table's own built-in export buttons only handle the currently displayed page reliably, which for several hundred thousand rows is almost never what you want.

---

## Exports

| Tab | Button | Contents |
|---|---|---|
| Metagene | Download Plot / Data (Excel) | PDF of the profile; per-bin mean values behind it |
| Gene Contributions | Download Heatmap / Curve (PDF) | PDF of the per-gene heatmap and the cumulative curve |
| Gene Contributions | Download Gene Table (Excel) | Per-gene mean and peak signal, share of total, cumulative share |
| Enrichment | Download Plot / Data (Excel) | PDF of the profile; per-position log ratios, gene counts, contributing genes |
| Region Browser | Download Plot | PDF of the current locus view |
| Differential Peaks | Download Plot / Data (Excel) | PDF of the per-mark bar chart; per-mark, per-direction region counts |
| Differential Peaks | Download Table (Excel) | Per-mark gained, lost, totals and all three percentages |
| Differential Peaks | Download Plot (PDF) ×2 | The filtered and genome-wide ranked figures |
| Gene Set Comparison | Download Plot / Data (Excel) | PDF of the box plot; statistics table and per-gene posteriors |
| Data Table | Download full table (CSV) | Every row of the per-bin table |

All PDF exports honour the current Plot Appearance settings.

---

## Troubleshooting

### A figure fails with an error mentioning a missing argument
If a plot area shows something like `argument "observed" is missing, with no default`, another package attached in your R session is masking a function the app uses. The known case is **randomForest**, which exports its own `margin()` and overrides the ggplot2 function of the same name.

The app guards against this specific case and reports the resolution in the startup banner. If the banner says `margin()` resolves to anything other than `ggplot2`, restart R and load fewer packages before starting the app. For any other masked function the error message names the failing call and the package it came from, and a full traceback is printed to the console.

### The genome-wide ranked figure shows a note instead of a plot
The object was loaded by an older build that discarded the `frequencies` slot. Re-upload the ChromstaR object; no other change is needed.

### A plot does not change when I move a slider
Expected behaviour — press the tab's Compute button. See [Controls that appear on every tab](#controls-that-appear-on-every-tab).

### Fewer genes are reported than I loaded
Usually chromosome naming. Genes on sequences absent from the ChromstaR bins — unplaced scaffolds are the common case — never overlap a bin and contribute nothing. Check that naming matches between the annotation and the object.

### Profiles look smeared or symmetric when they should not
Check that your gene table has a `strand` column. Without it every gene is treated as unstranded and minus-strand genes are read backwards.

### The app is slow
Whole-genome operations scale with the number of genes and bins. Narrow the gene scope, reduce the number of bins, or restrict the marks and conditions. The condition order box is a quick way to drop conditions you are not currently looking at.

---

## Author

**Janan Gawra** — IHPE (UMR 5244, CNRS – Université de Perpignan Via Domitia)

Application build `2026-09-21`.

Built for chromatin and epigenomics analysis by Janan Gawra.
[linkedin.com/in/janangawra](https://linkedin.com/in/janangawra)

---
