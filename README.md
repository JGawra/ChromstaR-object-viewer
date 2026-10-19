# ChromstaR Viewer

An interactive Shiny application for exploring chromatin state analyses of ChIP-seq data.

ChromstaR Viewer opens a ChromstaR combined model (the object produced by the [chromstaR](https://bioconductor.org/packages/3.20/bioc/html/chromstaR.html) Bioconductor package or by the Galaxy chromstaR workflow) and turns it into interactive figures and tables without writing R code. Everything is computed live from the object you upload. Nothing is precomputed or cached between sessions, except the optional saved RF + SHAP runs described below.

Developed for *Echinococcus multilocularis* epigenomics (8 histone modifications), but it works with any chromstaR combined model that uses the same mark set or a compatible one.

---

## Contents

- [What it does](#what-it-does)
- [Files you need](#files-you-need)
- [Installation](#installation)
- [Test data](#test-data)
- [Quick start](#quick-start)
- [Input file 1: the ChromstaR object](#input-file-1-the-chromstar-object)
- [Input file 2: the gene annotation table](#input-file-2-the-gene-annotation-table)
- [Input file 3: sample sheet and count tables](#input-file-3-sample-sheet-and-count-tables)
- [The Load Data tab](#the-load-data-tab)
- [Controls that appear on every tab](#controls-that-appear-on-every-tab)
- [Metagene Profile](#metagene-profile)
- [Gene Contributions](#gene-contributions)
- [Expression Deciles](#expression-deciles)
- [Enrichment Profile](#enrichment-profile)
- [Bivalent](#bivalent)
- [Region Browser](#region-browser)
- [Differential Peaks](#differential-peaks)
- [Gene Set Comparison](#gene-set-comparison)
- [Mark Combinations](#mark-combinations)
- [RF Region Model](#rf-region-model)
- [SHAP Mark x Region](#shap-mark-x-region)
- [RF Combined (Rmd)](#rf-combined-rmd)
- [SHAP Gene Explorer](#shap-gene-explorer)
- [Data Table](#data-table)
- [Exports](#exports)
- [Troubleshooting](#troubleshooting)

---

## What it does

The app is built around one idea: a chromatin state is a combination of histone marks present at a genomic position, and the interesting biology is in how those combinations differ between conditions. Each tab answers one question.

Two datasets can be loaded at once, called **Stage A** and **Stage B**, so two independent comparisons can sit side by side in one figure. For example, a protoscolex comparison (PNA vs PA) next to a metacestode comparison (TwoDO vs Meta). Every tab that produces a figure facets by stage automatically when both are loaded.

| Tab | What it answers |
|---|---|
| **Load Data** | Which objects, annotation and count tables am I working with, and in what condition order and colours? |
| **Metagene Profile** | Where does each mark sit relative to the TSS, the TES, or across the gene body? |
| **Gene Contributions** | Which individual genes produce the average metagene line, and how concentrated is the signal? |
| **Expression Deciles** | How does each mark's profile across the gene body change with expression level? |
| **Enrichment Profile** | How enriched or depleted is each mark around gene boundaries, relative to the genome-wide average? |
| **Bivalent** | Do two or more marks co-occur on the same bins, and which genes carry that combination in each condition? |
| **Region Browser** | What does the raw signal look like at one locus or chromosome? |
| **Differential Peaks** | Which marks change between conditions, in which direction, and in how many regions? |
| **Gene Set Comparison** | Does my gene list differ from a random background, a matched background, or a second list? |
| **Mark Combinations** | Which exact sets of marks do genes carry at their promoter or across the gene body? |
| **RF Region Model** | How well do the marks at eight positions along a gene predict its expression? |
| **SHAP Mark x Region** | Which mark at which position pushes expression up or down? |
| **RF Combined (Rmd)** | Cross-stage prediction, cross-validation, correlated features, and draft methods text for the same run. |
| **SHAP Gene Explorer** | Why does the model predict a given value for one gene? |
| **Data Table** | What are the underlying per-bin values, and can I export them? |
| **About / Help** | In-app reference for each tab. |

---

## Files you need

Keep these three files in the same folder:

| File | Role |
|---|---|
| `app_v29.R` | The Shiny app. Run this. |
| `rf_shap_tools.R` | RF + SHAP engine: region windows, feature table, models, SHAP, figures, output files. Sourced by the app at startup. |
| `rf_shap_extra.R` | Saving and loading RF + SHAP runs, cross-stage prediction, %IncMSE, correlated-pair tables, methods text, and extra figures. Sourced by the app at startup. |

Both `rf_shap_*.R` files are required. The app stops at startup if either is missing.

Optional, for testing: `make_test_data.R` (see [Test data](#test-data)).

---

## Installation

From CRAN:

```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr",
                   "tidyr", "scales", "shinyWidgets", "DT",
                   "magrittr", "openxlsx", "zoo", "patchwork"))
```

For the RF + SHAP tabs and the saved-run zip export:

```r
install.packages(c("randomForest", "caret", "treeshap", "ggbeeswarm",
                   "parallel", "zip"))
```

From Bioconductor:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "IRanges", "rtracklayer",
                       "GenomicFeatures", "S4Vectors"))
```

The RF + SHAP engine calls `randomForest`, `caret`, `treeshap`, `ggbeeswarm`, `parallel` and `zip` with `::` and does not attach them, so none of them can mask functions the app uses.

---

## Test data

A synthetic test set exercises every tab end to end. Generate it with `make_test_data.R` in the repository root:

```bash
Rscript make_test_data.R
```

This writes a folder `test_data/` containing:

| File | What it is |
|---|---|
| `stageA_protoscolex.rds` | ChromstaR-style object, conditions PNA and PA |
| `stageB_metacestode.rds` | ChromstaR-style object, conditions TwoDO and Meta |
| `genes.tsv` | Gene annotation with strand, for 4 placed chromosomes and 3 unplaced contigs |
| `normalised_counts.tsv` | Normalised expression, genes × samples |
| `colData.xlsx` | Sample sheet with `SampleName` and `Condition` |

The data is **synthetic**. It exercises every tab but does not represent real biology, so use it to check that the app runs and that figures look sensible, not to draw conclusions. In particular, Differential Peaks counts are not realistic estimates.

What the test data is designed to show:

- Metagene profiles with higher active marks at the TSS and across the gene body in expressed genes.
- RF + SHAP tabs in which active marks near the TSS and in the body contribute most, with a positive direction.
- The Bivalent tab recovering a planted set of promoters carrying H3K4me3 and H3K27me3 together in PNA, and a set gaining H3K4me3 in PA.
- Three unplaced contigs whose genes have no bins, so unmatched-gene behaviour can be checked.

---

## Quick start

```r
shiny::runApp("app_v29.R")
```

Keep `rf_shap_tools.R` and `rf_shap_extra.R` in the same folder as the app. The app sources both at startup, and it is started from that folder.

On startup a banner is printed to the R console. Check it whenever behaviour does not match this document. It confirms which build is running and which package `margin()` resolves to:

```
---------------------------------------------------------------
ChromstaR Viewer build: <date> / v29: RF + SHAP region model tabs
margin() resolves to: ggplot2   (must be 'ggplot2')
---------------------------------------------------------------
```

A first session using the test data:

1. Generate the test data with `Rscript make_test_data.R`.
2. Upload `test_data/stageA_protoscolex.rds` as **Stage A**, and `test_data/stageB_metacestode.rds` as **Stage B** (optional, but it exercises the stage features).
3. Upload `test_data/genes.tsv` as the gene annotation.
4. Upload `test_data/colData.xlsx` as the sample sheet, and `test_data/normalised_counts.tsv` as the count table for each stage. Set the count-table unit on the Load Data tab first.
5. Set the condition order to follow your biology, and set the condition colours.
6. Go to **Metagene Profile**, keep the defaults, press **Compute Profile**. A sensible TSS profile confirms the inputs are read correctly.
7. Go to **Bivalent**, tick H3K4me3 and H3K27me3, and press **Compute Bivalent Analysis**. The gene table should show the planted genes in PNA.
8. Go to **RF Region Model**, choose two conditions, set "Test genes to explain" to 2000, and press **Run RF + SHAP**.
9. Go to **Differential Peaks**, press **Compute Differential Peaks**, and read the per-mark summary table before interpreting any figure.

When you start with your own data, follow the same steps in the same order.

---

## Input file 1: the ChromstaR object

The main input. Without it the app can do nothing.

### Accepted formats

`.RData`, `.rda`, `.Rdata`, `.rds`. For an `.RData` file the app searches the saved environment for an object of class `combinedMultiHMM`, `multiHMM` or `uniHMM` and uses the first one found. If none matches it falls back to the first object in the file. A file containing several unrelated objects may load the wrong one.

Upload limit is **2 GB** per file.

### What the object must contain

To keep memory use manageable the app keeps five slots and discards the rest. If a slot is missing, the tabs that depend on it stop working but the rest of the app still runs.

| Slot | What it holds | Which tabs need it |
|---|---|---|
| `bins` | Genomic bins with RPKM signal, posteriors and per-condition state combination | Metagene, Gene Contributions, Expression Deciles, Enrichment, Bivalent, Region Browser, Data Table, Gene Set Comparison, Mark Combinations, RF Region Model |
| `segments` | Merged chromatin domains with `differential.score`, `width` and combination per condition | Differential Peaks |
| `frequencies` | Domain counts for each unique pair of state combinations | Genome-wide ranked figure on Differential Peaks |
| `info` | Sample metadata: mark, condition, replicate | Automatic detection of marks and conditions on load |
| `hmms` | Fallback metadata where `info` is absent | Detection fallback only |

> [!IMPORTANT]
> Earlier builds discarded `frequencies` when trimming the object. If you are continuing a session started with an older build, re-upload the object so that slot is kept. Otherwise the genome-wide ranked figure shows an explanatory note instead of a plot.

### Column naming the app depends on

Signal and state columns are located by name, so the object must follow the chromstaR convention. In `bins`:

```
counts.rpkm.<MARK>.<CONDITION>.<replicate>

e.g.  counts.rpkm.H3K4me3.Meta.rep1
```

Replicate suffixes are optional. Replicates are averaged per bin.

In `segments`, `frequencies`, and `bins` (for the Bivalent and Mark Combinations tabs):

```
combination.<CONDITION>

e.g.  combination.Meta   containing values like  [H3K4me3+H3K9ac]
```

Marks and conditions are detected from these names, so a mark spelled inconsistently is treated as a separate mark.

### Stage A and Stage B

Stage A is required, Stage B optional. Each has a free-text label used in figure legends and facet strips. Give them names that read well in a figure. The two objects need not share conditions, and every tab keeps them separate and labels them by stage.

---

## Input file 2: the gene annotation table

Supplies gene coordinates. Optional for the Region Browser and the Data Table. **Required** for every gene-anchored analysis (metagene, gene contributions, expression deciles, enrichment, bivalent, gene set comparison, mark combinations, RF Region Model).

### It must be a delimited table with a header, not a GTF

> [!WARNING]
> This is the most common source of confusion. The upload control is labelled "Upload genes.tsv file" and its file filter lists `.gtf` among accepted extensions, but the file is read with a plain table reader that expects a header row. **A real GTF file (nine unnamed columns with attributes packed into the ninth) will not load correctly.** Convert your annotation to a simple table first.

The separator is detected from the first line: tab if one is present, otherwise comma. Column names are lower-cased before matching.

### Required and optional columns

| Column | Required | Behaviour if absent |
|---|---|---|
| `chr` (or `seqnames`) | **Yes** | Loading fails. Values must match the chromosome names in the ChromstaR object. |
| `start` | **Yes** | Loading fails. |
| `end` | **Yes** | Loading fails. |
| `strand` | No | Every gene is treated as unstranded. See the warning below. |
| `gene_name` | No | Falls back to `gene_id`; if that is absent too, a label `chr:start-end` is generated. |
| `gene_id` | No | Falls back to `gene_name`. The RF + SHAP tabs match on `gene_id` first. |

> [!WARNING]
> **Include the strand column.** Without it every gene loads as unstranded, so the app cannot tell the TSS from the TES. Minus-strand genes are then read backwards, and every profile anchored on the TSS or TES (metagene, gene contributions, expression deciles, enrichment, bivalent, mark combinations, RF features) is smeared. The plots are still produced, with no warning.

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

After uploading, the panel beside the control reports how many genes were read. Check that against what you expect. A mismatch is usually a header or separator problem.

### Chromosome names must match

Gene coordinates are intersected with the ChromstaR bins by chromosome name. If the annotation calls a sequence `chr1` and the object calls it `1`, or if the annotation includes unplaced scaffolds that were never binned, those genes silently contribute nothing. Nothing errors. The gene counts in figure captions come out lower than the number of genes you loaded.

---

## Input file 3: sample sheet and count tables

Needed for the expression-based tabs: **Expression Deciles**, **Mark Combinations** (optional expression boxes), **RF Region Model**, **SHAP Mark x Region**, **RF Combined (Rmd)** and **SHAP Gene Explorer**.

### Sample sheet (colData)

A table with a sample column and a condition column. Accepted names (case-insensitive): `SampleName`, `sample`, `sample_name`, `sampleid` or `id` for samples; `Condition` or `group` for conditions. A `Replicate` column is allowed but not required. Accepted formats: `.xlsx`, `.xlsm`, `.csv`, `.tsv`, `.txt`.

Condition names must match the condition names in the ChromstaR object.

### Normalised count table

One row per gene, one numeric column per sample, first column gene IDs (its header may be blank). Comma-, tab- or semicolon-separated. Gene IDs must match `gene_id` or `gene_name` in the annotation.

The **unit** of the table is chosen on the Load Data tab before you upload the table: **FPKM** or **DESeq2 normalised counts**. The choice applies to tables loaded after it, so re-upload a table if you change the unit. The unit is printed in the RF output files and the methods text.

> [!NOTE]
> DESeq2 normalised counts are not divided by gene length, so long genes look more expressed. FPKM is the better target for the RF model. The app reports this in the count-table summary when DESeq2 counts are selected.

The app warns if every value is a whole number, which suggests raw counts. Raw counts are not comparable between samples.

---

## The Load Data tab

The only tab that changes global settings. Everything set here applies everywhere at once.

### Uploads and summaries

Each upload is followed by a summary panel reporting the object class, the number of bins and segments, the marks and conditions found, and the chromosomes present. Read these before moving on. A missing mark or condition is almost always a naming mismatch in the object's columns.

### Count-table unit

Choose FPKM or DESeq2 normalised counts **before** uploading the count tables.

### Condition / life-cycle order

Conditions are detected alphabetically, which rarely matches the biology. Drag them into a meaningful sequence, for example `TwoDO → Meta → PNA → PA`. That order controls legend order, facet order, the order of exported columns, and which colour each condition receives, everywhere in the app.

Removing a condition with its `×` hides it from every analysis. This is a quick way to focus a figure on two of four conditions without re-uploading. **Reset to detected order** brings them all back.

### Condition colours

One colour per condition. The colour is used for that condition in **every** figure the app produces, so a condition looks identical across every panel. **Reset colours** restores the defaults, assigned by position in the condition order.

In ranked figures each half of a bar means "the mark is present in this condition", so each half takes the colour of the condition in which the mark is present.

---

## Controls that appear on every tab

### The Compute button, and why nothing updates without it

Every analysis tab has a button such as **Compute Profile**, **Compute Enrichment**, **Load Region**, **Compute Bivalent Analysis**, **Compute Differential Peaks**, **Run RF + SHAP** or **Compute Combinations**. These analyses intersect whole-genome bin tables with gene annotations and can take from seconds to several minutes, so the app deliberately does not recompute when a setting changes.

> [!NOTE]
> Changing a setting and then looking at the plot shows you the **old** result. This includes the appearance sliders. After any change, press the tab's Compute button again.

### The Plot Appearance panel

Plotting tabs have a collapsible Plot Appearance panel. Each tab has its own independent set of controls.

| Control | Range | What it does |
|---|---|---|
| Title size | 10–20 | Point size of the plot title |
| Axis title size | 8–18 | Point size of the axis titles |
| Axis labels size | 6–16 | Point size of tick labels. Any rotation already applied is preserved |
| Legend text size | 8–14 | Point size of legend entries; the legend title scales with it |
| Line width | 0.3–2.5 | Thickness of lines in line-based plots |
| Point/marker size | 1–5 | Size of plotted points |
| Transparency | 0.1–1.0 | Alpha of shaded bands. Deliberately not applied to heatmap tiles, which it would wash out |
| Show gridlines | on/off | Removes major and minor gridlines, useful for print |

The Expression Deciles and Mark Combinations tabs also have a size panel (figure height and width on screen, and download width and height in inches). All PDF exports honour the current settings, so set sizes before exporting rather than scaling afterwards.

### Gene scope and stage selectors

Most tabs offer a choice between all genes in the annotation and a list you paste. Choosing the second reveals a gene box, with **Select all** and **Clear all** buttons. Where two datasets are loaded, a stage selector controls whether the figure shows Stage A, Stage B or both.

---

## Metagene Profile

Averages signal for each mark across selected genes, producing the classic profile around a reference point.

| Control | Purpose |
|---|---|
| Reference (one or more) | TSS, TES or gene body. Selecting several produces a panel each |
| Upstream / Downstream (bp) | Window either side of the reference point |
| Number of bins | How finely the window is divided. More bins give more detail and more noise |
| Gene scope | All genes, or only the genes you paste |
| Marks to display | Which marks to draw |
| Condition comparison | Side-by-side puts each condition in its own panel; Overlay draws them on shared axes |
| Smooth curve (LOESS) | Smooths the profile. Recommended when many genes are included |

The profile follows each gene's strand. For minus-strand genes the window is read from the 5′ end, so the TSS is always at the TSS end of the plot.

**Buttons:** *Download Plot* (PDF) · *Download Data (Excel)* (per-bin values behind it).

---

## Gene Contributions

A metagene line is an average, and averages hide their composition. This tab breaks one profile back into the genes behind it.

Choose one stage, mark and condition, then the same reference and window controls as the metagene tab. **Genes to show in heatmap** keeps the heatmap legible. **Log-scale the heatmap colours** compresses the range when a few genes dominate.

**Outputs:**

- A per-gene heatmap, one row per gene, sorted by total signal in the window. The column-wise average of this matrix is the metagene line from the previous tab.
- A cumulative contribution curve. A steep rise means a few genes carry the profile. A diagonal means every gene contributes roughly equally. A Gini coefficient and the top-1%, top-10% and top-50% shares are printed beneath it.
- A gene ranking table with coordinates, mean and peak signal, each gene's share of the total, and a running cumulative share.

**Buttons:** *Download Heatmap (PDF)* · *Download Curve (PDF)* · *Download Gene Table (Excel)*. The Excel file has two sheets: the ranking, and the per-gene matrix behind the heatmap.

---

## Expression Deciles

Splits genes into ten groups by mean normalised expression, then plots the gene-body profile of each mark for each decile. It shows how signal along the gene changes with expression level.

Requires the sample sheet and the count table for each stage used.

| Control | Purpose |
|---|---|
| Stages | Which loaded stages to include |
| Rank genes using these groups | Which colData groups define expression. Only samples present in the count table are used |
| Exclude genes with mean count below | Expression floor applied before ranking into deciles |
| Marks to display | Which marks to draw |
| Conditions to show | Which conditions to draw |
| Display | Separate panels (one per condition), or overlay (conditions in one panel, solid and dashed) |
| Upstream / Downstream flank (bp) | Size of the flanking windows |
| Gene body bins / Bins per flank | Resolution of the scaled gene body and the flanks |
| **Compute Profile** | Runs the analysis |

The gene-body window is scaled: flanks and body are each split into bins, so genes of different lengths contribute to the same positions. The caption lists the number of genes in each decile.

**Buttons:** *Download Plot (PNG)*, at the width and height set in the size panel.

---

## Enrichment Profile

Plots `log(observed/expected)` enrichment across a composite gene: upstream flank, gene body scaled 0–100%, downstream flank. Expected is the genome-wide mean signal for that mark and condition, so zero means "as expected by chance" and the sign gives enrichment or depletion.

| Control | Purpose |
|---|---|
| Upstream of TSS / Downstream of TES (bp) | Size of the two flanks |
| Bins per region | Resolution of the composite profile |
| Gene scope | All genes, or the list you paste |
| Marks to display / Stages to compare | Which marks and datasets to include |
| Panel layout | One panel per mark (conditions overlaid), or one panel per condition (marks overlaid, as in the Galaxy chromstaR output) |
| Ratio calculation | Ratio of means (recommended) or legacy per-bin log ratio |
| Smooth curve (LOESS) | Smooths each curve |

### Why the ratio calculation matters

**Ratio of means** averages the raw signal at each position, zero bins included, then takes a single log ratio against the genome-wide mean. This is what the Galaxy chromstaR profile does, and it is the correct default.

**Legacy** takes the log of each bin and averages the logs. Empty bins must be discarded, which biases every position upward, most severely where a mark is genuinely absent. Real depletion is flattened away. The option exists only to reproduce older figures.

### Reading the gene count in the caption

The subtitle reports the total number of distinct genes contributing anywhere in the profile, and the range of genes contributing per position. These differ because a gene is counted at a position only if one of its bins falls there. **Quote the total as your sample size**. The range tells a reviewer how evenly the profile is supported.

**Buttons:** *Download Plot* (PDF) · *Download Data (Excel)* (per-position log ratios, gene counts and contributing gene names).

---

## Bivalent

Treats two or more marks as one combination and compares it with each mark on its own, around the TSS, gene body and TES. Use it to ask whether a bivalent-like state (for example H3K4me3 with H3K27me3) is present, where, and in which conditions.

| Control | Purpose |
|---|---|
| Which combination | Tick the marks yourself, or pick a state that occurs in your data, ranked by bin count |
| Co-occurrence | Each line is the fraction of bins in the chromstaR state, as log(observed/expected). This mode shows real co-occurrence |
| Sum | Adds the RPKM of the chosen marks, either normalised to each mark's genome-wide mean or raw. A high sum can come from one mark alone, so it is not evidence of co-occurrence |
| Contains all / Exactly | Whether other marks may be present in the state, or the chosen marks must be the only ones |
| Gene call rule | A gene is called when any promoter bin is in the state, or when at least half of them are |
| Promoter window | Upstream and downstream of the TSS used for gene calls |
| Y axis | log(observed/expected), or the observed level. The observed axis shows losses directly |

**Outputs:** the profile figure, a summary of genes carrying the combination in each condition with gained, lost and kept counts between consecutive conditions in your life-cycle order, and a per-gene table you can filter to genes that gain, lose or carry the combination.

**Buttons:** *Download Plot (PDF)* · *Download Data (Excel)*, with three sheets: settings, profile and genes.

---

## Region Browser

Shows raw signal at a locus rather than an average. Use it to check whether a genome-wide pattern is real at a place you know.

| Control | Purpose |
|---|---|
| Chromosome scope | A single chromosome zoomed to a range, selected chromosomes at full length, or all chromosomes |
| Start / End (bp) | The visible window when a single chromosome is selected |
| Jump to gene + **Go** | Centres the view on a named gene plus a flanking distance |
| Flanking (bp) | How much context to show either side when jumping |
| Bin scope | All bins including intergenic, only bins inside a gene, or only intergenic bins. **The only place in the app that shows truly intergenic positions** |
| Marks to display / Stages to compare | Which marks and datasets to draw |
| Condition comparison | Side-by-side panels or overlaid conditions |
| **Load Region** | Runs the query. Nothing is drawn until this is pressed |

**Buttons:** *Download Plot* (PDF).

---

## Differential Peaks

The quantitative core: which marks change between conditions, in which direction, and in how many regions. Reproduces the Galaxy chromstaR differential analysis and extends it.

| Control | Purpose |
|---|---|
| Min differential score | chromstaR's confidence threshold. Galaxy default `0.9999` |
| Min merged region width (bp) | Discards very short domains. Galaxy default `300` |
| Gene scope | Whole genome, or only segments overlapping the genes you paste |
| Conditions to compare | 2 to 5 conditions, in your life-cycle order |
| Comparison mode | All pairwise combinations, or one reference against each of the others |
| Ranked figure style | Stacked (gained + lost) or Diverging (lost left of zero, gained right) |
| Genomic zones | Optionally restrict counts to chosen zones (needs the gene annotation) |
| **Compute Differential Peaks** | Runs the analysis |

### How the filtering works, and the two numbers

Filtering happens in two steps:

1. Segments are kept if their differential score **and** width pass the thresholds.
2. Of those, only segments whose chromatin state actually **differs** between the two conditions are retained. These are the *differential regions* everything on the tab is built from.

The two counts are usually almost identical but need not be. A segment can accumulate enough total difference to pass the score threshold without any single mark crossing the boundary that defines its state. Such a segment is correctly excluded. **The figure caption reports the number of differential regions, the number the bars are built from.**

### Outputs

- **Differential Peaks per Histone Mark**: grouped bars, one pair per mark, giving the number of regions where the mark is present in one condition and absent in the other, in each direction.
- **Per-mark summary table**: Gained, Lost and Total changed for each mark, each with its own percentage, plus the direction label.
- **Differential regions by genomic zone**: counts per zone (TSS, gene body, upstream, TES, intergenic).
- **Marks ranked by change, differential regions**: marks ranked by how many differential regions they change in, each bar split into gained and lost.
- **Marks ranked by change, genome-wide**: the same ranking across every domain in the genome, with no score or width filter, from the object's frequency table. Percentages here are of all domains and are **not** comparable with the filtered figure. This panel needs the `frequencies` slot.

### Two things to be careful about when quoting these numbers

> [!CAUTION]
> **The total is not a directional count.** "Total changed" for a mark is gained **plus** lost. The Direction label only says which is larger. If a mark shows 30,284 total changed regions labelled "Gained in Meta", that does **not** mean 30,284 regions gained the mark. Quote the Gained or Lost column when you mean a direction.

> [!CAUTION]
> **The percentages overlap and will not sum to 100%.** Most differential regions have several marks changing at once, and each region is counted once for every mark that changes in it. The note above the summary table reports the mean number of marks changing per region for your data.

**Buttons:** *Download Plot* (PDF, and Excel) · *Download Data (Excel)* · *Download Table (Excel)* · *Download all differential regions (CSV)*, one row per region, with coordinates, the state in each condition, the marks that changed, and the gene and zone when annotation is loaded.

---

## Gene Set Comparison

Tests whether a gene list of your own carries different promoter chromatin from a background. Use it when you arrive with a list from a differential expression analysis.

| Control | Purpose |
|---|---|
| Upstream / Downstream of TSS (bp) | Defines the promoter window scored for each gene |
| Stages to compare | Each stage runs independently and is tagged in the results |
| Replicates | Average across replicates, or keep each separate to check consistency |
| Summary statistic | Mean or median. Drives delta and the permutation test's central value. Both are always shown |
| Comparison type | (1) permutation against random genes from the whole annotation; (2) permutation against random draws from your own background list B; (3) a direct Wilcoxon test of A against B |
| Remove set A's genes from list B | Option for type 2, on by default |
| Gene set A / Gene set B | Paste gene identifiers, one per line. A match summary reports how many were found |
| Number of random draws | Permutation modes only |
| **Compare Gene Sets** | Runs the analysis |

**Outputs:** a statistics table ranked by FDR, a box plot of the distributions, and a per-gene table of promoter posterior probabilities for set A. In permutation mode the empirical p-value is the fraction of random draws as extreme as your real set.

**Background list B.** When the background is a matched set (for example all genes with expression similar to A), choose the second permutation option. Every random draw is then taken only from B, so expression level is controlled for. B must be larger than A, ideally several times larger. A warning appears when B is smaller than three times A.

> [!NOTE]
> Always read the match summary before interpreting anything. If only part of your identifiers matched the annotation, the test ran on that part. Identifier-format mismatches are the usual cause.

**Buttons:** *Download Plot* (PDF) · *Download Data (Excel)*. In permutation mode the Excel file has the settings, statistics, per-gene posteriors, permutation draws and set sizes.

---

## Mark Combinations

Counts genes by the **exact** set of marks they carry, in the UpSet style. Each gene is counted once, in the combination of exactly its marks.

| Control | Purpose |
|---|---|
| Stage / Condition | One stage and one condition at a time |
| Marks | Which marks to consider |
| Gene window | Promoter (TSS ± flanks), or the whole gene body |
| Upstream / Downstream flank (bp) | Size of the promoter window |
| Show expression boxes | Adds a top panel with the expression of the genes in each combination |
| Expression from these groups | Which colData groups define expression for those boxes |
| Combinations to show / Minimum genes per combination | Limits what is drawn |
| **Compute Combinations** | Runs the analysis |

**Outputs:** the intersection bars (number of genes per combination), the dot matrix showing which marks make up each combination, set sizes per mark (genes carrying the mark in any combination), and, optionally, expression per combination.

**Buttons:** *Download Plot (PNG)*, at the width and height set in the size panel. A text summary below the figure reports the number of genes with at least one mark and the set size of each mark.

---

## RF Region Model

These tabs ask a different question from the rest of the app: given the histone marks at eight positions along a gene, how well do they predict that gene's expression, and which mark at which position pushes expression up or down?

Each gene is described by **64 features**: eight marks × eight regions. The regions are:

| Region | Window relative to the TSS (gene orientation) |
|---|---|
| TSS −2 kb | −2000 to −1000 bp |
| TSS −1 kb | −1000 to 0 bp |
| TSS ± 200 bp | −200 to +200 bp |
| Body 5′ | First third of the gene body |
| Body mid | Middle third |
| Body 3′ | Last third |
| TES + 1 kb | 1000 bp past the 3′ end |
| TES + 2 kb | 1000 to 2000 bp past the 3′ end |

The eight marks are H3K23ac, H3K27me3, H3K36me3, H3K4me3, H3K9ac, H3K9me3, H4K12ac and H4K20me1. **Every stage compared must contain all eight marks.** The app stops with a message naming the missing marks if one does not.

Each region's signal is the mean RPKM of the bins under it, weighted by how many bases of the window each bin covers. A region is reported as missing (and the gene is dropped from the model) when less than half of its window is covered by bins.

The target is log₂(expression + 1), where expression is the mean over the selected samples of the count table. A random forest predicts the target from the 64 features, and TreeSHAP attributes each prediction to individual features.

### Requirements

The tabs use the objects already loaded on the Load Data tab: the ChromstaR object, the gene annotation, the sample sheet and the normalised count table for each stage. Nothing else needs to be uploaded. The count-table gene identifiers must match `gene_id` or `gene_name` in the annotation.

### Genes excluded by default

**Genes whose gene bodies overlap another gene's body are excluded** by default. Their region windows fall on the same bins, so the features are not independent. Flank overlap between neighbouring genes is common and is not used as a criterion. The option can be turned off.

Genes are also removed in two cases:

- genes with no bins under any region, and
- *outliers*: genes with zero signal for every feature but log₂(expression + 1) above 1.

The model needs at least 50 genes left after these steps. The numbers removed are reported in the performance table.

### RF Region Model tab

| Control | Purpose |
|---|---|
| Condition 1 / Condition 2 | Any two stage-and-condition pairs, taken from the same list as the rest of the app. Samples come from the sample sheet |
| Trees (combined / single region) | Forest size for the 64-feature model and for each single-region comparison |
| CV folds | Folds for the stability check on the training set |
| Cores for SHAP | Parallel workers for TreeSHAP (one on Windows) |
| Test genes to explain | 0 explains every held-out gene. A value such as 2000 gives a quick first run |
| Exclude genes whose gene bodies overlap another gene | On by default (see above) |
| **Run RF + SHAP** | Runs the pipeline. Nothing is computed until this is pressed |

The model is split **80% training / 20% held-out test** (caret partition, seeded).

### Saved runs

Every completed run is saved automatically to a folder (default `rf_shap_saved`, editable). Each run writes two files: an `.RData` file holding the results and a small `.meta.rds` file describing it.

- **If a saved run with the same settings exists**, pressing *Run RF + SHAP* loads it instead of recomputing, when *load it instead of recomputing* is ticked. A run matches only when the stages, conditions, samples, marks, number of bins and segments in each object, and every setting that changes the result are the same.
- **Load selected saved run** loads any run listed in the folder, newest first.
- **Refresh list** re-reads the folder.

### Outputs

- **Performance table**: correlation and RMSE on held-out genes, the number of training and test genes, outliers removed, genes dropped for lack of bins, and feature and mark counts.
- **Cross-validation table**: mean and SD of r, RMSE and R² over the folds.
- **Single vs all regions**: the correlation of each region's model alone, compared with all regions combined.
- **Feature correlation matrix**: Spearman correlations between the 64 features on the training set.

Correlated features share importance, so read the correlation matrix before interpreting SHAP values.

### Downloading everything

**Download all outputs (zip)** writes:

| File | Contents |
|---|---|
| `performance_combined_model.csv` | Held-out performance |
| `cv_<k>fold_results.csv` | Per-fold results |
| `genes_sharing_bins.csv` | Genes excluded for overlapping gene bodies (only if any) |
| `single_region_performance.csv` | Single-region correlations |
| `shap_feature_summary.csv` | Per-feature importance and direction |
| `feature_correlations.csv` | Full feature correlation table |
| `cross_stage_performance.csv` | Cross-stage prediction (see RF Combined) |
| `feature_correlated_pairs.csv` | Pairs with \|rho\| at or above the threshold |
| `cv_summary.csv` | Mean and SD across folds |
| `incmse_per_feature.csv` | Classic random forest importance |
| `mark_importance_summary.csv` | Importance summed over regions, per mark |
| `region_importance_summary.csv` | Importance summed over marks, per region |
| `methods_text.txt` | Draft methods paragraph filled from the run |
| `rf_shap_fits.rds` | The fitted models and their data |
| `Fig_*.pdf` | Correlation matrix, single vs combined, mark × region bubble map, expressed vs silent heatmap, one beeswarm per model, mark and region importance, %IncMSE, and stage comparison (when there are two stages) |

> [!NOTE]
> The analysis is a correlation and attribution study, not a causal one. A high SHAP value means the model relies on that feature to predict expression. It does not show that the mark regulates the gene.

---

## SHAP Mark x Region

Shows the attributions from the most recent run on the RF Region Model tab, or from a run loaded there.

- **Which mark matters where**: one bubble per mark and region. Size is importance (mean absolute SHAP). Colour is direction (Spearman correlation between the feature's value and its SHAP value).
- **Expressed vs silent genes**: average push of each mark at each position, separately for genes with expression at or above 1 (log₂ scale, as in FPKM ≥ 1) and below it.
- **Top features (beeswarm)**: one panel per model, showing the twenty features with the largest mean absolute SHAP. Each dot is a gene. Right of zero pushes predicted expression up. Red dots are genes where the mark is present, blue where it is absent.

---

## RF Combined (Rmd)

Shows the same run as the RF Region Model tab, arranged as the analyses of the reference Rmd script.

| Control | Purpose |
|---|---|
| Report feature pairs with \|rho\| at least | Threshold for the correlated-pairs table |

Outputs:

- **Cross-stage prediction**: each stage's model applied to its own held-out test genes, and to all genes of the other stage. Shows whether a model transfers.
- **Cross-validation summary**: mean and SD across folds, per stage.
- **Strongly correlated feature pairs**: a searchable table.
- **Which mark matters (summed over regions)**: SHAP importance per mark, with direction.
- **Classic random forest importance (%IncMSE)**: unscaled permutation importance summed over regions. It has no direction.
- **Importance per region**: SHAP importance summed over marks, per region.
- **Mark @ region importance, stage A vs stage B**: a scatter comparing each feature's importance across the two stages, with the Spearman correlation and the top ten features labelled. Needs two stages.
- **Methods (draft)**: a paragraph describing the analysis, filled in from the run's settings and counts. Check it against what you did before using it in a manuscript.

---

## SHAP Gene Explorer

Explains one gene's predicted value.

Choose a model (stage) and enter one or more gene identifiers. The first one found in that model's test set is drawn as a waterfall. It starts at the average predicted expression (the baseline) and shows how each of the twelve largest feature contributions, with the rest grouped as "other features", moves the prediction to the model's value for that gene. The subtitle reports the baseline, the predicted value and the measured value.

---

## Data Table

The per-bin values behind every other tab, as a searchable, sortable table. Use it to check a specific position or confirm a column exists under the name you expect.

Each bin is annotated with the gene and genomic zone it falls into, and with whether it overlaps a differential region when the Differential Peaks tab has been run.

> [!NOTE]
> Use **Download full table (CSV)** for a full export. The table's own export buttons handle only the displayed page reliably, which for several hundred thousand rows is almost never what you want.

---

## Exports

| Tab | Button | Contents |
|---|---|---|
| Metagene | Download Plot / Data (Excel) | PDF of the profile; per-bin mean values behind it |
| Gene Contributions | Download Heatmap / Curve (PDF) | Per-gene heatmap; cumulative curve |
| Gene Contributions | Download Gene Table (Excel) | Ranking and per-gene matrix |
| Expression Deciles | Download Plot (PNG) | Figure at the size set in the panel |
| Enrichment | Download Plot / Data (Excel) | PDF of the profile; per-position log ratios, gene counts, contributing genes |
| Bivalent | Download Plot (PDF) / Data (Excel) | PDF of the profile; settings, profile and gene tables |
| Region Browser | Download Plot | PDF of the current view |
| Differential Peaks | Download Plot / Data (Excel) | PDF of the per-mark bar chart; per-mark, per-direction counts |
| Differential Peaks | Download Table (Excel) | Per-mark gained, lost, totals and percentages |
| Differential Peaks | Download all differential regions (CSV) | One row per differential region |
| Differential Peaks | Download Plot (PDF) ×2 | Filtered and genome-wide ranked figures |
| Gene Set Comparison | Download Plot / Data (Excel) | PDF of the box plot; statistics, per-gene posteriors and settings |
| Mark Combinations | Download Plot (PNG) | UpSet-style figure |
| RF Region Model | Download all outputs (zip) | CSV tables, fitted models (RDS), methods text and PDF figures |
| Data Table | Download full table (CSV) | Every row of the per-bin table |

PDF exports honour the current Plot Appearance settings.

---

## Troubleshooting

### The app stops at startup, saying a file cannot be found
Both `rf_shap_tools.R` and `rf_shap_extra.R` must be in the same folder as `app_v29.R`, and the app must be started from that folder. Copy all three files together.

### A figure fails with an error mentioning a missing argument
If a plot area shows something like `argument "observed" is missing, with no default`, another package attached in your R session is masking a function the app uses. The known case is **randomForest**, which exports its own `margin()`.

The app guards against this case and reports the resolution in the startup banner. If the banner says `margin()` resolves to anything other than `ggplot2`, restart R, and load fewer packages before starting the app. The RF + SHAP engine calls randomForest, caret and treeshap with `::` and attaches none of them, so it does not cause this. For any other masked function, the error message names the failing call and the package it came from, and a full traceback is printed to the console.

### The RF tab says a stage has no counts.rpkm columns for some marks
The RF Region Model needs all eight marks in every stage it compares. The message names the missing marks. Either the object lacks those marks, or the mark names are spelled differently from the rest of the object.

### The RF tab says no genes are shared
The gene IDs in the count table do not match `gene_id` or `gene_name` in the annotation. Check the row names of the count table against the annotation's identifiers.

### The RF tab says fewer than 50 genes are left
Usually the count table or the annotation is wrong, or the expression floor is removing most genes. Check the count-table summary on the Load Data tab, and the chromosome naming in the annotation.

### The genome-wide ranked figure shows a note instead of a plot
The object was loaded by an older build that discarded the `frequencies` slot. Re-upload the ChromstaR object. No other change is needed.

### A plot does not change when I move a slider
Expected behaviour. Press the tab's Compute button. See [Controls that appear on every tab](#controls-that-appear-on-every-tab).

### Fewer genes are reported than I loaded
Usually chromosome naming. Genes on sequences absent from the ChromstaR bins, such as unplaced scaffolds, never overlap a bin and contribute nothing. Check that naming matches between the annotation and the object. The test data deliberately includes unplaced contigs to exercise this.

### Profiles look smeared or symmetric when they should not
Check that your gene table has a `strand` column. Without it every gene is treated as unstranded and minus-strand genes are read backwards.

### The app is slow
Whole-genome operations scale with the number of genes and bins. Narrow the gene scope, reduce the number of bins, or restrict the marks and conditions. The RF + SHAP tabs are the slowest. Use a smaller "Test genes to explain" value for the first run, and reuse saved runs where possible.

---

## Author

**Janan Gawra**, IHPE (UMR 5244, CNRS, Université de Perpignan Via Domitia)

Application build v29.

Built for chromatin and epigenomics analysis by Janan Gawra.
[linkedin.com/in/janangawra](https://linkedin.com/in/janangawra)
