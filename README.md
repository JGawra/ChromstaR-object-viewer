# About ChromstaR Viewer

## What this tool does

ChromstaR Viewer is a point-and-click app for exploring chromatin state and histone
modification signal from a ChromstaR `combinedMultiHMM` object, alongside a gene
annotation. It lets you compare experimental conditions or developmental stages
across multiple histone marks, browse specific genomic regions, quantify
differential chromatin states, and compare gene sets — all without writing any R
code.

Use the sidebar on the left to move between tabs. Each tab is independent, so you
can jump around in any order once data is loaded. Every tab has its own **Compute**
button: changing a setting does nothing until you click it.

---

## 1. Load Data tab

**Purpose:** upload your data. Nothing else works until this is done.

- **"Label for this dataset (used in plots)"** — the name this object gets in plot
  legends, facet strips and exported tables. Defaults to "Stage A".
- **"Upload ChromstaR .RData or .rds file" (Stage A)** — your ChromstaR
  `combinedMultiHMM` object. A summary appears once loaded (object class, bin
  count, and so on).
- **Stage B (optional)** — a second, independent ChromstaR object, e.g. a
  different developmental stage or a separate experiment. It does **not** need to
  share bin coordinates with Stage A: every analysis runs on each object
  separately and the results are combined and labelled by stage. Once Stage B is
  loaded, every analysis tab gains a **"Stages to compare"** control.
- **"Upload genes.tsv file"** — the gene annotation, shared by both stages. A
  plain table with columns `chr, start, end, gene_id, gene_name, strand`. Convert
  a GTF to this format first. Coordinates are read as 1-based inclusive (GTF
  convention); a BED file's 0-based starts will be off by one.
- **Detected Marks & Conditions panel** — fills in automatically and lists every
  histone mark, every condition, and the genome info found in each object. Check
  this before moving on: the chromosome names shown here must match the names in
  your gene table, or nothing will be annotated as genic.

### Condition / Life-cycle Order

The box at the bottom of the Load tab controls **both the order and the membership**
of conditions across the entire app.

- **Drag** the boxes to arrange conditions in a biologically meaningful order
  (e.g. `TwoDO → Meta → PNA → PA`, or `cerca → somula → mira → sp1 → adult`)
  instead of the default alphabetical order. That order is then applied
  everywhere at once: legend order, facet strips, the order of `rpkm_*` columns
  in the exported table, and which colour each condition gets.
- **Remove** a condition (the × on its box) and it is hidden from every plot,
  table and analysis — it is excluded from the computation itself, not filtered
  out afterwards.
- **"Reset to detected order"** restores every condition found in the loaded
  objects.
- A coloured preview strip below the box shows the current order and each
  condition's colour.

---

## 2. Metagene Profile tab

**Purpose:** plot average signal (RPKM) around a reference point — TSS, TES, or
across the gene body — for a chosen set of genes, one panel per histone mark.

- **"Reference" checkboxes** — tick TSS, TES and/or Gene body. More than one can
  be ticked; panels are laid out in biological order, **TSS → Gene body → TES**.
- **"Upstream (bp)" / "Downstream (bp)"** — shown when TSS or TES is ticked; sets
  the window before and after the reference point. A live summary line shows the
  exact window produced.
- **"Number of bins"** — how many segments the window or gene body is divided
  into. Choose this relative to your ChromstaR bin size: for a 4,000 bp window
  built from 250 bp bins, roughly 16–40 bins is right. Asking for many more bins
  than the data supports (e.g. 100 bins over 4 kb of 250 bp bins) produces a
  jagged line made of gaps, not finer resolution.
- **"Gene scope"** — "Only selected genes" uses the pasted list below; "All genes
  in GTF" ignores it and runs on every gene in the annotation.
- **"Select all" / "Clear all"** and the **gene textbox** — paste gene IDs, one
  per line (e.g. straight from an Excel column). A match count appears underneath
  showing how many were recognised and, if any failed, examples of the misses.
- **"Marks to display"** — which histone marks to include.
- **"Condition comparison"** — "Side-by-side" gives each condition its own facet;
  "Overlay" draws all conditions on the same axes in different colours.
- **"Smooth curve (LOESS)"** — draws a smoothed line instead of the raw signal;
  useful with many genes, but check the raw line too before interpreting a bump.
- **"Compute Profile"**, **"Download Plot"** (PDF), **"Download Data (Excel)"**
  (one row per bin/mark/condition).

**Strand.** Profiles are strand-aware: the window is placed relative to each
gene's own orientation and minus-strand genes are flipped before averaging, so
position 0% is always the 5′ end. Without that flip a promoter mark appears at
*both* ends of the gene body.

---

## 3. Enrichment Profile tab

**Purpose:** plot log(observed/expected) enrichment around gene boundaries — the
upstream flank, the gene body as a percentage of gene length, and the downstream
flank on one continuous axis.

- **"Upstream of TSS (bp)" / "Downstream of TES (bp)"** — how far before the start
  and after the end of each gene to include.
- **"Bins per region"** — resolution of the curve.
- **"Gene scope"**, **"Select all"/"Clear all"**, gene textbox, **"Marks to
  display"** — same behaviour as the Metagene tab.
- **"Panel layout"**
  - *One panel per mark (colour = condition)* — conditions overlaid inside each
    mark's panel. Best for asking how a single mark changes across the life cycle.
  - *One panel per condition (colour = mark)* — all marks overlaid inside each
    condition's panel, the layout used by the Galaxy chromstaR output. Best for
    asking which marks dominate at a given stage. Panels share one y-axis so
    conditions can be compared directly.
- **"Expected (baseline)"** — the denominator of log(observed/expected).
  - *Mean over covered bins (signal > 0)* — the default. The observed side has to
    drop zero-signal bins because log(0) is undefined, so this keeps numerator and
    denominator on the same set of bins.
  - *Mean over all genomic bins* — divides by a mean that includes every empty
    bin. This shifts each curve upward by log(mean covered / mean all); because
    coverage sparsity differs per mark and per condition, the shift differs per
    curve and conditions can no longer be compared by vertical position.
- **"Smooth curve (LOESS)"**, **"Compute Enrichment"**, **"Download Plot"**,
  **"Download Data (Excel)"**.

---

## 4. Region Browser tab

**Purpose:** a genome-browser-style view of raw signal across a chosen chromosome
and position range, or across whole chromosomes at once.

- **"Chromosome scope"** — "Single chromosome" to zoom into a bp range (the usual
  case); "Selected chromosomes" or "All chromosomes" to view entire chromosomes,
  useful as an overview but slower and more cluttered.
- **Chromosome dropdown** — which chromosome(s) to view (multi-select in
  "Selected chromosomes" mode).
- **"Start (bp)" / "End (bp)"** — the window, in single-chromosome mode.
- **"Jump to gene"** plus **"Flanking (bp)"** and **"Go"** — pick a gene by name,
  set the flanking sequence, and the Start/End boxes jump to that locus.
- **"Bin scope"** — all bins including intergenic, only bins inside a gene, or
  only intergenic bins, using the loaded annotation.
- **"Marks to display"**, **"Condition comparison"** — as in the Metagene tab.
- **"Load Region"**, **"Download Plot"** (PDF).

---

## 5. Differential Peaks tab

**Purpose:** counts, per histone mark, how many chromatin segments are confidently
present in one condition but not another (and vice versa) — a bar chart of which
marks change most between conditions. Any number of conditions is supported.

- **"Min differential score"** — keep only segments with a confidence score at or
  above this threshold (closer to 1 = stricter).
- **"Min merged region width (bp)"** — discard segments shorter than this.
- **"Gene scope"**
  - *Whole genome (all segments)* — the default, and what the Galaxy differential
    tool does.
  - *Only selected genes* — keeps only segments overlapping the genes you paste
    in, with the same Select all / Clear all buttons and live match counter as the
    other tabs. Use it to ask whether a specific gene set is remodelled between
    stages.
- **"Conditions to compare"** — pick 2 to 5 conditions from the loaded object.
  Drag to reorder, × to remove.
- **"Comparison mode"**
  - *All pairwise combinations* — one panel per pair (5 conditions → 10 panels).
  - *One reference vs the others* — every chosen condition against a reference you
    pick (5 conditions → 4 panels), usually what you want for a life-cycle
    baseline.
- **"Stages to compare"** — appears once Stage B is loaded; each stage is filtered
  and counted independently and gets its own row of panels.
- **"Compute Differential Peaks"**, **"Download Plot"**, **"Download Data
  (Excel)"**.

**How it works.** Segments are filtered once on `differential.score` and width.
Then, for each pair of conditions, only the segments whose combinatorial state
actually differs between those two conditions are kept, and each mark is counted
in both directions. Bars are coloured by the condition the mark is present in,
using the app-wide condition palette, and the plot canvas grows with the number of
panels.

---

## 6. Gene Set Comparison tab

**Purpose:** compares the average histone-mark posterior probability over the
promoter region of one gene list against either a second gene list or a random
background — for asking "is this curated gene set unusual for a given mark?" or
"do up- and down-regulated genes differ in chromatin state?"

- **"Upstream of TSS (bp)" / "Downstream of TSS (bp)"** — the promoter window
  averaged over for every gene.
- **"Replicates"** — "Average across replicates" merges rep1/rep2/… into one
  number per mark/condition; "Show each replicate separately" keeps them apart so
  you can check replicate consistency.
- **"Summary statistic"** — Mean or Median as the headline number (both are always
  shown in the results table).
- **"Comparison type"** — "Compare to a random background" runs a permutation test
  against many random gene sets of the same size (the statistically rigorous
  default); "Compare to a second gene list I provide" lets you paste an actual
  second list.
- **Gene set A textbox** — your main gene list, one ID per line.
- **"Number of random draws"** — random-background mode only. Higher gives more
  precise empirical p-values but is slower; 1,000 is a good working default,
  several thousand for a final figure.
- **Gene set B textbox** — manual mode only.
- **"Compare Gene Sets"** — fills in the results table, plot and per-gene table.
- **Results table** — one row per mark/condition, sorted by FDR, colour-coded by
  significance.
- **Boxplot** — the spread of individual gene values for set A vs set B (or vs one
  representative random draw), one panel per mark.
- **Per-gene table** — one row per gene in set A with its posterior value per
  mark/condition; sortable, filterable, exportable.
- **"Download Plot"** (PDF), **"Download Data (Excel)"** (multi-sheet).

---

## 7. Data Table tab

**Purpose:** the full per-bin chromatin state table, annotated with which gene (if
any) and which genomic zone each bin falls into — TSS ±200 bp, gene body thirds,
upstream/downstream flanks, or intergenic.

- **"Download full table (CSV)"** — exports every row. Use this rather than the
  table's own export, which only reliably handles small subsets.
- **Table itself** — scrollable, sortable, with per-column search boxes.
- `rpkm_<mark>_<condition>` columns are emitted in your chosen life-cycle order,
  and only for the conditions you kept in the Condition / Life-cycle Order box.

**Check the `genomic_zone` column.** If every bin says "Intergenic", the gene
annotation matched nothing — almost always a chromosome-naming mismatch between
the ChromstaR object and the gene table. The app now warns you when this happens
and prints both sets of names.

---

## Key concepts that apply across several tabs

- **Gene scope** (Metagene, Enrichment, Differential Peaks): "only selected genes"
  vs everything. On the Metagene and Enrichment tabs, both modes only ever include
  bins overlapping a gene — there is no "intergenic" position relative to a TSS.
- **Bin / chromosome scope** (Region Browser): every bin, only genic, or only
  intergenic; one chromosome zoomed in, or several at full length.
- **RPKM vs log(observed/expected)** — Metagene and Region Browser show raw mean
  RPKM. Enrichment Profile shows a log-ratio against a genome-wide baseline,
  better for comparing marks with very different baseline signal levels.
- **Posterior probability** (Gene Set Comparison) — ChromstaR's 0–1 confidence
  that a mark is genuinely present at a bin, distinct from the RPKM signal used
  elsewhere.
- **Stages vs conditions** — *conditions* live inside one ChromstaR object
  (jointly modelled by the HMM). *Stages* are two separate objects loaded side by
  side, computed independently and combined afterwards.

---

## Checking that a result is real

- **Set the bin count relative to your bin size.** More plotted points than the
  data supports gives noise, not resolution.
- **Split genes by strand** and run the metagene on `+` genes only, then `−` genes
  only. The two profiles must look the same. Mirror images mean orientation is
  broken.
- **Shift the annotation** by 50 kb and re-run: every profile should flatten.
  A peak that survives is an artifact of window geometry, not of the TSS.
- **Check a housekeeping gene** in the Region Browser — H3K4me3 should be a sharp
  peak at its TSS.
- **Check the mitochondrial and unplaced contigs** — they should carry essentially
  no real ChIP enrichment.
- **Check coverage per mark and condition** — the fraction of zero-signal bins. A
  mark with very few covered bins produces a metagene driven by a handful of
  regions.
- **Check replicate concordance** in Gene Set Comparison with "Show each replicate
  separately".
- **Check direction balance** in Differential Peaks — if nearly every change is
  gained in one condition for every mark, suspect sequencing depth or
  normalisation rather than biology.
- **Vary the thresholds.** Re-run at differential score 0.99 / 0.999 / 0.9999 and
  several widths; the ranking of marks should be stable.

---

## Performance tips

- Computations are vectorised and stay fast with all genes selected, but very
  large gene sets combined with many marks can take a few seconds.
- The Data Table's gene/zone annotation is computed once per session and cached;
  the first visit to that tab takes a little longer. Loading a new gene file
  clears the cache so the table is rebuilt against the new annotation.
- "All chromosomes" in the Region Browser, and thousands of permutations in Gene
  Set Comparison, are noticeably slower — preview with smaller values, then scale
  up for the final result.
- All-pairwise mode with 5 conditions draws 10 panels; the canvas grows to match,
  so expect a taller plot and a slightly longer render.

---

## Updates

### 09/09/2026

**Multi-condition support.** The app now works with ChromstaR objects containing
any number of conditions, not just two. Differential Peaks gained a condition
picker (2–5) and two comparison modes — all pairwise combinations, or one
reference vs the others.

**Condition / Life-cycle Order now filters as well as orders.** Removing a
condition hides it from every plot, table and analysis, and it is excluded from
the computation rather than filtered afterwards.

**Metagene profiles are strand-aware.** Minus-strand genes were being averaged in
backwards, which made promoter marks appear at both ends of the gene body and
symmetrised the TSS profile. Positions are now flipped for minus-strand genes, and
bins whose midpoint falls outside the window are dropped instead of being piled
into the first and last plot points. Reference panels are ordered TSS → Gene body
→ TES.

**Enrichment Profile: Galaxy-style layout.** A new panel layout puts one condition
per panel with all marks overlaid, alongside the existing one-panel-per-mark view.

**Enrichment Profile: corrected baseline.** The observed side drops zero-signal
bins, but the expected value was previously a mean over all bins including empty
ones, shifting each curve up by an amount that differed per mark and per
condition. The baseline now defaults to the mean over covered bins; the old
behaviour remains selectable.

**Differential Peaks: gene scope.** Restrict the analysis to a pasted gene list
instead of the whole genome.

**Bin annotation cache fixed.** The per-bin gene/zone annotation is now recomputed
when a new gene file is loaded, and the app warns when the annotation matches zero
bins (usually a chromosome-naming mismatch).

---

Built for chromatin and epigenomics analysis by Janan Gawra.
[linkedin.com/in/janangawra](https://linkedin.com/in/janangawra)
