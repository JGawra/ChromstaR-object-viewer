About ChromstaR Viewer

What this tool does
ChromstaR Viewer is a point-and-click app for exploring chromatin state and histone modification signal from a ChromstaR combinedMultiHMM object, alongside a gene annotation. It lets you compare experimental conditions (e.g. PA vs PNA) across multiple histone marks, browse specific genomic regions, and compare gene sets — all without writing any R code. Use the sidebar on the left to move between tabs; each tab is independent, so you can jump around in any order once data is loaded.

1. Load Data tab
Purpose: upload your data. Nothing else works until this is done.

"Upload ChromstaR .RData file" (left box) — click to browse for your ChromstaR combinedMultiHMM object (a .RData / .rda file). A text box below it shows a summary once loaded (object class, bin count, etc).
"Upload genes.tsv file" (right box) — click to browse for your gene annotation. This should be a simple table with columns chr, start, end, gene_id, gene_name, strand — convert a GTF to this format first if needed.
Detected Marks & Conditions panel at the bottom fills in automatically once both files are loaded — it lists every histone mark, every condition, and basic genome info found in your object. Use this to sanity-check the upload before moving to other tabs.
2. Metagene Profile tab
Purpose: plot average signal (RPKM) around a reference point (TSS, TES, or across the gene body) for a chosen set of genes, one panel per histone mark.

"Reference" checkboxes — tick TSS, TES, and/or Gene body. You can tick more than one to compare them side by side in the same plot.
"Upstream (bp)" / "Downstream (bp)" — only shown when TSS or TES is ticked; sets how far before/after the reference point to plot. A live summary line below shows exactly what window this produces.
"Number of bins" — how many segments to divide the window/gene body into; more bins = finer resolution but a noisier-looking line.
"Gene scope" — "Only selected genes" uses whatever you've pasted below; "All genes in GTF" ignores the pasted list and runs on every gene in the annotation.
"Select all" / "Clear all" buttons — fill or empty the gene textbox below with every gene name from the annotation.
Gene textbox — paste gene IDs here, one per line (e.g. copy-pasted from an Excel column). A match count appears underneath showing how many were recognised.
"Marks to display" checkboxes — tick which histone marks to include in the plot.
"Condition comparison" — "Side-by-side" puts each condition in its own facet panel; "Overlay" draws all conditions on the same axes with different colours.
"Smooth curve (LOESS)" checkbox — draws a smoothed line instead of the raw jagged signal; recommended when plotting many genes at once.
"Compute Profile" — runs the calculation and draws the plot. Nothing happens until you click this, even if you change settings above.
"Download Plot" — saves the current plot as a PDF.
"Download Data (Excel)" — saves the exact numbers behind the plot (one row per bin/mark/condition) as an .xlsx file.
3. Enrichment Profile tab
Purpose: plot log(observed/expected) enrichment of each mark around gene boundaries (TSS and TES together), one panel per mark — useful for seeing which marks are relatively enriched or depleted across a gene versus the genome-wide average.

"Upstream of TSS (bp)" / "Downstream of TES (bp)" — how far before the start and after the end of each gene to include.
"Bins per region" — resolution of the plot; more bins = finer detail.
"Gene scope", "Select all"/"Clear all", gene textbox, "Marks to display" — same behaviour as the equivalent controls in the Metagene Profile tab above.
"Smooth curve (LOESS)" checkbox — same as Metagene tab, smooths the line.
"Compute Enrichment" — runs the calculation and draws the plot.
"Download Plot" / "Download Data (Excel)" — same as Metagene tab: save the figure as PDF, or the underlying numbers as Excel.
4. Region Browser tab
Purpose: a genome-browser-style view of raw signal across a chosen chromosome and position range, or across whole chromosomes at once — for zooming into a specific locus.

"Chromosome scope" — "Single chromosome" lets you zoom into a specific bp range (most common use); "Selected chromosomes" or "All chromosomes" show entire chromosome(s) at full length instead, useful for a wide overview but can be slow/cluttered with many small contigs.
Chromosome dropdown — pick which chromosome(s) to view (becomes a multi-select box in "Selected chromosomes" mode).
"Start (bp)" / "End (bp)" — only shown in single-chromosome mode; defines the exact window to zoom into.
"Jump to gene" dropdown plus the "Flanking (bp)" box and the "Go" button — instead of typing coordinates by hand, pick a gene name here, set how much flanking sequence you want on each side, and click "Go" to jump the Start/End boxes straight to that gene's location.
"Bin scope" — "All bins, including intergenic" shows everything; "Only bins inside a gene" or "Only intergenic bins" filter the view using the loaded gene annotation.
"Marks to display" checkboxes — tick which histone marks to show as tracks.
"Condition comparison" — "Side-by-side" or "Overlay", same meaning as in the Metagene tab.
"Load Region" — fetches and draws the signal for the current settings.
"Download Plot" — saves the current view as a PDF.
5. Differential Peaks tab
Purpose: counts, per histone mark, how many chromatin segments are confidently present in one condition but not the other (and vice versa) — a bar chart summarising which marks change the most between conditions.

"Min differential score" — only keep chromatin segments with a confidence score at or above this threshold (closer to 1 = stricter, fewer but more confident segments).
"Min merged region width (bp)" — discard segments shorter than this, to avoid counting tiny noisy regions.
"Compute Differential Peaks" — runs the filtering and counting, then draws the bar chart.
"Download Plot" / "Download Data (Excel)" — save the figure as PDF, or the underlying per-mark counts as Excel.
6. Gene Set Comparison tab
Purpose: compares the average histone-mark "posterior probability" (confidence that a mark is present) over the promoter region of one gene list against either a second gene list you provide, or a random background — useful for asking "is this curated gene set unusual for a given mark?" or "do up- and down-regulated genes differ in chromatin state?"

"Upstream of TSS (bp)" / "Downstream of TSS (bp)" — defines the promoter window averaged over for every gene.
"Replicates" — "Average across replicates" merges rep1/rep2/etc into one number per mark/condition; "Show each replicate separately" keeps them apart so you can check replicate consistency.
"Summary statistic" — choose Mean or Median as the main number used to compare the two groups (both are always shown side by side in the results table regardless of this choice).
"Comparison type" — "Compare to a random background" runs a permutation test against many random gene sets of the same size (statistically rigorous, recommended default); "Compare to a second gene list I provide" lets you paste an actual second list instead (e.g. down-regulated genes).
Gene set A textbox — paste your main gene list here, one ID per line.
"Number of random draws" — only shown in random-background mode; how many random gene sets to sample for the statistical test. Higher = more precise p-values but slower (1000 is a good default; thousands for a final result).
Gene set B textbox — only shown in manual mode; paste your second gene list here.
"Compare Gene Sets" — runs the comparison and fills in the results table, plot, and per-gene table below.
Results table (top right) — one row per mark/condition, sorted by FDR (most significant first); colour-coded green/yellow by significance.
Boxplot — shows the spread of individual gene values for set A vs set B (or vs one representative random draw), one panel per mark.
"Download Plot" / "Download Data (Excel)" — save the boxplot as PDF, or the full results (stats + per-gene values) as a multi-sheet Excel file.
Per-gene table (bottom) — one row per gene in set A, showing its individual posterior value per mark/condition; sortable, filterable by column, and exportable directly as CSV/Excel using the buttons above the table.
7. Data Table tab
Purpose: the full per-bin chromatin state table, annotated with which gene (if any) and which genomic zone (TSS, gene body thirds, upstream/downstream flanks, or intergenic) each bin falls into.

"Download full table (CSV)" — exports every row (500k+) to CSV. Use this rather than the table's own export, which only reliably handles small subsets.
Table itself — scrollable and sortable by column; use the search boxes to filter by gene name, zone, or other fields.
Key concepts that apply across several tabs
Gene scope (Metagene/Enrichment tabs): "only selected genes" vs "all genes in GTF." Both modes only ever include bins overlapping a gene — there's no "intergenic" position relative to a TSS or gene body.
Bin/Chromosome scope (Region Browser tab): whether to see every bin (including intergenic), only genic bins, or only intergenic bins; and whether to zoom into one chromosome or view several/all at full length.
RPKM vs log(observed/expected) — Metagene and Region Browser show raw mean RPKM signal. Enrichment Profile shows a log-ratio against the genome-wide average for that mark/condition, better for comparing marks with very different baseline signal levels.
Posterior probability (Gene Set Comparison tab) — a 0-1 confidence score from ChromstaR that a given mark is genuinely present at a bin, distinct from the RPKM signal used elsewhere.
Performance tips
Computations are vectorised and should stay fast even with all genes selected, but very large gene sets combined with many marks can still take a few seconds.
The Data Table's gene/zone annotation is computed once per session and cached — the first visit to that tab may take a little longer.
Selecting "All chromosomes" in Region Browser, or running thousands of permutations in Gene Set Comparison, will be noticeably slower — start with smaller values to preview, then scale up for a final result.
Contact
Built for chromatin/epigenomics analysis in Echinococcus multilocularis by Janan Gawra.

 linkedin.com/in/janangawra
