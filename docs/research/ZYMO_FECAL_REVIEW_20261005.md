# Zymo fecal TaxoViz run: short research report

2026-10-05 · Descriptive review of saved outputs; no model rerun or taxonomy correction.

The maps recover several strong groups consistent with the added cultures, but also expose mixed bins, classifier disagreement, and strong read-length effects. **Version 2 is a different display, not demonstrated better classification.** All 290,227 plotted reads were matched by ID to the full LRBinner run. Neither the bin labels nor agreement colors are ground truth.

## Data and reproducibility

The subset is the owner's **intentional earliest-sequence 3% selection using a custom ontime-based sampler**, not a random sample. It contains 2.958 Gbp, median length 5,159 bp, and mean GC about 42.9%. The full FASTQ contains 9,674,223 reads and 90.827 Gbp of sequence; its physical size is 184.19 GB / 171.55 GiB. Thus 3% refers to reads, not sequenced bases or file bytes. Full-run and early-run proportions answer different questions.

Both HTML versions contain identical read IDs, lengths, labels and agreement counts. The requested notebook, HTMLs, four SVGs, pickle and log are in [the local run archive](../../runs/zymo_fecal_20261005/README.md). The remote pickle is named `binning_result.pkl`, singular. All 12 downloaded run/source/support files matched remote SHA-256 checksums. The newly saved notebook is retained separately from the stale initial snapshot.

| Setting | Version 1 | Version 2 |
|---|---:|---:|
| t-SNE perplexity | 250 | 100 |
| UMAP neighbors | 100 | 50 |
| UMAP min_dist | 0.30 | 0.15 |
| PCA provenance | Earlier fitted state not independently saved | 43 PCs, 55.0% variance; cap 50 |

The SVG titles establish the three display settings. The current notebook records canonical k4+k6 CLR features, pseudocount 1, GC, Manhattan distance, seed 42 and t-SNE 500 optimization iterations plus early exaggeration. We cannot establish the complete version-1 fitted state from overwritten notebook cells. Old 5,730,402-read and 9,993-read statistics still survive in unrelated cells; they do not describe these plots.

## What the structures mean

TaxoViz makes sequence features → PCA → t-SNE/UMAP, then colors reads by classifier agreement. It does not produce the LRBinner bins or correct labels. Agreement colors identify combinations of tools, not taxa. Genus/species panels share the same coordinates.

The [bin overlay](../../runs/zymo_fecal_20261005/derived/comparison_bins.png) and [taxon overlay](../../runs/zymo_fecal_20261005/derived/comparison_taxa.png) use the original HTML coordinates with every plotted read retained. Bin numbers below refer to LRBinner; taxonomic percentages use agreement between KrakenUniq and Centrifuger (K=C) among **all plotted reads in that bin**, including unresolved reads in the denominator.

| LRBinner bin | Full-run reads | Plotted reads | Main observation |
|---|---:|---:|---|
| 4 | 2,562,702 | 74,268 | Broad central/petal structure; 69.1% K=C *B. fragilis*, with other gut taxa mixed in |
| 5 | 1,959,140 | 53,988 | Large separate UMAP petal structure; 71.5% K=C *E. faecalis* |
| 3 | 797,077 | 23,854 | Well-separated island; 98.5% K=C *P. asaccharolytica* |
| 1 | 556,436 | 18,476 | Low-GC island; 89.0% K=C *F. nucleatum* |
| 2 | 416,364 | 11,466 | High-GC island; 94.2% K=C *K. pneumoniae* |
| 0 | 431,839 | 13,611 | Coherent island, but only 2.5% K=C *P. anaerobius*; strong tool disagreement |
| 7 | 476,742 | 16,557 | Isolated region with conflicting *E. coli*, lambda-phage and other labels |

Bins 18/20/21/22 are mixed by these diagnostic labels. Bin 22 alone includes K=C *B. fragilis* (32.1%), *P. vulgatus* (8.2%) and *Blautia wexlerae* (5.5%). Conversely, one taxon can span bins: *P. asaccharolytica* also appears in bin 10, *K. pneumoniae* in bin 6, and *E. faecalis* in bin 18. This supports neither “one island = one species” nor “one bin = one species.”

The [length/GC views](../../runs/zymo_fecal_20261005/derived/comparison_gradients.png) show long reads along many petals/loops and shorter reads toward dense cores. GC separates major groups: bin 1 averages 27.0%, bin 5 37.6%, bin 3 52.9%, bin 2 58.4%. These are observations. Length-dependent sparse k6 counts and count-pseudocount CLR, local genomic composition, and nonlinear projection are plausible contributors to the petals; errors, strains, plasmids or distinct species are not established by shape. Quality was not available for this assessment.

Version 2 has tighter UMAP groups and a rearranged t-SNE map. For 10,000 fixed query reads and 15 nearest displayed neighbors, same-bin fractions are 89.0%→89.2% for t-SNE and 85.5%→86.4% for UMAP. However, exact neighbor overlap across versions is only 47.1% and 11.7%, respectively. These are **2D concordance diagnostics**, not accuracy, high-dimensional neighborhood preservation, or proof of improvement. UMAP's authors also caution against interpreting its apparent density and separations as original-space structure ([documentation](https://umap-learn.readthedocs.io/en/latest/clustering.html)).

## How LRBinner differs

The inspected run used k4 composition plus 15-mer coverage profiles, an 8-dimensional VAE, then iterative density-valley clustering using cosine-derived distances. The final resumed run used exhaustive search (`-bit 0`) and minimum bin size 500. It found 23 seed clusters, then assigned **1,448,297 remaining reads (15.0%)** by scores from composition/coverage profiles. Every full-run read consequently has a bin; this is not evidence that every read fits its bin confidently. Read-mode clustering here is not HDBSCAN on the TaxoViz plot. Source snapshot: LRBinner `23478e23fd7d74b26a87958cbf1efbcc273009de`; see [method paper](https://doi.org/10.1186/s13015-022-00221-z).

The pickle covers all 9,674,223 positions exactly once and matches `bins.txt`. The exact-ID join has zero missing/duplicate matched IDs, and all full FASTQ lengths match `lengths.txt`. LRBinner's coverage and fitted population are full-run quantities; TaxoViz's representation is fitted to the early subset. Their agreement is useful cross-method evidence, not an independent truth benchmark.

## Zymo reference and qPCR

The [manufacturer's D6323 protocol](https://files.zymoresearch.com/protocols/d6323-zymobiomics_fecal_reference_protocol.pdf) describes pooled human stool and links public DNA/16S characterization reports. The workbook's 354 reference percentages match the public 16S report's 354 values to six decimal places; named leading taxa match too. This strongly identifies the profile's provenance, but does **not** establish the physical sample's catalog number or lot. Reference profiles are method-dependent measurements, not exact read-level truth for the diluted, supplemented sample.

For illustration, the manufacturer's DNA reference has *Bacteroides* 18.44%, *Faecalibacterium* 8.09%, *Blautia* 6.17%, and *B. fragilis* 0.1033% of its reported classified profile. The workbook 16S profile instead has *Blautia wexlerae* 7.542%, *F. prausnitzii* 7.303%, and *B. fragilis* 0.0263%. In this early ONT subset, K=C *B. fragilis* accounts for 21.03% of all reads / 34.14% of bases; *E. faecalis* for 14.64% / 27.15%. The dominance of added-culture candidates is qualitatively consistent with the workbook, not a discrepancy against an unmodified mock standard.

The workbook records 20 µL Zymo, a 10-fold dilution, a 1000-fold culture-sample dilution, and a 6:1 culture:sample mixture. Those notes are the available record, not independently verified preparation details. Its culture list includes *F. nucleatum*, *E. faecalis*, *P. anaerobius*, *P. gingivalis*, *P. asaccharolytica*, PKS-labelled *Klebsiella*, and BFT-labelled *B. fragilis*.

**qPCR cannot supply quantitative truth here.** Sheet 2 omits the ×6 culture contribution for *E. faecalis* and *F. nucleatum* while applying it to several other targets. “PEP stom” conflicts with the *P. anaerobius* culture label; the PKS description names *E. coli* while the culture list names *Klebsiella*. BFT/PKS marker copies are not interchangeable with universal 16S copies or organism abundance. N/A is not zero; very late Cq results need assay controls/limits. Formulas and cached values are preserved in [the qPCR audit](../../runs/zymo_fecal_20261005/derived/qpcr_formulas.csv).

## Most useful follow-ups

1. **Audit bin 0 and its classifier references.** KrakenUniq calls 13,368/13,611 reads *P. anaerobius*, MetaPhlAn calls 3,757, while Centrifuger spreads them across other taxa or leaves them unresolved. The culture record supports targeted validation, not automatic reassignment.
2. **Inspect bins 7/11/13 for mobile/control sequence.** Bin 7 has 4,792 Centrifuger lambda calls and almost no MetaPhlAn species calls. Shared/mobile sequence, control material and database/parser effects are hypotheses requiring alignment evidence; do not label this island simply *E. coli* or delete it.
3. **Fix agreement semantics before using colors for correction.** Of 43,794 species-level `K_only` reads, 37,289 have conflicting classified labels; 5,617/6,087 `M_only` reads also contain conflicts. “Only” currently hides disagreement. Genus-wide three-tool agreement is 20.94%, species-wide 11.07%; neither is calibrated confidence.
4. **Freeze complete per-version configuration and fitted transforms, then test length/GC effects.** Retain the intentional early-time design and compare the same IDs. Assess changes against independent sequence evidence and rare/unknown retention, not visual compactness. P0 correctness tasks remain open.

Detailed calculations, limits and provenance: [experiment record](../computing/experiments/EXP-20261005-ZYMO.md). No biological accuracy, strain identity, toxin carriage, or correction gain was validated.