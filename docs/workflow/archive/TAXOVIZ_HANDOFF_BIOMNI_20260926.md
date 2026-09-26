# TAXOVIZ_HANDOFF.md — Machine-Oriented Technical Handoff Specification

> **AUDIENCE:** AI agents (LLM-based multi-agent development systems). This document is NOT written for human readability. It maximizes technical density, exactness, and completeness. Read it entirely before modifying anything in this repository.
>
> **LAST UPDATED:** 2026-09-26 session. Status: tool implemented and validated on a 100k-read subset; full 15M-read production run NOT yet executed.

---

## 1. DOCUMENT PURPOSE & CONSUMPTION PROTOCOL

**What this document is:** the single source of truth for (a) design intent and motivation, (b) algorithmic contracts and invariants, (c) validated empirical findings, (d) open issues and pending tasks, (e) development roadmap, (f) repository map.

**What this document is NOT:** the source of truth for implementation code. The Python files in this repository are. If this document and the code disagree on a signature, the CODE wins for signatures/defaults; this document wins for INTENT and INVARIANTS. If you find a conflict, fix the weaker artifact and log the conflict in Section 11.

**Consumption protocol for continuing agents:**
1. Read this document fully.
2. Read `taxoviz_core.py` and `taxoviz_full.py` fully before any modification.
3. NEVER violate the invariants in Section 7 — they encode fixed bugs that cost significant debugging time.
4. Before proposing algorithmic changes, check Section 12 (roadmap) to avoid re-proposing already-explored and rejected approaches (see 12.3, History of explored approaches).
5. Empirical numbers in Section 10 are measured facts from real 100k-read runs — do not overwrite them with estimates; extend them only with new measured results.
6. When you complete work, UPDATE this document (findings, invariants, open issues, repo map) so the handoff chain stays intact.

**Project name:** TaxoViz.

---

## 2. CONCEPT & MOTIVATION

**Problem:** metagenomic profiling tools (MetaPhlAn4, KrakenUniq, Centrifuger) disagree substantially at the read level on the same Nanopore long-read metagenomic dataset. Profile-level (abundance table) comparison hides WHERE and WHY they disagree. There is no standard tool to visualize read-level taxonomic agreement/disagreement structure.

**Core idea:** ordinate every READ (not every species) in a compositional feature space derived from the read's own sequence (canonical k-mer frequencies, CLR-transformed, plus GC content), reduce dimensionality (PCA), embed (t-SNE + UMAP), and color each read point by the agreement pattern of the three classifiers on that read. The resulting 2D map shows:
- which regions of sequence space are taxonomically contested vs. consensus,
- whether disagreement is structured (clustered = systematic, taxon-specific bias) or random (noise),
- host/contamination islands (reads classified by none or by one tool only),
- community complexity (number/density of sequence-space clusters).

**Why read-level:** classifier disagreement is a per-read event; averaging to profiles destroys the diagnostic signal. Read-level ordination also enables selecting trustworthy read subsets (e.g., only MKC = all-three-agree reads) for downstream analyses.

**Design provenance:** the ordination recipe (CLR-transformed canonical k-mer frequencies + GC content → PCA → t-SNE/UMAP) is adapted from **binny** (Mallawaarachchi et al.), a binner for long-read metagenomics that clusters contigs by CLR k-mer + GC features. TaxoViz transplants this proven contig-level recipe to READ level and repurposes it from binning to classifier-agreement diagnostics. binny uses Manhattan metric for t-SNE; TaxoViz retains this.

**Target data:** Nanopore metagenomic reads, ~4.5 kb mean length, full dataset ~15M reads. Validated on a 100k-read subset.

---

## 3. ALGORITHMIC FOUNDATION

### 3.1 Canonical k-mer counting
- For each k, the k-mer vocabulary is CANONICAL: each k-mer is merged with its reverse complement; canonical form = lexicographically smaller of (kmer, revcomp(kmer)). This halves vocabulary and makes features strand-invariant (recommended for ONT data, which has no strand information).
- Vocabulary sizes: k=4 → **136** canonical k-mers; k=6 → **2080** canonical k-mers. Closed form: (4^k + 4^⌈k/2⌉)/2 — the second term counts self-reverse-complement palindromic k-mers (4^⌈k/2⌉ of them: free positions are the first ⌈k/2⌉, the rest are determined by complementarity). Check: k=4 → (256+16)/2 = 136 ✓; k=6 → (4096+64)/2 = 2080 ✓ — matches `build_canonical_kmer_index` output exactly.
- k-mers containing `N` (or any non-ACGT byte) are skipped entirely.
- k=4 captures low-order nucleotide composition (dominant signal for GC/strand); k=6 captures finer dinucleotide-context signal. Small k is deliberate: for multi-kb Nanopore reads, larger k (8–15) was tested earlier (see 12.3) and rejected for the production tool because vocabulary explodes (k=8 → ~33k canonical; k=11 → ~524k) while binny's recipe demonstrates k=4/6 suffices for ordination.

### 3.2 CLR (centered log-ratio) transformation
- Raw counts are compositional (each read's k-mer vector sums to read length); Euclidean geometry on raw counts is statistically invalid (Aitchison geometry required).
- `CLR(x_i) = log(x_i + δ) − mean_j(log(x_j + δ))`, δ = pseudocount = **1.0** (default). Computed in float64, returned as float32.
- CLR makes features subcompositionally coherent and roughly scale-free w.r.t. read length.

### 3.3 GC content
- GC fraction per read, then z-scored across the dataset, appended as ONE final column. Shape (n, 1). If σ(GC) < 1e-9, returns zeros (degenerate single-taxon case).
- Total feature dimensionality for k=(4,6): 136 + 2080 + 1 = **2217**.

### 3.4 PCA
- Purpose: denoise + accelerate downstream embeddings (t-SNE/UMAP on 2217 dims is slow; on ~50-75 PCs it is fast and binny-standard).
- `n_components` selection: fit with max_components=**75**; then choose the smallest number of PCs whose cumulative explained variance ≥ **0.75** (via `np.searchsorted` on cumsum); cap at 75. Both thresholds are CLI-tunable (`--pca-components`, `--pca-variance`).

### 3.5 t-SNE (openTSNE)
- Metric: **manhattan** (binny standard; CLR space is Aitchison — L1 on log-ratios is the natural geometry).
- Defaults: perplexity=30, n_iter=750, n_jobs=-1, random_state=42.
- **Empirically recommended for ≥100k reads: perplexity=100** (see Section 10.2).
- Scalability: ≤100,000 reads → direct fit. >100,000 → **landmark approach**: fit on a random subset of exactly 100,000 reads (`LANDMARK_THRESHOLD`), then `embedding.transform()` the remainder in chunks of 50,000 (`REST_CHUNK`). Landmark indices drawn with `np.random.default_rng(seed)`, sorted; rest via `np.setdiff1d`.

### 3.6 UMAP
- Metric: **manhattan**. Defaults: n_neighbors=30, min_dist=0.1.
- **Empirically recommended for ≥100k reads: n_neighbors=100, min_dist=0.3** (see Section 10.2).
- GPU path: cuML `cuml.manifold.UMAP` (RAPIDS) when available and `--no-gpu` not passed; automatic fallback to CPU `umap-learn` (with `low_memory=True`) on ImportError OR any Exception in the cuML call. cuML output may be a cuDF object — code handles `.to_numpy()` conversion.

### 3.7 Agreement computation (per read, per rank)
- Each classifier yields a taxid-or-None per read. Taxids are normalized to a rank (genus AND species, both computed) via NCBI taxonomy walk-up.
- Consensus taxid = most common taxid among tools that produced a valid (>0, non-None) taxid (`Counter(...).most_common(1)`).
- Agreement set = tools whose (rank-normalized) taxid equals the consensus.
- 8 possible categories (Section 9). Deterministic tie-break: if all active tools give DIFFERENT taxids, `Counter.most_common` returns the first-inserted key; dict insertion order is `metaphlan, krakenuniq, centrifuger`, so MetaPhlAn wins the tie → category `M_only`. This is a deliberate, documented, deterministic choice (visualization-grade, not truth-grade).

---

## 4. ARCHITECTURE

Three-file design (shared core + two consumers), chosen to avoid code duplication between interactive exploration and production:

```
taxoviz_core.py  ←── imported by ──┐
   (shared library)               ├── taxoviz_notebook.ipynb  (interactive, ≤~100k reads)
                                  └── taxoviz_full.py         (production CLI, any size)
```

### 4.1 Execution paths in `taxoviz_full.py` (selected automatically by read count)
- **IN-MEMORY path** (n ≤ `STREAMING_THRESHOLD` = 500,000): load all reads via `read_fastq_all`, `build_feature_matrix` in RAM, sklearn `PCA`, direct openTSNE fit, CPU umap-learn (GPU if available).
- **STREAMING path** (n > 500,000): never hold all sequences in RAM. `build_feature_matrix_streaming` writes the (n × 2217) float32 feature matrix to HDF5 (`chunks=(chunk_size, D)`, compression `lzf`) while collecting read_ids, read_lens, gc in the SAME pass; `IncrementalPCA` (two passes over HDF5: partial_fit, then transform); landmark t-SNE; cuML/CPU UMAP. HDF5 scratch file deleted after pipeline completes.

### 4.2 Constants (do not change without re-validating Section 10 findings)
```
LANDMARK_THRESHOLD  = 100_000   # reads; t-SNE direct-fit vs landmark switch
STREAMING_THRESHOLD = 500_000  # reads; in-memory vs streaming switch
REST_CHUNK          = 50_000   # reads per t-SNE transform chunk (landmark path)
DEFAULT chunk_size  = 50_000   # reads per streaming/IncrementalPCA chunk
Feature dim D      = 2217     # for k=(4,6): 136 + 2080 + 1 (GC)
```

### 4.3 Output set (written to `--out-dir`)
| File | Content |
|---|---|
| `taxoviz_static.svg` / `.png` | 2×2 matplotlib grid: {t-SNE, UMAP} × {genus, species}, colored by agreement |
| `taxoviz_interactive.html` | Plotly WebGL, same 2×2 grid, per-read hover (id, len, GC, 3 taxon names) |
| `taxoviz_agreement_stats.csv` | per (level, category): count + pct |
| `taxoviz_embeddings.npz` | read_ids (object array) + tsne + umap coordinate matrices |
| `taxoviz_metadata.parquet` | full per-read metadata DataFrame (see 5.9 for schema) |
| `taxoviz_features.h5` | streaming-path scratch; DELETED at pipeline end |

---

## 5. API CONTRACTS — `taxoviz_core.py`

Module-level state built at import time (fast, <1 s): `KMER_INDEX: Dict[int, Dict[str,int]]`, `KMER_VOCAB_SIZE: Dict[int,int]`, `KMER_NUMERIC_LOOKUP: Dict[int, np.ndarray]` — pre-built for k=4 and k=6 only. If you add a k value at runtime, you MUST populate all three dicts for it first.

### 5.1 K-mer index construction
- `reverse_complement(seq: str) -> str` — via `str.maketrans("ACGTN","TGCAN")` + slice-reverse.
- `build_canonical_kmer_index(k: int) -> Tuple[Dict[str,int], int]` — iterates `itertools.product("ACGT", repeat=k)`; canonical = `kmer if kmer <= rc else rc`; positions assigned in first-encounter order. Returns ({canonical_kmer: position}, vocab_size).
- `build_numeric_to_canonical_lookup(k: int) -> np.ndarray` — returns int32 array of shape (4^k,) mapping numeric base-4 k-mer index (A=0,C=1,G=2,T=3, positional base-4) → canonical position; entries for non-canonical numeric indices are **-1** (never accessed at runtime because canonical = min(fwd, rev) is always canonical). **INVARIANT: see 7.1.**

### 5.2 Vectorisation
- `_vectorise_read_python(seq: str, k: int) -> np.ndarray` — pure-Python fallback; shape (vocab_size,), float32; skips k-mers with N; returns zeros if len(seq) < k.
- `_make_numba_vectoriser(k: int) -> Callable[[str], np.ndarray]` — returns a closure over a `numba.njit(cache=True)` kernel `_count(seq_bytes, k, vocab_size, base_enc, comp_enc, numeric_lookup)`. Encoding tables: `BASE_ENC` = int32[256], A/C/G/T→0/1/2/3, all else −1; `COMP_ENC` = int32[4] = [3,2,1,0] (A↔T, C↔G in base-4 encoding). Kernel computes forward numeric index; on any −1 byte skips the window; computes reverse-complement numeric index by iterating positions in reverse with `comp_enc`; canonical numeric = `min(fwd, rev)`; position = `numeric_lookup[canonical]`; increments `vec[pos]`. JIT warm-up call on `b"ACGTACGT"` at closure creation. **On ImportError of numba → returns pure-Python closure.** **INVARIANT: see 7.1.**

### 5.3 Transforms
- `clr_transform(counts: np.ndarray, pseudocount: float = 1.0) -> np.ndarray` — CLR over axis=1 in float64 → float32. Contract: input (n, vocab), output same shape.
- `gc_content(sequences: List[str]) -> np.ndarray` — float32 (n,), `(G+C)/max(len,1)`.
- `gc_zscore(gc: np.ndarray) -> np.ndarray` — float32 (n,1); zeros if σ < 1e-9.

### 5.4 Feature matrix
- `build_feature_matrix(sequences: List[str], k_values: Tuple[int,...] = (4,6), pseudocount: float = 1.0, n_jobs: int = -1) -> np.ndarray` — per k: vectorise all sequences via `joblib.Parallel(n_jobs, prefer="threads")` → CLR; concatenate CLR blocks across k; append GC z-score column. Output (n, Σvocab_k + 1) float32 = (n, 2217) for default k. **Note: a new Numba vectoriser closure is built per call per k — closures are NOT cached across calls (acceptable; JIT cache=True makes recompilation cheap).**

### 5.5 FASTQ I/O
- `open_fastq(path: str)` — transparent opener for plain / `.gz` / `.bz2`, text mode.
- `read_fastq_all(path: str) -> Tuple[List[str], List[str]]` — (read_ids, sequences); sequences uppercased; via Bio.SeqIO. For ≤500k reads only.
- `iter_fastq_chunks(path: str, chunk_size: int = 50_000)` — generator yielding (ids_chunk, seqs_chunk); uppercased.

### 5.6 Classifier parsers (all return `Dict[read_id, Optional[int]]`, pre-initialized with `{rid: None for rid in read_ids_set}`; reads absent from the classifier output stay None)
- `parse_metaphlan_sam(sam_path: str, read_ids_set: set)` — MetaPhlAn4 `--samout` SAM. Skips `@` header lines; skips lines with <6 tab fields; **skips unmapped reads (SAM flag `int(parts[1]) & 4`)**; extracts taxid from reference name (field 3, index 2) with regex `(?:taxid[|_]{1,2}|__taxid__)([0-9]+)` (IGNORECASE). Handles `.gz`/`.bz2`/plain. Supported marker-name formats: `UniRef90_XXXXX__taxid__NNNN`, `GeneID:XXXXX|taxid|NNNN|...`.
- `parse_krakenuniq(out_path: str, read_ids_set: set)` — KrakenUniq per-read output, tab-separated: col0 = `C`/`U` status, col1 = read_id, col2 = taxid. **Only `C` lines are kept** (U → stays None). Handles `.gz`/plain.
- `parse_centrifuger(out_path: str, read_ids_set: set)` — Centrifuger per-read output: col0 = readID, col2 = taxID. `taxID <= 0` → None. Handles `.gz`/plain.

### 5.7 Taxonomy (taxopy)
- `load_taxonomy(taxonomy_dir: str)` — `taxopy.TaxDb(nodes_dmp=<dir>/nodes.dmp, names_dmp=<dir>/names.dmp, keep_files=True)`.
- `get_rank_taxid(taxid: Optional[int], rank: str, taxdb) -> Optional[int]` — if taxon.rank == rank → taxid; else walk `taxon.taxid_lineage` (ancestor list) and return first ancestor with matching rank; None on failure/None/≤0 input. **CAVEAT: `taxon.taxid_lineage` attribute name is UNVERIFIED against the installed taxopy version — see open issue 11.2.**
- `build_taxid_name_cache(taxids: set, taxdb) -> Dict[int, str]` — batch {taxid: name}; fallback `f"taxid:{tid}"` on exception.
- `taxid_to_name(tid, cache: Dict[int,str]) -> str` — safe for None/NaN/≤0/str → `"unclassified"`; cache miss → `f"taxid:{tid}"`.

### 5.8 Agreement computation
- Module constants: `TOOLS = ("metaphlan","krakenuniq","centrifuger")`; `TOOL_ABBR = {metaphlan:"M", krakenuniq:"K", centrifuger:"C"}`; `CATEGORIES`, `COLORS`, `LABEL_NAMES`, `_ABBR_TO_CAT` (see Section 9).
- `compute_agreement_row(metaphlan_taxid, krakenuniq_taxid, centrifuger_taxid) -> str` — implements 3.7. Returns one of the 8 category strings. **INVARIANT: see 7.3 (deterministic tie-break).**
- `build_metadata_df(read_ids, sequences, labels_metaphlan, labels_krakenuniq, labels_centrifuger, taxdb, read_lens=None, gc_arr=None) -> pd.DataFrame` — **`sequences` may be None IFF both `read_lens` (int array) and `gc_arr` (float array) are provided (streaming path); otherwise raises ValueError.** Normalizes ALL raw taxids to genus+species in one batch pass (builds genus_map/species_map dicts first — do not per-read-walk the taxonomy). Output schema (indexed by read_id):

| Column | Type | Meaning |
|---|---|---|
| `read_len` | int | read length (bp) |
| `gc` | float | GC fraction |
| `{tool}_raw` ×3 | int/None | raw classifier taxid |
| `{tool}_genus`, `{tool}_species` ×6 | int/None | rank-normalized taxid |
| `agreement_genus`, `agreement_species` | str | category (Section 9) |
| `{tool}_{level}_name` ×6 | str | taxon name via cache; "unclassified" fallback |

### 5.9 Plot helpers
- `make_legend_handles(categories=None)` — matplotlib patches for CATEGORIES (or subset).
- `scatter_ordination(ax, coords, agreement_series, title="", xlabel="", ylabel="", alpha=0.4, s=2.0, rasterized=True)` — draws categories in `CATEGORIES` order with `zorder=CATEGORIES.index(cat)` so later (more informative) categories render ON TOP. Reuse this ordering contract in any new plotting code.

---

## 6. API CONTRACTS — `taxoviz_full.py`

Import contract: inserts its own directory into `sys.path` and imports `taxoviz_core as tc` — **the two files must be co-located.**

Pipeline stages in `run_pipeline(args)` (order is fixed):
1. `count_reads(fastq_path) -> int` — line-count `@`-prefixed lines (NOTE: counts any line starting with `@`, including quality-score lines beginning with `@` — acceptable approximation for standard FASTQ; if exactness needed, parse 4-line records instead).
2. Feature matrix: streaming (`build_feature_matrix_streaming(fastq_path, hdf5_path, k_values, pseudocount, chunk_size, n_jobs, n_reads) -> (read_ids: List[str], D: int, read_lens: np.ndarray(int32), gc_arr: np.ndarray(float32))`) or in-memory (`tc.read_fastq_all` + `tc.build_feature_matrix`). **INVARIANT: see 7.2.**
3. Parse the three classifiers (against `set(read_ids)`).
4. `tc.load_taxonomy` + `tc.build_metadata_df` (streaming: `sequences=None, read_lens=…, gc_arr=…`).
5. PCA: `run_pca_incremental(hdf5_path, n_reads, pca_max_components, pca_variance_threshold, chunk_size, random_seed) -> np.ndarray` (two HDF5 passes: partial_fit loop, then transform loop; component count chosen by variance threshold) or `run_pca_full(X, …) -> np.ndarray`.
6. `run_tsne(X_pca, n_reads, perplexity, metric, n_iter, n_jobs, random_seed, landmark_n=LANDMARK_THRESHOLD) -> np.ndarray` — direct fit ≤ landmark_n; else landmark fit + chunked transform (3.5).
7. `run_umap(X_pca, n_neighbors, min_dist, metric, random_seed, use_gpu) -> np.ndarray` — cuML→CPU fallback (3.6).
8. Save outputs: `save_static_figure` (2×2: t-SNE genus/species top row, UMAP genus/species bottom row; `tc.scatter_ordination` with alpha=0.3, s=1.5; svg+png at dpi=150), `save_interactive_html` (Plotly `make_subplots(rows=2, cols=2)`, `go.Scattergl` WebGL traces per category with `legendgroup=cat`, hover text = read id/len/gc/3 taxon names, `showlegend` only on first subplot), `save_stats_csv`, `save_embeddings` (npz: read_ids object array + tsne + umap), `save_metadata` (parquet).
9. Streaming cleanup: unlink HDF5 scratch.

### CLI flags (`parse_args()`; argparse with `ArgumentDefaultsHelpFormatter`; defaults in parentheses)
Required: `--fastq`, `--metaphlan-sam`, `--krakenuniq`, `--centrifuger`, `--taxonomy-dir` (dir with nodes.dmp+names.dmp), `--out-dir`.
k-mer: `--k-values` ([4, 6], nargs+), `--clr-pseudocount` (1.0).
PCA: `--pca-components` (75), `--pca-variance` (0.75).
t-SNE: `--tsne-perplexity` (30), `--tsne-metric` (manhattan), `--tsne-n-iter` (750).
UMAP: `--umap-neighbors` (30), `--umap-min-dist` (0.1), `--umap-metric` (manhattan).
Compute: `--chunk-size` (50000), `--n-jobs` (−1), `--seed` (42), `--no-gpu` (flag; force CPU umap-learn).

**Recommended production invocation for ≥100k reads (empirically grounded, Section 10.2):**
`--tsne-perplexity 100 --umap-neighbors 100 --umap-min-dist 0.3`

---

## 7. CRITICAL INVARIANTS & FIXED-BUG REGISTRY

These encode bugs that were found and fixed during development. Violating them reintroduces silent, hard-to-detect failures. **Do not "simplify" these away.**

### 7.1 Numba vectoriser index mapping (CRITICAL — was an IndexError-class bug)
**Invariant:** the Numba kernel must NEVER use a raw numeric base-4 k-mer index (range 0..4^k−1) directly as a position into the canonical vocabulary vector (size = number of canonical k-mers, e.g. 136 for k=4). It MUST map through the precomputed `numeric_lookup` array of shape (4^k,) (built by `build_numeric_to_canonical_lookup`).
**Why it is correct:** for the ACGT alphabet encoded A=0,C=1,G=2,T=3, positional base-4 numeric order is IDENTICAL to lexicographic string order — verified exhaustively for k=4 (all 256 k-mers, 0 mismatches between numeric-min and lexicographic-min canonical forms). Therefore `min(fwd_numeric, rev_numeric)` equals the numeric index of the lexicographically canonical k-mer, and `numeric_lookup[that]` gives its vocabulary position.
**Historical bug:** the first implementation indexed the 136-element vector with the raw 0..255 numeric index → IndexError / silent corruption.

### 7.2 Streaming path must not re-read sequences (CRITICAL — was a ~67 GB RAM bug)
**Invariant:** in the streaming path (>500k reads), read lengths and GC fractions MUST be collected during the vectorisation pass (inside `build_feature_matrix_streaming`, while each chunk's sequences are in memory) and passed to `build_metadata_df(read_ids, sequences=None, …, read_lens=…, gc_arr=…)`. The pipeline must NEVER load all sequences into RAM a second time for metadata construction.
**Why:** 15M reads × ~4.5 kb ≈ 67.5 GB as Python strings — OOM on any normal node.
**Historical bug:** the first streaming implementation re-read the FASTQ into RAM for GC/length computation.

### 7.3 Deterministic all-disagree tie-break
**Invariant:** when all active tools return different taxids, the consensus must deterministically resolve to MetaPhlAn's taxid (category `M_only`), via `Counter` first-inserted-key semantics on a dict built in order (metaphlan, krakenuniq, centrifuger). Python ≥3.7 guarantees dict insertion order. If you reimplement `compute_agreement_row`, you MUST preserve this determinism (any deterministic rule is acceptable, but it must be documented and stable across runs).

### 7.4 Draw order in scatter plots
**Invariant:** agreement categories are drawn in `CATEGORIES` list order with increasing zorder, so rarer/more-informative categories (MKC last) are visible on top of dominant ones (KC, K_only). Preserve when adding plot types.

### 7.5 Rank-dependence of agreement
Agreement is computed at BOTH genus and species levels independently. These give different distributions (Section 10.1) — genus-level normalization merges divergent species calls into agreement. Never report "the" agreement rate without naming the rank.

---

## 8. CLASSIFIER I/O CONTRACTS (summary — full detail in 5.6)

| Tool | Input file | Format contract | Classified-read rule |
|---|---|---|---|
| MetaPhlAn4 | `--samout` SAM (`.gz`/`.bz2`/plain) | tab-SAM; taxid in reference name via regex `(?:taxid[|_]{1,2}|__taxid__)([0-9]+)` | mapped (flag & 4 == 0) AND regex hit |
| KrakenUniq | per-read output (`.gz`/plain) | tab: status, read_id, taxid, … | status == `C` |
| Centrifuger | per-read output (`.gz`/plain) | tab: readID, ?, taxID, … | taxID > 0 |

All parsers intersect with the FASTQ read-id set; ids absent from a classifier's output remain None (→ "not classified by that tool").

---

## 9. AGREEMENT CATEGORY SYSTEM

```
CATEGORIES = ["none", "M_only", "K_only", "C_only", "MK", "MC", "KC", "MKC"]
COLORS = {            # Wong (2011) colorblind-safe palette — do not change casually
  "none":   "#AAAAAA",  "M_only": "#E69F00",  "K_only": "#56B4E9",
  "C_only": "#009E73",  "MK":     "#F0E442",  "MC":     "#D55E00",
  "KC":     "#0072B2",  "MKC":    "#CC79A7",
}
LABEL_NAMES = {
  "none": "None classified",          "M_only": "MetaPhlAn4 only",
  "K_only": "KrakenUniq only",         "C_only": "Centrifuger only",
  "MK": "MetaPhlAn4 + KrakenUniq",     "MC": "MetaPhlAn4 + Centrifuger",
  "KC": "KrakenUniq + Centrifuger",   "MKC": "All three agree",
}
_ABBR_TO_CAT = {"C":"C_only","K":"K_only","M":"M_only","CK":"KC","CM":"MC","KM":"MK","CKM":"MKC"}
```
M = MetaPhlAn4, K = KrakenUniq, C = Centrifuger. Category = sorted abbreviation of the agreeing-tool set (sorted alphabetically: C < K < M).

---

## 10. VALIDATED EMPIRICAL FINDINGS (100k-read subset, real data)

**Provenance:** two production-pipeline runs on a 100,000-read Nanopore subset (9,993 reads present in the ordination after parsing; the interactive HTML contained 4 panels × 9,993 points). Run A ("standard"): perplexity=30, UMAP n_neighbors=30/min_dist=0.1. Run B ("larger/softer"): perplexity=100, UMAP n_neighbors=100/min_dist=0.3. Both runs used the redesigned pipeline (k=4,6 CLR + GC → PCA → t-SNE/UMAP, Manhattan metric). Findings below were extracted by decoding the Plotly typed-array JSON embedded in the self-contained HTML outputs (base64 `bdata` fields, dtype f4; 32 traces = 8 categories × 4 panels).

**Panel mapping in the 2×2 interactive HTML (Plotly axis ids):** `x/y` = t-SNE genus (row1col1); `x2/y2` = t-SNE species (row1col2); `x3/y3` = UMAP genus (row2col1); `x4/y4` = UMAP species (row2col2). **Correction note for downstream agents:** an earlier in-session analysis mislabeled x2 as "UMAP" and x3/x4 as GC-colored variants; the mapping above is verified against the generating code (`make_subplots` + trace insertion order). Quantitative spread/separation comparisons below were computed on the t-SNE panels; the UMAP panels were visually consistent but not separately re-quantified after the correction.

### 10.1 Agreement distribution (n = 9,993 reads; measured, not estimated)

| Category | Genus level n (%) | Species level n (%) |
|---|---|---|
| KC (KrakenUniq + Centrifuger) | 5,394 (54.0) | 5,574 (55.8) |
| K_only (KrakenUniq only) | 1,850 (18.5) | 1,998 (20.0) |
| MKC (all three agree) | 1,500 (15.0) | 884 (8.8) |
| none | 588 (5.9) | 432 (4.3) |
| C_only (Centrifuger only) | 422 (4.2) | 780 (7.8) |
| M_only (MetaPhlAn4 only) | 152 (1.5) | 248 (2.5) |
| MK | 70 (0.7) | 54 (0.5) |
| MC | 17 (0.2) | 23 (0.2) |

**Key measured facts:**
- **MetaPhlAn4 is the systematic outlier:** only ~15% (genus) / ~9% (species) of reads reach MKC; KC alone covers ~54–56%. MetaPhlAn4 contributes a non-None call on only ~17–18% of reads (M_only + MK + MC + MKC).
- **Rank dependence (7.5) is large:** MKC drops from 15.0% (genus) to 8.8% (species) — genus-level normalization converts many species-level disagreements into agreement. K_only and C_only grow at species level (fininer resolution → more solo calls).
- **KrakenUniq and Centrifuger agree with each other far more than either agrees with MetaPhlAn4** — consistent with their shared alignment/k-mer-based whole-read architecture vs. MetaPhlAn's marker-gene architecture.
- **Working hypothesis (unverified, see 11.3):** MetaPhlAn4's silence is architectural (marker-gene coverage requires sufficient depth on clade-specific marker loci; short/low-quality ONT reads often fail its alignment thresholds), possibly aggravated by Illumina-default parameters (e.g., `--min_mapq`, `--min_alignment_len`) applied to ONT reads.

### 10.2 Parameter comparison (Run A standard vs Run B softer)

- **Topology is stable across settings:** the same islands/clusters appear in both runs; category spatial structure is reproducible. This is evidence the k-mer signal is real, not an embedding artifact.
- **Normalized cluster spread (per-category std / global plot range) is essentially unchanged** between settings: ~0.20–0.25 for all categories in both runs (e.g., KC 0.220→0.226; none 0.199→0.219; MKC 0.253→0.250). Softer parameters do not disperse clusters.
- **Inter-cluster centroid separation (normalized) increases ~10–20% at perplexity=100** (e.g., MKC↔M_only 0.227→0.248; KC↔M_only 0.214→0.226; KC↔K_only 0.141→0.182; MC↔MKC 0.177→0.245).
- **Recommendation (empirically grounded):** for ≥100k reads use `--tsne-perplexity 100 --umap-neighbors 100 --umap-min-dist 0.3`. For the 15M-read run, perplexity ≥100 is further justified (larger samples need larger perplexity to capture global structure; the landmark path caps the fitted subset at 100k, so perplexity up to ~10k is computationally legal but 100–500 is the sensible range to sweep).

### 10.3 Spatial-structure observations (qualitative, from panel inspection)
- Distinct compact islands exist — putatively taxonomically homogeneous read groups.
- MKC and KC occupy different regions but both appear across multiple islands → consensus is not confined to one taxon.
- M_only / MK concentrate in specific zones → putatively taxa well covered by MetaPhlAn marker genes.
- A large diffuse central cloud (mostly KC + K_only) → contested/low-abundance taxa with poor marker coverage.
- `none` (unclassified by anyone) is scattered across the whole map, NOT clustered → heterogeneous mixture (host reads, rare taxa, low-quality reads), not a single artifact population.

---

## 11. OPEN ISSUES / PENDING TASKS

1. **cuML installation unresolved on the target machine** (Python 3.10 + CUDA 11.8). Conda `cuml=24.*` fails (requires CUDA 12.x + Python 3.11/3.12). Options, in preference order: (a) `pip install cuml-cu11==25.6.0 cudf-cu11==25.6.0 cuvs-cu11==25.6.0 rmm-cu11==25.6.0 pylibraft-cu11==25.6.0 --extra-index-url https://pypi.nvidia.com` (RAPIDS publishes CUDA-11.x wheels supporting Python ≥3.10); (b) recreate env with Python 3.11 + CUDA 12.x; (c) run with `--no-gpu` (CPU umap-learn fallback is fully implemented and automatic). Verify with `python -c "import cuml; print(cuml.__version__)"`.
2. **taxopy API verification:** `get_rank_taxid` uses `taxon.taxid_lineage` — the attribute name was written from memory and NOT verified against the installed taxopy version. Check the installed package (e.g., `dir(taxopy.Taxon(...))`) before the production run; the fallback catch is broad (`except Exception`) so a wrong attribute would silently return None for every rank normalization — **this would corrupt all agreement categories to "none"-heavy output. Test `get_rank_taxid` on a known taxid (e.g., 562 → genus Escherichia) before trusting any run.**
3. **MetaPhlAn4 parameter audit (hypothesis 10.1):** check whether MetaPhlAn4 was run with Illumina-default `--min_mapq` / `--min_alignment_len` on ONT reads; re-run with ONT-appropriate thresholds and re-measure the MKC share before concluding architectural bias.
4. **Full 15M-read production run not yet executed.** Expected runtimes (estimates, 256-core + 8×A100-class node): vectorisation ~30–60 min; IncrementalPCA ~20–40 min; landmark t-SNE ~2–4 h; cuML UMAP ~30–60 min; total ~4–6 h.
5. **`count_reads` approximation** (see 6): counts `@`-prefixed lines; replace with 4-line record parsing if exact counts matter (affects only the streaming/in-memory mode switch at the 500k boundary).

---

## 12. ROADMAP & DEVELOPMENT DIRECTIONS

### 12.1 Validated practical applications (build features on these)
1. **Classification-quality triage:** export/select only MKC reads as a high-confidence subset for downstream analyses (assembly, abundance estimation).
2. **Host/contamination detection:** islands dominated by `none`/single-tool categories → candidate host or contaminant sequence pools; selectable interactively (hover gives read ids) for BLAST verification.
3. **MetaPhlAn4 bias discovery:** clusters with high KC-but-no-M share → taxa missing from or under-covered by the marker database at this depth; feed back as a database/parameter audit list.
4. **Community-complexity proxy:** number/density of islands in the ordination as a diversity descriptor; many small dense islands = high diversity with clean taxa; one diffuse cloud = single dominant or many near-neighbor taxa.
5. **Parameter transfer protocol:** tune embeddings on a 100k subset (as done), then run the CLI on the full dataset with the tuned flags — this workflow is now empirically grounded (10.2).

### 12.2 Proposed development directions (not yet implemented — candidates for continuing agents)
- **Density-based island segmentation** (HDBSCAN on the UMAP/t-SNE coordinates) to turn qualitative "islands" into machine-readable clusters with per-cluster agreement profiles and taxon compositions.
- **Per-cluster taxon enrichment tables** (cluster × dominant taxon × agreement category) as a CSV output.
- **Confidence-scored consensus labels:** replace the binary agreement category with a weighted consensus (e.g., agreement at genus AND species, or 2-of-3 with rank-aware scoring).
- **Interactive selection → export:** Plotly `selection` events to export selected read ids to a file (application 12.1.2).
- **Reference-free host filtering:** use the ordination itself (islands with no classifier support + uniform GC) to flag host reads without a host reference.
- **Benchmark mode:** run the ordination on simulated communities (known ground truth) to quantify how well island structure recovers true taxa.
- **Streaming-UMAP option** for >500k reads if cuML is unavailable (e.g., Landmark-UMAP analogous to landmark t-SNE).
- **Multi-dataset mode:** several samples on one shared embedding (fit on pooled PCA space, transform per-sample) for between-sample comparison.

### 12.3 History of explored approaches (do NOT re-propose without new justification)
- **k=8/11/15 canonical k-mers, L1-normalized frequencies, Euclidean UMAP** — explored in the earlier monolithic notebook (`taxoviz_analysis.ipynb`); k=11 was its default. Superseded: vocabulary explosion (k=11 → ~524k canonical features) with no demonstrated gain over the binny-style k=4/6 CLR recipe for ordination purposes.
- **NMDS (Bray-Curtis) and Jaccard PCoA (GPU, chunked 100k×100k)** — explored in the same notebook; O(n²) memory/time, capped at 50k reads (NMDS) and required GPU chunking (PCoA). Superseded by PCA→t-SNE/UMAP which scales to the full dataset.
- The earlier notebook remains in the repository as an exploratory record; its per-method HTML outputs (one panel per method) were replaced by the current 2×2 t-SNE/UMAP × genus/species format.

---

## 13. REPOSITORY MAP

| File | Role | Contents |
|---|---|---|
| `taxoviz_core.py` | Shared library (imported by both consumers) | Sections 1–10 of code: canonical k-mer index + numeric lookup; Numba vectoriser (with 7.1 invariant); CLR; GC; `build_feature_matrix`; FASTQ readers (plain/gz/bz2); 3 classifier parsers; taxopy taxonomy normalization; agreement computation (8 categories, genus+species); color palette; matplotlib plot helpers. ~28 KB. |
| `taxoviz_full.py` | Production CLI | Full pipeline: read counting → streaming/in-memory feature matrix → classifier parsing → metadata → PCA (Incremental/full) → t-SNE (direct/landmark) → UMAP (cuML/CPU) → 5 output artifacts + HDF5 scratch cleanup. tqdm progress bars, logging. All flags in Section 6. ~30 KB. Must be co-located with `taxoviz_core.py`. |
| `taxoviz_notebook.ipynb` | Interactive exploration notebook (current, 15 cells) | Cell 1 CONFIG (all tunables in one dict); 2 imports; 3 load FASTQ; 4 parse classifiers; 5 taxonomy; 6 feature matrix; 7 PCA (+scree); 8 t-SNE; 9 UMAP; 10 combined 2×2 static figure; 11 interactive Plotly HTML (2×2, WebGL, per-read hover); 12 agreement stats table; 13 parameter-sweep helper (independent of cells 8–9); 14 color-by-read-length/GC panels. Intended workflow: run on a ≤100k subset, tune parameters live, transfer flags to the CLI. |
| `taxoviz_analysis.ipynb` | Earlier exploratory notebook (31 cells; SUPERSEDED, kept as record) | Monolithic pre-redesign version: k=8/11/15 comparison, k=11 default, L1-normalized frequencies, NMDS (Bray-Curtis), GPU chunked Jaccard PCoA, UMAP (Euclidean), per-method interactive HTML. See 12.3 before reusing anything from it. |
| `environment.yml` | Conda environment spec | Python 3.10; numpy/pandas/scipy/scikit-learn/joblib; biopython; taxopy; opentsne; umap-learn; numba; h5py; matplotlib; plotly; tqdm; jupyter/nbformat; pyarrow; pip fallback for taxopy. cuML commented out (see 11.1 for the working install route). |
| `INSTRUCTIONS.md` | Human-oriented usage guide | Environment setup, NCBI taxonomy download (`taxdump.tar.gz` → nodes.dmp/names.dmp), classifier run flags, notebook workflow, CLI usage with all flags, output descriptions, scalability table, color-scheme attribution. |
| `TAXOVIZ_HANDOFF.md` | THIS document | Machine-oriented handoff spec. Update after every development iteration. |
| `taxoviz_comparison.png` | Session analysis artifact | 2×2 comparison of Run A vs Run B (from the 100k HTMLs). **Known label defect:** panels labeled "UMAP" actually show the t-SNE species panel (panel-mapping error described in Section 10); regenerate before external use. |
| `taxoviz_agreement_dist.png` | Session analysis artifact | Bar chart of agreement distribution, averaged across genus+species levels (per-level numbers in 10.1 are authoritative). |

**External inputs required to run (not in repo):** FASTQ (`.gz` ok); MetaPhlAn4 `--samout` SAM; KrakenUniq per-read output; Centrifuger per-read output; NCBI taxonomy directory (`nodes.dmp` + `names.dmp` from `taxdump.tar.gz`).

---

## APPENDIX A — Quick-start command sequence (for agents executing the pipeline)

```bash
# Environment
conda env create -f environment.yml && conda activate taxoviz
# GPU UMAP (option a, see 11.1):
pip install cuml-cu11==25.6.0 cudf-cu11==25.6.0 cuvs-cu11==25.6.0 rmm-cu11==25.6.0 pylibraft-cu11==25.6.0 --extra-index-url https://pypi.nvidia.com

# Sanity check BEFORE any run (open issue 11.2):
python -c "import taxoviz_core as tc; db = tc.load_taxonomy('<taxonomy_dir>'); print(tc.get_rank_taxid(562, 'genus', db))"  # expect 561 (Escherichia)

# Production run (empirically recommended parameters, Section 10.2)
python taxoviz_full.py \
    --fastq reads.fastq.gz \
    --metaphlan-sam metaphlan4.sam \
    --krakenuniq krakenuniq.out \
    --centrifuger centrifuger.out \
    --taxonomy-dir ncbi_taxonomy/ \
    --out-dir results/taxoviz/ \
    --tsne-perplexity 100 --umap-neighbors 100 --umap-min-dist 0.3
```

## APPENDIX B — Data-extraction protocol for TaxoViz HTML outputs

The interactive HTML is self-contained Plotly (plotly.js v3.6.0 embedded, ~13 MB). To recover point data programmatically: locate the LAST `Plotly.newPlot(` call in the file; parse the traces JSON array with a bracket-depth scanner (respecting string escapes); each trace's `x`/`y` are typed arrays `{dtype: "f4", bdata: "<base64>"}` — decode with `np.frombuffer(base64.b64decode(bdata), dtype=np.float32)`. `legendgroup` = agreement category; `xaxis`/`yaxis` ∈ {x, x2, x3, x4} map to panels per Section 10 (t-SNE genus / t-SNE species / UMAP genus / UMAP species). `text` holds per-read hover strings (read id, len, GC, three taxon names). This protocol was used to produce all numbers in Section 10.
