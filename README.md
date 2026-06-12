# TaxoViz — Setup and Usage Guide

Read-level k-mer ordination coloured by classifier agreement across
MetaPhlAn4, KrakenUniq, and Centrifuger.

---

## File overview

| File | Purpose |
|------|---------|
| `taxoviz_core.py` | Shared library — k-mer vectorisation, CLR, parsers, taxonomy, agreement logic |
| `taxoviz_notebook.ipynb` | Interactive notebook for 10k-read subset exploration |
| `taxoviz_full.py` | Production CLI for full dataset (15M+ reads, streaming, GPU) |
| `environment.yml` | Conda environment specification |

All four files must be in the **same directory**.

---

## 1. Environment setup

```bash
mamba env create -f environment.yml
conda activate taxoviz
```

### GPU UMAP (optional but recommended for >500k reads)

Install RAPIDS cuML matching your CUDA version:
```bash
# CUDA 11.8 example — check https://rapids.ai/start.html for your version
mamba install -c rapidsai -c conda-forge cuml=24.04 cudatoolkit=11.8
```

If cuML is not installed, `taxoviz_full.py` automatically falls back to
CPU `umap-learn` (slower but correct).

### NCBI taxonomy files

Download once and reuse:
```bash
mkdir -p /data/ncbi_taxonomy
cd /data/ncbi_taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz   # extracts nodes.dmp, names.dmp, etc.
```

---

## 2. Classifier output requirements

### MetaPhlAn4
Run with `--samout` to produce a SAM file:
```bash
metaphlan sample.fastq.gz \
    --input_type fastq \
    --samout metaphlan4.sam \
    --read_min_len 100 \
    --min_alignment_len 100 \
    -o metaphlan4_profile.txt
```

### KrakenUniq
Run with per-read output enabled:
```bash
krakenuniq \
    --db /data/krakenuniq_db \
    --output krakenuniq.out \
    --report krakenuniq_report.txt \
    sample.fastq.gz
```

### Centrifuger
```bash
centrifuger \
    -x /data/centrifuger_nt \
    -U sample.fastq.gz \
    -S centrifuger.out \
    --report-file centrifuger_report.txt
```

---

## 3. Interactive notebook (10k subset)

```bash
jupyter notebook taxoviz_notebook.ipynb
```

**Workflow:**
1. **Cell 1 — CONFIG**: set paths to your FASTQ, classifier outputs, and taxonomy directory. Adjust hyperparameters as needed.
2. **Cells 2–6**: run sequentially to load data, parse classifiers, normalise taxonomy, and build the feature matrix.
3. **Cell 7 — PCA**: inspect the scree plot to see how many PCs are selected.
4. **Cell 8 — t-SNE**: run openTSNE. Tweak `tsne_perplexity` in CONFIG and re-run to compare.
5. **Cell 9 — UMAP**: run UMAP. Tweak `umap_n_neighbors` / `umap_min_dist` and re-run.
6. **Cell 10**: save the combined 2×2 static figure (SVG + PNG).
7. **Cell 11**: save the interactive HTML with per-read hover tooltips.
8. **Cell 12**: print and save the agreement statistics table.
9. **Cell 13 — Parameter sweep**: compare multiple perplexity/n_neighbors values side by side without re-running the full pipeline.
10. **Cell 14**: colour the ordination by read length and GC content to identify error-enriched diffuse clouds.

---

## 4. Full-dataset CLI (`taxoviz_full.py`)

### Basic usage

```bash
python taxoviz_full.py \
    --fastq         /data/sample.fastq.gz \
    --metaphlan-sam /data/metaphlan4.sam \
    --krakenuniq    /data/krakenuniq.out \
    --centrifuger   /data/centrifuger.out \
    --taxonomy-dir  /data/ncbi_taxonomy/ \
    --out-dir       /results/taxoviz/
```

### Passing hyperparameters tuned in the notebook

```bash
python taxoviz_full.py \
    --fastq         /data/sample.fastq.gz \
    --metaphlan-sam /data/metaphlan4.sam \
    --krakenuniq    /data/krakenuniq.out \
    --centrifuger   /data/centrifuger.out \
    --taxonomy-dir  /data/ncbi_taxonomy/ \
    --out-dir       /results/taxoviz/ \
    --tsne-perplexity 50 \
    --umap-neighbors  30 \
    --umap-min-dist   0.05 \
    --n-jobs          128
```

### All options

```
Required:
  --fastq             Input FASTQ (.gz ok)
  --metaphlan-sam     MetaPhlAn4 --samout SAM file
  --krakenuniq        KrakenUniq per-read output
  --centrifuger       Centrifuger per-read output
  --taxonomy-dir      Directory with nodes.dmp and names.dmp
  --out-dir           Output directory (created if absent)

k-mer features:
  --k-values          k values to concatenate (default: 4 6)
  --clr-pseudocount   CLR pseudocount δ (default: 1.0)

PCA:
  --pca-components    Max PCA components (default: 75)
  --pca-variance      Variance threshold for PC selection (default: 0.75)

t-SNE:
  --tsne-perplexity   (default: 30)
  --tsne-metric       (default: manhattan)
  --tsne-n-iter       (default: 750)

UMAP:
  --umap-neighbors    (default: 30)
  --umap-min-dist     (default: 0.1)
  --umap-metric       (default: manhattan)

Compute:
  --chunk-size        Reads per chunk for streaming (default: 50000)
  --n-jobs            Parallel jobs, -1 = all cores (default: -1)
  --seed              Random seed (default: 42)
  --no-gpu            Disable GPU UMAP, force CPU umap-learn
```

### Outputs

All files are written to `--out-dir`:

| File | Description |
|------|-------------|
| `taxoviz_static.svg` / `.png` | 2×2 grid: t-SNE + UMAP × genus + species |
| `taxoviz_interactive.html` | Plotly WebGL interactive plot with per-read hover |
| `taxoviz_agreement_stats.csv` | Per-category counts and percentages |
| `taxoviz_embeddings.npz` | Raw coordinates: `tsne`, `umap`, `read_ids` arrays |
| `taxoviz_metadata.parquet` | Full per-read metadata DataFrame |

Load embeddings later:
```python
import numpy as np
data = np.load("taxoviz_embeddings.npz", allow_pickle=True)
read_ids = data["read_ids"]
X_tsne   = data["tsne"]    # shape (n, 2)
X_umap   = data["umap"]    # shape (n, 2)
```

---

## 5. Scalability and expected runtimes

| Dataset size | Mode | t-SNE | UMAP | Total |
|---|---|---|---|---|
| 10k reads | In-memory | ~1 min | ~30 s | ~5 min |
| 500k reads | In-memory | ~30 min | ~10 min | ~1 h |
| 15M reads | Streaming | ~2–4 h | ~30–60 min (GPU) | ~4–6 h |

Streaming mode activates automatically above 500k reads.
Landmark t-SNE (100k landmarks) activates above 100k reads.

**Server specs assumed for 15M estimate:** 256 CPUs, 2 TB RAM, 8× NVIDIA A100 80 GB.

---

## 6. Agreement colour scheme

| Colour | Category | Meaning |
|--------|----------|---------|
| Grey | None | No tool classified this read |
| Orange | MetaPhlAn4 only | Only MetaPhlAn4 assigned a taxon |
| Sky blue | KrakenUniq only | Only KrakenUniq assigned a taxon |
| Green | Centrifuger only | Only Centrifuger assigned a taxon |
| Yellow | MetaPhlAn4 + KrakenUniq | Two tools agree, Centrifuger disagrees/unclassified |
| Vermilion | MetaPhlAn4 + Centrifuger | Two tools agree, KrakenUniq disagrees/unclassified |
| Blue | KrakenUniq + Centrifuger | Two tools agree, MetaPhlAn4 disagrees/unclassified |
| Pink | All three agree | Highest-confidence classification |

Palette: Wong (2011) colorblind-friendly.
