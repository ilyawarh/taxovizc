#!/usr/bin/env python3
"""
taxoviz_full.py
───────────────
Production CLI for full-dataset k-mer ordination + classifier agreement.

Usage
-----
    python taxoviz_full.py \\
        --fastq         /data/sample.fastq.gz \\
        --metaphlan-sam /data/metaphlan4.sam \\
        --krakenuniq    /data/krakenuniq.out \\
        --centrifuger   /data/centrifuger.out \\
        --taxonomy-dir  /data/ncbi_taxonomy/ \\
        --out-dir       /results/taxoviz/ \\
        [--k-values 4 6] \\
        [--pca-components 75] \\
        [--pca-variance 0.75] \\
        [--tsne-perplexity 30] \\
        [--tsne-metric manhattan] \\
        [--umap-neighbors 30] \\
        [--umap-min-dist 0.1] \\
        [--umap-metric manhattan] \\
        [--chunk-size 50000] \\
        [--n-jobs -1] \\
        [--seed 42] \\
        [--no-gpu]

Outputs (all written to --out-dir)
-------
    taxoviz_static.svg / .png   — 2×2 grid: t-SNE + UMAP × genus + species
    taxoviz_interactive.html    — Plotly WebGL interactive plot
    taxoviz_agreement_stats.csv — per-category counts and percentages
    taxoviz_embeddings.npz      — raw coordinates (tsne, umap) + read_ids
    taxoviz_metadata.parquet    — full per-read metadata DataFrame

Scalability
-----------
    ≤500k reads  : in-memory path (sklearn PCA, direct openTSNE, umap-learn CPU)
    >500k reads  : streaming path (HDF5 feature matrix, IncrementalPCA,
                   landmark openTSNE, cuML UMAP on GPU)
    Landmark t-SNE threshold: 100,000 landmarks
    GPU UMAP: requires cuML (RAPIDS); falls back to CPU umap-learn if unavailable
              or if --no-gpu is passed
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pickle
import bz2

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("taxoviz")

# ── taxoviz_core ─────────────────────────────────────────────────────────────
# taxoviz_core.py must be in the same directory as this script.
_HERE = Path(__file__).parent
sys.path.insert(0, str(_HERE))
import taxoviz_core as tc

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

LANDMARK_THRESHOLD = 100_000   # use landmark t-SNE above this many reads
STREAMING_THRESHOLD = 500_000  # use streaming/HDF5 path above this many reads


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: COUNT READS
# ─────────────────────────────────────────────────────────────────────────────

def count_reads(fastq_path: str) -> int:
    """Count reads in FASTQ without loading sequences."""
    log.info("Counting reads in %s ...", fastq_path)
    n = 0
    with tc.open_fastq(fastq_path) as fh:
        for line in fh:
            if line.startswith("@"):
                n += 1
    log.info("  Total reads: %s", f"{n:,}")
    return n


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: BUILD FEATURE MATRIX  (streaming → HDF5 for large datasets)
# ─────────────────────────────────────────────────────────────────────────────

def build_feature_matrix_streaming(
    fastq_path: str,
    hdf5_path: str,
    k_values: tuple,
    pseudocount: float,
    chunk_size: int,
    n_jobs: int,
    n_reads: int,
) -> tuple[list[str], int, np.ndarray, np.ndarray]:
    """
    Stream FASTQ in chunks, compute CLR k-mer + GC features per chunk,
    write to HDF5 dataset (n_reads × D).

    Also collects read_lens and gc_arr during the same pass so that
    sequences do NOT need to be re-read into RAM for metadata construction.
    (15M reads × 4.5 kb = ~67 GB RAM if held as strings — avoided here.)

    Returns (read_ids, feature_dim, read_lens, gc_arr).
    """
    from tqdm import tqdm

    # Feature dimension
    D = sum(tc.KMER_VOCAB_SIZE[k] for k in k_values) + 1  # +1 for GC

    log.info("Building feature matrix → %s  (D=%d)", hdf5_path, D)

    all_read_ids: list[str] = []
    all_read_lens: list[int] = []
    all_gc: list[float] = []
    n_written = 0

    with h5py.File(hdf5_path, "w") as hf:
        ds = hf.create_dataset(
            "features",
            shape=(n_reads, D),
            dtype=np.float32,
            chunks=(min(chunk_size, n_reads), D),
            compression="lzf",
        )

        pbar = tqdm(
            total=n_reads, unit="reads", desc="Vectorising",
            dynamic_ncols=True,
        )
        for chunk_ids, chunk_seqs in tc.iter_fastq_chunks(fastq_path, chunk_size):
            X_chunk = tc.build_feature_matrix(
                chunk_seqs, k_values=k_values,
                pseudocount=pseudocount, n_jobs=n_jobs,
            )
            end = n_written + len(chunk_ids)
            ds[n_written:end] = X_chunk
            all_read_ids.extend(chunk_ids)
            # Collect read lengths and GC while sequences are in memory
            all_read_lens.extend(len(s) for s in chunk_seqs)
            all_gc.extend(
                (s.count("G") + s.count("C")) / max(len(s), 1)
                for s in chunk_seqs
            )
            n_written = end
            pbar.update(len(chunk_ids))
        pbar.close()

    read_lens = np.array(all_read_lens, dtype=np.int32)
    gc_arr    = np.array(all_gc,        dtype=np.float32)
    log.info("  Feature matrix written: %s reads × %d features", f"{n_written:,}", D)
    return all_read_ids, D, read_lens, gc_arr


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: PCA
# ─────────────────────────────────────────────────────────────────────────────

def run_pca_incremental(
    hdf5_path: str,
    n_reads: int,
    pca_max_components: int,
    pca_variance_threshold: float,
    chunk_size: int,
    random_seed: int,
) -> np.ndarray:
    """
    IncrementalPCA from HDF5 feature matrix.
    Two passes: first to fit, second to transform.
    Returns X_pca (n_reads × n_components), float32.
    """
    from sklearn.decomposition import IncrementalPCA
    from tqdm import tqdm

    log.info("IncrementalPCA  (max_components=%d) ...", pca_max_components)
    ipca = IncrementalPCA(n_components=pca_max_components)

    with h5py.File(hdf5_path, "r") as hf:
        ds = hf["features"]

        # Pass 1: fit
        pbar = tqdm(range(0, n_reads, chunk_size), desc="PCA fit", unit="chunk",
                    dynamic_ncols=True)
        for start in pbar:
            end = min(start + chunk_size, n_reads)
            ipca.partial_fit(ds[start:end])
        pbar.close()

        # Determine number of components
        cumvar = np.cumsum(ipca.explained_variance_ratio_)
        n_components = int(np.searchsorted(cumvar, pca_variance_threshold) + 1)
        n_components = min(n_components, pca_max_components)
        log.info(
            "  %d PCs explain %.1f%% variance  (threshold %.0f%%)",
            n_components, 100 * cumvar[n_components - 1],
            100 * pca_variance_threshold,
        )

        # Pass 2: transform
        X_pca = np.empty((n_reads, n_components), dtype=np.float32)
        pbar = tqdm(range(0, n_reads, chunk_size), desc="PCA transform", unit="chunk",
                    dynamic_ncols=True)
        for start in pbar:
            end = min(start + chunk_size, n_reads)
            chunk_full = ipca.transform(ds[start:end])
            X_pca[start:end] = chunk_full[:, :n_components].astype(np.float32)
        pbar.close()

    log.info("  X_pca shape: %s", X_pca.shape)
    return X_pca


def run_pca_full(
    X: np.ndarray,
    pca_max_components: int,
    pca_variance_threshold: float,
    random_seed: int,
) -> np.ndarray:
    """Standard PCA for in-memory datasets (≤500k reads)."""
    from sklearn.decomposition import PCA

    log.info("PCA  (max_components=%d) ...", pca_max_components)
    pca = PCA(n_components=pca_max_components, random_state=random_seed)
    X_pca_full = pca.fit_transform(X)

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    n_components = int(np.searchsorted(cumvar, pca_variance_threshold) + 1)
    n_components = min(n_components, pca_max_components)
    log.info(
        "  %d PCs explain %.1f%% variance  (threshold %.0f%%)",
        n_components, 100 * cumvar[n_components - 1],
        100 * pca_variance_threshold,
    )
    return X_pca_full[:, :n_components].astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: t-SNE
# ─────────────────────────────────────────────────────────────────────────────

def run_tsne(
    X_pca: np.ndarray,
    n_reads: int,
    perplexity: int,
    metric: str,
    n_iter: int,
    n_jobs: int,
    random_seed: int,
    landmark_n: int = LANDMARK_THRESHOLD,
) -> np.ndarray:
    """
    openTSNE embedding.
    - ≤landmark_n reads: direct fit
    - >landmark_n reads: landmark fit on random subset, then transform remainder
    """
    from openTSNE import TSNE
    from tqdm import tqdm

    tsne = TSNE(
        perplexity=perplexity,
        metric=metric,
        n_iter=n_iter,
        n_jobs=n_jobs,
        random_state=random_seed,
        verbose=False,
    )

    if n_reads <= landmark_n:
        log.info("t-SNE: direct fit  (n=%s) ...", f"{n_reads:,}")
        t0 = time.time()
        X_tsne = np.array(tsne.fit(X_pca))
        log.info("  Done in %.1f min", (time.time() - t0) / 60)
        return X_tsne

    # Landmark approach
    log.info(
        "t-SNE: landmark fit  (landmarks=%s, total=%s) ...",
        f"{landmark_n:,}", f"{n_reads:,}",
    )
    rng = np.random.default_rng(random_seed)
    landmark_idx = rng.choice(n_reads, size=landmark_n, replace=False)
    landmark_idx.sort()
    rest_idx = np.setdiff1d(np.arange(n_reads), landmark_idx)

    t0 = time.time()
    log.info("  Fitting landmarks ...")
    embedding_landmarks = tsne.fit(X_pca[landmark_idx])
    log.info("  Landmark fit done in %.1f min", (time.time() - t0) / 60)

    # Transform remaining reads in chunks to show progress
    X_tsne = np.empty((n_reads, 2), dtype=np.float32)
    X_tsne[landmark_idx] = np.array(embedding_landmarks)

    REST_CHUNK = 50_000
    log.info("  Projecting %s remaining reads ...", f"{len(rest_idx):,}")
    pbar = tqdm(range(0, len(rest_idx), REST_CHUNK), desc="t-SNE project",
                unit="chunk", dynamic_ncols=True)
    for start in pbar:
        end = min(start + REST_CHUNK, len(rest_idx))
        idx_chunk = rest_idx[start:end]
        projected = embedding_landmarks.transform(X_pca[idx_chunk])
        X_tsne[idx_chunk] = np.array(projected)
    pbar.close()

    log.info("  t-SNE total time: %.1f min", (time.time() - t0) / 60)
    return X_tsne


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: UMAP
# ─────────────────────────────────────────────────────────────────────────────

def run_umap(
    X_pca: np.ndarray,
    n_neighbors: int,
    min_dist: float,
    metric: str,
    random_seed: int,
    use_gpu: bool,
) -> np.ndarray:
    """
    UMAP embedding.
    Tries cuML GPU UMAP first (if use_gpu=True), falls back to umap-learn CPU.
    """
    if use_gpu:
        try:
            from cuml.manifold import UMAP as cuUMAP
            log.info("UMAP: cuML GPU  (n_neighbors=%d, min_dist=%.2f) ...",
                     n_neighbors, min_dist)
            t0 = time.time()
            reducer = cuUMAP(
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                metric=metric,
                random_state=random_seed,
                verbose=True,
            )
            X_umap = reducer.fit_transform(X_pca)
            if hasattr(X_umap, "to_numpy"):
                X_umap = X_umap.to_numpy()
            log.info("  cuML UMAP done in %.1f min", (time.time() - t0) / 60)
            return X_umap.astype(np.float32)
        except ImportError:
            log.warning("cuML not available — falling back to CPU umap-learn")
        except Exception as e:
            log.warning("cuML UMAP failed (%s) — falling back to CPU", e)

    import umap as umap_learn
    log.info("UMAP: CPU umap-learn  (n_neighbors=%d, min_dist=%.2f) ...",
             n_neighbors, min_dist)
    t0 = time.time()
    reducer = umap_learn.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_seed,
        verbose=True,
        low_memory=True,
    )
    X_umap = reducer.fit_transform(X_pca).astype(np.float32)
    log.info("  CPU UMAP done in %.1f min", (time.time() - t0) / 60)
    return X_umap


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: SAVE OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────

def save_static_figure(
    out_dir: Path,
    X_tsne: np.ndarray,
    X_umap: np.ndarray,
    meta: pd.DataFrame,
    tsne_perplexity: int,
    umap_n_neighbors: int,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    log.info("Saving static figure ...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(
        "Read-level k-mer ordination coloured by classifier agreement",
        fontsize=13, fontweight="bold", y=1.01,
    )

    panels = [
        (axes[0, 0], X_tsne, "genus",
         f"t-SNE  (genus, perplexity={tsne_perplexity})"),
        (axes[0, 1], X_tsne, "species",
         f"t-SNE  (species, perplexity={tsne_perplexity})"),
        (axes[1, 0], X_umap, "genus",
         f"UMAP  (genus, n_neighbors={umap_n_neighbors})"),
        (axes[1, 1], X_umap, "species",
         f"UMAP  (species, n_neighbors={umap_n_neighbors})"),
    ]
    for ax, coords, level, title in panels:
        tc.scatter_ordination(
            ax, coords, meta[f"agreement_{level}"],
            title=title, xlabel="Dim 1", ylabel="Dim 2",
            alpha=0.3, s=1.5,
        )

    handles = tc.make_legend_handles()
    fig.legend(handles=handles, loc="lower center", ncol=4,
               fontsize=9, framealpha=0.8, bbox_to_anchor=(0.5, -0.04))
    plt.tight_layout()

    for ext in ("svg", "png"):
        path = out_dir / f"taxoviz_static.{ext}"
        fig.savefig(path, bbox_inches="tight", dpi=150)
        log.info("  Saved: %s", path)
    plt.close(fig)


def save_interactive_html(
    out_dir: Path,
    X_tsne: np.ndarray,
    X_umap: np.ndarray,
    meta: pd.DataFrame,
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    log.info("Saving interactive HTML ...")

    def _make_traces(coords, agreement_series, level, show_legend):
        traces = []
        for cat in tc.CATEGORIES:
            mask = (agreement_series == cat).values
            if mask.sum() == 0:
                continue
            rid_list = meta.index[mask].tolist()
            traces.append(go.Scattergl(
                x=coords[mask, 0], y=coords[mask, 1],
                mode="markers",
                marker=dict(size=2, color=tc.COLORS[cat], opacity=0.5),
                name=tc.LABEL_NAMES[cat],
                text=[
                    f"read: {rid}<br>len: {meta.loc[rid,'read_len']}<br>"
                    f"gc: {meta.loc[rid,'gc']:.2f}<br>"
                    f"M: {meta.loc[rid,f'metaphlan_{level}_name']}<br>"
                    f"K: {meta.loc[rid,f'krakenuniq_{level}_name']}<br>"
                    f"C: {meta.loc[rid,f'centrifuger_{level}_name']}"
                    for rid in rid_list
                ],
                hovertemplate="%{text}<extra></extra>",
                legendgroup=cat,
                showlegend=show_legend,
            ))
        return traces

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            "t-SNE genus", "t-SNE species",
            "UMAP genus",  "UMAP species",
        ],
    )
    for row, coords in enumerate([X_tsne, X_umap], start=1):
        for col, level in enumerate(["genus", "species"], start=1):
            show = (row == 1 and col == 1)
            for trace in _make_traces(coords, meta[f"agreement_{level}"], level, show):
                fig.add_trace(trace, row=row, col=col)

    fig.update_layout(
        title="TaxoViz — k-mer ordination × classifier agreement",
        height=900, width=1200,
        legend=dict(itemsizing="constant"),
    )
    path = out_dir / "taxoviz_interactive.html"
    fig.write_html(str(path))
    log.info("  Saved: %s", path)


def save_stats_csv(out_dir: Path, meta: pd.DataFrame) -> None:
    rows = []
    for level in ("genus", "species"):
        vc = meta[f"agreement_{level}"].value_counts()
        for cat in tc.CATEGORIES:
            n = vc.get(cat, 0)
            rows.append({
                "level":    level,
                "category": cat,
                "label":    tc.LABEL_NAMES[cat],
                "count":    n,
                "pct":      round(100 * n / len(meta), 2),
            })
    df = pd.DataFrame(rows)
    path = out_dir / "taxoviz_agreement_stats.csv"
    df.to_csv(path, index=False)
    log.info("  Saved: %s", path)


def save_embeddings(
    out_dir: Path,
    read_ids: list[str],
    X_tsne: np.ndarray,
    X_umap: np.ndarray,
) -> None:
    path = out_dir / "taxoviz_embeddings.npz"
    np.savez_compressed(
        path,
        read_ids=np.array(read_ids, dtype=object),
        tsne=X_tsne,
        umap=X_umap,
    )
    log.info("  Saved: %s", path)


def save_metadata(out_dir: Path, meta: pd.DataFrame) -> None:
    path = out_dir / "taxoviz_metadata.parquet"
    meta.to_parquet(path)
    log.info("  Saved: %s", path)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(args: argparse.Namespace) -> None:
    t_start = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("TaxoViz full-dataset pipeline")
    log.info("  FASTQ        : %s", args.fastq)
    log.info("  Output dir   : %s", out_dir)
    log.info("  k values     : %s", args.k_values)
    log.info("  GPU UMAP     : %s", not args.no_gpu)
    log.info("=" * 60)

    k_values = tuple(args.k_values)

    # ── Step 1: Count reads ──────────────────────────────────────────────────
    n_reads = count_reads(args.fastq)
    streaming = n_reads > STREAMING_THRESHOLD
    log.info("Mode: %s  (threshold %s)",
             "STREAMING" if streaming else "IN-MEMORY",
             f"{STREAMING_THRESHOLD:,}")

    # ── Step 2: Feature matrix ───────────────────────────────────────────────
    if streaming:
        hdf5_path = str(out_dir / "taxoviz_features.h5")
        read_ids, D, read_lens_arr, gc_arr = build_feature_matrix_streaming(
            args.fastq, hdf5_path, k_values,
            args.clr_pseudocount, args.chunk_size, args.n_jobs, n_reads,
        )
    else:
        log.info("Loading all reads into memory ...")
        from tqdm import tqdm
        read_ids, sequences = tc.read_fastq_all(args.fastq)
        log.info("  Vectorising %s reads ...", f"{len(read_ids):,}")
        X = tc.build_feature_matrix(
            sequences, k_values=k_values,
            pseudocount=args.clr_pseudocount, n_jobs=args.n_jobs,
        )
        log.info("  Feature matrix: %s  %.1f MB", X.shape, X.nbytes / 1e6)

    # ── Step 3: Parse classifiers ────────────────────────────────────────────
    log.info("Parsing classifier outputs ...")
    metphlan_taxdb = pickle.load(bz2.open('/mnt/raid0/Databases/DDTdb/MetaPhlAn4/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl', 'r'))
    mark2taxid = tc.build_marker_to_taxid(metphlan_taxdb)
    read_ids_set = set(read_ids)
    labels_m = tc.parse_metaphlan_sam(args.metaphlan_sam,  read_ids_set, mark2taxid)
    labels_k = tc.parse_krakenuniq(   args.krakenuniq,     read_ids_set)
    labels_c = tc.parse_centrifuger(  args.centrifuger,    read_ids_set)

    def _n_classified(d):
        return sum(1 for v in d.values() if v is not None and v > 0)

    log.info("  MetaPhlAn4  %s / %s", f"{_n_classified(labels_m):,}", f"{n_reads:,}")
    log.info("  KrakenUniq  %s / %s", f"{_n_classified(labels_k):,}", f"{n_reads:,}")
    log.info("  Centrifuger %s / %s", f"{_n_classified(labels_c):,}", f"{n_reads:,}")

    # ── Step 4: Taxonomy normalisation ──────────────────────────────────────
    log.info("Loading taxonomy from %s ...", args.taxonomy_dir)
    taxdb = tc.load_taxonomy(args.taxonomy_dir)

    log.info("Building metadata DataFrame (%s reads) ...", f"{n_reads:,}")
    if streaming:
        # Sequences are NOT re-read into RAM (would be ~67 GB for 15M reads).
        # read_lens_arr and gc_arr were collected during the vectorisation pass.
        meta = tc.build_metadata_df(
            read_ids, None,
            labels_m, labels_k, labels_c,
            taxdb,
            read_lens=read_lens_arr,
            gc_arr=gc_arr,
        )
    else:
        meta = tc.build_metadata_df(
            read_ids, sequences,
            labels_m, labels_k, labels_c,
            taxdb,
        )
    log.info("  Metadata shape: %s", meta.shape)

    # ── Step 5: PCA ──────────────────────────────────────────────────────────
    if streaming:
        X_pca = run_pca_incremental(
            hdf5_path, n_reads,
            args.pca_components, args.pca_variance,
            args.chunk_size, args.seed,
        )
    else:
        X_pca = run_pca_full(X, args.pca_components, args.pca_variance, args.seed)

    # ── Step 6: t-SNE ────────────────────────────────────────────────────────
    X_tsne = run_tsne(
        X_pca, n_reads,
        perplexity=args.tsne_perplexity,
        metric=args.tsne_metric,
        n_iter=args.tsne_n_iter,
        n_jobs=args.n_jobs,
        random_seed=args.seed,
        landmark_n=LANDMARK_THRESHOLD,
    )

    # ── Step 7: UMAP ─────────────────────────────────────────────────────────
    X_umap = run_umap(
        X_pca,
        n_neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        metric=args.umap_metric,
        random_seed=args.seed,
        use_gpu=not args.no_gpu,
    )

    # ── Step 8: Save outputs ─────────────────────────────────────────────────
    log.info("Saving outputs to %s ...", out_dir)
    save_static_figure(out_dir, X_tsne, X_umap, meta,
                       args.tsne_perplexity, args.umap_neighbors)
    save_interactive_html(out_dir, X_tsne, X_umap, meta)
    save_stats_csv(out_dir, meta)
    save_embeddings(out_dir, read_ids, X_tsne, X_umap)
    save_metadata(out_dir, meta)

    # Clean up HDF5 scratch file
    if streaming:
        Path(hdf5_path).unlink(missing_ok=True)
        log.info("  Removed HDF5 scratch: %s", hdf5_path)

    elapsed = (time.time() - t_start) / 60
    log.info("=" * 60)
    log.info("Pipeline complete in %.1f min", elapsed)
    log.info("Outputs in: %s", out_dir)
    log.info("=" * 60)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="TaxoViz: k-mer ordination + classifier agreement (full dataset)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required
    req = p.add_argument_group("Required")
    req.add_argument("--fastq",         required=True, help="Input FASTQ (.gz ok)")
    req.add_argument("--metaphlan-sam", required=True, dest="metaphlan_sam",
                     help="MetaPhlAn4 --samout SAM file")
    req.add_argument("--krakenuniq",    required=True,
                     help="KrakenUniq per-read output")
    req.add_argument("--centrifuger",   required=True,
                     help="Centrifuger per-read output")
    req.add_argument("--taxonomy-dir",  required=True, dest="taxonomy_dir",
                     help="Directory with nodes.dmp and names.dmp")
    req.add_argument("--out-dir",       required=True, dest="out_dir",
                     help="Output directory (created if absent)")

    # k-mer
    km = p.add_argument_group("k-mer features")
    km.add_argument("--k-values",       nargs="+", type=int, default=[4, 6],
                    dest="k_values",    help="k values (concatenated)")
    km.add_argument("--clr-pseudocount", type=float, default=1.0,
                    dest="clr_pseudocount", help="CLR pseudocount δ")

    # PCA
    pca = p.add_argument_group("PCA")
    pca.add_argument("--pca-components", type=int, default=75,
                     dest="pca_components", help="Max PCA components")
    pca.add_argument("--pca-variance",   type=float, default=0.75,
                     dest="pca_variance", help="Variance threshold for PC selection")

    # t-SNE
    tsne = p.add_argument_group("t-SNE")
    tsne.add_argument("--tsne-perplexity", type=int,   default=30,
                      dest="tsne_perplexity")
    tsne.add_argument("--tsne-metric",     type=str,   default="manhattan",
                      dest="tsne_metric")
    tsne.add_argument("--tsne-n-iter",     type=int,   default=750,
                      dest="tsne_n_iter")

    # UMAP
    um = p.add_argument_group("UMAP")
    um.add_argument("--umap-neighbors", type=int,   default=30,
                    dest="umap_neighbors")
    um.add_argument("--umap-min-dist",  type=float, default=0.1,
                    dest="umap_min_dist")
    um.add_argument("--umap-metric",    type=str,   default="manhattan",
                    dest="umap_metric")

    # Compute
    comp = p.add_argument_group("Compute")
    comp.add_argument("--chunk-size", type=int, default=50_000,
                      dest="chunk_size",
                      help="Reads per chunk for streaming/IncrementalPCA")
    comp.add_argument("--n-jobs",     type=int, default=-1,
                      dest="n_jobs",  help="Parallel jobs (-1 = all cores)")
    comp.add_argument("--seed",       type=int, default=42)
    comp.add_argument("--no-gpu",     action="store_true", dest="no_gpu",
                      help="Disable GPU UMAP (force CPU umap-learn)")

    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)
