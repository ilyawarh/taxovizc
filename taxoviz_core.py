"""
taxoviz_core.py
───────────────
Shared functions used by both the interactive notebook (taxoviz_notebook.ipynb)
and the full-dataset CLI tool (taxoviz_full.py).

Covers:
  - Canonical k-mer index building + numeric lookup arrays (Numba-compatible)
  - Per-read k-mer count vectorisation (Numba JIT with correct canonical mapping)
  - CLR transformation with pseudocount
  - GC-content computation
  - Feature matrix assembly (CLR k-mers + GC z-score)
  - Classifier output parsers (MetaPhlAn4 SAM, KrakenUniq, Centrifuger)
  - NCBI taxonomy normalisation (taxopy)
  - Agreement pattern computation (8 categories × genus/species)
  - Colour palette and plot helpers
"""

from __future__ import annotations

import bz2
import gzip
import re
import warnings
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  CANONICAL K-MER INDEX
# ─────────────────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def build_canonical_kmer_index(k: int) -> Tuple[Dict[str, int], int]:
    """
    Build {canonical_kmer_string: position} for all canonical k-mers of length k.
    Canonical = lexicographically smaller of (kmer, revcomp(kmer)).
    Positions are assigned in the order canonical k-mers are first encountered
    when iterating product("ACGT", repeat=k).

    Returns
    -------
    index : dict  {kmer_str -> int}
    vocab_size : int
    """
    index: Dict[str, int] = {}
    idx = 0
    for kmer_tuple in product("ACGT", repeat=k):
        kmer = "".join(kmer_tuple)
        rc = reverse_complement(kmer)
        canonical = kmer if kmer <= rc else rc
        if canonical not in index:
            index[canonical] = idx
            idx += 1
    return index, idx


def build_numeric_to_canonical_lookup(k: int) -> np.ndarray:
    """
    Build a lookup array of shape (4^k,) mapping numeric base4 k-mer index
    to canonical position (as defined by build_canonical_kmer_index).

    Numeric base4 index: A=0, C=1, G=2, T=3 → kmer encoded as base-4 number.
    This encoding preserves lexicographic order for the ACGT alphabet, so
    min(fwd_numeric, rc_numeric) == position of the lexicographically smaller
    string — verified exhaustively for k=4 and k=6.

    Positions for non-canonical k-mers are set to -1 (should never be accessed
    since we always take min(fwd, rc) before lookup).

    Returns
    -------
    lookup : np.ndarray  shape (4^k,), dtype int32
    """
    enc = {"A": 0, "C": 1, "G": 2, "T": 3}

    def base4(kmer: str) -> int:
        idx = 0
        for b in kmer:
            idx = idx * 4 + enc[b]
        return idx

    # Build canonical list in the same order as build_canonical_kmer_index
    canonical_numeric: List[int] = []
    seen: set = set()
    for kmer_tuple in product("ACGT", repeat=k):
        kmer = "".join(kmer_tuple)
        rc = reverse_complement(kmer)
        canonical = kmer if kmer <= rc else rc
        num = base4(canonical)
        if num not in seen:
            seen.add(num)
            canonical_numeric.append(num)

    lookup = np.full(4 ** k, -1, dtype=np.int32)
    for pos, num in enumerate(canonical_numeric):
        lookup[num] = pos
    return lookup


# Pre-build at import time (fast, <1 s for k=4,6)
KMER_INDEX: Dict[int, Dict[str, int]] = {}
KMER_VOCAB_SIZE: Dict[int, int] = {}
KMER_NUMERIC_LOOKUP: Dict[int, np.ndarray] = {}

for _k in (4, 6):
    KMER_INDEX[_k], KMER_VOCAB_SIZE[_k] = build_canonical_kmer_index(_k)
    KMER_NUMERIC_LOOKUP[_k] = build_numeric_to_canonical_lookup(_k)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  K-MER VECTORISATION
# ─────────────────────────────────────────────────────────────────────────────

def _vectorise_read_python(seq: str, k: int) -> np.ndarray:
    """
    Pure-Python canonical k-mer count vector for one read.
    Returns np.ndarray shape (vocab_size,), dtype float32.
    Skips k-mers containing N.
    """
    index = KMER_INDEX[k]
    vocab_size = KMER_VOCAB_SIZE[k]
    vec = np.zeros(vocab_size, dtype=np.float32)
    n = len(seq)
    if n < k:
        return vec
    for i in range(n - k + 1):
        kmer = seq[i: i + k]
        if "N" in kmer:
            continue
        rc = reverse_complement(kmer)
        canonical = kmer if kmer <= rc else rc
        idx = index.get(canonical)
        if idx is not None:
            vec[idx] += 1.0
    return vec


def _make_numba_vectoriser(k: int):
    """
    Return a Numba-JIT compiled vectoriser for a specific k.
    Falls back to pure Python if numba is not available.

    The Numba version uses:
      - BASE_ENC[byte] → 0..3 for ACGT, -1 for N/other
      - COMP_ENC[base] → complement in ACGT encoding (A=0↔T=3, C=1↔G=2)
      - NUMERIC_LOOKUP[canonical_numeric_idx] → canonical position in vocab
        (size 4^k; positions for non-canonical k-mers are -1 but never accessed)

    Correctness note: for the ACGT alphabet encoded as A=0,C=1,G=2,T=3,
    numeric base4 order is identical to lexicographic order, so
    min(fwd_numeric, rc_numeric) correctly identifies the canonical k-mer.
    """
    try:
        import numba as nb

        BASE_ENC = np.full(256, -1, dtype=np.int32)
        for byte_val, enc_val in zip(b"ACGT", [0, 1, 2, 3]):
            BASE_ENC[byte_val] = enc_val
        # Complement in ACGT encoding: A(0)↔T(3), C(1)↔G(2)
        COMP_ENC = np.array([3, 2, 1, 0], dtype=np.int32)

        vocab_size = KMER_VOCAB_SIZE[k]
        numeric_lookup = KMER_NUMERIC_LOOKUP[k].copy()  # shape (4^k,), int32

        @nb.njit(cache=True)
        def _count(
            seq_bytes: bytes,
            k: int,
            vocab_size: int,
            base_enc: np.ndarray,
            comp_enc: np.ndarray,
            numeric_lookup: np.ndarray,
        ) -> np.ndarray:
            vec = np.zeros(vocab_size, dtype=np.float32)
            n = len(seq_bytes)
            for i in range(n - k + 1):
                # Compute forward numeric index
                fwd = 0
                valid = True
                for j in range(k):
                    b = base_enc[seq_bytes[i + j]]
                    if b < 0:
                        valid = False
                        break
                    fwd = fwd * 4 + b
                if not valid:
                    continue

                # Compute reverse-complement numeric index (iterate in reverse)
                rev = 0
                for j in range(k - 1, -1, -1):
                    b = base_enc[seq_bytes[i + j]]
                    rev = rev * 4 + comp_enc[b]

                # Canonical = numerically smaller (equivalent to lexicographic min)
                canonical_num = fwd if fwd <= rev else rev
                pos = numeric_lookup[canonical_num]
                # pos == -1 only for non-canonical entries, which never occur here
                if pos >= 0:
                    vec[pos] += 1.0
            return vec

        # Warm-up JIT compilation
        _count(b"ACGTACGT", k, vocab_size, BASE_ENC, COMP_ENC, numeric_lookup)

        def vectorise_one(seq: str) -> np.ndarray:
            return _count(
                seq.encode("ascii"), k, vocab_size,
                BASE_ENC, COMP_ENC, numeric_lookup,
            )

        return vectorise_one

    except ImportError:
        def vectorise_one(seq: str) -> np.ndarray:
            return _vectorise_read_python(seq, k)
        return vectorise_one


# ─────────────────────────────────────────────────────────────────────────────
# 3.  CLR TRANSFORMATION
# ─────────────────────────────────────────────────────────────────────────────

def clr_transform(counts: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    """
    Centered log-ratio transformation with pseudocount.

    CLR(x_i) = log((x_i + δ) / geometric_mean(x + δ))

    Parameters
    ----------
    counts      : np.ndarray  shape (n, vocab_size), raw k-mer counts
    pseudocount : float       δ added before log (default 1.0)

    Returns
    -------
    np.ndarray  shape (n, vocab_size), dtype float32
    """
    x = counts.astype(np.float64) + pseudocount
    log_x = np.log(x)
    log_gm = log_x.mean(axis=1, keepdims=True)   # geometric mean in log space
    return (log_x - log_gm).astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 4.  GC-CONTENT
# ─────────────────────────────────────────────────────────────────────────────

def gc_content(sequences: List[str]) -> np.ndarray:
    """GC fraction per sequence. Returns shape (n,), dtype float32."""
    return np.array(
        [(s.count("G") + s.count("C")) / max(len(s), 1) for s in sequences],
        dtype=np.float32,
    )


def gc_zscore(gc: np.ndarray) -> np.ndarray:
    """Z-score normalise GC array. Returns float32 column vector (n, 1)."""
    mu, sigma = gc.mean(), gc.std()
    if sigma < 1e-9:
        return np.zeros((len(gc), 1), dtype=np.float32)
    return ((gc - mu) / sigma).reshape(-1, 1).astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  FULL FEATURE MATRIX  (CLR k-mers + GC z-score)
# ─────────────────────────────────────────────────────────────────────────────

def build_feature_matrix(
    sequences: List[str],
    k_values: Tuple[int, ...] = (4, 6),
    pseudocount: float = 1.0,
    n_jobs: int = -1,
) -> np.ndarray:
    """
    Build the full feature matrix for a list of sequences.

    Pipeline:
      For each k: canonical k-mer counts → CLR transform
      Concatenate CLR vectors across all k values
      Append GC z-score as final column

    Result shape: (n, sum(vocab_sizes_k) + 1)
    For k=(4,6): (n, 136 + 2080 + 1) = (n, 2217)

    Parameters
    ----------
    sequences   : list of str
    k_values    : tuple of int  (default (4, 6))
    pseudocount : float
    n_jobs      : int  (-1 = all cores via joblib)

    Returns
    -------
    np.ndarray  shape (n, D+1), dtype float32
    """
    import joblib

    clr_parts = []
    for k in k_values:
        vectorise_one = _make_numba_vectoriser(k)
        counts = np.array(
            joblib.Parallel(n_jobs=n_jobs, prefer="threads")(
                joblib.delayed(vectorise_one)(seq) for seq in sequences
            ),
            dtype=np.float32,
        )
        clr_parts.append(clr_transform(counts, pseudocount))

    clr_matrix = np.concatenate(clr_parts, axis=1)
    gc_z = gc_zscore(gc_content(sequences))
    return np.concatenate([clr_matrix, gc_z], axis=1)


# ─────────────────────────────────────────────────────────────────────────────
# 6.  FASTQ READER
# ─────────────────────────────────────────────────────────────────────────────

def open_fastq(path: str):
    """Open FASTQ, FASTQ.gz, or FASTQ.bz2 transparently."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    elif p.endswith(".bz2"):
        return bz2.open(p, "rt")
    return open(p, "r")


def read_fastq_all(path: str) -> Tuple[List[str], List[str]]:
    """
    Read all reads from a FASTQ file into memory.
    Returns (read_ids, sequences).
    Suitable for small datasets (≤500k reads).
    """
    read_ids, sequences = [], []
    with open_fastq(path) as fh:
        for record in SeqIO.parse(fh, "fastq"):
            read_ids.append(record.id)
            sequences.append(str(record.seq).upper())
    return read_ids, sequences


def iter_fastq_chunks(path: str, chunk_size: int = 50_000):
    """
    Generator: yield (read_ids_chunk, sequences_chunk) in chunks.
    Used for streaming large files without loading everything into RAM.
    """
    read_ids_chunk, seqs_chunk = [], []
    with open_fastq(path) as fh:
        for record in SeqIO.parse(fh, "fastq"):
            read_ids_chunk.append(record.id)
            seqs_chunk.append(str(record.seq).upper())
            if len(seqs_chunk) >= chunk_size:
                yield read_ids_chunk, seqs_chunk
                read_ids_chunk, seqs_chunk = [], []
    if seqs_chunk:
        yield read_ids_chunk, seqs_chunk


# ─────────────────────────────────────────────────────────────────────────────
# 7.  CLASSIFIER OUTPUT PARSERS
# ─────────────────────────────────────────────────────────────────────────────
import pickle
import bz2
from typing import Dict, Optional

def build_marker_to_taxid(metphlan_taxdb):
    marker_to_taxid = {}

    for lineage, info in metphlan_taxdb["taxonomy"].items():

        marker = lineage.split("|t__")[-1]   # SGB27704 или MVxxx...

        tax_path = info[0]
        taxid = int(tax_path.strip("|").split("|")[-1])

        marker_to_taxid[marker] = taxid

    return marker_to_taxid

def parse_metaphlan_sam(
    sam_path: str, read_ids_set: set,
    marker_to_taxid: dict,
) -> Dict[str, Optional[int]]:
    """
    Parse MetaPhlAn4 SAM output (--samout).
    Mapped reads → extract taxid from marker reference name.
    Unmapped reads (SAM flag & 4) → None.

    Supported marker name formats:
      'UniRef90_XXXXX__taxid__NNNN'
      'GeneID:XXXXX|taxid|NNNN|...'
    """
    result: Dict[str, Optional[int]] = {r: None for r in read_ids_set}

    sgb_pattern = re.compile(r"(SGB(\d+))$", re.IGNORECASE)

    p = str(sam_path)
    opener = bz2.open if p.endswith(".bz2") else (
        gzip.open if p.endswith(".gz") else open
    )
    with opener(p, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            read_id = parts[0]
            if read_id not in read_ids_set:
                continue
            if int(parts[1]) & 4:
                continue  # unmapped

            m = sgb_pattern.search(parts[2])
            if m:
                result[read_id] = m.group(1)
            
            taxname = result[read_id]
            taxid = marker_to_taxid.get(taxname)
            if taxid is not None:
                result[read_id] = int(taxid)
            else:
                result[read_id] = 0
    return result


def parse_krakenuniq(
    out_path: str, read_ids_set: set
) -> Dict[str, Optional[int]]:
    """
    Parse KrakenUniq per-read output.
    Format: C/U  read_id  taxid  length  kmer_hits
    """
    result: Dict[str, Optional[int]] = {r: None for r in read_ids_set}
    opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(out_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            status, read_id, taxid_str = parts[0], parts[1], parts[2]
            if read_id not in read_ids_set or status != "C":
                continue
            try:
                result[read_id] = int(taxid_str)
            except ValueError:
                pass
    return result


def parse_centrifuger(
    out_path: str, read_ids_set: set
) -> Dict[str, Optional[int]]:
    """
    Parse Centrifuger per-read output.
    Format: readID  seqID  taxID  score  2ndBestScore  hitLength  queryLength  numMatches
    taxID == 0 means unclassified → stored as None.
    """
    result: Dict[str, Optional[int]] = {r: None for r in read_ids_set}
    opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(out_path, "rt") as fh:
        header_skipped = False
        for line in fh:
            if not header_skipped:
                header_skipped = True
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            read_id, taxid_str = parts[0], parts[2]
            if read_id not in read_ids_set:
                continue
            try:
                taxid = int(taxid_str)
                result[read_id] = taxid if taxid > 0 else None
            except ValueError:
                pass
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 8.  TAXONOMY NORMALISATION
# ─────────────────────────────────────────────────────────────────────────────

def load_taxonomy(taxonomy_dir: str):
    """Load NCBI taxonomy via taxopy. Returns TaxDb object."""
    import taxopy
    return taxopy.TaxDb(
        nodes_dmp=str(Path(taxonomy_dir) / "nodes.dmp"),
        names_dmp=str(Path(taxonomy_dir) / "names.dmp"),
        keep_files=True,
    )


def get_rank_taxid(
    taxid: Optional[int], rank: str, taxdb
) -> Optional[int]:
    """
    Walk up the taxonomy tree from taxid to find the ancestor at `rank`.
    Returns taxid at that rank, or None if not found.
    """
    if taxid is None or taxid <= 0:
        return None
    import taxopy
    try:
        taxon = taxopy.Taxon(int(taxid), taxdb)
    except Exception:
        return None
    if taxon.rank == rank:
        return taxon.taxid
    for ancestor_id in taxon.taxid_lineage:
        try:
            ancestor = taxopy.Taxon(ancestor_id, taxdb)
            if ancestor.rank == rank:
                return ancestor.taxid
        except Exception:
            continue
    return None


def build_taxid_name_cache(taxids: set, taxdb) -> Dict[int, str]:
    """Build {taxid: name} cache for all unique taxids in one pass."""
    import taxopy
    cache: Dict[int, str] = {}
    for tid in taxids:
        if tid is None or tid <= 0:
            continue
        try:
            cache[int(tid)] = taxopy.Taxon(int(tid), taxdb).name
        except Exception:
            cache[int(tid)] = f"taxid:{tid}"
    return cache


def taxid_to_name(tid, cache: Dict[int, str]) -> str:
    """Resolve taxid to name using pre-built cache. Safe for None/NaN/0."""
    if tid is None:
        return "unclassified"
    try:
        tid_int = int(tid)
    except (ValueError, TypeError):
        return "unclassified"
    if tid_int <= 0:
        return "unclassified"
    return cache.get(tid_int, f"taxid:{tid_int}")


# ─────────────────────────────────────────────────────────────────────────────
# 9.  AGREEMENT PATTERN COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────

TOOLS = ("metaphlan", "krakenuniq", "centrifuger")
TOOL_ABBR = {"metaphlan": "M", "krakenuniq": "K", "centrifuger": "C"}

CATEGORIES = ["none", "M_only", "K_only", "C_only", "MK", "MC", "KC", "MKC"]

# Colorblind-friendly palette (Wong 2011)
COLORS = {
    "none":   "#AAAAAA",
    "M_only": "#E69F00",
    "K_only": "#56B4E9",
    "C_only": "#009E73",
    "MK":     "#F0E442",
    "MC":     "#D55E00",
    "KC":     "#0072B2",
    "MKC":    "#CC79A7",
}

LABEL_NAMES = {
    "none":   "None classified",
    "M_only": "MetaPhlAn4 only",
    "K_only": "KrakenUniq only",
    "C_only": "Centrifuger only",
    "MK":     "MetaPhlAn4 + KrakenUniq",
    "MC":     "MetaPhlAn4 + Centrifuger",
    "KC":     "KrakenUniq + Centrifuger",
    "MKC":    "All three agree",
}

_ABBR_TO_CAT = {
    "C": "C_only", "K": "K_only", "M": "M_only",
    "CK": "KC", "CM": "MC", "KM": "MK",
    "CKM": "MKC",
}


def compute_agreement_row(
    metaphlan_taxid: Optional[int],
    krakenuniq_taxid: Optional[int],
    centrifuger_taxid: Optional[int],
) -> str:
    """
    Compute agreement category for one read given three taxid values.

    Logic:
      1. Collect tools that gave a non-None, non-zero taxid
      2. Find consensus taxid (most common among active tools)
      3. Agreement = set of tools that match the consensus

    Edge case: if all three tools give different taxids, Counter picks the
    first-inserted key (Python 3.7+ dict order = metaphlan, krakenuniq,
    centrifuger), so metaphlan's taxid becomes the "consensus" and the
    result is "M_only". This is deterministic and acceptable for visualisation.
    """
    labels = {
        "metaphlan":   metaphlan_taxid,
        "krakenuniq":  krakenuniq_taxid,
        "centrifuger": centrifuger_taxid,
    }
    active = {t: v for t, v in labels.items() if v is not None and v > 0}
    if not active:
        return "none"
    consensus_taxid = Counter(active.values()).most_common(1)[0][0]
    agreeing = sorted(TOOL_ABBR[t] for t, v in active.items() if v == consensus_taxid)
    return _ABBR_TO_CAT.get("".join(agreeing), "none")


def build_metadata_df(
    read_ids: List[str],
    sequences: Optional[List[str]],
    labels_metaphlan: Dict[str, Optional[int]],
    labels_krakenuniq: Dict[str, Optional[int]],
    labels_centrifuger: Dict[str, Optional[int]],
    taxdb,
    read_lens: Optional[np.ndarray] = None,
    gc_arr: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Build the full metadata DataFrame with normalised taxids,
    agreement patterns at genus and species level, read lengths, and GC.

    Parameters
    ----------
    read_ids   : list of str
    sequences  : list of str, OR None if read_lens and gc_arr are provided
                 (streaming path: sequences are not held in RAM)
    read_lens  : optional np.ndarray shape (n,) — pre-computed read lengths
                 (used when sequences=None, e.g. streaming path)
    gc_arr     : optional np.ndarray shape (n,) — pre-computed GC fractions
                 (used when sequences=None, e.g. streaming path)

    Returns pd.DataFrame indexed by read_id.
    """
    if sequences is None and (read_lens is None or gc_arr is None):
        raise ValueError("Either sequences or both read_lens and gc_arr must be provided")

    # Collect all raw taxids for batch normalisation
    all_raw: set = set()
    for labels in (labels_metaphlan, labels_krakenuniq, labels_centrifuger):
        all_raw.update(v for v in labels.values() if v is not None and v > 0)

    # Normalise all taxids to genus and species in one pass
    genus_map: Dict[int, Optional[int]] = {}
    species_map: Dict[int, Optional[int]] = {}
    for tid in all_raw:
        genus_map[tid] = get_rank_taxid(tid, "genus", taxdb)
        species_map[tid] = get_rank_taxid(tid, "species", taxdb)

    def norm(tid, rank_map):
        if tid is None or tid <= 0:
            return None
        return rank_map.get(tid)

    rows = []
    _iter = (
        zip(read_ids, sequences)
        if sequences is not None
        else zip(read_ids, read_lens, gc_arr)
    )
    for i, read_id in enumerate(read_ids):
        if sequences is not None:
            seq = sequences[i]
            rl = len(seq)
            gc = (seq.count("G") + seq.count("C")) / max(rl, 1)
        else:
            rl = int(read_lens[i])
            gc = float(gc_arr[i])

        m_raw = labels_metaphlan.get(read_id)
        k_raw = labels_krakenuniq.get(read_id)
        c_raw = labels_centrifuger.get(read_id)

        m_g = norm(m_raw, genus_map)
        k_g = norm(k_raw, genus_map)
        c_g = norm(c_raw, genus_map)
        m_s = norm(m_raw, species_map)
        k_s = norm(k_raw, species_map)
        c_s = norm(c_raw, species_map)

        rows.append({
            "read_id":             read_id,
            "read_len":            rl,
            "gc":                  gc,
            "metaphlan_raw":       m_raw,
            "krakenuniq_raw":      k_raw,
            "centrifuger_raw":     c_raw,
            "metaphlan_genus":     m_g,
            "krakenuniq_genus":    k_g,
            "centrifuger_genus":   c_g,
            "metaphlan_species":   m_s,
            "krakenuniq_species":  k_s,
            "centrifuger_species": c_s,
            "agreement_genus":     compute_agreement_row(m_g, k_g, c_g),
            "agreement_species":   compute_agreement_row(m_s, k_s, c_s),
        })

    df = pd.DataFrame(rows).set_index("read_id")

    # Add human-readable taxon names
    all_norm_taxids: set = set()
    for col in ("metaphlan_genus", "krakenuniq_genus", "centrifuger_genus",
                "metaphlan_species", "krakenuniq_species", "centrifuger_species"):
        all_norm_taxids.update(
            int(v) for v in df[col].dropna() if pd.notna(v) and int(v) > 0
        )
    name_cache = build_taxid_name_cache(all_norm_taxids, taxdb)

    for tool in TOOLS:
        for level in ("genus", "species"):
            col = f"{tool}_{level}"
            df[f"{tool}_{level}_name"] = df[col].apply(
                lambda tid: taxid_to_name(tid, name_cache)
            )

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 10.  PLOT HELPERS
# ─────────────────────────────────────────────────────────────────────────────

import matplotlib.patches as mpatches


def make_legend_handles(categories=None):
    if categories is None:
        categories = CATEGORIES
    return [
        mpatches.Patch(color=COLORS[c], label=LABEL_NAMES[c])
        for c in categories
    ]


def scatter_ordination(
    ax, coords: np.ndarray, agreement_series: pd.Series,
    title: str = "", xlabel: str = "", ylabel: str = "",
    alpha: float = 0.4, s: float = 4.0, rasterized: bool = True,
):
    """
    Plot a 2D ordination scatter on ax, coloured by agreement pattern.
    Categories are drawn in CATEGORIES order so 'MKC' (most informative)
    is rendered on top.
    """
    for cat in CATEGORIES:
        mask = (agreement_series == cat).values
        if mask.sum() == 0:
            continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=COLORS[cat], s=s, alpha=alpha,
            linewidths=0, rasterized=rasterized,
            label=LABEL_NAMES[cat],
            zorder=CATEGORIES.index(cat),
        )
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
