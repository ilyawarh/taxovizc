"""Small audit probes; report observations, not an end-to-end accuracy benchmark.

Run from any directory with the project's CPU dependencies installed:
    python scripts/audit/reproduce_20260926.py > observations.json
The baseline bugs are recorded as observations, not encoded as desired behavior.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import tempfile

import h5py
import numpy as np
from sklearn.decomposition import IncrementalPCA

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import taxoviz_core as tc
import taxoviz_full as tf


def main():
    out = {
        "scope": "Synthetic audit probes; no real FASTQ/database or full pipeline run",
        "python": platform.python_version(),
        "versions": {p: importlib.metadata.version(p) for p in
                     ["numpy", "pandas", "scikit-learn", "biopython", "h5py", "numba", "taxopy"]},
        "sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest()
                   for f in ["taxoviz_core.py", "taxoviz_full.py"]},
    }
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        fq = folder / "tiny.fastq"
        seqs = ["AAAAAAAA", "AAAAAAAA", "GGGGGGGG", "GGGGGGGG"]
        fq.write_text("".join(f"@r{i}\n{s}\n+\n{'@' if i == 0 else 'I'}" +
                              "I" * (len(s) - 1) + "\n" for i, s in enumerate(seqs)))
        ids, parsed = tc.read_fastq_all(str(fq))
        counted = tf.count_reads(str(fq))
        out["read_count"] = {"parsed": len(ids), "count_reads": counted}
        full = tc.build_feature_matrix(parsed, n_jobs=1)
        h5 = folder / "features.h5"
        tf.build_feature_matrix_streaming(str(fq), str(h5), (4, 6), 1., 2, 1, counted)
        with h5py.File(h5) as handle:
            streamed = handle["features"][:]
        out["streaming"] = {
            "hdf5_rows": len(streamed), "real_ids": len(ids),
            "phantom_rows_all_zero": bool(np.all(streamed[len(ids):] == 0)),
            "full_gc_column": full[:, -1].tolist(),
            "stream_gc_column": streamed[:len(ids), -1].tolist(),
            "kmer_blocks_equal": bool(np.allclose(full[:, :-1], streamed[:len(ids), :-1])),
        }
        sam = folder / "tiny.sam"
        def line(read, marker, flag=0):
            return f"{read}\t{flag}\t{marker}\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        mapping = {"SGB1": 562, "MV001-c1": 123}
        sam.write_text(line("r1", "marker_SGB1") + line("r2", "MV001-c1"))
        out["marker_formats"] = tc.parse_metaphlan_sam(str(sam), {"r1", "r2"}, mapping)
        sam.write_text(line("r1", "marker_SGB1") + line("r1", "MV001-c1", 256))
        out["secondary_overwrite"] = tc.parse_metaphlan_sam(str(sam), {"r1"}, mapping)
        cent = folder / "cent.tsv"
        cent.write_text("r1\tref\t562\t10\n")
        out["headerless_centrifuger"] = tc.parse_centrifuger(str(cent), {"r1"})
        cent.write_text("readID\tseqID\ttaxID\tscore\nr1\ta\t562\t100\nr1\tb\t1280\t10\n")
        out["centrifuger_last_record"] = tc.parse_centrifuger(str(cent), {"r1"})
        # Minimal synthetic taxonomy, not a downloaded NCBI snapshot.
        nodes = [(1,1,"no rank"),(2,1,"superkingdom"),(561,2,"genus"),(562,561,"species")]
        (folder / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\n" for t,p,r in nodes))
        (folder / "names.dmp").write_text("".join(f"{t}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for t,name in
                                                        [(1,"root"),(2,"Bacteria"),(561,"Escherichia"),(562,"Escherichia coli")]))
        db = tc.load_taxonomy(str(folder))
        out["taxonomy"] = {"562_genus": tc.get_rank_taxid(562, "genus", db),
                           "562_species": tc.get_rank_taxid(562, "species", db)}
    out["agreement"] = {"one_only": tc.compute_agreement_row(1,None,None),
                        "three_conflict": tc.compute_agreement_row(1,2,3),
                        "two_conflict": tc.compute_agreement_row(None,2,3)}
    out["clr_length_effect"] = tc.clr_transform(np.array([[1,0,3],[10,0,30]])).tolist()
    rng = np.random.default_rng(42)
    checks = {}
    for k in (4,6):
        vectorizer = tc._make_numba_vectoriser(k)
        seqs = ["", "A", "NNNNNN", "ACGTACGTNNACGT"] + ["".join(rng.choice(list("ACGTN"),128)) for _ in range(20)]
        checks[str(k)] = all(np.array_equal(vectorizer(s), tc._vectorise_read_python(s,k)) and
                            np.array_equal(vectorizer(s), vectorizer(tc.reverse_complement(s))) for s in seqs)
    out["numba_python_and_strand_agree"] = checks
    for k in (3,5):
        out[f"canonical_vocab_k{k}"] = tc.build_canonical_kmer_index(k)[1]
    try:
        tc.build_feature_matrix(["ACGTACGT"], k_values=(3,), n_jobs=1)
        out["nondefault_k"] = "accepted"
    except Exception as exc:
        out["nondefault_k"] = f"{type(exc).__name__}: {exc}"
    try:
        ipca = IncrementalPCA(n_components=3)
        ipca.partial_fit(rng.normal(size=(8,5)))
        ipca.partial_fit(rng.normal(size=(2,5)))
        out["ipca_small_final_batch"] = "accepted in this installed version"
    except Exception as exc:
        out["ipca_small_final_batch"] = f"{type(exc).__name__}: {exc}"
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
