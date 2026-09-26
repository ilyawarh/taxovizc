# TAXOVIZ_HANDOFF — current agent entry point

Updated: 2026-09-26. Document version: 1.0. Baseline implementation SHA: 9528b62e5da4624384c9977561a6e8caa484702a. Workspace branch: docs/research-workspace-2026-09-26. This document does not imply that the branch has been merged.

## Read first

1. AGENTS.md and START_HERE.md.
2. docs/workflow/ROADMAP.md and DECISIONS.md.
3. docs/computing/AUDIT_20260926.md and evidence/audit_20260926.json.
4. docs/research/RESEARCH_BRIEF.md, METHODS_REVIEW.md and benchmarks/README.md.
5. Both Python source files before implementation changes.

## Intent

Develop assembly-free read-level taxonomic annotation refinement using sequence-derived neighborhoods plus noisy classifier evidence. Initial target: ONT long reads. Single-tool refinement must be possible; multiple tools are optional evidence. Standalone supervised classification and short-read support are later tracks. An unsupervised unnamed bin is not a species identification. Final outputs must preserve original calls, read identity, final rank, uncertainty and reasons for each action.

## Actual state

Existing implementation builds canonical k4/k6 CLR + GC features, PCA, t-SNE/UMAP and agreement-colored plots. It does not implement correction, cluster assignment or ground-truth F1 evaluation. The current documentation branch adds organization, research and reproducible audit probes without changing production behavior.

Verified by synthetic execution: count_reads overcounts @-prefixed quality lines and streaming preallocates phantom rows; GC normalization differs across chunks; MetaPhlAn misses MV-style references and a subsequent unrecognized record can erase a call; conflict categories collapse into tool-only labels; nondefault k=3 raises KeyError; Centrifuger last hit wins. Taxopy 0.14.0 correctly maps 562→561 genus on a mini fixture. Numba/Python counts and reverse complements agree on the sampled probes. See JSON for exact versions and hashes.

Stored notebook output (not rerun): 9,993 reads, 3,213 raw MetaPhlAn calls, 75 PCs retain 54.5% variance. A separate stale cell shows 5,730,402 reads. Do not infer complete provenance or a 100k validation. No actual input FASTQ/SAM/MetaPhlAn DB was supplied in this session.

## Immediate next task

Implement TV-001 (exact read accounting) and TV-002 (fit/transform preprocessing with global GC) in a focused fix branch based on the accepted workspace. Add regression tests for quality @, chunk size/order, duplicate IDs and matrix/metadata cardinality. Do not change the scientific model in that fix. Then implement TV-003/004: versioned marker mapping and explicit availability/conflict semantics. Obtain a small permitted real marker-format fixture when needed; never invent the actual database schema.

## Scientific direction

Benchmark first. Sequence features + PCA + sparse neighbor correction is the baseline. Compare VAE/contrastive learning only after input correctness and independent scoring. Never use agreement labels as truth. Separate read classification, clustering and abundance metrics. A 2D island is not necessarily a taxon; retain unknowns/coarser ranks and protect rare taxa. Closest analogs: Taxometer, LRBinner, MetaBCC-LR, TaxVAMB; TaxDistill is a 2026 preprint. Sources and concrete transfer limits are in METHODS_REVIEW.md.

## Known historical corrections

Original Biomni document is archived unchanged at docs/workflow/archive/TAXOVIZ_HANDOFF_BIOMNI_20260926.md. Its implementation signatures, universal canonical-count formula, Manhattan/Aitchison explanation, claimed input size and MetaPhlAn support inference are not authoritative. Preserve canonical-index correctness, no full-sequence reload in streaming and rank-specific reports. The old forced MetaPhlAn tie-break is legacy behavior, not a correctness invariant for the future classifier.

## Working arrangement and persistence

Research Lab: docs/research. Computing Lab: docs/computing, code, scripts/audit, benchmarks. Workflow Hub: docs/workflow plus this handoff. User selected one GitHub source of truth. Laptop holds code/docs/small fixtures; server is only for authorized compute jobs. No persistent server workspace or interface chats were created. Optional educational VAE tutorial is outside the repository.

Record work in LOG, update ROADMAP status, add decisions with evidence and update this handoff when next task or verified state changes. Include commit/config/input/db/taxonomy hashes and resource metrics for every scientific run. Full production data and large generated matrices are not Git artifacts.
