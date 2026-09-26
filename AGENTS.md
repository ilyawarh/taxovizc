# Agent instructions for TaxoViz

Read START_HERE.md, TAXOVIZ_HANDOFF.md, docs/workflow/ROADMAP.md, docs/workflow/DECISIONS.md and the relevant track documents before edits. Read both taxoviz_core.py and taxoviz_full.py before changing their contracts. No external agent service or server workspace is configured by these files.

## Authority and evidence

- Current user instructions have priority. Current code determines existing behavior; accepted decisions determine intended behavior. Archive documents are historical evidence, not immutable algorithmic requirements.
- Distinguish observed-in-code, reproduced-on-synthetic-input, stored-notebook-output, historical-report, hypothesis and proposed design. Never relabel a proposal as implemented.
- The original Biomni handoff is preserved in docs/workflow/archive. Its tie-break, geometric claims, sample size and parser contracts have corrections in the audit. Deterministic behavior does not justify incorrect semantics.
- Agreement among tools is not ground truth. Do not use the evaluated labels to define the truth, select thresholds or claim scientific improvement.
- Do not silently equate species, strain, genome, bin and individual organism. Preserve unresolved, unknown and mobile-element cases.

## Change discipline

- Work on a focused branch; separate correctness fixes, algorithm experiments and presentation changes. Keep main stable through review. Do not merge or release merely because documents describe a milestone.
- Preserve read identity and order. No read may disappear silently. Report duplicate identifiers and missing/malformed parser records explicitly.
- Separate representation used for inference from 2D display. Preserve original calls and expose every correction with a reason and calibrated score or explicitly uncalibrated score.
- Use small synthetic fixtures for parser and streaming bugs. Test invariants and scientific failure modes, not visual attractiveness. Run heavy benchmarks only as specifically authorized compute jobs.
- Every scientific run records commit, configuration, seeds, input/database checksums, taxonomy versions, environment, read counts, exclusions and resources.
- Keep FASTQ, databases, large matrices, credentials and restricted sample identifiers outside Git. Use synthetic/public fixtures. Check third-party license compatibility before copying code.
- Local laptop: source, documentation and small tests. Server: authorized job-specific checkout or source snapshot only, under the applicable institutional policy; no persistent project-management service.

## Context persistence

At completion of a meaningful task, update docs/workflow/LOG.md, the corresponding ROADMAP status, and TAXOVIZ_HANDOFF.md when state or next steps change. Record rationale in DECISIONS.md, not only chat. Use templates for experiments and decisions. Link exact commits and evidence. Keep negative results and rejected hypotheses.

Do not spawn subagents by default. Delegation requires an explicit task instruction. The optional VAE tutorial is outside this repository; method evaluation involving VAE belongs in Research Lab as normal research.
