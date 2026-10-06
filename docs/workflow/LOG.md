# Лог разработки

## 2026-09-26 — workspace bootstrap и ревизия baseline

- Владелец утвердил вариант 1: единый GitHub repository с локальной разработкой и server compute jobs. Учебный VAE-разбор исключён из repository deliverables.
- Зафиксирован baseline main: 9528b62e5da4624384c9977561a6e8caa484702a.
- Прочитаны Python core/CLI, README, environment и notebook; проведён read-only обзор первичной литературы до 2026-09-26.
- Выявлены ошибки подсчёта FASTQ/HDF5, chunk-wise GC scaling, MetaPhlAn marker/multi-hit handling, agreement semantics и ограничения adapters. Synthetic probes сохранены с версиями и hashes.
- Сохранённый notebook подтверждает 9,993 rows и achieved PCA variance 54.5%; handoff claims про 100k и MetaPhlAn support требуют исправления интерпретации.
- Созданы три трека, ROADMAP, benchmark protocol, decision log, templates, START_HERE, AGENTS и актуальный handoff. Исходный handoff архивирован без изменения. Навигация добавлена в README; создан .gitignore.
- Production Python-код, environment.yml и notebook не изменены. Нет full pipeline run, обучения VAE, measured F1 gain или 15M benchmark.
- Изменения подготовлены в docs/research-workspace-2026-09-26 для draft PR. Слияние и release не выполнялись.

Следующая задача: TV-001 + TV-002 (accounting и deterministic preprocessing), затем TV-003/TV-004 и чистый fixture smoke-run.


## 2026-10-05 - Zymo fecal saved-run review

- Completed owner-requested retrieval and descriptive comparison on branch `research/zymo-fecal-review-2026-10-05`, based on workspace HEAD `2e29b4b9f6b15631d2c588fbff708148efc1de6f`. Existing line-ending-only changes were preserved.
- Downloaded both Plotly HTMLs, four SVGs, updated server notebook (after owner saved session), singular `zymo_lrb/binning_result.pkl`, log and supporting code/order files into ignored `runs/zymo_fecal_20261005/`. All 12 remote checksums verified. Preserved stale initial notebook separately; copied qPCR workbook without editing it.
- Verified 290,227 identical plotted read IDs across versions and matched every ID to full-run LRBinner. Audited 9,674,223 FASTQ records/lengths and the complete pickle-to-bins.txt mapping: 23 bins, no missing or overlapping memberships.
- Owner clarified that the 3% subset deliberately selects earliest sequences via a custom ontime-based sampler. This is the intended experimental design, not a sampling bug. File-order concentration is retained as descriptive provenance; no random-sampling claim.
- Quantified bin/taxon/length/GC associations and display-neighborhood changes. Strong added-culture candidate groups coexist with mixed bins and disagreement; v2 improvement is not demonstrated. Flagged Peptostreptococcus bin 0, E. coli/lambda-labelled bins 7/11/13, and legacy agreement semantics.
- Retrieved public manufacturer DNA/16S references through protocol Appendix A. Workbook's 354 reference percentages match the public 16S values. Product/lot identity remains unknown. qPCR formula and target-label inconsistencies prevent exact mixture truth.
- Saved [research report](../research/ZYMO_FECAL_REVIEW_20261005.md), [experiment record](../computing/experiments/EXP-20261005-ZYMO.md), compact evidence and reproducible scripts. Source data, read IDs and large plots remain ignored. No production fixes, classifier reruns, model fitting, commit, merge or push.

Next: preserve per-version provenance and address P0 correctness; use independent alignment for bin 0 and bins 7/11/13 before any correction. See ADR-006.

## 2026-10-06 - Zymo report publication

Completed the saved-run review delivery: 14-page PDF, standalone aggregate JSON context, four original-coordinate comparison figures, supplemental length/rare-target evidence, and reproducible publication scripts. PDF rendered with Poppler and visually inspected; image scaling corrected. All 23 bins and all 39 qPCR target records preserved in the handoff. No model rerun or production-code change. Owner's intentional earliest-3% sampling design retained.

Workspace audit: Research/Computing/Workflow artifacts linked; large originals and read IDs remain under ignored runs; database/seed/provenance gaps explicitly documented rather than invented. Existing unrelated line-ending changes preserved. User authorized commit/push of research/zymo-fecal-review-2026-10-05; main remains unchanged. Delivery commit is discoverable by subject `Publish Zymo experiment analysis and full-context report`.

Artifacts: [PDF](../research/reports/Zymo_Experiment_20261006.pdf), [agent context](../research/ZYMO_AGENT_CONTEXT_20261006.json), [publication manifest](../computing/evidence/zymo_20261005/publication_manifest.json). See ADR-007. P0 correctness and independent biological validation remain open.

Synchronization receipt: research delivery [659dffdebef981c62d832f2f78ef6dd9d68cfe67](https://github.com/ilyawarh/taxovizc/commit/659dffdebef981c62d832f2f78ef6dd9d68cfe67) pushed successfully. HTTPS had no configured write credentials; existing GitHub SSH identity authenticated as the repository owner. Origin fetch URL remains HTTPS; push URL now uses SSH to the same repository. No credentials were copied. Main was not merged.
