# Журнал решений

## ADR-001 — единый GitHub workspace

Дата: 2026-09-26. Статус: accepted by owner.
Один taxovizc содержит Research Lab, Computing Lab и Workflow Hub. Основная рабочая копия на ноутбуке; сервер для конкретных разрешённых вычислений. Причина: код и решения переносятся вместе между агентами. Последствие: не хранить здесь full data и restricted inputs. Отдельный учебный VAE-разбор находится вне репозитория по уточнению владельца.

## ADR-002 — текущий статус прототипа

Дата: 2026-09-26. Статус: adopted working assessment, supported by code audit.
Текущий инструмент — визуализатор agreement, не валидированный корректор. Исторический handoff архивируется; точные прошлые числа не перезаписываются, но их интерпретация и уровень доказательности уточняются текущим аудитом. 15M performance и improved F1 остаются неподтверждёнными.

## ADR-003 — inference отделён от 2D display

Дата: 2026-09-26. Статус: proposed architecture, default research direction.
Начать с PCA + neighborhood correction; 2D map используется для диагностики. VAE/contrastive/UMAP representations конкурируют в benchmark. Причина: сохранение соседства и biological assignment нельзя оценить привлекательностью рисунка. Пересмотреть при независимом доказательстве лучшего 2D варианта.

## ADR-004 — исправление не равно принудительному назначению

Дата: 2026-09-26. Статус: proposed model contract.
Хранить исходные labels и reason каждого изменения; поддерживать abstention, coarser rank и unnamed clusters. Species/genome/strain не отождествлять. Причина: unknown taxa, horizontal transfer и ошибки anchors. Проверка: open-set и rare-taxon benchmark.

## ADR-005 — порядок разработки

Дата: 2026-09-26. Статус: adopted work plan, numeric acceptance thresholds pending.
P0 correctness → P1 benchmark → P2 simple baseline → P3 learned representations → P4 scaling → P5 generalization. Текущая ветка содержит документы и audit probes; fixes идут следующим отдельным PR. Причина: не смешивать организацию проекта с невидимым изменением уже апробированной реализации.


## ADR-006 - interpret the Zymo early-time run as exploratory evidence

Date: 2026-10-05.
Status: accepted scope based on owner clarification; biological follow-ups proposed, not implemented.
Context: requested review of two saved TaxoViz maps, full-run LRBinner and third-party qPCR workbook.
Options and evidence: exact-ID join permits comparison of full-run bins on identical plotted reads. Workbook reference percentages match manufacturer 16S characterization, but mixed-sample qPCR formulas/labels are inconsistent. Version changes alter local display neighborhoods without independent truth.
Decision: retain the owner's custom ontime-based earliest-sequence 3% sampling design. Do not call it an error or assume it is representative random sampling. Treat product identity as unconfirmed, manufacturer profiles as method-dependent reference evidence, and K=C/bin agreement as diagnostics rather than ground truth. Preserve original annotations and ambiguous/mobile cases. Store large originals and per-read outputs under ignored runs; compact evidence and report belong to Computing/Research tracks.
Trade-offs: this yields interpretable saved-run evidence but cannot establish correction accuracy, causal parameter effects, absolute mixture abundance, or v2 superiority.
Validation required: frozen per-version config/transforms; correct parser and agreement semantics; independent sequence-level checks for Peptostreptococcus and E. coli/lambda candidates; length/GC-controlled experiments.
Conditions for reversal: independent truth, exact sample preparation/assay documentation or validated reference mapping resolving the ambiguities. Do not request missing preparation details again without new evidence: owner stated the supplied workbook is all that is available.
Related experiment/commit: [EXP-20261005-ZYMO](../computing/experiments/EXP-20261005-ZYMO.md), workspace baseline `2e29b4b9f6b15631d2c588fbff708148efc1de6f`.
Supersedes: no previous architecture decision; sampling clarification overrides the initial provisional concern raised during exploration.

## ADR-007 - publish aggregate Zymo evidence with reproducible figures

Date: 2026-10-06. Status: accepted delivery decision under owner's report/Git request.
Context: the completed analysis needed a readable PDF and self-contained context for other agents.
Options and evidence: PDF provides annotated figures and interpretation; JSON retains complete compact tables and source narratives without per-read identifiers.
Decision: publish the PDF, four figures, aggregate JSON, evidence and scripts on the focused research branch. Keep raw HTML, notebook snapshots, sequences, identifiers and copied upstream source in ignored runs. Synchronize the branch without merging main.
Trade-offs: publication is portable; recomputing figures still requires local ignored inputs. PDF is descriptive, not validation of TaxoViz correction.
Validation: 14 PDF pages rendered and inspected; JSON parses; all 23 bins and 39 qPCR targets retained; checksums recorded in publication manifest. Full run/database provenance gaps remain explicit.
Conditions for reversal: new biological evidence or corrected source data require a versioned analysis update, not silent reinterpretation.
Related experiment/commit: EXP-20261005-ZYMO; baseline 2e29b4b9f6b15631d2c588fbff708148efc1de6f; delivery commit subject `Publish Zymo experiment analysis and full-context report`.
Supersedes: no scientific/architecture decision; supplements ADR-006.
