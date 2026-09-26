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
