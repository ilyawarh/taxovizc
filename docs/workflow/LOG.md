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
