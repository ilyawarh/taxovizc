# Computing Lab

Текущая база: commit 9528b62e5da4624384c9977561a6e8caa484702a. Python-ядро, CLI, README, environment.yml и все исходные ячейки notebook прочитаны при аудите 2026-09-26.

- [Аудит](AUDIT_20260926.md): приоритеты ошибок и проверенные наблюдения.
- [Машинные результаты](evidence/audit_20260926.json): версии окружения, хеши кода, синтетические примеры.
- [Воспроизведение](../../scripts/audit/reproduce_20260926.py): лёгкие проверки без реальной базы и полного embedding.
- [Бенчмарки](../../benchmarks/README.md): план научной проверки.
- [План](../workflow/ROADMAP.md): сначала P0, затем baseline коррекции.

Рекомендуемая будущая структура пакета: io/adapters, taxonomy, features, representations, neighbors, refinement, evaluation, visualization. Переносить существующий код поэтапно после фиксации корректности. Добавить pyproject.toml и отдельные CPU/GPU зависимости; CLI разделить на validate, features, fit, transform, refine, evaluate, plot. Это проектное решение, перечисленных подкоманд пока нет.

Нормализованная запись свидетельства должна хранить sample_id, read_id, mate/window, tool, database_version, taxonomy_namespace, raw_taxid, candidate_taxids, raw_score, score_type, support_span и parse_status. Итоговая запись: original_taxid, final_taxid, final_rank, action, reason, confidence, confidence_type, neighborhood_support, uncertainty, model_id.

Принцип адаптеров: tool output бывает read-level, marker-level, contig-level и profile-level. Их нельзя приводить к одинаковым голосам без указания семантики. Sourmash gather/profile не превращается автоматически в таблицу назначений всех ридов; MetaPhlAn SAM содержит маркерные выравнивания, а не полный результат индивидуальной классификации каждого рида.
