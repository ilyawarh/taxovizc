# TaxoViz: план разработки

Версия документа: 1.0, 2026-09-26. Исходный код: 9528b62e5da4624384c9977561a6e8caa484702a. Статус: план принят как рабочее направление; научные архитектуры — кандидаты до бенчмарка. Обновлять статус задач и handoff после каждой итерации.

## Последовательность и контрольные точки

| Этап | Статус | Результат | Условие перехода |
|---|---|---|---|
| P0. Восстановить корректность | TODO, высший приоритет | Воспроизводимые данные и честная диагностика | Закрыты критические ошибки аудита; read accounting и streaming equivalence проверены |
| P1. Зафиксировать оценку | TODO; спецификация подготовлена | Версионированный benchmark и базовые метрики | Truth из независимого источника; split и endpoints заморожены |
| P2. Простой корректор | TODO | PCA + локальная коррекция с отказом | Добавочная ценность сверх исходного tool и ансамбля на development; test не использован для выбора |
| P3. Представления | TODO | Сравнение k-mer/PCA, VAE, contrastive и optional UMAP | Победитель по качеству, сохранению редких таксонов и ресурсам, не по картинке |
| P4. Масштабирование | TODO | Контроль ресурсов, checkpoints, ограниченная визуализация | Измеренные 100k → 1M → полный разрешённый набор; нет потери идентичности |
| P5. Обобщение | TODO | Новые tools, short reads, внешние embeddings | Отдельная валидация каждого режима |

P0 и подготовка P1 могут идти параллельно как задачи; тяжёлые вычисления пока не запускались. P4 начинается с проектирования в P0, но полномасштабный benchmark выполняется после корректности. Никакие оценки времени из старого README не считать измеренным SLA.

## P0 — задачи для ближайшего coding-сеанса

- TV-001: точный FASTQ accounting. Парсить записи, проверять malformed/duplicate ids. После записи матрицы требовать равенство n_parsed = n_written = len(ids) = feature_rows = metadata_rows = prediction_rows. Избыточный preallocation не допускается без resize. Проверить quality-строки с @, gzip/bz2, empty и chunk boundaries.
- TV-002: общая модель preprocessing. Вычислять GC statistics глобально либо fit на явно сохранённой обучающей выборке, затем transform всех chunks. Не нормировать GC отдельно по блоку. Проверить инвариантность признаков к перестановке/разбиению; PCA/IPCA сравнивать по подпространству и downstream neighbors, а не требовать побитового равенства.
- TV-003: MetaPhlAn adapter. CLI/config для точной версии marker DB; явное разрешение reference→marker→clade→taxonomy; SGB и MV fixtures; статус unresolved_marker; многократные SAM hits агрегировать по явной политике с flag/MAPQ/alignment span. Формат реальной pkl ещё нужно проверить на небольшом разрешённом примере. Убрать зависимость от абсолютного пути.
- TV-004: разделить availability, conflicts и agreement. Хранить число ответивших tools, исходные labels, множество конфликтов, rank coverage. Полное разногласие не превращать в M_only. Legacy column при необходимости сохранить под явно устаревшим именем; миграцию задокументировать.
- TV-005: остальные parsers/taxonomy. Header-aware Centrifuger, multi-hit aggregation, malformed counters, нормализация read/mate IDs, merged/deleted taxid, namespaces и версии. Не превращать любой сбой в biological unclassified.
- TV-006: CLI validation и CPU окружение. Проверять k, положительность pseudocount, число PCs/reads/features, first IPCA batch, параметры embeddings, output paths. Разделить core/visualization/GPU, убрать обязательный rapids из CPU установки; сохранить точные версии проверенного окружения. Warning suppression сделать локальным.
- TV-007: очистить и перепроверить notebook. Сохранить исторические outputs отдельно перед очисткой; единый CONFIG; убрать двойную загрузку полного FASTQ; исполнить сверху вниз на разрешённом fixture. Показать фактически достигнутую PCA variance и предупреждение о недостигнутом threshold.

Definition of done P0: synthetic regression tests, один согласованный notebook/CLI smoke-run, отчёт обо всех входных/выходных строках; новые тесты не требуют многогигабайтных баз. Ошибки P0 исправлять отдельным PR от текущего организационного PR.

## P1 — baseline и критерий успеха

Спецификация: ../../benchmarks/README.md. Включить исходные классификаторы по отдельности, majority/LCA/weighted ensemble, confidence thresholding без геометрии, затем варианты коррекции. Обязательна абляция «только sequence», «только labels», «labels + geometry». Сравнение Bracken делать отдельно на abundance endpoint.

На development выбрать primary endpoint: read-level macro-F1 на species с учётом abstentions, плюс constraint по rare-taxon recall и wrong-correction rate. Зафиксировать точные допустимые изменения до закрытого test. Если не достигнуты, результат считается отрицательным, а не скрывается выбором другого rank или seed. Числовые пороги пока не назначены: их нужно обосновать базовой вариабельностью и прикладной целью.

## P2 — минимальная модель коррекции

Вход: валидированные признаки X, optional external Z, таблица evidence E, taxonomy T. Начать с 1–2 исходных инструментов; архитектура не должна требовать строго M/K/C.

1. Построить representation независимо от предсказаний для первого baseline; fit transforms фиксировать.
2. Вычислить kNN/ANN и расстояния; проверить recall ANN на малом exact поднаборе. Mutual-kNN и порог расстояния — кандидаты, но они могут изолировать редкие группы.
3. Получить исходные мягкие распределения/веса из E. Некалиброванные scores разных tools не складывать как одинаковые вероятности. Корреляцию ошибок учитывать через ablation/tool dropout, а не предположение независимости голосов.
4. Ограниченно агрегировать поддержку соседей. Baseline: q_i = (1-lambda_i)*p_i + lambda_i*sum_j(w_ij*p_j)/sum_j(w_ij). Это правило сглаживания, не доказанная Bayesian posterior. Для нулевой поддержки сохранить неопределённость. Не делать обязательный global n_reads × n_taxa dense массив: sparse top candidates.
5. Выбирать наиболее конкретный поддержанный узел taxonomy с порогами confidence, margin, neighborhood support и OOD. Если уверенности недостаточно — сохранить вызов, поднять rank или abstain согласно заранее зафиксированной политике. Unknown не должен означать «назначить ближайший известный вид».
6. Сохранять original/final/reason, а также изменения correct→wrong. Точки без labels могут образовать безымянную группу; автоматическое именование запрещено без свидетельства.

Итеративное label propagation, графовые нейросети и whole-cluster reassignment — отдельные следующие эксперименты, поскольку сильнее распространяют систематические ошибки.

## P3 — контролируемая матрица экспериментов

| Компонент | Сначала | Затем | Пока отложить |
|---|---|---|---|
| Sequence features | k3/k4 frequency/Hellinger/CLR, quality/length audit | k4+k6 block weighting, coverage histograms, minimizer evidence | Dense high-k vocabulary |
| Representation | PCA 16/32/64 | VAE 16/32/64; contrastive; denoising AE; UMAP 10–50 | Foundation model training from scratch |
| Grouping | Neighbors without hard bins | HDBSCAN, Leiden, mixture/prototypes | Mandatory 2D island = species rule |
| Labels | Confidence thresholds + hierarchy | Learned calibration, multimodal loss, open-set head | Uncalibrated all-tool agreement as truth |
| Heavy embedding | Frozen DNABERT-S small subset if simpler methods fail | Distillation after demonstrated gain | Full 15M transformer inference as first baseline |

Числа размерностей — проектная сетка, не «оптимальные параметры». Подбирать на development по независимой taxonomic truth. Для external embeddings записывать producer/version/metric/dimension/fit population и mapping read_id; измерять качество независимо от внутреннего encoder.

## P4 — ресурсы и инженерия

- Rolling 2-bit canonical k-mer counter и пакетный compiled loop; benchmark threads/processes/nogil, не полагаться на n_jobs=-1. Сравнить с корректным Python reference.
- Уменьшить HDF5 chunk относительно нынешних 443.4 MB uncompressed при 50k×2217×4. Chunk byte target выбирать измерением I/O.
- Features/latents/metadata — потоковые columnar arrays/Parquet или memmap, taxonomy/name dictionaries отдельно; отказаться от списка миллионов Python row-dicts и повторных hover strings.
- Сохранять encoder/PCA/scaler и manifest; atomic stages, resume и checksums. Удалять промежуточный файл только по явной cleanup policy после успешного завершения всех потребителей.
- UMAP/GPU fallback должен быть видимым и управляемым; выбирать CPU/GPU явно для сравнения. Наличие восьми GPU не означает multi-GPU выполнение текущего cuML вызова.
- Визуализация: стратифицированная выборка и density raster; полный read-level результат отдельно. Для редких групп показывать oversampling weights, не выдавать плотность визуальной выборки за abundance.
- Размеры до запуска: 15M×2217 float32 ≈133.02 GB; 15M×75 ≈4.5 GB; 15M×32 ≈1.92 GB. Граф k=30 содержит 450M направленных рёбер: около 3.6 GB только int32 indices + float32 weights, без offsets, symmetrization и ANN overhead.

## P5 — обобщение и продуктовый интерфейс

Проверить Kraken family, Centrifuge/Centrifuger, Kaiju и MetaPhlAn по отдельным адаптерам. Для Sourmash явно ограничить поддержку доступным типом output. Затем Illumina paired-end, HiFi, несколько samples и внешние embeddings. Cross-sample coverage для сырых ридов требует собственного способа установления соответствий; нельзя копировать contig abundance matrix без её построения.

Дальнейший standalone classifier требует reference/label training и отдельного benchmark held-out species/genus. Чистый unsupervised binning может оставаться полезным режимом, но не называть его полностью таксономическим классификатором.


## 2026-10-05 - completed exploratory Zymo review (TV-R001)

Status: DONE for artifact retrieval, saved-output interpretation and exact-ID comparison; does not complete P0/P1 or validate correction.
Evidence: [research report](../research/ZYMO_FECAL_REVIEW_20261005.md), [experiment](../computing/experiments/EXP-20261005-ZYMO.md).

- Intentional earliest-time 3% subset: 290,227 reads, matched to all 23 full-run bins.
- Real-run evidence reinforces TV-004 (hidden conflicts in -only- colors) and TV-007 (stale notebook outputs, missing independent v1 fitted state). TV-001/002/003/004/007 remain open; this review changed no production behavior.
- P1 preparation advanced: public Zymo DNA/16S references located and qPCR limitations audited. No independent read-level truth or frozen evaluation split yet.
- Proposed research follow-up: verify bin 0 and mobile/control candidates in bins 7/11/13 by independent alignment; assess length/GC effects on a fixed read set after P0. Maintain early-time sampling as a user-selected design.

## Zymo delivery - 2026-10-06

TV-R001 publication: DONE. [PDF](../research/reports/Zymo_Experiment_20261006.pdf) and [standalone agent context](../research/ZYMO_AGENT_CONTEXT_20261006.json) complete the exploratory review. No change to P0 acceptance: parser/agreement fixes and independent validation are still required. Follow-ups: bin 0 references, bins 7/11/13 mobile/control candidates, frozen transforms and length/GC-controlled analysis.
