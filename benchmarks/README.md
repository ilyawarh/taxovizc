# Benchmark protocol — proposed, not executed

Дата: 2026-09-26. До запуска закрепить manifest, development/test границу и численные критерии приёмки. Маленькие probes в scripts/audit не измеряют биологическую точность.

## Данные и разбиения

1. Synthetic ONT/HiFi communities с точной read→source genome истиной; несколько независимых communities и seeds; разные error/length distributions, abundance skew и близкие штаммы.
2. Reference-held-out сценарии: известные genomes; новый strain; исключённый species; исключённый genus. Фиксировать состав исходных classifier DB, версии и отсутствие запрещённых тестовых референсов. Ground truth и lookup для scoring не доступны алгоритму.
3. Public mock communities (например Zymo) для внешней проверки, с отдельным учётом неопределённых shared sequences. Реальные клинические данные без истины — проверка согласованности и targeted validation, не измерение F1.
4. Отдельные strata: rare taxa, human/host, fungi/archaea, viruses, plasmids, repeats, low-complexity, chimeras, low-quality. Нельзя объявлять unclassified «host» по одному GC.
5. Illumina paired-end — отдельный benchmark. Обе mates, окна одного рида, фрагменты одного генома/strain не должны утекать через supervised split. Финальный перенос проверять на независимых communities/genomes.

Transductive режим может использовать все последовательности и исходные noisy labels текущего test sample, если это явно часть задачи, но не test truth. Hyperparameters/thresholds выбираются на development communities. Inductive fit→transform режим оценивать отдельно. Случайный split ридов одного генома недостаточен для доказательства переносимости на новые виды.

## Сравнения

Один и тот же input и taxonomy reconciliation: original tool; threshold-only; majority/LCA/weighted ensemble; PCA-neighbor correction; альтернативные representations. Taxometer/TaxVAMB на сборке — отдельный assembly-assisted comparator с учётом стоимости сборки, не нечестная замена read-level истины. LRBinner/MetaBCC-LR — clustering baselines; переименование bins через truth допустимо лишь как явно помеченный oracle clustering upper bound, не как classifier score.

## Метрики

- Для каждого rank: per-taxon precision/recall/F1, macro/micro summaries и числа TP/FP/FN. True known read с abstention учитывается как FN его класса; ошибочный вызов даёт FP предсказанному и FN истинному классу. Зафиксировать множество taxa для macro averaging, не исключать таксоны с нулём предсказаний.
- Coverage = доля ридов с назначением на целевом rank. Precision среди назначенных всегда рядом с coverage. Ошибки на unknown taxa и доля корректных отказов — отдельные endpoints; coarse rank не засчитывать species TP.
- Матрица переходов: correct→correct/wrong/abstain; wrong→correct/wrong/abstain; abstain→correct/wrong/abstain. Отдельно correction precision и rare-taxon recall.
- Clustering: ARI/AMI, homogeneity/completeness или B-cubed, purity вместе с completeness и fraction unbinned. Singleton bins не должны обеспечивать иллюзию успеха.
- Profiles: taxon detection F1, L1/Bray–Curtis abundance error, ложная abundance вне community. Разделять read fraction, nucleotide fraction и organism abundance: длина генома и bias секвенирования мешают их отождествлению.
- Calibration: Brier/reliability и risk–coverage curves; отдельно raw confidence и calibrated probability. Variance VAE сама по себе не калибровка taxid.
- Cost: wall time, CPU-hours/GPU-hours, peak RSS/VRAM, temporary disk, output size. Несколько seeds; paired bootstrap по независимым samples/communities или source genomes согласно endpoint, не миллионы зависимых ридов как независимые реплики.

## Абляции и отрицательные контроли

Отключить geometry; отключить labels; отключить GC; k3/k4/k6; разная нормировка/длина; shuffled labels; shuffled neighborhoods; bias всех tools к одному ложному виду; редкий близкий вид рядом с abundant видом; удалённый из DB таксон; отсутствие MetaPhlAn marker. Разные chunk sizes и порядок FASTQ не должны менять preprocessing semantics.

## Gates

P0: точная идентичность и корректный parser. P1: frozen endpoint и независимый test. P2: статистически и практически обоснованный выигрыш сверх baseline без неприемлемой потери rare/unknown; конкретные допуски записать до test. P4: отдельная лестница масштаба; оценки 15M из handoff не являются доказательством.
