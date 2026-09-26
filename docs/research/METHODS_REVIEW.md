# Research Lab: обзор методов и конкурирующих решений

Поиск: 2026-09-26. Это целевой обзор ближайших алгоритмических семейств, не исчерпывающий systematic review и не патентный поиск. Основания — первичные статьи и официальные репозитории/документация. Числа из чужих бенчмарков не переносим на TaxoViz.

## Ближайшие аналоги

| Метод | Объект и задача | Что применимо в TaxoViz | Ограничение переноса |
|---|---|---|---|
| MetaBCC-LR, 2020 [1] | Биннинг сырых длинных ридов по композиции и k-mer coverage | Независимый от сборки контекст abundance; baseline read-level clustering | Кластеры не дают taxid; subsampling и неравномерная представленность |
| LRBinner, 2021/2022 [2] | VAE объединяет композицию и гистограммы встречаемости k-mer; кластеризуются длинные риды | Наиболее прямой read-level VAE baseline; выделение групп разного размера | Не готовый корректор профайлеров; исторические ONT ошибки отличаются от современных |
| Taxometer, 2024 [3] | Коррекция таксономии контигов нейросетью по TNF и abundance, иерархический loss | Ближайший концептуальный конкурент; обучение на частично размеченных фрагментах, отказ по порогу | Контиги и coverage надёжнее сырых коротких ридов; labels остаются шумными |
| TaxVAMB, 2026 [4] | Полуобучаемый мультимодальный VAE: композиция/abundance + таксономия контигов | Пример объединения геометрии и таксономических свидетельств | MAG recovery не равен read-label F1; использовать как сравнительную архитектуру |
| TaxDistill, preprint 2026 [5] | Дистилляция GenomeOcean в корректор таксономии контигов | Снижение переобучения на шумных метках, мягкие targets, отказ | Препринт; дорогой teacher, ground-truth-free обучение не делает псевдометки истинными |

Taxometer опровергает широкую формулировку «коррекция классификаторов по композиционному сходству впервые». LRBinner опровергает «VAE для кластеризации метагеномных ридов впервые». Научная ниша должна быть конкретнее этих формулировок.

## Методы, из которых стоит заимствовать компоненты

| Метод | Подход | Приоритет |
|---|---|---|
| VAMB, 2021 [6] | VAE объединяет композицию и coabundance контигов; затем отдельная кластеризация | Высокий как архитектурный baseline, не как доказательство превосходства VAE |
| binny, 2022 [7] | Итеративная ординация/кластеризация контигов с оценкой по маркерам | Источник текущей идеи визуализации; качество островов проверять независимо |
| SemiBin / SemiBin2, 2022/2023 [8] | Обучение представлениям через ограничения/самообучение; contig-level binning | Высокий: reverse-complement и фрагменты одного длинного рида как осторожные positive pairs |
| COMEBin, 2024 [9] | Multi-view contrastive learning, признаки последовательности и coverage, Leiden | Высокий: сравнить с VAE; отличать augmentation от размножения независимых наблюдений |
| WEVOTE, 2016 [10] | Ансамбль read-level классификаторов с таксономическим голосованием | Обязательный простой baseline: доказывать добавочную ценность геометрии сверх ансамбля |
| MetaMaps, 2019 [11] | Вероятностное назначение длинных ридов и оценка состава по референсам | Baseline ONT; модель нескольких кандидатов и sample-level контекста |
| Bracken, 2017 [12] | Перерасчёт abundance по Kraken и структуре базы | Сравнивать профильные оценки; не считать прямым аналогом коррекции каждого read_id |
| Metabuli, 2024 [13] | Совместный анализ DNA/AA для классификации | Сильный исходный классификатор для проверки универсальности корректора |
| DNABERT-S, 2025 [14] | Видоориентированное представление ДНК | Второй эшелон: frozen embeddings на небольшом стратифицированном наборе; учитывать предобучение и расход |

Не переносить механически must-link/cannot-link из contig-binning. Разные исходные taxid могут быть ошибкой, а два разных рида — одним организмом. Random negatives в контрастивном обучении могут стать ложными отрицательными парами. Mobile/chimeric reads нарушают даже некоторые positive-pair предположения.

## Краткая траектория области

Композиционные сигнатуры и coverage были основой биннинга задолго до нейросетей. Затем появились direct long-read binning и обучаемые совместные представления. В 2021–2024 годах VAE и contrastive learning стали важными семействами MAG-биннинга, но большей частью работали с контигами. Следующий шаг — использование этих признаков для исправления самих таксономических аннотаций (Taxometer), затем явная интеграция таксономии с VAE (TaxVAMB). Foundation-model embeddings и дистилляция — перспективное, но вычислительно более тяжёлое направление. Ни одна стадия не снимает проблемы неполной базы, шумных labels и редких организмов.

## Выбор алгоритмов для первых экспериментов

1. Композиция: k=3/4 с L1 или Hellinger, отдельно CLR с контролем pseudocount; k=6 как абляция. При 150 bp максимум 145 окон k=6 на 2080 canonical признаков: разреженность делает длину и сглаживание особенно существенными.
2. PCA (например 16/32/64 координаты; это диапазон эксперимента) + Euclidean neighbors — прозрачный и дешёвый baseline. Долю дисперсии логировать, а не считать критерием таксономического качества.
3. Разреженный kNN/mutual-kNN граф + ограниченное локальное распространение вероятностей. Якоря не считать безошибочными; избегать неограниченной рекурсии новых псевдометок.
4. HDBSCAN на PCA/latent coordinates, Leiden на графе, mixture/prototype assignment — альтернативы. HDBSCAN даёт noise и группы разной формы, но чувствителен к неодинаковой плотности; Leiden может делить континуум, а его resolution не является числом видов. Mixture models удобны для мягких вероятностей, но требуют проверки предположений о форме групп.
5. VAE и denoising autoencoder сравнить с contrastive representations; потери реконструкции взвешивать по модальностям. VAE не превращает слабые признаки в отсутствующую информацию о происхождении.
6. UMAP в 10–50 измерениях — отдельный вариант; 2D UMAP/t-SNE — преимущественно диагностика. HDBSCAN по 2D разрешён как сравнительный эксперимент, но не как основа по умолчанию: UMAP может искусственно разрывать группы и искажать плотность [15].

Нецелесообразны для полного 15M набора: полная матрица попарных расстояний, стандартные NMDS/PCoA на ней и обязательное отображение всех ридов в четырёх WebGL-панелях. Нецелесообразно сейчас: foundation model с нуля, автоматическое удаление всех «неподдержанных» островов, безусловное majority relabeling целых кластеров.

## Источники

1. Wickramarachchi et al. MetaBCC-LR, Bioinformatics 2020. https://doi.org/10.1093/bioinformatics/btaa441 ; https://github.com/anuradhawick/MetaBCC-LR
2. Wickramarachchi & Lin. LRBinner, Algorithms for Molecular Biology 2022; WABI 2021. https://doi.org/10.1186/s13015-022-00221-z ; https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2021.11
3. Kutuzova et al. Taxometer, Nature Communications 2024. https://doi.org/10.1038/s41467-024-52771-y
4. Kutuzova et al. TaxVAMB, Nature Biotechnology, 27 April 2026. https://doi.org/10.1038/s41587-026-03098-0 ; https://github.com/RasmussenLab/vamb
5. Ye et al. TaxDistill, arXiv v1, 22 May 2026, preprint. https://arxiv.org/abs/2605.28868 ; full text https://arxiv.org/html/2605.28868v1
6. Nissen et al. VAMB, Nature Biotechnology 2021. https://doi.org/10.1038/s41587-020-00777-4
7. Hickl et al. binny, Briefings in Bioinformatics 2022. https://doi.org/10.1093/bib/bbac431 ; https://github.com/a-h-b/binny
8. SemiBin: https://doi.org/10.1038/s41467-022-29843-y ; SemiBin2: https://academic.oup.com/bioinformatics/article/39/Supplement_1/i21/7210480 ; https://github.com/BigDataBiology/SemiBin
9. Wang et al. COMEBin, Nature Communications 2024. https://doi.org/10.1038/s41467-023-44290-z
10. WEVOTE, 2016. https://pmc.ncbi.nlm.nih.gov/articles/PMC5040256/
11. Dilthey et al. MetaMaps, Nature Communications 2019. https://doi.org/10.1038/s41467-019-10934-2
12. Lu et al. Bracken, PeerJ Computer Science 2017. https://doi.org/10.7717/peerj-cs.104 ; https://github.com/jenniferlu717/Bracken
13. Kim & Steinegger. Metabuli, Nature Methods 2024. https://doi.org/10.1038/s41592-024-02273-y
14. DNABERT-S, Bioinformatics 2025. https://academic.oup.com/bioinformatics/article/41/Supplement_1/i255/8199363 ; https://github.com/MAGICS-LAB/DNABERT_S
15. UMAP authors' documentation: https://umap-learn.readthedocs.io/en/latest/clustering.html
16. Sourmash authors' documentation: https://sourmash.readthedocs.io/en/stable/classifying-signatures.html

Статьи и программные интерфейсы проверялись выборочно по доступным первичным материалам; полный воспроизводимый запуск конкурирующих инструментов в этом сеансе не выполнен. Перед интеграцией кода фиксировать release/commit и лицензию конкретного компонента.
