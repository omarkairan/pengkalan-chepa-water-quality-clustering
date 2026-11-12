# Pengkalan Chepa water-quality clustering

This package accompanies the JMCS manuscript by Omar Kairan, Fukayhah@Fatiha Suhami and Nurul Hafizah Azizan. The analysis was executed on the supplied DOE workbook, not simulated data.

## Quick start in Jupyter (Bahasa Melayu)

1. Extract fail ZIP ke satu folder.
2. Letakkan fail asal `LAMPIRAN A - KELANTAN.xlsx` dalam folder yang sama dengan notebook dan skrip Python. Fail asal tidak disertakan semula dalam pakej ini.
3. Gunakan Python 3.12. Pasang pakej dalam `requirements.txt`. Dalam Jupyter, gunakan `%pip install -r requirements.txt`, kemudian restart kernel sebelum import. Pakej rasmi `sdoclust==0.5.1` diperlukan untuk menghasilkan semula semakan parameter major revision.
4. Buka `Pengkalan_Chepa_Analysis.ipynb`. Semak laluan `EXCEL` pada sel konfigurasi.
5. Run sel mengikut susunan. Sel analisis penuh menjalankan 200 subsampel tarikh dan 100 percubaan nilai bawah had pelaporan. Tempoh pelaksanaan bergantung pada komputer.
6. Keputusan baharu disimpan dalam `results_rerun`. Keputusan asal untuk manuskrip sudah tersedia dalam `results`.

Jangan salin bahagian `argparse` ke dalam sel Jupyter atau gunakan `%run` tanpa laluan fail Excel. Notebook menggunakan `run_analysis(EXCEL, OUT)` secara terus untuk mengelakkan masalah argumen kernel `-f`.

## Command line

```bash
python -m pip install -r requirements.txt
python pengkalan_chepa_clustering.py "LAMPIRAN A - KELANTAN.xlsx" --out results_rerun
```

For a quick installation check only, add `--repeats 5 --draws 3`. Those reduced settings do not reproduce the manuscript stability summaries. Use the defaults (200 and 100) for manuscript reproduction.

## Analysis steps

| Step | Implementation and output |
|---|---|
| Audit | Read five annual sheets; retain source strings, censoring flags, reporting limits, Excel row and record ID; reject unknown non-numeric strings and duplicate station–date records. |
| Missingness | One missing NH3-N value at station 4DPCH004, 23 August 2020; 179 complete records in the primary analysis. |
| Censoring | Half-limit substitution for the primary Euclidean analysis; audit retains the fact that these are censored observations. This is not censored-likelihood modelling. |
| Transformation | Natural log of one plus BOD, COD, SS and NH3-N numerical concentrations (mg/L); DO and pH unchanged; all six columns standardised. |
| Descriptive analysis | Overall and station-specific statistics, Spearman rank correlations. No independence-based significance tests. |
| PCA | All six transformed variables; eigenvectors, scores and explained variance exported. PCA does not reduce the clustering inputs. |
| Clustering | SDOclust default settings when `sdoclust` is installed; otherwise a documented one-group fallback matching the manuscript interpretation. K-means and Ward are evaluated for k = 2 through 6. K-means uses 50 starts. |
| Number of groups | Gap statistic for k = 1 through 6 using 200 PCA-aligned uniform reference datasets and the one-standard-error rule; 500 repeated Hopkins statistics are exported as a supplementary clustering-tendency diagnostic. |
| Selection | Maximum silhouette within each conventional method conditional on k = 2 through 6. Detailed profiles use the full-coverage, non-degenerate solution with greatest silhouette. DBI is a separate diagnostic. |
| Cluster-specific separation | Mean, median, range and negative silhouette counts for every cluster in the two-group K-means and Ward solutions. |
| SDOclust robustness | Eleven one-factor-at-a-time parameter settings around the defaults, each across 30 seeds (330 runs). |
| Stability | 200 subsamples of 24/30 dates without replacement; retain all available stations per chosen date; refit transformation/scaling and reselect k. ARI compares shared records. |
| Cluster recovery | Maximum Jaccard overlap against a matching refitted group, on shared observations. |
| Sensitivity | Zero and limit/sqrt(2) substitution, 100 uniform censoring draws, DO saturation, no log transformation, median imputation of missing NH3-N, 30 seeds and omission of each station. |
| Interpretation | Cluster medians/IQR and station/year sample proportions. Labels are not regulatory water classes or continuous-time occupancy estimates. |

## Main results and their limits

- Default SDOclust: one group, 179 observations; silhouette and DBI are undefined.
- Gap statistic: k = 1 under the one-standard-error rule. Gap values for k = 1,...,6 are 0.843938, 0.846581, 0.849602, 0.841371, 0.886294 and 0.903051.
- K-means: two groups, silhouette 0.231094, DBI 1.614006; C1 n=83 and C2 n=96.
- K-means cluster-specific mean silhouettes: C1 0.210949 and C2 0.248510.
- Ward: two groups, silhouette 0.221775, DBI 1.646614.
- PCA PC1 + PC2: 66.3786% of variation in the transformed, standardised data.
- Median date-block ARI: K-means 0.644232; Ward 0.344551.
- K-means reselected two groups in 113/200 subsamples. The partition should be treated as exploratory.
- SDOclust returned one group in 187/200 date subsamples and 27/30 full-data seeds. A perfect ARI comparing two one-group solutions is mathematically possible but is not evidence of useful cluster recovery. Such values are excluded from the main stability summary. Raw sensitivity CSV files retain mathematical ARI values for transparency; read them alongside the group counts.
- Across the major-revision SDOclust parameter grid, one group occurred in 300/330 runs. Every setting had a median of one group.
- Without log transformation, K-means selected five groups and differed substantially from the primary result. Do not report the clusters as independent of preprocessing.
- C1 and C2 are ordered by median NH3-N. Algorithm-native labels in the method columns are zero-based and arbitrary; use the `cluster` column for the manuscript display labels.

## Output dictionary

`audit_all_records.csv` retains raw and parsed values for all 180 records. Numeric parsed censored entries are initially blank by design; use the `_censored`, `_limit` and `_raw` columns to distinguish them from genuinely missing values. The processed numeric values used in clustering are in `cluster_assignments.csv`.

`candidate_metrics.csv` contains all candidate k values, not only the selected results. `selected_metrics.csv` contains the primary method comparison. `cluster_assignments.csv` contains record IDs, stations, dates, original-unit working values, method labels, display clusters and PCA scores.

`stability_replicates.csv` contains all 200 replicates for each method. `stability_summary.csv` reports empirical spread (not confidence intervals). `cluster_jaccard.csv` summarises label-specific recovery. For this run, K-means native label 0 is C1 and label 1 is C2.

`gap_statistic.csv`, `per_cluster_silhouette.csv`, `hopkins_replicates.csv`, `hopkins_summary.csv`, `sdoclust_parameter_runs.csv` and `sdoclust_parameter_summary.csv` contain the major-revision analyses. `sensitivity.csv`, `censor_uniform_draws.csv`, `seed_sensitivity.csv` and `station_omission.csv` preserve the earlier checks. Missing-median sensitivity has 180 fitted observations but agreement is assessed on the original 179 complete records.

`figure1_correlation.png`, `figure2_pca.png` and `figure3_composition.png` are 800-dpi manuscript figures. `figureS1_ward.png` is a supplementary truncated dendrogram. Other CSV files provide descriptive summaries, PCA coefficients, station information and cluster compositions.

`run_metadata.json` records software versions, the seed, replicate counts and the original workbook SHA-256 hash. Floating-point rounding and library versions can cause small differences. The tested scikit-learn 1.5.2 pin is retained for reproducibility. The script no longer requires `sdoclust` at import time; install the official implementation only if it is available in your environment.

## Manuscript information requiring author confirmation

The author names and emails follow the user's instructions. The affiliation is currently recorded at institution/branch level as Universiti Teknologi MARA Cawangan Kelantan, Malaysia. Add the correct department, campus and postal address for each author if required. Confirm the exact preferred spelling and spacing of `Fukayhah@Fatiha Suhami`.

JMCS requests author contributions and a conflict-of-interest statement. These cannot be inferred from author order, so the draft explicitly marks them for confirmation rather than inventing declarations. Funding is stated as none, and the requested acknowledgements to UiTM Cawangan Kelantan and DOE are included.

Before submission, confirm DOE's data-use terms and the laboratory/quality-control information if available. Where related manuscripts reuse these monitoring records, disclose that overlap and identify the distinct research questions and outputs. The present manuscript does not reproduce a paired upstream–downstream analysis or a station vulnerability ranking.

## Key source links

- SDOclust recent evaluation: https://doi.org/10.1016/j.is.2025.102562
- Original SDOclust: https://doi.org/10.1007/978-3-031-46994-7_16
- Official implementation: https://github.com/CN-TU/pysdoclust
- Censoring caveats: https://doi.org/10.1016/j.chemosphere.2006.04.051
- PCA: https://doi.org/10.1098/rsta.2015.0202
- Cluster stability: https://doi.org/10.1016/j.csda.2006.11.025
- ARI: https://doi.org/10.1007/BF01908075

Before resubmission, replace the highlighted repository placeholder in the revised manuscript and response letter with a public URL or DOI. DOE data can be requested using Form DA/1/08, subject to DOE approval and conditions.
