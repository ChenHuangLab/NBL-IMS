# Neuroblastoma IMS analysis scripts

English-language copies of 14 R scripts from seven selected modules of the neuroblastoma (NBL) analysis project. The scripts cover DNA methylation, Cox regression and NMF, three-group IMS characterization, single-cell processing, cell subgroups and communication, HRC deconvolution, and machine-learning feature selection.

## Repository structure

| Folder | Scripts | Purpose |
| --- | ---: | --- |
| `01_methylation_analysis` | 5 | Package setup, beta-matrix loading and quality control, differential methylation, and CpG identifier recovery |
| `02_cox_regression_and_nmf` | 3 | Historical Cox analysis, univariable Cox screening, and NMF clustering |
| `03_ims_three_group_characterization` | 1 | Three-group differential-expression visualization |
| `04_single_cell_preprocessing_and_clustering` | 1 | Single-cell import, doublet detection, quality control, SCTransform integration, clustering, and marker exploration |
| `05_cell_subgroups_and_communication` | 2 | Risk-state annotation, CellChat analysis, and T-cell/myeloid subclustering |
| `06_hrc_deconvolution_of_ims` | 1 | Reference expression profiles and CIBERSORT deconvolution |
| `07_machine_learning` | 1 | Iterative feature selection using random forest, logistic regression, and XGBoost |

The folder numbers indicate conceptual order. These are research-session scripts, with exploratory alternatives and historical examples retained; they are not a single automated workflow. In particular, the NMF script also retains an APA/breast-cancer example from the supplied source.

## Inputs and environment

Use R with the packages required by the relevant module. Package versions from the original analyses were not supplied and have not been reconstructed. The scripts contain existing Windows and Linux paths; review and adapt all input paths, output paths, and `setwd()` calls before execution.

| Module | Main requirements and inputs |
| --- | --- |
| Methylation | ChAMP, minfi, wateRmelon, impute, limma, GEOquery, Biobase, data.table, plotting packages; phenotype `group.txt`, methylation `data.txt`, and prepared R workspaces |
| Cox/NMF | survival (including `Surv()` and `coxph()`), NMF, Rmpi/doMPI, pheatmap; survival/expression objects and prepared Cox results |
| IMS groups | pacman, tidyverse, ggrepel, reshape2, rio, ggsci; `NMF_3grp.matrix.rds` containing `NMF_3grp` |
| Single-cell processing | Seurat, DropletUtils, scDblFinder, scuttle, sctransform, glmGamPoi, harmony, clustree, irlba, stringr, ggplot2; raw/prepared objects from GSE137804, GSE140819, GSE192906, and the GOSH cohort |
| Subgroups/communication | Seurat, CellChat, harmony, dplyr, patchwork and the other imported single-cell packages; `sce_all6.rds` with the expected metadata |
| HRC deconvolution | Seurat, CIBERSORT, bseqsc, e1071, preprocessCore, purrr, IRanges, GenomicRanges, plotting packages; prepared reference, bulk-expression, and group objects |
| Machine learning | randomForest, pROC, Matrix, xgboost; pre-existing `h2omoedl_save`, `ooutre_sorted`, `model`, `exp`, and `nmfg` objects |

Only the selected scripts are included. Data files, models, manuscript files, analysis results, and scripts from other modules are not included. The methylation scripts refer to `functions_DEM.R`, which was not present in the selected folders and is not bundled. Obtain required data and auxiliary functions separately. Some files use `load()` for files ending in `.rds`; retain the original serialization convention or adapt it to the actual input file format.

The English copies translate three path components: `generated_data`, `catecholamine_related`, and `three_group_markers`. Point these paths to the corresponding local directories when running the scripts.

## Running the scripts

1. Prepare the required input data and R packages for the module of interest.
2. Inspect the script and update local paths. Restore or create its required workspace objects and metadata columns.
3. Run the applicable analysis sections in R or RStudio. Review exploratory sections before running an entire file, particularly those that clear the workspace, install packages, change directories, or write outputs.

Several scripts rely on objects created in earlier interactive sessions. Examples retained from the source include `ne2` in the deconvolution script and `meta$celltype` for CellChat group sizes while communication groups use `newthree`. Reconcile these with the actual objects and metadata before execution. The feature-selection script retains the original repeated use of the test partition and binomial-GLM blocks; review their suitability for the intended three-class analysis before treating the resulting performance as an independent validation.

## Translation and validation

- File and module names, Chinese comments, and Chinese plot labels were translated into English.
- Analysis expressions, model settings, thresholds, cluster assignments, and exploratory code were preserved.
- Two source syntax issues were repaired in the English copies: copied R-console `>` prompts were removed from the Cox history script, and a free-text marker note was converted to a comment in the T-cell/myeloid script.
- All 14 English scripts passed syntax parsing using R 4.5.2. Parsed expressions were compared with the supplied source after accounting for the documented string translations and syntax repairs; all 14 comparisons matched.
- Full biological analyses were not rerun. Successful parsing does not validate data availability, package compatibility, or analytical results.

`SCRIPT_MANIFEST.tsv` records SHA-256 hashes for the supplied source files and the corresponding English copies. Original local files were preserved.

## Attribution

The supplied methylation tutorial scripts retain their original author headers crediting Jianming Zeng and Bioinformatics Skill Tree, including their original tutorial links. Other source comments and code were retained where applicable. No new blanket license is asserted for third-party tutorial material.
