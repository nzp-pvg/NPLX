# 🫀 NPLX Reproducibility Hub

This repository provides a reproducible analytical framework for identifying and validating an **NPL-centered plaque-core macrophage program** associated with foam remodeling in human atherosclerosis. The workflow integrates bulk transcriptomics, DNA methylation, single-cell RNA sequencing, cross-cohort validation, macrophage-state inference, trajectory analysis, virtual gene perturbation, cell–cell communication analysis, and exploratory translational validation.

The repository is organized for traceability, parameter transparency, modular execution, and figure-level regeneration.

![Project](https://img.shields.io/badge/project-NPL--centered%20atherosclerosis-005493)
![Analysis](https://img.shields.io/badge/analysis-multi--omics-712238)
![Reproducibility](https://img.shields.io/badge/workflow-reproducible-F5BD4D)
![Language](https://img.shields.io/badge/code-R%20%7C%20Python%20%7C%20Bash-8F979A)

> 🧠 **Core question:** Is the NPL-centered program reproducibly localized to plaque-core macrophages, associated with foam remodeling, supported across molecular layers, and detectable in independent translational datasets?

---

## 🎯 What you can reproduce here

- ✅ Bulk atherosclerosis transcriptomic preprocessing and limma differential analysis
- ✅ Cross-cohort consensus candidate construction across public human datasets
- ✅ Mechanism gene-set construction and bulk mechanism scoring
- ✅ DNA methylation preprocessing and methylation differential analysis
- ✅ Multi-omics candidate integration across transcriptomic and methylation evidence
- ✅ Single-cell cell-type localization in plaque-core and adjacent tissue contexts
- ✅ Macrophage-state scoring, foam/TREM2-like versus C1Q-like comparisons, and pseudobulk analysis
- ✅ Macrophage continuum and pseudotime-style transition analysis
- ✅ NPL compact-module construction, cross-layer scoring, and bulk coherence assessment
- ✅ Virtual gene perturbation using scTenifoldKnk-compatible workflows
- ✅ Cell–cell communication analysis centered on NPL-module-high macrophages
- ✅ Blood-axis and translational signature analysis
- ✅ Cross-cohort robustness checks and figure-level regeneration

⚠️ Virtual perturbation, trajectory inference, cell–cell communication, and translational classifier outputs are computational evidence. They should not be interpreted as direct experimental perturbation, causal proof, or clinical validation.

---

## 🚀 Quick start

### 1. Project paths

The main project code is organized under:

```text
script/R/       # Main R analysis pipeline
script/python/  # Metadata and cohort preparation helpers
script/bash/    # Environment, download, and pipeline entrypoints
scripts/figures/ # Figure-specific plotting scripts
```

Run commands from the project root:

```bash
cd ~
```

### 2. Bootstrap and data preparation

```bash
bash script/bash/00_env_snapshot.sh
bash script/bash/01_tree_snapshot.sh
bash script/bash/02_download_geo_soft.sh
bash script/bash/03_extract_soft_metadata.sh
bash script/bash/04_download_series_matrix.sh
bash script/bash/05_download_platform_tables.sh
bash script/bash/run_main_analysis.sh
```

### 3. Recommended analytical order

```text
01–03  Platform annotation and bulk-expression preparation
04–07  Bulk differential analysis, consensus candidates, and mechanism scores
08–09  GSE46394 methylation preparation and limma analysis
10–11  GSE221911 expression preparation and blood-axis analysis
12–16  Single-cell marker sets, localization, support, and candidate shortlist
17–21  Macrophage-state definition, refinement, and transition scoring
22–28  NPL module construction, cross-layer scoring, coherence, nulls, and virtual perturbation
29–31  Macrophage continuum, pseudobulk, and cell-type-restricted program analysis
32–35  Blood translation, methylation integration, translational support, and robustness
36–37  Single-cell state validation and cell–cell communication analysis
```

---

## 🧬 Data and cohorts

| Data layer | Dataset or source | Analytical role |
|---|---|---|
| Bulk plaque transcriptomics | `GSE43292`, `GSE28829`, `GSE100927` | Discovery, external validation, and plaque-context comparison |
| Single-cell plaque data | `GSE159677` | Cell-type localization and macrophage-state analysis |
| Single-cell support data | `GSE155512` | Independent cell-state and macrophage support |
| DNA methylation | `GSE46394` | Methylation-level support for integrated candidates |
| Blood transcriptomics | `GSE221911` | Translational blood-axis and classifier analysis |
| Translational support | `GSE21545` | Additional cross-cohort support |
| Cell–cell communication | CellChatDB human reference | Ligand–receptor communication context |

Public GEO inputs and processed cohort manifests are stored under `data/`. Data provenance is recorded in `data/metadata/` and `Log/`.

---

## 🧹 Core analytical modules

### 🧬 Bulk transcriptomics

- Platform-specific annotation tables are prepared before expression harmonization.
- Probe-level matrices are converted to gene-level matrices.
- Cohort-specific differential expression is performed with `limma`.
- Consensus candidates are constructed across independent atherosclerosis datasets.
- Mechanism gene sets are scored across bulk expression layers.

### 🧪 DNA methylation

- `GSE46394` beta-value matrices and phenotype data are prepared separately from expression data.
- Methylation-level differential analysis is performed with `limma`.
- Methylation evidence is integrated with transcriptomic candidate scores.

### 🔬 Single-cell localization

- `GSE159677` is used for plaque-core and adjacent-tissue cell-type localization.
- Major vascular cell types include macrophages, endothelial cells, smooth muscle cells, fibroblast/mesenchymal cells, T cells, B cells, and mast cells.
- The NPL compact module is evaluated across cell types and macrophage states.
- `GSE155512` provides independent single-cell support.

### 🫧 Macrophage foam-remodeling axis

- Macrophage states include C1Q-like, inflammatory, IFN-like, and FOAM/TREM2-like programs.
- NPL-centered modules are constructed from core genes and foam-remodeling neighbors.
- Macrophage continuum, pseudobulk, patient-pair, and state-level analyses are retained as separate evidence layers.

### 🧠 Virtual perturbation and communication

- NPL-module-high and NPL-module-low macrophage groups are compared computationally.
- Virtual gene perturbation is performed using the project scTenifoldKnk-compatible workflow.
- Cell–cell communication is summarized using a cached human CellChat reference and sample-aware macrophage grouping.

---

## ⚙️ Key parameters and interpretation rules

- **NPL compact module:** `NPL`, `FABP5`, `GPNMB`, `APOC1`, `PLA2G7`, `SPP1`, `CD36`, `CYP27A1`, and `APOE`.
- **Primary macrophage comparison:** FOAM/TREM2-like versus C1Q-like macrophage states.
- **External validation principle:** independent cohorts are used as support layers rather than being treated as interchangeable replicas.
- **Patient-pair analyses:** paired core-versus-adjacent comparisons are preserved where sample metadata allow.
- **Null and robustness checks:** module null benchmarking, cross-cohort robustness, and perturbation sensitivity are retained in the result tables.
- **Statistical caution:** descriptive single-cell and communication results should be interpreted with their cohort size, sampling structure, and cell-level dependence in mind.

---

## 🗂️ Repository layout

```text
data/
├── raw/                         # GEO SOFT, series matrices, platform tables, references
├── processed/                   # Bulk, methylation, phenotype, annotation, and single-cell inputs
└── metadata/                    # Dataset, platform, sample, and study-design manifests

res/
├── qc/                          # QC and pipeline summaries
├── tables/bulk/                 # Bulk DEG, signature, and classifier tables
└── tables/mechanism/            # Methylation, single-cell, NPL, state, and communication tables

script/
├── R/                           # Numbered analysis scripts and project configuration
├── python/                      # Metadata and phenotype preparation helpers
└── bash/                        # Bootstrap, download, and execution entrypoints

scripts/figures/                 # Figure-specific regeneration scripts
figures/                         # Final figure exports
results/                         # Figure-specific and validation outputs
Log/                             # Run history and environment records
manuscript/                      # Notes, proposal files, and manuscript materials
```

---

## 🎨 Figure regeneration

The main figure-specific script currently available is:

```bash
Rscript scripts/figures/plot_Figure3_single_cell_localization_compact_NPL_module.R
```

This script regenerates the NPL compact-module single-cell localization figure and exports panel-level data to `results/figure3/`. Other numbered R scripts write intermediate tables and figures to `res/`, `results/`, and `figures/` according to their module-specific configuration.

---

## 📦 Important outputs

- `res/tables/bulk/` — bulk differential-expression and translational classifier tables
- `res/tables/mechanism/` — methylation, NPL-module, macrophage-state, pseudobulk, and communication tables
- `res/qc/` — dataset, matrix, methylation, single-cell, and robustness summaries
- `results/figure3/` — panel-level exports and key numbers for the compact NPL localization figure
- `figures/main/` — main figure exports

---

## ♻️ Reproducibility and limitations

This project preserves raw public inputs, processed matrices, phenotype manifests, QC summaries, numbered scripts, and run history where available. Results should be interpreted alongside the following limitations:

- public cohorts differ in tissue source, platform, phenotype definition, and sample size;
- single-cell analyses are sensitive to cell-state annotation, sampling structure, and cell-level dependence;
- methylation and transcriptomic evidence are integrated as complementary association layers, not as proof of regulatory causality;
- trajectory inference and virtual perturbation are computational hypotheses that require experimental validation;
- CellChat-like communication scores depend on the reference database and sample grouping rules;
- translational blood analyses provide peripheral support rather than direct plaque-tissue confirmation.

Historical outputs and run records should be preserved when regenerating figures. New figure revisions should use a new version suffix rather than overwriting prior outputs.

---

## 📖 Citation and license

If you use this repository, please cite the associated manuscript once the final publication reference is available.

The repository currently does not specify a software license. Please contact the authors before reuse or redistribution until a license is added.
