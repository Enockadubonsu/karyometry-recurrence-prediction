# From Nuclei to Patients: A Systematic Evaluation of Aggregation
Strategies for Bladder Cancer Recurrence Prediction

**Author:** Enock Adu Bonsu  
**Course:** BIOS 648 — High-Dimensional Health Data Analysis
and Machine Learning  
**Institution:** University of Arizona, Department of
Epidemiology and Biostatistics  


---

## Overview

This repository contains the complete analysis code, manuscript,
and figures for a systematic evaluation of patient-level
aggregation strategies for predicting bladder cancer recurrence
from karyometric nuclear morphometry data.

The central finding: **aggregation strategy matters more than
classifier choice** in multiple-instance learning problems
with hierarchical biomedical data.

---

## Key Results

| Model | MCE | Sensitivity | Specificity |
|-------|-----|-------------|-------------|
| Majority class baseline (M0) | 0.464 | 0.000 | 1.000 |
| Best individual model (M5: Lasso, Rep B) | 0.429 | 0.513 | 0.622 |
| Stacking meta-ensemble (META) | **0.321** | **0.641** | **0.711** |

**Leakage finding:** Nucleus-level cross-validation inflated
apparent performance by 33.9 percentage points (56% overoptimism)
relative to the valid patient-level estimate.

---

## Repository Structure

```
karyometry-recurrence-prediction/
│
├── code/
│   ├── karyometry_analysis.R       # Full LOOCV pipeline (Layers 1-3)
│   └── figures.R                   # All 7 figures
│
├── manuscript/
│   ├── BIOS648_Final_Manuscript.tex    # LaTeX source
│   └── BIOS648_Final_Manuscript.pdf    # Compiled PDF
│
├── figures/
│   ├── Fig1_MCE_forestplot.png
│   ├── Fig2_SensSPec.png
│   ├── Fig6_ConfusionMatrices.png
│   ├── FigS1_FeatureImportance.png
│   ├── FigS2_Leakage.png
│   └── FigS3_Heatmap_MCE.png
│
└── data/
    └── README.md                   # Data source and access instructions
```
---

## Methods Summary

Three patient-level representations:
- **Rep A:** Mean + SD + Q10/Q50/Q90 per feature (460-dim)
- **Rep B:** Score-guided softmax-weighted mean (92-dim)
- **Rep C:** Probability aggregation — mean, max, Q90 (3-dim)

Four classifiers: Lasso logistic regression, Random Forest,
PCA + LDA, Single hidden-layer neural network

Evaluation: Patient-level leave-one-out cross-validation (n=84)
with 95% Clopper-Pearson confidence intervals and McNemar's
test with Bonferroni correction.

---

## Requirements

```r
# R packages required
install.packages(c("glmnet", "randomForest", "MASS", "nnet"))
```

R version 4.3 or later.

---

## Citation

If you use this code or framework, please cite:

> Adu Bonsu, E. (2026). From Nuclei to Patients: A Systematic
> Evaluation of Aggregation Strategies for Bladder Cancer
> Recurrence Prediction Using Karyometric Morphometry Data.
> BIOS 648 Final Project, University of Arizona.