## Pleiotropy Correction in pyVIPER

In transcriptional regulatory network analysis, pleiotropy occurs when a target gene is regulated by multiple transcription factors (TFs). This shared regulation can artificially inflate TF activity scores, because the same gene contributes signal to several regulators. As a result, TFs may appear active even when their inferred activity is driven by pleiotropic targets rather than true regulator-specific effects.

Pleiotropy correction addresses this issue by down-weighting or penalizing targets that are shared across multiple regulons, thereby reducing spurious enrichment and improving the specificity of inferred TF activity scores.

## What This Repository Implements

This repository contains a faithful Python re-implementation of the original shadowRegulon pleiotropy-correction algorithm from the R VIPER package, adapted for use in pyVIPER.
Specifically, the implementation:
- Reproduces the sample-wise pleiotropy correction logic used in VIPER (R)
- Penalizes target contributions based on regulon overlap (pleiotropy)
- Preserves normalized enrichment score (NES) behavior
- Produces numerically consistent results with the original R implementation when starting from the same z-scored expression data and regulons

The goal was not to approximate the method, but to replicate the exact algorithmic behavior of VIPER’s shadowRegulon in Python.

## Why This Matters

pyVIPER previously lacked a native pleiotropy-correction step, leading to:
- Overestimation of TF activity
- Inflated NES values for highly connected regulators
- Reduced agreement with VIPER (R) results
By implementing pleiotropy correction directly in Python:

TF activity inference becomes more biologically specific

Results are directly comparable between R VIPER and pyVIPER

Large-scale and HPC-based workflows can remain fully Python-native
