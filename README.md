# CAMDA Challenge: Predict drug toxicity
An ensemble learning approach for modeling the systems biology of drug-induced injury in human liver

## Getting Started

These instructions will guide you through the code used to participate at the CAMDA CMap drug safety challenge. All the scripts are written in R and most of them are R Markdowns that improve the readability of the code.

### Prerequisites

What packages you need to install. Some instructions are given at the R Markdowns.

```
cmapR
ComplexHeatmap
caret
```

### Outline

1. Analyze the compounds
2. Identify genes that separate better Most/Less vs. DILI drugs
3. Classify drugs

### 1. Analyze the compounds

In the script `analysis_of_compounds.Rmd`, we first analyzed the different DILIRank categories for each compound, removing the ambiguous drugs. Then, we analyzed the number of gene expression samples for each compound depending on the conditions of cell line, dose and time. We written a summary of the analysis in the table `summary_compounds.tsv`.