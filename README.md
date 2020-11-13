# CAMDA Challenge: Predict drug toxicity
An ensemble learning approach for modeling the systems biology of drug-induced injury in human liver

## Getting Started

These instructions will guide you through the code used to participate at the CAMDA CMap drug safety challenge. All the scripts are written in R.

### Prerequisites

What packages you need to install. Some instructions are given at the R Markdowns.

```
cmapR
ComplexHeatmap
caret
RxnSim
```

### Outline

1. Gold standard analysis
2. Data collection
3. Identification of gene signatures to separate DILI and no-DILI drugs
4. Machine learning classification
5. Plots


### 1. Gold standard analysis

Script: `analysis_of_compounds.Rmd`.
We first analyzed the different DILIRank categories for each compound, removing the ambiguous drugs. Then, we analyzed the number of gene expression samples for each compound depending on the conditions of cell line, concentration and time. We written a summary of the analysis in the table `summary_compounds.tsv`.


### 2. Data collection

#### 2.1. CMap gene expression

We retrieved the gene expression data from the GCT object `CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda`, provided by CAMDA.

#### 2.2. Drug chemical structure

Script: `calculate_tanimoto.R`.
We calculated the similarity between all compounds, creating a matrix of chemical similarity. We used the R package RxnSim to calculate the similarity matrix using the Tanimoto distance.
The resulting table is called `tanimoto_smiles.tsv` (Supplementary Table 3).

#### 2.3. Drug target

The targets of the compounds (drugs) considered in the study were retrieved from three different databases: DGIdb, HitPick and SEA. 
* **DGIdb**: The compound names were used to retrieve the targets from the web server.
* **SEA and HitPick**: The SMILES of the compounds were used to get predictions in batch mode from the web servers.
The summary table is called `targets_dgidb_hitpick_sea.tsv` (Supplementary Table 4).

### 3. Identification of gene signatures to separate DILI and no-DILI drugs

#### 3.1. DisGeNET phenotype-gene associations
We compiled manually a list of 15 phenotypes related with DILI from DisGeNET.
Using the web server of DisGeNET, we retrieved phenotype-gene associations from expertly curated repositories (UniProt, the Comparative Toxicogenomics Database (CTD), ORPHANET, the Clinical Genome Resource (CLINGEN), the Genomics England PanelApp, the Cancer Genome Interpreter (CGI) and PsyGeNET).
We analyzed the genes of the phenotypes, checked which ones are redundant and put them in a table with the script `get_genes_from_disgenet_guildify.R`.
The output table of phenotype-gene associations is: `disease2gene_disgenet_guildify.tsv`.
The output table of redundant phenotypes is: `redundant_phenotypes.tsv`

#### 3.2. GUILDify expansions of phenotype-gene associations
Script: `gene_sets_phenotypes.R`.
We used the R package `guildifyR` to make the expansions of the genes associated to the 15 phenotypes. 

#### 3.3. Non-parametric test
Script: `gene_signature.R`.
First, we retrieved the gene expression samples from cell line PHH (liver primary cell), dose concentration 10 µM and time 24 h. 
Then, for each landmark gene, we calculated a Wilcoxon test between the gene expression values of the DILI and no-DILI drugs. We written a summary of the output values in the table `reverse_signature_phh_notcorrected_info.txt` (Supplementary Table 7).

#### 3.4. Correlated samples
Script: `analyse_top_correlated_samples.R` and `get_signature_from_top_correlated_samples.R`.
First, we calculate the Pearson's correlation between all the samples.
Then, we select the samples with a correlation threshold above 0.5, or otherwise the two most correlated samples.


### 4. Machine learning classification

We used different scripts to classify the drugs depending on the feature:

* **Landmark genes (CMap expression)**: `classify_by_landmark.R`
* **DisGeNET genes (CMap expression)**: `classify_by_disgenet.R`
* **GUILDify genes (CMap expression)**: `classify_by_guildify.R`
* **DILI landmark signature (CMap expression)**: `classify_by_signature.Rmd`
* **SMILES**: `classify_by_smiles.R`
* **Drug targets**: `classify_by_targets.R`
* **Combining DILI landmark signature & SMILES**: `classify_by_signature_smiles.Rmd`
* **Combining DisGeNET signature & SMILES**: `classify_by_disgenet_smiles.Rmd`
* **Combining GUILDify signature & SMILES**: `classify_by_guildify_smiles.Rmd`

The validations are in the folder `outputs/validations`.


### 5. Plots and Tables

* **Figure 1 & Supplementary Figure 8**: `check_means_results_with_F1_and_MCC.R`
* **Figure 2**: PowerPoint.
* **Figure 3 & Supplementary Figures 4 & 9**: `check_means_by_phenotype.R`
* **Figure 4**: `plot_heatmap_signature.R`
* **Supplementary Figure 1**: PowerPoint.
* **Supplementary Figure 2**: `plot_KM.R`
* **Supplementary Figure 3**: `plot_heatmap_phenotypes.R`
* **Supplementary Figure 5**: `plot_heatmap_smiles.R`
* **Supplementary Figure 6**: `plot_targets_percentage.R`
* **Supplementary Figure 7**: `check_means_results_correlated.R`
* **Table 1**: Manual.
* **Table 2**: `compare_genes.R`
* **Supplementary Table 1**: `create_data_files.R`
* **Supplementary Table 2**: `create_data_files.R`
* **Supplementary Table 3**: `calculate_tanimoto.R`
* **Supplementary Table 4**: Information retrieved using the web servers.
* **Supplementary Table 5**: `analysis_of_types_of_compounds.R`
* **Supplementary Table 6**: `comparison_genes.R`
* **Supplementary Table 7**: `gene_signature.R`


