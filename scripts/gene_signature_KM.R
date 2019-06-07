
### Load packages ###
# For cmap analysis
library(cmapR)
library(ggplot2)
library(caret)
# For KMeans
library(stringr)
library(FactoMineR)
library(factoextra)
library(dplyr)

### Define files ###
#setwd("/sbi/users/interchange/emre/quim/camda")
drugs_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-pert_iname.rda"
drug_info_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-info.rda"
dilirank_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda"
expression_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda"
landmark_genes_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt"
gene_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_gene_info.txt"
cell_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_cell_info.txt"

### Load files ###
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR

### Functions ###
get.expression.matrix<-function(gct_subset, selected_genes, drank.sel) {
  id_subset_sig <- which(gct_subset@rdesc$id %in% selected_genes)
  gct_subset_sig <- subset.gct(gct_subset, rid=id_subset_sig)
  # Get matrix in form of data frame
  sigmatrix <- data.frame(t(gct_subset_sig@mat), check.names = F, stringsAsFactors = F)
  sigmatrix$pert_iname <- gct_subset_sig@cdesc$pert_iname
  sigmatrix <- sigmatrix[sigmatrix$pert_iname %in% drank.sel$pert_iname,]
  sigmatrix$dili <- drank.sel$vDILIConcern[match(sigmatrix$pert_iname, drank.sel$pert_iname)]
  return(sigmatrix);
}

### Get landmark genes ###
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]

### Subset drugs ###
drank.sel$DILIConcern[drank.sel$DILIConcern == "No-DILI-concern"] <- "No-DILI-Concern" # Correct the 4 entries with lowercase concern
# We remove the first letter from vDILIConcern entries (except Ambiguous)
drank.sel$vDILIConcern[drank.sel$vDILIConcern != "Ambiguous DILI-concern"] <- sub('.', '', drank.sel$vDILIConcern[drank.sel$vDILIConcern != "Ambiguous DILI-concern"])
# We subset the table by those drugs that have the same value in vDILIConcern and DILIConcern columns
dilirank_df <- drank.sel[drank.sel$vDILIConcern == drank.sel$DILIConcern,]
dili_concern_drugs <- dilirank_df$pert_iname[dilirank_df$DILIConcern == "Most-DILI-Concern" | dilirank_df$DILIConcern == "Less-DILI-Concern" ]
no_concern_drugs <- dilirank_df$pert_iname[dilirank_df$DILIConcern == "No-DILI-Concern"]

### Subset GCT ###
# Subset the GCT object by cell ID = PHH, dose concentration 10 µM and time 24h
id_subset <- which(gct@cdesc$pert_idose == "10 µM" & gct@cdesc$pert_itime == "24 h" & gct@cdesc$cell_id == "PHH")
gct_subset <- subset.gct(gct, cid=id_subset)
# Subset the GCT object by DILI-Concern drugs
id_dili <- which(gct_subset@cdesc$pert_iname %in% dili_concern_drugs)
gct_dili <- subset.gct(gct_subset, cid=id_dili)
# Subset the GCT object by No-DILI-Concern drugs
id_no <- which(gct_subset@cdesc$pert_iname %in% no_concern_drugs)
gct_no <- subset.gct(gct_subset, cid=id_no)

### Load wilcoxon test analysis ###
wilcox_file <- "/home/quim/PHD/Projects/camda/results/gene_test_landmark_phh_10_24_mixlessmost.tsv"
wilcox_df = read.csv(wilcox_file, header = TRUE, sep = "\t")
selected_genes1 <- unique(wilcox_df$gene_id[wilcox_df$p.value<0.05 & wilcox_df$landmark==TRUE])
selected_genes2 <- unique(wilcox_df$gene_id[wilcox_df$p.value<0.1 & wilcox_df$landmark==TRUE])

### Load wilcoxon test analysis removing outliers ###
wilcox_file <- "/home/quim/PHD/Projects/camda/results/gene_test_landmark_phh_10_24_noout_mixlessmost.tsv"
wilcoxout_df = read.csv(wilcox_file, header = TRUE, sep = "\t")
selected_genes_out1 <- unique(wilcoxout_df$gene_id[wilcoxout_df$p.value<0.05 & wilcoxout_df$landmark==TRUE])
selected_genes_out2 <- unique(wilcoxout_df$gene_id[wilcoxout_df$p.value<0.1 & wilcoxout_df$landmark==TRUE])

### Get signature matrices for plotting KMeans ###
sigmatrix1 <- get.expression.matrix(gct_subset, selected_genes1, drank.sel)
sigmatrix2 <- get.expression.matrix(gct_subset, selected_genes2, drank.sel)
sigmatrixout1 <- get.expression.matrix(gct_subset, selected_genes_out1, drank.sel)
sigmatrixout2 <- get.expression.matrix(gct_subset, selected_genes_out2, drank.sel)

### Plot KMeans ###
plots <- list()
matrices <- list(sigmatrix1, sigmatrix2, sigmatrixout1, sigmatrixout2)
titles <- list("Gene signature (p.val<0.05)", "Gene signature (p.val<0.1)", "Gene signature (no outliers, p.val<0.05)", "Gene signature (no outliers, p.val<0.1)")
for(i in 1:length(matrices)){
  sigmatrix <- matrices[[i]]
  expression <- sigmatrix[,!colnames(sigmatrix) %in% c("pert_iname", "dili")]
  info <- sigmatrix[,colnames(sigmatrix) %in% c("pert_iname", "dili")]
  k3 <- kmeans(expression, 3)
  p3 <- fviz_cluster(k3, data = expression, geom = "point") + ggtitle("K = 3")
  data.plot <- p3
  coord_clusters <- cbind(data.plot$data, info)
  plots[[i]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(color = dili)) + #geom_text(label = coord_clusters$compound) +
    ggtitle(titles[[i]])
}

pdf("/home/quim/PHD/Projects/camda/results/kmeans/plots_phh_10_24.pdf", onefile = T)
for(i in 1:length(plots)){
  print(plots[[i]])
}
dev.off()

