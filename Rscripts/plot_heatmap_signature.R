### Load packages ###
#install.packages('circlize', dependencies = TRUE)
#BiocManager::install('ComplexHeatmap')
library(Cairo)
library(circlize) # For Heatmap()
library(cmapR)
library(ComplexHeatmap) # For Heatmap()
library(ggplot2)


### Define variables ###
place = "home" #home or work
remove.outliers = TRUE
outliers = c('daunorubicin', 'vorinostat')
cell_id = "PHH"
pert_idose="10 µM"
pert_itime="24 h"

if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}


### Define files ###
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected.txt", sep="/")
#wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_noout_notcorrected.txt", sep="/")
output_table <- paste(main_directory, "results/reverse_engineering/gene_test_phh_noout.tsv", sep="/")
#output_plot_pdf <- paste(main_directory, "results/plots/heatmap_gene_signature.pdf", sep="/")
#output_plot_png <- paste(main_directory, "results/plots/heatmap_gene_signature.png", sep="/")
output_plot_pdf <- paste(main_directory, "results/plots/heatmap_gene_signature_v2.pdf", sep="/")
output_plot_png <- paste(main_directory, "results/plots/heatmap_gene_signature_v2.png", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


### Get landmark genes ###
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]


### Subset drugs ###
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Subset by landmark genes, cell ID, dose and time ###
id_selected_samples <- which(gct@cdesc$cell_id == cell_id & gct@cdesc$pert_idose == pert_idose & gct@cdesc$pert_itime == pert_itime)
id_selected_genes <- which(gct@rid %in% landmark_genes)
gct_subset <- subset.gct(gct, cid=id_selected_samples, rid=id_selected_genes)


### Subset by DILI/NO-DILI drugs ###
# Subset DILI drugs
dili_concern_drugs <- union(drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs)
id_subset_dili <- which(gct_subset@cdesc$pert_iname %in% dili_concern_drugs)
gct_subset_dili <- subset.gct(gct_subset, cid=id_subset_dili)
# Subset NO-DILI drugs
id_subset_nodili <- which(gct_subset@cdesc$pert_iname %in% drug.dataset$no_concern_drugs)
gct_subset_nodili <- subset.gct(gct_subset, cid=id_subset_nodili)


### Load wilcoxon test analysis ###
heat_df = read.csv(output_table, header = TRUE, sep = "\t")
wilcox_df = read.csv(wilcox_file, header = FALSE, sep = "\t")
selected_genes <- unique(wilcox_df[,1])


### Get a vector of gene symbols ordered equal to the vector of gene IDs ###
genes_df <- data.frame(1:length(selected_genes), selected_genes)
colnames(genes_df)<- c("id", "pr_gene_id")
genes_df <- merge(genes_df, gene_info_df[,c("pr_gene_id", "pr_gene_symbol")], by="pr_gene_id")
genes_df <- genes_df[order(genes_df$id), ]
genes_df$id <- NULL


### Create new dataframe to store the matrix of expression ###
#cols <- c(drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs)
cols <- c(drug.dataset$most_concern_drugs, drug.dataset$no_concern_drugs)
rows <- genes_df$pr_gene_symbol
heat_df <- data.frame(matrix(ncol = length(cols), nrow=length(rows)))
colnames(heat_df) <- cols
rownames(heat_df) <- rows
gene_ids <- genes_df$pr_gene_id


### Introduce the expression in the matrix ###
##### First the DILI drugs
for (drug in drug.dataset$most_concern_drugs){
  id_drug <- which(gct_subset_dili@cdesc$pert_iname == drug)
  id_genes <- which(gct_subset_dili@rid %in% gene_ids)
  gct_drug <- subset.gct(gct_subset_dili, cid=id_drug, rid=id_genes)
  #means<- rowMeans(gct_drug@mat)
  medians <- apply(gct_drug@mat, 1, median) 
  heat_df[drug] <- medians
}
# for (drug in drug.dataset$less_concern_drugs){
#   id_drug <- which(gct_subset_dili@cdesc$pert_iname == drug)
#   id_genes <- which(gct_subset_dili@rid %in% gene_ids)
#   gct_drug <- subset.gct(gct_subset_dili, cid=id_drug, rid=id_genes)
#   #means<- rowMeans(gct_drug@mat)
#   medians <- apply(gct_drug@mat, 1, median) 
#   heat_df[drug] <- medians
# }
##### Then the No-DILI drugs
for (drug in drug.dataset$no_concern_drugs){
  id_drug <- which(gct_subset_nodili@cdesc$pert_iname == drug)
  id_genes <- which(gct_subset_nodili@rid %in% gene_ids)
  gct_drug <- subset.gct(gct_subset_nodili, cid=id_drug, rid=id_genes)
  #means<- rowMeans(gct_drug@mat)
  medians <- apply(gct_drug@mat, 1, median) 
  heat_df[drug] <- medians
}


### ComplexHeatmap removing expression lower than abs(1.5) ###
heat_df_for_plot <- apply(heat_df, c(1,2),  function(x) {ifelse(abs(x)>1.5, x, 0)})
#type <- c(rep("Most-DILI-Concern", length(drug.dataset$most_concern_drugs)), rep("Less-DILI-Concern", times=length(drug.dataset$less_concern_drugs)), rep("No-DILI-Concern", times=length(drug.dataset$no_concern_drugs)))
type <- c(rep("Most-DILI-Concern", length(drug.dataset$most_concern_drugs)), rep("No-DILI-Concern", times=length(drug.dataset$no_concern_drugs)))
#annotation <- HeatmapAnnotation(df = data.frame(type = type))
##### Plot PDF
#Cairo::CairoPDF(output_plot_pdf, width = 12, height = 8) # Save in PDF
#Cairo::CairoPDF(output_plot_pdf, width = 14, height = 10) # Save in PDF
Cairo::CairoPDF(output_plot_pdf, width = 18, height = 15) # Save in PDF
Heatmap(heat_df_for_plot, name = "MODZscore", km = 0, 
        col = colorRamp2(c(min(heat_df), 0, max(heat_df)), c("blue", "white", "red")), 
        #top_annotation = annotation, 
        cluster_rows=TRUE, cluster_columns=FALSE, 
        #row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 4),  
        #row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),  
        row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(fontsize = 14),  
        #column_split = c(rep("Most-DILI-Concern", length(drug.dataset$most_concern_drugs)), rep("Less-DILI-Concern", length(drug.dataset$less_concern_drugs)), rep("No-DILI-Concern", length(drug.dataset$no_concern_drugs))),
        column_split = c(rep("Most-DILI-Concern", length(drug.dataset$most_concern_drugs)), rep("No-DILI-Concern", length(drug.dataset$no_concern_drugs))),
        #column_title_gp = gpar(fontsize = 12),
        #column_title_gp = gpar(fill = c("#FFB448", "#FF6C65", "#99EE99"), fontsize = 12),
        column_title_gp = gpar(fill = c("#FFB448", "#FF6C65", "#99EE99"), fontsize = 16),
        column_gap = unit(5, "mm"),
        border = TRUE
)
dev.off()
##### Plot PNG
#Cairo::CairoPNG(output_plot_png, dpi=300, width = 10, height = 8, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
Cairo::CairoPNG(output_plot_png, dpi=300, width = 18, height = 15, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
Heatmap(heat_df_for_plot, name = "MODZscore", km = 0, 
        col = colorRamp2(c(min(heat_df), 0, max(heat_df)), c("blue", "white", "red")), 
        #top_annotation = annotation, 
        cluster_rows=TRUE, cluster_columns=FALSE, 
        #row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 4),  
        row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(fontsize = 14),  
        #column_split = c(rep("Most-DILI-Concern", length(drug.dataset$most_concern_drugs)), rep("Less-DILI-Concern", length(drug.dataset$less_concern_drugs)), rep("No-DILI-Concern", length(drug.dataset$no_concern_drugs))),
        column_split = c(rep("Most-DILI-Concern", length(drug.dataset$most_concern_drugs)), rep("No-DILI-Concern", length(drug.dataset$no_concern_drugs))),
        #column_title_gp = gpar(fontsize = 12),
        column_title_gp = gpar(fill = c("#FFB448", "#FF6C65", "#99EE99"), fontsize = 16),
        column_gap = unit(5, "mm"),
        border = TRUE
)
dev.off()
