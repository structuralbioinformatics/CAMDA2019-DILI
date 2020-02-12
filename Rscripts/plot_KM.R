### Load packages ###
# For cmap analysis
library(cmapR)
library(ggplot2)
library(caret)
# For KMeans
library(ggfortify)
library(stringr)
library(FactoMineR)
library(factoextra)
library(dplyr)


### Define variables ###
place = "home" #home or work
remove.outliers = FALSE
outliers = c('daunorubicin', 'vorinostat')


if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PhD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/Projects/camda"
}


### Define files ###
#setwd("/sbi/users/interchange/emre/quim/camda")
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
drug_info_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-info.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
landmark_genes_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
# Output files
output_plot_pdf <- paste(main_directory, "results/plots/km_landmark.pdf", sep="/")
output_plot_png <- paste(main_directory, "results/plots/km_landmark.png", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


### Get landmark genes ###
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Subset GCT object ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_df <- subset.expression(gct, landmark_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_ind_df <- subset.expression(gct, landmark_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)


### Plot KMeans ###
expression_df$DILIrank <- expression_df$dilirank
expression_df$dilirank <- NULL
info_fields <- c("pert_iname", "severity", "DILIrank")
expression <- expression_df[,!colnames(expression_df) %in% info_fields]
info <- expression_df[,colnames(expression_df) %in% info_fields]
k3 <- kmeans(expression, 3)
p3 <- fviz_cluster(k3, data = expression, geom = "point") + ggtitle("K = 3")
data.plot <- p3
coord_clusters <- cbind(data.plot$data, info)
ggplot(data = coord_clusters, aes(x, y, color = DILIrank)) + 
  geom_point(aes(shape = DILIrank, color = DILIrank)) + 
  stat_ellipse() +
  labs(y=data.plot$labels$y, x = data.plot$labels$x)
  #geom_text(label = coord_clusters$pert_iname) +

### Plot KMeans (with ggfortify) ###
Cairo::CairoPNG(output_plot_png, dpi=300, width = 6, height = 5, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
autoplot(kmeans(expression,3), data = expression_df, colour = 'DILIrank', frame = TRUE)
dev.off()


