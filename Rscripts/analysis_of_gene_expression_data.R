#### Load packages ####
# For cmap analysis
library(cmapR)
library(ggplot2)


### Define variables ###
place = "home" #home or work
remove.outliers = FALSE
outliers = c('daunorubicin', 'vorinostat')
number.cv = 10
number.repetitions = 10
fraction_train = 0.7


#### Define files ####
# Data files
if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
# Output files
output_table_num_samples <- paste(main_directory, "camda_data/summary_num_samples.tsv", sep="/")


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
all_drugs <- unique(c(drug.dataset$drugs, drug.dataset$independent_drugs))


### Analyze number of samples for each concentration and time for all cell lines
cols <- c("cell_id", "pert_idose", "pert_itime", "num_drugs", "num_samples")
samples_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(samples_df) <- cols
cell_ids <- unique(gct@cdesc$cell_id)
cell_ids <- cell_ids[order(cell_ids)]

for (cell_id in cell_ids){
  idx <- which(gct@cdesc$cell_id==cell_id & gct@cdesc$pert_iname %in% all_drugs)
  id_selected_genes <- which(gct@rid %in% landmark_genes)
  x <- subset.gct(gct, cid=idx, rid=id_selected_genes) # Subset by drug
  concentrations <- unique(x@cdesc$pert_idose)
  concentrations <- concentrations[order(concentrations)]
  for (concentration in concentrations){
    idx_c <- which(x@cdesc$pert_idose==concentration)
    id_selected_genes <- which(x@rid %in% landmark_genes)
    x_c <- subset.gct(x, cid=idx_c, rid=id_selected_genes)
    times <- unique(x_c@cdesc$pert_itime)
    times <- times[order(times)]
    for (time in times){
      idx_t <- which(x_c@cdesc$pert_itime==time)
      id_selected_genes <- which(x_c@rid %in% landmark_genes)
      x_t <- subset.gct(x_c, cid=idx_t, rid=id_selected_genes)
      num_drugs <- length(unique(x_t@cdesc$pert_iname))
      num_doses <- length(x_t@cid)
      print(c(cell_id, concentration, time, num_drugs, num_doses))
      samples_df[nrow(samples_df)+1,] <- c(cell_id, concentration, time, num_drugs, num_doses)
    }
  }
}
write.table(samples_df, file = output_table_num_samples,row.names=FALSE, na="-",col.names=TRUE, sep="\t")
head(samples_df)


### Analyze number of samples for each concentration and time for cell line PHH ###
cell_id = "PHH"
pert_idose="10 µM"
pert_itime="24 h"
id_selected_samples <- which(gct@cdesc$cell_id == cell_id & gct@cdesc$pert_iname %in% all_drugs)
id_selected_genes <- which(gct@rid %in% landmark_genes)
#gct_subset <- subset.gct(gct, cid=id_selected_samples)
gct_subset <- subset.gct(gct, cid=id_selected_samples, rid=id_selected_genes)
table(gct_subset@cdesc$pert_idose, gct_subset@cdesc$pert_itime)


