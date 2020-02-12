### Load packages ###
library(cmapR)
library(ggplot2)
library(caret)


### Define variables ###
place = "home" #home or work
remove.outliers = TRUE
outliers = c('daunorubicin', 'vorinostat')
number.cv = 10
number.repetitions = 10
fraction_train = 0.7
# For 1 drug or more
#cutoff_perc_targets = 0
#width_pdf = 200
#height_pdf = 7
#name_plot_pdf = "barplot_targets_percentage_complete.pdf"
#name_plot_png = "barplot_targets_percentage_complete.png"
# For 20 drugs or more
cutoff_perc_targets = 5
cutoff_targets = 20
width_pdf = 10
height_pdf = 5
name_plot_pdf = "barplot_targets_percentage.pdf"
name_plot_png = "barplot_targets_percentage.pdf"


if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PhD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/Projects/camda"
}

### Define files ###
# Data files
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
drug_info_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-info.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
landmark_genes_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
# Specific data files
targets_file <- paste(main_directory, "additional_data/targets/targets_dgidb_hitpick_sea.tsv", sep="/")
# Output files
output_plot_pdf <- paste(main_directory, "results/plots", name_plot_pdf, sep="/")
output_plot_png <- paste(main_directory, "results/plots", name_plot_png, sep="/")
output_venn_pdf <- paste(main_directory, "results/plots/venn_targets.pdf", sep="/")
output_venn_png <- paste(main_directory, "results/plots/venn_targets.png", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Prepare targets info ###
# Read targets and map them to dilirank info (the dilirank info on targets is wrong!)
targets_df <- read.csv(targets_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(targets_df)[match("drug", colnames(targets_df))] <- "pert_iname" # Change name to pert_iname
targets_ind_df <- targets_df[targets_df$pert_iname %in% drug.dataset$independent_drugs,] # Get independent set
targets_df <- merge(x = targets_df, y = drug.dataset$dilirank_df[c("pert_iname", "DILIConcern", "Severity.Class")], by = "pert_iname")
targets_df$DILI <- NULL
targets_df$severity <- NULL
targets_ind_df$DILI <- NULL
targets_ind_df$severity <- NULL
colnames(targets_df)[match("DILIConcern", colnames(targets_df))] <- "dilirank"
colnames(targets_df)[match("Severity.Class", colnames(targets_df))] <- "severity"
targets_df <- targets_df[targets_df$pert_iname %in% drug.dataset$drugs,]


### Prepare balanced machine learning datasets ###
datasets.list <- prepare.balanced.datasets(targets_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)


### Prepare the data for the plot ###
targets <- colnames(targets_df)[!(colnames(targets_df) %in% c("pert_iname", "dilirank", "severity"))]
plot_df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(plot_df) <- c("target","perc_dili", "perc_nodili", "perc_ind", "perc_total")
for (target in targets){
  tar_dili <- targets_df[targets_df$pert_iname %in% c(drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs), colnames(targets_df)==target]
  tar_nodili <- targets_df[targets_df$pert_iname %in% drug.dataset$no_concern_drugs, colnames(targets_df)==target]
  tar_ind <- targets_ind_df[targets_ind_df$pert_iname %in% drug.dataset$independent_drugs, colnames(targets_ind_df)==target]
  tar_total <- c(tar_dili, tar_nodili, tar_ind)
  perc_dili <- length(which(tar_dili == 1)) / length(c(drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs)) * 100
  perc_nodili <- length(which(tar_nodili == 1)) / length(drug.dataset$no_concern_drugs)  * 100
  perc_ind <- length(which(tar_ind == 1)) / length(drug.dataset$independent_drugs) * 100
  perc_total <- (length(which(tar_dili == 1))+length(which(tar_nodili == 1))+length(which(tar_ind == 1))) / length(c(drug.dataset$drugs, drug.dataset$independent_drugs)) * 100
  plot_df[nrow(plot_df)+1,] <- c(target, as.numeric(perc_dili), as.numeric(perc_nodili), as.numeric(perc_ind), as.numeric(perc_total))
}
plot_df <- transform(plot_df, perc_dili = as.numeric(perc_dili), perc_nodili = as.numeric(perc_nodili), perc_ind = as.numeric(perc_ind), perc_total = as.numeric(perc_total))
plot_df <- plot_df[order(plot_df$perc_ind, plot_df$perc_total, plot_df$perc_dili, decreasing = T), ] # Order the data
#plot_filt_df <- plot_df[plot_df$perc_total>cutoff_perc_targets,]
plot_filt_df <- plot_df[c(1:cutoff_targets),]
plot_filt_df$target <- factor(plot_filt_df$target, levels = plot_filt_df$target)  # to retain the order in plot.
targets_order <- plot_filt_df$target


### combine the num columns into a single column with separate rows for each type of drug; assign to new vector ###
library(tidyr)
plot_filt_df$perc_total <- NULL
plot_reordered_df <- gather(plot_filt_df,type_drug,perc_drugs,-target)
plot_reordered_df$type_drug <- factor(plot_reordered_df$type_drug, levels = c("perc_dili", "perc_nodili", "perc_ind")) # Change the order of the stacked groups

# Format numbers as percentages
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
#plot_reordered_df$perc_drugs <- percent(plot_reordered_df$perc_drugs, digits = 2)

### Draw plot ###
Cairo::CairoPDF(output_plot_pdf, width = width_pdf, height = height_pdf) # Save in PDF
theme_set(theme_bw())
# Stacked barplot
ggplot(plot_reordered_df, aes(fill=type_drug, y=perc_drugs, x=target)) + 
  #geom_bar(position="stack", stat="identity") +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_text(aes(label=perc_drugs), size = 3, position = position_stack(vjust = 0.5))+
  #labs(title="Number of drugs per target", x="Targets", y="Number of drugs") + 
  labs(x="Targets", y="% of drugs") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #scale_fill_discrete(name = "Number of drugs")
  scale_fill_discrete(name = "Number of drugs", labels = c("% DILI drugs", "% No-DILI drugs", "% independent drugs"))
dev.off()


