### Load packages ###
#install.packages('circlize', dependencies = TRUE)
#BiocManager::install('ComplexHeatmap')
library(Cairo)
library(circlize) # For Heatmap()
library(ComplexHeatmap) # For Heatmap()

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
tanimoto_file <- paste(main_directory, "additional_data/tanimoto_smiles.tsv", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
output_plot_pdf <- paste(main_directory, "results/plots/heatmap_smiles.pdf", sep="/")
output_plot_png <- paste(main_directory, "results/plots/heatmap_smiles.png", sep="/")
output_hist_pdf <- paste(main_directory, "results/plots/histogram_smiles.pdf", sep="/")
output_hist_png <- paste(main_directory, "results/plots/histogram_smiles.png", sep="/")
output_hist2_png <- paste(main_directory, "results/plots/histogram_smiles2.png", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(dilirank_file)

#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)

### Prepare SMILES data ###
tanimoto_df <- read.csv(tanimoto_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tanimoto_df$DILIConcern[tanimoto_df$DILIConcern == "No-DILI-concern"] <- "No-DILI-Concern" # Correct the 4 entries with lowercase concern
tanimoto_df$pert_iname = rownames(tanimoto_df)
colnames(tanimoto_df)[match("DILIConcern", colnames(tanimoto_df))] <- "dilirank"
tanimoto_ind_df <- tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$independent_drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]
tanimoto_df<-tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]

### Create a Heatmap with the results ###
heat_df <- as.matrix(tanimoto_df[,1:(ncol(tanimoto_df)-4)]) # Matrix needed to insert content into Heatmap

Cairo::CairoPDF(output_plot_pdf, width = 15, height = 15) # Save in PDF
Heatmap(heat_df, name = "results", km = 0, 
        col = colorRamp2(c(min(heat_df), 0.25, 1), c("white", "yellow", "red")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 6), column_names_gp =  gpar(fontsize = 6), column_names_rot = 45,
        column_title_gp = gpar(fontsize = 16),
        row_split = c(rep("Most-DILI-Concern", length(tanimoto_df$pert_iname[tanimoto_df$dilirank == "Most-DILI-Concern"])), rep("Less-DILI-Concern", length(tanimoto_df$pert_iname[tanimoto_df$dilirank == "Less-DILI-Concern"])), rep("No-DILI-Concern", length(tanimoto_df$pert_iname[tanimoto_df$dilirank == "No-DILI-Concern"]))),
        column_split = c(rep("Most-DILI-Concern", length(tanimoto_df$pert_iname[tanimoto_df$dilirank == "Most-DILI-Concern"])), rep("Less-DILI-Concern", length(tanimoto_df$pert_iname[tanimoto_df$dilirank == "Less-DILI-Concern"])), rep("No-DILI-Concern", length(tanimoto_df$pert_iname[tanimoto_df$dilirank == "No-DILI-Concern"]))),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm")
)
dev.off()


### Calculate means between groups ###
tanimoto_most_df <- tanimoto_df[rownames(tanimoto_df) %in% drug.dataset$most_concern_drugs, colnames(tanimoto_df) %in% drug.dataset$most_concern_drugs]
tanimoto_less_df <- tanimoto_df[rownames(tanimoto_df) %in% drug.dataset$less_concern_drugs, colnames(tanimoto_df) %in% drug.dataset$less_concern_drugs]
tanimoto_no_df <- tanimoto_df[rownames(tanimoto_df) %in% drug.dataset$no_concern_drugs, colnames(tanimoto_df) %in% drug.dataset$no_concern_drugs]
tanimoto_most_no_df <- tanimoto_df[rownames(tanimoto_df) %in% drug.dataset$most_concern_drugs, colnames(tanimoto_df) %in% drug.dataset$no_concern_drugs]
tanimoto_most_less_df <- tanimoto_df[rownames(tanimoto_df) %in% drug.dataset$most_concern_drugs, colnames(tanimoto_df) %in% drug.dataset$less_concern_drugs]
tanimoto_less_no_df <- tanimoto_df[rownames(tanimoto_df) %in% drug.dataset$less_concern_drugs, colnames(tanimoto_df) %in% drug.dataset$no_concern_drugs]

tanimoto_most_val <- tanimoto_most_df[upper.tri(tanimoto_most_df, diag = FALSE)]
tanimoto_less_val <- tanimoto_less_df[upper.tri(tanimoto_less_df, diag = FALSE)]
tanimoto_no_val <- tanimoto_no_df[upper.tri(tanimoto_no_df, diag = FALSE)]
tanimoto_most_no_val <- tanimoto_most_no_df[upper.tri(tanimoto_most_no_df, diag = FALSE)]
tanimoto_most_less_val <- tanimoto_most_less_df[upper.tri(tanimoto_most_less_df, diag = FALSE)]
tanimoto_less_no_val <- tanimoto_less_no_df[upper.tri(tanimoto_less_no_df, diag = FALSE)]

mean(tanimoto_most_val)
mean(tanimoto_less_val)
mean(tanimoto_no_val)
mean(tanimoto_most_no_val)
mean(tanimoto_most_less_val)
mean(tanimoto_less_no_val)

most_df <- data.frame(Tanimoto=tanimoto_most_val, DILIrank="Most vs. Most")
less_df <- data.frame(Tanimoto=tanimoto_less_val, DILIrank="Less vs. Less")
no_df <- data.frame(Tanimoto=tanimoto_no_val, DILIrank="No vs. No")
most_no_df <- data.frame(Tanimoto=tanimoto_most_no_val, DILIrank="Most vs. No")
most_less_df <- data.frame(Tanimoto=tanimoto_most_less_val, DILIrank="Most vs. Less")
less_no_df <- data.frame(Tanimoto=tanimoto_less_no_val, DILIrank="Less vs. No")
all_df <- rbind(most_df, less_df, no_df)
most_no_comparison_df <- rbind(most_df, most_no_df, no_df)

Cairo::CairoPNG(output_hist_png, dpi=300, width = 5, height = 4, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
ggplot(all_df, aes(x=Tanimoto, fill=DILIrank)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, binwidth=0.05, color='black', position="identity") +
  scale_fill_manual(values=c("#FF6C65", "#FFB448", "#99EE99")) +
  labs(x="Tanimoto index", y="Density") +
  geom_density(alpha=0.4)
dev.off()

Cairo::CairoPNG(output_hist2_png, dpi=300, width = 5, height = 4, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
ggplot(most_no_comparison_df, aes(x=Tanimoto, fill=DILIrank)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, binwidth=0.05, color='black', position="identity") +
  scale_fill_manual(values=c("#FF6C65", "#A3DCF8", "#99EE99")) +
  labs(x="Tanimoto index", y="Density") +
  geom_density(alpha=0.4)
dev.off()



