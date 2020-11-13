### Load packages ###
library(ComplexHeatmap) # For Heatmap()
library(circlize) # For Heatmap()

### Define variables ###
place = "home" #home or work
#feature_names <- c("disgenet","guildify","disgenet_smiles","guildify_smiles")
#plot_names <- c("DisGeNET","GUILDify","DisGeNET + SMILES","GUILDify + SMILES")
feature_names <- c("disgenet")
plot_names <- c("DisGeNET")
type_genes = 'curated' # curated or all
ml_algorithms <- c("rf", "gbm") #c("rf", "gbm")
fontsize <- 10 # 8 for both rf and gbm

if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PhD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/Projects/camda"
}


### Define directories and files ###
output_dir <- paste(main_directory, "results/crossvalidation", sep="/")
phenotype2gene_file <- paste(main_directory, "guildify_data/phenotype2gene.tsv", sep="/")
disease2gene_file <- paste(main_directory, "guildify_data/disease2gene_disgenet_guildify.tsv", sep="/")
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")
guildify_dir <- paste(main_directory, "guildify_data/guildify_results/", sep="/")
validation_file <- paste(main_directory, "results/validations/3rd_round_validations.tsv", sep="/")


### Read DisGeNET and GUILDify phenotypes ###
disgenet2gene <- read.csv(disease2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
redundantphenotypes_df <- read.csv(redundantphenotypes_file, header=TRUE, sep="\t", stringsAsFactors = F)
categories <- c(paste("disgenet.",type_genes, sep = ""), paste("guildify.",type_genes, sep = ""))
redundantphenotypes_df <- redundantphenotypes_df[redundantphenotypes_df$category %in% categories,]
disgenet2gene <- disgenet2gene[!(disgenet2gene$diseaseid %in% redundantphenotypes_df$redundant.phenotype),]
phenotypes <- unique(disgenet2gene[,c('diseaseid', 'diseaseterm')])
phenotypes <- phenotypes[order(phenotypes$diseaseterm),]


### Read cross-validation results and calculate the mean of all results ###
#cols <- c("name","Accuracy", "Precision", "Sensitivity", "Specificity")
cols <- c("name","Accuracy", "Precision", "Sensitivity", "Specificity", "F1", "MCC")
results_phenotypes <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(results_phenotypes) <- cols

for (i in 1:nrow(phenotypes)){
  diseaseid <- phenotypes[i,c("diseaseid")]
  diseaseterm <- phenotypes[i,c("diseaseterm")]
  for (j in 1:length(feature_names)){
    feature_name <- feature_names[j]
    formal_feature_name <- plot_names[j]
    database <- strsplit(feature_name, "_")[[1]][1]
    type_analysis <- paste(database, type_genes, sep='.')
    disgenet2gene_df <- disgenet2gene[c("geneid", "diseaseid", "diseaseterm", type_analysis)]
    colnames(disgenet2gene_df) <- c("geneid", "diseaseid", "diseaseterm", "type_analysis")
    genes <- unique(disgenet2gene_df$geneid[disgenet2gene_df$type_analysis==1 & disgenet2gene_df$diseaseid==diseaseid])
    if (length(genes)>=10) {
      for (ml_algorithm in ml_algorithms){
        file_name <- paste("cv", feature_name, ml_algorithm, paste(diseaseid, ".txt", sep=""), sep="_")
        feature_file <- paste(output_dir, feature_name, file_name, sep="/")
        if (file.exists(feature_file)){
          feature_df <- read.csv(feature_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
          if (ncol(feature_df) == 5){
            feature_df$f1 <- apply(feature_df[,2:5], 1, function(y) calculate.f1.from.metrics(precision = y['precision'], sensitivity = y['sensitivity']))
            feature_df$mcc <- apply(feature_df[,2:5], 1, function(y) calculate.mcc.from.metrics(accuracy = y['accuracy'], precision = y['precision'], sensitivity = y['sensitivity'], specificity = y['specificity']))
          }
          feature_data <- feature_df[!feature_df$model=="pred.comb",]
          res <- apply(feature_data[,2:length(cols)], 2, mean)
          #plot_name <- paste(diseaseterm, " (",formal_feature_name,", ", toupper(ml_algorithm), ")", sep="")
          if (length(ml_algorithms) > 1){
            plot_name <- paste(diseaseterm, " (", toupper(ml_algorithm), ")", sep="")
          } else {
            plot_name <- paste(diseaseterm, sep="")
          }
          results_phenotypes[nrow(results_phenotypes)+1,] <- c(plot_name, res)
        }
      }
    }
  }
}
res_phen <- apply(results_phenotypes[,2:length(cols)], 2, as.numeric) # Convert into numeric the content of the matrix
res_mean <- apply(res_phen, 2, mean)
results_phenotypes[nrow(results_phenotypes)+1,] <- c("Mean", res_mean)



### Define output plots ###
feature_name <- feature_names[1]
if(length(ml_algorithms) == 1){
  ml_algorithm <- ml_algorithms[1]
  output_plot_pdf <- paste(main_directory, sprintf("results/plots/heatmap_results_phenotypes_%s_%s.pdf", feature_name, ml_algorithm), sep="/")
  output_plot_png <- paste(main_directory, sprintf("results/plots/heatmap_results_phenotypes_%s_%s.png", feature_name, ml_algorithm), sep="/")
  plot_width <- 8
  plot_height <- 5
} else {
  ml_algorithm <- paste(ml_algorithms, collapse = "_")
  output_plot_pdf <- paste(main_directory, sprintf("results/plots/heatmap_results_phenotypes_%s_%s.pdf", feature_name, ml_algorithm), sep="/")
  output_plot_png <- paste(main_directory, sprintf("results/plots/heatmap_results_phenotypes_%s_%s.png", feature_name, ml_algorithm), sep="/")
  plot_width <- 8
  plot_height <- 8
}



### Create a Heatmap with the results ###

# heat_df <- as.matrix(results_phenotypes[,2:5]) # Matrix needed to insert content into Heatmap
# colnames(heat_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity")
# heat_df <- apply(heat_df, 2, as.numeric) # Convert into numeric the content of the matrix
# rownames(heat_df) <- results_phenotypes[,1]
# Cairo::CairoPDF(output_plot_pdf, width = 6, height = 5) # Save in PDF
# Heatmap(heat_df, name = "results", km = 0, 
#         col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
#         na_col = "grey", # NA color
#         cluster_rows=FALSE, cluster_columns=FALSE, 
#         row_names_gp = gpar(fontsize = fontsize), column_names_gp =  gpar(fontsize = fontsize), column_names_rot = 45,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = fontsize))
#         }
# )
# dev.off()
# Cairo::CairoPNG(output_plot_png, dpi=300, width = 6, height = 5, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
# Heatmap(heat_df, name = "results", km = 0, 
#         col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
#         na_col = "grey", # NA color
#         cluster_rows=FALSE, cluster_columns=FALSE, 
#         row_names_gp = gpar(fontsize = fontsize), column_names_gp =  gpar(fontsize = fontsize), column_names_rot = 60,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = fontsize))
#         }
# )
# dev.off()

heat_df <- as.matrix(results_phenotypes[,2:6]) # Matrix needed to insert content into Heatmap
colnames(heat_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score")
heat_df <- apply(heat_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(heat_df) <- results_phenotypes[,1]
Cairo::CairoPDF(output_plot_pdf, width = plot_width, height = plot_height) # Save in PDF
h1 <- Heatmap(heat_df, name = "results", km = 0, 
              col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
              na_col = "grey", # NA color
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = fontsize), column_names_gp =  gpar(fontsize = fontsize), column_names_rot = 45,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = fontsize))
              }
)
mcc_df <- data.frame(matrix(ncol = 1, nrow = nrow(results_phenotypes)))
colnames(mcc_df) <- c("MCC")
mcc_df$MCC <- results_phenotypes[,7] # Matrix needed to insert content into Heatmap
mcc_df <- apply(mcc_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(mcc_df) <- results_phenotypes[,1]
h2 <- Heatmap(mcc_df, name = "MCC", km = 0, 
              #col = colorRamp2(c(min(mcc_df, na.rm = TRUE), 0, max(mcc_df, na.rm = TRUE)), c("yellow", "white", "#00b8ff")),
              col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 45,
              column_title_gp = gpar(fontsize = 10),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mcc_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
ht_list = h1 + h2
draw(ht_list, ht_gap = unit(2, "mm"))
dev.off()


Cairo::CairoPNG(output_plot_png, dpi=300, width = plot_width, height = plot_height, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
Heatmap(heat_df, name = "results", km = 0, 
        col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = fontsize), column_names_gp =  gpar(fontsize = fontsize), column_names_rot = 60,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = fontsize))
        }
)
h2 <- Heatmap(mcc_df, name = "MCC", km = 0, 
              #col = colorRamp2(c(min(mcc_df, na.rm = TRUE), 0, max(mcc_df, na.rm = TRUE)), c("yellow", "white", "#00b8ff")),
              col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 10),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mcc_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
ht_list = h1 + h2
draw(ht_list, ht_gap = unit(2, "mm"))
dev.off()

