### Load packages ###
library(ComplexHeatmap) # For Heatmap()
library(circlize) # For Heatmap()

### Define variables ###
place = "home" #home or work
type_genes = 'curated' # curated or all
ml_algorithms <- c("rf") #c("rf", "gbm")
disgenet_phenotype <- "C0023892" #FALSE
if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}
feature_names <- c("landmark",
                   "disgenet",
                   "guildify",
                   "signature",
                   #"signature_by_tas",
                   #"peng",
                   "targets",
                   "smiles",
                   "disgenet_smiles",
                   "guildify_smiles",
                   "signature_smiles"
)
plot_names <- c("Landmark",
                "DisGeNET",
                "GUILDify",
                "DILI landmark",
                #"Peng",
                #"Reverse eng. TAS",
                "Targets",
                "SMILES",
                "DisGeNET + SMILES",
                "GUILDify + SMILES",
                "DILI landmark + SMILES"
)


### Define directories and files ###
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
phenotype2gene_file <- paste(main_directory, "guildify_data/phenotype2gene.tsv", sep="/")
disease2gene_file <- paste(main_directory, "guildify_data/disease2gene_disgenet_guildify.tsv", sep="/")
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")
guildify_dir <- paste(main_directory, "guildify_data/guildify_results/", sep="/")
validation_file <- paste(main_directory, "results/validations/3rd_round_validations.tsv", sep="/")
output_dir <- paste(main_directory, "results/crossvalidation", sep="/")


### Load files ###
source(functions_file)


### Read DisGeNET and GUILDify phenotypes ###
disgenet2gene <- read.csv(disease2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
redundantphenotypes_df <- read.csv(redundantphenotypes_file, header=TRUE, sep="\t", stringsAsFactors = F)
categories <- c(paste("disgenet.",type_genes, sep = ""), paste("guildify.",type_genes, sep = ""))
redundantphenotypes_df <- redundantphenotypes_df[redundantphenotypes_df$category %in% categories,]
disgenet2gene <- disgenet2gene[!(disgenet2gene$diseaseid %in% redundantphenotypes_df$redundant.phenotype),]
phenotypes <- unique(disgenet2gene[,c('diseaseid', 'diseaseterm')])
phenotypes <- phenotypes[order(phenotypes$diseaseterm),]


### Read the independent validation results ###
vm = read.csv(validation_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)


### Read cross-validation results and calculate the mean of all results ###
cols <- c("name","Accuracy", "Precision", "Sensitivity", "Specificity", "F1", "MCC", "Accuracy", "Precision", "Sensitivity", "Specificity", "F1", "MCC")
results <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(results) <- cols

for (i in 1:length(feature_names)){
  feature_name <- feature_names[i]
  formal_feature_name <- plot_names[i]
  for (ml_algorithm in ml_algorithms){
    # Get the independent validation results
    val_file_name <- paste("val.JAguirre_predictions", feature_name, paste(ml_algorithm, ".txt", sep=""), sep="_")
    valres <- vm[vm["Name"] == val_file_name,][,2:7]
    if (nrow(valres)==0){
      valres <- rep(NA, 6)
    }
    
    # Get the cross-validation results
    if (feature_name %in% c("disgenet", "disgenet_smiles", "guildify", "guildify_smiles")){
      cols <- c("name","accuracy", "precision", "sensitivity", "specificity", "f1", "mcc")
      results_phenotypes <- data.frame(matrix(ncol = length(cols), nrow = 0))
      colnames(results_phenotypes) <- cols
      if (length(ml_algorithms) > 1){
        plot_name <- paste(formal_feature_name, " (", toupper(ml_algorithm), ")", sep="")
      } else {
        plot_name <- paste(formal_feature_name, sep="")
      }
      if (disgenet_phenotype != FALSE){
        val_file_name <- paste("val.JAguirre_predictions", feature_name, ml_algorithm, paste(disgenet_phenotype, ".txt", sep=""), sep="_")
        print(val_file_name)
        valres <- vm[vm["Name"] == val_file_name,][,2:7]
        if (nrow(valres)==0){
          valres <- rep(NA, 6)
        }
        file_name <- paste("cv", feature_name, ml_algorithm, paste(disgenet_phenotype, ".txt", sep=""), sep="_")
        feature_file <- paste(output_dir, feature_name, file_name, sep="/")
        print(feature_file)
        if (file.exists(feature_file)){
          feature_df <- read.csv(feature_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
          if (ncol(feature_df) == 5){
            feature_df$f1 <- apply(feature_df[,2:5], 1, function(y) calculate.f1.from.metrics(precision = y['precision'], sensitivity = y['sensitivity']))
            feature_df$mcc <- apply(feature_df[,2:5], 1, function(y) calculate.mcc.from.metrics(accuracy = y['accuracy'], precision = y['precision'], sensitivity = y['sensitivity'], specificity = y['specificity']))
          }
          feature_data <- feature_df[!feature_df$model=="pred.comb",]
          res <- apply(feature_df[,2:7], 2, mean)
          results[nrow(results)+1,] <- c(plot_name, res, valres)
        } else {
          valres <- rep(NA, 6)
          results[nrow(results)+1,] <- c(plot_name, res, valres)
        }
      } else {
        for (j in 1:nrow(phenotypes)){
          diseaseid <- phenotypes[j,c("diseaseid")]
          diseaseterm <- phenotypes[j,c("diseaseterm")]
          database <- strsplit(feature_name, "_")[[1]][1]
          type_analysis <- paste(database, type_genes, sep='.')
          disgenet2gene_df <- disgenet2gene[c("geneid", "diseaseid", "diseaseterm", type_analysis)]
          colnames(disgenet2gene_df) <- c("geneid", "diseaseid", "diseaseterm", "type_analysis")
          genes <- unique(disgenet2gene_df$geneid[disgenet2gene_df$type_analysis==1 & disgenet2gene_df$diseaseid==diseaseid])
          if (length(genes)>=10) {
            file_name <- paste("cv", feature_name, ml_algorithm, paste(diseaseid, ".txt", sep=""), sep="_")
            feature_file <- paste(output_dir, feature_name, file_name, sep="/")
            if (file.exists(feature_file)){
              feature_df <- read.csv(feature_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
              if (ncol(feature_df) == 5){
                feature_df$f1 <- apply(feature_df[,2:5], 1, function(y) calculate.f1.from.metrics(precision = y['precision'], sensitivity = y['sensitivity']))
                feature_df$mcc <- apply(feature_df[,2:5], 1, function(y) calculate.mcc.from.metrics(accuracy = y['accuracy'], precision = y['precision'], sensitivity = y['sensitivity'], specificity = y['specificity']))
              }
              feature_data <- feature_df[!feature_df$model=="pred.comb",]
              res <- apply(feature_data[,2:7], 2, mean)
              results_phenotypes[nrow(results_phenotypes)+1,] <- c(diseaseid, res)
            }
          }
        }
        res <- apply(results_phenotypes[,2:7], 2, as.numeric) # Convert into numeric the content of the matrix
        res_mean <- apply(res, 2, mean)
        results[nrow(results)+1,] <- c(plot_name, res_mean, valres)
      }
    } else {
      file_name <- paste("cv", feature_name, paste(ml_algorithm, ".txt", sep=""), sep="_")
      feature_file <- paste(output_dir, file_name, sep="/")
      cm = read.csv(feature_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      if (ncol(cm) == 5){
        cm$f1 <- apply(cm[,2:5], 1, function(y) calculate.f1.from.metrics(precision = y['precision'], sensitivity = y['sensitivity']))
        cm$mcc <- apply(cm[,2:5], 1, function(y) calculate.mcc.from.metrics(accuracy = y['accuracy'], precision = y['precision'], sensitivity = y['sensitivity'], specificity = y['specificity']))
      }
      res <- apply(cm[1:10,2:7], 2, mean)
      #print(file_name)
      #print(res)
      if (length(ml_algorithms) > 1){
        plot_name <- paste(formal_feature_name, " (", toupper(ml_algorithm), ")", sep="")
      } else {
        plot_name <- paste(formal_feature_name, sep="")
      }
      results[nrow(results)+1,] <- c(plot_name, res, valres)
    }
  }  
}


### Create a Heatmap with the results ###

# heat_df <- as.matrix(results[,2:13]) # Matrix needed to insert content into Heatmap
# colnames(heat_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score", "MCC", "Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score", "MCC")
# rownames(heat_df) <- results[,1]
# Cairo::CairoPDF(output_plot_pdf) # Save in PDF
# Heatmap(heat_df, name = "results", km = 0, 
#         col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
#         na_col = "grey", # NA color
#         cluster_rows=FALSE, cluster_columns=FALSE, 
#         row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
#         column_title_gp = gpar(fontsize = 12),
#         column_split = c(rep("Testing set", ncol(heat_df)/2), rep("Validation set", ncol(heat_df)/2)),
#         #column_title = c("Testing set", "Independent set"),
#         column_gap = unit(5, "mm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 10))
#         }
# )
# dev.off()


# Cairo::CairoPNG(output_plot_png, dpi=300, width = 8, height = 8, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
# Heatmap(heat_df, name = "results", km = 0, 
#         col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
#         na_col = "grey", # NA color
#         cluster_rows=FALSE, cluster_columns=FALSE, 
#         row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
#         column_title_gp = gpar(fontsize = 12),
#         column_split = c(rep("Testing set", 6), rep("Validation set", 6)),
#         column_gap = unit(5, "mm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 10))
#         }
# )
# dev.off()


### Define output plots ###
if(length(ml_algorithms) == 1){
  output_plot_pdf <- paste(main_directory, sprintf("results/plots/heatmap_results_crossvalidation_with_F1_MCC_%s.pdf", ml_algorithms[1]), sep="/")
  output_plot_png <- paste(main_directory, sprintf("results/plots/heatmap_results_crossvalidation_with_F1_MCC_%s.png", ml_algorithms[1]), sep="/")
  plot_width <- 8
  plot_height <- 7
} else {
  ml_algorithm <- paste(ml_algorithms, collapse = "_")
  output_plot_pdf <- paste(main_directory, sprintf("results/plots/heatmap_results_crossvalidation_with_F1_MCC_%s.pdf", ml_algorithm), sep="/")
  output_plot_png <- paste(main_directory, sprintf("results/plots/heatmap_results_crossvalidation_with_F1_MCC_%s.png", ml_algorithm), sep="/")
  plot_width <- 9
  plot_height <- 9
}

# ### Define table for plot ###
# # Define the main heatmap
# heat_df <- as.matrix(results[,c(2:6,8:12)]) # Matrix needed to insert content into Heatmap
# colnames(heat_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score", "Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score")
# heat_df <- apply(heat_df, 2, as.numeric) # Convert into numeric the content of the matrix
# rownames(heat_df) <- results[,1]
# # Define the mcc heatmap
# mcc_df <- as.matrix(results[,c(7,13)]) # Matrix needed to insert content into Heatmap
# colnames(mcc_df) <- c("MCC testing", "MCC hold-out")
# mcc_df <- apply(mcc_df, 2, as.numeric) # Convert into numeric the content of the matrix
# rownames(mcc_df) <- results[,1]
# 
# ### Plot PDF ###
# Cairo::CairoPDF(output_plot_pdf, width = plot_width, height = plot_height) # Save in PDF
# h1 <- Heatmap(heat_df, name = "results", km = 0, 
#               col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
#               na_col = "grey", # NA color
#               cluster_rows=FALSE, cluster_columns=FALSE, 
#               row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
#               column_title_gp = gpar(fontsize = 12),
#               column_split = factor(c(rep("Testing set", ncol(heat_df)/2), rep("Independent hold-out test set", ncol(heat_df)/2)), levels = c("Testing set", "Independent hold-out test set")),
#               #column_title = c("Testing set", "Independent set"),
#               column_gap = unit(5, "mm"),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 10))
#               }
# )
# h2 <- Heatmap(mcc_df, name = "MCC", km = 0, 
#               #col = colorRamp2(c(min(mcc_df, na.rm = TRUE), 0, max(mcc_df, na.rm = TRUE)), c("yellow", "white", "#00b8ff")),
#               col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
#               cluster_rows=FALSE, cluster_columns=FALSE, 
#               row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
#               column_title_gp = gpar(fontsize = 10),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(sprintf("%.2f", mcc_df[i, j]), x, y, gp = gpar(fontsize = 10))
#               }
# )
# ht_list = h1 + h2
# draw(ht_list, ht_gap = unit(5, "mm"))
# dev.off()
# 
# ### Plot PNG ###
# Cairo::CairoPNG(output_plot_png, dpi=300, width = plot_width, height = plot_height, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
# h1 <- Heatmap(heat_df, name = "results", km = 0, 
#               col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
#               na_col = "grey", # NA color
#               cluster_rows=FALSE, cluster_columns=FALSE, 
#               row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
#               column_title_gp = gpar(fontsize = 12),
#               column_split = factor(c(rep("Testing set", ncol(heat_df)/2), rep("Independent hold-out test set", ncol(heat_df)/2)), levels = c("Testing set", "Independent hold-out test set")),
#               cluster_column_slices=FALSE,
#               #column_title = c("Testing set", "Independent set"),
#               column_gap = unit(5, "mm"),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 10))
#               }
# )
# h2 <- Heatmap(mcc_df, name = "MCC", km = 0, 
#               #col = colorRamp2(c(min(mcc_df, na.rm = TRUE), 0, max(mcc_df, na.rm = TRUE)), c("yellow", "white", "#00b8ff")),
#               col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
#               cluster_rows=FALSE, cluster_columns=FALSE, 
#               row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
#               column_title_gp = gpar(fontsize = 10),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(sprintf("%.2f", mcc_df[i, j]), x, y, gp = gpar(fontsize = 10))
#               }
# )
# ht_list = h1 + h2
# draw(ht_list, ht_gap = unit(5, "mm"), merge_legend = TRUE)
# dev.off()


### Plot PNG with each MCC separated ###
# Define the heatmap for testing set
heat_df <- as.matrix(results[,c(2:6,8:12)]) # Matrix needed to insert content into Heatmap
heat_test_df <- as.matrix(results[,c(2:6)]) # Matrix needed to insert content into Heatmap
colnames(heat_test_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score")
heat_test_df <- apply(heat_test_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(heat_test_df) <- results[,1]
# Define the mcc for testing set
mcc_df <- as.matrix(results[,c(7,13)]) # Matrix needed to insert content into Heatmap
mcc_test_df <- as.matrix(results[,c(7)]) # Matrix needed to insert content into Heatmap
colnames(mcc_test_df) <- c("MCC")
mcc_test_df <- apply(mcc_test_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(mcc_test_df) <- results[,1]
# Define the heatmap for independent set
heat_ind_df <- as.matrix(results[,c(8:12)]) # Matrix needed to insert content into Heatmap
colnames(heat_ind_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "F1-score")
heat_ind_df <- apply(heat_ind_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(heat_ind_df) <- results[,1]
# Define the mcc for independent set
mcc_ind_df <- as.matrix(results[,c(13)]) # Matrix needed to insert content into Heatmap
colnames(mcc_ind_df) <- c("MCC")
mcc_ind_df <- apply(mcc_ind_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(mcc_ind_df) <- results[,1]

Cairo::CairoPNG(output_plot_png, dpi=300, width = plot_width, height = plot_height, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
h1 <- Heatmap(heat_test_df, name = "results", km = 0, 
              col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
              na_col = "grey", # NA color
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 12),
              column_title = "          Testing set",
              cluster_column_slices=FALSE,
              show_heatmap_legend = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", heat_test_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
h2 <- Heatmap(mcc_test_df, name = "MCC", km = 0, 
              col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 10),
              show_heatmap_legend = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mcc_test_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
h3 <- Heatmap(heat_ind_df, name = "results2", km = 0, 
              col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
              na_col = "grey", # NA color
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 12),
              column_title = "          Independent hold-out test set",
              cluster_column_slices=FALSE,
              show_heatmap_legend = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", heat_ind_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
h4 <- Heatmap(mcc_ind_df, name = "MCC2", km = 0, 
              col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 10),
              show_heatmap_legend = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mcc_ind_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
ht_list = h1 + h2 + h3 + h4
draw(ht_list, ht_gap = unit(c(2,5,2,5), "mm"), merge_legend = TRUE)
dev.off()



Cairo::CairoPDF(output_plot_pdf, width = plot_width, height = plot_height) # Save in PDF
h1 <- Heatmap(heat_test_df, name = "results", km = 0, 
              col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
              na_col = "grey", # NA color
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 12),
              column_title = "          Testing set",
              cluster_column_slices=FALSE,
              show_heatmap_legend = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", heat_test_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
h2 <- Heatmap(mcc_test_df, name = "MCC", km = 0, 
              col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 10),
              show_heatmap_legend = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mcc_test_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
h3 <- Heatmap(heat_ind_df, name = "results2", km = 0, 
              col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
              na_col = "grey", # NA color
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 12),
              column_title = "          Independent hold-out test set",
              cluster_column_slices=FALSE,
              show_heatmap_legend = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", heat_ind_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
h4 <- Heatmap(mcc_ind_df, name = "MCC2", km = 0, 
              col = colorRamp2(c(-1, 0, 1), c("yellow", "white", "#00b8ff")),
              cluster_rows=FALSE, cluster_columns=FALSE, 
              row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
              column_title_gp = gpar(fontsize = 10),
              show_heatmap_legend = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mcc_ind_df[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)
ht_list = h1 + h2 + h3 + h4
draw(ht_list, ht_gap = unit(c(2,5,2,5), "mm"), merge_legend = TRUE)
dev.off()


