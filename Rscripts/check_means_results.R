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
  main_directory = "/Users/quim/Dropbox/UPF/PhD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/Projects/camda"
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
output_dir <- paste(main_directory, "results/crossvalidation", sep="/")
output_plot_pdf <- paste(main_directory, "results/plots/heatmap_results_crossvalidation_rf.pdf", sep="/")
output_plot_png <- paste(main_directory, "results/plots/heatmap_results_crossvalidation_rf.png", sep="/")
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


### Read the independent validation results ###
vm = read.csv(validation_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)


### Read cross-validation results and calculate the mean of all results ###
results <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(results) <- c("name","Accuracy", "Precision", "Sensitivity", "Specificity", "Accuracy", "Precision", "Sensitivity", "Specificity")

for (i in 1:length(feature_names)){
  feature_name <- feature_names[i]
  formal_feature_name <- plot_names[i]
  for (ml_algorithm in ml_algorithms){
    # Get the independent validation results
    val_file_name <- paste("val.JAguirre_predictions", feature_name, paste(ml_algorithm, ".txt", sep=""), sep="_")
    valres <- vm[vm["Name"] == val_file_name,][,2:5]
    if (nrow(valres)==0){
      valres <- c(NA, NA, NA, NA)
    }

    # Get the cross-validation results
    if (feature_name %in% c("disgenet", "disgenet_smiles", "guildify", "guildify_smiles")){
      results_phenotypes <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(results_phenotypes) <- c("name","accuracy", "precision", "sensitivity", "specificity")
      if (length(ml_algorithms) > 1){
        plot_name <- paste(formal_feature_name, " (", toupper(ml_algorithm), ")", sep="")
      } else {
        plot_name <- paste(formal_feature_name, sep="")
      }
      if (disgenet_phenotype != FALSE){
        val_file_name <- paste("val.JAguirre_predictions", feature_name, ml_algorithm, paste(disgenet_phenotype, ".txt", sep=""), sep="_")
        valres <- vm[vm["Name"] == val_file_name,][,2:5]
        if (nrow(valres)==0){
          valres <- c(NA, NA, NA, NA)
        }
        file_name <- paste("cv", feature_name, ml_algorithm, paste(disgenet_phenotype, ".txt", sep=""), sep="_")
        feature_file <- paste(output_dir, feature_name, file_name, sep="/")
        if (file.exists(feature_file)){
          feature_df <- read.csv(feature_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
          feature_data <- feature_df[!feature_df$model=="pred.comb",]
          res <- apply(feature_data[,2:5], 2, mean)
          results[nrow(results)+1,] <- c(plot_name, res, valres)
        } else {
          res <- c(NA, NA, NA, NA)
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
              feature_data <- feature_df[!feature_df$model=="pred.comb",]
              res <- apply(feature_data[,2:5], 2, mean)
              results_phenotypes[nrow(results_phenotypes)+1,] <- c(diseaseid, res)
            }
          }
        }
        res <- apply(results_phenotypes[,2:5], 2, as.numeric) # Convert into numeric the content of the matrix
        res_mean <- apply(res, 2, mean)
        results[nrow(results)+1,] <- c(plot_name, res_mean, valres)
      }
    } else {
      file_name <- paste("cv", feature_name, paste(ml_algorithm, ".txt", sep=""), sep="_")
      feature_file <- paste(output_dir, file_name, sep="/")
      cm = read.csv(feature_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      res <- apply(cm[1:10,2:5], 2, mean)
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
heat_df <- as.matrix(results[,2:9]) # Matrix needed to insert content into Heatmap
colnames(heat_df) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "Accuracy", "Precision", "Sensitivity", "Specificity")
heat_df <- apply(heat_df, 2, as.numeric) # Convert into numeric the content of the matrix
rownames(heat_df) <- results[,1]
Cairo::CairoPDF(output_plot_pdf) # Save in PDF
Heatmap(heat_df, name = "results", km = 0, 
        col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
        column_title_gp = gpar(fontsize = 12),
        column_split = c(rep("Testing set", 4), rep("Validation set", 4)),
        #column_title = c("Testing set", "Independent set"),
        column_gap = unit(5, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()
Cairo::CairoPNG(output_plot_png, dpi=300, width = 6, height = 6, units = "in") # Resolution taken from: https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
Heatmap(heat_df, name = "results", km = 0, 
        col = colorRamp2(c(min(heat_df, na.rm = TRUE), 0.5, max(heat_df, na.rm = TRUE)), c("red", "white", "green")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 10), column_names_gp =  gpar(fontsize = 10), column_names_rot = 60,
        column_title_gp = gpar(fontsize = 12),
        column_split = c(rep("Testing set", 4), rep("Validation set", 4)),
        column_gap = unit(5, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()
