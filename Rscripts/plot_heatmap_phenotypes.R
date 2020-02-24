### Load packages ###
#install.packages('circlize', dependencies = TRUE)
#install.packages('bayesbio', dependencies = TRUE)
#BiocManager::install('ComplexHeatmap')
#BiocManager::install('jaccard')
library(Cairo)
library(circlize) # For Heatmap()
library(ComplexHeatmap) # For Heatmap()
library(bayesbio)

### Define variables ###
place = "home" #home or work

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
disease2gene_file <- paste(main_directory, "guildify_data/disease2gene_disgenet_guildify.tsv", sep="/")
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")
output_plot_pdf <- paste(main_directory, "results/plots/heatmap_phenotypes.pdf", sep="/")
output_plot_png <- paste(main_directory, "results/plots/heatmap_phenotypes.png", sep="/")
output_plot_inter_pdf <- paste(main_directory, "results/plots/heatmap_phenotypes_intersection.pdf", sep="/")
output_plot_inter_png <- paste(main_directory, "results/plots/heatmap_phenotypes_intersection.png", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(dilirank_file)

#### Obtain genes from phenotypes ####
disgenet2gene <- read.csv(disease2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
redundantphenotypes_df <- read.csv(redundantphenotypes_file, header=TRUE, sep="\t", stringsAsFactors = F)
type_genes = 'curated' # curated or all
database = 'disgenet'
type_analysis <- paste(database,type_genes,sep=".")
redundantphenotypes_df <- redundantphenotypes_df[redundantphenotypes_df$category == type_analysis,]
disgenet2gene <- disgenet2gene[!(disgenet2gene$diseaseid %in% redundantphenotypes_df$redundant.phenotype),]
disgenet2gene_df <- disgenet2gene[c("geneid", "diseaseid", "diseaseterm", type_analysis)]
colnames(disgenet2gene_df) <- c("geneid", "diseaseid", "diseaseterm", "type_analysis")
phenotypes <- unique(disgenet2gene[,c('diseaseid', 'diseaseterm')])
phenotypes <- phenotypes[order(phenotypes$diseaseterm),]


### Obtain the phenotypes with 10 genes or more ###
phenotypes_df <- c()
phenotypes_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(phenotypes_df) <- c("diseaseid","diseaseterm", "abbrev")

for (i in 1:nrow(phenotypes)){
  diseaseid <- phenotypes[i,c("diseaseid")]
  diseaseterm <- phenotypes[i,c("diseaseterm")]
  genes <- unique(disgenet2gene_df$geneid[disgenet2gene_df$type_analysis==1 & disgenet2gene_df$diseaseid==diseaseid])
  if (length(genes)>=10) {
    abbrev <- diseaseterm
    if (nchar(diseaseterm)>=20){
      abbrev <- paste(substr(diseaseterm, 1, 19), '.', sep = "")
    }
    phenotypes_df[nrow(phenotypes_df)+1,] <- c(diseaseid, diseaseterm, abbrev)
  }
}


### Calculate the Jaccard and nº common for each pair of phenotypes ###
jaccard_df <- data.frame(matrix(ncol = nrow(phenotypes_df), nrow = nrow(phenotypes_df)))
colnames(jaccard_df) <- phenotypes_df$abbrev
rownames(jaccard_df) <- phenotypes_df$abbrev
intersect_df <- data.frame(matrix(ncol = nrow(phenotypes_df), nrow = nrow(phenotypes_df)))
colnames(intersect_df) <- phenotypes_df$abbrev
rownames(intersect_df) <- phenotypes_df$abbrev

for (i in 1:nrow(jaccard_df)){
  phen1 <- rownames(jaccard_df)[i]
  phenid1 <- phenotypes_df$diseaseid[phenotypes_df$abbrev==phen1]
  genes1 <- unique(disgenet2gene_df$geneid[disgenet2gene_df$type_analysis==1 & disgenet2gene_df$diseaseid==phenid1])
  for (j in 1:ncol(jaccard_df)){
    phen2 <- colnames(jaccard_df)[j]
    phenid2 <- phenotypes_df$diseaseid[phenotypes_df$abbrev==phen2]
    genes2 <- unique(disgenet2gene_df$geneid[disgenet2gene_df$type_analysis==1 & disgenet2gene_df$diseaseid==phenid2])
    jacc <- jaccardSets(genes1, genes2)
    jaccard_df[i,j] <- jacc # [num_row, num_column]
    inter <- length(intersect(genes1, genes2))
    intersect_df[i,j] <- inter # [num_row, num_column]
  }
}


### Create a Heatmap with the results of Jaccard ###
heat_df <- as.matrix(jaccard_df) # Matrix needed to insert content into Heatmap

Cairo::CairoPDF(output_plot_pdf, width = 6, height = 5) # Save in PDF
Heatmap(heat_df, name = "Jaccard Index", km = 0, 
        col = colorRamp2(c(min(heat_df), 0.15, 1), c("white", "yellow", "red")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 8), column_names_gp =  gpar(fontsize = 8), column_names_rot = 80,
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)
dev.off()
Cairo::CairoPNG(output_plot_png, dpi=300, width = 6, height = 5, units = "in") # Save in PDF
Heatmap(heat_df, name = "Jaccard Index", km = 0, 
        col = colorRamp2(c(min(heat_df), 0.15, 1), c("white", "yellow", "red")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 8), column_names_gp =  gpar(fontsize = 8), column_names_rot = 80,
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", heat_df[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)
dev.off()


### Create a Heatmap with the results of Intersection ###
heat_inter_df <- as.matrix(intersect_df) # Matrix needed to insert content into Heatmap

Cairo::CairoPDF(output_plot_inter_pdf, width = 6, height = 5) # Save in PDF
Heatmap(heat_inter_df, name = "Common genes", km = 0, 
        col = colorRamp2(c(0, 50, 300), c("white", "yellow", "red")), 
        na_col = "grey", # NA color 
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 8), column_names_gp =  gpar(fontsize = 8), column_names_rot = 80,
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%d", heat_inter_df[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)
dev.off()
Cairo::CairoPNG(output_plot_inter_png, dpi=300, width = 6, height = 5, units = "in") # Save in PDF
Heatmap(heat_inter_df, name = "Common genes", km = 0, 
        col = colorRamp2(c(0, 50, 300), c("white", "yellow", "red")), 
        na_col = "grey", # NA color
        cluster_rows=FALSE, cluster_columns=FALSE, 
        row_names_gp = gpar(fontsize = 8), column_names_gp =  gpar(fontsize = 8), column_names_rot = 80,
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%d", heat_inter_df[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)
dev.off()



