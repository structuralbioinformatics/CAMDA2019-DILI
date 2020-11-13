### Define variables ###
place = "home" #home or work

if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PhD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}

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
disease2gene_file <- paste(main_directory, "guildify_data/disease2gene_disgenet_guildify.tsv", sep="/")
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")
wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected.txt", sep="/")
tanimoto_file <- paste(main_directory, "/additional_data/tanimoto_smiles.tsv", sep="/")
targets_file <- paste(main_directory, "additional_data/targets/targets_dgidb_hitpick_sea.tsv", sep="/")
# Output files
output_smiles_file <- paste(main_directory, "outputs/files/drug_smiles.tsv", sep="/")
output_targets_file <- paste(main_directory, "outputs/files/drug_targets.tsv", sep="/")
output_dili_landmark_file <- paste(main_directory, "outputs/files/dili_landmark.tsv", sep="/")
output_disgenet_guildify_file <- paste(main_directory, "outputs/files/disgenet_guildify_genes.tsv", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)
all_drugs <- unique(c(drug.dataset$drugs, drug.dataset$independent_drugs))


#### Get landmark genes ####
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]


### Get genes associated to phenotypes from DisGeNET and write output file ###
disgenet2gene <- read.csv(disease2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
redundantphenotypes_df <- read.csv(redundantphenotypes_file, header=TRUE, sep="\t", stringsAsFactors = F)
types_analysis <- c('disgenet.curated', 'guildify.curated')
redundantphenotypes_df <- redundantphenotypes_df[redundantphenotypes_df$category %in% types_analysis,]
disgenet2gene <- disgenet2gene[!(disgenet2gene$diseaseid %in% redundantphenotypes_df$redundant.phenotype),]
disgenet2gene_df <- disgenet2gene[c("geneid", "diseaseid", "diseaseterm", types_analysis)]
disgenet2gene_filt_df <- disgenet2gene_df[(disgenet2gene_df$disgenet.curated == 1) | (disgenet2gene_df$guildify.curated == 1),]
gene_symbols_disgenet_df <- gene_info_df[c("pr_gene_id", "pr_gene_symbol")][gene_info_df$pr_gene_id %in% disgenet2gene_df$geneid,]
rownames(gene_symbols_disgenet_df) <- NULL
colnames(gene_symbols_disgenet_df) <- c("geneid", "genesymbol")
disgenet_df <- merge(x = disgenet2gene_filt_df, y = gene_symbols_disgenet_df, by = "geneid")
disgenet_df <- disgenet_df[c("diseaseterm", "diseaseid", "genesymbol", "geneid", "disgenet.curated", "guildify.curated")][order(disgenet_df$diseaseterm),]
colnames(disgenet_df) <- c("Disease term", "Disease ID (UMLS)", "Gene symbol", "Gene ID (Entrez)", "DisGeNET", "GUILDify")
rownames(disgenet_df) <- NULL
write.table(disgenet_df, file = output_disgenet_guildify_file, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)


### Create SMILES file ###
smiles_dilrank_df <- drug.dataset$dilirank_df[drug.dataset$dilirank_df$pert_iname %in% drug.dataset$drugs,][c("pert_iname", "SMILES", "vDILIConcern")]
smiles_ind_df <- drank.sel[drank.sel$pert_iname %in% drug.dataset$independent_drugs,][c("pert_iname", "SMILES", "vDILIConcern")]
smiles_df <- rbind(smiles_dilrank_df, smiles_ind_df)
colnames(smiles_df) <- c("Drug name", "SMILES", "DILIConcern")
write.table(smiles_df, file = output_smiles_file, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)


### Load wilcoxon test analysis ###
wilcox_df = read.csv(wilcox_file, header = FALSE, sep = "\t")
#selected_genes <- unique(wilcox_df$gene_id[wilcox_df$p.value<0.05 & wilcox_df$landmark==TRUE])
selected_genes <- unique(wilcox_df[,1])


### Subset GCT object by landmark genes ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_df <- subset.expression(gct, landmark_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_ind_df <- subset.expression(gct, landmark_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)


### Subset GCT object by wilcoxon genes ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_wilcox_df <- subset.expression(gct, selected_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_wilcox_ind_df <- subset.expression(gct, selected_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)


### Create DILI landmark file ###
### Get a vector of gene symbols ordered equal to the vector of gene IDs ###
genes_df <- data.frame(1:length(selected_genes), selected_genes)
colnames(genes_df)<- c("id", "pr_gene_id")
genes_df <- merge(genes_df, gene_info_df[,c("pr_gene_id", "pr_gene_symbol")], by="pr_gene_id")
genes_df <- genes_df[order(genes_df$id), ]
genes_df$id <- NULL

