### Load packages ###
#install.packages('RxnSim', dependencies = TRUE)
library(RxnSim)


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
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
output_table <- paste(main_directory, "additional_data/tanimoto_smiles.tsv", sep="/")


### Load R object with compound information ###
load(dilirank_file)


### Calculate Tanimoto similarity matrix ###
tanimoto_matrix <- ms.compute.sim.matrix (drank.sel$SMILES, format = 'smiles', standardize = T, explicitH = F,
                                          sim.method = 'tanimoto', fp.type = 'extended', fp.mode = 'bit', fp.depth = 6,
                                          fp.size = 1024, clearCache = T)
tanimoto_matrix <- data.frame(tanimoto_matrix)
colnames(tanimoto_matrix) <- drugs
rownames(tanimoto_matrix) <- drugs
tanimoto_matrix$DILIConcern <- drank.sel$DILIConcern
tanimoto_matrix$vDILIConcern <- drank.sel$vDILIConcern
tanimoto_matrix$severity <- drank.sel$Severity.Class
write.table(tanimoto_matrix, file = output_table,row.names=TRUE, na='-',col.names=TRUE, sep="\t")



