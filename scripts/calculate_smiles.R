#install.packages('RxnSim', dependencies = TRUE)

library(RxnSim)

# Load R object with compound information
dilirank_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda"
load(dilirank_file)

tanimoto_matrix <- ms.compute.sim.matrix (drank.sel$SMILES, format = 'smiles', standardize = T, explicitH = F,
                                          sim.method = 'tanimoto', fp.type = 'extended', fp.mode = 'bit', fp.depth = 6,
                                          fp.size = 1024, clearCache = T)
tanimoto_matrix <- data.frame(tanimoto_matrix)
colnames(tanimoto_matrix) <- drugs
rownames(tanimoto_matrix) <- drugs
tanimoto_matrix$DILIConcern <- drank.sel$DILIConcern
tanimoto_matrix$vDILIConcern <- drank.sel$vDILIConcern
tanimoto_matrix$severity <- drank.sel$Severity.Class

output_table <- "/home/quim/PHD/Projects/camda/additional_data/tanimoto_smiles.tsv"
write.table(tanimoto_matrix, file = output_table,row.names=TRUE, na='-',col.names=TRUE, sep="\t")
