### Load packages ###
library(caret)
library(mltools)

### Define variables ###
place = "home" #home or work

if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}

### Define directories and files ###
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
validation_dir <- paste(main_directory, "results/validations/3rd_round", sep="/")
output_table <- paste(main_directory, "results/validations/3rd_round_validations.tsv", sep="/")
file_names <- c(
  "val.JAguirre_predictions_mix_smiles_continuous_gbm.txt",
  "val.JAguirre_predictions_mix_smiles_continuous_rf.txt",
  "val.JAguirre_predictions_mix_smiles_discrete_gbm.txt",
  "val.JAguirre_predictions_mix_smiles_discrete_rf.txt",
  "val.JAguirre_predictions_mix_wilcox_continuous_gbm.txt",
  "val.JAguirre_predictions_mix_wilcox_continuous_rf.txt",
  "val.JAguirre_predictions_mix_wilcox_discrete_gbm.txt",
  "val.JAguirre_predictions_mix_wilcox_discrete_rf.txt",
  "val.JAguirre_predictions_mix_wilcoxsmiles_continuous_gbm.txt",
  "val.JAguirre_predictions_mix_wilcoxsmiles_continuous_rf.txt",
  "val.JAguirre_predictions_mix_wilcoxsmiles_discrete_gbm.txt",
  "val.JAguirre_predictions_mix_wilcoxsmiles_discrete_rf.txt",
  "val.GRIB_predictions_mix_disgenet_discrete_rf.txt",
  "val.GRIB_predictions_mix_disgenet_discrete_gbm.txt",
  "val.GRIB_predictions_mix_targets_discrete_rf.txt",
  "val.GRIB_predictions_mix_targets_discrete_gbm.txt",
  "val.GRIB_predictions_mix_wilcox_discrete_gbm.txt",
  "val.GRIB_predictions_mix_wilcox_discrete_rf.txt",
  "val.GRIB_predictions_mix_hub_discrete_gbm.txt",
  "val.GRIB_predictions_mix_hub_discrete_rf.txt",
  "val.GRIB_predictions_mix_wilcoxtargets_discrete_rf.txt",
  "val.GRIB_predictions_mix_wilcoxtargets_discrete_gbm.txt",
  "val.GRIB_predictions_mix_wilcoxtargetsmiles_discrete_rf.txt",
  "val.GRIB_predictions_mix_wilcoxtargetsmiles_discrete_gbm.txt"
)

file_names <- c("val.JCKevin_predictions_mulliner.csv",
                "val.JCKevin_predictions_dilirank.csv",
                "val.JAguirre_predictions_targets_rf.txt",
                "val.JAguirre_predictions_targets_gbm.txt",
                "val.JAguirre_predictions_smiles_rf.txt",
                "val.JAguirre_predictions_smiles_gbm.txt",
                "val.JAguirre_predictions_signature_smiles_rf.txt",
                "val.JAguirre_predictions_signature_smiles_gbm.txt",
                "val.JAguirre_predictions_signature_rf.txt",
                "val.JAguirre_predictions_signature_gbm.txt",
                "val.JAguirre_predictions_signature_glm.txt",
                "val.JAguirre_predictions_hub_smiles_gbm.txt",
                "val.JAguirre_predictions_hub_smiles_rf.txt",
                "val.JAguirre_predictions_disgenet_gbm_C0023892.txt",
                "val.JAguirre_predictions_disgenet_rf_C0023892.txt",
                "val.JAguirre_predictions_disgenet_smiles_gbm_C0023892.txt",
                "val.JAguirre_predictions_disgenet_smiles_rf_C0023892.txt",
                "val.JAguirre_predictions_guildify_gbm_C0023892.txt",
                "val.JAguirre_predictions_guildify_rf_C0023892.txt",
                "val.JAguirre_predictions_guildify_smiles_gbm_C0023892.txt",
                "val.JAguirre_predictions_guildify_smiles_rf_C0023892.txt",
                "val.JAguirre_predictions_landmark_gbm.txt",
                "val.JAguirre_predictions_landmark_rf.txt",
                "val.JAguirre_predictions_peng_gbm.txt",
                "val.JAguirre_predictions_peng_rf.txt"
)


### Load files ###
source(functions_file)


### Read validation results ###
sum <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(sum) <- c("Name", "Accuracy", "Precision", "Recall", "Specificity", "F1", "MCC")
for (file_name in file_names){
  file <- file.path(validation_dir, file_name)
  # Code from: https://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html
  cm = as.matrix(read.csv(file, header=TRUE, sep=",", stringsAsFactors=FALSE))
  # 0 = NO-DILI
  # 1 = DILI
  #  X0 X1 <== reference (cols)
  #0 18  6
  #1 14 11
  # <== predictions values (rows)
  #n = sum(cm) # number of instances
  #nc = nrow(cm) # number of classes
  #diag = diag(cm) # number of correctly classified instances per class 
  #rowsums = apply(cm, 1, sum) # number of predictions per class
  #colsums = apply(cm, 2, sum) # number of instances per class
  #p = colsums / n # distribution of instances over the actual classes
  #q = rowsums / n # distribution of instances over the predicted classes
  #accuracy = sum(diag) / n # Fraction of the instances that are correctly classified
  #precision = diag / rowsums # Fraction of correct predictions for certain class 
  #recall = diag / colsums # Fraction of instances of certain class correctly predicted
  #res <-  data.frame(precision, recall, f1) 
  TP <- cm[2,2]
  FP <- cm[2,1]
  TN <- cm[1,1]
  FN <- cm[1,2]
  #print(file_name)
  #print(cm)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- TP / (TP + FP)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  f1.from.rates <- calculate.f1.from.rates(TP=TP, FP=FP, FN=FN)
  f1.from.metrics = calculate.f1.from.metrics(precision = precision, sensitivity = sensitivity)
  mcc.from.function <- mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  mcc.from.rates <- calculate.mcc.from.rates(tp=TP, tn=TN, fp=FP, fn=FN)
  mcc.from.metrics <- calculate.mcc.from.metrics(accuracy = accuracy, precision = precision, sensitivity = sensitivity, specificity = specificity)
  print(paste("TP: ", TP))
  print(paste("TN: ", TN))
  print(paste("FP: ", FP))
  print(paste("FN: ", FN))
  print(paste("accuracy: ", accuracy))
  print(paste("precision: ", precision))
  print(paste("sensitivity: ", sensitivity))
  print(paste("specificity: ", specificity))
  print(paste("f1-score with tp/tn/fp/fn: ", f1.from.rates))
  print(paste("f1-score with a/p/n/s: ", f1.from.metrics))
  print(paste("mcc with r function: ", mcc.from.function))
  print(paste("mcc with tp/tn/fp/fn: ", mcc.from.rates))
  print(paste("mcc with a/p/n/s: ", mcc.from.metrics))
  sum[nrow(sum)+1,] <- c(file_name, accuracy, precision, sensitivity, specificity, f1.from.metrics, mcc.from.metrics)
  #print(res)
}

write.table(sum, file = output_table,row.names=FALSE, na="-",col.names=TRUE, sep="\t")