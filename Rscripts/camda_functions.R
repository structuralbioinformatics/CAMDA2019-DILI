#####################################################################
#####################################################################
############################# FUNCTIONS #############################
#####################################################################
#####################################################################


#####################################################################
######################### GENERAL FUNCTIONS #########################
#####################################################################


#####################################################################
######################## subset.drug.dataset ########################

#### Function to subset the drugs used in the analysis.
subset.drug.dataset<-function(drank.sel, outliers=c('daunorubicin', 'vorinostat'), remove.outliers=TRUE) {
  
  # Get independent dataset (Ambiguous drugs)
  independent_df <- drank.sel[drank.sel$vDILIConcern == "Ambiguous DILI-concern",]
  independent_drugs <- independent_df$pert_iname
  # Correct the 4 entries with lowercase concern
  drank.sel$DILIConcern[drank.sel$DILIConcern == "No-DILI-concern"] <- "No-DILI-Concern" 
  # We remove the first letter from vDILIConcern entries (except Ambiguous)
  drank.sel$vDILIConcern[drank.sel$vDILIConcern != "Ambiguous DILI-concern"] <- sub('.', '', drank.sel$vDILIConcern[drank.sel$vDILIConcern != "Ambiguous DILI-concern"])
  # We subset the table by those drugs that have the same value in vDILIConcern and DILIConcern columns
  dilirank_all_df <- drank.sel[drank.sel$vDILIConcern == drank.sel$DILIConcern,]
  # We remove the outlier drugs if necessary
  if (remove.outliers==TRUE){
    dilirank_df <- dilirank_all_df[!dilirank_all_df$pert_iname %in% outliers,]
  } else {
    dilirank_df <- dilirank_all_df
  }
  drugs <- dilirank_df$pert_iname
  most_concern_drugs <- dilirank_df$pert_iname[dilirank_df$DILIConcern == "Most-DILI-Concern"]
  less_concern_drugs <- dilirank_df$pert_iname[dilirank_df$DILIConcern == "Less-DILI-Concern"]
  #dili_concern_drugs <- dilirank_df$pert_iname[dilirank_df$DILIConcern == "Most-DILI-Concern" | dilirank_df$DILIConcern == "Less-DILI-Concern" ]
  no_concern_drugs <- dilirank_df$pert_iname[dilirank_df$DILIConcern == "No-DILI-Concern"]
  
  return(list(dilirank_df=dilirank_df, independent_df=independent_df, drugs=drugs, most_concern_drugs=most_concern_drugs, less_concern_drugs=less_concern_drugs, no_concern_drugs=no_concern_drugs, independent_drugs=independent_drugs));
}
#####################################################################
#####################################################################


#####################################################################
######################### subset.expression #########################

#### Function to subset the expression of a set of genes from samples
#### of a specific cell ID, dose and time, and merging the repeated
#### samples.
subset.expression<-function(gct, selected_genes, drugs, drugs_df, cell_id="PHH", pert_idose="10 ??M", pert_itime="24 h", merge_samples=FALSE) {
  
  # Subset the gct object by cell line, dose, time and genes
  id_selected_samples <- which(gct@cdesc$cell_id == cell_id & gct@cdesc$pert_idose == pert_idose & gct@cdesc$pert_itime == pert_itime)
  id_selected_genes <- which(gct@rid %in% selected_genes)
  gct_subset <- subset.gct(gct, cid=id_selected_samples, rid=id_selected_genes)
  
  # Subset by drugs of interest
  id_drugs <- which(gct_subset@cdesc$pert_iname %in% drugs)
  gct_drugs <- subset.gct(gct_subset, cid=id_drugs)
  
  # Create table of expression
  expression_df <- data.frame(t(gct_drugs@mat))
  
  # Get info of DILIrank
  mapping_drugs <- data.frame(gct_drugs@cid, gct_drugs@cdesc$pert_id, gct_drugs@cdesc$pert_iname)
  colnames(mapping_drugs) <- c("cid", "pert_id", "pert_iname")
  mapping_drugs <- merge(x = mapping_drugs, y = drugs_df[c("pert_iname", "Severity.Class", "DILIConcern")], by = "pert_iname")

  # Get severity values and DILIrank categories and map them to the table of expression
  severity <- c()
  dilirank <- c()
  pert_inames <- c()
  for (cid in rownames(expression_df)){
    row <- mapping_drugs[mapping_drugs$cid==cid,]
    severity[length(severity)+1] <- row$Severity.Class
    dilirank[length(dilirank)+1] <- row$DILIConcern
    pert_inames[length(pert_inames)+1] <- as.character(row$pert_iname)
  }
  expression_df$pert_iname <- pert_inames
  expression_df$severity <- severity
  expression_df$dilirank <- dilirank
  if (merge_samples==TRUE){
    #expression_df <- aggregate(.~pert_iname+severity+dilirank,expression_df,median) # Join samples from same drug
    expression_df <- aggregate(.~pert_iname+severity+dilirank,expression_df,function(val){y <- val[which.max( abs(val) )]; return(y)}) # Join samples from same drug
  }
  return(expression_df);
}
#####################################################################
#####################################################################


#####################################################################
##################### prepare.balanced.datasets #####################

#### Prepare balanced training and testing datasets having the same
#### number of DILI and no-DILI drugs, and the same proportion of
#### Most and Less DILI drugs.
prepare.balanced.datasets<-function(expression_data, num_datasets, most_concern_drugs, less_concern_drugs, no_concern_drugs, type_analysis="discrete", fraction_train=0.7) {

  fraction_test <- 1 - fraction_train
  
  # Make sure that the drugs considered are the ones that are in expression_data
  most_concern_drugs <- most_concern_drugs[most_concern_drugs %in% expression_data$pert_iname]
  less_concern_drugs <- less_concern_drugs[less_concern_drugs %in% expression_data$pert_iname]
  no_concern_drugs <- no_concern_drugs[no_concern_drugs %in% expression_data$pert_iname]
  
  # Calculate the number of less/most drugs taking into account the
  # proportions of less/most drugs and that the sum must be equal to
  # the number of no-concern drugs
  most_less_concern_drugs <- union(most_concern_drugs, less_concern_drugs)
  min_number_drugs <- min(length(most_less_concern_drugs), length(no_concern_drugs))
  proportion_most = length(most_concern_drugs) / length(most_less_concern_drugs)
  proportion_less = length(less_concern_drugs) / length(most_less_concern_drugs)
  num_most = round(proportion_most * min_number_drugs)
  num_less = round(proportion_less * min_number_drugs)
  
  # Get the drugs for the testing dataset
  drugs_testing <- c(sample(less_concern_drugs, round(num_less*fraction_test), replace=F), sample(most_concern_drugs, round(num_most*fraction_test), replace=F), sample(no_concern_drugs, round(min_number_drugs*fraction_test), replace=F))
  testing <- expression_data[expression_data$pert_iname %in% drugs_testing,]
  testing_with_names <- expression_data[expression_data$pert_iname %in% drugs_testing,]
  testing$pert_iname <- NULL
  if (type_analysis == "discrete"){
    testing$severity <- NULL
    colnames(testing)[match("dilirank", colnames(testing))] <- "cat"
    testing$cat[testing$cat == "No-DILI-Concern"] <- "vNo-DILI-Concern" # Rename No records
    testing$cat[testing$cat == "Less-DILI-Concern" | testing$cat == "Most-DILI-Concern"] <- "vMost-DILI-Concern / vLess-DILI-Concern" # Rename Most/Less records
  } else {
    testing$dilirank <- NULL
    colnames(testing)[match("severity", colnames(testing))] <- "cat"
  }
  
  # Get the drugs for the training datasets
  training_datasets <- list()
  training_datasets_with_names <- list()
  for (i in 1:num_datasets){
    # Get drugs
    drugs_training <- c(sample(less_concern_drugs[!less_concern_drugs%in%drugs_testing], round(num_less*fraction_train), replace=F), sample(most_concern_drugs[!most_concern_drugs%in%drugs_testing], round(num_most*fraction_train), replace=F), sample(no_concern_drugs[!no_concern_drugs%in%drugs_testing], round(min_number_drugs*fraction_train), replace=F))
    # Get gene expression
    training <- expression_data[expression_data$pert_iname %in% drugs_training,]
    training_with_names <- expression_data[expression_data$pert_iname %in% drugs_training,]
    # Prepare table
    training$pert_iname <- NULL
    if (type_analysis == "discrete"){
      training$severity <- NULL
      colnames(training)[match("dilirank", colnames(training))] <- "cat"
      training$cat[training$cat == "No-DILI-Concern"] <- "vNo-DILI-Concern" # Rename No records
      training$cat[training$cat == "Less-DILI-Concern" | training$cat == "Most-DILI-Concern"] <- "vMost-DILI-Concern / vLess-DILI-Concern" # Rename Most/Less records
    } else{
      training$dilirank <- NULL
      colnames(training)[match("severity", colnames(training))] <- "cat"
    }
    training_datasets[[length(training_datasets)+1]] <- training
    training_datasets_with_names[[length(training_datasets_with_names)+1]] <- training_with_names
  }
  return(list(testing=testing, training_datasets=training_datasets, testing_with_names=testing_with_names, training_datasets_with_names=training_datasets_with_names));
}
#####################################################################
#####################################################################


#####################################################################
############################ train.model ############################

#### Train a machine learning model using a cross-validation.
train.model<-function(training, testing, type_analysis="discrete", ml.method="rf", number.cv=10) {
  
  # Define parameters
  ctrl = trainControl(method = "cv", number = number.cv)
  
  # Train
  if (ml.method=="glm"){
    modFit = train(cat ~ ., data = training, method = "glm", family=binomial(), trControl = ctrl) 
  } else {
    modFit = train(cat ~ ., data = training, method = ml.method, trControl = ctrl, verbose=F) 
  }
  pred = predict(modFit, testing)
  if(type_analysis == "discrete") { 
    cor = NA
    a = confusionMatrix(table(pred, testing$cat))
    #print('CONFUSION MATRIX:')
    #print(a)
  } else {
    # Calculate correlation
    pred <- as.numeric(pred)
    testing$cat <- as.numeric(testing$cat)
    cor = cor(as.numeric(pred), as.numeric(testing$cat))
    # Assign category to testing and predictions
    pred_df<-data.frame(pred)
    pred_df$DILI <- "no"
    if (length(pred_df[pred_df$pred>=2,]$DILI)>0){
      pred_df[pred_df$pred>=2,]$DILI <- "less/most"
    }
    testing$DILI <- "no"
    if (length(testing$cat>=2)){
      testing[testing$cat>=2,]$DILI <- "less/most"
    }
    # Create confusion matrix based on category
    u <- union(pred_df$DILI, testing$DILI)
    t <- table(factor(pred_df$DILI, u), factor(testing$DILI, u))
    a = confusionMatrix(t)
  }
  return(list(modFit=modFit, a=a, cor=cor));
}
#####################################################################
#####################################################################


#####################################################################
##################### train.combine.classifiers #####################

#### Create n different models using a cross-validation and combine
#### all the models in a final model.
train.combine.classifiers<-function(datasets.list, ml.method="rf", type_analysis="discrete", number.cv=10, number.repetitions=10) {

  # Train different models n times by cross-validation
  train.list <- list()
  for (i in 1:number.repetitions){
    train.result <- train.model(training = datasets.list$training_datasets[[i]], testing = datasets.list$testing, type_analysis = type_analysis, ml.method = ml.method, number.cv = number.cv)
    train.list[[length(train.list)+1]] <- train.result
  }
  
  # Combine them
  pred.list <- list()
  pred.comb <- data.frame(matrix(ncol = 0, nrow = length(datasets.list$testing$cat)), stringsAsFactors = F)
  for (i in 1:length(train.list)){
    pred.result = predict(train.list[[i]]$modFit, datasets.list$testing)
    pred.list[[length(pred.list)+1]] <- pred.result
    Newcolname <- paste("pred.", i, sep = "") 
    pred.comb[[Newcolname]] <- pred.result
  }
  pred.comb[["cat"]] <- datasets.list$testing$cat

  if (type_analysis=="discrete"){
    modFit = train(cat ~ ., data=pred.comb, method = "rf")
    pred = predict(modFit)
    a = confusionMatrix(table(pred, datasets.list$testing$cat))
  } else{
    modFit = train(cat ~ ., data=pred.comb, method = "gam")
    pred = predict(modFit)
    pred_df<-data.frame(pred=pred, cat=datasets.list$testing$cat)
    pred_df$DILI <- "no"
    if (length(pred_df[pred_df$pred>=2,]$DILI)>0){
      pred_df[pred_df$pred>=2,]$DILI <- "less/most"
    }
    datasets.list$testing$DILI <- "no"
    if (length(datasets.list$testing$cat[datasets.list$testing$cat>=2])>0){
      datasets.list$testing[datasets.list$testing$cat>=2,]$DILI <- "less/most"
    }
    a = confusionMatrix(table(pred_df$DILI, datasets.list$testing$DILI))
  }
  return(list(modFit.comb=modFit, train.list=train.list, a=a));
}
#####################################################################
#####################################################################


#####################################################################
######################## validate.classifier ########################

#### Test the model on an independent dataset.
validate.classifier<-function(independent_data, output_file, modFit.comb, train.list, type_analysis="discrete") {
  
  # Redefine dataframe of independent data
  val.data <- independent_data
  val.data$pert_iname <- NULL
  val.data$severity <- NULL
  val.data$dilirank <- NULL

  # Make predictions using individual models
  pred.list <- list()
  pred.comb <- data.frame(matrix(ncol = 0, nrow = nrow(val.data)), stringsAsFactors = F)
  for (i in 1:length(train.list)){
    pred.result = predict(train.list[[i]]$modFit, val.data)
    pred.list[[length(pred.list)+1]] <- pred.result
    Newcolname <- paste("pred.", i, sep = "") 
    pred.comb[[Newcolname]] <- pred.result
  }
  # Make predictions using combined model
  pred = predict(modFit.comb, pred.comb)
  
  # Write predictions
  if (type_analysis=="discrete"){
    pred_df <- data.frame(Compound.Name=independent_data$pert_iname, Predicted.Label=pred)
    write.table(pred_df, file = output_file, row.names=FALSE, na="-", col.names=TRUE, sep=",", quote=FALSE)
  } else {
    pred_df <- data.frame(Compound.Name=independent_data$pert_iname, pred=pred)
    pred_df$Predicted.Label <- "vNo-DILI-Concern"
    if (length(pred_df[pred_df$pred>=2,]$Predicted.Label)>0){
      pred_df[pred_df$pred>=2,]$Predicted.Label <- "vMost-DILI-Concern / vLess-DILI-Concern"
    }
    pred_df$pred <- NULL
    write.table(pred_df, file = output_file, row.names=FALSE, na="-", col.names=TRUE, sep=",", quote=FALSE)
  }
  return(pred_df);
}
#####################################################################
#####################################################################


#####################################################################
################## write.machine.learning.results ###################

#### Write the machine learning results in a file.
write.machine.learning.results<-function(ml.results, output.file) {
  # Define table
  results.summary <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(results.summary) <- c("model", "accuracy", "precision", "sensitivity", "specificity")
  # Annotate single model results
  for (i in 1:length(ml.results$train.list)){
    result <- data.frame(ml.results$train.list[[i]]$a$byClass)
    model.name <- paste("pred.", i, sep = "") 
    summary <- c(model.name, result["Balanced Accuracy",], result["Precision",], result["Sensitivity",], result["Specificity",])
    results.summary[i,] <- summary
  }
  # Annotate combined model results
  result.comb <- data.frame(ml.results$a$byClass)
  summary.comb <- c("pred.comb", result.comb["Balanced Accuracy",], result.comb["Precision",], result.comb["Sensitivity",], result.comb["Specificity",])
  results.summary[nrow(results.summary)+1,] <- summary.comb
  # Write output file
  write.table(results.summary, file = output.file, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)
}
#####################################################################
#####################################################################


#####################################################################
####################### get.expression.matrix #######################

#### Get the expression matrix from a GCT object
get.expression.matrix<-function(gct_subset, selected_genes, dilirank_table, merge_samples=TRUE) {
  id_subset_sig <- which(gct_subset@rdesc$id %in% selected_genes)
  gct_subset_sig <- subset.gct(gct_subset, rid=id_subset_sig)
  # Get matrix in form of data frame
  sigmatrix <- data.frame(t(gct_subset_sig@mat), check.names = F, stringsAsFactors = F)
  sigmatrix$pert_iname <- gct_subset_sig@cdesc$pert_iname
  sigmatrix$cell_id <- gct_subset_sig@cdesc$cell_id
  sigmatrix$pert_idose <- gct_subset_sig@cdesc$pert_idose
  sigmatrix$pert_itime <- gct_subset_sig@cdesc$pert_itime
  sigmatrix <- sigmatrix[sigmatrix$pert_iname %in% dilirank_table$pert_iname,]
  sigmatrix$dili <- dilirank_table$vDILIConcern[match(sigmatrix$pert_iname, dilirank_table$pert_iname)]
  if(merge_samples==TRUE){
    sigmatrix <- aggregate(.~pert_iname+cell_id+pert_idose+pert_itime+dili,sigmatrix,function(val){y <- val[which.max( abs(val) )]; return(y)}) # Join replicates from same drug/cell/dose/time/dili
  }
  return(sigmatrix);
}
#####################################################################
#####################################################################


#####################################################################
########################### calculate.pca ###########################

#### Calculate the PCA of a signature matrix and plot the two main
#### components.
calculate.pca<-function(sigmatrix, analysis_name, merge_samples=TRUE) {
  
  # merge samples from same drug
  if(merge_samples==TRUE){
    sigmatrix <- sigmatrix[,!colnames(sigmatrix) %in% c("cell_id", "pert_idose", "pert_itime")]
    sigmatrix <- aggregate(.~pert_iname+dili,sigmatrix,function(val){y <- val[which.max( abs(val) )]; return(y)}) # Join samples from same drug
  }
  # get expression matrix
  info_fields <- c("pert_iname", "cell_id", "pert_idose", "pert_itime", "dili")
  expression <- sigmatrix[,!colnames(sigmatrix) %in% info_fields]
  info <- sigmatrix[,colnames(sigmatrix) %in% info_fields]
  # make pca
  pca <- prcomp(expression, scale=TRUE) # samples to be rows and genes to be columns
  # make a scree plot
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # computes the percentage of variation of each PC
  scree.plot <- barplot(pca.var.per, main=analysis_name, xlab="Principal Component", ylab="Percent Variation")
  # now make a fancy looking plot that shows the PCs and the variation:
  pca.data <- data.frame(Sample=rownames(data.frame(pca$x)),
                         X=pca$x[,1],
                         Y=pca$x[,2])
  pca.data <- cbind(pca.data, info)
  if(merge_samples==TRUE){
    pca.plot <- ggplot(data=pca.data, aes(x=X, y=Y, color=dili, shape = dili)) +
      geom_point() +
      xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
      ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
      theme_bw() +
      ggtitle(analysis_name)
  } else {
    pca.plot <- ggplot(data=pca.data, aes(x=X, y=Y, color=cell_id, shape = dili)) +
      geom_point() +
      xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
      ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
      theme_bw() +
      ggtitle(analysis_name)
  }
  return(list(pca.plot=pca.plot, scree.plot=scree.plot));
}
#####################################################################
#####################################################################


#####################################################################
#####################################################################
#####################################################################


#####################################################################
######################### SPECIFIC FUNCTIONS ########################
#####################################################################


#####################################################################
####################### prepare.balanced.drugs ######################

#### Prepare balanced sets of drugs having the same number of DILI
#### and no-DILI drugs, and the same proportion of Most and Less DILI
#### drugs.
prepare.balanced.drugs<-function(num_datasets, most_concern_drugs, less_concern_drugs, no_concern_drugs, fraction_train=0.7) {
  
  fraction_test <- 1 - fraction_train
  
  # Calculate the number of less/most drugs taking into account the
  # proportions of less/most drugs and that the sum must be equal to
  # the category with minimum number of drugs (DILI or NO-DILI)
  most_less_concern_drugs <- union(most_concern_drugs, less_concern_drugs)
  min_number_drugs <- min(length(most_less_concern_drugs), length(no_concern_drugs))
  proportion_most = length(most_concern_drugs) / length(most_less_concern_drugs)
  proportion_less = length(less_concern_drugs) / length(most_less_concern_drugs)
  num_most = round(proportion_most * min_number_drugs)
  num_less = round(proportion_less * min_number_drugs)
  
  # Get the drugs for the testing dataset
  drugs_testing <- c(sample(less_concern_drugs, round(num_less*fraction_test), replace=F), sample(most_concern_drugs, round(num_most*fraction_test), replace=F), sample(no_concern_drugs, round(min_number_drugs*fraction_test), replace=F))
  
  # Get the drugs for the training datasets
  drugs_training_list <- list()
  for (i in 1:num_datasets){
    # Get drugs
    drugs_training <- c(sample(less_concern_drugs[!less_concern_drugs%in%drugs_testing], round(num_less*fraction_train), replace=F), sample(most_concern_drugs[!most_concern_drugs%in%drugs_testing], round(num_most*fraction_train), replace=F), sample(no_concern_drugs[!no_concern_drugs%in%drugs_testing], round(min_number_drugs*fraction_train), replace=F))
    drugs_training_list[[length(drugs_training_list)+1]] <- drugs_training
  }
  return(list(drugs_testing=drugs_testing, drugs_training_list=drugs_training_list));
}
#####################################################################
#####################################################################


#####################################################################
################ prepare.datasets.from.balanced.drugs ###############

#### Prepare balanced sets of drugs having the same number of DILI
#### and no-DILI drugs, and the same proportion of Most and Less DILI
#### drugs.
prepare.datasets.from.balanced.drugs<-function(dataset, drugs_testing, drugs_training_list, type_analysis="discrete") {
  
  # Get the drugs for the testing dataset
  testing <- dataset[dataset$pert_iname %in% drugs_testing,]
  testing$pert_iname <- NULL
  if (type_analysis == "discrete"){
    testing$severity <- NULL
    colnames(testing)[match("dilirank", colnames(testing))] <- "cat"
    testing$cat[testing$cat == "No-DILI-Concern"] <- "vNo-DILI-Concern" # Rename No records
    testing$cat[testing$cat == "Less-DILI-Concern" | testing$cat == "Most-DILI-Concern"] <- "vMost-DILI-Concern / vLess-DILI-Concern" # Rename Most/Less records
  } else {
    testing$dilirank <- NULL
    colnames(testing)[match("severity", colnames(testing))] <- "cat"
  }
  
  # Get the features of drugs for the training datasets
  training_datasets <- list()
  for (i in 1:length(drugs_training_list)){
    # Get features
    drugs_training <- drugs_training_list[[i]]
    training <- dataset[dataset$pert_iname %in% drugs_training,]
    # Prepare data
    training$pert_iname <- NULL
    if (type_analysis == "discrete"){
      training$severity <- NULL
      colnames(training)[match("dilirank", colnames(training))] <- "cat"
      training$cat[training$cat == "No-DILI-Concern"] <- "vNo-DILI-Concern" # Rename No records
      training$cat[training$cat == "Less-DILI-Concern" | training$cat == "Most-DILI-Concern"] <- "vMost-DILI-Concern / vLess-DILI-Concern" # Rename Most/Less records
    } else{
      training$dilirank <- NULL
      colnames(training)[match("severity", colnames(training))] <- "cat"
    }
    training_datasets[[i]] <- training
  }
  return(list(testing=testing, training_datasets=training_datasets));
}
#####################################################################
#####################################################################


#####################################################################
####################### combine.disgenet.models #####################

#### Combine all the single models, and also the combined ones from
#### the DisGeNET analysis.
combine.disgenet.models<-function(disgenet.datasets.list, disgenet.results.list, output.file) {

  # Define table
  pred.single.models <- data.frame(matrix(ncol = 0, nrow = length(disgenet.datasets.list[[1]]$testing$cat)), stringsAsFactors = F)
  pred.combined.models <- data.frame(matrix(ncol = 0, nrow = length(disgenet.datasets.list[[1]]$testing$cat)), stringsAsFactors = F)
  # Loop over all phenotypes
  for (i in 1:length(disgenet.datasets.list)){
    # Get dataset and results of the phenotype
    datasets.list <- disgenet.datasets.list[[i]]
    phenotype.results <- disgenet.results.list[[i]]
    # Get prediction of the combined model
    pred.combined.result = predict(phenotype.results$modFit.comb)
    comb.model.name <- paste("pred", i, sep = ".") 
    pred.combined.models[[comb.model.name]] <- pred.combined.result
    # Loop over all sigle models of the phenotype
    for (j in 1:length(phenotype.results$train.list)){
      # Get prediction of single model of the phenotype
      pred.single.result = predict(phenotype.results$train.list[[j]]$modFit, datasets.list$testing)
      single.model.name <- paste(i, "pred", j, sep = ".") 
      pred.single.models[[single.model.name]] <- pred.single.result
    }
  }
  pred.combined.models[["cat"]] <- disgenet.datasets.list[[1]]$testing$cat
  pred.single.models[["cat"]] <- disgenet.datasets.list[[1]]$testing$cat
  # Train combining the combined models
  modFit.combined = train(cat ~ ., data=pred.combined.models, method = "rf")
  prediction.combined = predict(modFit.combined)
  a.combined = confusionMatrix(table(prediction.combined, disgenet.datasets.list[[1]]$testing$cat))
  # Train combining the single models
  modFit.single= train(cat ~ ., data=pred.single.models, method = "rf")
  prediction.single = predict(modFit.single)
  a.single = confusionMatrix(table(prediction.single, disgenet.datasets.list[[1]]$testing$cat))
  # Write results
  results.summary <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(results.summary) <- c("model", "accuracy", "precision", "sensitivity", "specificity")
  result.single <- data.frame(a.single$byClass)
  summary.single <- c("single models", result.single["Balanced Accuracy",], result.single["Precision",], result.single["Sensitivity",], result.single["Specificity",])
  result.combined <- data.frame(a.combined$byClass)
  summary.combined <- c("combined models", result.combined["Balanced Accuracy",], result.combined["Precision",], result.combined["Sensitivity",], result.combined["Specificity",])
  results.summary[1,] <- summary.single
  results.summary[2,] <- summary.combined
  write.table(results.summary, file = output.file, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)
  return(list(a.combined=a.combined, a.single=a.single, modFit.combined=modFit.combined, modFit.single=modFit.single));
}
#####################################################################
#####################################################################


#####################################################################
################# validate.disgenet.combined.models #################

#### Test the model on an independent dataset.
validate.disgenet.combined.models<-function(independent_data, output_file_single, output_file_comb, modFit.single, modFit.comb, disgenet.results.list) {
  
  # Make predictions using individual models
  pred.single <- data.frame(matrix(ncol = 0, nrow = nrow(independent_data[[1]])), stringsAsFactors = F) # The nrow is the number of drugs in the validation dataset
  pred.comb <- data.frame(matrix(ncol = 0, nrow = nrow(independent_data[[1]])), stringsAsFactors = F)
  for (i in 1:length(disgenet.results.list)){

    # Redefine dataframe of independent data
    val.data <- independent_data[[i]]
    val.data$pert_iname <- NULL
    val.data$severity <- NULL
    val.data$dilirank <- NULL
    
    # Get dataset and results of the phenotype
    phenotype.results <- disgenet.results.list[[i]]
    
    pred.phenotype.comb <- data.frame(matrix(ncol = 0, nrow = nrow(val.data)), stringsAsFactors = F)
    for (j in 1:length(phenotype.results$train.list)){
      pred.result = predict(phenotype.results$train.list[[j]]$modFit, val.data)
      single.model.name <- paste(i, "pred", j, sep = ".")
      comb.model.name <- paste("pred", j, sep = ".")
      pred.single[[single.model.name]] <- pred.result
      pred.phenotype.comb[[comb.model.name]] <- pred.result
    }
    pred.comb.result = predict(phenotype.results$modFit.comb, pred.phenotype.comb)
    comb.final.model.name <- paste("pred", i, sep = ".")
    pred.comb[[comb.final.model.name]] <- pred.comb.result
    
  }
  # Make predictions using combined model
  prediction.single = predict(modFit.single, pred.single)
  prediction.comb = predict(modFit.comb, pred.comb)
  
  # Write predictions
  pred_single_df <- data.frame(Compound.Name=independent_data[[1]]$pert_iname, Predicted.Label=prediction.single)
  write.table(pred_single_df, file = output_file_single, row.names=FALSE, na="-", col.names=TRUE, sep=",", quote=FALSE)
  pred_comb_df <- data.frame(Compound.Name=independent_data[[1]]$pert_iname, Predicted.Label=prediction.comb)
  write.table(pred_comb_df, file = output_file_comb, row.names=FALSE, na="-", col.names=TRUE, sep=",", quote=FALSE)
  return(list(pred_single_df=pred_single_df, pred_comb_df=pred_comb_df));
}
#####################################################################
#####################################################################


#####################################################################
################### subset.expression.by.samples ####################

#### Function to subset the expression of a set of genes from 
#### specific sample IDs.
subset.expression.by.samples<-function(gct, selected_genes, selected_samples, drugs, drugs_df, merge_samples=TRUE, merge_criteria="median") {
  
  # Subset the gct object by samples
  id_selected_samples <- which(gct@cid %in% selected_samples)
  id_selected_genes <- which(gct@rid %in% selected_genes)
  gct_subset <- subset.gct(gct, cid=id_selected_samples, rid=id_selected_genes)

  # Subset by drugs of interest
  id_drugs <- which(gct_subset@cdesc$pert_iname %in% drugs)
  gct_drugs <- subset.gct(gct_subset, cid=id_drugs)
  
  # Create table of expression
  expression_df <- data.frame(t(gct_drugs@mat))
  
  # Get info of DILIrank
  mapping_drugs <- data.frame(gct_drugs@cid, gct_drugs@cdesc$pert_id, gct_drugs@cdesc$pert_iname)
  colnames(mapping_drugs) <- c("cid", "pert_id", "pert_iname")
  mapping_drugs <- merge(x = mapping_drugs, y = drugs_df[c("pert_iname", "Severity.Class", "DILIConcern")], by = "pert_iname")
  
  # Get severity values and DILIrank categories and map them to the table of expression
  severity <- c()
  dilirank <- c()
  pert_inames <- c()
  for (cid in rownames(expression_df)){
    row <- mapping_drugs[mapping_drugs$cid==cid,]
    severity[length(severity)+1] <- row$Severity.Class
    dilirank[length(dilirank)+1] <- row$DILIConcern
    pert_inames[length(pert_inames)+1] <- as.character(row$pert_iname)
  }
  expression_df$pert_iname <- pert_inames
  expression_df$severity <- severity
  expression_df$dilirank <- dilirank
  if (merge_samples==TRUE){
    if (merge_criteria=="median"){
      expression_df <- aggregate(.~pert_iname+severity+dilirank,expression_df,median) # Join samples from same drug by median
    } else if (merge_criteria=="mean"){
      expression_df <- aggregate(.~pert_iname+severity+dilirank,expression_df,mean) # Join samples from same drug by mean
    } else {
      expression_df <- aggregate(.~pert_iname+severity+dilirank,expression_df,function(val){y <- val[which.max( abs(val) )]; return(y)}) # Join samples from same drug by highest value
    }
  }
  return(expression_df);
}
#####################################################################
#####################################################################


#####################################################################
#################### calculate.mcc.from.metrics #####################
#####################################################################
#####################################################################

#### Function to calculate MCC from accuracy/precision/sensitivity/
#### specificity
calculate.mcc.from.metrics<-function(accuracy, precision, sensitivity, specificity) {
  A <- accuracy
  P <- precision
  N <- sensitivity
  S <- specificity
  mcc.1st <- (P/(1-P)) * (S/(1-S)) - (P/(A*(1-P))) - (S/(A*(1-S))) + 1 + (P/(1-P)) + (S/(1-S))
  mcc.2nd <- ((P/(1-P)) + 1) * ( (P/(1-P)) + (P/(A*(1-P))) + (S/(A*(1-S))) - 1 - (P/(1-P)) - (S/(1-S)) ) * ((S/(1-S)) + 1) * ( (S/(1-S)) + (P/(A*(1-P))) + (S/(A*(1-S))) - 1 - (P/(1-P)) - (S/(1-S)) )
  mcc <- mcc.1st / sqrt(mcc.2nd)
  return(mcc);
}
#####################################################################
#####################################################################


#####################################################################
##################### calculate.mcc.from.rates ######################
#####################################################################
#####################################################################

#### Function to calculate MCC from TP/TN/FP/FN rates
calculate.mcc.from.rates<-function(tp, tn, fp, fn) {
  mcc <- ( tp*tn - fp*fn ) / sqrt( (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) )
  return(mcc);
}
#####################################################################
#####################################################################


#####################################################################
#####################################################################
#####################################################################
