rm(list = ls())
library(guildifyR)

# list of curated phenotypes
list_of_diseases<- c("C0023892", "C1262760,C0860207,C3658290", "C0220994", "C0023890,C0239946", "C0023891", "C0086565,C0023895", "C0162557", "C3241937,C0400966", "C2711227,C0015695")


# call to the disgenet API to retrieve genes associated to the list of liver phenotypes in disgenet CURATED
list_of_genes <- list()
for (l in list_of_diseases) {
  url <- paste0("https://www.disgenet.org/api/gda/disease/", l , "?source=CURATED&format=tsv")
  print(url)
  dataTsv <- RCurl::getURLContent( url     )
  myTextConnection <- textConnection( dataTsv )
  result <- read.csv( myTextConnection, header = TRUE, sep = "\t", colClasses=c("character"))
  close(myTextConnection)
  list_of_genes[[l]] <-  unique(result$gene_symbol)
} 



# call to the guildify server to expand genes associated to the list of liver phenotypes using PPI data
species="9606"
tissue="all"
network.source="BIANA"

scoring.options = list(netscore=T, repetitionSelector=3, iterationSelector=2) # NetScore
list_of_job_ids <- list()

for (l in list_of_diseases)   {
  keywords = paste0(list_of_genes[[l]], collapse = ";")
  print(keywords)
  result.table = query(keywords, species, tissue, network.source)
  jobid <- submit.job(result.table, species, tissue, network.source, scoring.options)
  list_of_job_ids[[as.character(l)]] <- jobid
}


# retrieve results
list_of_job_results <- list()
# names(list_of_job_results)
#for (l in setdiff(list_of_diseases, names(list_of_job_results)) ) {
for (l in list_of_diseases ) {
  result = retrieve.job(list_of_job_ids[[l]])
  list_of_job_results[[as.character(l)]] <- result@scores
}
