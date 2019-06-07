##CAMDA challenge
library(stringr)
library(FactoMineR)
library(factoextra)
library(dplyr)

#Load files
setwd("E:/Dropbox/TransQST2/CAMDA2019/")
load("CAMDA_l1000_1314compounds-dilirank.v2.rda")
camda.all <- read.delim("matrix_compounds_morethan3.tsv", sep = "\t", header = T,
                        row.names = 1, check.names = F, stringsAsFactors = F)
landmark.genes <- read.delim("landmark_genes.tsv", row.names = 1, header = T, sep = "\t")

#Process the file into a nice df
camda.tall <- data.frame(t(camda.all[,-c(1,2)]), check.names = F, stringsAsFactors = F)
camda.tall.split <- data.frame(str_split_fixed(row.names(camda.tall), "_", 5), 
                               stringsAsFactors = F)
colnames(camda.tall.split) <- c("compound", "cell_line", "dose", "time", "replicate")
##Add Dili class to the processed set
camda.tall.split$dili <- drank.sel$vDILIConcern[match(camda.tall.split$compound, drank.sel$Compound.Name)]
camda.proc <- cbind(camda.tall.split, camda.tall)

##Have a subset of only landmark genes
camda.landmark <- camda.proc[,colnames(camda.proc) %in% row.names(landmark.genes)]
camda.landmark <- cbind(camda.proc[,1:6], camda.landmark)


##PCA analyses##########
##Select dataset 
dataset <- camda.landmark

##Selection based on DILI risk
#dataset <- droplevels(dataset[dataset$dili %in% c("vMost-DILI-Concern", "vNo-DILI-Concern"),])
  
#Run exploratory PCA  
res.pca <- PCA(dataset[,7:ncol(dataset)], scale.unit = T, graph = T)
var <- get_pca_var(res.pca)
fviz_pca_ind(res.pca, addEllipses = T, col.ind = dataset$dili)
fviz_pca_ind(res.pca, addEllipses = T, col.ind = dataset$cell_line)
fviz_pca_ind(res.pca, addEllipses = T, col.ind = dataset$dose)
fviz_pca_ind(res.pca, addEllipses = T, col.ind = dataset$time)

##Kmeans with the dataset]
cellline <- "PHH"
dataset1 <- dataset[dataset$cell_line == cellline,]
dataset <- dataset1[,7:ncol(dataset1)]
#k2 <- kmeans(dataset[,7:ncol(dataset)], 2)
aggregate(dataset[,7:ncol(dataset)], by = list(k2$cluster), FUN = mean)
mydata <- data.frame(dataset[,1:7], k2$cluster)
summarise(group_by(mydata, dili, fit.cluster), count = n())

fviz_cluster(fit, data = dataset[,7:ncol(dataset)], geom = "point")
dataset %>% as_tibble() %>% mutate(cluster = fit$cluster,
                                   compound = row.names(dataset[,7:ncol(dataset)])) %>%

##Compare several clusters
#k2 <- kmeans(dataset, 2)
k3 <- kmeans(dataset, 3)
k4 <- kmeans(dataset, 4)
k5 <- kmeans(dataset, 5)

#p2 <- fviz_cluster(k2, data = dataset, geom = "point") + ggtitle("K = 2")
p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
p4 <- fviz_cluster(k4, data = dataset, geom = "point") + ggtitle("K = 4")
p5 <- fviz_cluster(k5, data = dataset, geom = "point") + ggtitle("K = 5")
library(gridExtra)
grid.arrange(p3, p4, p5, nrow =2)

##Determine the optimal number of clusters
fviz_nbclust(dataset, kmeans, method = "wss", k.max = 48)

##Summarize based on DILI classes
mydata <- data.frame(dataset1[,1:7], k3$cluster)
summarise(group_by(mydata, dili, fit.cluster), count = n())

##Plot the results
##Coordinates from the set 
cellline <- "HCC515"
dataset1 <- camda.landmark
#dataset1 <- dataset[dataset$cell_line == cellline,]
dataset <- dataset1[,7:ncol(dataset1)]

k3 <- kmeans(dataset, 3)
p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
data.plot <- p3
coord_clusters <- cbind(data.plot$data, dataset1[,1:6])

ggplot(data = coord_clusters, aes(x, y)) + 
  geom_point(aes(color=dili, shape = cluster, size =10)) + ggtitle(cellline) #+
  #geom_label_repel(label = coord_clusters$compound)
  #geom_text(label = coord_clusters$compound)

summarise(group_by(coord_clusters, compound, cluster), count = n())


##Removing daunorubicin and vorinostat
##Plot the results
##Coordinates from the set 

compounds <- c("daunorubicin", "vorinostat")
dataset2 <- camda.landmark[!camda.landmark$compound %in% compounds,]
celllines <- unique(camda.landmark$cell_line)
plots <- list()
for(i in 1:length(celllines)){
dataset1 <- camda.landmark[camda.landmark$cell_line == celllines[[i]],]
#dataset1 <- dataset2[dataset2$cell_line == celllines[[i]],]
dataset <- dataset1[,7:ncol(dataset1)]
k3 <- kmeans(dataset, 3)
p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
data.plot <- p3
coord_clusters <- cbind(data.plot$data, dataset1[,1:6])
plots[[celllines[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
  geom_point(aes(shape = dili, color = time)) + #geom_text(label = coord_clusters$compound) +
    ggtitle(celllines[[i]])}#+
#geom_label_repel(label = coord_clusters$compound)
#geom_text(label = coord_clusters$compound)

pdf("plots_cell_lines_all_time.pdf", onefile = T)
for(i in 1:length(plots)){
  print(plots[[i]])
}
dev.off()

k3 <- kmeans(dataset, 3)
p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
k4 <- kmeans(dataset, 4)
p4 <- fviz_cluster(k4, data = dataset, geom = "point") + ggtitle("K = 4")
k5 <- kmeans(dataset, 4)
p5 <- fviz_cluster(k5, data = dataset, geom = "point") + ggtitle("K = 5")
library(gridExtra)
grid.arrange(p3, p4, p5, nrow =2)
fviz_nbclust(dataset, kmeans, method = "wss", k.max = 100)


summarise(group_by(coord_clusters, compound, cluster), count = n())



###Add more info to the sets
add1 <- read.delim("additional_info/GSE92742_Broad_LINCS_cell_info.txt", header = T, check.names = F, sep = "\t", stringsAsFactors = F)
add1.subset <- add1[,c("cell_id", "cell_type", "sample_type", "primary_site")]
add2 <- read.delim("additional_info/GSE92742_Broad_LINCS_gene_info.txt", header = T, check.names = F, sep = "\t")
add3 <- read.delim("additional_info/GSE92742_Broad_LINCS_pert_info.txt", header = T, check.names = F, sep = "\t")
add4 <- read.delim("additional_info/GSE92742_Broad_LINCS_sig_info.txt", header = T, check.names = F, sep = "\t")
hubs <- read.delim("additional_info/PHH_hubgenes.csv", header = T,
                   stringsAsFactors = F, check.names = F, sep = ";")

##Using the sets and selecting by gene hubs
compounds <- c("daunorubicin", "vorinostat")
dataset2 <- camda.landmark[!camda.landmark$compound %in% compounds,]
dataset3 <- cbind(merge(dataset2[,1:6], add1.subset, by.x = "cell_line", by.y = "cell_id"))
dataset4 <- cbind(dataset3, dataset2[,7:ncol(dataset2)])
##Only 24h
dataset4 <- dataset4[dataset4$time =="24h",]

##Separate per organ type
features <- unique(dataset4$primary_site)
plots <- list()
clusters.plots <- list()
plots.time <- list()
plots.compound <- list()
plots.compound2 <- list()
for(i in 1:length(features)){
  dataset1 <- dataset4[dataset4$primary_site == features[[i]],]
  #dataset1 <- dataset2[dataset2$cell_line == celllines[[i]],]
  dataset <- dataset1[,10:ncol(dataset1)]
  k3 <- kmeans(dataset, 3)
  p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
  data.plot <- p3
  coord_clusters <- cbind(data.plot$data, dataset1[,1:9])
  clusters.plots[[i]] <- p3
  plots[[features[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = dili, color = sample_type)) + #geom_text(label = coord_clusters$compound) +
    ggtitle(features[[i]])
  plots.time[[features[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = time)) #+ geom_text(label = coord_clusters$compound) +
    #theme(legend.position = "none")
  plots.compound[[features[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = compound)) + #geom_text(label = coord_clusters$compound) +
    theme(legend.position = "none")
  plots.compound2[[features[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = compound)) 
  }

pdf("plots_features_all.pdf", onefile = T)
for(i in 1:length(plots)){
  grid.arrange(clusters.plots[[i]],plots[[i]],  plots.time[[i]], plots.compound[[i]], ncol = 2)
}
dev.off()

##Subsetting using hub genes
dataset5 <- cbind(dataset4[,1:9], dataset4[,hubs$entrezID %in% colnames(dataset4)])


plots.hubs <- list()
plots <- list()
clusters.plots <- list()
plots.time <- list()
plots.compound <- list()
plots.compound2 <- list()
for(i in 1:length(celllines)){
  dataset1 <- dataset5[dataset5$cell_line == celllines[[i]],]
  #dataset1 <- dataset2[dataset2$cell_line == celllines[[i]],]
  dataset <- dataset1[,10:ncol(dataset1)]
  k3 <- kmeans(dataset, 3)
  p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
  data.plot <- p3
  coord_clusters <- cbind(data.plot$data, dataset1[,1:9])
  clusters.plots.hubs[[i]] <- p3
  plots[[celllines[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = dili, color = sample_type)) + #geom_text(label = coord_clusters$compound) +
    ggtitle(celllines[[i]])
  plots.time[[celllines[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = time)) #+ geom_text(label = coord_clusters$compound) +
  #theme(legend.position = "none")
  plots.compound[[celllines[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = compound)) + #geom_text(label = coord_clusters$compound) +
    theme(legend.position = "none")
  plots.compound2[[celllines[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = compound)) 
  plots.hubs[[celllines[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = dili, color = compound)) + #geom_text(label = coord_clusters$compound) +
    ggtitle(celllines[[i]])}#+
#geom_label_repel(label = coord_clusters$compound)
#geom_text(label = coord_clusters$compound)

pdf("hub_genes_subset.pdf", onefile = T)
for(i in 1:length(plots.hubs)){
  grid.arrange(clusters.plots.hubs[[i]],plots[[i]],  plots.time[[i]], plots.compound[[i]], ncol = 2)
}
dev.off()

### Liver diseases DisGenet
diseases <- unique(disgenet$name)


for(i in 1:length(diseases)){
  genes.disease <- disgenet$geneid[disgenet$name == diseases[[i]]]
  dataset1 <- cbind(dataset4[,1:9], dataset4[,colnames(dataset4) %in% genes.disease])
  #dataset1 <- dataset2[dataset2$cell_line == celllines[[i]],]
  if(ncol(dataset1) > 10){
  dataset <- dataset1[,10:ncol(dataset1)]
  k3 <- kmeans(dataset, 3)
  p3 <- fviz_cluster(k3, data = dataset, geom = "point") + ggtitle("K = 3")
  data.plot <- p3
  coord_clusters <- cbind(data.plot$data, dataset1[,1:9])
  clusters.plots[[i]] <- p3
  plots[[diseases[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = dili, color = sample_type)) + #geom_text(label = coord_clusters$compound) +
    ggtitle(diseases[[i]])
  plots.compound[[diseases[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = time, colour = cell_line)) + #geom_text(label = coord_clusters$compound) +
    theme(legend.position = "none")
  plots.compound2[[diseases[[i]]]] <- ggplot(data = coord_clusters, aes(x, y)) + 
    geom_point(aes(shape = cell_line, colour = compound)) 
  } else {i <- i +1}
  }


pdf("disgenet_genes_subset.pdf", onefile = T)
for(i in 1:length(plots)){
  if(!is.null(clusters.plots[[i]]) && !is.null(plots[[i]]) && !is.null(plots.compound[[i]])){
  grid.arrange(clusters.plots[[i]],plots[[i]], plots.compound[[i]], ncol = 2)
} else {i <- i + 1}}
dev.off()

