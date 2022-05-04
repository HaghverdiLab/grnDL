setwd("/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_GRN_construction_DataDriven/grnDL/Scripts/")

source("functions.R")
source("plots.R")
library(igraph)
library(batchelor)
library(network)
library(ggraph)
library(STRINGdb)


args = commandArgs(trailingOnly=TRUE)
folder_train <- args[1]
folder_test <- args[2]
out <- args[3]
print(args)

score_threshold <- ifelse(length(args) < 6, as.numeric(100),  as.numeric(args[6]))
edges <- ifelse(length(args) < 7, 3000, as.numeric(args[7]))
species <- ifelse(length(args) < 8, 9606, as.numeric(args[8]))

hvgs <- 1000
min_cells <- 10
min_connections <- 3
name <- paste(basename(folder_train), paste(sep="_",args[4:7], collapse = "_"), sep="_")
if(!dir.exists(paste(sep="/", out, name))) dir.create(paste(sep="/",out, name))


print(paste("Starting with", name))
start_time <- Sys.time()
data_train_raw <- readData(paste(sep="/", folder_train, "data_train.txt"))
data_test_raw <-  readData(paste(sep="/", folder_test, "data_test.txt"))
meta_train <- getMeta(paste(sep="/", folder_train, "meta_train.txt"))
meta_test  <- getMeta(paste(sep="/", folder_test, "meta_test.txt"))


data_train_raw <- data_train_raw[rownames(data_train_raw) %in% rownames(data_test_raw),]


if ( score_threshold > 0){
    string_db <- STRINGdb$new(version="11.5", species=species, score_threshold=score_threshold)
    genes <- rownames(data_test_raw)
    genes_translation <- as.data.frame(genes)
    genes_translation <- string_db$map(my_data_frame = genes_translation , my_data_frame_id_col_names="genes")
    genes_translation <- genes_translation[!is.na(genes_translation$STRING_id),]
    data_train_raw <- data_train_raw[rownames(data_train_raw) %in% genes_translation$genes,]
}

data <- data_train_raw[!(rowSums(data_train_raw != 0) < min_cells),]
meta <- meta_train[meta_train$id %in% colnames(data),]
rownames(meta) <- meta$id


data_train <- select_hvg(data, hvgs)


# Cell wise (Result: Sum of each row = 1)
 data_train <- apply(data_train, 2, normalize_l1)
# Gene wise (Variance of each column is 1, mean = 0) (Colum should be Gene)
data_train <- t(scale(t(data_train), center = TRUE, scale = TRUE))

AM <- getPartialCor(data_train)
if(score_threshold > 0) AM <- getScores(AM, string_db)

print(paste("Filtering", edges, min_connections))
AM <- getGRN(AM,edges, min_connections)


end_time <- Sys.time()
runtime <- end_time - start_time

print("Save results")
printProperties(as.matrix(AM),out, args[4:7], runtime, basename(folder_train))

grns <- showResults(AM, data_train, meta, "GRN")
ggsave(paste( sep="/",out, name, "grns.png"), plot = grns, device = "png", dpi=350,
       width=178, height = 178, units = "mm")


folder <- paste( sep="/",out, name)
trainingData <- getDataFile(data_train_raw, AM, folder, meta_train, "train")

testData <- getDataFile(data_test_raw, AM, folder, meta_test,"test")

write.table(rownames(trainingData), paste(folder, "feature_names.txt", sep="/"),
            sep=",", quote=F,row.names = F, col.names = F)

write.table(AM, paste( sep="/",folder, "adj.txt"), sep=",", col.names=T, row.names=F, quote=F, append=F)
write.table(data_train, paste( sep="/",folder, "expression.txt"), sep=",", col.names=T, row.names=T, quote=F, append=F)


u = as.data.frame(umap::umap(t(data_train))$layout)
colnames(u) <- c("x", "y")
u$id <- rownames(u)
u <- merge(u, meta, by="id")
umap <- ggplot(u, aes(x,y, color=class_))+geom_point()+theme_minimal()
ggsave(paste(sep="/",out, name, "umap.png"), plot = umap, device = "png", dpi=350,
       width=178, height = 178, units = "mm")

