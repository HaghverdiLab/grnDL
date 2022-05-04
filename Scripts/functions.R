 
readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ",", verbose = F))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}

getMeta <- function(metafile, presentTypes=NULL){
  meta <- read.csv(metafile)
  meta <- meta[order(meta$class_),]
  if(!is.null(presentTypes)){
    meta <- rbind(meta[meta$class_ %in% presentTypes,],
                  meta[!(meta$class_ %in% presentTypes),] )
  }
  return(meta)
}

printProperties <- function(AM, out, settings, runtime, data ){
   G <- igraph::graph_from_adjacency_matrix(AM, mode = "undirected", weighted = T)
   info <- c(data, settings, runtime, length(V(G)), gsize(G), edge_density(G),
             min(AM[AM > 0]), max(AM), mean(AM[AM > 0]), mean(degree(G)))
   info <- t(as.data.frame(info))
    print(info)
    write.table(info, file = paste(sep="/", out, "grn_info.csv"), append = T, col.names=F, row.names = F, quote = F, sep=",")
}

pearson_residuals <- function(counts, theta){
    counts_sum1 = rowSums(counts)
    counts_sum0 = colSums(counts)
    counts_sum  = sum(counts)

    #get residuals
    mu = (counts_sum1  %*% t(counts_sum0)) / counts_sum
    z = (counts - mu) / sqrt(mu + mu**2/theta)

    #clip to sqrt(n)
    n = ncol(counts)
    z[z >  sqrt(n)] = sqrt(n)
    z[z < -sqrt(n)] = -sqrt(n)
    return(as.matrix(z))
}

select_hvg <- function(data, hvgs){
  residuals = pearson_residuals(data,100)
  residual_var = rowVars(residuals)
  names(residual_var) <- rownames(residuals)
  residual_var <- names(residual_var[order(residual_var, decreasing = T)])
  data = data[head(residual_var, hvgs), ]
  return(data)
}


normalize_l1 <- function(x) {
    return (log(x+0.001) / (sum(log(x+0.001))))
   
}


getPartialCor <- function(data){
    cor <- ppcor::pcor(t(data),method = "spearman")$estimate
    cor <- abs(cor)
    cor[is.na(cor)] <- 0
    diag(cor)     <- 0
    cor[cor < 0] <- 0
    colnames(cor) <- rownames(data)
    rownames(cor) <- rownames(data)
    return(cor)
}



getScores <- function(adj, string_db){
    adj <- as.data.frame(adj)
    adj$id <- rownames(adj)
    predicted_interactions <- reshape2::melt(adj)
    predicted_interactions <- translateGenes(predicted_interactions, "variable", "to", string_db)
    predicted_interactions <- translateGenes(predicted_interactions, "id", "from", string_db)
    string_interactions <- string_db$get_interactions(unique(c(predicted_interactions$from, predicted_interactions$to)))
    
    interactions <- merge(predicted_interactions, string_interactions, by = c("to", "from"), all=F)     
    interactions$combined_score <- interactions$combined_score  / 1000
    
    interactions$score <- interactions$value * interactions$combined_score
    interactions <- interactions[, c("id", "variable", "score")]
    interactions_copy <- interactions
    colnames(interactions_copy) <- c("variable", "id", "score")
    interactions <- rbind(interactions, interactions_copy)
    AM <- reshape2::acast(interactions, list(names(interactions)[1], names(interactions)[2]),
                          value.var="score", fun.aggregate = max,fill = 0)

    return(AM)       
}

translateGenes <- function(data, col, name, string_db){
    data <- string_db$map(my_data_frame = data, my_data_frame_id_col_names=col)
    colnames(data)[colnames(data) == "STRING_id"] <- name
    return(data)
}


getGRN <- function(AM, edges, minCon){
    print("Filtering....")
    AM <- as.data.frame(filter_interactions(as.matrix(AM), edges, minCon)) 
    print(paste(min(AM), max(AM)))
    print("Symmetricize...")
    AM <- symmetricize(AM)
    print(paste(min(AM), max(AM)))
    print("Largest component..")
    AM <- getLargestComponent(AM)
    print(paste(min(AM), max(AM)))
    
    return(AM)
}


filter_interactions <- function(cor, edges=2000, minCon=0){
  
    data <- reshape2::melt(cor)
    data <- data[order(data$value, decreasing=T),]
    data <- head(data, edges*2)
    data_copy <- data
    colnames(data_copy) <- c("Var2", "Var1", "value")
    data <- rbind(data, data_copy)
    data$value[data$value < 0] <- 0
    cor <- reshape2::acast(data, list(names(data)[1], names(data)[2]), value.var="value", fun.aggregate=max)
    cor[is.na(cor)]<- 0
    cor[cor < 0 ]<- 0
    genes <- unique(colnames(cor))
    if (minCon > 0){  #(length(genes) > 500){
        keep <- unlist(lapply(genes, function(g) sum(cor[,g] > 0))) > 3           
        cor <- cor[keep, keep] 
    }                          
    return(cor)
}  
                       

getLargestComponent <- function(AM){
  G = igraph::graph_from_adjacency_matrix(as.matrix(AM), weighted = T)
  if (igraph::is.connected(G, "weak")) return(AM)
  else{
    components <- igraph::clusters(G)
    biggest_cluster_id <- which.max(components$csize)
    vert_ids <- igraph::V(G)[components$membership == biggest_cluster_id]
    largeComponent <- igraph::induced_subgraph(G, vert_ids)
    AM <- as.matrix(as_adjacency_matrix(largeComponent, sparse=F, attr = "weight"))
    return(AM)
  }
}
                          

symmetricize <- function(x){

    #retrieve the (row,col) indices for the lower-diagonal entries in the matrix
    ldi <- matrix(c(row(x)[lower.tri(x)], col(x)[lower.tri(x)]), ncol = 2)
     #already column-major, so we can leave this

    #retrieve the (row,col) indices for the upper-diagonal entries in the matrix
    udi <- matrix(c(row(x)[upper.tri(x)], col(x)[upper.tri(x)]), ncol = 2)
    #in order for these to be in the symmetrical order as ldi, we need to sort.
    udi <- udi[order(udi[,1], udi[,2]),]

    #extract the upper and lower diagonal elements in a way that's symmetrical in their indexces
    ud <- x[udi]
    ld <- x[ldi]
 
    #replace with either the min, max, or mean, depending on the selection
    x[ldi] <- apply(rbind(ld,ud),2,max);
    x[udi] <- apply(rbind(ld,ud),2,max);

    
    return(x) 
  }
                              
getDataFile <- function(data, grn, output, meta, set){
  meta <- meta[meta$id %in% colnames(data),]
  data <- data[colnames(grn), meta$id]
  write.table(data, getFileName(output, "data", set), sep=",", quote=F,
              row.names = F, col.names = T)
  write.table(meta, getFileName(output, "classes", set), sep=",", quote=F,
              row.names = F, col.names = T)
  return(data)
}
                              
getFileName <- function(folder, type, set){
  name <- paste(folder, paste(type, set, sep="_"), sep="/")
  return(paste(name, "txt", sep="."))
}
 