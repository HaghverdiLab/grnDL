get_class_info <- function(expression, meta, celltype){
    
    e <- as.matrix(expression[, colnames(expression) %in%
                             meta$id[meta$class_ == celltype]])
    
    means <- rowMeans(e)
    variances <- as.numeric(matrixStats::rowVars(e)) #/ as.numeric(means)
    x <- data.frame(means,variances)
    x$variances[x$means == 0] <- 0
    x$celltype = celltype
    x$gene     = rownames(e)
    return(x)
}

get_GRN_graph <- function(adj){
    graph <- graph_from_adjacency_matrix(as.matrix(adj),mode = "undirected",  diag=F, weighted =TRUE )
    graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = igraph_opt("edge.attr.comb"))
    V(graph)$name <- colnames(adj)
    graph = igraph::delete.edges(graph, which(E(graph)$weight==0))
    return(graph)
}

visualize_GRN <- function(graph, data, color, celltype = NULL,
                          size=NULL, layout="nicely", rangeMean,
                          rangeVar, type="Expression"){
    
    if(!is.null(celltype)){
        x <- data[data$celltype==celltype,]
        x[,color] <- as.numeric(x[,color])
        rownames(x) <- x$gene
        x <- x[names(V(graph)),]
        o <- order(x[,color],na.last = FALSE)
        x <- x[order(x[,color],na.last = FALSE),]
        graph <- permute(graph, o)
    } else x <- data
    
    adj <- as_adjacency_matrix(graph)
    graph <- as.network(adj)
    graph %v% "value" = as.numeric(x[,color])
    graph %v% "var" = as.numeric(x[,size])
    plot <- ggraph(graph, layout =  layout) + coord_fixed() +
            geom_edge_link(color="lightgrey") +
            theme_light()+ ggtitle(celltype)  +
            theme(axis.text=element_blank(), axis.title=element_blank(),
                  legend.text=element_text(size=15), legend.title=element_text(size=15),
                  plot.title=element_text(size=20),  legend.position = "right")
    
    if(!is.null(size)){
        plot <- plot +  geom_node_point(aes(colour = value, size= var))+
                labs(size="Variance") + scale_size(range = rangeVar )
    } else{
         plot <- plot +  guides(size="none")+
                 geom_node_point(aes(colour = value), size= 3)

    }
    if(type=="Expression"){
        plot <- plot +labs(color="Mean") + scale_colour_gradient2(low = "darkgreen", mid = "white",
                                       high = "darkred", space = "Lab",
                                       na.value = "grey50", guide = "colourbar",
                                       aesthetics = "colour", midpoint=0)
    
    } else if (type == "Attention"){
        plot <- plot +  labs(color="Attention") +
                scale_colour_gradient2(low = "green", mid = "grey",
                                       high = "red", space = "Lab",
                                       na.value = "grey50", guide = "colourbar",
                                       aesthetics = "colour", midpoint=0)
    }
    #plot <- plot + guides(color= guide_colourbar(title.position = "top",
     #                           title.hjust = .5))
   return(plot)

}

summarize_data <- function(data, label=NULL, steps=NULL){
    if(!is.null(label)) data <- data[data$Version == label,]
    if(!is.null(label)) data <- data[data$CellsPerCelltype %in% steps,]
    
    data_melt <- reshape2::melt(data, value.name = "Accuracy",
                                id.vars = c("Method", "Version", "CellsPerCelltype", "SetNr"))
    
    colnames(data_melt)[colnames(data_melt) == "variable"] <- "CellType"
    data_summary <- data_melt %>%
                    group_by(Method, Version, CellsPerCelltype, CellType) %>%
                    summarise(mean = mean(Accuracy),
                              sd= sd(Accuracy),
                              q25= quantile(Accuracy, probs=0.25),
                              q75= quantile(Accuracy, probs=0.75),
                              median = median(Accuracy))
    
    return(as.data.frame(data_summary))   
}

showResults <- function(AM, dataTrain, meta_train, showPlot="Hist", celltypes=NULL){
    expression <- dataTrain[rownames(dataTrain) %in% colnames(AM),]
    print(paste("data", ncol(expression), nrow(expression)))
    if(nrow(AM)==0) return()
    
    graph <- get_GRN_graph(AM)
 
    if(is.null(celltypes))celltypes <- unique(meta_train$class_)

    sets <- lapply(celltypes, function (celltype) 
                   get_class_info(as.data.frame(expression), as.data.frame(meta_train), celltype))               
    set <- do.call(rbind, sets)  
    genes <- unique(set$gene)
    x <- unlist(lapply(genes, function(g) min(set$means[set$gene == g])))

    if(showPlot == "Hist" & nrow(AM) > 90 ){
        values <- reshape2::melt(AM)
        return(ggplot(values,aes(x=value))+ geom_histogram()+
               theme_minimal()+scale_x_continuous(breaks = seq(0,max(values$value), by = 0.01)))  }
    if(showPlot=="GRN" & nrow(AM) > 90 &  length(E(graph)) < 10000){
        #layout <- umap::umap(AM, n_components = 2, n_components = 10)$layout  
        layout <- get_grn_layout(graph) 


        graphs_train <- lapply(celltypes, function(type) visualize_GRN(graph = graph, data = set, color = "means",
                                                                       celltype=type,  layout = layout,
                                                                       rangeMean = c(min(set$means), max(set$means)), 
                                                                       type="Expression"))

    names(graphs_train) <- celltypes

    plot <- ggpubr::ggarrange(plotlist = graphs_train, common.legend = T, legend="bottom",
                              ncol=sqrt(length(celltypes)), nrow= sqrt(length(celltypes)))                   
    return(plot) 
    }

}
                    
get_grn_layout <-  function(graph, cluster=NULL){
    if(is.null(cluster)) cluster <- cluster_louvain(graph)
    weigths <- ifelse(igraph::crossing(cluster, graph), 1 * E(graph)$weight, 100 * E(graph)$weight)
    layout<- layout_with_fr(graph, weights = weigths)
    return(layout)                         
}