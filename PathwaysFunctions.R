read.phensim.file <- function(file)
{
  pathways.data <- as.data.frame(fread(file))
  if("Activity Score" %in% colnames(pathways.data)) {
    organism <- map.organism[sapply(map.organism$Code,grepl,tolower(pathways.data[1,1])),"Organism"]
    #organism <- map.organism[map.organism$Code==tolower(unlist(strsplit(pathways.data[1,1],":|-"))[2]),"Organism"]
    pathways.data <- unique(pathways.data[,c("Node Id","Activity Score")])
  } else if("Perturbation" %in% colnames(pathways.data)) {
    organism <- map.organism[sapply(map.organism$Code,grepl,tolower(pathways.data[1,1])),"Organism"]
    #organism <- map.organism[map.organism$Code==tolower(unlist(strsplit(pathways.data[1,1],":|-"))[2]),"Organism"]
    pathways.data <- unique(pathways.data[,c("Gene Id","Perturbation")])
  } else {
    organism.name <- names(fread(file,nrows=1))
    organism <- as.character(map.organism[map.organism$Organism==organism.name,"Organism"])
    pathways.data <- unique(pathways.data)
  }
  colnames(pathways.data) <- c("node","activity")
  final.data <- list(organism=organism,data=pathways.data)
  return(final.data)
}

build.pathway.net <- function(data.list,metapathway.list,pathway.list,ortho.list,
                              networks,pathway,gene,hide.elements)
{
  multilayer.net <- list()
  
  #Get reference genes for gene-centric visualization
  ref.gene <- strsplit(gene,"\n")[[1]][1]
  if(ref.gene!="All") {
    ref.net <- strsplit(gene,"\n")[[1]][2]
    ref.organism <- as.character(data.list[[ref.net]]$organism)
    ref.pathways <- pathway.list[[ref.organism]]
    ref.id <- c(ref.pathways[ref.pathways$pathwayName==pathway & ref.pathways$nodeName==ref.gene,"node"])
    ortho.nodes <- c(ref.id)
    for(net in networks) {
      if(net!=ref.net) {
        net.organism <- as.character(data.list[[net]]$organism)
        if(ref.organism<=net.organism) {
          ref.ortho <- ortho.list[[paste0(ref.organism,"-",net.organism)]]
        } else {
          ref.ortho <- ortho.list[[paste0(net.organism,"-",ref.organism)]]
        }
        ortho.nodes <- c(ortho.nodes,ref.ortho[ref.ortho[,paste0(ref.organism," id")]==ref.id,paste0(net.organism," id")])
      }
    }
  }
  id.count <- 1
  
  for(net in networks)
  {
    net.node.data <- data.list[[net]]
    pathway.info <- pathway.list[[net.node.data$organism]]
    metapathway.info <- metapathway.list[[net.node.data$organism]]
    
    #Retrieve pathway nodes
    pathway.info <- pathway.info[pathway.info$pathwayName==pathway,c("node","nodeName","endpoint")]
    pathway.nodes.info <- merge(pathway.info,net.node.data$data,all.x=T)
    pathway.nodes.info[is.na(pathway.nodes.info$activity),"activity"] <- 0
    #Hide extra element, if needed
    if("Hide chemical entities" %in% hide.elements){
      pathway.nodes.info <- pathway.nodes.info[!startsWith(pathway.nodes.info$node,"chebi:"),]
      pathway.nodes.info <- pathway.nodes.info[!startsWith(pathway.nodes.info$node,"cpd:"),]
      pathway.nodes.info <- pathway.nodes.info[!startsWith(pathway.nodes.info$node,"gl:"),]
    }
    if("Hide drugs" %in% hide.elements){
      pathway.nodes.info <- pathway.nodes.info[!startsWith(pathway.nodes.info$node,"dr:"),]
    }
    if("Hide miRNAs" %in% hide.elements){
      pathway.nodes.info <- pathway.nodes.info[!grepl("-miR",pathway.nodes.info$node),]
      pathway.nodes.info <- pathway.nodes.info[!grepl("-let",pathway.nodes.info$node),]
    }
    pathway.nodes.info$layer <- net
    
    #Retrieve pathway edges
    pathway.edges.info <- metapathway.info[metapathway.info$source %in% pathway.nodes.info$node & metapathway.info$target %in% pathway.nodes.info$node,]
    if(ref.gene!="All") {
      #gene.id <- pathway.nodes.info[pathway.nodes.info$nodeName==gene,"node"]
      gene.edges <- pathway.edges.info[pathway.edges.info$source %in% ortho.nodes | pathway.edges.info$target %in% ortho.nodes,]
      sub.nodes.list <- unique(c(gene.edges$source,gene.edges$target))
      pathway.edges.info <- pathway.edges.info[pathway.edges.info$source %in% sub.nodes.list & pathway.edges.info$target %in% sub.nodes.list,]
    }
    if(nrow(pathway.edges.info)>0) {
      pathway.nodes.ids <- unique(c(pathway.edges.info$source,pathway.edges.info$target))
      pathway.nodes.info <- unique(pathway.nodes.info[pathway.nodes.info$node %in% pathway.nodes.ids,])
    } else if(ref.gene!="All") {
      pathway.nodes.info <- pathway.nodes.info[pathway.nodes.info$node %in% ortho.nodes,]
    }
    
    #Re-map ids
    if(nrow(pathway.nodes.info)>0) {
      pathway.nodes.info$id <- id.count:(id.count+nrow(pathway.nodes.info)-1)
    } else {
      pathway.nodes.info$id <- numeric(0)
    }
    sub.nodes.info <- pathway.nodes.info[,c("node","id")]
    pathway.edges.info <- merge(pathway.edges.info,sub.nodes.info,by.x="source",by.y="node")
    pathway.edges.info <- merge(pathway.edges.info,sub.nodes.info,by.x="target",by.y="node")
    pathway.edges.info <- pathway.edges.info[,c(6,2,3,7,1,4,5)]
    colnames(pathway.edges.info) <- c("sourceId","source","sourceName","targetId","target","targetName","weight")
    id.count <- id.count+nrow(pathway.nodes.info)
    
    #Add node and edge data to the final network
    multilayer.net[[net]] <- list(nodes=pathway.nodes.info,edges=pathway.edges.info)
    
  }
  
  #Merge nodes and edges in a unique multilayer network
  multilayer.nodes <- lapply(multilayer.net,function(el){el$nodes})
  multilayer.nodes <- Reduce(rbind,multilayer.nodes)
  multilayer.nodes <- multilayer.nodes[!duplicated(multilayer.nodes$id),]
  multilayer.edges <- lapply(multilayer.net,function(el){el$edges})
  multilayer.edges <- Reduce(rbind,multilayer.edges)
  if(nrow(multilayer.edges)>0) {
    multilayer.edges$type <- "intra"
  } else {
    multilayer.edges$type <- character(0)
  }
    
  #Create inter-layer connections
  if(length(multilayer.net)>1) {
    org.pairs <- t(combn(names(multilayer.net),2))
    for(i in 1:nrow(org.pairs)) {
      first.org <- data.list[[org.pairs[i,1]]]$organism
      second.org <- data.list[[org.pairs[i,2]]]$organism
      ortho.data <- ortho.list[[paste0(first.org,"-",second.org)]]
      ortho.data <- ortho.data[ortho.data[,paste0(first.org," id")] %in% multilayer.nodes$node & 
                               ortho.data[,paste0(second.org," id")] %in% multilayer.nodes$node,]
      if(nrow(ortho.data)>0) {
        ortho.data <- ortho.data[,c(paste0(first.org," id"),paste0(first.org," symbol"),
                                    paste0(second.org," id"),paste0(second.org," symbol"))]
        colnames(ortho.data) <- c("source","sourceName","target","targetName")
        ortho.data <- merge(ortho.data,multilayer.net[[org.pairs[i,1]]]$nodes[,c("node","id")],by.x="source",by.y="node")
        ortho.data <- merge(ortho.data,multilayer.net[[org.pairs[i,2]]]$nodes[,c("node","id")],by.x="target",by.y="node")
        ortho.data <- ortho.data[,c(5,2,3,6,1,4)]
        colnames(ortho.data) <- c("sourceId","source","sourceName","targetId","target","targetName")
        ortho.data$weight <- 0
        ortho.data$type <- "inter"
        multilayer.edges <- rbind(multilayer.edges,ortho.data)
        multilayer.edges <- multilayer.edges[!duplicated(multilayer.edges[,c("sourceId","targetId")]),]
      }
    }
  }
  
  
  return(list(nodes=multilayer.nodes,edges=multilayer.edges))
  
}

plot.pathway <- function(multilayer.nodes,multilayer.edges,pathway)
{
  if(nrow(multilayer.nodes)>0)
  {
    #Set node layer for hierarchical layout
    freq.layer <- table(multilayer.nodes$layer)
    multilayer.nodes$level <- unlist(lapply(1:length(freq.layer),function(x){rep(x,freq.layer[x])}))
    
    #Set node labels
    colnames(multilayer.nodes)[2] <- "label"
  
    #Set node colors
    if(all(multilayer.nodes$activity==0)) {
      multilayer.nodes$color.background <- "grey"
    } else {
      max.val <- ceiling(max(abs(min(multilayer.nodes$activity)),
                         abs(max(multilayer.nodes$activity))))               
      min.val <- -max.val
      activity.scores.normalized <- (multilayer.nodes$activity-min.val)/(max.val-min.val)
      activity.colors <- colorRamp(c("blue","grey","red"))(activity.scores.normalized)
      activity.colors <- apply(activity.colors, 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
      multilayer.nodes$color.background <- activity.colors
    }
  
    #Set border color for orthologous nodes
    border.colors <- rep("brown",nrow(multilayer.nodes))
    border.colors[multilayer.nodes$level==2] <- "darkgreen"
    border.colors[multilayer.nodes$level==3] <- "gold"
    multilayer.nodes$color.border <- border.colors
  
    #Set node border width
    multilayer.nodes$borderWidth <- 3
  
    #Set a square shape for endpoint nodes
    nodes.shape <- rep("dot",nrow(multilayer.nodes))
    nodes.shape[multilayer.nodes$endpoint==T] <- "square"
    multilayer.nodes$shape <- nodes.shape
    
    #Set tooltip string for nodes
    multilayer.nodes$title <- paste0("Score: ",round(multilayer.nodes$activity,3))
  }
  
  if(nrow(multilayer.edges)>0)
  {
    #Set arrow type for inhibition and activation edges
    arrow.type <- rep("arrow",nrow(multilayer.edges))
    arrow.type[multilayer.edges$weight<0] <- "bar"
    arrow.type[multilayer.edges$type=="inter"] <- NA 
    multilayer.edges$arrows.to.type <- arrow.type
  
    #Set edge style for intra- and inter-layer edges
    multilayer.edges$dashes <- F
    multilayer.edges[multilayer.edges$type=="inter","dashes"] <- T
  }
  
  #Set legend for nodes and edges
  net.names <- unique(multilayer.nodes$layer)
  num.layers <- length(net.names)
  col.borders <- c("brown","darkgreen","gold")
  legend.nodes <- data.frame(label = c(net.names,"Endpoint", "Score"),
                             shape = c(rep("dot",num.layers),"square","image"),
                             color.background = c(rep("grey",num.layers),"white","white"),
                             color.border = c(col.borders[1:num.layers], "black","black"),
                             font.size = c(rep(16,num.layers),16,16),
                             borderWidth = c(rep(3,num.layers),3,3),
                             image=c(rep("",num.layers),"","https://alpha.dmi.unict.it/~gmicale/Documents/activity2.png"))
  
  #Set edge weights for plotting network
  edge.plot.weigths <- rep(1,nrow(multilayer.edges))
  edge.plot.weigths[multilayer.edges$weight==0] <- 2
  
  #Plot pathway
  #colnames(multilayer.nodes)[1] <- "id"
  colnames(multilayer.edges)[1] <- "from"
  colnames(multilayer.edges)[4] <- "to"
  multilayer.plot <- visNetwork(multilayer.nodes, multilayer.edges, main=pathway) %>%
    visPhysics(enabled = F) %>%
    visEdges(arrows="to",color = "black") %>%
    visOptions(clickToUse = F) %>%
    visInteraction(hover = TRUE) %>%
    visLegend(addNodes = legend.nodes, useGroups = FALSE, width=0.25, position="right", zoom=F) %>%
    visNodes(font=list(color="black", size=20))
  if(num.layers==1) {
    multilayer.plot <- visIgraphLayout(multilayer.plot,layout="layout_with_kk", type = "full", weights=edge.plot.weigths)
  } else {
    multilayer.plot <- visIgraphLayout(multilayer.plot,layout="layout_with_sugiyama", layers = multilayer.nodes$level, type="full")
    window <- 1
    span <- 0.25
    net.space <- window+span
    lower.limit <- -(net.space*num.layers/2)
    multilayer.plot$x$nodes[multilayer.plot$x$nodes$level==1,"y"] <- runif(sum(multilayer.plot$x$nodes$level==1),lower.limit,lower.limit+window)
    tmp <- multilayer.edges[multilayer.edges$type=="inter",]
    for(i in 2:num.layers) {
      sub.tmp <- tmp[tmp$from %in% multilayer.nodes[multilayer.nodes$level==(i-1),"id"] & tmp$to %in% multilayer.nodes[multilayer.nodes$level==i,"id"],]
      list.ids <- sub.tmp[match(multilayer.plot$x$nodes[multilayer.plot$x$nodes$level==i & multilayer.plot$x$nodes$id %in% multilayer.edges[multilayer.edges$type=="inter","to"],"id"],sub.tmp$to),"from"]
      multilayer.plot$x$nodes[multilayer.plot$x$nodes$level==i & multilayer.plot$x$nodes$id %in% multilayer.edges[multilayer.edges$type=="inter","to"],"y"][!is.na(list.ids)] <- net.space+multilayer.plot$x$nodes[match(list.ids,multilayer.plot$x$nodes$id),"y"][!is.na(list.ids)]
    }
  }
  
  return(multilayer.plot)
  
}

get.list.selectable.nodes <- function(pathway,list.networks,hide.elements)
{
  list.options <- list()
  list.options$All <- c("All")
  for(network in list.networks) {
    organism <- data.list[[network]]$organism
    pathway.nodes <- pathway.list[[organism]]
    list.nodes <- pathway.nodes[pathway.nodes$pathwayName==pathway,c("node","nodeName")]
    if("Hide chemical entities" %in% hide.elements) {
      list.nodes <- list.nodes[!startsWith(list.nodes$node,"chebi:"),]
      list.nodes <- list.nodes[!startsWith(list.nodes$node,"cpd:"),]
      list.nodes <- list.nodes[!startsWith(list.nodes$node,"gl:"),]
    }
    if("Hide drugs" %in% hide.elements) {
      list.nodes <- list.nodes[!startsWith(list.nodes$node,"dr:"),]
    }
    if("Hide miRNAs" %in% hide.elements) {
      list.nodes <- list.nodes[!grepl("-miR",list.nodes$node),]
      list.nodes <- list.nodes[!grepl("-let",list.nodes$node),]
    }
    list.nodes <- sort(list.nodes$nodeName)
    final.options <- paste0(list.nodes,"\n",network)
    names(final.options) <- list.nodes
    list.options[[network]] <- final.options
  }
  return(list.options)
}

metapathway.list <- list()
pathway.list <- list()
ortho.list <- list()
data.list <- list()
map.organism <- read.table("Data/mapOrganism.txt",header=T,sep="\t")
