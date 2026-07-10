classify.node.type <- function(node.ids)
{
  node.type <- rep("gene", length(node.ids))
  node.type[grepl("^(chebi:|cpd:|gl:)", node.ids, ignore.case = TRUE)] <- "chemical"
  node.type[startsWith(node.ids, "dr:")] <- "drug"
  is.mirna <- grepl("-mir|-let|mir-|let-", node.ids, ignore.case = TRUE) & node.type == "gene"
  node.type[is.mirna] <- "mirna"
  return(node.type)
}

apply.hide.filters <- function(df, node.col, hide.elements)
{
  if(nrow(df) == 0) return(df)
  if(!"node.type" %in% colnames(df)) {
    df$node.type <- classify.node.type(df[[node.col]])
  }
  hide.types <- c()
  if("Hide chemical entities" %in% hide.elements) hide.types <- c(hide.types, "chemical")
  if("Hide drugs" %in% hide.elements)             hide.types <- c(hide.types, "drug")
  if("Hide miRNAs" %in% hide.elements)            hide.types <- c(hide.types, "mirna")
  if(length(hide.types) == 0) return(df)
  return(df[!df$node.type %in% hide.types, , drop = FALSE])
}

resolve.organism.by.code <- function(id.string)
{
  is.match <- sapply(map.organism$Code, function(code) grepl(code, id.string, fixed = TRUE))
  matches <- map.organism$Code[is.match]
  if(length(matches) == 0) return(NA_character_)
  as.character(map.organism[map.organism$Code == matches[1], "Organism"][1])
}

read.phensim.file <- function(file)
{
  #Read file header
  con <- file(file,"r")
  file.header <- readLines(con,n=1)
  close(con)
  #Read data according to the file format
  if(grepl("Activity Score",file.header,fixed=T)) {
    pathways.data <- as.data.frame(fread(file))
    organism <- resolve.organism.by.code(tolower(pathways.data[1,1]))
    pathways.data <- unique(pathways.data[,c("Node Id","Activity Score")])
  } else if(grepl("Perturbation",file.header,fixed=T)) {
    pathways.data <- as.data.frame(fread(file))
    organism <- resolve.organism.by.code(tolower(pathways.data[1,1]))
    pathways.data <- unique(pathways.data[,c("Gene Id","Perturbation")])
  } else {
    pathways.data <- as.data.frame(fread(file,skip=1))
    pathways.data <- unique(pathways.data)
    organism.name <- file.header
    organism.row <- map.organism[map.organism$Organism==organism.name,"Organism"]
    organism <- if(length(organism.row)==0) NA_character_ else as.character(organism.row[1])
  }
  colnames(pathways.data) <- c("node","activity")
  final.data <- list(organism=organism,data=pathways.data)
  return(final.data)
}

get.ancestor.tree.edges <- function(edges.df, target.nodes, max.hops = Inf)
{
  empty.result <- list(edges = edges.df[0,], hops = data.frame(node = character(0), hop = numeric(0), stringsAsFactors = FALSE))
  if(nrow(edges.df) == 0) return(empty.result)
  all.nodes <- unique(c(edges.df$source, edges.df$target))
  target.nodes <- intersect(target.nodes, all.nodes)
  if(length(target.nodes) == 0) return(empty.result)
  
  g <- igraph::graph_from_data_frame(unique(edges.df[,c("source","target")]), directed = TRUE,
                                     vertices = data.frame(name = all.nodes, stringsAsFactors = FALSE))
  g.super.sink <- igraph::add_vertices(g, 1, name = "__SUPERSINK__")
  super.sink.id <- which(igraph::V(g.super.sink)$name == "__SUPERSINK__")
  target.ids <- match(target.nodes, igraph::V(g.super.sink)$name)
  g.super.sink <- igraph::add_edges(g.super.sink, as.vector(rbind(target.ids, rep(super.sink.id, length(target.ids)))))
  
  bfs.res <- igraph::bfs(g.super.sink, root = super.sink.id, mode = "in", father = TRUE, dist = TRUE)
  father <- as.integer(bfs.res$father)
  dist <- as.integer(bfs.res$dist)
  
  discovered <- setdiff(which(!is.na(father)), super.sink.id)
  if(length(discovered) == 0) return(empty.result)
  if(is.finite(max.hops)) {
    discovered <- discovered[dist[discovered] <= (max.hops + 1)]
  }
  if(length(discovered) == 0) return(empty.result)
  
  node.names <- igraph::V(g.super.sink)$name[discovered]
  hops.df <- data.frame(node = node.names, hop = dist[discovered] - 1, stringsAsFactors = FALSE)
  
  edge.src <- node.names
  edge.tgt <- igraph::V(g.super.sink)$name[father[discovered]]
  keep <- edge.tgt != "__SUPERSINK__"  #scarta gli "archi" verso il nodo virtuale
  tree.pairs <- paste(edge.src[keep], edge.tgt[keep], sep = "->")
  
  tree.edges <- edges.df[paste(edges.df$source, edges.df$target, sep = "->") %in% tree.pairs, , drop = FALSE]
  list(edges = tree.edges, hops = hops.df)
}

build.pathway.net <- function(data.list,metapathway.list,pathway.list,ortho.list,
                              networks,pathways,genes,hide.elements,view.mode="neighbors",max.hops=Inf)
{
  multilayer.net <- list()
  
  #Get reference genes for gene-centric visualization
  ref.genes <- sapply(strsplit(genes,"\n"),function(x){x[1]})
  ortho.nodes <- c()
  if(!"All" %in% ref.genes) {
    ref.nets <- unique(sapply(strsplit(genes,"\n"),function(x){x[2]}))
    ref.nets <- ref.nets[ref.nets %in% names(data.list)]
    if(length(ref.nets) > 0) {
      ref.organisms <- unname(sapply(data.list[ref.nets],function(x){x$organism}))
      ref.organisms <- unique(ref.organisms[!is.na(ref.organisms) & ref.organisms %in% names(pathway.list)])
      if(length(ref.organisms) > 0) {
        ref.pathways <- pathway.list[ref.organisms]
        ref.ids <- lapply(ref.pathways,function(x){unique(x[x$nodeName %in% ref.genes,"node"])})
        ortho.nodes <- c(unname(unlist(ref.ids)))
        for(ref.organism in names(ref.ids)) {
          for(net.organism in names(pathway.list)) {
            if(net.organism!=ref.organism) {
              if(ref.organism<=net.organism) {
                ref.ortho <- ortho.list[[paste0(ref.organism,"-",net.organism)]]
              } else {
                ref.ortho <- ortho.list[[paste0(net.organism,"-",ref.organism)]]
              }
              if(!is.null(ref.ortho) && nrow(ref.ortho) > 0) {
                ortho.nodes <- c(ortho.nodes,ref.ortho[ref.ortho[,paste0(ref.organism," id")] %in% ref.ids[[ref.organism]],paste0(net.organism," id")])
              }
            }
          }
        }
      }
    }
  }
  id.count <- 1
  
  for(net in networks) {
    net.node.data <- data.list[[net]]
    if(is.null(net.node.data)) next
    pathway.info <- pathway.list[[net.node.data$organism]]
    metapathway.info <- metapathway.list[[net.node.data$organism]]
    
    #Retrieve pathway nodes
    pathway.info <- unique(pathway.info[pathway.info$pathwayName %in% pathways,c("node","nodeName","endpoint","node.type")])
    pathway.nodes.info <- merge(pathway.info,net.node.data$data,all.x=T)
    pathway.nodes.info[is.na(pathway.nodes.info$activity),"activity"] <- 0
    
    #Hide extra elements, if needed (single shared implementation, see apply.hide.filters)
    pathway.nodes.info <- apply.hide.filters(pathway.nodes.info,"node",hide.elements)
    pathway.nodes.info$layer <- net
    
    #Retrieve pathway edges
    pathway.edges.info <- metapathway.info[metapathway.info$source %in% pathway.nodes.info$node & metapathway.info$target %in% pathway.nodes.info$node,]
    node.hops <- NULL
    if(!"All" %in% ref.genes) {
      if(view.mode=="paths") {
        ancestor.result <- get.ancestor.tree.edges(pathway.edges.info, ortho.nodes, max.hops = max.hops)
        pathway.edges.info <- ancestor.result$edges
        node.hops <- ancestor.result$hops
      } else {
        gene.edges <- pathway.edges.info[pathway.edges.info$source %in% ortho.nodes | pathway.edges.info$target %in% ortho.nodes,]
        sub.nodes.list <- unique(c(gene.edges$source,gene.edges$target))
        pathway.edges.info <- pathway.edges.info[pathway.edges.info$source %in% sub.nodes.list & pathway.edges.info$target %in% sub.nodes.list,]
      }
    }
    if(nrow(pathway.edges.info)>0) {
      pathway.nodes.ids <- unique(c(pathway.edges.info$source,pathway.edges.info$target))
      pathway.nodes.info <- unique(pathway.nodes.info[pathway.nodes.info$node %in% pathway.nodes.ids | pathway.nodes.info$node %in% ortho.nodes,])
    } else if(!"All" %in% ref.genes) {
      pathway.nodes.info <- pathway.nodes.info[pathway.nodes.info$node %in% ortho.nodes,]
    }
    
    if(!is.null(node.hops) && nrow(node.hops)>0) {
      pathway.nodes.info <- merge(pathway.nodes.info, node.hops, by="node", all.x=TRUE)
    } else {
      pathway.nodes.info$hop <- NA_real_
    }
    
    #Re-map ids
    if(nrow(pathway.nodes.info)>0) {
      pathway.nodes.info$id <- id.count:(id.count+nrow(pathway.nodes.info)-1)
    } else {
      pathway.nodes.info$id <- numeric(0)
    }
    sub.nodes.info <- pathway.nodes.info[,c("node","id")]
    pathway.edges.info <- merge(pathway.edges.info,sub.nodes.info,by.x="source",by.y="node")
    colnames(pathway.edges.info)[colnames(pathway.edges.info)=="id"] <- "sourceId"
    pathway.edges.info <- merge(pathway.edges.info,sub.nodes.info,by.x="target",by.y="node")
    colnames(pathway.edges.info)[colnames(pathway.edges.info)=="id"] <- "targetId"
    pathway.edges.info <- pathway.edges.info[,c("sourceId","source","sourceName","targetId","target","targetName","weight")]
    id.count <- id.count+nrow(pathway.nodes.info)
    
    #Add node and edge data to the final network
    multilayer.net[[net]] <- list(nodes=pathway.nodes.info,edges=pathway.edges.info)
    
  }
  
  if(length(multilayer.net) == 0) {
    empty.nodes <- data.frame(node=character(),nodeName=character(),endpoint=logical(),
                              node.type=character(),activity=numeric(),layer=character(),
                              id=integer(),stringsAsFactors=FALSE)
    empty.edges <- data.frame(sourceId=integer(),source=character(),sourceName=character(),
                              targetId=integer(),target=character(),targetName=character(),
                              weight=numeric(),type=character(),stringsAsFactors=FALSE)
    return(list(nodes=empty.nodes,edges=empty.edges))
  }
  
  #Merge nodes and edges in a unique multilayer network
  multilayer.nodes <- lapply(multilayer.net,function(el){el$nodes})
  multilayer.nodes <- do.call(rbind,multilayer.nodes)
  multilayer.nodes <- multilayer.nodes[!duplicated(multilayer.nodes$id),]
  multilayer.edges <- lapply(multilayer.net,function(el){el$edges})
  multilayer.edges <- do.call(rbind,multilayer.edges)
  if(nrow(multilayer.edges)>0) {
    multilayer.edges$type <- "intra"
  } else {
    multilayer.edges$type <- character(0)
  }
  
  #Create inter-layer connections
  if(length(multilayer.net)>1) {
    net.pairs <- t(combn(names(multilayer.net),2))
    for(i in 1:nrow(net.pairs)) {
      net1 <- net.pairs[i,1]
      net2 <- net.pairs[i,2]
      org1 <- data.list[[net1]]$organism
      org2 <- data.list[[net2]]$organism
      
      sorted.orgs <- sort(c(org1,org2))
      first.org <- sorted.orgs[1]
      second.org <- sorted.orgs[2]
      net.for.first  <- if(org1==first.org) net1 else net2
      net.for.second <- if(org1==first.org) net2 else net1
      
      ortho.data <- ortho.list[[paste0(first.org,"-",second.org)]]
      if(is.null(ortho.data) || nrow(ortho.data)==0) next
      
      ortho.data <- ortho.data[ortho.data[,paste0(first.org," id")] %in% multilayer.nodes$node &
                                 ortho.data[,paste0(second.org," id")] %in% multilayer.nodes$node,]
      if(nrow(ortho.data)>0) {
        ortho.data <- ortho.data[,c(paste0(first.org," id"),paste0(first.org," symbol"),
                                    paste0(second.org," id"),paste0(second.org," symbol"))]
        colnames(ortho.data) <- c("source","sourceName","target","targetName")
        ortho.data <- merge(ortho.data,multilayer.net[[net.for.first]]$nodes[,c("node","id")],by.x="source",by.y="node")
        colnames(ortho.data)[colnames(ortho.data)=="id"] <- "sourceId"
        ortho.data <- merge(ortho.data,multilayer.net[[net.for.second]]$nodes[,c("node","id")],by.x="target",by.y="node")
        colnames(ortho.data)[colnames(ortho.data)=="id"] <- "targetId"
        ortho.data <- ortho.data[,c("sourceId","source","sourceName","targetId","target","targetName")]
        ortho.data$weight <- 0
        ortho.data$type <- "inter"
        multilayer.edges <- rbind(multilayer.edges,ortho.data)
        multilayer.edges <- multilayer.edges[!duplicated(multilayer.edges[,c("sourceId","targetId")]),]
      }
    }
  }
  
  return(list(nodes=multilayer.nodes,edges=multilayer.edges))
  
}

plot.pathway <- function(multilayer.nodes,multilayer.edges,background="white",view.mode="neighbors")
{
  if(nrow(multilayer.nodes)==0) {
    return(visNetwork(data.frame(id=integer(),label=character()),
                      data.frame(from=integer(),to=integer()),
                      background = background) %>%
             visPhysics(enabled = F))
  }
  
  if(nrow(multilayer.nodes)>0)
  {
    layer.order <- unique(multilayer.nodes$layer)
    multilayer.nodes$level <- match(multilayer.nodes$layer, layer.order)
    
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
    
    #Set border color for orthologous nodes.
    border.colors <- rep("#0072B2",nrow(multilayer.nodes))
    border.colors[multilayer.nodes$level==2] <- "#D55E00"
    border.colors[multilayer.nodes$level==3] <- "#009E73"
    multilayer.nodes$color.border <- border.colors
    
    #Set node border width
    multilayer.nodes$borderWidth <- 3
    
    #Set a square shape for endpoint nodes
    nodes.shape <- rep("dot",nrow(multilayer.nodes))
    nodes.shape[multilayer.nodes$endpoint==T] <- "square"
    multilayer.nodes$shape <- nodes.shape
    
    #Set node size by hop distance
    default.size <- 25
    if("hop" %in% colnames(multilayer.nodes) && any(!is.na(multilayer.nodes$hop))) {
      multilayer.nodes$size <- ifelse(is.na(multilayer.nodes$hop), default.size,
                                      pmax(40 - multilayer.nodes$hop*8, 12))
    } else {
      multilayer.nodes$size <- default.size
    }
    
    #Set tooltip string for nodes
    multilayer.nodes$title <- paste0("Score: ",round(multilayer.nodes$activity,3))
    if("hop" %in% colnames(multilayer.nodes)) {
      multilayer.nodes$title <- ifelse(is.na(multilayer.nodes$hop), multilayer.nodes$title,
                                       paste0(multilayer.nodes$title,"<br>Hop: ",multilayer.nodes$hop))
    }
  }
  
  if(nrow(multilayer.edges)>0)
  {
    #Set arrow type for inhibition and activation edges
    arrow.type <- rep("arrow",nrow(multilayer.edges))
    arrow.type[multilayer.edges$weight<0] <- "bar"
    arrow.type[multilayer.edges$type=="intra" & multilayer.edges$weight==0] <- NA
    arrow.type[multilayer.edges$type=="inter"] <- NA
    multilayer.edges$arrows.to.type <- arrow.type
    
    #Set edge style for intra- and inter-layer edges
    multilayer.edges$dashes <- F
    multilayer.edges[multilayer.edges$type=="inter","dashes"] <- T
  }
  
  #Set legend for nodes and edges
  net.names <- unique(multilayer.nodes$layer)
  num.layers <- length(net.names)
  col.borders <- c("#0072B2","#D55E00","#009E73")
  legend.nodes <- data.frame(label = c(net.names,"Endpoint", "Score"),
                             shape = c(rep("dot",num.layers),"square","image"),
                             color.background = c(rep("grey",num.layers),"white","white"),
                             color.border = c(col.borders[1:num.layers], "black","black"),
                             font.size = c(rep(16,num.layers),16,16),
                             borderWidth = c(rep(3,num.layers),3,3),
                             image=c(rep("",num.layers),"","Icons/ScoreBar.png"))
  
  #Set edge weights for plotting network
  edge.plot.weigths <- rep(1,nrow(multilayer.edges))
  edge.plot.weigths[multilayer.edges$weight==0] <- 2
  
  #Plot pathway
  colnames(multilayer.edges)[1] <- "from"
  colnames(multilayer.edges)[4] <- "to"
  multilayer.plot <- visNetwork(multilayer.nodes, multilayer.edges, background = background) %>%
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

get.list.selectable.nodes <- function(list.pathways,list.networks,hide.elements,data.list,pathway.list)
{
  list.options <- list()
  for(network in list.networks) {
    organism <- data.list[[network]]$organism
    pathway.nodes <- pathway.list[[organism]]
    list.nodes <- unique(pathway.nodes[pathway.nodes$pathwayName %in% list.pathways,c("node","nodeName","node.type")])
    list.nodes <- apply.hide.filters(list.nodes,"node",hide.elements)
    list.nodes <- list.nodes[order(list.nodes$nodeName),]
    final.options <- paste0(list.nodes$nodeName,"\n",network,"\n",list.nodes$node.type)
    names(final.options) <- list.nodes$nodeName
    list.options[[network]] <- final.options
  }
  return(list.options)
}

get.list.selectable.pathways <- function(list.nodes,data.list,pathway.list)
{
  #Group nodes by network
  nodes.info <- do.call(rbind,strsplit(list.nodes, "\n"))
  nodes.info <- split(nodes.info[,1],nodes.info[,2])
  
  #For each organism, get list of selectable pathways
  list.options <- unique(unlist(sapply(names(nodes.info),function(network){
    nodes.list <- nodes.info[[network]]
    organism <- data.list[[network]]$organism
    pathway.info <- pathway.list[[organism]]
    sel.pathways <- unique(pathway.info[pathway.info$nodeName %in% nodes.list,"pathwayName"])
    return(sel.pathways)
  })))
  return(list.options)
}

filter.selectable.nodes <- function(starting.list.nodes,list.networks,hide.elements)
{
  list.selectable.nodes <- starting.list.nodes[list.networks]
  list.selectable.nodes <- lapply(list.selectable.nodes,function(el){
    if(length(el)==0) return(el)
    parts <- do.call(rbind,strsplit(el,"\n"))
    el.genes <- parts[,1]
    ref.net  <- unique(parts[,2])
    el.types <- if(ncol(parts)>=3) parts[,3] else classify.node.type(el.genes)
    
    hide.types <- c()
    if("Hide chemical entities" %in% hide.elements) hide.types <- c(hide.types,"chemical")
    if("Hide drugs" %in% hide.elements)             hide.types <- c(hide.types,"drug")
    if("Hide miRNAs" %in% hide.elements)            hide.types <- c(hide.types,"mirna")
    
    keep <- !el.types %in% hide.types
    el.genes <- el.genes[keep]
    
    final.options <- paste0(el.genes,"\n",ref.net)
    names(final.options) <- el.genes
    return(final.options)
  })
  return(list.selectable.nodes)
}

metapathway.list <- list()
pathway.list <- list()
ortho.list <- list()
map.organism <- read.table("Data/mapOrganism.txt",header=T,sep="\t")
