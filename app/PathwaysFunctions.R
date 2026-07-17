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

expand.to.ortho.nodes <- function(genes, data.list, pathway.list, ortho.list)
{
  ref.genes <- if(length(genes)==0) character(0) else sapply(strsplit(genes,"\n"),function(x){x[1]})
  ortho.nodes <- c()
  if(length(ref.genes)==0 || "All" %in% ref.genes) return(list(ref.genes=ref.genes, ortho.nodes=ortho.nodes))
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
  list(ref.genes=ref.genes, ortho.nodes=unique(ortho.nodes))
}

get.paths.between.edges <- function(edges.df, source.nodes, dest.nodes, max.length = Inf)
{
  empty.result <- list(edges = edges.df[0,], roles = data.frame(node = character(0), path.role = character(0), stringsAsFactors = FALSE))
  if(nrow(edges.df) == 0) return(empty.result)
  all.nodes <- unique(c(edges.df$source, edges.df$target))
  source.nodes <- intersect(source.nodes, all.nodes)
  dest.nodes <- intersect(dest.nodes, all.nodes)
  if(length(source.nodes) == 0 || length(dest.nodes) == 0) return(empty.result)
  
  g <- igraph::graph_from_data_frame(unique(edges.df[,c("source","target")]), directed = TRUE,
                                     vertices = data.frame(name = all.nodes, stringsAsFactors = FALSE))
  
  dist.from <- function(root.name, mode) {
    root.id <- match(root.name, igraph::V(g)$name)
    bfs.res <- igraph::bfs(g, root = root.id, mode = mode, father = TRUE, dist = TRUE, unreachable = FALSE)
    father <- as.integer(bfs.res$father)
    dist <- as.integer(bfs.res$dist)
    reached <- (seq_along(father) == root.id) | !is.na(father)
    out <- rep(NA_integer_, length(dist))
    out[reached] <- dist[reached]
    names(out) <- igraph::V(g)$name
    out
  }
  
  min.forward <- setNames(rep(NA_integer_, length(igraph::V(g))), igraph::V(g)$name)
  for(src in source.nodes) min.forward <- pmin(min.forward, dist.from(src, "out"), na.rm = TRUE)
  min.backward <- setNames(rep(NA_integer_, length(igraph::V(g))), igraph::V(g)$name)
  for(dst in dest.nodes) min.backward <- pmin(min.backward, dist.from(dst, "in"), na.rm = TRUE)
  
  relevant.nodes <- names(min.forward)[!is.na(min.forward) & !is.na(min.backward)]
  if(length(relevant.nodes) == 0) return(empty.result)
  if(is.finite(max.length)) {
    total.dist <- min.forward[relevant.nodes] + min.backward[relevant.nodes]
    relevant.nodes <- relevant.nodes[total.dist <= max.length]
  }
  if(length(relevant.nodes) == 0) return(empty.result)
  
  relevant.edges <- edges.df[edges.df$source %in% relevant.nodes & edges.df$target %in% relevant.nodes, , drop = FALSE]
  
  is.source <- relevant.nodes %in% source.nodes
  is.dest <- relevant.nodes %in% dest.nodes
  roles.df <- data.frame(node = relevant.nodes,
                         path.role = ifelse(is.source & is.dest, "both", ifelse(is.source, "source", ifelse(is.dest, "destination", NA_character_))),
                         stringsAsFactors = FALSE)
  
  list(edges = relevant.edges, roles = roles.df)
}

get.ego.network.edges <- function(edges.df, ref.nodes, max.hops = 1)
{
  empty.result <- list(edges = edges.df[0,], roles = data.frame(node = character(0), path.role = character(0), stringsAsFactors = FALSE))
  if(nrow(edges.df) == 0) return(empty.result)
  all.nodes <- unique(c(edges.df$source, edges.df$target))
  ref.nodes <- intersect(ref.nodes, all.nodes)
  if(length(ref.nodes) == 0) return(empty.result)
  if(is.null(max.hops) || is.na(max.hops) || max.hops < 0) max.hops <- 1
  
  g <- igraph::graph_from_data_frame(unique(edges.df[,c("source","target")]), directed = TRUE,
                                     vertices = data.frame(name = all.nodes, stringsAsFactors = FALSE))
  
  neighborhood.nodes <- character(0)
  for(ref in ref.nodes) {
    ref.id <- match(ref, igraph::V(g)$name)
    bfs.res <- igraph::bfs(g, root = ref.id, mode = "all", father = TRUE, dist = TRUE, unreachable = FALSE)
    father <- as.integer(bfs.res$father)
    dist <- as.integer(bfs.res$dist)
    reached <- (seq_along(father) == ref.id) | !is.na(father)
    within.hops <- reached & !is.na(dist) & dist <= max.hops
    neighborhood.nodes <- c(neighborhood.nodes, igraph::V(g)$name[which(within.hops)])
  }
  neighborhood.nodes <- unique(neighborhood.nodes)
  if(length(neighborhood.nodes) == 0) return(empty.result)
  
  relevant.edges <- edges.df[edges.df$source %in% neighborhood.nodes & edges.df$target %in% neighborhood.nodes, , drop = FALSE]
  roles.df <- data.frame(node = ref.nodes, path.role = "ego", stringsAsFactors = FALSE)
  
  list(edges = relevant.edges, roles = roles.df)
}

build.pathway.net <- function(data.list,metapathway.list,pathway.list,ortho.list,
                              networks,pathways,genes,hide.elements,view.mode="neighbors",
                              source.genes=NULL,dest.genes=NULL,max.hops=1,max.length=Inf)
{
  multilayer.net <- list()
  
  if(view.mode=="paths") {
    src.expansion <- expand.to.ortho.nodes(source.genes, data.list, pathway.list, ortho.list)
    dst.expansion <- expand.to.ortho.nodes(dest.genes, data.list, pathway.list, ortho.list)
    ref.genes <- unique(c(src.expansion$ref.genes, dst.expansion$ref.genes))
    source.ortho.nodes <- src.expansion$ortho.nodes
    dest.ortho.nodes <- dst.expansion$ortho.nodes
  } else {
    gene.expansion <- expand.to.ortho.nodes(genes, data.list, pathway.list, ortho.list)
    ref.genes <- gene.expansion$ref.genes
    ortho.nodes <- gene.expansion$ortho.nodes
  }
  id.count <- 1
  
  for(net in networks) {
    net.node.data <- data.list[[net]]
    if(is.null(net.node.data)) next
    pathway.info <- pathway.list[[net.node.data$organism]]
    metapathway.info <- metapathway.list[[net.node.data$organism]]
    
    #Retrieve pathway nodes
    pathway.info <- pathway.info[pathway.info$pathwayName %in% pathways,c("node","nodeName","endpoint","node.type")]
    pathway.info <- aggregate(endpoint ~ node + nodeName + node.type, data = pathway.info, FUN = any)
    pathway.nodes.info <- merge(pathway.info,net.node.data$data,all.x=T)
    pathway.nodes.info[is.na(pathway.nodes.info$activity),"activity"] <- 0
    
    #Hide extra elements, if needed
    always.allowed <- if(view.mode=="paths") unique(c(source.ortho.nodes, dest.ortho.nodes)) else ortho.nodes
    pathway.nodes.filtered <- apply.hide.filters(pathway.nodes.info,"node",hide.elements)
    pathway.nodes.info <- pathway.nodes.info[pathway.nodes.info$node %in% pathway.nodes.filtered$node | pathway.nodes.info$node %in% always.allowed,]
    pathway.nodes.info$layer <- rep(net, nrow(pathway.nodes.info))
    
    #Retrieve pathway edges
    pathway.edges.info <- metapathway.info[metapathway.info$source %in% pathway.nodes.info$node & metapathway.info$target %in% pathway.nodes.info$node,]
    node.roles <- NULL
    if(view.mode=="paths") {
      if(length(source.genes)>0 && length(dest.genes)>0) {
        path.result <- get.paths.between.edges(pathway.edges.info, source.ortho.nodes, dest.ortho.nodes, max.length = max.length)
        pathway.edges.info <- path.result$edges
        node.roles <- path.result$roles
      } else {
        pathway.edges.info <- pathway.edges.info[0,]
      }
    } else if(!"All" %in% ref.genes) {
      ego.result <- get.ego.network.edges(pathway.edges.info, ortho.nodes, max.hops = max.hops)
      pathway.edges.info <- ego.result$edges
      node.roles <- ego.result$roles
    }
    if(nrow(pathway.edges.info)>0) {
      pathway.nodes.ids <- unique(c(pathway.edges.info$source,pathway.edges.info$target))
      if(view.mode!="paths" && !"All" %in% ref.genes) pathway.nodes.ids <- unique(c(pathway.nodes.ids, ortho.nodes))
      pathway.nodes.info <- unique(pathway.nodes.info[pathway.nodes.info$node %in% pathway.nodes.ids,])
    } else if(view.mode=="paths") {
      pathway.nodes.info <- pathway.nodes.info[0,]
    } else if(!"All" %in% ref.genes) {
      pathway.nodes.info <- pathway.nodes.info[pathway.nodes.info$node %in% ortho.nodes,]
    }
    
    if(!is.null(node.roles) && nrow(node.roles)>0) {
      pathway.nodes.info <- merge(pathway.nodes.info, node.roles, by="node", all.x=TRUE)
    } else {
      pathway.nodes.info$path.role <- rep(NA_character_, nrow(pathway.nodes.info))
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

get.legend.info <- function(multilayer.nodes)
{
  net.names <- unique(multilayer.nodes$layer)
  num.layers <- length(net.names)
  col.borders <- c("#0072B2","#D55E00","#009E73")
  list(net.names = net.names, colors = col.borders[seq_len(num.layers)])
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
    
    #Node size
    multilayer.nodes$size <- 25
    
    #Set tooltip string for nodes
    multilayer.nodes$title <- paste0("Score: ",round(multilayer.nodes$activity,3))
    
    #Highlight source/destination/reference nodes
    if("path.role" %in% colnames(multilayer.nodes)) {
      is.source <- !is.na(multilayer.nodes$path.role) & multilayer.nodes$path.role %in% c("source","both")
      is.dest   <- !is.na(multilayer.nodes$path.role) & multilayer.nodes$path.role %in% c("destination","both")
      is.ego    <- !is.na(multilayer.nodes$path.role) & multilayer.nodes$path.role == "ego"
      multilayer.nodes$color.border[is.source & !is.dest] <- "#00A651"
      multilayer.nodes$color.border[is.dest & !is.source] <- "#ED1C24"
      multilayer.nodes$color.border[is.source & is.dest]  <- "#F7941D"
      multilayer.nodes$color.border[is.ego] <- "#8E44AD"
      multilayer.nodes$borderWidth[is.source | is.dest | is.ego] <- 6
      role.label <- ifelse(is.source & is.dest, "Source &amp; destination node",
                           ifelse(is.source, "Source node",
                                  ifelse(is.dest, "Destination node",
                                         ifelse(is.ego, "Reference node", NA_character_))))
      multilayer.nodes$title <- ifelse(is.na(role.label), multilayer.nodes$title,
                                       paste0(multilayer.nodes$title,"<br>",role.label))
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

truncate.label <- function(x, max.chars = 60)
{
  ifelse(nchar(x) > max.chars, paste0(substr(x, 1, max.chars - 1), "\u2026"), x)
}

get.list.selectable.nodes <- function(list.pathways,list.networks,data.list,pathway.list)
{
  list.options <- list()
  for(network in list.networks) {
    organism <- data.list[[network]]$organism
    pathway.nodes <- pathway.list[[organism]]
    list.nodes <- unique(pathway.nodes[pathway.nodes$pathwayName %in% list.pathways,c("node","nodeName","node.type")])
    list.nodes <- list.nodes[order(list.nodes$nodeName),]
    final.options <- paste0(list.nodes$nodeName,"\n",network,"\n",list.nodes$node.type)
    names(final.options) <- truncate.label(list.nodes$nodeName)
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
  names(list.options) <- truncate.label(list.options)
  return(list.options)
}

filter.selectable.nodes <- function(starting.list.nodes,list.networks)
{
  starting.list.nodes[list.networks]
}

metapathway.list <- list()
pathway.list <- list()
ortho.list <- list()
map.organism <- read.table("Data/mapOrganism.txt",header=T,sep="\t")
