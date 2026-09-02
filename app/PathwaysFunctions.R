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

#-----------------------------------------------------------------------------
# Dato un vettore di riferimenti nodo nel formato "nodeName\nnetwork" (come
# restituito dai picker pathwaySel/geneSel/sourceSel/destSel), espande la
# selezione ai nodi ortologhi corrispondenti in TUTTI gli organismi
# rappresentati in pathway.list, cosi' che lo stesso gene selezionato in un
# solo layer/organismo venga comunque riconosciuto negli altri layer
# caricati per organismi diversi. Estratta come funzione a se' perche' la
# modalita' "paths" (cammini da sorgenti a destinazioni) deve ripetere
# questa stessa espansione due volte, una per i nodi sorgente e una per i
# nodi destinazione.
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# Dato un elenco di archi diretti (colonne "source","target") e due insiemi
# di nodi (source.nodes, dest.nodes), restituisce SOLO i nodi/archi che
# giacciono su ALMENO UN cammino diretto da un nodo sorgente a un nodo
# destinazione - cioe' l'intersezione tra "raggiungibile in avanti da
# qualche sorgente" e "puo' raggiungere all'indietro qualche destinazione".
# Un nodo che soddisfa entrambe le condizioni e' per definizione
# attraversato da almeno un cammino sorgente->destinazione: qualunque arco
# tra due nodi di questo insieme e' quindi utilizzabile su un simile
# cammino, e viene incluso nel risultato (sottografo indotto).
#
# max.length (default Inf, nessun limite) scarta i nodi la cui distanza
# combinata (dalla sorgente piu' vicina + verso la destinazione piu'
# vicina) supera il limite - un'approssimazione pratica e standard quando
# ci sono piu' sorgenti/destinazioni contemporaneamente (il cammino piu'
# breve complessivo che attraversa il nodo potrebbe usare una sorgente per
# il primo tratto e una destinazione diversa per il secondo).
#
# Restituisce anche un data.frame "roles" (node, path.role) che marca
# ciascun nodo incluso come "source", "destination", o "both" (se lo
# stesso nodo e' presente in entrambi gli insiemi) - usato per evidenziare
# visivamente sorgenti e destinazioni nel grafico.
#-----------------------------------------------------------------------------
get.paths.between.edges <- function(edges.df, source.nodes, dest.nodes, max.length = Inf, perturbed.nodes = NULL)
{
  empty.result <- list(edges = edges.df[0,], roles = data.frame(node = character(0), path.role = character(0), stringsAsFactors = FALSE))
  if(nrow(edges.df) == 0) return(empty.result)
  all.nodes <- unique(c(edges.df$source, edges.df$target))
  source.nodes <- intersect(source.nodes, all.nodes)
  dest.nodes <- intersect(dest.nodes, all.nodes)
  if(length(source.nodes) == 0 || length(dest.nodes) == 0) return(empty.result)
  
  g <- igraph::graph_from_data_frame(unique(edges.df[,c("source","target")]), directed = TRUE,
                                     vertices = data.frame(name = all.nodes, stringsAsFactors = FALSE))
  
  #BUG FIX (segnalato in produzione: "sembra che la ricerca dei cammini non
  #tenga conto della direzione degli archi"): le pathway biologiche
  #contengono spesso cicli/anelli di retroazione. Con una BFS "normale"
  #(illimitata), la ricerca in avanti dalle sorgenti puo' PROSEGUIRE OLTRE
  #una destinazione gia' raggiunta seguendo un ciclo, "rientrando" poi da
  #un'altra direzione e facendo risultare raggiungibili all'indietro anche
  #nodi che in realta' si trovano solo "a valle" di una destinazione -
  #cioe' su un cammino che passa GIA' per la destinazione e prosegue oltre,
  #il che non ha senso per una query "mostrami i cammini da X a Y". La BFS
  #qui sotto si comporta come una normale BFS (stessa distanza calcolata
  #per ogni nodo raggiunto), con l'unica differenza che, una volta
  #raggiunto un nodo "terminale" (le destinazioni per la ricerca in avanti,
  #le sorgenti per quella all'indietro), non prosegue l'esplorazione OLTRE
  #quel nodo - lo marca comunque come raggiunto, semplicemente non lo usa
  #per scoprire ulteriori nodi.
  bfs.stop.at <- function(start.node, mode, terminal.nodes) {
    node.names <- igraph::V(g)$name
    dist <- setNames(rep(NA_integer_, length(node.names)), node.names)
    start.id <- match(start.node, node.names)
    if(is.na(start.id)) return(dist)
    dist[start.node] <- 0L
    frontier <- start.id
    #BUG FIX (segnalato in produzione: selezionando un nodo sia come
    #sorgente sia come una delle destinazioni - es. ISG15 sorgente, e
    #ISG15 stesso incluso tra le destinazioni insieme a UBB - la ricerca
    #non trovava PIU' ALCUN cammino, nemmeno quello verso UBB che esisteva
    #prima di aggiungere ISG15 alle destinazioni). Il nodo di PARTENZA di
    #ciascuna ricerca va sempre escluso dal proprio insieme "terminale":
    #la regola "non proseguire oltre un terminale gia' raggiunto" si
    #applica ai nodi INCONTRATI durante l'esplorazione (per non
    #proseguire oltre una destinazione gia' raggiunta), non al punto di
    #partenza stesso - da cui dobbiamo sempre poter esplorare, anche se e'
    #esso stesso "terminale" per un ALTRO ruolo (qui: sorgente=ISG15 e
    #ISG15 e' ANCHE tra le destinazioni).
    terminal.set <- setdiff(intersect(terminal.nodes, node.names), start.node)
    d <- 0L
    repeat {
      expandable <- frontier[!(node.names[frontier] %in% terminal.set)]
      if(length(expandable) == 0) break
      neighbor.ids <- unique(unlist(igraph::adjacent_vertices(g, expandable, mode = mode)))
      neighbor.ids <- neighbor.ids[is.na(dist[node.names[neighbor.ids]])]
      if(length(neighbor.ids) == 0) break
      d <- d + 1L
      dist[node.names[neighbor.ids]] <- d
      frontier <- neighbor.ids
    }
    dist
  }
  
  #Distanza minima da QUALUNQUE sorgente (in avanti, fermandosi alle
  #destinazioni) e verso QUALUNQUE destinazione (all'indietro, fermandosi
  #alle sorgenti), per ciascun nodo del grafo.
  min.forward <- setNames(rep(NA_integer_, length(igraph::V(g))), igraph::V(g)$name)
  for(src in source.nodes) min.forward <- pmin(min.forward, bfs.stop.at(src, "out", dest.nodes), na.rm = TRUE)
  min.backward <- setNames(rep(NA_integer_, length(igraph::V(g))), igraph::V(g)$name)
  for(dst in dest.nodes) min.backward <- pmin(min.backward, bfs.stop.at(dst, "in", source.nodes), na.rm = TRUE)
  
  relevant.nodes <- names(min.forward)[!is.na(min.forward) & !is.na(min.backward)]
  if(length(relevant.nodes) == 0) return(empty.result)
  if(is.finite(max.length)) {
    total.dist <- min.forward[relevant.nodes] + min.backward[relevant.nodes]
    relevant.nodes <- relevant.nodes[total.dist <= max.length]
  }
  if(length(relevant.nodes) == 0) return(empty.result)
  
  #Filtro opzionale: mantieni solo i nodi che giacciono su un cammino che
  #attraversa ALMENO UN nodo perturbato (score diverso da zero). Un nodo e'
  #"utile" per questo filtro se e' esso stesso perturbato, oppure se puo'
  #raggiungere un nodo perturbato restando comunque in grado di
  #raggiungere una destinazione (e' "prima" di un nodo perturbato su
  #qualche cammino), oppure se e' raggiungibile da un nodo perturbato
  #restando comunque raggiungibile da una sorgente (e' "dopo" un nodo
  #perturbato su qualche cammino).
  #BUG FIX: NULL (filtro disattivato) e un vettore VUOTO ma non-NULL
  #(filtro attivo, ma nessun nodo perturbato esiste affatto nella pathway)
  #sono due situazioni DIVERSE - la seconda deve dare un risultato vuoto
  #(nessun cammino passa per un nodo perturbato se non ce n'e' nessuno),
  #non essere trattata come "filtro assente". Il controllo va quindi fatto
  #solo su is.null(), non sulla lunghezza.
  if(!is.null(perturbed.nodes)) {
    perturbed.relevant <- intersect(perturbed.nodes, relevant.nodes)
    if(length(perturbed.relevant) == 0) return(empty.result)
    
    reaches.perturbed <- setNames(rep(FALSE, length(relevant.nodes)), relevant.nodes)
    reached.from.perturbed <- setNames(rep(FALSE, length(relevant.nodes)), relevant.nodes)
    for(p in perturbed.relevant) {
      d.fwd <- bfs.stop.at(p, "out", dest.nodes)
      reached <- intersect(names(d.fwd)[!is.na(d.fwd)], relevant.nodes)
      reached.from.perturbed[reached] <- TRUE
      d.bwd <- bfs.stop.at(p, "in", source.nodes)
      reaches <- intersect(names(d.bwd)[!is.na(d.bwd)], relevant.nodes)
      reaches.perturbed[reaches] <- TRUE
    }
    useful <- (relevant.nodes %in% perturbed.relevant) | reaches.perturbed[relevant.nodes] | reached.from.perturbed[relevant.nodes]
    relevant.nodes <- relevant.nodes[useful]
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

#-----------------------------------------------------------------------------
# Dato un elenco di archi diretti e un insieme di nodi di riferimento,
# restituisce il sottografo indotto su tutti i nodi entro max.hops passi
# (in ENTRAMBE le direzioni, cioe' ignorando il verso degli archi - a
# differenza di get.paths.between.edges, qui non c'e' una direzione
# sorgente->destinazione, solo "vicinanza" al nodo di riferimento) da
# ALMENO UNO dei nodi di riferimento. Con max.hops=1 riproduce esattamente
# il comportamento storico della modalita' Ego-network (nodo selezionato +
# vicini immediati, con gli archi tra i vicini stessi inclusi).
#
# Restituisce anche un data.frame "roles" (node, path.role="ego") che
# marca i nodi di riferimento stessi, per evidenziarli visivamente nel
# grafico.
#-----------------------------------------------------------------------------
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
    #mode="all": ignora il verso degli archi, cosi' da includere sia
    #predecessori che successori del nodo di riferimento, come nella
    #logica storica del "vicinato" (non direzionale).
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
                              source.genes=NULL,dest.genes=NULL,max.hops=1,max.length=Inf,
                              only.perturbed.paths=FALSE)
{
  multilayer.net <- list()
  
  #Get reference genes for gene-centric visualization. In modalita' "paths"
  #(cammini da sorgenti a destinazioni) l'espansione agli ortologhi va
  #fatta separatamente per i nodi sorgente e per quelli destinazione;
  #altrimenti (modalita' "neighbors"/ego-network) si usa il singolo
  #elenco "genes" come prima.
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
    
    #Elenco dei nodi perturbati (score diverso da zero), catturato PRIMA
    #del filtro "Show ..." qui sotto - cosi' un nodo perturbato ma di tipo
    #nascosto (es. un composto chimico con "Show compounds" non spuntato)
    #conta comunque come "perturbato" ai fini del filtro opzionale "solo
    #cammini che passano per un nodo perturbato", coerentemente con il
    #principio gia' seguito altrove: i filtri "Show ..." controllano solo
    #cosa viene DISEGNATO, non cosa viene cercato/calcolato.
    all.perturbed.nodes <- pathway.nodes.info[!is.na(pathway.nodes.info$activity) & pathway.nodes.info$activity != 0, "node"]
    
    #Hide extra elements, if needed (single shared implementation, see apply.hide.filters)
    #BUG FIX (segnalato in produzione per la modalita' "paths": selezionando
    #ad es. un miRNA come sorgente e un gene come destinazione, senza "Show
    #miRNAs" spuntato, non veniva mostrato nulla). I nodi sorgente/
    #destinazione (modalita' "paths") o di riferimento (modalita'
    #"neighbors") non devono essere esclusi dal filtro "Show ...", anche se
    #sono di un tipo nascosto - solo i nodi INTERMEDI/di contorno devono
    #rispettarlo.
    always.allowed <- if(view.mode=="paths") unique(c(source.ortho.nodes, dest.ortho.nodes)) else ortho.nodes
    pathway.nodes.filtered <- apply.hide.filters(pathway.nodes.info,"node",hide.elements)
    pathway.nodes.info <- pathway.nodes.info[pathway.nodes.info$node %in% pathway.nodes.filtered$node | pathway.nodes.info$node %in% always.allowed,]
    pathway.nodes.info$layer <- rep(net, nrow(pathway.nodes.info))
    
    #Retrieve pathway edges
    pathway.edges.info <- metapathway.info[metapathway.info$source %in% pathway.nodes.info$node & metapathway.info$target %in% pathway.nodes.info$node,]
    node.roles <- NULL
    if(view.mode=="paths") {
      # "All paths": tutti i nodi/archi che giacciono su almeno un cammino
      # diretto da un nodo sorgente a un nodo destinazione (fino a
      # max.length, se specificato), nella rete GIA' RIDOTTA DAI FILTRI
      # HIDE (stessa convenzione della modalita' "neighbors": se ad esempio
      # i miRNA sono nascosti, la ricerca dei cammini avviene nella rete
      # senza miRNA, salvo le sorgenti/destinazioni stesse).
      if(length(source.genes)>0 && length(dest.genes)>0) {
        #Nuova opzione: mostra solo i cammini che passano per ALMENO UN
        #nodo perturbato (score diverso da zero) - indipendente dal limite
        #sulla lunghezza dei cammini (max.length), i due filtri si
        #combinano.
        perturbed.nodes <- if(isTRUE(only.perturbed.paths)) all.perturbed.nodes else NULL
        path.result <- get.paths.between.edges(pathway.edges.info, source.ortho.nodes, dest.ortho.nodes,
                                               max.length = max.length, perturbed.nodes = perturbed.nodes)
        pathway.edges.info <- path.result$edges
        node.roles <- path.result$roles
      } else {
        pathway.edges.info <- pathway.edges.info[0,]
      }
    } else if(!"All" %in% ref.genes) {
      # "Ego-network": nodi entro max.hops passi (in entrambe le direzioni)
      # dal/dai nodo/i di riferimento selezionato/i.
      ego.result <- get.ego.network.edges(pathway.edges.info, ortho.nodes, max.hops = max.hops)
      pathway.edges.info <- ego.result$edges
      node.roles <- ego.result$roles
    }
    if(nrow(pathway.edges.info)>0) {
      pathway.nodes.ids <- unique(c(pathway.edges.info$source,pathway.edges.info$target))
      if(view.mode!="paths" && !"All" %in% ref.genes) pathway.nodes.ids <- unique(c(pathway.nodes.ids, ortho.nodes))
      pathway.nodes.info <- unique(pathway.nodes.info[pathway.nodes.info$node %in% pathway.nodes.ids,])
    } else if(view.mode=="paths") {
      #Nessun cammino trovato tra le sorgenti e le destinazioni selezionate.
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
    
    #Node size: valore fisso di default (non piu' variabile in base alla
    #distanza/"hop", concetto rimosso insieme a "Max hops" con la nuova
    #logica sorgente->destinazione dei cammini).
    multilayer.nodes$size <- 25
    
    #Set tooltip string for nodes
    multilayer.nodes$title <- paste0("Score: ",round(multilayer.nodes$activity,3))
    
    #Evidenzia i nodi sorgente/destinazione dei cammini (modalita' "All
    #paths") o di riferimento (modalita' "Ego-network"): un bordo distinto
    #e piu' spesso, che sovrascrive - solo per questi nodi specifici - il
    #colore usato altrove per identificare il layer/organismo. L'
    #informazione sul layer resta comunque visibile su tutti gli altri
    #nodi del grafico, quindi non se ne perde la lettura complessiva.
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

#-----------------------------------------------------------------------------
# Elenco dei nodi selezionabili (per la ricerca "Search by Node"/"Source
# node"/"Destination node") per ciascun network/layer, dato un insieme di
# pathway. Include SEMPRE tutti i tipi di nodo (geni, miRNA, farmaci,
# composti chimici): a differenza di prima, i filtri "Show ..." si
# applicano solo alla VISUALIZZAZIONE del grafico finale, non a cosa e'
# ricercabile/selezionabile come nodo di partenza.
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Abbrevia un'etichetta troppo lunga con puntini di sospensione, per la
# visualizzazione nei menu a tendina. Alcuni nomi di composti chimici sono
# estremamente lunghi e rendevano il menu molto largo indipendentemente
# dalle regole CSS di larghezza massima/troncamento (probabilmente perche'
# bootstrap-select calcola la propria larghezza "auto" misurando il testo
# per intero) - troncare l'etichetta VISUALIZZATA direttamente qui, lato
# dati, risolve il problema a monte, senza dipendere da alcun
# comportamento CSS/JS difficile da verificare senza un browser reale.
# Il VALORE sottostante (usato per la selezione/il confronto con
# pathway.list$nodeName) NON viene mai troncato: solo l'etichetta mostrata.
#-----------------------------------------------------------------------------
truncate.label <- function(x, max.chars = 40)
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
  #Etichetta abbreviata per i nomi di pathway molto lunghi (stesso motivo
  #e stessa tecnica usata per i nomi di nodo in get.list.selectable.nodes:
  #il VALORE selezionabile resta il nome completo, solo il testo mostrato
  #nel menu viene troncato).
  names(list.options) <- truncate.label(list.options)
  return(list.options)
}

#-----------------------------------------------------------------------------
# Restringe l'elenco completo dei nodi selezionabili (starting.list.nodes,
# una lista per network/layer) ai soli network attualmente selezionati -
# non filtra piu' per tipo di nodo (vedi nota sopra su
# get.list.selectable.nodes).
#-----------------------------------------------------------------------------
filter.selectable.nodes <- function(starting.list.nodes,list.networks)
{
  starting.list.nodes[list.networks]
}

#-----------------------------------------------------------------------------
# Costruisce la tabella "larga" nodo x network per la scheda "Node
# comparison": una riga per ciascun nodeName distinto tra i network
# caricati, una colonna per network con il relativo score, e una colonna
# "delta" (max-min tra i valori NON mancanti di quella riga - NA se un
# network non ha quel gene/nodo per il proprio organismo, es. perche' e'
# di un organismo diverso). Ordinata per delta decrescente (i nodi piu'
# "diversi" tra i file per primi). Tutti i nodi vengono inclusi, non solo
# quelli con score diverso da zero.
#-----------------------------------------------------------------------------
build.node.comparison.table <- function(data.list, pathway.list, pathways = NULL)
{
  empty.result <- data.frame(nodeName=character(0), delta=numeric(0), pathwaysList=character(0), stringsAsFactors = FALSE)
  if(length(data.list) == 0) return(empty.result)
  
  organisms.present <- unique(as.character(sapply(data.list, function(d) d$organism)))
  organisms.present <- organisms.present[!is.na(organisms.present)]
  organisms.present <- intersect(organisms.present, names(pathway.list))
  
  #Elenco COMPLETO (mai filtrato dal parametro "pathways") delle pathway a
  #cui ciascun nodeName appartiene, aggregato su tutti gli organismi
  #attualmente rilevanti - usato per il tooltip "a quali pathway
  #appartiene questo nodo" al passaggio del mouse sul nome, che deve
  #restare informativo indipendentemente dal filtro pathway attivo.
  pathway.membership.agg <- data.frame(nodeName=character(0), pathwaysList=character(0), stringsAsFactors = FALSE)
  if(length(organisms.present) > 0) {
    pathway.membership <- do.call(rbind, lapply(organisms.present, function(org) {
      unique(pathway.list[[org]][,c("nodeName","pathwayName")])
    }))
    pathway.membership.agg <- aggregate(pathwayName ~ nodeName, data = pathway.membership,
                                        FUN = function(x) paste(sort(unique(x)), collapse = "|||"))
    colnames(pathway.membership.agg)[2] <- "pathwaysList"
  }
  
  per.network <- lapply(names(data.list), function(net) {
    organism <- data.list[[net]]$organism
    if(is.na(organism) || !(organism %in% names(pathway.list))) return(NULL)
    pw.data <- pathway.list[[organism]]
    #Filtro opzionale per pathway selezionate: un nodo resta nella tabella
    #se appartiene ad ALMENO UNA delle pathway scelte (unione, non
    #intersezione) - NULL o vettore vuoto significa "nessun filtro",
    #comportamento identico a prima.
    if(!is.null(pathways) && length(pathways) > 0) {
      pw.data <- pw.data[pw.data$pathwayName %in% pathways, ]
    }
    node.names <- unique(pw.data[,c("node","nodeName")])
    if(nrow(node.names) == 0) return(NULL)
    merged <- merge(node.names, data.list[[net]]$data, by="node", all.x = FALSE)
    if(nrow(merged) == 0) return(NULL)
    merged <- merged[!duplicated(merged$nodeName), c("nodeName","activity")]
    colnames(merged)[2] <- net
    merged
  })
  per.network <- Filter(Negate(is.null), per.network)
  if(length(per.network) == 0) return(empty.result)
  
  wide <- Reduce(function(a,b) merge(a, b, by = "nodeName", all = TRUE), per.network)
  
  score.cols <- setdiff(colnames(wide), "nodeName")
  wide$delta <- apply(wide[,score.cols,drop=FALSE], 1, function(row) {
    vals <- row[!is.na(row)]
    if(length(vals) == 0) return(NA_real_)
    max(vals) - min(vals)
  })
  
  wide <- merge(wide, pathway.membership.agg, by = "nodeName", all.x = TRUE)
  wide$pathwaysList[is.na(wide$pathwaysList)] <- ""
  
  wide <- wide[order(-wide$delta), ]
  rownames(wide) <- NULL
  wide
}

#-----------------------------------------------------------------------------
# Massimo assoluto (simmetrico, arrotondato per eccesso) su cui basare la
# scala colore della tabella "Node comparison" - stessa identica logica
# gia' usata per colorare i nodi nel grafico della rete (vedi
# build.pathway.net -> "Set node colors"), ma calcolata su TUTTI i valori
# della tabella nodo x network (non sul solo sottoinsieme mostrato in una
# pagina), cosi' che uno stesso valore abbia sempre lo stesso colore in
# qualunque pagina/ricerca ci si trovi.
#-----------------------------------------------------------------------------
compute.abs.max <- function(node.comparison.table)
{
  score.cols <- setdiff(colnames(node.comparison.table), c("nodeName","delta","pathwaysList"))
  if(length(score.cols) == 0) return(1)
  all.vals <- unlist(node.comparison.table[,score.cols,drop=FALSE])
  all.vals <- all.vals[!is.na(all.vals)]
  if(length(all.vals) == 0 || all(all.vals == 0)) return(1)
  ceiling(max(abs(min(all.vals)), abs(max(all.vals))))
}

metapathway.list <- list()
pathway.list <- list()
ortho.list <- list()
map.organism <- read.table("Data/mapOrganism.txt",header=T,sep="\t")
