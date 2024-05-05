#organism <- "Human"
#organism <- "Mouse"
organism <- "Rat"
gene <- "ISG15"

#INPUT DATA
#input.data <- as.data.frame(fread(paste0("Data/Input/Phensim/",gene,"-up-",tolower(organism),".tsv"),sep="\t",header=T))
input.data <- as.data.frame(fread(paste0("C:/Ricerca/Phensim/ISG15_up/",organism,"/PhensimSimulation.txt"),sep="\t",header=T))
organism <- map.organism[map.organism$Code==tolower(unlist(strsplit(input.data[1,1],":|-"))[2]),"Organism"]
input.data <- as.data.frame(unique(input.data[,c("Node Id","Activity Score")]))
input.data$Organism <- organism
input.data <- input.data[,c(3,1,2)]
colnames(input.data) <- c("Organism","Node","Score")
write.table(input.data,paste0("Data/Input/Custom/",gene,"-up-",tolower(organism),".txt"),sep="\t",row.names = F,quote=F)


#METAPATHWAY
meta.data <- as.data.frame(fread(paste0("C:/Ricerca/Phensim/ISG15_up/",organism,"/PhensimMetapathway.txt"),sep="\t",header=T))
meta.data <- meta.data[,c("source","source_name","target","target_name","edge_weight")]
colnames(meta.data) <- c("source","sourceName","target","targetName","weight")
#write.table(meta.data,paste0("Data/Metapathways/",organism,".txt"),row.names = F,sep="\t",quote=F)
saveRDS(meta.data,paste0("Data/Metapathways/",organism,".rds"))

#PATHWAYS
path.data <- as.data.frame(fread(paste0("C:/Ricerca/Phensim/ISG15_up/",organism,"/PhensimSimulation.txt"),sep="\t",header=T))
path.data <- path.data[,c("Pathway Id","Pathway Name","Node Id","Node Name","Is Endpoint")]
colnames(path.data) <- c("pathway","pathwayName","node","nodeName","endpoint")
path.data$endpoint[path.data$endpoint=="Yes"] <- T
path.data$endpoint[path.data$endpoint=="No"] <- F
path.data$endpoint <- as.logical(path.data$endpoint)
saveRDS(path.data,paste0("Data/Pathways/",organism,".rds"))

#ORTHOLOGS - SAME ORGANISM
organism <- "Rat"
homo.data <- as.data.frame(fread(paste0("C:/Ricerca/Phensim/ISG15_up/",organism,"/PhensimSimulation.txt"),sep="\t",header=T))
homo.data <- unique(homo.data[,c("Node Name","Node Id")])
homo.data <- homo.data[!startsWith(homo.data$`Node Id`,"chebi:"),]
homo.data <- homo.data[!startsWith(homo.data$`Node Id`,"cpd:"),]
homo.data <- homo.data[!startsWith(homo.data$`Node Id`,"gl:"),]
homo.data <- homo.data[!startsWith(homo.data$`Node Id`,"dr:"),]
colnames(homo.data) <- c(paste0(organism," symbol"),paste0(organism," id"))
ortho.data <- readRDS("Data/Orthologs/Rat-Rat.rds")
sub.data <- ortho.data[,c(paste0(organism," symbol"),paste0(organism," id"))]
homo.data <- unique(rbind(sub.data,homo.data))
saveRDS(homo.data,paste0("Data/Orthologs/",organism,"-",organism,".rds"))

#ORTHOLOGS - DIFFERENT ORGANISM

organism <- "Rat"

ortho.data <- read.table("C:/Ricerca/Bio/Orthologs/RatGenome/ORTHOLOGS_RAT.txt",sep="\t",header=T,skip=15,
                         blank.lines.skip = F)
sub.data <- ortho.data[,c("RAT_GENE_SYMBOL","RAT_GENE_NCBI_GENE_ID",
                          paste0(toupper(organism),"_GENE_SYMBOL"),paste0(toupper(organism),"_GENE_NCBI_GENE_ID"))]
sub.data <- sub.data[!is.na(sub.data[,paste0(toupper(organism),"_GENE_NCBI_GENE_ID")]),]
if(organism<"Rat") {
  sub.data <- sub.data[,c(3,4,1,2)]
  colnames(sub.data) <- c(paste0(organism, " symbol"),paste0(organism," id"),"Rat symbol","Rat id")
  saveRDS(sub.data,paste0("Data/Orthologs/",organism,"-Rat.rds"))
} else {
  colnames(sub.data) <- c("Rat symbol","Rat id",paste0(organism, " symbol"),paste0(organism," id"))
  saveRDS(sub.data,paste0("Data/Orthologs/Rat-",organism,".rds"))
}

ortho.data <- readRDS("Data/Orthologs/Human-Mouse.rds")
#saveRDS(ortho.data,"Data/Orthologs/Human-Mouse.rds")
sub.data <- ortho.data[,c("Human symbol","Human id")]
sub.data <- unique(sub.data)
saveRDS(sub.data,"Data/Orthologs/Human-Human.rds")
sub.data <- ortho.data[,c("Mouse symbol","Mouse id")]
sub.data <- unique(sub.data)
saveRDS(sub.data,"Data/Orthologs/Mouse-Mouse.rds")
