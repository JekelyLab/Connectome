rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#load RColorBrewer for color palette 
library(RColorBrewer)

#set working directory
setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Figures/Figure_celltypes/')

# catmaid connection, needs username, password AND token - weird!
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info


annotation_celltypelist = list()
number_of_neuron_types=0

#first we read all celltypes from 1-182 and all annotations
for (i in c(1:182)){
  number_of_neuron_types <- number_of_neuron_types + 1
  annotation = paste("annotation:^celltype", i, "$", sep="")
    #read presyn neuron group by annotation, also retrieve all annotation
  #celltypelist[[i]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = FALSE) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
  }
#we read all non-neuronal celltypes from 1-90 and all annotations
for (i in c(1:90)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #read presyn neuron group by annotation, also retrieve all annotation
  #celltypelist[[i+number_of_neuron_types]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = FALSE) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i+number_of_neuron_types]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

length(annotation_celltypelist)

#define empty synapse list with the right dimensions
synapse_list <- vector("list", length(annotation_celltypelist)*length(annotation_celltypelist))
#we iterate through the celltype lists and body regions and retrieve connectivity between pre- and postsyn sets of skids

#retrieve unique skids for a celltype from the annotation_celltype list
unique(annotation_celltypelist[[1]]$skid)

list_position=0;cycle=0
#this can take a while...
for (celltype_skids_pre in annotation_celltypelist){
  presyn_skids <- unique(celltype_skids_pre$skid)
  for (celltype_skids_post in annotation_celltypelist){
    postsyn_skids <- unique(celltype_skids_post$skid)
    assign("celltype_conn", NULL, envir = .GlobalEnv)  #in every iteration we empty the connectivity list  
    # get connectors between neurons of interest
    celltype_conn = catmaid_get_connectors_between(pre=presyn_skids, post=postsyn_skids, pid=11)
      N_synapses=nrow(celltype_conn)
      if(length(celltype_conn) == 0) {N_synapses <- 0}
        list_position <- list_position + 1
      synapse_list[[list_position]] <- N_synapses
    }
  cycle <- cycle + 1;print (cycle)
}

synapse_list[73984] <- 0
#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=length(annotation_celltypelist) )

#retrieve one skid and cell name per cell type
skids1per <- unlist(lapply(annotation_celltypelist,function(x) x$skid[1]))
names1per <- catmaid_get_neuronnames(skids1per, pid=11)

#use regex to truncate names at first _ which should give the name of the cell type
celltype_names=(gsub("(\\_|\\;).+", "", names1per))
names1per
synapse_matrix = as.data.frame(synapse_matrix)

#assign column names to matrix
synapse_matrix=setNames(synapse_matrix, as.character(celltype_names))

#assign row names to matrix
rownames(synapse_matrix) <- as.character(celltype_names)
synapse_matrix = as.matrix(synapse_matrix)

###############################
#graph visualisation
##############################

# Load igraph
library(igraph)
#https://rdrr.io/cran/igraph/man/
library(networkD3)
library(htmlwidgets)

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
celltype_conn_graph <- graph_from_adjacency_matrix(synapse_matrix,
    mode = c("directed"),
    weighted = NULL,  diag = TRUE, add.colnames = NULL, add.rownames = NA)
celltype_conn_graph


wc <- cluster_walktrap(celltype_conn_graph)
members <- membership(wc)

# Convert to object suitable for networkD3
celltype_conn_graph_d3 <- igraph_to_networkD3(celltype_conn_graph, group = members)

#The NodeGroup vector in the Nodes data frame needs to be non-numeric so we convert it to character
celltype_conn_graph_d3$nodes$group <- as.character(celltype_conn_graph_d3$nodes$group)

#need to add a value to each link for it to work
celltype_conn_graph_d3$links$value <- 1

#assign names
celltype_conn_graph_d3$nodes$name <- row.names(synapse_matrix)

# Plot
sankeyNetwork(Links = celltype_conn_graph_d3$links, Nodes = celltype_conn_graph_d3$nodes, Source = "source",
              Target = "target",  NodeID = "name",Value = "value",
              LinkGroup = NULL, units = "", NodeGroup = "group",
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"), fontSize = 22,
              fontFamily = "sans", nodeWidth = 20, nodePadding = 10, margin = 0,
              height = NULL, width = NULL, iterations = 500, sinksRight = F)




####################
#Saving files
####################

#save the matrix as pdf
pdf('Celltype_connectivity_matrix.pdf', width=7.5, height=7)
heatmap((synapse_matrix.no0),   #show matrix 
        Rowv=NA, Colv=NA,
        cexRow = 0.07, cexCol = 0.07, revC=T, 
        scale = 'none',
        col=hcl.colors(300, "Oslo", alpha = 1, rev = FALSE, fixup = TRUE), #col=brewer.pal(9, 'YlOrRd'), #col=terrain.colors(500, rev=T),
        symm=F, margins= c(3,2))
dev.off()


write.csv(synapse_matrix.no0, file = "Celltype_connectivity_matrix.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")
