#R code to generate the Sankey graph in Figure_segmental_connections Panel D in the Verasztó et al 2021 Platynereis connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely Feb 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#load RColorBrewer for color palette 
library(RColorBrewer)

#set working directory
setwd('/working_dir/')

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info

#retrieve all annotations for the cells in the connectome
connectome_cells <- catmaid_get_annotations_for_skeletons("annotation:^connectome$", pid = 11)


#define the six body regions, matching the catmaid annotations
body_regions <- c('episphere','segment_0', 'segment_1', 'segment_2', 'segment_3', 'pygidium')
#define the cell categories, matching the catmaid annotations
cell_category <- c('Sensory neuron','interneuron', 'motorneuron', 'effector')

list_position=0; synapse_list=list(); matrix_lists=list()
#these iterated loops will query the connectome_cells data based on annotations and retrieve the connectivity between sets of cells (defined by their skids)
for (i in c(1:6)){  #iterate through the six body regions
  body_region_i = body_regions[i]
  cycle=0
  for (j in c(1:4)){  #iterate through the cell categories
    cell_category_j = cell_category[j]
    print (body_regions[i]); print (cell_category[j])
    presyn_skids_list1 <- list(); presyn_skids_list2 <- list(); counter_k1=0; counter_k2=0
    
      for (k in c(1:length(connectome_cells$annotation))){
        if (connectome_cells$annotation[k] == body_region_i)
        {#here we collect the presyn skids that match the first annotation in the nested loops
          counter_k1 <- counter_k1+1; presyn_skids_list1[[counter_k1]] <- connectome_cells$skid[k]}
        if (connectome_cells$annotation[k] == cell_category_j)
        {#here we collect the presyn skids that match the second annotation in the nested loops
          counter_k2 <- counter_k2+1; presyn_skids_list2[[counter_k2]] <- connectome_cells$skid[k]}
        #find skids that are shared between the two lists
        }
      presyn_skids <- intersect(presyn_skids_list1,presyn_skids_list2)

      #we do these nested iterations to get the postsyn skids per category
      for (l in c(1:6)){  #iterate through the six body regions
        body_region_l = body_regions[l]
        for (m in c(1:4)){  #iterate through the cell categories
          cell_category_m = cell_category[m]
          postsyn_skids_list1 <- list(); postsyn_skids_list2 <- list(); counter_n1=0; counter_n2=0  
          
          for (n in c(1:length(connectome_cells$annotation))){
            if (connectome_cells$annotation[n] == body_region_l)
            {#here we collect the presyn skids that match the first annotation in the nested loops
              counter_n1 <- counter_n1+1; postsyn_skids_list1[[counter_n1]] <- connectome_cells$skid[n]}
            if (connectome_cells$annotation[n] == cell_category_m)
            {#here we collect the presyn skids that match the second annotation in the nested loops
              counter_n2 <- counter_n2+1; postsyn_skids_list2[[counter_n2]] <- connectome_cells$skid[n]}
            #find skids that are shared between the two lists
          }
          postsyn_skids <- intersect(postsyn_skids_list1,postsyn_skids_list2)
          cellgroup_conn=list()
          # get connectors between neurons of interest
          if (length(presyn_skids)==0 | length(postsyn_skids)==0) {
            N_synapses <- 0
          } else {
            cellgroup_conn = catmaid_get_connectors_between(pre=presyn_skids, post=postsyn_skids, pid=11)
            N_synapses=nrow(cellgroup_conn)
            if(length(cellgroup_conn) == 0) {N_synapses <- 0}
                }
          list_position <- list_position + 1; print (list_position); print(N_synapses)
          synapse_list [[list_position]] <- N_synapses
          
          }
      }
    }
}

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=24 )

#define the names for plotting of the six body regions and cell categories
body_regions_name <- c('head','sg0', 'sg1', 'sg2', 'sg3', 'pyg')
#define the cell categories, matching the catmaid annotations
cell_category_name <- c('SN','IN', 'MN', 'eff')

#we make a cellgroup name list 
cell_group_names=list();counter=0
for (i in c(1:6)){  #iterate through the six body regions
  body_region_i = body_regions_name[i]
  for (j in c(1:4)){  #iterate through the cell categories
    counter <- counter+1
    cell_group_names[[counter]] <- paste(cell_category_name[j], body_regions_name[i], sep="-")
}}  

synapse_matrix = as.data.frame(synapse_matrix)

#assign column names to matrix
synapse_matrix=setNames(synapse_matrix, as.character(cell_group_names))

#assign row names to matrix
rownames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix = as.matrix(synapse_matrix)
synapse_matrix

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
    weighted = T,  diag = TRUE, add.colnames = NULL, add.rownames = NA)

wc <- cluster_walktrap(celltype_conn_graph)
members <- membership(wc)

# Convert to object suitable for networkD3
celltype_conn_graph_d3 <- igraph_to_networkD3(celltype_conn_graph, group = members)

#The NodeGroup vector in the Nodes data frame needs to be non-numeric so we convert it to character
celltype_conn_graph_d3$nodes$group <- as.character(celltype_conn_graph_d3$nodes$group)

# Plot
sankeyNetwork(Links = celltype_conn_graph_d3$links, Nodes = celltype_conn_graph_d3$nodes, Source = "source",
                    Target = "target",  NodeID = "name",Value = "value",
                    LinkGroup = NULL, units = "", NodeGroup = "group",
                    colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"), fontSize = 66,
                    fontFamily = "sans", nodeWidth = 200, nodePadding = 10, margin = 0,
                    height = NULL, width = NULL, iterations = 500, sinksRight = T)

#to save as pdf open in the browser (Show in new window icon in the Viewer window) and save as pdf

#filtered by edge strength (optional)
link_strength <- celltype_conn_graph_d3_thr$links$value
hist(link_strength, col='blue',breaks=100)

syn_threshold=25
celltype_conn_graph_d3_thr <- celltype_conn_graph_d3

for (i in seq_along(link_strength)) {
  if (link_strength[i]<syn_threshold){
        link_strength[i] <- 0}
}
link_strength
hist(link_strength, col='blue',breaks=100)

celltype_conn_graph_d3_thr$links$value <- link_strength

# Plot
sankeyNetwork(Links = celltype_conn_graph_d3_thr$links, Nodes = celltype_conn_graph_d3_thr$nodes, Source = "source",
              Target = "target",  NodeID = "name",Value = "value",
              LinkGroup = NULL, units = "", NodeGroup = "group",
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"), fontSize = 16,
              fontFamily = "sans", nodeWidth = 30, nodePadding = 5, margin = NULL,
              height = NULL, width = NULL, iterations = 300, sinksRight = T)

