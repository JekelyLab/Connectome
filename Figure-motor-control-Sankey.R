#R code to generate the Sankey graph in Figure_motor_control
#in the Verasztó et al 2021 Platynereis connectome paper
#Uses a graph file generated in Catmaid as input
#Gaspar Jekely Feb 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

#set working directory
setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Figures/Figure-motor-control/')

IN_MN_connectivity <- read.csv('adjacency_matrix.csv', header = T, row.names = 1)

dim(IN_MN_connectivity)
IN_MN_connectivity <- as.matrix(IN_MN_connectivity)
###############################
#graph visualisation
##############################

# Load igraph
library(igraph)
#https://rdrr.io/cran/igraph/man/
library(networkD3)
library(htmlwidgets)

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
celltype_conn_graph <- graph_from_adjacency_matrix(IN_MN_connectivity,
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
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"), fontSize = 12,
              fontFamily = "sans", nodeWidth = 20, nodePadding = 10, margin = 0,
              height = NULL, width = NULL, iterations = 500, sinksRight = T)

#to save as pdf open in the browser (Show in new window icon in the Viewer window) and save as pdf

library(heatmaply)

heatmaply(IN_MN_connectivity[1:18,],
          colors = hcl.colors(40, palette='Oslo'),
          column_text_angle = 90,row_text_angle = 0,
          Rowv=T,Colv=T,
          dist_method="euclidean",
          hclustfun = hclust,
          hclust_method='ward.D2' ,
          dendrogram = c("both"),
          show_dendrogram = c(T,T),
          grid_gap = 1,
          hide_colorbar = F,
          plot_method = c("plotly"),
          fontsize_row = 8,
          fontsize_col = 8,
          revC=F
)


#filtered by edge strength (optional)
link_strength <- celltype_conn_graph_d3$links$value
hist(link_strength, col='blue',breaks=100)

syn_threshold=5
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
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"), fontSize = 12,
              fontFamily = "sans", nodeWidth = 20, nodePadding = 5, margin = NULL,
              height = NULL, width = NULL, iterations = 300, sinksRight = T)

