rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)
#set working directory
setwd('/Users/gaspar/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/statistics/Khulood/')

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")

#read all connectome cells
connectome <- read.neurons.catmaid('^connectome$', pid = 11)

#get the list of skids
skids <- names(connectome)

length(skids)

#retrieve the list of neuron names
connectome_names <- catmaid_get_neuronnames(skids, pid=11)
length(connectome_names)
connectome_names[1:2]

#retrieve the annotations for all cells
annotations <- catmaid_get_annotations_for_skeletons(skids, pid=11)

connections <- catmaid_get_connectors_between(pre_skids = skids, post_skids = skids, pid=11, get_names = TRUE)

# use table() to cross-tabulate number of connections between skids
connectome_conn_table = table(connections[,1:2])

#check if names are unique
which(duplicated(connectome_names))


#create connectivity matrix
conn_mat = matrix(0,nrow=length(connectome),ncol=length(connectome))
# name rows and cols for skids
colnames(conn_mat) = names(connectome)
rownames(conn_mat) = names(connectome)
# populate connectivity matrix based on connectome_conn_table
for (x in rownames(connectome_conn_table)){
  for (y in colnames(connectome_conn_table)){
    conn_mat[x,y] = connectome_conn_table[x,y]
  }
}

conn_mat[1:20,1:10]
dim(conn_mat)


write.csv(conn_mat, file = "Connectome_matrix.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")

write.csv(connectome_names, file = "Connectome_skids_names.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")

write.csv(annotations, file = "Connectome_annotations.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")
