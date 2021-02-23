#This code was used to generate the cell type connectivity matrix of the 3 day old Platynereis larva 
#described in Veraszto et al. 2021
#Gaspar Jekely 2021 Feb

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#set working directory
setwd('/Users/gaspar/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Figures/Figure_celltypes')

# catmaid connection, needs username, password AND token
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
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
  }

#we read all non-neuronal celltypes from 1-90 and all annotations
for (i in c(1:90)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i+number_of_neuron_types]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

#define empty synapse list with the right dimensions
synapse_list <- vector("list", length(annotation_celltypelist)*length(annotation_celltypelist))

#define empty cells_per_celltype list with the right dimensions
cells_per_celltype <- vector("list", length(annotation_celltypelist))

#we iterate through the celltype lists and body regions and retrieve connectivity between pre- and postsyn sets of skids
list_position=0
cycle=0
for (df1 in annotation_celltypelist){    #iterate through the celltype list  (all presyn cells in groups)
    presyn_skids <- df1$skid
    cycle <- cycle + 1
    print (cycle)
    cells_per_celltype[cycle] <- length(presyn_skids)
    if (cycle>182){next} #we only collect the connectivity for presyn cells 1-182 (these are the neurons)
    
        for (df2 in annotation_celltypelist){  #nested iteration through the celltype list (postsyn cells)
          postsyn_skids <- df2$skid
        assign("cellgroup_conn", NULL, envir = .GlobalEnv)  #in every iteration we empty the connectivity list
        # get connectors betwwen neurons of interest
        if (length(presyn_skids)==0 | length(postsyn_skids)==0) {
          N_synapses <- 0
      } else {
          cellgroup_conn = catmaid_get_connectors_between(pre=presyn_skids, post=postsyn_skids, pid=11)
          N_synapses=nrow(cellgroup_conn)
         if(length(cellgroup_conn) == 0) {N_synapses <- 0}
        }
        list_position <- list_position + 1
        synapse_list [[list_position]] <- N_synapses
        }
}



#convert synapse list into a matrix of appropriate dimensions

synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=182 )

dim(synapse_matrix)

#we make a celltype name list 
id_pre_left <- c(paste('left',(1:(nrow(synapse_matrix)/2)),sep='_'))
id_pre_right <- c(paste('right',(1:(nrow(synapse_matrix)/2)),sep='_'))
id_pre <- c(id_pre_left,id_pre_right)

presyn_name <- paste('celltype',1:(nrow(synapse_matrix)),sep="_")
presyn_name
id_post_neu <- c(paste('neuronal',(1:(nrow(synapse_matrix))),sep='_'))
id_post_non_neu <- c(paste('non_neuronal',(1:90),sep='_'))
id_post <- c(id_post_neu,id_post_non_neu)
postsyn_name <- paste('celltype',id_post,sep="_")
postsyn_name

synapse_matrix = as.data.frame(synapse_matrix)
#assign column names to matrix
synapse_matrix=setNames(synapse_matrix, as.character(postsyn_name))
#assign row names to matrix
rownames(synapse_matrix) <- as.character(presyn_name)
synapse_matrix = as.matrix(synapse_matrix)


#plot heatmap with heatmaply
library(heatmaply)
library(RColorBrewer)

#save the matrix as pdf
pdf('Celltype_connectivity_matrix.pdf', width=7.5, height=7)
heatmaply(sqrt(synapse_matrix),
              column_text_angle = 90,row_text_angle = 0,
              fontsize_row = 3,
              fontsize_col = 3,
              Rowv=F, Colv=F,
              show_dendrogram = c(F, F),
              grid_gap = 0,
              hide_colorbar = F,
              plot_method = c("ggplot"),
              revC=F,
             #col=brewer.pal(9, 'Spectral'),
              col=c('grey20','cyan','orange'),
            #col=hcl.colors(20, "Oslo", alpha = 1, rev = F, fixup = F),
              xlab = "postsynaptic celltypes, neuronal 1-182 and non-neuronal 1-90",
              ylab = "presynaptic neuronal celltypes 1-182",
              main = "connectivity matrix of celltypes (sqrt of summed synapses)",
              
    )
dev.off()

write.csv(synapse_matrix, file = "Left_right_celltype_connectivity_matrix.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")
