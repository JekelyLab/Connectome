#This code was used to compare the connectivity of cell types on the left and right body sides of the 3 day old Platynereis larva described in Veraszto et al. 2021
#Gaspar Jekely 2021 Feb

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#set working directory
setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Figures/Figure-left-right-connectivity')

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


#define the six body regions, matching the catmaid annotations
body_sides <- c('left_side','right_side')

#define empty synapse list with the right dimensions
synapse_list <- vector("list", length(annotation_celltypelist)*length(body_sides)*length(annotation_celltypelist))

#define empty cells_per_celltype list with the right dimensions
cells_per_celltype <- vector("list", length(annotation_celltypelist))

#we iterate through the celltype lists and body regions and retrieve connectivity between pre- and postsyn sets of skids
list_position=0
for (i in c(1:2)){  #iterate through the two body sides
  body_side_i = body_sides[i]
  cycle=0
  print (body_sides[i])
  for (df1 in annotation_celltypelist){    #iterate through the celltype list  (all presyn cells per body region)
    presyn_skids <- df1[df1$annotation == body_side_i,1]
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
print(list_position)
}

#convert synapse list into a matrix of appropriate dimensions
length(synapse_list)/(182*2)
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=182*length(body_sides) )

dim(synapse_matrix)

#we make a celltype name list 
id_pre_left <- c(paste('left',(1:(nrow(synapse_matrix)/2)),sep='_'))
id_pre_right <- c(paste('right',(1:(nrow(synapse_matrix)/2)),sep='_'))
id_pre <- c(id_pre_left,id_pre_right)

presyn_name <- paste('celltype',id_pre,sep="_")
presyn_name
id_post_neu <- c(paste('neuronal',(1:(nrow(synapse_matrix)/2)),sep='_'))
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
pdf('Left_right_celltype_connectivity_matrix.pdf', width=7.5, height=7)
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
              ylab = "presynaptic neuronal celltypes, left side 1-182, right side 1-182",
              main = "left-right connectivity matrix of celltypes (sqrt of summed synapses)",
              
    )
dev.off()

#plot the same with ggplot
synapse_matrix
#plot as dot plot
synapse_matrix_sqr_tb <- as.data.frame((sqrt(synapse_matrix))) %>%
  rownames_to_column(var = "presyn")%>%
  pivot_longer(-presyn, names_to = "postsyn", values_to = "synapses")%>%
  group_by(postsyn)%>%
  mutate(percent_presyn =  synapses / sum(synapses, na.rm = TRUE))



#alternative visualisation with ggplot (not used in figure)  
library(ggplot2)
ggplot(synapse_matrix_sqr_tb) +
  aes(x = postsyn, y = presyn, colour = percent_presyn) +
  geom_point(size=0.3) +
  scale_color_gradient2(low = "white",
                        mid = "#56B4E9",
                        high = "#D55E00",
                        midpoint = 0.5,
                       na.value = "white",
                       guide = "colourbar") +
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x='postsynaptic cell types', y='presynaptic cell types')


##########################################
#correlation calculations
##########################################


#calculate pearson between top half (left as presyn) and bottom half (right as presyn) of the synaptic matrix
dim(synapse_matrix)
synapse_matrix[1:(nrow(synapse_matrix)/2),]
synapse_matrix[((nrow(synapse_matrix)/2)+1):(nrow(synapse_matrix)),]
(nrow(synapse_matrix)/2)

pearson_left_right <- cor.test(synapse_matrix[1:(nrow(synapse_matrix)/2),]
,synapse_matrix[((nrow(synapse_matrix)/2)+1):(nrow(synapse_matrix)),], 
method = "pearson")

pearson_left_right


#calculate pearson correlation for each left-right row pair (presyn pairs)
pearson <- matrix(nrow=(nrow(synapse_matrix)/2), ncol=(nrow(synapse_matrix)/2))

for (i in c(1:(nrow(synapse_matrix)/2))){
  for (j in c((nrow(synapse_matrix)/2)+1):nrow(synapse_matrix)){
  #calculate correlation, use $estimate to extract the correlation value from the cor.test output
  pearson[i,(j-nrow(synapse_matrix)/2)] <- cor.test(synapse_matrix[i,],synapse_matrix[j,], method = "pearson")$estimate
  }
}

dim(pearson)
pearson[1:10,1:10]
#plot heatmap with heatmaply
library(heatmaply)

heatmaply_cor(pearson,
              column_text_angle = 90,row_text_angle = 0,
              fontsize_row = 4,
              fontsize_col = 4,
              Rowv=NA, Colv=NA,
              show_dendrogram = c(F, F),
              hide_colorbar = F,
              plot_method = c("ggplot"),
              revC=F,
              na.value='grey0',
              col=c('black','grey20','white','#1098CD','red'),
              #file = "Left_right_celltype_connectivity_pearson.pdf"
              )
           

######################################################   
#divide all rows by the number of presyn cells in that celltype
#calculate pearson
cells_per_celltype[1:182]
cells_per_celltype_left_right=c(cells_per_celltype[1:182],cells_per_celltype[1:182])
cells_per_celltype_mat <- matrix(nrow=364, ncol=272)              
for (i in c(1:272)){cells_per_celltype_mat[,i] <- unlist(cells_per_celltype[1:182])}
cells_per_celltype_mat[1:13,1:3]

synapse_matrix_div <-  synapse_matrix/cells_per_celltype_mat
dim(synapse_matrix)
dim(cells_per_celltype_mat)
synapse_matrix_div[1:10,1:10]
synapse_matrix_div[is.na(synapse_matrix_div)] <- 0
synapse_matrix_div[is.infinite(synapse_matrix_div)] <- 0

#heatmap of divided synapse matrix
heatmaply(sqrt(synapse_matrix_div),
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
          col=c('grey20','cyan','orange','orange','orange','red','red','red'),
          #col=hcl.colors(20, "Oslo", alpha = 1, rev = F, fixup = F),
          xlab = "postsynaptic celltypes, neuronal 1-182 and non-neuronal 1-90",
          ylab = "presynaptic neuronal celltypes, left side 1-182, right side 1-182",
          main = "left-right connectivity matrix of celltypes (sqrt of summed synapses)",
          
)



#calculate pearson between top half (left as presyn) and bottom half (right as presyn) of the synaptic matrix

pearson_left_right_div <- cor.test(synapse_matrix_div[1:(nrow(synapse_matrix_div)/2),]
                               ,synapse_matrix_div[((nrow(synapse_matrix_div)/2)+1):(nrow(synapse_matrix_div)),], 
                               method = "pearson")

pearson_left_right_div


#calculate pearson correlation for each left-right row pair (presyn pairs)
pearson_div <- matrix(nrow=(nrow(synapse_matrix_div)/2), ncol=(nrow(synapse_matrix_div)/2))

for (i in c(1:(nrow(synapse_matrix_div)/2))){
  for (j in c((nrow(synapse_matrix_div)/2)+1):nrow(synapse_matrix_div)){
    #calculate correlation, use $estimate to extract the correlation value from the cor.test output
    pearson_div[i,(j-nrow(synapse_matrix_div)/2)] <- cor.test(synapse_matrix_div[i,],synapse_matrix_div[j,], method = "pearson")$estimate
  }
}

dim(pearson_div)
pearson_div[1:10,1:10]
#plot heatmap with heatmaply
library(heatmaply)

heatmaply_cor(pearson_div,
              column_text_angle = 90,row_text_angle = 0,
              fontsize_row = 4,
              fontsize_col = 4,
              Rowv=NA, Colv=NA,
              show_dendrogram = c(F, F),
              hide_colorbar = F,
              plot_method = c("ggplot"),
              revC=F,
              na.value='grey0',
              col=c('black','grey20','white','#1098CD','red'),
              #file = "Left_right_celltype_connectivity_pearson.pdf"
)


              