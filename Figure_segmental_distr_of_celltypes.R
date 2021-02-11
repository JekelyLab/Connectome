#R code to generate the graph in Figure_segmental_connections panels H and I in the Verasztó et al 2021 Platynereis connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely Feb 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#set working directory
setwd('/working_dir/')

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info

annotation_neuronal_celltypelist = list()
annotation_non_neuronal_celltypelist = list()

#first we read all celltypes from 1-182 and all annotations
for (i in c(1:182)){
  annotation = paste("annotation:^celltype", i, "$", sep="")
    #read presyn neuron group by annotation, also retrieve all annotation
  #celltypelist[[i]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = FALSE) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
  }

#we read all non-neuronal celltypes from 1-90 and all annotations
for (i in c(1:90)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #read presyn neuron group by annotation, also retrieve all annotation
  #celltypelist[[i+number_of_neuron_types]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = FALSE) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_non_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

#define the six body regions, matching the catmaid annotations
body_regions <- c('episphere','segment_0', 'segment_1', 'segment_2', 'segment_3', 'pygidium')

#define empty matrix of presence-absence of celltype per segment
segmental_neuronal_celltype_complement <- matrix(,nrow=length(body_regions), ncol=length(annotation_neuronal_celltypelist), byrow=T)
segmental_non_neuronal_celltype_complement <- matrix(,nrow=length(body_regions), ncol=length(annotation_non_neuronal_celltypelist), byrow=T)

for (i in c(1:length(annotation_neuronal_celltypelist))){    #iterate through the celltype list
  skids <- unique(annotation_neuronal_celltypelist[[i]]$skid) #retrieve unique skids for a celltype from the annotation_celltype list
    for (j in c(1:length(body_regions))){
      print (body_regions[j])
      if (sum(match(annotation_neuronal_celltypelist[[i]]$annotation, body_regions[j],nomatch = 0))>0){
        print ("match")
        segmental_neuronal_celltype_complement[j,i] <- 1
        #check for sensory neurons etc and change value to use for histogram coloring
        if (sum(match(annotation_neuronal_celltypelist[[i]]$annotation, "interneuron",nomatch = 0))>0){
          segmental_neuronal_celltype_complement[j,i] <- 4}
        if (sum(match(annotation_neuronal_celltypelist[[i]]$annotation, "Sensory neuron",nomatch = 0))>0){
          segmental_neuronal_celltype_complement[j,i] <- 2}
        if (sum(match(annotation_neuronal_celltypelist[[i]]$annotation, "motorneuron",nomatch = 0))>0){
          segmental_neuronal_celltype_complement[j,i] <- 6}
        }
      
      else( segmental_neuronal_celltype_complement[j,i] <- 0)
    }
}


for (i in c(1:length(annotation_non_neuronal_celltypelist))){    #iterate through the celltype list
  skids <- unique(annotation_non_neuronal_celltypelist[[i]]$skid) #retrieve unique skids for a celltype from the annotation_celltype list
  for (j in c(1:length(body_regions))){
    print (body_regions[j])
    if (sum(match(annotation_non_neuronal_celltypelist[[i]]$annotation, body_regions[j],nomatch = 0))>0){
      print ("match")
      segmental_non_neuronal_celltype_complement[j,i] <- 1
      #check for muscle etc and change value to ease histogram coloring
      if (sum(match(annotation_non_neuronal_celltypelist[[i]]$annotation, "muscle",nomatch = 0))>0){
        segmental_non_neuronal_celltype_complement[j,i] <- 3}
      if (sum(match(annotation_non_neuronal_celltypelist[[i]]$annotation, "ciliated cell",nomatch = 0))>0){
        segmental_non_neuronal_celltype_complement[j,i] <- 2}
      if (sum(match(annotation_non_neuronal_celltypelist[[i]]$annotation, "gland cell",nomatch = 0))>0){
        segmental_non_neuronal_celltype_complement[j,i] <- 4}
      if (sum(match(annotation_non_neuronal_celltypelist[[i]]$annotation, "pigment cell",nomatch = 0))>0){
        segmental_non_neuronal_celltype_complement[j,i] <- 5}
    }
    else( segmental_non_neuronal_celltype_complement[j,i] <- 0)
  }
}


#celltype name lists
non_neuronal_celltype_names = c("akrotroch", "crescent cell", "prototroch", "nuchal cilia", "metatroch", "paratroch", "spinGland", "covercell", "ciliatedGland", "eyespot pigment cell", "pigment cell AE", "bright droplets parapodial", "macrophage", "yolk cover cell", "flat glia", "EC rad glia ", "MVGland", "microvillarCell", "protonephridium", "nephridium", "nephridiumTip", "chaeta", "aciculoblast", "circumacicular", "hemichaetal", "ER circumchaetal", "noER circumchaetal", "EC circumchaetal", "HeadGland", "InterparaGland", "spinMicroGland", "CB pigment", "vacuolar cell_head", "Glia pigmented", "pygidial pigment", "meso", "MUSac_notA", "MUSac_notP", "MUSac_notM", "MUSac_neuAV", "MUSac_neuPD", "MUSac_neuPV", "MUSac_neuDy", "MUSac_neuDx", "MUSac_neuDach", "MUSac_neure", "MUSac_i", "MUSob-ant_re", "MUSob-ant_arc", "MUSob-ant_m-pp", "MUSob-ant_ml-pp", "MUSob-ant_l-pp", "MUSob-ant_trans", "MUSob-post_notD", "MUSob-post_neuDlong", "MUSob-post_neuDprox", "MUSob-post_neuDdist", "MUSob-post_neuV", "MUSob-post_notV", "MUSob-postM", "MUSob-post_noty", "MUSob-post_i", "MUSchae_notDob", "MUSchae_notD", "MUSchae_notDn", "MUSchae_notA", "MUSchae_notAac", "MUSchae_notAre", "MUSchae_neuVob", "MUSchae_neuDac", "MUSchae_neuAVo", "MUSchae_neuAVt", "MUSchae_Are", "MUStrans_pyg", "MUSlong_D", "MUSlong_V", "MUSax", "MUSring", "MUSph", "MUSll", "MUSant", "MUSly", "MUSpl", "MUSci", "MUSch", "MUSpx", "MUSpr", "MUStri", "MUSmed_head", "EC")
neuronal_celltype_names =c("PRC", "IN1", "INton", "INpreSer", "cPRC", "INRGWa1", "IN_NOS", "Ser-h1", "MC", "cMNPDF", "cMNATO", "cMNdc", "SNnuch", "SNnuchNS", "SNnuchBx", "SN_LHSO2golden", "SNhorn", "SNhook", "MNant", "SNlasso", "INlasso", "INdecussPre", "INpreMN", "SNMIP-vc", "SNmus", "SNantlerPDF", "SNPDF-dc", "SN_DLSO1.2", "SN_DLSO", "SN_DLSO1.2", "SN_DLSO2", "SN_LHSO2", "eyespot_PRC_R3", "eyespot_PRC_R1", "MS1", "MS4", "MS3", "MS5", "SN_ASTC", "SN_DLSO3_PDF", "SN_DLSO3", "MC2biax", "SN_IRP2-burs", "SN_IRP2-FMRF", "SN_MIP4", "SN_MIP1", "MNspinning", "SN_WLD1", "SN_DSO17", "SNbicil", "SNasym", "SN47Ach", "INCrossbow", "pygPBunp", "INarc1", "INarc2", "INsn1", "INrope", "Loop", "INCM", "MNspider-ant", "MNspider-post", "MNbiramous", "MNring", "MNcrab", "MNhose", "MNbow", "MNwave", "MNacicX", "INchaeMech", "chaeMech", "INFVa_pyg", "INsplitCR ", "INsplitCRATO", "INMC3", "INbiax", "INsplitPB_RF/Ya", "INsplitPBant", "INsplitPB", "ventralParaPB2", "interparaPM1", "cioMNcover", "INsplitPUh", "vMN3", "vMN1", "vMN2", "MC3cover", "antPUc1", "spinPU1", "VentraltrunkPUunp", "ventralParaPU1", "hPUc1", "pygCirrusPU", "hCirrusPUc1", "ventralpostParaPUc2", "hPUc3", "dsoPUc1unp", "pygCirrusPU2", "hPU2l_asymPDF", "hCR2", "VentralpygCR2", "doCRunp", "SN_DLSO1.0", "INdc", "INbackcross3", "MNmouth", "MNob-ipsi", "INbiaxLeg", "INbiaxH", "SN_NS1", "SN_NS3", "SN_NS4", "SN_NS15", "SNaant_3", "SNPDF_pyg", "INfoot", "INW", "INmidL", "INhook", "INlat1", "INlat2", "INSturn", "INcross", "INpro", "INproT2", "INpear", "INZ", "MNpostv", "SN_NS22", "SN_NS5", "SN_NS17", "SN_YF5cil", "SN_NS19", "SN_NS20", "SN_NS6", "SN_NS29", "SN_NS27", "SN_NS18", "SN_NS16", "INbigloop", "SNadNS22", "Ser-tr1", "SNFVa", "INUturn", "INATO_pyg", "INcomm-Upstairs", "INasc_pyg", "SNblunt", "INsplitBronto", "INleucoPU", "MNladder", "INdecussfoot", "INdecusshook", "INpreLadder", "INcomm-Upcross", "MNacic", "INsplitVent", "INcomm-DownL", "MNarm", "INsnl", "MNantacic", "MNpostacic", "SN_DLSO1.3_PDF", "MNakro", "MNche", "MNgland_head", "SNstiff", "SNbronto", "INcommascFV", "SNpygM", "INR1", "SN_DLSO1.1_1l_NP", "SN_DLSO1.1", "INmus", "INdescLuqinPDF", "INcosplit", "INcommdescFVa", "MNheadV", "INsqNSasym", "INbiaxHmid", "MNantelope", "MNsmile")
#check if there are duplicates
duplicated(non_neuronal_celltype_names)
length(non_neuronal_celltype_names)
length(neuronal_celltype_names)

#assign column names to matrix
colnames(segmental_neuronal_celltype_complement) <- as.character(neuronal_celltype_names)
colnames(segmental_non_neuronal_celltype_complement) <- as.character(non_neuronal_celltype_names)

#assign row names to matrix
#define the six body region names
body_region_names <- c('head','sg0', 'sg1', 'sg2', 'sg3', 'pyg')

rownames(segmental_neuronal_celltype_complement) <- as.character(body_region_names)
rownames(segmental_non_neuronal_celltype_complement) <- as.character(body_region_names)

#plot heatmap with heatmaply
library(heatmaply)

hm_neuronal <- heatmaply(segmental_neuronal_celltype_complement,
          colors = c("white","#F28D24","#AEC7E8","#2078B5"),
          column_text_angle = 90,row_text_angle = 0,
          Rowv=NA,
          dist_method="euclidean",
          hclustfun = hclust,
          hclust_method='ward.D' ,
          dendrogram = c("col"),
          show_dendrogram = c(F, F),
          grid_gap = 1,
          hide_colorbar = T,
          plot_method = c("plotly"),
          fontsize_row = 30,
          fontsize_col = 10,
          revC=T
          )
hm_neuronal

hm_non_neuronal <- heatmaply(segmental_non_neuronal_celltype_complement,
          colors = c("white","grey90","#CA1523","#A373AB","grey60"),
          column_text_angle = 90,row_text_angle = 0,
          Rowv=NA,
          dist_method="euclidean",
          hclustfun = hclust,
          hclust_method = "ward.D",
          dendrogram = c("col"),
          show_dendrogram = c(F, F),
          grid_gap = 1,
          hide_colorbar = T,
          plot_method = c("plotly"),
          fontsize_row = 30,
          fontsize_col = 10,
          revC=T
)
hm_non_neuronal



####################
#Saving files
####################

#save the heatmaps as webpage in the Viewer window -> Export -> Save as Web page
