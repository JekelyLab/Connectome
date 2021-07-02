#This code was used to carry out the sholl analysis of Suppl_Figure_Sholl_left_right of the 3 day old Platynereis larva 
#described in Veraszto et al. 2021
#Gaspar Jekely 2021 Feb

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 


#select some colorblind friendly color combinations
# Taken mostly from https://tradeblotter.wordpress.com/2013/02/28/the-paul-tol-21-color-salute/
#Color blind friendly palettes
########## DEFINE Color palettes #########
#From Paul Tol: https://personal.sron.nl/~pault/
Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933',
               '#CC6677', '#882255', '#AA4499', '#DDDDDD')
Tol_light <- c('#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', 
               '#99DDFF', '#44BB99', '#DDDDDD')
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
########################## SHOW as multiple PIE Charts ###################
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
pie(rep(1,length(Tol_bright)), col=Tol_bright, Tol_bright, main='Tol bright')
pie(rep(1,length(Tol_muted)), col=Tol_muted, Tol_muted, main='Tol muted')
pie(rep(1,length(Okabe_Ito)), col=Okabe_Ito, Okabe_Ito, main='Tol Okabe Ito')
pie(rep(1,length(Tol_light)), col=Tol_light, Tol_light, main='Tol light')




# load nat and all associated packages, incl catmaid
library(natverse)

#set working directory
setwd('.')

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info

#demonstration of sholl analysis and averaging for a celltype on an example
neurons = nlapply(read.neurons.catmaid("celltype1$", pid=11, fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

sholl <- sholl_analysis(
  neurons,
  start = soma(neurons),
  starting.radius = 1000,
  ending.radius = 180000,
  radius.step = 1000
)

sholl_intersections <- lapply(sholl, function(x) (x$intersections))
length(sholl_intersections)
sholl_means <- rowMeans(as.data.frame(sholl_intersections[1:length(sholl_intersections)]))

annotation_celltypelist = list()

#first we read all celltypes from 1-182 and all annotations
for (i in c(1:182)){
  annotation = paste("annotation:^celltype", i, "$", sep="")
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
  }

#define empty cells_per_celltype list with the right dimensions
#cells_per_celltype <- vector("list", length(annotation_celltypelist))


#we iterate through the celltype lists and body regions and retrieve connectivity between pre- and postsyn sets of skids
list_position=0
cycle <- 0
sholl_left_mat <- matrix(data=NA, nrow=182, ncol=180)
sholl_right_mat <- matrix(data=NA, nrow=182, ncol=180)

for (df1 in annotation_celltypelist){    #iterate through the celltype list  (all presyn cells per body region)
    left_skids <- df1[df1$annotation == 'left_side',1]
    right_skids <- df1[df1$annotation == 'right_side',1]
    cycle <- cycle + 1
    
    #skip asymmetric or midline cells (e.g. MC, pygPBunp)
    if (length(left_skids)==0 | length(right_skids)==0){sholl_left_mat[cycle,] <- 0
    sholl_right_mat[cycle,] <- 0; next}
    
    left_neurons = nlapply(read.neurons.catmaid(left_skids, pid=11, fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    right_neurons = nlapply(read.neurons.catmaid(right_skids, pid=11, fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    
    print(left_neurons)
    
    #do sholl analysis on left and right cells separately
    sholl_left <- sholl_analysis(
      left_neurons,
      start = soma(left_neurons),
      starting.radius = 1000,
      ending.radius = 180000,
      radius.step = 1000
    )

    print(right_neurons)
    
    sholl_right <- sholl_analysis(
      right_neurons,
      start = soma(right_neurons),
      starting.radius = 1000,
      ending.radius = 180000,
      radius.step = 1000
    )
    
    sholl_intersections_left <- lapply(sholl_left, function(x) (x$intersections))
    sholl_means_left <- rowMeans(as.data.frame(sholl_intersections_left[1:length(sholl_intersections_left)]))
    
    sholl_intersections_right <- lapply(sholl_right, function(x) (x$intersections))
    sholl_means_right <- rowMeans(as.data.frame(sholl_intersections_right[1:length(sholl_intersections_right)]))
    
    print (cycle)
    sholl_left_mat[cycle,] <- sholl_means_left
    sholl_right_mat[cycle,] <- sholl_means_right
}   


#plot all sholl plots
matplot(t(sholl_right_mat), type = "l")

#define smoothing function
lowess_fun <- function(x,f=0){
  return(lowess(x, f=0.04, iter=100))
}

lowess_sholl_left <- apply(sholl_left_mat, 1, lowess_fun,f=0.04)
lowess_sholl_right <- apply(sholl_right_mat, 1, lowess_fun,f=0.04)

#test plot
plot(lowess_sholl_left[[2]])

#create smoothed matrix
sholl_left_smoothed_mat <- sholl_left_mat
sholl_right_smoothed_mat <- sholl_right_mat

#add smoothed values
for (i in c(1:length(lowess_sholl_left))){
  sholl_left_smoothed_mat[i,] <- lowess_sholl_left[[i]]$y
  sholl_right_smoothed_mat[i,] <- lowess_sholl_right[[i]]$y
}

#plot the smoothed sholl diagrams
matplot(t(sholl_left_smoothed_mat), type = "l")
sholl_left_smoothed_mat
#we make a celltype name list 
celltype_names <- c("PRC", "IN1", "INton", "INpreSer", "cPRC", "INRGWa1", "IN_NOS", "Ser-h1", "MC", "cMNPDF", "cMNATO", "cMNdc", "SNnuch", "SNnuchNS", "SNnuchBx", "SN_LHSO2golden", "SNhorn", "SNhook", "MNant", "SNlasso", "INlasso", "INdecussPre", "INpreMN", "SNMIP-vc", "SNmus", "SNantlerPDF", "SNPDF-dc", "SN_DLSO1.2_PDF", "SN_DLSO", "SN_DLSO1.2", "SN_DLSO2", "SN_LHSO2", "eyespot_PRC_R3", "eyespot_PRC_R1", "MS1", "MS4", "MS3", "MS5", "SN_ASTC", "SN_DLSO3_PDF", "SN_DLSO3", "MC2biax", "SN_IRP2-burs", "SN_IRP2-FMRF", "SN_MIP4", "SN_MIP1", "MNspinning", "SN_WLD1", "SN_DSO17", "SNbicil", "SNasym", "SN47Ach", "INCrossbow", "pygPBunp", "INarc1", "INarc2", "INsn1", "INrope", "Loop", "INCM", "MNspider-ant", "MNspider-post", "MNbiramous", "MNring", "MNcrab", "MNhose", "MNbow", "MNwave", "MNacicX", "INchaeMech", "chaeMech", "INFVa_pyg", "INsplitCR ", "INsplitCRATO", "INMC3", "INbiax", "INsplitPB_RF/Ya", "INsplitPBant", "INsplitPB", "ventralParaPB2", "interparaPM1", "cioMNcover", "INsplitPUh", "vMN3", "vMN1", "vMN2", "MC3cover", "antPUc1", "spinPU1", "VentraltrunkPUunp", "ventralParaPU1", "hPUc1", "pygCirrusPU", "hCirrusPUc1", "ventralpostParaPUc2", "hPUc3", "dsoPUc1unp", "pygCirrusPU2", "hPU2l_asymPDF", "hCR2", "VentralpygCR2", "doCRunp", "SN_DLSO1.0", "INdc", "INbackcross3", "MNmouth", "MNob-ipsi", "INbiaxLeg", "INbiaxH", "SN_NS1", "SN_NS3", "SN_NS4", "SN_NS15", "SNaant_3", "SNPDF_pyg", "INfoot", "INW", "INmidL", "INhook", "INlat1", "INlat2", "INSturn", "INcross", "INpro", "INproT2", "INpear", "INZ", "MNpostv", "SN_NS22", "SN_NS5", "SN_NS17", "SN_YF5cil", "SN_NS19", "SN_NS20", "SN_NS6", "SN_NS29", "SN_NS27", "SN_NS18", "SN_NS16", "INbigloop", "SNadNS22", "Ser-tr1", "SNFVa", "INUturn", "INATO_pyg", "INcomm-Upstairs", "INasc_pyg", "SNblunt", "INsplitBronto", "INleucoPU", "MNladder", "INdecussfoot", "INdecusshook", "INpreLadder", "INcomm-Upcross", "MNacic", "INsplitVent", "INcomm-DownL", "MNarm", "INsnl", "MNantacic", "MNpostacic", "SN_DLSO1.3_PDF", "MNakro", "MNche", "MNgland_head", "SNstiff", "SNbronto", "INcommascFV", "SNpygM", "INR1", "SN_DLSO1.1_1l_NP", "SN_DLSO1.1", "INmus", "INdescLuqinPDF", "INcosplit", "INcommdescFVa", "MNheadV", "INsqNSasym", "INbiaxHmid", "MNantelope", "MNsmile")
#check if there are duplicates
duplicated(celltype_names)

#assign rownames
rownames(sholl_left_smoothed_mat) <- celltype_names
rownames(sholl_right_smoothed_mat) <- celltype_names

#remove all-zero rows
sholl_left_smoothed_mat.no0 = sholl_left_smoothed_mat[ rowSums(sholl_left_smoothed_mat)!=0, ] 
sholl_right_smoothed_mat.no0 = sholl_right_smoothed_mat[ rowSums(sholl_right_smoothed_mat)!=0, ] 

#plot heatmap with heatmaply
library(heatmaply)

#plot as heatmap
heatmaply(sholl_left_smoothed_mat.no0,
              column_text_angle = 90,row_text_angle = 0,
              fontsize_row = 3,
              fontsize_col = 3,
              Rowv=F, Colv=F,
              show_dendrogram = c(F, F),
              grid_gap = 0,
              hide_colorbar = F,
              plot_method = c("ggplot"),
              revC=F,
              col=c('white','#0072B2','#E69F00','#E69F00','#E69F00','#D55E00','#D55E00','#D55E00'),
              xlab = "distance from soma (µm)",
              ylab = "Cell types",
              main = "left side",
          key.title="Sholl value"
          )

heatmaply(sholl_right_smoothed_mat.no0,
          column_text_angle = 90,row_text_angle = 0,
          fontsize_row = 3,
          fontsize_col = 3,
          Rowv=F, Colv=F,
          show_dendrogram = c(F, F),
          grid_gap = 0,
          hide_colorbar = F,
          plot_method = c("ggplot"),
          revC=T,
          col=c('white','#316DB4','orange','orange','orange','red','red','red'),
          xlab = "distance from soma (µm)",
          ylab = "Cell types",
          main = "right side",
          key.title="Sholl value",
)

