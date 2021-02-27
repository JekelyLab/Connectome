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
presyn_name <- c("PRC", "IN1", "INton", "INpreSer", "cPRC", "INRGWa1", "IN_NOS", "Ser-h1", "MC", "cMNPDF", "cMNATO", "cMNdc", "SNnuch", "SNnuchNS", "SNnuchBx", "SN_LHSO2golden", "SNhorn", "SNhook", "MNant", "SNlasso", "INlasso", "INdecussPre", "INpreMN", "SNMIP-vc", "SNmus", "SNantlerPDF", "SNPDF-dc", "SN_DLSO1.2_PDF", "SN_DLSO", "SN_DLSO1.2", "SN_DLSO2", "SN_LHSO2", "eyespot_PRC_R3", "eyespot_PRC_R1", "MS1", "MS4", "MS3", "MS5", "SN_ASTC", "SN_DLSO3_PDF", "SN_DLSO3", "MC2biax", "SN_IRP2-burs", "SN_IRP2-FMRF", "SN_MIP4", "SN_MIP1", "MNspinning", "SN_WLD1", "SN_DSO17", "SNbicil", "SNasym", "SN47Ach", "INCrossbow", "pygPBunp", "INarc1", "INarc2", "INsn1", "INrope", "Loop", "INCM", "MNspider-ant", "MNspider-post", "MNbiramous", "MNring", "MNcrab", "MNhose", "MNbow", "MNwave", "MNacicX", "INchaeMech", "chaeMech", "INFVa_pyg", "INsplitCR ", "INsplitCRATO", "INMC3", "INbiax", "INsplitPB_RF/Ya", "INsplitPBant", "INsplitPB", "ventralParaPB2", "interparaPM1", "cioMNcover", "INsplitPUh", "vMN3", "vMN1", "vMN2", "MC3cover", "antPUc1", "spinPU1", "VentraltrunkPUunp", "ventralParaPU1", "hPUc1", "pygCirrusPU", "hCirrusPUc1", "ventralpostParaPUc2", "hPUc3", "dsoPUc1unp", "pygCirrusPU2", "hPU2l_asymPDF", "hCR2", "VentralpygCR2", "doCRunp", "SN_DLSO1.0", "INdc", "INbackcross3", "MNmouth", "MNob-ipsi", "INbiaxLeg", "INbiaxH", "SN_NS1", "SN_NS3", "SN_NS4", "SN_NS15", "SNaant_3", "SNPDF_pyg", "INfoot", "INW", "INmidL", "INhook", "INlat1", "INlat2", "INSturn", "INcross", "INpro", "INproT2", "INpear", "INZ", "MNpostv", "SN_NS22", "SN_NS5", "SN_NS17", "SN_YF5cil", "SN_NS19", "SN_NS20", "SN_NS6", "SN_NS29", "SN_NS27", "SN_NS18", "SN_NS16", "INbigloop", "SNadNS22", "Ser-tr1", "SNFVa", "INUturn", "INATO_pyg", "INcomm-Upstairs", "INasc_pyg", "SNblunt", "INsplitBronto", "INleucoPU", "MNladder", "INdecussfoot", "INdecusshook", "INpreLadder", "INcomm-Upcross", "MNacic", "INsplitVent", "INcomm-DownL", "MNarm", "INsnl", "MNantacic", "MNpostacic", "SN_DLSO1.3_PDF", "MNakro", "MNche", "MNgland_head", "SNstiff", "SNbronto", "INcommascFV", "SNpygM", "INR1", "SN_DLSO1.1_1l_NP", "SN_DLSO1.1", "INmus", "INdescLuqinPDF", "INcosplit", "INcommdescFVa", "MNheadV", "INsqNSasym", "INbiaxHmid", "MNantelope", "MNsmile")

postsyn_name <- c("PRC", "IN1", "INton", "INpreSer", "cPRC", "INRGWa1", "IN_NOS", "Ser-h1", "MC", "cMNPDF", "cMNATO", "cMNdc", "SNnuch", "SNnuchNS", "SNnuchBx", "SN_LHSO2golden", "SNhorn", "SNhook", "MNant", "SNlasso", "INlasso", "INdecussPre", "INpreMN", "SNMIP-vc", "SNmus", "SNantlerPDF", "SNPDF-dc", "SN_DLSO1.2_PDF", "SN_DLSO", "SN_DLSO1.2", "SN_DLSO2", "SN_LHSO2", "eyespot_PRC_R3", "eyespot_PRC_R1", "MS1", "MS4", "MS3", "MS5", "SN_ASTC", "SN_DLSO3_PDF", "SN_DLSO3", "MC2biax", "SN_IRP2-burs", "SN_IRP2-FMRF", "SN_MIP4", "SN_MIP1", "MNspinning", "SN_WLD1", "SN_DSO17", "SNbicil", "SNasym", "SN47Ach", "INCrossbow", "pygPBunp", "INarc1", "INarc2", "INsn1", "INrope", "Loop", "INCM", "MNspider-ant", "MNspider-post", "MNbiramous", "MNring", "MNcrab", "MNhose", "MNbow", "MNwave", "MNacicX", "INchaeMech", "chaeMech", "INFVa_pyg", "INsplitCR ", "INsplitCRATO", "INMC3", "INbiax", "INsplitPB_RF/Ya", "INsplitPBant", "INsplitPB", "ventralParaPB2", "interparaPM1", "cioMNcover", "INsplitPUh", "vMN3", "vMN1", "vMN2", "MC3cover", "antPUc1", "spinPU1", "VentraltrunkPUunp", "ventralParaPU1", "hPUc1", "pygCirrusPU", "hCirrusPUc1", "ventralpostParaPUc2", "hPUc3", "dsoPUc1unp", "pygCirrusPU2", "hPU2l_asymPDF", "hCR2", "VentralpygCR2", "doCRunp", "SN_DLSO1.0", "INdc", "INbackcross3", "MNmouth", "MNob-ipsi", "INbiaxLeg", "INbiaxH", "SN_NS1", "SN_NS3", "SN_NS4", "SN_NS15", "SNaant_3", "SNPDF_pyg", "INfoot", "INW", "INmidL", "INhook", "INlat1", "INlat2", "INSturn", "INcross", "INpro", "INproT2", "INpear", "INZ", "MNpostv", "SN_NS22", "SN_NS5", "SN_NS17", "SN_YF5cil", "SN_NS19", "SN_NS20", "SN_NS6", "SN_NS29", "SN_NS27", "SN_NS18", "SN_NS16", "INbigloop", "SNadNS22", "Ser-tr1", "SNFVa", "INUturn", "INATO_pyg", "INcomm-Upstairs", "INasc_pyg", "SNblunt", "INsplitBronto", "INleucoPU", "MNladder", "INdecussfoot", "INdecusshook", "INpreLadder", "INcomm-Upcross", "MNacic", "INsplitVent", "INcomm-DownL", "MNarm", "INsnl", "MNantacic", "MNpostacic", "SN_DLSO1.3_PDF", "MNakro", "MNche", "MNgland_head", "SNstiff", "SNbronto", "INcommascFV", "SNpygM", "INR1", "SN_DLSO1.1_1l_NP", "SN_DLSO1.1", "INmus", "INdescLuqinPDF", "INcosplit", "INcommdescFVa", "MNheadV", "INsqNSasym", "INbiaxHmid", "MNantelope",
                  "MNsmile", "akrotroch", "crescent cell", "prototroch", "nuchal cilia", "metatroch", "paratroch", "spinGland", "covercell", "ciliatedGland", "eyespot pigment cell", "pigment cell AE", "bright droplets parapodial", "macrophage", "yolk cover cell", "flat glia", "EC rad glia ", "MVGland", "microvillarCell", "protonephridium", "nephridium", "nephridiumTip", "chaeta", "aciculoblast", "circumacicular", "hemichaetal", "ER circumchaetal", "noER circumchaetal", "EC circumchaetal", "HeadGland", "InterparaGland", "spinMicroGland", "CB pigment", "vacuolar cell_head", "Glia pigmented", "pygidial pigment", "meso", "MUSac_notA", "MUSac_notP", "MUSac_notM", "MUSac_neuAV", "MUSac_neuPD", "MUSac_neuPV", "MUSac_neuDy", "MUSac_neuDx", "MUSac_neuDach", "MUSac_neure", "MUSac_i", "MUSob-ant_re", "MUSob-ant_arc", "MUSob-ant_m-pp", "MUSob-ant_ml-pp", "MUSob-ant_l-pp", "MUSob-ant_trans", "MUSob-post_notD", "MUSob-post_neuDlong", "MUSob-post_neuDprox", "MUSob-post_neuDdist", "MUSob-post_neuV", "MUSob-post_notV", "MUSob-postM", "MUSob-post_noty", "MUSob-post_i", "MUSchae_notDob", "MUSchae_notD", "MUSchae_notDn", "MUSchae_notA", "MUSchae_notAac", "MUSchae_notAre", "MUSchae_neuVob", "MUSchae_neuDac", "MUSchae_neuAVo", "MUSchae_neuAVt", "MUSchae_Are", "MUStrans_pyg", "MUSlong_D", "MUSlong_V", "MUSax", "MUSring", "MUSph", "MUSll", "MUSant", "MUSly", "MUSpl", "MUSci", "MUSch", "MUSpx", "MUSpr", "MUStri", "MUSmed_head", "EC")

duplicated(postsyn_name)
postsyn_name[28]

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

write.csv(synapse_matrix, file = "Celltype_connectivity_matrix.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")
