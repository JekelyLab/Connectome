#R code to generate the anatomical overview images of the Platynereis connectome in Figure 1 of the Veraszto et al 2021 connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely Feb 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)

# load nat and all associated packages, incl catmaid
library(natverse)
library(magick)
library(nat)
options(nat.plotengine = 'rgl')
require("graphics")
# catmaid connection, needs username, password AND token
# load nat and all associated packages, incl catmaid

# catmaid connection, needs username, password AND token - weird!
# can run this separate file using source function
conn <- source("~/R/conn.R")

#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))


#see the available volumes
catmaid_get_volumelist(conn = NULL, pid = 11)

#read volumes
{
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)
}


#catmaid_get_connector_table("^connectome$", pid= 11, direction = "incoming", conn = conn_http1)

#read cells
{
neurons = nlapply(read.neurons.catmaid("^connectome_neuron$", pid=11, conn = conn_http1,
                                       fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
connectome = nlapply(read.neurons.catmaid("^connectome$", pid=11, conn = conn_http1,
                                          fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
Sensoryneuron = nlapply(read.neurons.catmaid("^connectome_Sensory_neuron$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
Motorneuron = nlapply(read.neurons.catmaid("^connectome_motorneuron$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
Interneuron = nlapply(read.neurons.catmaid("^connectome_interneuron$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
chaeta = nlapply(read.neurons.catmaid("^chaeta$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
muscle = nlapply(read.neurons.catmaid("^muscle$", pid=11, conn = conn_http1,
                                      fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
endoderm = nlapply(read.neurons.catmaid("^endoderm$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
epithelia = nlapply(read.neurons.catmaid("^epithelia_cell$", pid=11, conn = conn_http1,
                                         fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
gland = nlapply(read.neurons.catmaid("^gland cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
Ciliary_band_cell = nlapply(read.neurons.catmaid("^Ciliary_band_cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
glia = nlapply(read.neurons.catmaid("^glia cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
               function(x) smooth_neuron(x, sigma=6000))
pnb = nlapply(read.neurons.catmaid("^pnb$", pid=11, conn = conn_http1,
                                   fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
notopodium = nlapply(read.neurons.catmaid("^notopodium$", pid=11, conn = conn_http1,
                                   fetch.annotations = FALSE),
              function(x) smooth_neuron(x, sigma=6000))
neuropodium = nlapply(read.neurons.catmaid("^neuropodium$", pid=11, conn = conn_http1,
                                   fetch.annotations = FALSE),
              function(x) smooth_neuron(x, sigma=6000))

#these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
}

#check if there are any cells with two or more tagged somas
sum = summary(epithelia)
attributes(sum)
sum[sum$nsoma!=1,]
as.numeric(rownames(sum[sum$nsoma==2,]))


##########################################
# 3d plotting
#########################################

plot_background <- function(x){
nopen3d() # opens a pannable 3d window
par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
       col="white") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.3,
       col="black")
plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.5,
       col="grey") 
par3d(zoom=0.48)
}


#extract connectors to be able to plot them by unique colours
{
  SN_conn <- connectors(Sensoryneuron)
str(SN_conn)
presyn_SN_conn <- SN_conn[SN_conn$prepost == 0,]
postsyn_SN_conn <- SN_conn[SN_conn$prepost == 1,]

IN_conn <- connectors(Interneuron)
presyn_IN_conn <- IN_conn[IN_conn$prepost == 0,]
postsyn_IN_conn <- IN_conn[IN_conn$prepost == 1,]

#we can also subset with the subset shorthand function  
MN_conn <- connectors(Motorneuron)
presyn_MN_conn <- subset(MN_conn, prepost == 0)
postsyn_MN_conn <- subset(MN_conn, prepost == 1)
}

#plot only the synapses coloured by SN MN IN
#cb friendly colour codes interneuron = "#CC79A7", motoneuron = "#0072B2",  `sensory neuron` = "#E69F00"
{plot_background()
#plot only the presyn connectors
plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z, size=5, alpha=0.5, col="#E69F00", add=T)
plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z, size=5, alpha=0.5, col="#0072B2", add=T)
plot3d(presyn_IN_conn$x+1, presyn_IN_conn$y, presyn_IN_conn$z, size=5, alpha=0.5, col="#CC79A7", add=T)
}
rgl.snapshot("pictures/connectome_SN_IN_MN_synapses_ventral.png")

#we define a z clipping plane for the frontal view
{
par3d(windowRect = c(20, 30, 800, 800)) #resize for frontal view
nview3d("frontal", extramat=rotationMatrix(0.4, 1, 0.1, 0))
clipplanes3d(0, 0, -1, 60000)
par3d(zoom=0.63)
}
rgl.snapshot("pictures/connectome_SN_IN_MN_synapses_frontal.png")

#plot lateral view
{
  plot_background()
  #plot only the presyn connectors
  plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z, size=5, alpha=0.5, col="#E69F00", add=T)
  plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z, size=5, alpha=0.5, col="#0072B2", add=T)
  plot3d(presyn_IN_conn$x+1, presyn_IN_conn$y, presyn_IN_conn$z, size=5, alpha=0.5, col="#CC79A7", add=T)
  nview3d("left", extramat=rotationMatrix(-pi/2, pi, -0.2, 0))
  clipplanes3d(1, 0, 0.16, -75700)
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.07,
       col="#E2E2E2") 
  par3d(zoom=0.52)
}

rgl.snapshot("pictures/connectome_SN_IN_MN_synapses_left.png")

close3d()



#plot the cells with soma coloured by SN MN IN
{
plot_background()
plot3d(Sensoryneuron, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=0.6, col="#E69F00")
plot3d(Motorneuron, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=0.6, col="#0072B2",)
plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=0.6, col="#CC79A7")
}
rgl.snapshot("pictures/connectome_SN_IN_MN_cells_ventral.png")

#we define a z clipping plane for the frontal view
{
nview3d("frontal", extramat=rotationMatrix(0.4, 1, 0.1, 0))
clipplanes3d(0, 0, -1, 60000)
par3d(windowRect = c(20, 30, 800, 800))
par3d(zoom=0.63)
}
rgl.snapshot("pictures/connectome_SN_IN_MN_cells_frontal.png")
close3d()

#plot the same cells without a soma
{
  plot_background()
  plot3d(Sensoryneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, alpha=0.6, col="#E69F00")
  plot3d(Motorneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=0.6, col="#0072B2",)
  plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=0.6, col="#CC79A7")
}
rgl.snapshot("pictures/connectome_SN_IN_MN_cells_nosoma_ventral.png")

#we define a z clipping plane for the frontal view
{
  nview3d("frontal", extramat=rotationMatrix(0.4, 1, 0.1, 0))
  clipplanes3d(0, 0, -1, 60000)
  par3d(windowRect = c(20, 30, 800, 800))
  par3d(zoom=0.63)
}
rgl.snapshot("pictures/connectome_SN_IN_MN_cells_nosoma_frontal.png")
close3d()

#plot effectors only
{
plot_background()
plot3d(gland, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#0072B2")
plot3d(Ciliary_band_cell, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="#E69F00") 
plot3d(muscle, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col=hcl.colors(1200, palette='Reds'))
}
rgl.snapshot("pictures/connectome_effectors_frontal.png")

#add all other cells
{
plot3d(Sensoryneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.6, col="#E69F00")
plot3d(Motorneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=0.6, col="#0072B2",)
plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=0.6, col="#CC79A7")
plot3d(glia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(95, palette='Peach'))
plot3d(pnb, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.2,
       col=hcl.colors(2200, palette='Mint', rev=T))
plot3d(epithelia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.3,
       col=hcl.colors(1172, palette='Blues')) 
}

rgl.snapshot("pictures/connectome_body_all_cells_ventral.png")



##############################################
#Export images for Video1
##############################################
#extract connectors for left only to be able to plot them by unique colours
{
  skids_to_plot_left <- skids_by_2annotations('Sensory neuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  SN_conn_left <- connectors(skeletons_to_plot_left)
  str(SN_conn_left)
  presyn_SN_conn_left <- SN_conn_left[SN_conn_left$prepost == 0,]
  postsyn_SN_conn_left <- SN_conn_left[SN_conn_left$prepost == 1,]
  
  skids_to_plot_left <- skids_by_2annotations('interneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  
  IN_conn_left <- connectors(skeletons_to_plot_left)
  presyn_IN_conn_left <- IN_conn_left[IN_conn_left$prepost == 0,]
  postsyn_IN_conn_left <- IN_conn_left[IN_conn_left$prepost == 1,]
  
  skids_to_plot_left <- skids_by_2annotations('motorneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  #we can also subset with the subset shorthand function  
  MN_conn_left <- connectors(skeletons_to_plot_left)
  presyn_MN_conn_left <- subset(MN_conn_left, prepost == 0)
  postsyn_MN_conn_left <- subset(MN_conn_left, prepost == 1)
}

#define two windows, plot background
{
nopen3d() # opens a pannable 3d window
mfrow3d(1, 2)  #defines the two scenes
par3d(windowRect = c(20, 30, 1200, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
par3d(zoom=0.48)
next3d(clear=F)
nview3d("left", extramat=rotationMatrix(-pi/2, pi, -0.2, 0))
#clipplanes3d(1, 0, 0.16, -75700)
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
       col="white") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
       col="#E2E2E2") 
par3d(zoom=0.48)
}

rgl.snapshot("pictures/Video1B_1A.png")

#function to retrieve skids based on two annotations
skids_by_2annotations <- function(annotation1,annotation2){
  annotations_cells = list()
  annotation1 <- paste("^", annotation1, "$", sep="")
  annotations_cells[[1]] <- catmaid_get_annotations_for_skeletons(annotation1, pid = 11)
  #we retrieve those skeletons that are also annotated with right_side
  return(unlist(lapply(annotations_cells,function(x) x[x$annotation==annotation2,1])))
}


#plot ventral views and side views of left side cells only
{
next3d(clear=F)
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
       col="grey85")
skids_to_plot_left <- skids_by_2annotations('acicula','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="grey55")
rgl.snapshot("pictures/Video1B_1B.png")


next3d(clear=F)
plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
       col="grey90") 
skids_to_plot_left <- skids_by_2annotations('chaeta','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=1, col="grey90")
rgl.snapshot("pictures/Video1B_1C.png")


next3d(clear=F)
plot3d(endoderm, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="grey55")
skids_to_plot_left <- skids_by_2annotations('endoderm','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="grey55")
rgl.snapshot("pictures/Video1B_2.png")

next3d(clear=F)
plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#4477AA")
skids_to_plot_left <- skids_by_2annotations('stomodeum','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="#4477AA")
rgl.snapshot("pictures/Video1B_3.png")
}


#plot synapses
{
next3d(clear=F)
plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z, size=5, alpha=1, col="#E69F00", add=T)
next3d(clear=F)
plot3d(presyn_SN_conn_left$x, presyn_SN_conn_left$y, presyn_SN_conn_left$z, size=5, alpha=1, col="#E69F00", add=T)
rgl.snapshot("pictures/Video1B_4.png")

next3d(clear=F)
plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z, size=5, alpha=1, col="#0072B2", add=T)
next3d(clear=F)
plot3d(presyn_MN_conn_left$x, presyn_MN_conn_left$y, presyn_MN_conn_left$z, size=5, alpha=1, col="#0072B2", add=T)
rgl.snapshot("pictures/Video1B_5.png")

next3d(clear=F)
plot3d(presyn_IN_conn$x+1, presyn_IN_conn$y, presyn_IN_conn$z, size=5, alpha=1, col="#CC79A7", add=T)
next3d(clear=F)
plot3d(presyn_IN_conn_left$x+1, presyn_IN_conn_left$y, presyn_IN_conn_left$z, size=5, alpha=1, col="#CC79A7", add=T)
rgl.snapshot("pictures/Video1B_6.png")
}

#continue plotting cells
{
next3d(clear=F)
plot3d(Sensoryneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=1, col="#E69F00")
skids_to_plot_left <- skids_by_2annotations('Sensory neuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=1, col="#E69F00")
par3d(zoom=0.55)
rgl.snapshot("pictures/Video1B_7.png")

next3d(clear=F)
plot3d(Motorneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#0072B2",)
skids_to_plot_left <- skids_by_2annotations('motorneuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#0072B2",)
rgl.snapshot("pictures/Video1B_8.png")

next3d(clear=F)
plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#CC79A7")
skids_to_plot_left <- skids_by_2annotations('interneuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#CC79A7")
rgl.snapshot("pictures/Video1B_9.png")


next3d(clear=F)
plot3d(gland, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="#0072B2")
skids_to_plot_left <- skids_by_2annotations('gland cell','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#0072B2")
rgl.snapshot("pictures/Video1B_10.png")

next3d(clear=F)
plot3d(Ciliary_band_cell, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="#E69F00") 
skids_to_plot_left <- skids_by_2annotations('Ciliary_band_cell','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#E69F00") 
rgl.snapshot("pictures/Video1B_11.png")

next3d(clear=F)
plot3d(muscle, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(1200, palette='Reds'))
skids_to_plot_left <- skids_by_2annotations('muscle','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(1200, palette='Reds'))
rgl.snapshot("pictures/Video1B_12.png")

next3d(clear=F)
plot3d(glia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(95, palette='Peach'))
skids_to_plot_left <- skids_by_2annotations('glia cell','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(95, palette='Peach'))
rgl.snapshot("pictures/Video1B_13.png")

next3d(clear=F)
plot3d(pnb, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(2200, palette='Mint', rev=T))
skids_to_plot_left <- skids_by_2annotations('pnb','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(2200, palette='Mint', rev=T))
rgl.snapshot("pictures/Video1B_14.png")

next3d(clear=F)
plot3d(epithelia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(1290, palette='Blues')) 
skids_to_plot_left <- skids_by_2annotations('epithelia_cell','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(1290, palette='Blues')) 
rgl.snapshot("pictures/Video1B_15.png")

next3d(clear=F)
plot3d(notopodium, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="#009E73")
skids_to_plot_left <- skids_by_2annotations('notopodium','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="#009E73")
rgl.snapshot("pictures/Video1B_16.png")

next3d(clear=F)
plot3d(neuropodium, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="#CC79A7")
skids_to_plot_left <- skids_by_2annotations('neuropodium','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="#CC79A7")
rgl.snapshot("pictures/Video1B_17.png")
#rotation

next3d(clear=F)
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
next3d(clear=F)
nview3d("left", extramat=rotationMatrix(-pi/2, pi, -0.2, 0))
next3d(clear=F)
}


#export rotation by frame for video
for (i in 1:240){
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 2)
  next3d(clear=F)
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 2)
  next3d(clear=F)
  print (i)
  #save a snapshot
  filename <- paste("pictures/Video1B_spin", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
}

############################
#plot SN IN and MN only on left side with soma

#plot one panel background
{
  nopen3d() # opens a pannable 3d window
  mfrow3d(1, 1)  #defines one scene
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.3,
         col="black")
  plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.5,
         col="grey") 
  par3d(zoom=0.48)
}

#plot SN IN MN on left side only
#cb friendly colour codes interneuron = "#CC79A7", motoneuron = "#0072B2",  `sensory neuron` = "#E69F00"
{
  skids_to_plot_left <- skids_by_2annotations('motorneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, alpha=0.6, col="#0072B2")
  rgl.snapshot("pictures/MN_left.png")
  
skids_to_plot_left <- skids_by_2annotations('Sensory neuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       add=T, alpha=0.6, col="#E69F00")
par3d(zoom=0.55)
rgl.snapshot("pictures/SN_MN_left.png")
}

#plot one panel background
{
  nopen3d() # opens a pannable 3d window
  mfrow3d(1, 1)  #defines one scene
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.3,
         col="black")
  plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.5,
         col="grey") 
  par3d(zoom=0.48)
}

{
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=0.6, col="#E69F00")
rgl.snapshot("pictures/SN_left.png")
skids_to_plot_left <- skids_by_2annotations('interneuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=0.6, col="#CC79A7")
rgl.snapshot("pictures/SN_IN_left.png")
}

#########################################################
#Make multi-panel figure
#########################################################


library(cowplot)
library(ggplot2)
library(png)

#read png
{
img1 <- readPNG("pictures/Platynereis_SEM_inverted.png")
img2 <- readPNG("pictures/connectome_SN_IN_MN_synapses_ventral.png")
img3 <- readPNG("pictures/connectome_SN_IN_MN_synapses_frontal.png")
img4 <- readPNG("pictures/connectome_SN_IN_MN_cells_ventral.png")
img5 <- readPNG("pictures/connectome_SN_IN_MN_cells_nosoma_ventral.png")
img6 <- readPNG("pictures/connectome_SN_IN_MN_cells_nosoma_frontal.png")
img7 <- readPNG("pictures/connectome_SN_IN_MN_cells_frontal.png")
img8 <- readPNG("pictures/connectome_effectors_frontal.png")
img9 <- readPNG("pictures/connectome_body_all_cells_ventral.png")
img10 <- readPNG("pictures/SN_left.png")
img11 <- readPNG("pictures/SN_MN_left.png")
img12 <- readPNG("pictures/SN_IN_left.png")
img13 <- readPNG("pictures/connectome_SN_IN_MN_synapses_left.png")
}
#convert png to image panel
{
panel1 <- ggdraw() + draw_image(img1, scale = 1)
panel2 <- ggdraw() + draw_image(img2, scale = 1) + 
  draw_label("synapses", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("ventral", x = 0.9, y = 0.1,  hjust=1,fontfamily = "sans", fontface = "plain",
             color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1)
panel3 <- ggdraw() + draw_image(img3, scale = 1.2)  +
  draw_label("frontal", x = 0.9, y = 0.1,  hjust=1,fontfamily = "sans", fontface = "plain",
             color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1)
panel4 <- ggdraw() + draw_image(img4, scale = 1)
panel5 <- ggdraw() + draw_image(img5, scale = 1) + draw_label("SN", x = 0.9, y = 0.95, fontfamily = "sans", fontface = "plain",
  color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("IN", x = 0.9, y = 0.9, fontfamily = "sans", fontface = "plain",
               color = "#CC79A7", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("MN", x = 0.9, y = 0.85, fontfamily = "sans", fontface = "plain",
             color = "#0072B2", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("neurons (no soma)", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("ventral", x = 0.9, y = 0.1,  hjust=1,fontfamily = "sans", fontface = "plain",
             color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1)
panel6 <- ggdraw() + draw_image(img6, scale = 1)
panel7 <- ggdraw() + draw_image(img7, scale = 1)
panel8 <- ggdraw() + draw_image(img8, scale = 1) + 
  draw_label("effectors", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1)
panel9 <- ggdraw() + draw_image(img9, scale = 1) + 
  draw_label("all cells", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1)
panel10 <- ggdraw() + draw_image(img10, scale = 1) + draw_label("SN", x = 0.93, y = 0.95, fontfamily = "sans", fontface = "plain",
                                                                color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("left", x = 0.93, y = 0.9, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1)  +
  draw_label("neurons with soma", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1)
panel11 <- ggdraw() + draw_image(img11, scale = 1) + draw_label("SN", x = 0.93, y = 0.95, fontfamily = "sans", fontface = "plain",
                                                                color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("MN", x = 0.93, y = 0.9, fontfamily = "sans", fontface = "plain",
             color = "#0072B2", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("left", x = 0.93, y = 0.85, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1) 
panel12 <- ggdraw() + draw_image(img12, scale = 1) + draw_label("SN", x = 0.93, y = 0.95, fontfamily = "sans", fontface = "plain",
                                                                color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("IN", x = 0.93, y = 0.9, fontfamily = "sans", fontface = "plain",
             color = "#CC79A7", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("left", x = 0.93, y = 0.85, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1) 
panel13 <- ggdraw() + draw_image(img13, scale = 1) + draw_label("SN", x = 0.93, y = 0.95, fontfamily = "sans", fontface = "plain",
                                                                color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("IN", x = 0.93, y = 0.9, fontfamily = "sans", fontface = "plain",
             color = "#CC79A7", size = 10, angle = 0, lineheight = 0.9, alpha = 1) +
  draw_label("MN", x = 0.93, y = 0.85, fontfamily = "sans", fontface = "plain",
             color = "#0072B2", size = 10, angle = 0, lineheight = 0.9, alpha = 1)  +
  draw_label("left", x = 0.9, y = 0.1, hjust=1,fontfamily = "sans", fontface = "plain",
             color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1)
}
#cb friendly colour codes interneuron = "#CC79A7", motoneuron = "#0072B2",  `sensory neuron` = "#E69F00"




######################################
#code for SN IN MN statistics plots
######################################


#This code was used to generate Figure 1 in the Veraszto et al. 2021 connectome paper to plot the ratio of input to output synapses (I-O)/(I+O) 
#for sensory, inter and motoneurons, as a function of cable length
#the data were downloaded from the Catmaid measurements widget

library(tidyverse)
library(ggplot2)

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


#read data from data folder
SN <- read.csv2('data/SN_skeleton_measurements.csv', sep = ",")
IN <- read.csv2('data/IN_skeleton_measurements.csv', sep = ",")
MN <- read.csv2('data/MN_skeleton_measurements.csv', sep = ",")

#convert to tibble
SN <- as_tibble(SN)
IN <- as_tibble(IN)
MN <- as_tibble(MN)

#add neuron type description as a new column
SN <- SN %>% 
  mutate(neuron_type = 'sensory neuron')
IN <- IN %>% 
  mutate(neuron_type = 'interneuron')
MN <- MN %>% 
  mutate(neuron_type = 'motoneuron')

#combine tibbles
SN_IN <- full_join(SN,IN)
SN_IN_MN <- full_join(SN_IN,MN)

#add a column with input to output ratio
SN_IN_MN <- mutate(SN_IN_MN, in_out_ratio = (N.inputs-N.outputs)/(N.inputs+N.outputs))

#add a column with ratio of number of branches and length ('branchiness')
SN_IN_MN <- mutate(SN_IN_MN, branchpoints_length_ratio = (N.branch.nodes)/(Smooth.cable..nm.))

#add a column with sum of input and output synapses ('total synapses')
SN_IN_MN <- mutate(SN_IN_MN, total_synapses = (N.inputs+N.outputs))

#start plotting with ggplot
library(ggplot2)
{
  p1 <- ggplot(SN_IN_MN) +
    aes(x = Smooth.cable..nm./1000, y = in_out_ratio, colour = neuron_type,
        shape=neuron_type, size = N.inputs, alpha=neuron_type) +
    # geom_point(shape = "circle") + 
    geom_jitter(stroke=0, width=0, height = 0.01)+
    scale_color_manual(values = list(interneuron = "#CC79A7", motoneuron = "#0072B2", 
                                     `sensory neuron` = "#E69F00")) +
    scale_alpha_manual(values=c(0.5,0.8,1)) +
    scale_size_area(max_size=3)+
    labs(x='Cable length (µm)', y='Ratio (I-O) / (I+O)')+
    scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,1000,2000))+   #change "log10" to "identity" to remove log scale
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))+
    geom_text( 
      data=SN_IN_MN %>% filter(N.inputs>100), # Filter data first
      aes(label=Neuron), size=2.5, alpha=0.7, nudge_x =-0.1, check_overlap = TRUE, col='black')+                                            # Apply guides function
    guides(size = guide_legend('postsynapses'), colour="none", shape="none", alpha="none")
  
  
  p2 <- ggplot(SN_IN_MN) +
    aes(x = Smooth.cable..nm./1000, y = in_out_ratio, colour = neuron_type,
        shape=neuron_type, size = N.outputs, alpha=neuron_type) +
    # geom_point(shape = "circle") + 
    geom_jitter(stroke=0, width=0, height = 0.01)+
    scale_color_manual(values = list(interneuron = "#CC79A7", motoneuron = "#0072B2", 
                                     `sensory neuron` = "#E69F00")) +
    scale_alpha_manual(values=c(0.5,0.8,1)) +
    scale_size_area(max_size=3)+
    labs(x='Cable length (µm)', y='Ratio (I-O) / (I+O)')+
    scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,1000,2000))+   #change "log10" to "identity" to remove log scale
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))+
    geom_text( 
      data=SN_IN_MN %>% filter(N.outputs>130), # Filter data first
      aes(label=Neuron), size=2.5, alpha=0.7, check_overlap = TRUE, col='black')+                                            # Apply guides function
    guides(size = guide_legend('presynapses'), colour="none", shape="none", alpha="none")
  
  
  #plot branch numbers
  p3 <- ggplot(SN_IN_MN) +
    aes(x = Smooth.cable..nm./1000, y = N.branch.nodes, colour = neuron_type,
        shape=neuron_type, size = total_synapses, alpha=neuron_type) +
    # geom_point(shape = "circle") + 
    geom_jitter(stroke=0, width=0, height = 0.01)+
    scale_color_manual(values = list(interneuron = "#CC79A7", motoneuron = "#0072B2", 
                                     `sensory neuron` = "#E69F00")) +
    scale_alpha_manual(values=c(0.5,0.8,1)) +
    scale_size_area(max_size=3)+
    ylim(0,400)+
    labs(x='Cable length (µm)', y='number of branch nodes')+
    scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,500,1000,2000))+   #change "log10" to "identity" to remove log scale
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))+
    geom_text(data=SN_IN_MN %>% filter(total_synapses>200), # Filter data first
              aes(label=Neuron), size=2.5, alpha=0.7, nudge_x =-0.1, check_overlap = TRUE, col='black')+                                            # Apply guides function
    guides(size = guide_legend("all synapses"), shape = guide_legend("neuron type"), 
           alpha = guide_legend("neuron type"), 
           colour = guide_legend(override.aes = list(size = 4), 'neuron type'))
  
  #plot branch numbers
  p4 <- ggplot(SN_IN_MN) +
    aes(x = Smooth.cable..nm./1000, y = total_synapses, colour = neuron_type,
        shape=neuron_type, size = total_synapses, alpha=neuron_type) +
    # geom_point(shape = "circle") + 
    geom_jitter(stroke=0, width=0, height = 0.01)+
    scale_color_manual(values = list(interneuron = "#CC79A7", motoneuron = "#0072B2", 
                                     `sensory neuron` = "#E69F00")) +
    scale_alpha_manual(values=c(0.5,0.8,1)) +
    scale_size_area(max_size=3)+
    ylim(0,475)+
    labs(x='Cable length (µm)', y='number of synapses')+
    scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,500,1000,2000))+   #change "log10" to "identity" to remove log scale
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))+
    geom_text(data=SN_IN_MN %>% filter(total_synapses>200), # Filter data first
              aes(label=Neuron), size=2.5, alpha=0.7, nudge_x =-0.1, check_overlap = TRUE, col='black')+                                            # Apply guides function
    guides(size = guide_legend('all synapses'), colour="none", shape="none", alpha="none")
  
  
  # Saving plots with R ggsave Function
  ggsave("plots/SN_IN_MN_A.pdf", width = 15, height = 10, limitsize = FALSE, 
         units = c("cm"), p1)
  ggsave("plots/SN_IN_MN_B.pdf", width = 15, height = 10, limitsize = FALSE, 
         units = c("cm"), p2)
  ggsave("plots/SN_IN_MN_C.pdf", width = 15, height = 10, limitsize = FALSE, 
         units = c("cm"), p3)
  ggsave("plots/SN_IN_MN_D.pdf", width = 15, height = 10, limitsize = FALSE, 
         units = c("cm"), p4)
}




#read data for synapse distribution plots
{
  SN_in <- read.csv2('data/SN_Radial_density_of_input_synapses.csv', sep = ",")
  SN_out <- read.csv2('data/SN_Radial_density_of_output_synapses.csv', sep = ",")
  IN_in <- read.csv2('data/IN_Radial_density_of_input_synapses.csv', sep = ",")
  IN_out <- read.csv2('data/IN_Radial_density_of_output_synapses.csv', sep = ",")
  MN_in <- read.csv2('data/MN_Radial_density_of_input_synapses.csv', sep = ",")
  MN_out <- read.csv2('data/MN_Radial_density_of_output_synapses.csv', sep = ",")
  
  
  #combine data into a dat.frame and add NA values to the end of the rows to bring them to the same 175 length
  SN_IN_MN_in_out <- data.frame( MN_in = c(colMeans(MN_in[,-1]), rep(NA, 175 - length(colMeans(MN_in[,-1])) )),
                                 MN_out = c(colMeans(MN_out[,-1]), rep(NA, 175 - length(colMeans(MN_out[,-1])))),
                                 IN_in = c(colMeans(IN_in[,-1]), rep(NA, 175 - length(colMeans(IN_in[,-1])))),
                                 IN_out = c(colMeans(IN_out[,-1]), rep(NA, 175 - length(colMeans(IN_out[,-1])))),
                                 SN_in = c(colMeans(SN_in[,-1]), rep(NA, 175 - length(colMeans(SN_in[,-1])))),
                                 SN_out = c(colMeans(SN_out[,-1]), rep(NA, 175 - length(colMeans(SN_out[,-1])))))
  
  SN_IN_MN_in_out_tb <- tibble(SN_IN_MN_in_out)
  
  #add a column with um from soma
  SN_IN_MN_in_out_tb <- mutate(SN_IN_MN_in_out_tb, distance_from_soma = c(0:174))
  
  SN_IN_MN_in_out_tb
}


#generate synapse distribution plots
{
  syn1 <- ggplot(data=SN_IN_MN_in_out_tb)+
    stat_smooth(mapping = aes(x=distance_from_soma, y=IN_in), color='grey95', alpha=0.3, span=0.1)+
    geom_line(mapping = aes(x=distance_from_soma, y=IN_in), color='#E69F00', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
    geom_point(mapping = aes(x=distance_from_soma, y=IN_in), color='#E69F00', alpha=0.8, shape=17)+
    stat_smooth(mapping = aes(x=distance_from_soma, y=IN_out), color='grey95', alpha=0.3, span=0.1)+
    geom_line(mapping = aes(x=distance_from_soma, y=IN_out), color='#0072B2', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
    geom_point(mapping = aes(x=distance_from_soma, y=IN_out), color='#0072B2', alpha=0.5, shape=1)+
    theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
    ylim(0,0.55)+
    labs(x='distance from soma (µm)', y='mean synapse number')+
    draw_plot_label(label = "interneurons", size = 9, x = 20, y=0.55, fontface = "plain")+
    draw_plot_label(label = "incoming", size = 9, x = 3, y=0.4, fontface = "plain", color='#E69F00', alpha=0.9)+
    draw_plot_label(label = "outgoing", size = 9, x = 15, y=0.3, fontface = "plain", color='#0072B2', alpha=0.6)+
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))
  
  syn2 <- ggplot(data=SN_IN_MN_in_out_tb)+
    stat_smooth(mapping = aes(x=distance_from_soma, y=MN_in), color='grey95', alpha=0.3, span=0.1)+
    geom_line(mapping = aes(x=distance_from_soma, y=MN_in), color='#E69F00', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
    geom_point(mapping = aes(x=distance_from_soma, y=MN_in), color='#E69F00', alpha=0.8, shape=17)+
    stat_smooth(mapping = aes(x=distance_from_soma, y=MN_out), color='grey95', alpha=0.3, span=0.1)+
    geom_line(mapping = aes(x=distance_from_soma, y=MN_out), color='#0072B2', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
    geom_point(mapping = aes(x=distance_from_soma, y=MN_out), color='#0072B2', alpha=0.5, shape=1)+
    theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
    ylim(0,0.55)+
    labs(x='distance from soma (µm)', y='mean synapse number')+
    draw_plot_label(label = "motoneurons", size = 9, x = 35, y=0.55, fontface = "plain")+
    draw_plot_label(label = "incoming", size = 9, x = 25, y=0.2, fontface = "plain", color='#E69F00', alpha=0.9)+
    draw_plot_label(label = "outgoing", size = 9, x = 76, y=0.3, fontface = "plain", color='#0072B2', alpha=0.6)+
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))
  
  
  syn3 <- ggplot(data=SN_IN_MN_in_out_tb)+
    stat_smooth(mapping = aes(x=distance_from_soma, y=SN_in), color='grey95', alpha=0.3, span=0.1)+
    geom_line(mapping = aes(x=distance_from_soma, y=SN_in), color='#E69F00', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
    geom_point(mapping = aes(x=distance_from_soma, y=SN_in), color='#E69F00', alpha=0.8, shape=17)+
    stat_smooth(mapping = aes(x=distance_from_soma, y=SN_out), color='grey95', alpha=0.3, span=0.1)+
    geom_line(mapping = aes(x=distance_from_soma, y=SN_out), color='#0072B2', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
    geom_point(mapping = aes(x=distance_from_soma, y=SN_out), color='#0072B2', alpha=0.5, shape=1)+
    theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
    ylim(0,0.55)+
    labs(x='distance from soma (µm)', y='mean synapse number')+
    draw_plot_label(label = "sensory neurons", size = 9, x = 2, y=0.55, fontface = "plain")+
    draw_plot_label(label = "incoming", size = 9, x = -2, y=0.03, fontface = "plain", color='#E69F00', alpha=0.9)+
    draw_plot_label(label = "outgoing", size = 9, x = 24, y=0.3, fontface = "plain", color='#0072B2', alpha=0.6)+
    theme_minimal()+
    theme(axis.title.x = element_text(size=9))+
    theme(axis.title.y = element_text(size=9))+
    theme(axis.text = element_text(size=7))+
    theme(legend.text = element_text(size=7))+
    theme(legend.title = element_text(size=8))+
    theme(legend.key.size = unit(3, "mm"))
  
  ggsave("plots/SN_IN_MN_C.pdf", width = 10, height = 10, limitsize = FALSE, 
         units = c("cm"), syn1)
  ggsave("plots/SN_IN_MN_D.pdf", width = 10, height = 10, limitsize = FALSE, 
         units = c("cm"), syn2)
  ggsave("plots/SN_IN_MN_E.pdf", width = 10, height = 10, limitsize = FALSE, 
         units = c("cm"), syn3)
}
#combine bottom multi-panels



#plot in a multi panel figure
Fig1top <- plot_grid(panel1, panel5, panel2, panel3, panel13,
                  NULL,NULL,NULL,NULL,NULL,
                  panel8, panel9, panel10, panel11, panel12,
          ncol=5,
          align="h",
          # A negative rel_height shrinks space between elements
          rel_widths = c(1, 1, 1, 1, 1),
          rel_heights = c(1, 0.03, 1),
          label_size = 12,
          label_y = 1.02,
          label_x = 0,
          label_fontfamily = "sans", label_fontface = "plain",
          labels=c("A","B","C","D","E", 
                   "", "", "", "", "",
                   "F","G", "H", "I", "J"))

Fig1mid <- plot_grid(p2, p1, p3,
                     ncol=3,
                     align="h",
                     rel_widths = c(1, 1, 1),
                     label_size = 12,
                     label_y = 1.1,
                     label_x = 0,
                     label_fontfamily = "sans", label_fontface = "plain",
                     labels=c("K","L","M"))

Fig1bottom <- plot_grid(p4, syn1, syn2, syn3,
                        ncol=4,
                        align="h",
                        rel_widths = c(1, 0.7, 0.7,0.7),
                        label_size = 12,
                        label_y = 1.1,
                        label_x = 0,
                        label_fontfamily = "sans", label_fontface = "plain",
                        labels=c("N","O","P", "Q"))


Fig1combined <- plot_grid(Fig1top, NULL, Fig1mid, NULL, Fig1bottom,
                  ncol=1,
                  rel_heights = c(1, 0.03, 0.4, 0.05,0.4)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), units = "pt"))

ggsave("figures/Figure1.pdf", limitsize = FALSE, 
       units = c("cm"), Fig1combined, width = 23.9, height = 27)

ggsave("figures/Figure1.png", limitsize = FALSE, 
       units = c("cm"), Fig1combined, width = 23.9, height = 27)


