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