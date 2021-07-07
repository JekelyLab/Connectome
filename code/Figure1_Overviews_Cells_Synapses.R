#R code to generate the anatomical overview images of the Platynereis connectome in Figure 1 of the Veraszto et al 2021 connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely Feb 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=80000000000)

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
workdir <- "/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Figures/Figure1_overview/"
setwd(workdir)

#see the available volumes
catmaid_get_volumelist(conn = NULL, pid = 11)

#read volumes
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1.1))

#catmaid_get_connector_table("^connectome$", pid= 11, direction = "incoming", conn = conn_http1)

#read cells
neurons = nlapply(read.neurons.catmaid("^connectome_neuron$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000), conn = conn_http1)

connectome = nlapply(read.neurons.catmaid("^connectome$", pid=11, conn = conn_http1),
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
muscle = nlapply(read.neurons.catmaid("^muscle$", pid=11),
                function(x) smooth_neuron(x, sigma=6000), conn = conn_http1)
endoderm = nlapply(read.neurons.catmaid("^endoderm$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
epithelia = nlapply(read.neurons.catmaid("^epithelia_cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000), conn = conn_http1)
gland = nlapply(read.neurons.catmaid("^gland cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
Ciliary_band_cell = nlapply(read.neurons.catmaid("^Ciliary_band_cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
glia = nlapply(read.neurons.catmaid("^glia cell$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))

pnb = nlapply(read.neurons.catmaid("^pnb$", pid=11, 
                                    fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

#these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11, 
                                                                   fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

#check if there are any cells with two or more tagged somas
sum = summary(gland)
sum[sum$nsoma!=1,]
as.numeric(rownames(sum[sum$nsoma==2,]))

##########################################
# 3d plotting
#########################################

nopen3d() # opens apannable 3d window
mfrow3d(1, 1)  #defines the two scenes
par3d(windowRect = c(20, 30, 800, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
par3d(zoom=0.72)



plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="white"
) 


plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2"
) 

rgl.snapshot("connectome_body_1.png")

plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="black"
) 
rgl.snapshot("connectome_body_2.png")

plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey"
) 
rgl.snapshot("connectome_body_3.png")

plot3d(gland, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey"
) 
rgl.snapshot("connectome_body_3b.png")

plot3d(Ciliary_band_cell, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(333, palette='RedOr', rev = T)
) 
rgl.snapshot("connectome_body_4.png")

plot3d(muscle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
     rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
    col=hcl.colors(1200, palette='Reds')
   )

#plot3d(hckcs_mus, k=500, db=muscle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
#       col=hcl.colors(900, palette='Reds'))

rgl.snapshot("connectome_body_5.png")

plot3d(hckcs, k=500, db=neurons, WithConnectors = T, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0,
) 
rgl.snapshot("connectome_all_neurons_6.png")

plot3d(glia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(95, palette='Peach'))
rgl.snapshot("connectome_body_7.png")

plot3d(pnb, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.2,
       col=hcl.colors(2200, palette='Mint', rev=T)
) 
rgl.snapshot("connectome_body_8.png")

plot3d(epithelia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.3,
       col=hcl.colors(1172, palette='Blues')
) 
rgl.snapshot("connectome_body_9.png")







#plot3d(desmosome_connectome_new_non_muscle, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
#       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.4,
#       col=hcl.colors(2000, palette='Blues')) 


for (i in c(89:94)){
  nview3d("ventral", extramat=rotationMatrix((pi/180)*2*i, 0, 0, 1))
  Sys.sleep(15)
  plot3d(epithelia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.3,
         col=hcl.colors(1172, palette='Blues')
  ) 
  filename <- paste("test_movie", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
}
(3.1414/180)*188

#zoom in 
nview3d("ventral", extramat=rotationMatrix(pi+0.15, 0, 0, 1))
for (i in c(0:60)){
  par3d(zoom=0.72-i/200)
  filename <- paste("test_movie_zoom", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
}


#to save frames for a movie
movie3d(spin3d(axis = c(0, 0, 10), rpm = 1), duration = 1, dir=workdir, convert = NULL, clean = NULL)
