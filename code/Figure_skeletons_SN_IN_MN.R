#This code was used to generate Figure supplement in the Veraszto et al. 2021 connectome paper to plot the retio of input to output synapses (I-O)/(I+O) 
#for sensory, inter and motoneurons, as a function of cable length
#the data were downloaded from the Catmaid measurements widget

library(tidyverse)
library(ggplot2)
setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Figures/Figure_skeleton_stats_SN_IN_MN/')

SN <- read.csv2('SN_skeleton_measurements.csv', sep = ",")
IN <- read.csv2('IN_skeleton_measurements.csv', sep = ",")
MN <- read.csv2('MN_skeleton_measurements.csv', sep = ",")

SN <- as_tibble(SN)
IN <- as_tibble(IN)
MN <- as_tibble(MN)

SN <- SN %>% 
   mutate(neuron_type = 'sensory neuron')
IN <- IN %>% 
  mutate(neuron_type = 'interneuron')
MN <- MN %>% 
  mutate(neuron_type = 'motoneuron')

SN_IN <- full_join(SN,IN)
SN_IN_MN <- full_join(SN_IN,MN)

#add a column with input to output ratio
SN_IN_MN <- mutate(SN_IN_MN, in_out_ratio = (N.inputs-N.outputs)/(N.inputs+N.outputs))


SN_IN_MN

{
p1 <- ggplot(data=SN_IN_MN, mapping = aes(x=Smooth.cable..nm./1000, y=in_out_ratio, color=neuron_type, shape=neuron_type, size=N.inputs)) +
  geom_jitter(stroke=0, alpha=0.6, width=0, height = 0.01)+
  scale_x_continuous(trans = "log10",  limits=c(50,3000), breaks=c(100,1000,3000))+   #change "log10" to "identity" to remove log scale
  theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
  labs(x='Cable length (µm)', y='Ratio (I-O) / (I+O)')+
  scale_size_area(max_size=6)

p2 <- ggplot(data=SN_IN_MN, mapping = aes(x=Smooth.cable..nm./1000, y=in_out_ratio, color=neuron_type, shape=neuron_type, size=N.outputs)) +
  geom_jitter(stroke=0, alpha=0.6, width=0, height = 0.01)+
  scale_x_continuous(trans = "log10",  limits=c(50,3000), breaks=c(100,1000,3000))+   #change "log10" to "identity" to remove log scale
  theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
  labs(x='Cable length (µm)', y='Ratio (I-O) / (I+O)')+
  scale_size_area(max_size=6) 

#to plot multiple plots use
library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1, labels  = "AUTO", hjust = c(0, 0), font.label = list(size = 14, face = "plain", family = 'sans'))

arrange <- ggarrange(p1, p2, ncol = 2, nrow = 1, labels  = "AUTO", hjust = c(0, 0), font.label = list(size = 14, face = "plain", family = 'sans'))
# Saving R ggplot with R ggsave Function
ggsave("SN_IN_MN.pdf", width = 30, height = 10, limitsize = FALSE, 
       units = c("cm"), arrange)
}






#synapse distribution plots

SN_in <- read.csv2('SN_Radial_density_of_input_synapses.csv', sep = ",")
SN_out <- read.csv2('SN_Radial_density_of_output_synapses.csv', sep = ",")
IN_in <- read.csv2('IN_Radial_density_of_input_synapses.csv', sep = ",")
IN_out <- read.csv2('IN_Radial_density_of_output_synapses.csv', sep = ",")
MN_in <- read.csv2('MN_Radial_density_of_input_synapses.csv', sep = ",")
MN_out <- read.csv2('MN_Radial_density_of_output_synapses.csv', sep = ",")


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

{
syn1 <- ggplot(data=SN_IN_MN_in_out_tb)+
  stat_smooth(mapping = aes(x=distance_from_soma, y=IN_in), color='grey95', alpha=0.3, span=0.1)+
  geom_line(mapping = aes(x=distance_from_soma, y=IN_in), color='red', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
  geom_point(mapping = aes(x=distance_from_soma, y=IN_in), color='red', alpha=0.5, shape=17)+
  stat_smooth(mapping = aes(x=distance_from_soma, y=IN_out), color='grey95', alpha=0.3, span=0.1)+
  geom_line(mapping = aes(x=distance_from_soma, y=IN_out), color='blue', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
  geom_point(mapping = aes(x=distance_from_soma, y=IN_out), color='blue', alpha=0.5, shape=1)+
  theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
  ylim(0,0.55)+
  labs(x='distance from soma (µm)', y='mean synapse number')+
  draw_plot_label(label = "interneurons", size = 12, x = 20, y=0.55, fontface = "plain")+
  draw_plot_label(label = "incoming", size = 12, x = 3, y=0.4, fontface = "plain", color='red', alpha=0.6)+
  draw_plot_label(label = "outgoing", size = 12, x = 15, y=0.3, fontface = "plain", color='blue', alpha=0.6)


syn2 <- ggplot(data=SN_IN_MN_in_out_tb)+
  stat_smooth(mapping = aes(x=distance_from_soma, y=MN_in), color='grey95', alpha=0.3, span=0.1)+
  geom_line(mapping = aes(x=distance_from_soma, y=MN_in), color='red', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
  geom_point(mapping = aes(x=distance_from_soma, y=MN_in), color='red', alpha=0.5, shape=17)+
  stat_smooth(mapping = aes(x=distance_from_soma, y=MN_out), color='grey95', alpha=0.3, span=0.1)+
  geom_line(mapping = aes(x=distance_from_soma, y=MN_out), color='blue', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
  geom_point(mapping = aes(x=distance_from_soma, y=MN_out), color='blue', alpha=0.5, shape=1)+
  theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
  ylim(0,0.55)+
  labs(x='distance from soma (µm)', y='mean synapse number')+
  draw_plot_label(label = "motoneurons", size = 12, x = 35, y=0.55, fontface = "plain")+
  draw_plot_label(label = "incoming", size = 12, x = 25, y=0.2, fontface = "plain", color='red', alpha=0.6)+
  draw_plot_label(label = "outgoing", size = 12, x = 76, y=0.3, fontface = "plain", color='blue', alpha=0.6)


syn3 <- ggplot(data=SN_IN_MN_in_out_tb)+
  stat_smooth(mapping = aes(x=distance_from_soma, y=SN_in), color='grey95', alpha=0.3, span=0.1)+
  geom_line(mapping = aes(x=distance_from_soma, y=SN_in), color='red', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
  geom_point(mapping = aes(x=distance_from_soma, y=SN_in), color='red', alpha=0.5, shape=17)+
  stat_smooth(mapping = aes(x=distance_from_soma, y=SN_out), color='grey95', alpha=0.3, span=0.1)+
  geom_line(mapping = aes(x=distance_from_soma, y=SN_out), color='blue', shape=2,stat="smooth",span=0.1, alpha=0.3,size=1)+
  geom_point(mapping = aes(x=distance_from_soma, y=SN_out), color='blue', alpha=0.5, shape=1)+
  theme(panel.background = element_rect(fill = "grey95", color = "grey"))+
  ylim(0,0.55)+
  labs(x='distance from soma (µm)', y='mean synapse number')+
  draw_plot_label(label = "sensory neurons", size = 12, x = 2, y=0.55, fontface = "plain")+
  draw_plot_label(label = "incoming", size = 12, x = -2, y=0.03, fontface = "plain", color='red', alpha=0.6)+
  draw_plot_label(label = "outgoing", size = 12, x = 24, y=0.3, fontface = "plain", color='blue', alpha=0.6)



library("cowplot")
arrange2 <- 
  ggdraw() +
  draw_plot(p1, x = 0, y = 0.5, width = .5, height = .5) +
  draw_plot(p2, x = .5, y = 0.5, width = .5, height = .5) +
  draw_plot(syn1, x = 0, y = -0.02, width = .3, height = 0.5) +
  draw_plot(syn2, x = 0.33, y = -0.02, width = .3, height = 0.5) +
  draw_plot(syn3, x = 0.66, y = -0.02, width = .3, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, 0.5, 0,0.33,0.66), y = c(1.01, 1.01, 0.49,0.49,0.49), fontface = "plain")+
  theme(plot.margin = unit(c(1,0,2,0), "mm")) #set margins, top, right, bottom, left


# Saving R ggplot with R ggsave Function
ggsave("SN_IN_MN_synapses.pdf", width = 30, height = 19, limitsize = FALSE, 
       units = c("cm"), arrange2)
}
}
  