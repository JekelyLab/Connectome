#This code was used to generate Figure supplement in the Veraszto et al. 2021 connectome paper to plot the retio of input to output synapses (I-O)/(I+O) 
#for sensory, inter and motoneurons, as a function of cable length
#the data were downloaded from the Catmaid measurements widget

library(tidyverse)
library(ggplot2)
setwd('.')

#select some colorblind friendly color combinations
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(8, 'Dark2')
brewer.pal(8, 'Dark2')
Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
pie(rep(1,8), col=Tol_bright, Tol_bright)

SN <- read.csv2('data/SN_skeleton_measurements.csv', sep = ",")
IN <- read.csv2('data/IN_skeleton_measurements.csv', sep = ",")
MN <- read.csv2('data/MN_skeleton_measurements.csv', sep = ",")

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

#add a column with ratio of number of branches and length ('branchiness')
SN_IN_MN <- mutate(SN_IN_MN, branchpoints_length_ratio = (N.branch.nodes)/(Smooth.cable..nm.))

#add a column with sum of input and output synapses ('total synapses')
SN_IN_MN <- mutate(SN_IN_MN, total_synapses = (N.inputs+N.outputs))


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
    scale_size_area(max_size=6)+
    labs(x='Cable length (µm)', y='Ratio (I-O) / (I+O)')+
    scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,1000,2000))+   #change "log10" to "identity" to remove log scale
    theme_minimal()+
    geom_text( 
      data=SN_IN_MN %>% filter(N.branch.nodes>210 | Smooth.cable..nm./1000>800), # Filter data first
      aes(label=Neuron), size=3, alpha=0.7, check_overlap = TRUE, col='black'
    )

p2 <- ggplot(SN_IN_MN) +
  aes(x = Smooth.cable..nm./1000, y = in_out_ratio, colour = neuron_type,
      shape=neuron_type, size = N.outputs, alpha=neuron_type) +
  # geom_point(shape = "circle") + 
  geom_jitter(stroke=0, width=0, height = 0.01)+
  scale_color_manual(values = list(interneuron = "#CC79A7", motoneuron = "#0072B2", 
                                   `sensory neuron` = "#E69F00")) +
  scale_alpha_manual(values=c(0.5,0.8,1)) +
  scale_size_area(max_size=6)+
  labs(x='Cable length (µm)', y='Ratio (I-O) / (I+O)')+
  scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,1000,2000))+   #change "log10" to "identity" to remove log scale
  theme_minimal()+
  geom_text( 
    data=SN_IN_MN %>% filter(N.branch.nodes>210 | Smooth.cable..nm./1000>800), # Filter data first
    aes(label=Neuron), size=3, alpha=0.7, check_overlap = TRUE, col='black'
  )

#to plot multiple plots use
library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1, labels  = "AUTO", hjust = c(0, 0), font.label = list(size = 14, face = "plain", family = 'sans'))

arrange <- ggarrange(p1, p2, ncol = 2, nrow = 1, labels  = "AUTO", hjust = c(0, 0), font.label = list(size = 14, face = "plain", family = 'sans'))
# Saving R ggplot with R ggsave Function
ggsave("plots/SN_IN_MN_A.pdf", width = 15, height = 10, limitsize = FALSE, 
       units = c("cm"), p1)
ggsave("plots/SN_IN_MN_B.pdf", width = 15, height = 10, limitsize = FALSE, 
       units = c("cm"), p2)
ggsave("plots/SN_IN_MN.pdf", width = 30, height = 10, limitsize = FALSE, 
       units = c("cm"), arrange)
}





#synapse distribution plots
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
  draw_plot_label(label = "interneurons", size = 12, x = 20, y=0.55, fontface = "plain")+
  draw_plot_label(label = "incoming", size = 12, x = 3, y=0.4, fontface = "plain", color='#E69F00', alpha=0.9)+
  draw_plot_label(label = "outgoing", size = 12, x = 15, y=0.3, fontface = "plain", color='#0072B2', alpha=0.6)+
    theme_minimal()


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
  draw_plot_label(label = "motoneurons", size = 12, x = 35, y=0.55, fontface = "plain")+
  draw_plot_label(label = "incoming", size = 12, x = 25, y=0.2, fontface = "plain", color='#E69F00', alpha=0.9)+
  draw_plot_label(label = "outgoing", size = 12, x = 76, y=0.3, fontface = "plain", color='#0072B2', alpha=0.6)+
  theme_minimal()


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
  draw_plot_label(label = "sensory neurons", size = 12, x = 2, y=0.55, fontface = "plain")+
  draw_plot_label(label = "incoming", size = 12, x = -2, y=0.03, fontface = "plain", color='#E69F00', alpha=0.9)+
  draw_plot_label(label = "outgoing", size = 12, x = 24, y=0.3, fontface = "plain", color='#0072B2', alpha=0.6)+
  theme_minimal()

ggsave("plots/SN_IN_MN_C.pdf", width = 10, height = 10, limitsize = FALSE, 
       units = c("cm"), syn1)
ggsave("plots/SN_IN_MN_D.pdf", width = 10, height = 10, limitsize = FALSE, 
       units = c("cm"), syn2)
ggsave("plots/SN_IN_MN_E.pdf", width = 10, height = 10, limitsize = FALSE, 
                                    units = c("cm"), syn3)

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
ggsave("figures/SN_IN_MN_synapses.pdf", width = 30, height = 19, limitsize = FALSE, 
       units = c("cm"), arrange2)
}



#########################################################
#plot branch numbers

SN_IN_MN$Neuron[1:10]

ggplot(SN_IN_MN) +
  aes(x = Smooth.cable..nm./1000, y = N.branch.nodes, colour = neuron_type,
      shape=neuron_type, size = total_synapses, alpha=neuron_type) +
  # geom_point(shape = "circle") + 
  geom_jitter(stroke=0, width=0, height = 0.01)+
  scale_color_manual(values = list(interneuron = "#CC79A7", motoneuron = "#0072B2", 
                                   `sensory neuron` = "#E69F00")) +
  scale_alpha_manual(values=c(0.5,0.8,1)) +
  scale_size_area(max_size=6)+
  ylim(0,400)+
  labs(x='Cable length (µm)', y='number of branch nodes')+
  scale_x_continuous(trans = "log10",  limits=c(50,2178), breaks=c(100,500,1000,2000))+   #change "log10" to "identity" to remove log scale
  theme_minimal()+
  geom_text(data=SN_IN_MN %>% filter(N.branch.nodes>210 | Smooth.cable..nm./1000>800), # Filter data first
    aes(label=Neuron), size=3, alpha=0.7, check_overlap = TRUE, col='black')
