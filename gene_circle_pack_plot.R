
library(ggraph)
library(igraph)
library(tidyverse)


# read the mean gene abundance data
dt <- read.csv("gene abundance_dsc.csv")
dt[1:5,]   # genes already sorted by the abundance
#    drug_class gene_name abundance
# 1 Glycopeptide     vanRO 20.570277
# 2    Rifamycin      rpoB  6.155165
# 3 Glycopeptide     vanWI  5.148098
# 4    Rifamycin     rpoB2  2.628936
# 5 Glycopeptide     vanSO  1.429324


drugs <- unique(dt$drug_class)


edge1<- rbind(data.frame(from= rep("ARG", length(drugs)), to=drugs),
  data.frame(from = dt$drug_class, to = dt$gene_name))


vertice1 <- data.frame(name=dt$gene_name, size=dt$abundance, 
                       top10=c(dt$gene_name[1:10], 
                                        rep(NA, nrow(dt)-10)),
                       top20=c(dt$gene_name[1:20], 
                               rep(NA, nrow(dt)-20)),   
                       top30=c(dt$gene_name[1:30], 
                               rep(NA, nrow(dt)-30)))

vertice2 <- data.frame(name= c("ARG",drugs), size=30,
                       top10=c(NA,drugs),
                       top20=c(NA,drugs),top30=c(NA,drugs))

vertice_all <- rbind( vertice2, vertice1)

# Rebuild the graph object
mygraph2 <- graph_from_data_frame( edge1, vertices=vertice_all)

# 
ggraph(mygraph2, layout = 'circlepack', weight=size ) + 
  geom_node_circle(aes(fill = depth)) +
  geom_node_text( aes(label=top10)) +
  #geom_node_text( aes(label=top10, filter=leaf)) +
  theme_void() + 
  theme(legend.position="FALSE") + 
  scale_fill_viridis()



