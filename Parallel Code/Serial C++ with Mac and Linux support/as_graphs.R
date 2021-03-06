# David G. Anderson
#
# Plotting instructions for
# Autonomous Systems graph sparsification
# example
#
# 2016
#

library(igraph)

# import data
g <- read.graph("as.txt", format="edgelist")
p <- read.table("p.txt",header=FALSE)

# create minimum spanning tree from selected edges
nE <- dim(p)[1]-1
E(g)$color="#55555533"
E(g)$color[p[1:nE,1]]="red"
g_min <- subgraph.edges(g, p[1:nE,1])

# a force-directed algorithm for graph plotting
coords_min <- layout.fruchterman.reingold(g_min)

# plot
E(g)$color[p[1:nE,1]]="blue"
plot(g,layout=coords_min,vertex.label=NA,vertex.size=0,edge.arrow.size=0) + 
  title("Minimum spanning tree of Autonomous Systems data")
# zoom graph for more clarity
