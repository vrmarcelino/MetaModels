#### Differential Network Analyses, showing the differences in H2S consumption/production in health and CD

# Identify the H2S producers that INCREASED H2S PRODUCTION in CD,
# and the consumers that REDUCED H2S CONSUMPTION in CD

# Note

#### created:18 - Apr - 2022
## Updated 5 Oct 2022

rm(list=ls(all=TRUE))

### Load packages
library (igraph)
library(plyr) # for mapvalues
library(wesanderson)
library(dplyr) # for case_when function
library(MetBrewer)
library(ggplot2)


### Read files
edges_healthy <- read.csv("0_input/He_healthy_edges.csv", header = T, as.is = T)
nodes_healthy <- read.csv("0_input/He_healthy_nodes.csv", header = T, as.is = T)

edges_CD <- read.csv("0_input/He_CD_edges.csv", header = T, as.is = T)
nodes_CD <- read.csv("0_input/He_CD_nodes.csv", header = T, as.is = T)

nodes_healthy["diagnosis"] <- "healthy"
nodes_CD["diagnosis"] <- "CD"

nodes_healthy["flux_weighted_sum_healthy"] <- nodes_healthy$flux_weighted_sum
nodes_healthy["flux_weighted_sum_CD"] <- 0
nodes_CD["flux_weighted_sum_CD"] <- nodes_CD$flux_weighted_sum
nodes_CD["flux_weighted_sum_healthy"] <- 0
nodes_healthy["occurrences_healthy"] <- nodes_healthy$occurrences
nodes_healthy["occurrences_CD"] <- 0
nodes_CD["occurrences_CD"] <- nodes_CD$occurrences
nodes_CD["occurrences_healthy"] <- 0


merged_nodes <- bind_rows(nodes_healthy, nodes_CD)

merged_nodes_sum <- merged_nodes %>% group_by(node) %>% summarize(
  lineage=first(lineage),
  species=first(species),
  prod_cons=first(prod_cons),
  health_disease = paste(paste(diagnosis, collapse = "_and_")),
  occurrences_total=sum(occurrences),
  occurrences_healthy=sum(occurrences_healthy),
  occurrences_CD=sum(occurrences_CD),
  flux_weighted_sum_healthy=sum(flux_weighted_sum_healthy),
  flux_weighted_sum_CD=sum(flux_weighted_sum_CD),
  flux_weighted_diff=sum(flux_weighted_sum_healthy) - sum(flux_weighted_sum_CD))

# add direction of change in difference between health and CD:
# negative means that the production or consumption was higher in CD than in healthy. 
merged_nodes_sum <- merged_nodes_sum %>% mutate(direction_diff = case_when(
  flux_weighted_diff > 0 ~ "reduced_in_CD",
  flux_weighted_diff < 0 ~ "increased_in_CD",
  flux_weighted_diff == 0 ~ "no_difference"))

#transform differences into absolute values
merged_nodes_sum$flux_weighted_diff <- abs(merged_nodes_sum$flux_weighted_diff)

#write.csv(merged_nodes_sum, "2.2_Sup_TableX.csv")


# species of interest (producers that increased production, and consumers that decreased consumption in CD)
wanted_producers_nodes <- merged_nodes_sum[with(merged_nodes_sum, prod_cons == "p" & direction_diff == "increased_in_CD"), ]
wanted_consumer_nodes <- merged_nodes_sum[with(merged_nodes_sum, prod_cons == "c" & direction_diff == "reduced_in_CD"), ]

wanted_nodes <- rbind(merged_nodes_sum[1,],wanted_producers_nodes, wanted_consumer_nodes)


#### Combine edges
edges_healthy["flux_weighted_sum_healthy"] <- edges_healthy$flux_weighted_sum
edges_healthy["flux_weighted_sum_CD"] <- 0
edges_CD["flux_weighted_sum_CD"] <- edges_CD$flux_weighted_sum
edges_CD["flux_weighted_sum_healthy"] <- 0


merged_edges <- bind_rows(edges_healthy, edges_CD)

merged_edges_sum <- merged_edges %>% group_by(source,target) %>% summarize(
  occurrences_sum=sum(occurrences),
  flux_weighted_diff=sum(flux_weighted_sum_healthy) - sum(flux_weighted_sum_CD))

# add direction to edges:
merged_edges_sum <- merged_edges_sum %>% mutate(direction = case_when(
  endsWith(source, "_p") ~ "export",
  endsWith(target, "_c") ~ "import"))

#write.csv(merged_edges_sum, "check_edges.csv")

# species of interest 
wanted_producers_edges <- merged_edges_sum[with(merged_edges_sum, direction == "export" & flux_weighted_diff < 0), ]
wanted_consumer_edges <- merged_edges_sum[with(merged_edges_sum, direction == "import" & flux_weighted_diff > 0), ]

wanted_edges <- rbind(wanted_producers_edges, wanted_consumer_edges)

#transform differences into absolute values
wanted_edges$flux_weighted_diff <- abs(wanted_edges$flux_weighted_diff)


#write.csv(wanted_edges, "wanted_edges.csv")
#write.csv(wanted_nodes, "wanted_nodes.csv")

# Turning files into an igraph object:
net <- graph_from_data_frame(d = wanted_edges, vertices = wanted_nodes, directed = TRUE)

###############################
#### set plot attributes ######
###############################


##### Edges:

## set width of the arrows based on the weighted sum of fluxes 

NewMax = 15
NewMin = 0.2
NewRange = (NewMax - NewMin)

old_min_prod = min(E(net)$flux_weighted_diff[E(net)$direction == "export"])
old_min_cons = min(E(net)$flux_weighted_diff[E(net)$direction == "import"])
OldRange_producers = max(E(net)$flux_weighted_diff[E(net)$direction == "export"]) - old_min_prod
OldRange_consumers = max(E(net)$flux_weighted_diff[E(net)$direction == "import"]) - old_min_cons

new_range_producers = (((E(net)$flux_weighted_diff[E(net)$direction == "export"] - old_min_prod) * NewRange) / OldRange_producers) + NewMin
new_range_consumers = (((E(net)$flux_weighted_diff[E(net)$direction == "import"] - old_min_cons) * NewRange) / OldRange_consumers) + NewMin

E(net)$width_factor <- c(new_range_producers,new_range_consumers)
  
# color edges
#pal = wes_palette("Darjeeling1", 5, type = "discrete") # check colors
colrs_opaque <- met.brewer("Troy", 2)
colrs <- grDevices::adjustcolor(colrs_opaque, alpha = .75) # make it transparent
colrs_expanded = mapvalues(E(net)$direction, from=c("export", "import"), to=colrs)
E(net)$color <- colrs_expanded



##### Nodes:

# set different colors for  producers, metabolites and consumers
colrs_opaque <- c("#000000",met.brewer("Troy", 2))
#colrs <- grDevices::adjustcolor(colrs_opaque, alpha = 0.95) # make it transparent

colrs_expanded = mapvalues(V(net)$prod_cons, from=c("metab","p", "c"), to=colrs_opaque)
V(net)$color <- colrs_expanded

## set size of the balls based on weighted sum of fluxes

# convert to a more readable scale:
NewMax = 15
NewMin = 0.1

old_min_prod = min(V(net)$flux_weighted_diff [V(net)$prod_cons == "p"])
old_min_cons = min(V(net)$flux_weighted_diff [V(net)$prod_cons == "c"])

OldRange_producers = max(V(net)$flux_weighted_diff[V(net)$prod_cons == "p"]) - old_min_prod
OldRange_consumers = max(V(net)$flux_weighted_diff[V(net)$prod_cons == "c"]) - old_min_cons

new_range_producers = (((V(net)$flux_weighted_diff[V(net)$prod_cons == "p"] - old_min_prod) * NewRange) / OldRange_producers) + NewMin
new_range_consumers = (((V(net)$flux_weighted_diff[V(net)$prod_cons == "c"] - old_min_cons) * NewRange) / OldRange_consumers) + NewMin

V(net)$size <- c(30, new_range_producers,new_range_consumers) # added 30 for H2S (difference is zero)



top_x = 15 # add +1 to account for H2S

nodes_ordered <- V(net)$size[order(V(net)$size,decreasing=TRUE)]
threshold <- nodes_ordered[top_x]



##############################
####### PLOT #######
##############################


# sequence of numbers:
metab_coord <- cbind(2,110)

number_of_prod = length(grep("p",V(net)$prod_cons))
seq_x_p <- seq(from=number_of_prod, to=1)
prod_coord <- cbind(1, seq_x_p)

number_of_cons = length(grep("c",V(net)$prod_cons))
seq_x_c <- seq(from=200, to=21, length.out = number_of_cons)
cons_coord <- cbind(3,seq_x_c)


# But I want a semicircle instead of a flat edge:
## Thanks Jamie!!
r <- sqrt((LO[2, 1] - LO[1, 1])^2 + (LO[2, 2] - LO[1, 2])^2)
new.x <- sqrt(r^2 - (LO[, 2] - LO[1, 2])^2) + 2
new.x[1] <- LO[1, 1]
new.x[LO[, 1] == 1] <- -new.x[LO[, 1] == 1] # convert producers to a negative value, so they are in the other side
LO2 <- LO
LO2[, 1] <- new.x

pdf(file = "2.3_DiNA_tripartite_graph.pdf", width = 5, height = 5)
plot(net, edge.width=E(net)$width_factor, edge.arrow.size=0, vertex.label.cex=0.4,
     vertex.frame.color = NA, vertex.label.color="black",
     display.isolates=F,vertex.label.dist=0,
     vertex.label=ifelse(V(net)$size >= threshold, V(net)$species, NA),
     layout=LO2)

dev.off()



