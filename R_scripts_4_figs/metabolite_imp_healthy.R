#### Produce a graph of average metabolite exchange score in healthy individuals

#### created: 15 - Sep - 2021

rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
library(dplyr)


####################################################
### Input and parse tables

input_exchanges_fp <- "0_input/producers_consumers.csv"
metadata_fp <- "0_input/metadata_rewiring_microbiome.csv"
metab_names_fp <- "0_input/bigg_models_w_classes_curated_simplified_names.tsv"

met_exchanges <- read.csv(input_exchanges_fp)
met_exchanges$sample_id <- gsub('_cat', '', met_exchanges$sample_id) 

metad <- read.csv(metadata_fp)
rownames(metad) <- metad[,1]

metab_names <- read.csv(metab_names_fp, sep = "\t")
rownames(metab_names) <- metab_names[,1]

####################################################
### Calculate metabolite exchange Score
# where exchange score == (2 x  ((n_produc * n_cons)/(n_produc+n_cons)))

met_exchanges['exchange_score'] <- 2*((met_exchanges['n_producers']*met_exchanges['n_consumers']) / (met_exchanges['n_producers']+met_exchanges['n_consumers']))
dim(met_exchanges)
write.csv(met_exchanges, "2.1_met_exchanges_exchange_score.csv")



####################################################
### Calculate average exchange score in healthy individuals

metad_healthy <- metad[metad$HD == "healthy",]

#keep only the metabolic exchanges found in healthy samples:
met_exchanges_healthy <- subset(met_exchanges, sample_id %in% metad_healthy$Sample)
dim(met_exchanges_healthy)



########################################################
############## summary and error bars  #################
summary_exchange  <- met_exchanges_healthy %>% 
  group_by(metabolite) %>%   # the grouping variable
  summarise(mean_exchange = mean(exchange_score),  # calculates the mean of each group
            sd_exchange = sd(exchange_score)) # calculates the standard deviation of each group

summ_exchange <- as.data.frame(summary_exchange)
rownames(summ_exchange) <- summ_exchange[,1]

# add metabolite names:
summ_exchange_named <- merge(summ_exchange, metab_names[,c("bigg_id","name_simple")], by="row.names")
summ_exchange_named$bigg_id <- NULL
summ_exchange_named$Row.names <- NULL

# order by metabolite exchange:
summ_exchange_named_order <- summ_exchange_named[order(summ_exchange_named$mean_exchange, decreasing = TRUE),]
rownames(summ_exchange_named_order) <- summ_exchange_named_order$metabolite
write.csv(summ_exchange_named_order, "2.2_summ_exchange_ordered.csv")



########################################################
########### Plot the top xxx metabolites  ###############

## define number of metabolites to plot:
n_wanted = 30

top_metab <- summ_exchange_named_order[1:n_wanted,]

ggplot(top_metab, aes(x=reorder(name_simple, mean_exchange), y=mean_exchange)) + 
  geom_bar(stat="identity", fill='steelblue') + coord_flip() +
  geom_errorbar(aes(ymin=mean_exchange-sd_exchange, ymax=mean_exchange+sd_exchange), width=.2) +
  theme_classic() + labs(y= "Mean exchange Score", x = "Metabolite Exchanged")
  





