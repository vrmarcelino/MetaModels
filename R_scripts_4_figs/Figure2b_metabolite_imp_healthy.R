#### Produce a graph of average metabolite exchange score in healthy individuals

#### created: 31 - Jul - 2023

rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
library(dplyr)
library(unikn) #colors
library(MetBrewer)#colors

####################################################
### Input and parse tables

input_exchanges_fp <- "0_input/producers_consumers.csv"
metadata_fp <- "0_input/metadata_rewiring_microbiome.csv"
metab_names_fp <- "0_input/bigg_models_w_classes_curated_simplified_names.tsv"

met_exchanges <- read.csv(input_exchanges_fp)
met_exchanges$sample_id <- gsub('_cat', '', met_exchanges$sample_id) # remove the "_cat" from sample names

metad <- read.csv(metadata_fp)
rownames(metad) <- metad[,1]

metab_names <- read.csv(metab_names_fp, sep = "\t")
rownames(metab_names) <- metab_names[,1]

####################################################
### Calculate metabolite Importance Score
# where importance score == (2 x  ((n_produc * n_cons)/(n_produc+n_cons)))

met_exchanges['importance_score'] <- 2*((met_exchanges['n_producers']*met_exchanges['n_consumers']) / (met_exchanges['n_producers']+met_exchanges['n_consumers']))
dim(met_exchanges)
#write.csv(met_exchanges, "2.1_met_exchanges_importance_score.csv")



####################################################
### Calculate average Importance score in healthy individuals

metad_healthy <- metad[metad$HD == "healthy",]

#keep only the metabolic exchanges found in healthy samples:
met_exchanges_healthy <- subset(met_exchanges, sample_id %in% metad_healthy$Sample)
dim(met_exchanges_healthy)



########################################################
############## summary and error bars  #################
summary_importance  <- met_exchanges_healthy %>% 
  group_by(metabolite) %>%   # the grouping variable
  summarise(mean_importance = mean(importance_score),  # calculates the mean of each group
            sd_importance = sd(importance_score)) # calculates the standard deviation of each group

summ_importance <- as.data.frame(summary_importance)
rownames(summ_importance) <- summ_importance[,1]

# add metabolite names:
summ_importance_named <- merge(summ_importance, metab_names[,c("bigg_id","name_simple", "SubClass")], by="row.names")
summ_importance_named$bigg_id <- NULL
summ_importance_named$Row.names <- NULL

# order by metabolite importance:
summ_importance_named_order <- summ_importance_named[order(summ_importance_named$mean_importance, decreasing = TRUE),]
rownames(summ_importance_named_order) <- summ_importance_named_order$metabolite
#write.csv(summ_importance_named_order, "2.2_summ_importance_ordered.csv")

# Remove O2 and H2O to improve visualization
summ_importance_named_order <-summ_importance_named_order[!summ_importance_named_order$metabolite == "o2_e", ]
summ_importance_named_order <- summ_importance_named_order[!summ_importance_named_order$metabolite == "h2o_e", ]


########################################################
########### Select the top xxx metabolites  ############

## define number of metabolites to plot:
n_wanted = 15

top_metab <- summ_importance_named_order[1:n_wanted,]




##############################################################################################
########### Get the whole dataset for the top X metabolites (to plot all samples)  ############


top_metabs <- met_exchanges_healthy[met_exchanges_healthy$metabolite %in% top_metab$metabolite,]
metab_names$metabolite <- rownames(metab_names)

top_metabs_named <- merge(top_metabs, metab_names[,c("metabolite","bigg_id","name_simple", "SubClass")], by="metabolite")


#convert to factor with specified order
top_metabs_named$name_simple<- factor(top_metabs_named$name_simple,levels = c(rev(top_metab$name_simple)))


#colors - same as the ones used for panel 2d:
# n_colors <- length(unique(top_metabs_named$SubClass))
#wanted_colors <- met.brewer("Juarez", n=9)


#match to the subclasses used here:
unique(top_metabs_named$SubClass)
unique(top_metabs_named$name_simple)

top_metabs_named$SubClass <- factor(top_metabs_named$SubClass, levels = c(
  "Alcohols and polyols", "Amines", "Amino acids, peptides, and analogues",
  "Beta hydroxy acids and derivatives", "Carbohydrates and carbohydrate conjugates",
  "Carbonyl compounds", "Fatty acids and conjugates", "Non-metal phosphates",
  "Other non-metal sulfides", "Purine 2'-deoxyribonucleosides", "Pyrimidines and pyrimidine derivatives", 
  "Quaternary ammonium salts","Tricarboxylic acids and derivatives", "no sub class"))


wanted_colors_all <- met.brewer("Juarez", n=14)
wanted_colors <- c(wanted_colors_all[2],
                   wanted_colors_all[4],
                   wanted_colors_all[5],
                   wanted_colors_all[6],
                   wanted_colors_all[8],
                   wanted_colors_all[9],
                   wanted_colors_all[11],
                   wanted_colors_all[12],
                   wanted_colors_all[14])



##############################################################################
#### plot and save:


pdf(file = "MESSI_healthy.pdf", width = 5, height = 10)

top_metabs_named %>% ggplot(aes(x=name_simple, y=importance_score))+
  geom_violin(show.legend = F, aes(fill=SubClass)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=.4,show.legend = FALSE, size = 0.3) +
  coord_flip() + labs(y = "MESSI", x = "Metabolite") +
  scale_fill_manual(values=wanted_colors) + theme_classic(base_size =18)

dev.off()


## Legend:
pdf(file = "MESSI_healthy_legend.pdf", width = 4, height = 4)

top_metabs_named %>% ggplot(aes(x=name_simple, y=importance_score))+
  geom_violin(show.legend = T, aes(fill=SubClass)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=.4,show.legend = F, size = 0.3) +
  coord_flip() +
  scale_fill_manual(values=wanted_colors) + theme_classic(base_size =18)

dev.off()



