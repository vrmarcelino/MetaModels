#### Invetsigating H2S redundancy loss in IBD samples
#### Check ratio of H2S producers to consumers

#### created: 25 - Sep - 2021
#### updated: 13 - Apr - 2022

rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)

####################################################
### Input and parse table
input_parsed_table_fp <- "0_input/metad_and_metab_h2s.csv"
metad_and_metab_study <- read.csv(input_parsed_table_fp)
rownames(metad_and_metab_study) <- metad_and_metab_study[,1]



###########################################
############## summary   #################
# note that this analyses are done prior log2 transformation for easier interpretation
summary_prod_cons  <- metad_and_metab_study %>% # the names of the new data frame and the data frame to be summarised
  group_by(Diagnosis) %>%   # the grouping variable
  summarise(mean_n_prod = mean(n_producers),  # calculates the mean of each group
            sd_n_prod = sd(n_producers), # calculates the standard deviation of each group
            
            mean_n_cons = mean(n_consumers),
            sd_n_cons = sd(n_consumers),
            
            mean_total_prod = mean(h2s_production),
            sd_total_prod = sd(h2s_production),
            
            mean_total_cons = mean(h2s_consumption),
            sd_total_cons = sd(h2s_consumption))

write.table(summary_prod_cons , file = "1.3_summary_prod_cons_prior_log2.csv", sep = ",")



###########################################################################
################# Format tables for plots and analyses  ####################

#convert consumer fluxes to absolute abundances
metad_and_metab_study['flux_consumers'] <- abs(metad_and_metab_study['flux_consumers'])

# get prodcer/consumer ratio (log2 - transformed, but calculate using raw fluxes)
metad_and_metab_study['n_prod_cons_ratio'] <- log2(metad_and_metab_study['n_producers'] / metad_and_metab_study['n_consumers'])
metad_and_metab_study['total_flux_prod_cons_ratio'] <- log2(metad_and_metab_study['h2s_production'] / metad_and_metab_study['h2s_consumption'])
metad_and_metab_study['flux_prod_cons_ratio'] <- log2(metad_and_metab_study['flux_producers'] / metad_and_metab_study['flux_consumers'])


## log2 transform fluxes
metad_and_metab_study['h2s_production'] <- log2(metad_and_metab_study['h2s_production'])
metad_and_metab_study['h2s_consumption'] <- log2(metad_and_metab_study['h2s_consumption'])
metad_and_metab_study['flux_producers'] <- log2(metad_and_metab_study['flux_producers'])
metad_and_metab_study['flux_consumers'] <- log2(metad_and_metab_study['flux_consumers'])


#melt to plot boxplots
metad_and_metab_study_melt <- melt(metad_and_metab_study)

melt_prod_n <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'n_producers',]
melt_cons_n <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'n_consumers',]
melt_total_prod <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'h2s_production',]
melt_total_cons <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'h2s_consumption',]
melt_prod_flux <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'flux_producers',]
melt_cons_flux <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'flux_consumers',]


div_table <- rbind(melt_prod_n,melt_cons_n)
flux_total_table <- rbind(melt_total_prod,melt_total_cons) # takes abundance into consideration
flux_efficiency_table <-rbind(melt_prod_flux,melt_cons_flux) # does not that spp abundances in consideration

names(div_table)[5] <- "category_prod_cons"
names(flux_total_table)[5] <- "category_prod_cons"
names(flux_efficiency_table)[5] <- "category_prod_cons"


####################################################
################### Plot it!! ######################
wanted_colors <- c("#3B9AB2","#E1AF00") #blue and yellow

##### DIVERSITY of producers and consumers and their ratio #####

# Removing the outliers of the box plot with outlier.shape = NA, as the jitter will add them again.
pdiv <- ggplot(div_table, aes(x=category_prod_cons, y=value, fill=Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=wanted_colors) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.15), alpha=0.5) +
  theme_classic() + labs(y = "Number of species", x = "") 


### stats for Producers:
prod_table <- div_table[div_table$category_prod_cons == "n_producers",]
kruskal.test(prod_table$value ~ prod_table$Diagnosis)
#Kruskal-Wallis chi-squared = 21.099, df = 1, p-value = 4.363e-06

### stats for Consumers:
cons_table <- div_table[div_table$category_prod_cons == "n_consumers",]
kruskal.test(cons_table$value ~ cons_table$Diagnosis)
#Kruskal-Wallis chi-squared = 34.153, df = 1, p-value = 5.094e-09


### Ratio of producers and consumers diversity
pdiv_ratio <- ggplot(metad_and_metab_study, aes(x=Diagnosis, y=n_prod_cons_ratio, fill=Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=wanted_colors) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.15), alpha=0.5) +
  theme_classic() + labs(y = "H2S producer to consumer species ratio (log2)", x = "") 

# test (note that log2 transformation does not affect the results)
kruskal.test(metad_and_metab_study$n_prod_cons_ratio ~ metad_and_metab_study$Diagnosis)
# Kruskal-Wallis chi-squared = 19.432, df = 1, p-value = 1.042e-05


##### FLUX of producers and consumers, corrected for their absolute spp abundance #####
pflux_total <- ggplot(flux_total_table, aes(x=category_prod_cons, y=value, fill=Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=wanted_colors) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.15), alpha=0.5) +
  theme_classic() + labs(y = "Total flux of H2S exchanged (log2)", x = "") 


### stats for Producers:
prod_table <- flux_total_table[flux_total_table$category_prod_cons == "h2s_production",]
kruskal.test(prod_table$value ~ prod_table$Diagnosis)
#Kruskal-Wallis chi-squared = 1.8415, df = 1, p-value = 0.1748

### stats for Consumers:
cons_table <- flux_total_table[flux_total_table$category_prod_cons == "h2s_consumption",]
kruskal.test(cons_table$value ~ cons_table$Diagnosis)
#Kruskal-Wallis chi-squared = 18.376, df = 1, p-value = 1.813e-05

### Ratio of producers and consumers flux - w/ ssp. abundances
pflux_total_ratio <- ggplot(metad_and_metab_study, aes(x=Diagnosis, y=total_flux_prod_cons_ratio, fill=Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=wanted_colors) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.15), alpha=0.5) +
  theme_classic() + labs(y = "Total H2S production to consumption ratio (log2)", x = "") 


# test (note that log2 transformation does not affect the results)
kruskal.test(metad_and_metab_study$total_flux_prod_cons_ratio ~ metad_and_metab_study$Diagnosis)
# Kruskal-Wallis chi-squared = 15.922, df = 1, p-value = 6.602e-05


#### Plot/save 4 graphs:
pdf(file = "1.3_sulphur_stats.pdf", width = 8, height = 10)
ggarrange(pdiv, pdiv_ratio, pflux_total, pflux_total_ratio,
          labels = c("A)", "B)","C)","D)"),
          ncol = 2, nrow = 2)

dev.off()



