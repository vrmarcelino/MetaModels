#### Invetsigating H2S redundancy loss in IBD samples
#### Check if exchanges with media were also affected

#### created: 03 - Jan 203

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
            sd_total_cons = sd(h2s_consumption),
            
            mean_h2s_in_media = mean(h2s_in_media),
            sd_h2s_in_media = sd(h2s_in_media))

#write.table(summary_prod_cons , file = "1.3.2_summary_prod_cons_prior_log2_w_h2s.csv", sep = ",")



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

melt_h2s_in_media <- metad_and_metab_study_melt[metad_and_metab_study_melt$variable == 'h2s_in_media',]



### stats for H2S in media:
kruskal.test(melt_h2s_in_media$value ~ melt_h2s_in_media$Diagnosis)


