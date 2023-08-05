#### Produce a graph showing the difference in metabolite importance in healthy and CD samples 
# focusing on He te al only!!
# only testing it for metabolites present in > half of the samples

#### created: 30 - Sep - 2021
#### updated: 6 - Apr - 2022

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
met_exchanges$sample_id <- gsub('_cat', '', met_exchanges$sample_id) # remove the "_cat" from sample names

metad <- read.csv(metadata_fp)
rownames(metad) <- metad[,1]

metab_names <- read.csv(metab_names_fp, sep = "\t")
rownames(metab_names) <- metab_names[,1]

####################################################
### Calculate metabolite Importance Score

met_exchanges['importance_score'] <- 2*((met_exchanges['n_producers']*met_exchanges['n_consumers']) / (met_exchanges['n_producers']+met_exchanges['n_consumers']))
dim(met_exchanges)
#write.csv(met_exchanges, "2.1_met_exchanges_importance_score.csv")



####################################################
### Calculate average Importance score in healthy and CD from He et al cohhort:

metad_He2017 <- metad[metad$Author_Year == "He_2017",]
metad_He2017_healthy <- metad_He2017[metad_He2017$HD == "healthy",]

#keep only the metabolic exchanges found in healthy samples:
met_exchanges_healthy <- subset(met_exchanges, sample_id %in% metad_He2017_healthy$Sample)
dim(met_exchanges_healthy)

metad_He2017_CD <- metad_He2017[metad_He2017$HD == "diseased",]
met_exchanges_CD <- subset(met_exchanges, sample_id %in% metad_He2017_CD$Sample)

# total of 84 samples (46 disease and 38 healthy)


########################################################
############## summary and error bars  #################
summary_importance_healthy  <- met_exchanges_healthy %>% 
  group_by(metabolite) %>%   # the grouping variable
  summarise(mean_importance_healthy = mean(importance_score),  # calculates the mean of each group
            sd_importance_healthy = sd(importance_score)) # calculates the standard deviation of each group

summ_importance_healthy <- as.data.frame(summary_importance_healthy)
rownames(summ_importance_healthy) <- summ_importance_healthy[,1]



summary_importance_CD  <- met_exchanges_CD %>% 
  group_by(metabolite) %>%   # the grouping variable
  summarise(mean_importance_CD = mean(importance_score),  # calculates the mean of each group
            sd_importance_CD = sd(importance_score)) # calculates the standard deviation of each group

summ_importance_CD <- as.data.frame(summary_importance_CD)
rownames(summ_importance_CD) <- summ_importance_CD[,1]


# merge tables
summ_importance <- merge(summ_importance_healthy, summ_importance_CD, by='row.names', all=T) #
summ_importance[is.na(summ_importance)] <- 0 
rownames(summ_importance) <- summ_importance[,1]
summ_importance$Row.names <- NULL
summ_importance$metabolite.x <- NULL
summ_importance$metabolite.y <- NULL


# add metabolite names:
summ_importance_named <- merge(summ_importance, metab_names[,c("bigg_id","name_simple")], by="row.names")
rownames(summ_importance_named) <- summ_importance_named[,1]
summ_importance_named$bigg_id <- NULL
summ_importance_named$Row.names <- NULL


# calculate difference in Metabolite Importance Scores:
summ_importance_named['diff_MIS'] <- summ_importance_named['mean_importance_healthy'] - summ_importance_named['mean_importance_CD']


######################################################
#### Test if the differences are significant  ########
######################################################

# merge healthy and CD dataframes
met_exchanges_healthy['HD'] = 'healthy'
met_exchanges_CD['HD'] = 'CD'
He_met_exchanges <- rbind(met_exchanges_healthy, met_exchanges_CD)

# get a list of unique metabolites
unique_metabolites <- unique(He_met_exchanges$metabolite)

summ_importance_named['p_value'] = "NaN" # column to store p-values
count_tests = 0
for (met in unique_metabolites){
  wanted_met_df <-  He_met_exchanges[He_met_exchanges$metabolite == met,]
  
  # only test if metabolite is present in > half of the samples (42 samples):
  if (dim(wanted_met_df)[1] > 42) {
    wanted_abundant_met_df <- wanted_met_df
    count_tests = count_tests + 1 # keep track of how many tests are being performed
    kw <- kruskal.test(wanted_abundant_met_df$importance_score ~ wanted_abundant_met_df$HD)
    summ_importance_named[met, 'p_value'] <- kw$p.value
  }
}

# calculate the bonferroni-corrected p-value:
bon_p <- 0.05 /count_tests
bon_p 

# add significance if kruskal-wallis p is < than bin_p:
summ_importance_named$sig_after_bonferroni <- "ns"
summ_importance_named$sig_after_bonferroni <- ifelse(as.numeric(summ_importance_named$p_value) < bon_p, "*", summ_importance_named$sig_after_bonferroni)


# order by metabolite importance difference:
summ_importance_named_order <- summ_importance_named[order(summ_importance_named$diff_MIS, decreasing = TRUE),]

# remove NAs (generated when importance was zero)
summ_importance_named_order <- na.omit(summ_importance_named_order)

# save it
write.csv(summ_importance_named_order, "1.2_mean_importance_diff_ordered.csv")

