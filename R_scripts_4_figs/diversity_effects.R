#### Test if the loss of species diversity affects H2S producers and consumers differently

## created: 28 - March - 2022
## updated: 13 - Apr - 2022

rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(lsmeans)


####################################################
### Input and parse tables


# This table was generated with 0_make_table.R script, includes the total + net production & consumption of H2S
input_parsed_table_fp <- "0_input/metad_and_metab_h2s.csv"
metad_and_metab_study <- read.csv(input_parsed_table_fp)
rownames(metad_and_metab_study) <- metad_and_metab_study[,1]

spp_composition_fp <- "0_input/merged_kma_res.csv"
spp_conposition <- read.csv(spp_composition_fp)
names(spp_conposition) <- gsub('_cat', '', names(spp_conposition)) # remove the "_cat" from sample names
spp_conposition["spp_diversity",] <- colSums(spp_conposition != 0) # calculate the number of spp in each sample
spp_div <- t(spp_conposition["spp_diversity",])
mode(spp_div) = "numeric"


metad_and_metab <- cbind(metad_and_metab_study,spp_div[match(rownames(metad_and_metab_study), rownames(spp_div))])
names(metad_and_metab)[length(metad_and_metab)] <- "Species_diversity"


#####################################################
# convert total h2s production and consumption to log2
metad_and_metab['h2s_production'] <-log2(metad_and_metab['h2s_production'])
metad_and_metab['h2s_consumption'] <-log2(metad_and_metab['h2s_consumption']) # already in absolute terms



################################################
############## Test and Plot!  #################

###### NUMBER of producers and consumers ######

#format table
metad_and_metab_study_melt<-melt(metad_and_metab,measure.vars = c("n_producers", "n_consumers"))
names(metad_and_metab_study_melt)[13] <- "category_prod_cons"
names(metad_and_metab_study_melt)[14] <- "number_prod_or_cons"

# visualize
wanted_colors <- c("#d95f02","#009E73")

p1 <- ggplot(metad_and_metab_study_melt, aes(Species_diversity, number_prod_or_cons, color = category_prod_cons)) +
  geom_point() + scale_color_manual(values =wanted_colors) +
  geom_smooth(method = "lm", se = T) +
  ggpubr::theme_pubr() + labs(y = "Number of producers or consumers", x = "Species diversity") 



m.interaction <- lm(number_prod_or_cons ~ Species_diversity*category_prod_cons, data = metad_and_metab_study_melt)
summary(m.interaction)

# Obtain slopes
m.lst <- lstrends(m.interaction, "category_prod_cons", var="Species_diversity")
m.lst


# Compare slopes (compare the Estimated marginal means)
pairs(m.lst)

############################################################################
###### FLUX of producers and consumers, accounting for their abundance ######

#format table
metad_and_metab_study_melt<-melt(metad_and_metab,measure.vars = c("h2s_production", "h2s_consumption"))
names(metad_and_metab_study_melt)[13] <- "category_prod_cons"
names(metad_and_metab_study_melt)[14] <- "flux_prod_or_cons_abund"

# visualize
p2 <- ggplot(metad_and_metab_study_melt, aes(Species_diversity, flux_prod_or_cons_abund, color = category_prod_cons)) +
  geom_point() + scale_color_manual(values =wanted_colors) +
  geom_smooth(method = "lm", se = T) +
  ggpubr::theme_pubr() +
  labs(y = "Total H2S production or consumption (log2)",x = "Species diversity")


m.interaction <- lm(flux_prod_or_cons_abund ~ Species_diversity*category_prod_cons, data = metad_and_metab_study_melt)
summary(m.interaction)
m.lst <- lstrends(m.interaction, "category_prod_cons", var="Species_diversity")
m.lst
pairs(m.lst)



#### Plot/save both graphs:
pdf(file = "1.4_diversity_effect.pdf", width = 8, height = 5)
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()


