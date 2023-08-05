#### Test if different the microbiome associated with diseases have a significantly different alpha diversity

rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan) # shannon div
library("wesanderson")

####################################################
### Input and parse tables

metadata_fp <- "0_input_data/metadata_rewiring_microbiome.csv"
spp_composition_fp <- "0_input_data/merged_kma_res.csv"

metad <- read.csv(metadata_fp)
rownames(metad) <- metad[,1]

spp_conposition <- read.csv(spp_composition_fp)
names(spp_conposition) <- gsub('_cat', '', names(spp_conposition)) # remove the "_cat" from sample names

samples_that_worked <- read.csv("0_input_data/samples_that_worked_w_MICOM.tsv", sep = '\t')
rownames(samples_that_worked) <- gsub('_cat', '', samples_that_worked[,1]) # remove the "_cat" from sample names
wanted_samples <- rownames(samples_that_worked)



####################################################
### Calculate species diversity

#species richness
spp_richness_pre <- spp_conposition
spp_richness_pre["numb_species",] <- colSums(spp_richness_pre != 0) # calculate the number of spp in each sample
spp_richness <- t(spp_richness_pre["numb_species",])
mode(spp_richness) = "numeric"


# shannon (info about both richness and evenness)
shannon_pre <- spp_conposition[1:ncol(spp_conposition)-1] #remove tax
rownames(shannon_pre) <- shannon_pre[,1]
shannon_pre <- shannon_pre[-1] #remove first column
shannon_pre <- t(shannon_pre)
dim(shannon_pre) # 1697, 955

shannon_div <- as.data.frame(diversity(shannon_pre, index = "shannon"))
names(shannon_div) <- "shannon"

# keep only the wanted samples (rmeov ethe ones that failed MICOM optimization):
shannon_div <- shannon_div %>% filter(row.names(shannon_div) %in% wanted_samples)

dim(shannon_div)
# 1661 samples

####################################################
### Merge with Metadata:
metad_and_div_pre <- cbind(spp_richness, metad[,"Diagnosis"][match(rownames(spp_richness), rownames(metad))])
metad_and_div <- merge(metad_and_div_pre, shannon_div, by="row.names")
names(metad_and_div) <- c("sample", "spp_richness", "Diagnosis", "shannon")
rownames(metad_and_div) <- metad_and_div[,1]
metad_and_div$spp_richness <- as.numeric(metad_and_div$spp_richness)


####################################################
###  Wilcox stats 
# note = Default for p.adjust.method = "holm".
richness_stats <- compare_means(spp_richness ~ Diagnosis, metad_and_div, ref.group = "healthy")
write.csv(richness_stats,file = "1.1_wilcox_stats_healthy_vs_others_richness.csv")

shannon_stats <- compare_means(shannon ~ Diagnosis, metad_and_div, ref.group = "healthy")
#write.csv(shannon_stats,file = "1.2_wilcox_stats_healthy_vs_others_shannon.csv")


####################################################
### plot box plots with jitter dots
### Figure #: comparison across disease categories

## order by species div median (first calculate the mean Spp div per disease, then order based on it)
summary_Diagnosis  <- metad_and_div %>% # the names of the new data frame and the data frame to be summarised
  group_by(Diagnosis) %>%   # the grouping variable
  summarise(median_richness = median(spp_richness),  # calculates the median of each group
            sd_richness = sd(spp_richness)) # calculates the standard deviation of each group

summary_Diagnosis_order <- as.data.frame(summary_Diagnosis[order(summary_Diagnosis$median_richness),])

#convert to factor with specified order
metad_and_div$Diagnosis <- factor(metad_and_div$Diagnosis,levels = c(summary_Diagnosis_order$Diagnosis))

# colors:
wanted_colors <- wes_palette("Zissou1", 5, type = "discrete")
wanted_colors <- c(wanted_colors[4],wanted_colors[4],wanted_colors[4],wanted_colors[4],
                   wanted_colors[4],wanted_colors[4],wanted_colors[4],wanted_colors[1],
                   wanted_colors[4],wanted_colors[4],wanted_colors[4],wanted_colors[4])



# plot it:
shan_p <- ggboxplot(metad_and_div, x="Diagnosis", y="shannon") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=.9, aes(colour=Diagnosis), show.legend = FALSE) +
  coord_flip() + labs(y = "Shannon Ind", x = "Diagnosis") +
  scale_color_manual(values=wanted_colors) + theme_classic(base_size =18)

rich_p <- ggboxplot(metad_and_div, x="Diagnosis", y="spp_richness") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=.9, aes(colour=Diagnosis), show.legend = FALSE) +
  coord_flip() + labs(y = "Species richness", x = "Diagnosis") + 
  scale_color_manual(values=wanted_colors) + theme_classic(base_size =18)

## save it
pdf(file = "1.3_alpha_diversity_per_disease.pdf", width = 10, height = 6)
ggarrange(rich_p, shan_p,
          labels = c("A)", "B)"),
          ncol = 2, nrow = 1)
dev.off()

