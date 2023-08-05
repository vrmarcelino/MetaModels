#### Produce a beta diversity plot highlighting healthy and diseased samples

#### created: 13 - Jun - 2023

rm(list=ls(all=TRUE))

library(mixOmics)
library(dplyr)

####################################################
### Input and parse tables

metadata_fp <- "0_input_data/metadata_rewiring_microbiome.csv"
spp_composition_fp <- "0_input_data/merged_kma_res.csv"

metad <- read.csv(metadata_fp)
rownames(metad) <- metad[,1]

spp_conposition <- read.csv(spp_composition_fp, header=T)
names(spp_conposition) <- gsub('_cat', '', names(spp_conposition)) # remove the "_cat" from sample names
rownames(spp_conposition) <- spp_conposition[,1]

spp_names <- read.csv("0_input_data/wanted_spp_classification_simp.tsv", sep='\t')
rownames(spp_names) <- spp_names[,1]


# select the samples from HE et al (balanced dataset)
CD_samples <- read.table('0_input_data/He_CD_prefixes_rarefied.txt', header = F)
Healthy_samples <- read.table('0_input_data/He_healthy_prefixes_all_that_worked.txt', header = F)
wanted_samples <- rbind(CD_samples,Healthy_samples)

metad_wanted_samples <- metad[wanted_samples$V1,]
spp_wanted_samples <- spp_conposition[,wanted_samples$V1]
spp_wanted_samples <- t(spp_wanted_samples)

# remove MAGs with zero counts
spp_wanted_samples <- spp_wanted_samples[,colSums(spp_wanted_samples[])>0]
dim(spp_wanted_samples) # 617 species remained


###################
## pre-process the data:

# add in value to remove zeroes:
min_depth <- min(spp_wanted_samples[spp_wanted_samples > 0])
data.offset <- spp_wanted_samples + min_depth

# filtering is not needed here (KMA is quite restrictive and we are already working with MAGs)


###################
## do CLR and PCA:

pca.result <- pca(data.offset, logratio = 'CLR') 


## add diagnosis:
samples_ordered <- pca.result$names$sample
diag <- cbind(samples_ordered,metad_wanted_samples[,'Diagnosis'][match(samples_ordered, rownames(metad_wanted_samples))])
pca.result$diag <- diag[,2]


## Plot,  color by abundance in healthy vs disease:
pdf(file = "2.1_PCA_BetaDiv.pdf", width = 6, height = 4)

plotIndiv(pca.result,  # plot samples
          group = pca.result$diag,
          title = 'PCA Comps 1&2',
          pch = 16, cex = 4,
          col = c("#3B9AB2","#E1AF00"),
          legend = TRUE)

dev.off()





