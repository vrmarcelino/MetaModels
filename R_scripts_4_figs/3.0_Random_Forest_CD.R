#### Test if we can identify the same 'potentially therapeutic' bacteria for CD
# using Random Forest

rm(list=ls(all=TRUE))

library(randomForest)
library(caret)
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

merged_df <- cbind(spp_wanted_samples,metad_wanted_samples[,'Diagnosis'][match(rownames(spp_wanted_samples), rownames(metad_wanted_samples))])
colnames(merged_df)[ncol(merged_df)] <- "Diagnosis"
merged_df<- as.data.frame(merged_df)

# Diagnosis should be a factor variable.
merged_df$Diagnosis <- as.factor(merged_df$Diagnosis)
table(merged_df$Diagnosis) #sample number of samples in both categories

# convert abundances to numeric:
sapply(merged_df, class) # check that they are character
cols.name <- colnames(merged_df)[1:(ncol(merged_df)-1)] # all MAG names, minus Diagnosis
merged_df[cols.name] <- sapply(merged_df[cols.name],as.numeric) 
sapply(merged_df, class)



###################################################
### Random Forest

# Partition into Train and Test data
set.seed(200)
ind <- sample(2, nrow(merged_df), replace = TRUE, prob = c(0.7, 0.3))
train <- merged_df[ind==1,]
test <- merged_df[ind==2,]


# Random Forest 
rf <- randomForest(Diagnosis~., data=train, proximity=TRUE, importance=TRUE)
print(rf)


# Prediction & Confusion Matrix – train data
p1 <- predict(rf, train)
confusionMatrix(p1, train$Diagnosis)
# Train data accuracy is 100% that indicates all the values classified correctly.


# Prediction & Confusion Matrix – test data
p2 <- predict(rf, test)
confusionMatrix(p2, test$Diagnosis)

# Error rate of Random Forest
plot(rf)


####################################################
### Calculate the most 'important' variables and plot it:
importance(rf)

# make dataframe from importance() output
feat_imp_df <- importance(rf) %>% 
  data.frame() %>% 
  mutate(feature = row.names(.)) 

# add spp names
feat_imp_df_named <- cbind(feat_imp_df, spp_names[,"spp_name"][match(rownames(feat_imp_df), rownames(spp_names))])

# keep only the top 30
feat_imp_df_ordered <- feat_imp_df_named[order(feat_imp_df_named$MeanDecreaseGini,decreasing=TRUE),]
names(feat_imp_df_ordered)[6] <- "species"
feat_imp_top <- feat_imp_df_ordered[1:30,]


# plot dataframe
ggplot(feat_imp_top, aes(x = reorder(species, MeanDecreaseGini), 
                        y = MeanDecreaseGini)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_classic() +
  labs(
    x     = "Feature",
    y     = "Importance",
    title = "Feature Importance: <Model>"
  )


# Color by whether the abundance is increase or decreasing in CD


## calculate the mean abundance of each MAG in health and disease
## also cal their prevalence
merged_df
agg_df_mean <- aggregate(merged_df[cols.name], by=list(merged_df$Diagnosis), FUN=mean)
agg_df_mean_t <- as.data.frame(t(agg_df_mean))
names(agg_df_mean_t) <- agg_df_mean_t[1,]
agg_df_mean_t <- agg_df_mean_t[-1,] #remove first row


# test if mean abundance in healthy is greater than in IBD
agg_df_mean_t <- agg_df_mean_t %>%
  mutate(more_abund_in = if_else(healthy > IBD, 'healthy', 'ibd'))


## merge with the plotting df:
feat_imp_top_n_abund <- merge(feat_imp_top, agg_df_mean_t['more_abund_in'], by="row.names")

wanted_colors <- c("#3B9AB2","#E1AF00") #blue and yellow


# plot dataframe, color by abundance in healthy vs disease:
pdf(file = "3.1_Random_forest.pdf", width = 8, height = 8)
ggplot(feat_imp_top_n_abund, aes(x = reorder(species, MeanDecreaseGini), 
                         y = MeanDecreaseGini, fill=more_abund_in)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=wanted_colors)+
  theme_classic() +
  labs(
    x     = "Feature",
    y     = "Importance",
    title = "Feature Importance: <Model>"
  )
dev.off()



