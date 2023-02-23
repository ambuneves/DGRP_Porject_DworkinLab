################################################################################
# Script to compute gene expression means for each line and to compute shape   #
# means for each line. Write these as their own .csv files                     #
# The files created in this script are:                                        #
# dgrp_shape_means.csv and dgrp_gene_means.csv                                 #
################################################################################

#### Libraries
library(geomorph)
library(ggfortify)

##### Reading in data ####
### set up the working directory
setwd("/Users/amandan/Desktop/Dworkin/dgrp/data")

### Load in metadata
sample_info <- read.csv("dgrp-sample_info.csv")

#order sample_info so that it matches quant_files
sample_info <- sample_info[order(sample_info$fastqFile),]

#add a column for "Line" so that lines can be matched with the shape data Line column
sample_info$LineNum <- gsub("_|_1|_2|-1|-2|_T1|T2", "", x = sample_info$lineName)

### Load in gene expression data
cts <- read.csv("lengthscaled_counts.csv", row.names = 1)

### Load in shape data
raw_data <- read.csv("BothLabs_Wings_28Oct.csv")

plotAllSpecimens(arrayspecs(shape_data[,10:105], p = 48, k = 2)) #plot the landmark coordinates. Does this look like a wing to you? (yes)

#### Wrangling the data ####

# Subset data to only the Dworkin lab. Do not include the Houle lab wings. Both 
# labs took same landmarks but rearing conditions and wing imaging process differed,
# this showed to make size vary a lot between the two labs. This will screw up
# any sort of allometric relationships. To avoid this, use only the wings collected
# by the dworkin lab. 

#Only the dworkin lab
shape_data <- raw_data[grep("Dwo", raw_data$Lab),]
dim(shape_data) #15749

#Combine $Line and $Ind into one column for easier identification
shape_data$ID <- paste(shape_data$Line, shape_data$Ind, sep = "_")

# Multiply the centroid size by the scale. This won't give us different results
# for the analysis but it will make results in the proper scale (a more accurate
# measure) 

shape_data$Csize <- shape_data$Csize * shape_data$Scale

#double check that all the scale is the same
unique(shape_data$Scale)

# Checking how many replicates we have for each sample, might want to 'take out'(?) samples with lower reps
table(shape_data$Rep)
table(shape_data$LineNum, shape_data$Rep)  #Why would some lines have no value for rep1 or rep2 but have values for rep3 and/or rep4? -> because blocked design
table(shape_data$LineNum, shape_data$Sex) 
table(shape_data[grep("256|306|887|757", shape_data$LineNum),]$LineNum, 
      shape_data[grep("256|306|887|757", shape_data$LineNum),]$Sex) #lines 256,306,887 are imbalanced between males and females and 
#have low rep numbers so I'll remove those... should I remove 757? Will for now. Not even sure if we have rna for it anyways 

#256 is in the RNA sample, 306,757, and 887 are not 

#Is there a way to quantify the number of males and females in a sample compared to the whole and
#identify the outliers? 

summary(as.data.frame(table(shape_data$LineNum, shape_data$Sex)))

#### Figuring out if it would be better to compute a M/F average and then average those,
# or take the mean shape regardless of sex 

#line means to compare 
line_means <- aggregate(shape_data[,10:106], by = list(shape_data$LineNum), FUN = "mean")
colnames(line_means)[1] <- "LineNum" #rename the first column
dim(line_means)
#write.csv(line_means, "dgrp_shape_csize_means.csv")

#include centroid size where it is an average of the m and f averages 
c_means <- aggregate(shape_data[,10:106], by = list(shape_data$LineNum, shape_data$Sex), FUN = "mean")
dim(c_means)
csize_sex_avg <- aggregate(c_means[,3:99], by = list(c_means$Group.1), FUN = "mean")
colnames(csize_sex_avg)[1] <- "LineNum"
#write.csv(x = csize_sex_avg, "dgrp_shape_csize_mf-avg.csv")

#first, sex and line
sex_line_m <- aggregate(shape_data[,10:105], by = list(shape_data$LineNum, shape_data$Sex), FUN = "mean")
colnames(sex_line_m)[1] <- "LineNum"
colnames(sex_line_m)[2] <- "Sex"

#then, average the sexes 
sex_avg_means <- aggregate(sex_line_m[,3:98], by = list(sex_line_m$LineNum), FUN = "mean")
colnames(sex_avg_means)[1] <- "LineNum"

#One which includes centroid size, cuz that would be helpful
sex_line_csize <- aggregate(shape_data[,10:106], by = list(shape_data$LineNum, shape_data$Sex), FUN = "mean")
colnames(sex_line_csize)[1] <- "LineNum"
colnames(sex_line_csize)[2] <- "Sex"
#Write this one to a .csv file
#write.csv(sex_line_csize, "sex_line_csize_shape_means.csv")

#create a new matrix with both of these, plus a value for combined or not
df_1 <- cbind.data.frame(rep("combined_mean", times = length(line_means$LineNum)),line_means)
colnames(df_1)[1] <- "method"
df_2 <- cbind.data.frame(rep("averaged_mean", times = length(sex_avg_means$LineNum)), sex_avg_means)
colnames(df_2)[1] <- "method"
check_means <- rbind.data.frame(df_1, df_2)

autoplot(prcomp(check_means[,3:98]), data = check_means, colour = "method", label = FALSE, alpha = 0.7) #makes sense that sexes separate

check_means[check_means$LineNum=="256",] #In the case with 20M and 0F, the male means were used 
# (combined_mean and averaged_mean are the same)

# Removing the "problematic" samples --> Not going to do this for now 
# dim(shape_data) #15749
# shape_data <- shape_data[-grep("256|306|887|757", shape_data$LineNum),]
# dim(shape_data) #15673

# mean per line, per rep
rep_mean <- aggregate(shape_data[,10:105], by = list(shape_data$LineNum, shape_data$Sex, shape_data$Rep), FUN = "mean")
colnames(rep_mean)[1:3] <- c("LineNum", "Sex", "Rep")


#### Plot PCAs to compare the different mean groups to see which are reasonable to use for the analysis 
#using geomorph prcomp

autoplot(prcomp(sex_line_m[,3:98]), data = sex_line_m, colour = "Sex", label = FALSE) #makes sense that sexes separate
autoplot(prcomp(rep_mean[,4:99]), data = rep_mean, colour = "Rep", label = FALSE) #would be bad if replicates separated

# Plot PCA with a distinct colour combination for each DGRP (rep_mean)
# Looks... not helpful
autoplot(prcomp(rep_mean[,4:99]), data = rep_mean, colour = "LineNum", label = FALSE) #would be bad if replicates separated

# To model the sources of variation in the data, I'm going to create a linear model that includes 
# sex, rep, idk as terms in the model... maybe I won't include gene expression (for now)
# use the rep_mean dataset to do this 
table(rep_mean$Sex) #pretty even between the sexes (301F and 302M)
table(rep_mean$Rep) #most have 3 or 4 reps, some have only 1 or two
lm.test <- lm(as.matrix(rep_mean[,4:99]) ~ 1 + Rep + Sex, data = rep_mean) #See a model for each coordinate... but don't I want the whole thing?
summary(lm.test) #model for each coord... not exactly what I want? or is it?

coords <- rep_mean[,4:99]
gdf <- geomorph.data.frame(rep_mean[,4:99], line = rep_mean$LineNum, sex = rep_mean$Sex, rep = rep_mean$Rep)
fit1 <- procD.lm(coords ~ sex + rep, data = gdf) #sex and replicate both seem to have significant effects though R^2 is not large especially for rep
summary(fit1)
fit2 <- procD.lm(coords ~ rep + sex, data = gdf, iter = 99, SS.type = 2) #type II sum of squares
fit(fit2)

# Since the males and females are pretty even, is it alright if I just go ahead with using the mixed line means?
# Larval sample is mixed anyways..... *****

#Write line means to a .csv file
#write.csv(line_means, "dgrp_shape_means.csv")

#### Mean gene expression per gene, per line ####

genes_agg <- cbind.data.frame(data.frame(sample_info)$LineNum, t(cts)) #combine transpose of cts with LineNum identifier 
colnames(genes_agg)[1] <- "LineNum"
dim(genes_agg) #108, 13702
gene_means <- aggregate(genes_agg[,2:13702], by = list(genes_agg$LineNum), FUN = mean) #aggregate by mean 
dim(gene_means) #94, 13702 ... where [,2:13702] are genes and col 1 is the line names 
colnames(gene_means)[1] <- "LineNum"

gene_means[1:5,1:5] #Just a look at the first five genes, first five lines

#write as a .csv
#write.csv(gene_means, "dgrp_gene_abd_means.csv")
write.csv(gene_means, "dgrp_gene_cts_means.csv")

#### Data set matching #### *** Maybe do this after the gene means 
#Since the exact lines measured for the RNA seq and shape analysis might not be the same,
#I want to subset both datasets so that only individuals from matching lines are kept (otherwise, they are not really comparable)

gene_sample <- sample_info[sample_info$LineNum %in% line_means$LineNum,]
dim(gene_sample) #only contains the 82 matching lines (out of 149)

shape_sample <- line_means[line_means$LineNum %in% gene_sample$LineNum,]
dim(shape_sample) #now contains 82 shape measurements, lower than above because i didn't deal with the gene reps yet 
