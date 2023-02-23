#### Script to run linear model and generate vectors as the output

##### Make sure to run the shapeRsq() and shapePRsq() from WRP_FUNCTIONS.R

# Read in data
setwd("/Users/amandan/Desktop/Dworkin/dgrp/data")
covariates <- read.csv(file = "DGRP_covariates_final.csv")
gene_sample <- read.csv(file = "DGRP_genes_final.csv")
shape_coords <- read.csv(file = "DGRP_shape_final.csv")

#Grab the vector for each gene in the dataset:

gene_vectors <- matrix(data = NA, ncol = ncol(shape_coords), nrow = ncol(gene_sample))

for (i in 1:ncol(gene_sample)){
  gene_id <- colnames(gene_sample)[i]
  gene_exp <- gene_sample[,i] 
  model_fit <- lm(as.matrix(shape_coords) ~ 1 + as.matrix(covariates$log_csize_cent) + as.matrix(covariates[,2:5]) + gene_exp + as.matrix(covariates$pop_sub...2.) + as.matrix(covariates$pop_sub...3.) + as.matrix(covariates$wolbachia))
  vec <- coef(model_fit)[7,]
  gene_vectors[i,] <- vec 
} 

gene_vectors <- data.frame(colnames(gene_sample), gene_vectors)
colnames(gene_vectors)[1] <- "gene_id"

write.csv(gene_vectors, "gene_vectors_final.csv", row.names = FALSE)
