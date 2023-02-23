# This is the script that calculates the vector correlations for the random (naive) GO term group
# and the correlations for the matched gene expression correlation group (these groups were composed in a different script,
# called compose_gene-exp-corr_groups.R

#### First, calculate the actual correlations of vectors in the GO term groups ####

#### "True" pairwise gene correlations 

# Initialize an empty matrix to put the values in
pathway_mat <- matrix(data = NA, nrow = 9, ncol = 4)
pathway_mat <- as.data.frame(pathway_mat)
colnames(pathway_mat) <- c("pathway", "path_length", "mean_cor", "cor_sd")

# Do the loop!

for (k in 1:length(my_pathways)){ # For each pathway, create an nxn matrix where n represents each gene in that pathway
  pathway1 <- my_pathways[[k]]
  pathway_mat$pathway[k] <- names(my_pathways)[k]
  pathway_mat$path_length[k] <- nrow(pathway1)
  pathway_comp <- matrix(data = NA, ncol = nrow(pathway1), nrow = nrow(pathway1))
  for (i in 1:nrow(pathway1)){ # For each gene in the pathway, get the vector of expression
    gene_i <- which(gene_ids == pathway1$genes[i]) # index the gene in the expression matrix 
    vec1 <- gene_vecs[gene_i,] - mean(gene_vecs[gene_i,]) #the first vector, mean center 
    for (j in 1:nrow(pathway1)){ # then for each gene in the pathway (vec 2), compute the correlation with vec 1
      gene_j <- which(gene_ids == pathway1$genes[j])
      vec2 <- gene_vecs[gene_j,] - mean(gene_vecs[gene_j,])
      #vec.cor <- abs((t(vec1) %*% vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2))))
      vec.cor <- cor(vec1,vec2)
      pathway_comp[i,j] <- vec.cor
    }
  }
  
  pathway_mat$mean_cor[k] <- round(mean(pathway_comp[upper.tri(pathway_comp)]), 2)
  pathway_mat$cor_sd[k] <- round(sd(pathway_comp[upper.tri(pathway_comp)]), 2)
}

all_pathways <- data.frame(pathway_mat)

#write.csv(all_pathways, "/Users/amandan/Desktop/Dworkin/dgrp/data/true_pathway_correlation_values.csv", row.names = FALSE)


#### Second, calculate correlations between random groups of genes matched in length to the real GO term groups ####
#### "Naive" Permutation Values 
# Making the assumption that the null is a "pathway" of random genes matching the various lengths of the "real" pathways we used above

datalist = list() # initiate a datalist 

iterations <- 1000 # set the number of iterations

for (n in 1: length(my_pathways)){
   pathway1 <- my_pathways[[n]]
   number_of_genes <- nrow(pathway1)
   gene_vecs_subset <- gene_vecs[-which(gene_ids %in% pathway1$genes),] # remove the genes that are in one of our known pathways 
   pathway_comp_iter <- matrix(data = NA, ncol = number_of_genes, nrow = number_of_genes) # initialize the matrix for storing values 
   iteration_values <- matrix(data = NA, ncol=2, nrow = iterations)
   #iteration_values <- as.data.frame(iteration_values)
   naive_iteration_vals[n,1] <- nrow(pathway1)
   for (k in 1:iterations){
     random_pathway <- data.matrix(gene_vecs_subset[sample(nrow(gene_vecs_subset), size = number_of_genes),]) # create a random pathway of length of the first pathway in our list 
     for (i in 1:number_of_genes){ # grab the first vector from that pathway
       vec1 <- random_pathway[i,] - mean(random_pathway[i,]) #the first vector
       for (j in 1:number_of_genes){
         vec2 <- random_pathway[j,] - mean(random_pathway[j,])
         vec.cor <- abs((t(vec1) %*% vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2))))
         pathway_comp_iter[i,j] <- vec.cor
       }
     }
     # store the mean and sd values for this iteration
     iteration_values[k,1] <-  mean(pathway_comp_iter[upper.tri(pathway_comp_iter)])
     iteration_values[k,2] <- sd(pathway_comp_iter[upper.tri(pathway_comp_iter)])
   }

   datalist[[n]] <- iteration_values
 }


naive_pathway_iterations = do.call(rbind, datalist) 
naive_pathway_iterations <- as.data.frame(naive_pathway_iterations)
colnames(naive_pathway_iterations) <- c("correlation_mean", "correlation_sd")

naive_pathway_iterations$pathway <- rep(names(my_pathways), each = iterations) # add a column indicating which pathway the values come from 

# Extract 2.5% and 97.5% quantiles 
iteration_means <- do.call("rbind", tapply(naive_pathway_iterations$correlation_mean, naive_pathway_iterations$pathway, 
                                           quantile, c(0.025, 0.975)))

iteration_sds <- do.call("rbind", tapply(naive_pathway_iterations$correlation_sd, naive_pathway_iterations$pathway, 
                                         quantile, c(0.025, 0.975)))

path_vals <- cbind.data.frame(iteration_means, iteration_sds)
colnames(path_vals) <- c("mean_2.5", "mean_97.5", "sd_2.5", "sd_97.5")

# write output (May 16, 2022)
#write.csv(path_vals, "/Users/amandan/Desktop/Dworkin/dgrp/data/naive_1000_iteration_values.csv")


#### Put the information so far into a plotting dataframe ####

plot_this <- cbind.data.frame(path_vals[,1:2], all_pathways[1:3])

plot_this$x_labels <- paste(plot_this$pathway, " ", "(",plot_this$path_length,")", sep = "")

plot_this$x_labels <- factor(plot_this$x_labels, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", 
                                      "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))


#### Third, calculate the vector correlations for the matched gene expression GO term groups, which I read in here ####

# Because I saved the groups in a very silly and bad file format, use the below code to rescue the groups:
a <- readLines("/Users/amandan/Desktop/Dworkin/dgrp/data/all_gene_groups.csv")
b <- grep("^ *\\[\\d+", unlist(strsplit(a, ' *" *"?')), value = TRUE, invert = TRUE)

all_gene_groups <- lapply(unname(split(b,cumsum(grepl("^\\[\\[\\d+\\]\\]$",b)))), function(x)
  lapply(unname(split(x[-1], cumsum(grepl("^\\[", x[-1])))), `[`,-1))


cor_list <- list()
full_cor_list <- list()


# Do the loop!

Sys.time()
for (h in 1:length(all_gene_groups)){
  group1 <- all_gene_groups[[h]]
  for (i in 1:length(group1)){
    path1 <- group1[[i]]
    pathway_comp <- matrix(data = NA, ncol = length(path1), nrow = length(path1))
    
    for (k in 1:length(path1)){
      gene_i <- which(gene_ids == path1[k])
      vec1 <- gene_vecs[gene_i,] - mean(gene_vecs[gene_i,])
      
      for (j in 1:length(path1)){
        gene_j <- which(gene_ids == path1[j])
        vec2 <- gene_vecs[gene_j,] - mean(gene_vecs[gene_j,])
        #vec.cor <- abs((t(vec1) %*% vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2))))
        vec.cor <- cor(vec1,vec2)
        pathway_comp[k,j] <- vec.cor
      }
    }
    
    
    cor_list[[i]] <- data.frame(mean_correlation = round(mean(pathway_comp[upper.tri(pathway_comp)]), 2),
                                pathway_length = length(path1))
    
  }
  cor_df <- do.call("rbind", cor_list)
  full_cor_list[[h]] <- cor_df
}
Sys.time()


length(full_cor_list)
save_this <- full_cor_list

full_cor_list <- do.call("rbind", full_cor_list)
# write to output
write.csv(x = full_cor_list, "/Users/amandan/Desktop/Dworkin/dgrp/data/matched_exp_correlations_May252022.csv")

full_cor_values <- do.call("rbind", tapply(full_cor_list$mean_correlation, full_cor_list$pathway_length, 
                                         quantile, c(0.025, 0.975)))


#### Fourth, plot the figure again and this time add all of the information together ####

# Order the df correctly first
full_cor_values <- full_cor_values[c(4,1,2,5,8,6,7,3,9),]

plot_this$matched_2.5 <- full_cor_values[,1]
plot_this$matched_97.5 <- full_cor_values[,2]

# Write this dataframe so I can plot it later and not have to run a bunch of stuff
# write.csv(x = plot_this, "/Users/amandan/Desktop/Dworkin/dgrp/data/typical_figure_plot_df.csv")

#### plot pairwise vector correlations ####

plot_this <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/typical_figure_plot_df.csv")
plot_this$x_labels <- factor(plot_this$x_labels, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))

pcols <- c("Matched gene expression correlation group" = "lightblue",
           "Random gene group" = "pink",
           "GO term group mean" = "black" )

ggplot() + 
  geom_crossbar(data = plot_this, aes(x = factor(x_labels), y = 1.5, ymin = mean_2.5, ymax = mean_97.5, colour = "Random gene group", fill = "pink"), fill = "pink", alpha = 0.3) + 
  coord_cartesian(ylim = c(0, 1)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_crossbar(data = plot_this, aes(x = factor(x_labels), y = 1.5, ymin = matched_2.5, ymax = matched_97.5, colour= "Matched gene expression correlation group", fill = "lightblue"),fill = "lightblue", alpha = 0.3) + 
  geom_crossbar(data = plot_this, aes(x = factor(x_labels), y = mean_cor, ymin = mean_cor, ymax = mean_cor, colour = "GO term group mean", fill = "black")) +
  ylab(expression(atop("Pairwise correlation", paste("in direction of effect")))) + 
  xlab("GO term group (number of genes in group)") + 
  scale_colour_manual(name="",values=pcols) + 
  scale_fill_manual(name="",values=pcols) +
  theme(legend.position="bottom", legend.direction="horizontal")

ggsave("pathway_vec_correlations.png", 
       path = "/Users/amandan/Desktop/Dworkin/dgrp/outputs",
       width = 7,
       height = 4)
       
       
       
 #### I can also do this for the vector magnitudes ####
 
 #### Vector Magnitudes ####
# I would like to compare the vector magnitudes... for the true groups and for the 1000 matched expression groups I computed above. So... compute the l2 norm!
# In my mind, how this works is that I would compute the l2 norm of the vector... which would be the estimates from the model... which I do have saved within gene_vecs.. cool.. probably still I should check this with Ian?

# 1. Compute the vector magnitudes for the actual groups
# 2. Compute the vector magnitude for each of the matched gene expression groups

# "True" group vector magnitudes 

# Initialize an empty matrix to put the values in
mag_mat <- matrix(data = NA, nrow = 9, ncol = 3)
mag_mat <- as.data.frame(mag_mat)
colnames(mag_mat) <- c("pathway", "path_length", "vec_mag")

# Do the loop!

for (k in 1:length(my_pathways)){ # For each pathway, create an nxn matrix where n represents each gene in that pathway
  pathway1 <- my_pathways[[k]]
  mag_mat$pathway[k] <- names(my_pathways)[k]
  mag_mat$path_length[k] <- nrow(pathway1)
  mag_comp <- matrix(data = NA, ncol = 1, nrow = nrow(pathway1))
  for (i in 1:nrow(pathway1)){ # For each gene in the pathway, get the vector of expression
    gene_i <- which(gene_ids == pathway1$genes[i]) # index the gene in the expression matrix 
    gene_i_vec <- gene_vecs[gene_i,]
    mag_i <- sqrt(sum((gene_i_vec)^2)) # calculate the magnitude for this gene  
    mag_comp[i] <- mag_i
  }
  
  mag_mat$vec_mag[k] <- round(mean(mag_comp), 5)
}

mag_pathways <- data.frame(mag_mat) # these all seem quite low.. oh well

# Do this for the naive groups

datalist = list() # initiate a datalist 

iterations <- 1000 # set the number of iterations

for (n in 1: length(my_pathways)){
  pathway1 <- my_pathways[[n]]
  number_of_genes <- nrow(pathway1)
  gene_vecs_subset <- gene_vecs[-which(gene_ids %in% pathway1$genes),] # remove the genes that are in one of our known pathways 
  iteration_values <- matrix(data = NA, ncol=2, nrow = iterations)
  naive_mat <- matrix(data = NA, ncol = 1, nrow = number_of_genes)
  for (k in 1:iterations){
    random_pathway <- data.matrix(gene_vecs_subset[sample(nrow(gene_vecs_subset), size = number_of_genes),]) # create a random pathway of length of the first pathway in our list 
    for (i in 1:number_of_genes){ # grab the first vector from that pathway
      naive_mat[i] <- sqrt(sum((random_pathway[i,])^2)) # compute the vector magnitude 
    }
    # store the mean and sd values for this iteration
    iteration_values[k,1] <- mean(naive_mat)
    iteration_values[k,2] <- number_of_genes
  }
  
  datalist[[n]] <- iteration_values
}


mag_naive_iterations = do.call(rbind, datalist) 
mag_naive_iterations <- as.data.frame(mag_naive_iterations)
colnames(mag_naive_iterations) <- c("mean_vec_mag", "pathway_length")


# Extract 2.5% and 97.5% quantiles 
mag_means <- do.call("rbind", tapply(mag_naive_iterations$mean_vec_mag, mag_naive_iterations$pathway_length, 
                                           quantile, c(0.025, 0.975)))




# Finally, do this with the matched correlation in gene expression groups

vec_mag_here <- matrix(data = NA, ncol = 3, nrow = length(all_gene_groups))

for (i in 1:length(all_gene_groups)){
  gene_group_all <- all_gene_groups[[i]]
  gene_group_mean <- matrix(data = NA, nrow = length(gene_group_all))
  for (k in 1:length(gene_group_all)){
    gene_group <- unlist(gene_group_all[k])
    vec_means <- matrix(data = NA, ncol = 1, nrow = length(gene_group))
    for (j in 1:length(gene_group)){
      gene_i <- which(gene_ids == gene_group[j]) # index the gene in the expression matrix 
      gene_i_vec <- gene_vecs[gene_i,]
      mag_i <- sqrt(sum((gene_i_vec)^2)) # calculate the magnitude for this gene  
      vec_means[j] <- mag_i
    }
    gene_group_mean[k] <- mean(vec_means)
  }
  vec_mag_here[i,1] <- length(gene_group)
  vec_mag_here[i,2] <- quantile(gene_group_mean, probs = c(0.025))
  vec_mag_here[i,3] <- quantile(gene_group_mean, probs = c(0.975))
}

#### Plot this Magnitudes ####

# Put it all together and make it into a plot!

plot_mag <- mag_pathways
mag_means <- mag_means[c(4,1,2,5,8,6,7,3,9),]
plot_mag$naive_2.5 <- mag_means[,1]
plot_mag$naive_97.5 <- mag_means[,2]
plot_mag$matched_2.5 <- vec_mag_here[,2]
plot_mag$matched_97.5 <- vec_mag_here[,3]
plot_mag$x_labels <- paste(plot_mag$pathway, " ", "(", plot_mag$path_length, ")", sep = "")

plot_mag$x_labels <- factor(plot_mag$x_labels, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))

#write.csv(plot_mag, "/Users/amandan/Desktop/Dworkin/dgrp/data/plot_mag.csv")

pcols <- c("Matched gene expression correlation group" = "lightblue",
           "Random gene group" = "pink",
           "GO term group mean" = "black" )

ggplot() + 
  geom_crossbar(data = plot_mag, aes(x = x_labels, y = 1, ymin = naive_2.5, ymax = naive_97.5, colour = "Random gene group", fill = "Random gene group"), alpha = 0.3) +
  geom_crossbar(data = plot_mag, aes(x = x_labels, y = 1, ymin = matched_2.5, ymax = matched_97.5, colour = "Matched gene expression correlation group", fill = "Matched gene expression correlation group"), alpha = 0.3) + 
  geom_crossbar(data = plot_mag, aes(x = x_labels, y = vec_mag, ymin = vec_mag, ymax = vec_mag, colour = "GO term group mean", fill = "GO term group mean")) + 
  coord_cartesian(ylim = c(0.001, 0.004)) + 
  theme_classic() + 
  labs(y = "Group mean vector magnitude", x = "Gene group") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(name="",values=pcols) + 
  scale_colour_manual(name="",values=pcols) +
  theme(legend.position="bottom", legend.direction="horizontal")

ggsave("pathway_vec_mags.png", 
       path = "/Users/amandan/Desktop/Dworkin/dgrp/outputs",
       width = 7,
       height = 4)
