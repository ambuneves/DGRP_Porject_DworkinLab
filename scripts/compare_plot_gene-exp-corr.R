#### Script to look at the correlation of gene expression of the genes in different pathways 
# in our dataset, this will be handy because for the "full" permutations
# I am grabbing the genes with the highest correlations, but maybe it would make 
# more sense to get the ones with a degree of correlation that's similar
# to what we see for the "actual" pathways?

# Libaries:

library(ggplot2)

# Read in the pathway data
all_pathways <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/pathways/all_pathways_processed.csv")

my_pathways <- split(all_pathways, f = all_pathways$pathway)

# Read in the gene expression data 

genes <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/dgrp_gene_abd_means.csv")

# Initiate an empty dataframe to put the final values in 
cor_values <- matrix(data = NA, nrow = 9, ncol = 4)
cor_values <- as.data.frame(cor_values)
colnames(cor_values) <- c("pathway_length", "pathway_name", "mean_correlation", "pathway_sd")

for (k in 1:length(my_pathways)){
  path1 <- my_pathways[[k]]$genes # get just the genes 
  path_name <- my_pathways[[k]]$pathway # get the pathway name 
  path_name <- path_name[1]
  path_genes <- genes[,colnames(genes) %in% path1] # get the expression values for the genes in the pathway
  path_cor <- cor(path_genes) # make the correlation matrix 
  unique_cors <- path_cor[upper.tri(path_cor, diag = F)] #exclude the diagonal and the repeated genes
  path_mean <- mean(abs(unique_cors))
  path_sd <- sd(unique_cors)
  
  cor_values[k,1] <- length(path_genes)
  cor_values[k,2] <- path_name
  cor_values[k,3] <- path_mean
  cor_values[k,4] <- path_sd 
}

cor_values$type <- c("blue", "red", "red", "blue", "blue", "blue", "blue", "red", "blue")

#Write this to output

#write.csv(cor_values, "/Users/amandan/Desktop/Dworkin/dgrp/data/pathway_gene-exp_correlation_values.csv")


cor_values <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/pathway_gene-exp_correlation_values.csv")
 

# Initiate an empty dataframe to put the final values in 

path_cor_list <- list()

for (k in 1:length(my_pathways)){
  path1 <- my_pathways[[k]]$genes # get just the genes 
  path_genes <- genes[,colnames(genes) %in% path1] # get the expression values for the genes in the pathway
  path_cor <- cor(path_genes) # make the correlation matrix 
  unique_cors <- path_cor[upper.tri(path_cor, diag = F)] #exclude the diagonal and the repeated genes
  
  path_cor_list[[k]] <- unique_cors

}

names(path_cor_list) <- names(my_pathways)

group_values <- c(abs(path_cor_list[["bmp"]]), abs(path_cor_list[["cell_elong"]]), abs(path_cor_list[["disc_size"]]), abs(path_cor_list[["egfr"]]), abs(path_cor_list[["hh"]]), abs(path_cor_list[["hippo"]]), abs(path_cor_list[["insulin"]]), abs(path_cor_list[["wing_vein"]]), abs(path_cor_list[["wnt"]]))

group_val_names <- c(rep("bmp (57)", times = 1596), rep("cell_elong (11)", times = 55), rep("disc_size (13)", times = 78), rep("egfr (66)", times = 2145), rep("hh (89)", times = 3916), rep("hippo (75)", times = 2775), rep("insulin (76)", times = 2850), rep("wing_vein (44)", times = 946), rep("wnt (107)", times = 5671))

spread_cor_df <- data.frame(group_values = group_values, group_name = group_val_names)


#write.csv(spread_cor_df, "/Users/amandan/Desktop/Dworkin/dgrp/data/plotting_gene-exp_boxplot_df.csv")

spread_cor_df <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/plotting_gene-exp_boxplot_df.csv")
spread_cor_df$group_name <- factor(spread_cor_df$group_name, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))

#### Plot expression box plots ####

boxplot1 <- ggplot() + 
geom_boxplot(data = spread_cor_df, aes(x = group_name, y = group_values)) + 
theme_classic() + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
ylab("Pairwise gene expression correlation") + 
xlab("GO term group (number of genes in group)")

# save to output
ggsave("gene_exp_cor_boxplot.png", 
       path = "/Users/amandan/Desktop/Dworkin/dgrp/outputs",
       width = 7,
       height = 4)

# same thing but violin plot instead of boxplot
violin1 <- ggplot() + 
geom_violin(data = spread_cor_df, 
aes(x = factor(group_name, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)")), y = group_values)) + 
theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
ylab(expression(atop("Pairwise correlation in", paste("gene expression")))) + 
xlab("GO term group (number of genes in group)")

# and then overlay the cross bars to the violin plot to compare the distribution of gene expression correltions
plot_df_cor <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/plot_df_cor.csv")

# Add a column that includes the pathway names 
plot_df_cor <- plot_df_cor[order(plot_df_cor$pathway_length),]
plot_df_cor$path_length_name <- c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)")
plot_df_cor$path_length_names <- factor(plot_df_cor$path_length_name, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))

pcols <- c("Matched gene expression correlation group" = "lightblue",
           "Random gene group" = "pink",
           "GO term group mean" = "black" )

# Plot for thesis

violin1 + geom_crossbar(data = plot_df_cor, aes(x = factor(path_length_names), y = 1.5, ymin = X2.5, ymax = X97.5, colour = "Matched gene expression correlation group", fill = "lightblue"), fill = "lightblue", alpha = 0.3) + 
  coord_cartesian(ylim = c(0, 1)) + 
  theme_classic() +
  geom_crossbar(data = plot_df_cor, aes(x = factor(path_length_names), y = 1.5, ymin = naive_2.5, ymax = naive_97.5, colour = "Random gene group", fill = "pink"), fill = "pink", alpha = 0.3) + coord_cartesian(ylim = c(0,1)) + geom_crossbar(data = plot_df_cor, aes(x = path_length_names, y = group_correlation, ymin = group_correlation, ymax = group_correlation, colour = "GO term group mean", fill = "black")) + scale_colour_manual(name="",values=pcols) + scale_fill_manual(name="",values=pcols) + theme(legend.position="bottom", legend.direction="horizontal") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# save to output
ggsave("gene_exp_cor_violin_overlay.png", 
       path = "/Users/amandan/Desktop/Dworkin/dgrp/outputs",
       width = 7,
       height = 4)
