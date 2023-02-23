# Read in the correlation values
gene_exp_cor <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/pathway_gene-exp_correlation_values.csv", row.names = 1)

#### Get matched correlation values and pathways ####
# Make sure no gene appears in a pathway twice
# Grab each pathway

iterations = 1000
threshold = 0
all_cor_vals <- matrix(data = NA, nrow = nrow(gene_exp_cor), ncol = 4)
all_cor_vals <- data.frame(all_cor_vals)
colnames(all_cor_vals) <- c("pathway_length", "mean_correlation", "97.5", "2.5")
gene_list <- list()
all_gene_groups <- list()

for (j in 1:nrow(gene_exp_cor)){
  target_cor <- gene_exp_cor$mean_correlation[j]
  sim_gene_cor_2 <- matrix(data = NA, nrow = iterations)
  all_cor_vals$pathway_length[j] <- gene_exp_cor$pathway_length[j]
  
  for ( i in 1:iterations){
    a = 1
    comp_pathway <- matrix(data = NA, nrow = nrow(my_pathways[[j]]))
    comp_pathway <- as.data.frame(comp_pathway)
    gene_n_id <- sample(x = 1:ncol(gene_cor_mat_short), size = 1) # start the "seed"
    gene_n_id <- colnames(gene_cor_mat_short)[gene_n_id]
    while( a < nrow(my_pathways[[j]]) + 1){
      comp_pathway[a,1] <- gene_n_id
      gene_n_id <- which(colnames(gene_cor_mat_short) == gene_n_id)
      gene_1 <- gene_cor_mat_short[gene_n_id,]
      if (length(which(abs(gene_1) >= (target_cor - threshold))) < 5) {
        gene_n_id <- sample(names(potential_genes), size = 1)
        next
      }
      potential_genes <- which(abs(gene_1) >= (target_cor - threshold))
      gene_next <- sample(names(potential_genes), size = 1)
      if (gene_next == gene_n_id | gene_next %in% comp_pathway[,1]) {
        gene_n_id <- sample(names(potential_genes), size = 1)
        next # ensure that we do not have the same gene in a pathway twice
      }
      gene_n_id <- gene_next
      a = a + 1
    }
    
    test <- gene_exp_short[,colnames(gene_exp_short) %in% comp_pathway$V1]
    
    test_cor <- cor(test)
    
    sim_gene_cor_2[i,1] <- mean(abs(test_cor[upper.tri(test_cor)]))
    
    gene_list[[i]] <- comp_pathway$V1
    
  }
  
  all_cor_vals$mean_correlation[j] <- mean(sim_gene_cor_2)
  all_cor_vals$`97.5`[j] <- quantile(sim_gene_cor_2, probs = 0.975)
  all_cor_vals$`2.5`[j] <- quantile(sim_gene_cor_2, probs = 0.025)
  all_gene_groups[[j]] <- gene_list
  
}


# Here is where I made a mistake, I saved the output of this as a .csv
# DON'T EVER DO THIS... save as an rdata file in the future
write.table(all_gene_groups, "/Users/amandan/Desktop/Dworkin/dgrp/data/all_gene_groups.csv")

# Here is the code I used to revcover the information from this file
a <- readLines("/Users/amandan/Desktop/Dworkin/dgrp/data/all_gene_groups.csv")
b <- grep("^ *\\[\\d+", unlist(strsplit(a, ' *" *"?')), value = TRUE, invert = TRUE)

all_gene_groups <- lapply(unname(split(b,cumsum(grepl("^\\[\\[\\d+\\]\\]$",b)))), function(x)
  lapply(unname(split(x[-1], cumsum(grepl("^\\[", x[-1])))), `[`,-1))
