#### Script to compare the vectors of shape (PC1, PC2, PC3) to the vectors of the expression of some genes 

#### Amanda Neves June 6, 2022

##########################

# Packages 

library("geomorph")
library("ggplot2")

# Load in data
# shape data
shape <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/BothLabs_Wings_28Oct.csv")
dim(shape)

# gene vector data 
# all of the gene vectors
gene_vectors <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/gene_vectors_final.csv")
gene_ids <- gene_vectors$gene_id # for easier indexing 
gene_vecs <- as.matrix(gene_vectors[,-1])

# Load in the dictionaries to convert from gene name to fbgn and back
gene_dic <- read.delim("/Users/amandan/Desktop/Dworkin/background_effects/data/FlyGeneDictionary.txt",
                       header = TRUE,
                       sep = "\t",
                       dec = ".")


##### Cleaning shape data ####

# Before computing the shape PCs, need to subset the data to just contain the Dworkin wings
dat <- shape[shape$Lab == "Dwo",]

# Multiple Csize by scale
dat$Csize <- dat$Csize*dat$Scale

# Aggregate by Line number and sex 
sex_line <- aggregate(dat[,10:106], by = list(dat$Line, dat$Sex), FUN = "mean")
colnames(sex_line)[1] <- "LineNum"
colnames(sex_line)[2] <- "Sex"

# Then average between the sexes
shape_df <- aggregate(sex_line[,3:99], by = list(sex_line$LineNum), FUN = "mean")

#### regress out the effect of size so that I can compare this to the gene vectors later

### Linear model to regress size out of shape

lm.1 <- lm(as.matrix(shape_df[,2:97]) ~ as.matrix(shape_df$Csize))
shape_resid <- lm.1$residuals

log2cs <- log2(shape_df$Csize)
lm.1 <- lm(as.matrix(shape_df[,2:97]) ~ as.matrix(log2cs))

shape_resid <- lm.1$residuals

shape_pc <- prcomp(t(shape_resid))
pc_vals <- shape_pc$x
pc_vals[,1]
pc_vals[,2]
pc_vals[,3]

# Read in the pathway data
all_pathways <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/pathways/all_pathways_processed.csv")

my_pathways <- split(all_pathways, f = all_pathways$pathway) # split the pathways into a list by pathway name

pathlist <- list()

for (j in 1:length(my_pathways)){
  datalist <- list()
  path_1 <- as.data.frame(my_pathways[j])
  for (i in 1:nrow(path_1)){
    pc_cors <- list()
    gene_i <- path_1[i,1]
    vec_i <- gene_vecs[which(gene_ids == gene_i),]
    vec_i <- vec_i - mean(vec_i)
    # Loop to compare to the first 5 PCs
    for (n in 1:n_PC){
      pc_cors[[n]] <- c(gene_i, path_1[i,2], cor(vec_i, pc_vals[,n]), n)
    }
    
    pc_cors <- do.call(rbind, pc_cors)
    datalist[[i]] <- pc_cors
  }
  datalist <- do.call(rbind, datalist)
  pathlist[[j]] <- datalist
}

# Make it into a dataframe
path_pc <- as.data.frame(do.call(rbind, pathlist))
colnames(path_pc) <- c("gene_id", "pathway", "value", "PC")
path_pc$value <- as.numeric(path_pc$value)
path_pc$abs <- abs(path_pc$value) # a column for the absolute value cuz why not


#### naive distribution for all the go term gene groups plus facet wrap plot

iterations = 1000 # set the number of iterations

# okay cool!

all_paths <- list()

for (k in 1:length(my_pathways)){
  gene_vecs_sub <- gene_vecs[which(gene_ids %in% my_pathways[[k]]$genes == FALSE),]
  naive_vecs <- list()
  
  for (i in 1:iterations){
    datalist <- list()
    random_path <- gene_vecs_sub[sample(x = 1:nrow(gene_vecs_sub), size = nrow(my_pathways[[k]]), replace = FALSE),]
    for (j in 1:nrow(random_path)){
      pc_cors <- list()
      vec_j <- random_path[j,]
      vec_j <- vec_j - mean(vec_j)
      # Loop to compare to the first 5 PCs
      for (n in 1:n_PC){
        pc_cors[[n]] <- c(abs(cor(vec_j, pc_vals[,n])), n, i, j)
      }
      
      pc_cors <- do.call(rbind, pc_cors)
      datalist[[j]] <- pc_cors
    }
    
    datalist <- do.call(rbind, datalist)
    naive_vecs[[i]] <- datalist
    
  }
  
  temp_df <- as.data.frame(do.call(rbind, naive_vecs))
  colnames(temp_df) <- c("value", "PC", "iteration", "gene")
  temp_df$value <- as.numeric(temp_df$value)
  
  temp_df <- do.call("rbind", tapply(temp_df$value, temp_df$PC, 
                                 quantile, c(0.025, 0.975)))
  
  temp_df <- as.data.frame(temp_df)
  
  temp_df$PC <- 1:5
  
  temp_df$pathway <- names(my_pathways)[k]
  
  colnames(temp_df) <- c("low", "high", "PC", "pathway")
  
  all_paths[[k]] <- temp_df
  
}


# Then mush that dataframe up together....

c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)")

all_paths_df <- as.data.frame(do.call(rbind, all_paths))

path_length <- rep(c("(57)", "(11)", "(13)", "(66)", "(89)", "(75)", "(76)", "(44)", "(107)"), each = 5)

all_paths_df$path_length <- paste(all_paths_df$pathway, path_length, sep = " ")

all_paths_df$path_length <- factor(all_paths_df$path_length, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))

pcols <- c("Matched gene expression correlation group" = "lightblue",
           "Random gene group" = "pink",
           "GO term group mean" = "black" )


# Need to add column to path_pc which has the correct facet wrap title names 

facet_names <- list(
  'bmp'="bmp (57)",
  'cell_elong'="cell_elong (11)",
  'disc_size'="disc_size (13)",
  'egfr'="egfr (66)",
  'hh'="hh (89)",
  'hippo'="hippo (75)",
  'wnt'="wnt (107)",
  'insulin'="insulin (76)",
  'wing_vein'="wing_vein (44)"
)

path_pc$plot_label <- "hello"

for (i in 1:length(path_pc$pathway)){
  this_word <- as.character(path_pc$pathway[i])
  path_pc$plot_label[i] <- as.character(facet_names[this_word])
}


path_pc$plot_label <- factor(path_pc$plot_label, levels = c("cell_elong (11)", "disc_size (13)", "wing_vein (44)", "bmp (57)", "egfr (66)", "hippo (75)", "insulin (76)", "hh (89)", "wnt (107)"))

# Plot for thesis 

ggplot(data = path_pc, aes(x = PC, y = abs)) + 
  geom_boxplot() + 
  geom_crossbar(data = all_paths_df, aes(x = factor(PC), y = 1.5, ymin = low, ymax = high, colour = "Random gene group", alpha = 0.3),fill = "pink", alpha = 0.3) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0, 1)) + 
  theme_classic() + 
  scale_color_manual(name = "", values = c("Random gene group" = "pink")) + 
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  facet_wrap(~plot_label) +
  ylab(expression(atop("Absolute value of correlation", paste("with each principle component")))) + 
  theme(legend.position="bottom", legend.direction="horizontal")

ggsave("pathway_pc_direction_comp_log2.png", 
       path = "/Users/amandan/Desktop/Dworkin/dgrp/outputs",
       width = 7,
       height = 7)
