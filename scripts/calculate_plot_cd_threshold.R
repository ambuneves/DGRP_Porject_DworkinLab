#### R script to combine results of the 5 different permutation tests and
#### get the thresholds for each metric!!!

library(tidyr)
library(ggplot2)

### Load in files and combine into one dataframe

# I ran the permutation tests in separate runs to speed things up, 
# here I combine those files into one 
setwd("/Users/amandan/Desktop/Dworkin/dgrp/permutation_tests/results")

my_files <- list.files(pattern = ".csv")

results <- lapply(my_files, function(i){
  x = read.csv(i)
  x
})

results[[1]]
results <- do.call("rbind.data.frame", results)

# calculate the significance threshold, which is the 95th percentile
quantile(results$pRsq, 0.95) #0.08948
quantile(results$Rsq, 0.95) #0.20912
quantile(results$Vec, 0.95) #0.00599

#Visualize the results

# these are the "actual" values 
full_model <- read.csv("/Users/amandan/Desktop/Dworkin/dgrp/data/model_results_nov1.csv")

#### 

results$Vec[results$Vec >= quantile(results$Vec, probs = 0.95)]

hist(full_model$partialR2, col = "red", xlim = c(0,0.1))
hist(results$pRsq, col = "blue", add = TRUE)
legend("topright", legend = c("linear model results", "permutation values"), fill = c("red", "blue"), cex = 0.5)

hist(full_model$vector_magnitude, col = "red", xlim = c(0,0.006))
hist(results$Vec, col = "blue", add = TRUE)
legend("topright", legend = c("linear model results", "permutation values"), fill = c("red", "blue"), cex = 0.7)

plot(density(full_model$partialR2), col = "red", xlim =c(0,0.12) )
lines(density(results$pRsq), col = "blue")

ggplot() + geom_density(aes(x=partialR2), fill = "red", data = full_model, alpha = 0.5) + 
  geom_density(aes(x=pRsq), fill = "blue", data = results, alpha=0.5)

quantile(results$pRsq, probs = c(0.01, 0.05, 0.1, 0.2))
quantile(full_model$partialR2, probs = c(0.9,0.95,1))

# Make plot for thesis 
library(stringr)
cols <- c("All observed values" = "lightblue", "Highest values from each\n of 1000 shuffled permutations" = "grey")

ggplot() + 
  stat_density(geom = "area", aes(x=vector_magnitude, fill = "All observed values"), data = full_model, alpha = 0.5) + 
  stat_density(geom = "area", aes(x=Vec, fill = "Highest values from each\n of 1000 shuffled permutations"), data = results, alpha=0.5) + 
  geom_vline(xintercept = quantile(results$Vec, 0.95), linetype="dashed", 
             color = "black", size=0.3) + 
  theme_classic() +
  xlab("Magnitude of effect vector") +
  scale_fill_manual(name = "",values = cols)
  

ggsave("mag_vec_density.png", 
       path = "/Users/amandan/Desktop/Dworkin/dgrp/outputs",
       width = 8,
       height = 4)

quantile(results$Vec, probs = c(0.01, 0.05, 0.1, 0.2))
quantile(full_model$vector_magnitude, probs = c(0.9,0.95,1))

max(full_model$vector_magnitude)
