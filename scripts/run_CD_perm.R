##### R script to be used in batch mode to run the Churchill and Doerge permutation test

#Load in the functions for full R2 and partial R2

shapeRsq <- function( model ){
  fitted.variance <- sum(diag(var(model$fitted)))
  total.variance	<- sum(diag(var(model$fitted + model$resid)))
  fitted.variance / total.variance
}

comment(shapeRsq) <- "this function takes a multivariate model object and returns the fitted variance / total variance: equivalent to an Rsquared"


shapePRsq <- function( model ){
  # Based on the derivation from page 269 of Kutner et. al. (Applied linear statistical models edition 5.)
  residual.variance <- var(model$resid)
  variables <- attr(terms(model), "term.labels")
  model.length <- length(variables)
  variable.name <- rep(NA, model.length )
  partial.Rsq <- rep(NA, model.length )
  
  for (i in 1:model.length){
    variable.name[i] <- variables[i]
    drop <- parse( text=variables[i] )
    new.formula <- as.formula( paste( ".~.-", variables[i], sep=""))
    new.model <- update(model, new.formula )
    partial.Rsq[i] <- (sum ( diag( var(new.model$resid))) - sum( diag( residual.variance)) ) / sum( diag( var(new.model$resid)))
  }
  R2 <- shapeRsq( model )
  list(Rsquared=R2, partials=data.frame( cbind( variable.name, partial.Rsq ))	)
}

#### Load in data
covariates <- read.csv(file = "DGRP_covariates_final.csv")
gene_sample <- read.csv(file = "DGRP_genes_final.csv")
shape_coords <- read.csv(file = "DGRP_shape_final.csv")

##### Run the permutation #####

n <- 200 # how many iterations of the whole permutation

#partial_vals_max <- rep(NA, length = n) # to store the maximum values from each iteration.
partial_vals_max <- matrix(data = NA, nrow = n, ncol = 3)
colnames(partial_vals_max) <- c("pRsq", "Rsq", "Vec")

for (j in 1:n){
  
  # shuffle the rows of gene expression
  shuffle_index <- sample(nrow(shape_coords)) # number of strains
  gene_exp <- gene_sample[shuffle_index,] 
  #partial_vals <- rep(NA, length = ncol(gene_exp))
  partial_vals <- matrix(data = NA, nrow = ncol(gene_exp), ncol = 3)
  
  for (i in 1:ncol(gene_exp)){
    
    model_fit <- lm(as.matrix(shape_coords) ~ 1 + as.matrix(covariates$log_csize_cent) + as.matrix(covariates[,2:5]) + gene_exp[,i] + as.matrix(covariates$pop_sub...2.) + as.matrix(covariates$pop_sub...3.) + as.matrix(covariates$wolbachia))
    get_prsq <- shapePRsq(model_fit)
    partial_vals[i,1] <- as.numeric(get_prsq$partials[3,2]) #partial R2
    partial_vals[i,2] <- shapeRsq(model_fit) #full R2
    partial_vals[i,3] <- sqrt(sum((coef(model_fit)[7,])^2)) #magnitude of effect size vector 
    
  } 
  
  partial_vals_max[j,1] <- max(partial_vals[,1])
  partial_vals_max[j,2] <- max(partial_vals[,2])
  partial_vals_max[j,3] <- max(partial_vals[,3])
}

#My awful workaround to the overwriting data issue and ugly output issue in bash is to use 5
#separate R scripts.... one for each iteration of 200 permutations

write.csv(partial_vals_max, "perm_max_vals_2.csv", row.names = FALSE)
