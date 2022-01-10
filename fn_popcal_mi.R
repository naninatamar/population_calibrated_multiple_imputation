#' Function to multiply impute missing variable entries using 
#' population calibrated multiple imutation. 
#'
#' @param data Data containing all variables (including the missing variable) that should be used in the imputation model. The missing variable has to be of class factor.
#' @param missing_covariate The name of the missing variable (which should be imputed) in the data. Must be a string/character. 
#' @param true_distribution A vector giving the true proportions of each level of the missing variable known on aggregate on the population level. The first entry has to correspond to 
#' the proportion for the 1st level of the missing variable, the second entry of the vector to the second level of the missing variable etc. 
#' @param true_distribution_n Only used if the external distribution is estimated ("known with uncertainty"). Corresponds to the sample size true_distribution_n of the sample that was used to estimate the external_proportion_true.  
#' @param m The number of imputed datasets
#' @return imp = a list containing the m imputation of the missing observations. Ordered such that the for each entry in the list, the first value corresponds to the first missing value, the second to the second etc. 
#' all_converged = if TRUE then convergence was met in all 50 imputations. 

fn_popcal_mi <- function(data, missing_covariate, true_distribution, true_distribution_n = NULL, m=50){
imp <- list(NULL)
if(!is.null(true_distribution_n)){
  varcov.true_distribution <- -true_distribution %*% t(true_distribution) * true_distribution_n
  diag(varcov.true_distribution) <- true_distribution * (1-true_distribution) * true_distribution_n
  exp_mu = true_distribution * true_distribution_n
  true_distribution <- mvrnorm(n=1, mu = exp_mu, Sigma = varcov.true_distribution)
  true_distribution <- true_distribution/sum(true_distribution)
}
true_distribution = true_distribution[-1]
missing_cov <- data %>% dplyr::select(missing_covariate) %>% as_vector()
levels_missing <- levels(missing_cov)
n_levels <- length(levels(missing_cov))
which_col <- which(names(data) == missing_covariate)

data_mis <- data[is.na(data[,which_col]), ]
data_nonmis <- data[!is.na(data[, which_col]),]
nmis <- nrow(data_mis)
nnonmis <- nrow(data_nonmis)

data.popcal.impute <- matrix(rep(missing_cov, times=m), ncol=m)
converged <- rep(NA, m)

formula_reg <- as.formula(paste(missing_covariate, "~", "."))

for(i in 1:m){
  # 1st step fit multinomial regression for missing variable based on non-missing observations
  fit.nonmis <- multinom(formula_reg, data=data_nonmis)
  coef.nonmis <- summary(fit.nonmis)$coefficients
  var.nonmis <- summary(fit.nonmis)$standard.errors^2

  non_converged <- TRUE
  for(j in 1:1000){
    while(non_converged){
      # Step 2 random draw of regression coefficients
      coeff.r <- matrix(NA, nrow=nrow(coef.nonmis), ncol=ncol(coef.nonmis))
      for(k in 1:(n_levels-1)){
        coeff.r[k, ] <- mvrnorm(1, mu = coef.nonmis[k,], Sigma = diag(var.nonmis[k,]))
      }
      
      # Step 3: Random draw of probability of observing incomplete variable  
      prop.r <-  nnonmis/(nmis+nnonmis)
      prop.r.draw <- rnorm(1, mean = prop.r, sd = sqrt(prop.r*(1-prop.r)/(nmis+nnonmis)))
      
      # Step 4: probability of observing level=l of incomplete variable 
      design_mat <- model.matrix(fit.nonmis)
      ## linear predictor
      linpred <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat))
      for(k in 1:(n_levels-1)){
      linpred[,k] <- design_mat %*% coeff.r[k,]
      }
      ## probabilities
      prob <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat))
      for(k in 1:(n_levels-1)){
        prob[,k] <- apply(linpred, 1, function(x){exp(x[k])/(1+sum(exp(x)))})
      }
      
      # mean probability
      mean.prob <- 1/nnonmis * apply(prob, 2, sum)
     
      # calculate linear predictor and probabilites for missing data based on covariates:
      design_mat_mis <- model.matrix(~ . , data=data_mis[, -which_col])
      ## linear predictor missing data
      linpred_mis <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat_mis))
      for(k in 1:(n_levels-1)){
        linpred_mis[,k] <- design_mat_mis %*% coeff.r[k,]
      }
      explinpred_mis <- exp(linpred_mis)
      ## probabilities
      prob_mis <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat_mis))
      for(k in 1:(n_levels-1)){
        prob_mis[,k] <- apply(linpred_mis, 1, function(x){exp(x[k])/(1+sum(exp(x)))})
      }
      
      p_missing <- rep(NA, n_levels-1)
      for(k in 1:(n_levels-1)){
        p_missing[k] <- (true_distribution[k] - mean.prob[k]*prop.r.draw)/(1 - prop.r.draw)
      }
      
      # step 5; root solve 
      fun_root <- function(x){
        jj <- 1:nmis
        r <- rep(NA, length(x))
        for(ii in 1:length(x)){
        r[ii] <- p_missing[ii] - 1/nmis * sum((exp(x[ii])*explinpred_mis[jj,ii]) /
                                                      (1+rowSums(exp(x)*explinpred_mis[jj,]))) 
        }
        r
      }
      p0 <- c(0,0,0)
      rootsovle <- BBsolve(par=p0, fn = fun_root)
      non_converged <- rootsovle$message != "Successful convergence"
    }}
  converged[i] <- rootsovle$message
  delta <- rootsovle$par
  
  ## new random draw of regression coefficients
    coeff.r.new <- matrix(NA, nrow=nrow(coef.nonmis), ncol=ncol(coef.nonmis))
  for(k in 1:(n_levels-1)){
    coeff.r.new[k, ] <- mvrnorm(1, mu = coef.nonmis[k,], Sigma = diag(var.nonmis[k,]))
  }
  
  ## linear predictors for missing variables including the delta
  design_mat_mis
  linpred_mis_new <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat_mis))
  for(k in 1:(n_levels-1)){
    linpred_mis_new[,k] <- design_mat_mis %*% coeff.r.new[k,] + delta[k]
  }
  ## corresponding probabilities
  prob_post <- data.frame(matrix(NA, ncol=(n_levels), nrow=nrow(design_mat_mis)))
  for(k in 1:(n_levels)){
    if (k==1){
      prob_post[,k] <- apply(linpred_mis_new, 1, function(x){1/(1+sum(exp(x)))})
    } else{
    prob_post[,k] <- apply(linpred_mis_new, 1, function(x){exp(x[k-1])/(1+sum(exp(x)))})
    }
  }
  
  nc <- ncol(prob_post)
  un <- rep(runif(nmis), each = nc)
  draws <- un > apply(prob_post, 1, cumsum)
  idx <- 1+apply(draws, 2, sum)
  ret.lev <- levels_missing[idx]
  ret.lev <- as.factor(ret.lev)
  imp[[i]] <- ret.lev
}
all_converged = all(converged == "Successful convergence")
  return(list(imp = imp, all_converged = all_converged))
}
  