#' Function to multiply impute missing variable entries for a binomial incomplete variable using 
#' population calibrated multiple imputation. 
#'
#' @param data Data containing all variables (including the incomplete variable) that should be used in the imputation model. The incomplete variable has to be numeric with entries 0 and 1. 
#' @param missing_covariate The name of the missing variable (which should be imputed) in the data. Must be a string/character. 
#' @param true_prop A numeric giving the true proportion of incomplete variable == 1 known on aggregate on the population level.  
#' @param true_prop_n Only used if the external distribution is estimated ("known with uncertainty"). Corresponds to the sample size true_distribution_n of the sample that was used to estimate the true_prop.  
#' @param m The number of imputed datasets
#' @return imp = a list containing the m imputation of the missing observations. Ordered such that the for each entry in the list, the first value corresponds to the first missing value, the second to the second etc. 
#' all_converged = if TRUE then convergence was met in all 50 imputations. 

fn_popcal_mi_bin <- function(data, missing_covariate, true_prop, true_prop_n = NULL, m=50){
imp <- list(NULL)
if(!is.null(true_prop_n)){
  true_prop <- rnorm(1, mean = true_prop, sd = sqrt(true_prop*(1-true_prop)/true_prop_n))
  }
missing_cov <- data %>% dplyr::select(missing_covariate) %>% as_vector()
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
  fit.nonmis <- glm(formula_reg, data=data_nonmis, family = "binomial")
  coef.nonmis <- summary(fit.nonmis)$coefficients[,1]
  var.nonmis <- vcov(fit.nonmis)

  non_converged <- TRUE
  for(j in 1:1000){
    while(non_converged){
      # Step 2 random draw of regression coefficients
      coeff.r <- mvrnorm(n = 1, mu=coef.nonmis, Sigma = var.nonmis)
      
      # Step 3: Random draw of probability of observing incomplete variable  
      prop.r <-  nnonmis/(nmis+nnonmis)
      prop.r.draw <- rnorm(1, mean = prop.r, sd = sqrt(prop.r*(1-prop.r)/(nmis+nnonmis)))
      
      # Step 4: draw probability of observing incomplete variable =1 
      predict.probs <- predict(fit.nonmis, data=data_nonmis, type="response")
      
      # design_mat <- model.matrix(fit.nonmis)
      # ## linear predictor
      # linpred <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat))
      # for(k in 1:(n_levels-1)){
      # linpred[,k] <- design_mat %*% coeff.r[k,]
      # }
      # ## probabilities
      # prob <- matrix(NA, ncol=(n_levels-1), nrow=nrow(design_mat))
      # for(k in 1:(n_levels-1)){
      #   prob[,k] <- apply(linpred, 1, function(x){exp(x[k])/(1+sum(exp(x)))})
      # }
      # 
      # mean probability random draw
      mean.prob <- mean(predict.probs)
      mean.prob.draw <- rnorm(n=1, mean.prob, sd=sqrt(mean.prob*(1-mean.prob)/nnonmis))
     
       # calculate linear predictor and probabilites for missing data based on covariates:
      design_mat_mis <- model.matrix(~ . , data=data_mis[, -which_col])
      ## linear predictor missing data
      linpred_mis <- design_mat_mis %*% coeff.r
      explinpred_mis <- exp(linpred_mis)
      ## probabilities
      prob_mis <- 1/(1+exp(-linpred_mis))
      
      p_missing <- (true_prop - mean.prob.draw*prop.r.draw)/(1-prop.r.draw)
      # step 5; root solve 
      fun_root <- function(x){
        jj <- 1:nmis
        r <- abs(p_missing - 1/nmis * sum(1/(1+exp(-x)*exp(-linpred_mis[jj]))))
        r
      }
      r0 <- 0
      rootsolve <- optim(r0, fn = fun_root, method ="BFGS")
      non_converged <- rootsolve$convergence != 0
    }}
  converged[i] <- rootsolve$convergence
  delta <- rootsolve$par
  
  ## new random draw of regression coefficients
  coeff.r.new <- mvrnorm(n = 1, mu=coef.nonmis, Sigma = var.nonmis)
  
  ## linear predictors for missing variables including the delta
  linpred_mis_new <- design_mat_mis %*% coeff.r.new + delta
  ## corresponding probabilities
  prob_post <- 1/(1+exp(-linpred_mis_new))
  vec <- (runif(nrow(prob_post)) <= prob_post)
  vec[vec] <- 1
  imp[[i]] <- vec
}
all_converged = all(converged == 0)
  return(list(imp = imp, all_converged = all_converged))
}
  