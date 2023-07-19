# QUES 1: GIBBS SAMPLER

prec_data <- readRDS("precipitation.rds")

log_data <- log(prec_data)

par(mfrow = c(1,2))
hist(log_data, breaks = 40, prob = TRUE)
hist(prec_data, breaks = 40, prob = TRUE)

#PART A

# setting up the initial values of the priors

n = length(log_data)

# for mu
mu0 = mean(log_data)
tau20 = 1
# for sigma
sig20 = var(log_data)
nu0 = 1

set.seed(1234)

#Gibbs Sampling function
Gibbs_sampler <- function(steps = 500, 
                          mu_0 = mu0, 
                          tau2_0 = tau20, 
                          sig2_0 = sig20, 
                          nu_0 = nu0,
                          n = length(log_data)){
  
  #initial value for gibbs sampler
  mu_i = rnorm(1, mean = mu_0, sd = sqrt(tau2_0))
  sigma2_i =  rinvchisq(1, df = nu_0, scale = sig2_0)
  
  #initialize result DF
  sample_results = data.frame("mu" = mu_i,
                              "sigma2" = sigma2_i)
  
  #initialize loop over specified gibbs steps
  for (i in 1:steps) {
    
    #gibbs step for mu
    w = n/sigma2_i / (n/sigma2_i + 1/tau2_0)
    mu_n = w * mean(log_data) + (1-w) * mu_0
    tau2_n = (n/sigma2_i + 1/tau2_0)^-1
    mu_i = rnorm(1, mean = mu_n, sd = sqrt(tau2_n))
    
    #gibbs step for sigma
    nu_n = nu_0 + n
    scale_term = (nu_0 * sigma2_i + sum((log_data - mu_i)^2))/nu_n
    sigma_i = rinvchisq(1, df = nu_n, scale = scale_term)
    
    #IF
    #IF_mu = 1 + 2*(sum(cor(sample_results[i,1], mu_i)))
    
    #save step values to DF
    sample_results = rbind.data.frame(sample_results, c(mu_i, sigma_i))
    
  }
  
  return(sample_results)

}

joint_posterior <- Gibbs_sampler(steps = 500)

#autocorrelation
par(mfrow = c(2,1))
mu_Gibbs <- acf(joint_posterior[,1])
sig2_Gibbs <- acf(joint_posterior[,2])

#inefficiency factor
IF_Gibbs_mu <- 1+2*sum(mu_Gibbs$acf[-1])
IF_Gibbs_sig2 <- 1+2*sum(sig2_Gibbs$acf[-1])

cat("\n", "IF of mu:", round(IF_Gibbs_mu,3))
cat("\n", "IF of sigma2:", round(IF_Gibbs_sig2,3))

#trace plot
par(mfrow = c(2,1))
plot(1:length(joint_posterior$mu), joint_posterior[,1], type = "l",col="blue", xlab = "steps", ylab = "mu")
plot(1:length(joint_posterior$sigma2), joint_posterior[,2], type = "l",col="blue", xlab = "steps", ylab = "sigma2")

# PART B

set.seed(1234)
post_pred <- c()
post_burnremoved <- joint_posterior[51:length(joint_posterior$mu),]
for (i in 1:(length(post_burnremoved$mu))){
  post_pred[i] = rnorm(1, mean = post_burnremoved[i,1], sd = sqrt(post_burnremoved[i,2]))
}

e_post_draws <- exp(post_pred)


plot(density(prec_data), 
     col = "blue", 
     ylim = c(0,0.2), 
     xlim = c(-10,100),
     main = "Posterior Prediction vs Y")
lines(density(e_post_draws), 
      col = "red")

# QUES 2: METROPOLIS RW

rm(list = ls())

data <- read.table("eBayNumberOfBidderData.dat", header = T)


# PART A

#drop intercept
data_nointercept <- data[,-2]

glm_model <- glm(nBids ~ ., family = poisson, data = data_nointercept)
summary(glm_model)

# PART B
#split response and features
response <- as.matrix(data$nBids)
covariates <- as.matrix(data[,2:10])

#log posterior

logPost <- function(beta, 
                    X = covariates, 
                    Y = response){
  
  loglik <- sum(Y * (X %*% beta) - exp(X %*% beta))
  logprior <- dmvnorm(t(beta), mean = matrix(0, nrow = ncol(X)), sigma = (100*(solve(t(X) %*% X))), log = T)
  return(loglik + logprior)

}

#optimize log posterior

opt_res <- optim(par = matrix(1, nrow = 9),
                 fn = logPost,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

beta_mode <- opt_res$par
inv_jacobian <- -solve(opt_res$hessian)

rownames(beta_mode) <- colnames(covariates)

t(beta_mode)
glm_model$coefficients

# PART C

# Metropolis Random Walk function
RMW_func <- function(target_density,
                     c,
                     theta_i_1,
                     sigma_proposal,
                     steps,
                     X,
                     Y){
  
  set.seed(12345)
  
  # initialize result matrix
  result <- matrix(t(theta_i_1), ncol = 9)
  
  accept = 0
  
  for (i in 1:steps) {
    
    #sample from the proposal distribution
    theta_p <- rmvnorm(1,
                       mean = as.vector(theta_i_1), 
                       sigma = c*sigma_proposal)
    
    #calculate the ratio for the acceptance probability
    ratio = ((target_density(as.vector(theta_p), X, Y)) 
             - (target_density(as.vector(theta_i_1), X, Y)))
    
    #since target & proposal is defined in log form, we exponentiate to revert
    ratio = exp(ratio)
    
    #calculate alpha
    alpha = min(1,ratio)
    
    #draw from uniform
    u <- runif(1)
    
    #run test
    if(u < alpha){
      accept = accept + 1
      theta_i_1 = theta_p
      result <- rbind(result, theta_i_1)
    }
    else {
      result = rbind(result, as.vector(theta_i_1))
    }
    
  }
  
  return(list(result = result, acceptance = accept/steps))
  
}

# choose same params for the prior of target density
# choose data observed at mode for theta(i-1) & sigma for proposal density
sigma_proposal <- inv_jacobian
#theta_init <- beta_mode
theta_init <- matrix(rep(0.5,9), nrow = 9)

test_mrw <- RMW_func(target_density = logPost,
                     c = 0.6,
                     theta_i_1 = theta_init,
                     sigma_proposal = sigma_proposal,
                     steps = 5000,
                     X = covariates,
                     Y = response)

MRW_coeff <- test_mrw$result
MRW_coeff_means <- apply(MRW_coeff, 2, mean)
names(MRW_coeff_means) <- rownames(beta_mode)
MRW_coeff_means

colnames(MRW_coeff) <- rownames(beta_mode)

par(mfrow = c(3,3))

plot(MRW_coeff[,1], type = 'l', ylab = "Constant", xlab = "steps")
plot(MRW_coeff[,2], type = 'l', ylab = "PowerSeller", xlab = "steps")
plot(MRW_coeff[,3], type = 'l', ylab = "VerifyID", xlab = "steps")
plot(MRW_coeff[,4], type = 'l', ylab = "Sealed", xlab = "steps")
plot(MRW_coeff[,5], type = 'l', ylab = "Minblem", xlab = "steps")
plot(MRW_coeff[,6], type = 'l', ylab = "Majblem", xlab = "steps")
plot(MRW_coeff[,7], type = 'l', ylab = "LargNeg", xlab = "steps")
plot(MRW_coeff[,8], type = 'l', ylab = "LogBook", xlab = "steps")
plot(MRW_coeff[,9], type = 'l', ylab = "MinBidShare", xlab = "steps")

# PART D

sample <- as.matrix(c(1, 1, 0, 1, 0, 1, 0, 1.2, 0.8), nrow = 9)
names(sample) <- rownames(beta_mode)
posterior_samples <- MRW_coeff[301:length(MRW_coeff[,1]),]

linear_predictions <- posterior_samples %*% sample
lambda_predictons <- exp(linear_predictions)

sample_predictions <- c()

set.seed(12345)
for (i in 1:length(lambda_predictons)){
  sample_predictions[i] <- rpois(1, lambda_predictons[i])
}

hist(sample_predictions, breaks = 20)

cat("\n", "Probability of no bids on the test data:",length(which(sample_predictions == 0))/length(sample_predictions))


# QUES 3: TIME SERIES

rm(list = ls())

AR <- function(mu = 13, 
               phi = c(-0.95, -0.5, 0, 0.25, 0.75, 0.95), 
               sigma2 = 3, 
               T = 300) {
  
  chains = data.frame(0, ncol = length(phi))
  
  for (i in 1:length(phi)){
    
    x_t <- mu
    chains[1,i] = x_t
    
    for (j in 2:T){
      
      x_t = mu + phi[i] * (x_t - mu) + rnorm(1, 0 , sqrt(sigma2))
      chains[j,i] = x_t
      
    }
    
  }
  
  colnames(chains) <- paste("phi_",phi)
  
  return(chains)
  
}

test_AR1 <- AR()

par(mfrow = c(3,3))
plot(y = test_AR1[,1], x = c(1:300), type = 'l', ylab = "phi = -0.95", xlab = "T")
plot(y = test_AR1[,2], x = c(1:300), type = 'l', ylab = "phi = -0.5", xlab = "T")
plot(y = test_AR1[,3], x = c(1:300), type = 'l', ylab = "phi = 0", xlab = "T")
plot(y = test_AR1[,4], x = c(1:300), type = 'l', ylab = "phi = 0.25", xlab = "T")
plot(y = test_AR1[,5], x = c(1:300), type = 'l', ylab = "phi = 0.75", xlab = "T")
plot(y = test_AR1[,6], x = c(1:300), type = 'l', ylab = "phi = 0.95", xlab = "T")

par(mfrow = c(2,3))
sapply(test_AR1, function(x) acf(x, main=""))

# PART B

model_AR <- AR(phi = c(0.2, 0.95))

StanModel = '
data {
  int<lower=0> T; // Number of observations
  vector[T] x;  // indicating the chain x over T obs
}
parameters {
  real mu;
  real<lower = 0> sigma2;
  real<lower = -1, upper = 1> phi;
}
model {
  mu ~ normal(0,50); // Normal with mean 0, st.dev. 50
  sigma2 ~ scaled_inv_chi_square(1,10); // Scaled-inv-chi2 with nu 1,sigma 10
  phi ~ uniform(-1,1);
  
  for(i in 2:T){
    x[i] ~ normal(mu + phi * (x[i-1] - mu), sqrt(sigma2));
  }
}'

#fit model for phi 0.2
fit_0.2 = stan(model_code = StanModel,
               data = list(x = model_AR$`phi_ 0.2`, T = 300),
               warmup = 1000,
               iter = 2000,
               chains = 4)

#fit model for phi 0.95
fit_0.95 = stan(model_code = StanModel,
                data = list(x = model_AR$`phi_ 0.95`, T = 300),
                warmup = 1000,
                iter = 2000,
                chains = 4)

#PART B.i
#get posterior samples and means
post_samples_0.2 <- extract(fit_0.2)
post_mean_0.2 <- get_posterior_mean(fit_0.2)

post_samples_0.95 <- extract(fit_0.95)
post_mean_0.95 <- get_posterior_mean(fit_0.95)

cat("\n")
print("Posterior Means when phi = 0.2 :")
cat("\n")
print(post_mean_0.2)


cat("\n")
print("Posterior Means when phi = 0.95 :")
cat("\n")
print(post_mean_0.95)



post_samples_0.2_df <- data.frame(mu = post_samples_0.2$mu,
                                  sigma2 = post_samples_0.2$sigma2,
                                  phi = post_samples_0.2$phi)

CI_0.2 <- sapply(post_samples_0.2_df, function(x) quantile(x, probs=c(0.025, 0.975)))

post_samples_0.95_df <- data.frame(mu = post_samples_0.95$mu,
                                   sigma2 = post_samples_0.95$sigma2,
                                   phi = post_samples_0.95$phi)

CI_0.95 <- sapply(post_samples_0.95_df, function(x) quantile(x, probs=c(0.025, 0.975)))

cat("\n")
print("Posterior 95% CI when phi = 0.2 :")
cat("\n")
print(CI_0.2)

cat("\n")
print("Posterior 95% CI when phi = 0.95 :")
cat("\n")
print(CI_0.95)

#PART B.ii.

#plotting the two models
par(mfrow = c(2,1))
plot(y = model_AR[,1], x = c(1:300), type = 'l', ylab = "phi = 0.2", xlab = "T")
plot(y = model_AR[,2], x = c(1:300), type = 'l', ylab = "phi = 0.95", xlab = "T")

#checking convergence of the samplers
par(mfrow = c(3,2))
plot(post_samples_0.2$mu, type = 'l', ylab = "posterior Mu", main = "phi = 0.2")
plot(post_samples_0.95$mu, type = 'l', ylab = "posterior Mu", main = "phi = 0.95")
plot(post_samples_0.2$sigma2, type = 'l', ylab = "posterior Sigma2", main = "phi = 0.2")
plot(post_samples_0.95$sigma2, type = 'l', ylab = "posterior Sigma2", main = "phi = 0.95")
plot(post_samples_0.2$phi, type = 'l', ylab = "posterior phi", main = "phi = 0.2")
plot(post_samples_0.95$phi, type = 'l', ylab = "posterior phi", main = "phi = 0.95")

#joint posterior
par(mfrow = c(1,2))
plot(x = post_samples_0.2_df$mu,
     y = post_samples_0.2_df$phi,
     main = "phi = 0.2",
     xlab = "mu",
     ylab = "phi")
plot(x = post_samples_0.95_df$mu,
     y = post_samples_0.95_df$phi,
     main = "phi = 0.95",
     xlab = "mu",
     ylab = "phi")

