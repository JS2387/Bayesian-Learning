# Question 1

# PART A

#read the temperature data
temp_lamb <- read.table("TempLambohov.txt", header = TRUE)
plot(temp_lamb$time, temp_lamb$temp)

#set the starting parameters

mu_0 <- c(-10,100,-100)
omega_0 <- 0.02 * diag(3)
nu_0 <- 3
sigma2_0 <- 2

Y <- temp_lamb$temp
X <- cbind(1, temp_lamb$time, temp_lamb$time^2)

#prior draw

prior_fn <- function(covariates = X, 
                     mu = mu_0, 
                     omega = omega_0, 
                     sigma2 = sigma2_0, 
                     nu = nu_0) {
  
  #Inverse chi sq draw for variance
  #prior_sigma2 <- rchisq(1, df = nu) #draw from chi-sq with nu df
  #prior_sigma2 <- nu * sigma2/prior_sigma2 #compute for inv-chi-sq draws
  
  prior_sigma2 <- rinvchisq(n = 1, df = nu, scale = sigma2)
  
  prior_beta <-rmvnorm(1, mean = mu, sigma = prior_sigma2 * solve(omega))
  
  prior_final <- covariates %*% t(prior_beta)
  
  return(prior_final)
  
}

plot(temp_lamb$time, temp_lamb$temp)
for (i in 1:100) {
  lines(x = temp_lamb$time, y = prior_fn(), col = 'red')
}

plot(temp_lamb$time, temp_lamb$temp)
for (i in 1:100) {
  lines(x = temp_lamb$time, 
        y = prior_fn(mu= c(-20,145,-135), 
                     omega = 0.07 * diag(3), 
                     nu = 5, 
                     sigma2 = 0.15), 
        col = 'red')
}

#setting new initial params

mu_set= c(-20,145,-135)
omega_set = 0.07 * diag(3)
nu_set = 5
sigma2_set = 0.15

#PART B.1.

#Joint Posterior Distribution

joint_posterior <- function(covariates = X, 
                            responses = Y, 
                            omega = omega_set, 
                            mu = mu_set, 
                            nu = nu_set, 
                            sigma2 = sigma2_set) {
  
  #Setup new params for Posterior
  
  beta_hat <- solve(t(covariates) %*% covariates) %*% t(covariates) %*% responses
  mu_n <- (solve(t(covariates) %*% covariates + omega)) %*% (t(covariates) %*% covariates %*% beta_hat + omega %*% mu)
  omega_n <- t(covariates) %*% covariates + omega
  nu_n <- nu + length(responses)
  sigma2_n <- 1/nu_n * (nu * sigma2 + t(responses) %*% responses + t(mu) %*% omega %*% mu - t(mu_n) %*% omega_n %*% mu_n)
  
  
  #marginal posteriors of betas and variance
  
  sigma_posterior <- rinvchisq(n = 10000, df = nu_n, scale = as.numeric(sigma2_n))
  
  beta_posterior <- rmvnorm(n = 10000, mean = mu_n, sigma = as.numeric(sigma2_n) * solve(omega_n))
  
  return(cbind.data.frame(sigma_posterior, beta_posterior))
  
}

posterior <- joint_posterior()

hist(posterior[,1], breaks = 100, probability = TRUE, xlab = "sigma2", main = NULL)
hist(posterior[,2], breaks = 100, probability = TRUE, xlab = "beta0", main = NULL)
hist(posterior[,3], breaks = 100, probability = TRUE, xlab = "beta1", main = NULL)
hist(posterior[,4], breaks = 100, probability = TRUE, xlab = "beta2", main = NULL)

#PART B.2.

beta_post <- posterior[,-1]

pred_all <- as.matrix(beta_post) %*% t(X)

pred_median <- apply(pred_all, 2, FUN = median)

temp_lamb_pred <- temp_lamb
temp_lamb_pred$pred_median <- pred_median

pred_pi <- data.frame(lower = double(), upper = double())

for (i in 1:length(pred_median)) {
  pred_pi[i,] <- quantile(pred_all[,i], c(0.025,0.975)) 
}

temp_lamb_pred <- cbind.data.frame(temp_lamb_pred, pred_pi)

gp <- ggplot(data = temp_lamb_pred, aes(x = time)) + 
  geom_point(aes(y = temp), size = 0.5) +
  geom_line(aes(y = pred_median, color = "red")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.05, 
              color = "blue", 
              fill = "blue", 
              linetype = "dashed")

gp  

#PART C

max_preds <- apply(pred_all, 2, max)

temp_lamb_pred$max_pred <- max_preds

ggplot(data = temp_lamb_pred, aes(x = time)) + 
  geom_point(aes(y = temp), size = 0.5) +
  geom_line(aes(y = pred_median, color = "red")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.05, 
              color = "blue", 
              fill = "blue", 
              linetype = "dashed") +
  geom_line(aes(y = max_pred), color = "green")


# Question 2: Logisitic Regression

rm(list = ls())

# PART A

women_work <- read.table("WomenAtWork.dat", header = TRUE)

#setting up data & parameters for prior

#data
X <- as.matrix(women_work[,2:8])
Y <- as.matrix(women_work[,1], NCOL = 1)

#prior
nfeatures <- dim(X)[2]
tau <- 5
mu <- as.matrix(rep(0,nfeatures)) # Prior mean vector
sigma <- (tau^2)*diag(nfeatures) # Prior covariance matrix

#loglikelihood function

logPosterior <- function(betas, Y, X, mu, sigma){
  LL <- sum((X %*% betas) * Y - log(1 + exp(X %*% betas)))
  LPrior <- dmvnorm(betas, mu, sigma, log=TRUE)
  return(LL + LPrior)
}

#initialize beta for optim
beta_init <- matrix(0, nfeatures, ncol = 1)

#maximize logPosterior
OptimRes <- optim(beta_init,
                  logPosterior,
                  gr=NULL,
                  Y,X,mu,sigma,
                  method=c("BFGS"),
                  control=list(fnscale=-1),
                  hessian=TRUE)

#posterior mode
beta_mode <- OptimRes$par
names(beta_mode) <- colnames(X)

#Hessian is -J(beta)
beta_jacobian <- -OptimRes$hessian

#diagonal of hessian is the second derivates corresponding to the covariance
beta_inverse_jacobian <- solve(beta_jacobian)
colnames(beta_inverse_jacobian) <- colnames(X)

cat("\n", "The posterior mode is: ", "\n")
print(beta_mode)

cat("\n", "The Inverse Jacobian is: ", "\n")
print(beta_inverse_jacobian)

#draw large sample for beta corresponding to NSmallChild i.e. beta6
beta_sim <- rmvnorm(n = 1000, mean = beta_mode, sigma = solve(beta_jacobian))

#isolate draws for beta6
beta6_sim <- beta_sim[,6]

#95% PI for beta6
beta6_pi <- quantile(beta6_sim, prob = c(0.025,0.975))
names(beta6_pi) <- c("Lower", "Upper")

cat("\n", "The PI for Beta_nSmallChild:", "\n")
print(beta6_pi)

#comparison with glm output

glmModel <- glm(Work ~ 0 + ., data = women_work, family = binomial)

summary(glmModel)

#Part B

#Create matrix of sample features
X_data <- as.matrix(c(1, 20, 12, 8, 43, 0, 2))

#reuse the betas simulated earlier
#use sigmoid function to model probabilities of predictions

post_pred <- function(X, beta) {
  
  pred <- beta %*% X #linear predictions
  
  pred_logit <- exp(pred)/(1 + exp(pred)) #sigmoid function to model probabilities
  
  return(data.frame(Probability = pred_logit))
}

#function call on sample data
post_pred_sim <- post_pred(X_data, beta_sim)

#plotting the posterior prediction density & histogram
ggplot(data = post_pred_sim, aes(x = Probability)) + xlim(c(0,1)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 100, 
                 color = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.1, fill = "green") + 
  geom_vline(xintercept = 0.5, alpha = 0.8, linetype = "dashed", color = "red")

#PART C

post_pred_binom <- function(X, beta) {
  
  pred <- beta %*% X_data #linear predictions
  
  pred_logit <- 1/(1 + exp(-pred)) #sigmoid function to model probabilities
  
  trials <- c()
  
  for (i in 1:length(pred_logit)) {
    trials[i] <- rbinom(n = 1, size = 11, prob = pred_logit[i])
  }
  
  return(trials)
}

test_trials <- post_pred_binom(X_data, beta_sim)
hist(test_trials)