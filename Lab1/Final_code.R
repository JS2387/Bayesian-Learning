# Question 1 : Daniel Bernoulli

rm(list = ls())

#Alpha and Beta values for prior
a0 = 5
b0 = 5

#sample data
n = 50
s = 13

#Alpha and beta values for posterior
a1 <- a0 + s
b1 <- b0 + n - s

#PART A

#true mean and sd
true_mean <- a1/(a1+b1)
true_sd <- sqrt(a1*b1/((a1+b1)^2 * (a1+b1+1)))

#initialize vectors for posterior draws, mean and sd
post_mean <- c()
post_sd <- c()
posterior <- c()

#initialize vector for simulating number of draws
nDraws <- seq(10, 10000, by = 10)

#loop to draw random numbers and calculate each mean and sd
for(i in 1:length(nDraws)){
  posterior <- rbeta(nDraws[i], a1, b1)
  post_mean[i] <- mean(posterior)
  post_sd[i] <- sd(posterior)
}

#plot the means
plot(x = nDraws, y = post_mean, type = "l")
abline(h = true_mean, col = "red")

#plot the sd
plot(x = nDraws, y = post_sd, type = "l")
abline(h = true_sd, col = "red")


#PART B

exact_posterior <- pbeta(0.3, a1, b1)

compute_posterior <- rbeta(10000, a1, b1)

compute_est_posterior <- length(compute_posterior[which(compute_posterior < 0.3)])/10000

cat("\n", "Estimated Posterior Probability: ", round(compute_est_posterior,4))

cat("\n", "Exact Posterior Probability: ", round(exact_posterior,4))

#PART C

phi <- log(compute_posterior/(1 - compute_posterior))

hist(phi, probability = TRUE)

lines(density(phi), col="red", ylim = 1.5)

# Question 2: Log-normal Distribution and Gini Coefficient

rm(list = ls())

y_i <- c(33, 24, 48, 32, 55, 74, 23, 76, 17)
mu = 3.5

#PART A

nDraws = 10000
tau_sq <- sum((log(y_i) - mu)^2)/length(y_i) #given expression for tau^2
X <- rchisq(nDraws, length(y_i)) #draw x from chi-sq with n df
sig_sq <- length(y_i) * tau_sq/X #compute for inv-chi-sq draws

hist(sig_sq, probability = TRUE, breaks = 100)
lines(density(sig_sq), col = "red", ylim = 10)

#PART B

#CDF of std normal with mean = 0, sd = 1, data = sigma^2
phi <- pnorm(sqrt(sig_sq/2))

#Gini Coefficient using given expression 
G <- 2*phi - 1

hist(G, breaks = 100, probability = TRUE)
lines(density(G), col = "red")

#PART C

#ETI cutoffs for 2.5% and 97.5% cutoffs
margin <- quantile(G, probs = c(0.025, 0.975))

hist(G, breaks = 100, probability = TRUE)
lines(density(G), col = "red")
abline(v = margin, col = "blue")

#PART D
G_density <- density(G) #convert G into density
G_df <- data.frame(x = G_density$x, y = G_density$y) #create df of density
G_df <- G_df[order(G_df$y, decreasing = TRUE),] #order df to find cutoff
G_df$cum_y <- cumsum(G_df$y) #add column of cum density in ordered df
G_df_hdi <- G_df[which(G_df$cum_y < 0.95*sum(G_df$y)),] #subset density samples outside HPDI
hdi_low <- min(G_df_hdi$x) #lower HPDI cutoff
hdi_max <- max(G_df_hdi$x) #upper HPDI cutoff

hist(G, breaks = 100, probability = TRUE)
lines(density(G), col = "red")
abline(v = margin, col = "blue")
abline(h = G_df_hdi[G_df_hdi$x == hdi_low,2], col = "green")
abline(v = hdi_low, col = "green", lty  = 2)
abline(v = hdi_max, col = "green", lty  = 2)


# Question 3: Bayesian inference for the concentration parameter in the von Mises distribution

rm(list = ls())

y <- c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)
mu = 2.51

kappa <- seq(0.01,10, by = 0.001)

posterior <- function(k, y, mu){
  post_pdf <- exp(k * (sum(cos(y - mu)) - 1)) / (besselI(k, nu = 0)^length(y))
  return(post_pdf)
}

#draws from the posterior
posterior_test <- posterior(kappa, y, mu)

#normalize the posterior
posterior_density <- posterior_test/sum(kappa*posterior_test)

plot(kappa, posterior_density, type = 'l')

#find the kappa value of the posterior mode
k_mode = kappa[which.max(posterior_density)]

#find the posterior mode
posterior_k_mode <- posterior_density[which.max(posterior_density)]

cat("Approximate Posterior Mode of kappa: ", round(posterior_k_mode,6), "\n")
cat("Kappa value for posterior mode:", round(k_mode,3), "\n")

plot(kappa, posterior_density, type = 'l')
abline(v = k_mode)
