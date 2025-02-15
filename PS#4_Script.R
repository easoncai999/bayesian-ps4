### PART 2 ###
set.seed(107)

log_posterior <- function(F) {
  if (F <= 0 || F >= 182) {
    return(-Inf)
  }
  return(dnorm(F, mean = 34, sd = 20, log = TRUE))
}

n <- 10000
F_current <- 100  
chain <- numeric(n)
chain[1] <- F_current
current_logp <- log_posterior(F_current)
proposal_sd <- 50
accepts <- 0
for (i in 2:n) {
  F_proposed <- rnorm(1, mean = F_current, sd = proposal_sd)
  proposed_logp <- log_posterior(F_proposed)
  log_acceptance <- proposed_logp - current_logp
  
  if (log(runif(1)) < log_acceptance) {
    F_current <- F_proposed
    current_logp <- proposed_logp
    accepts <- accepts + 1
  }
  chain[i] <- F_current
}

accept_rate <- (accepts / n) * 100
print(paste("Acceptance Rate of the Metropolis Hastings Algorithm:", 
            round(accept_rate, 2), "%"))

plot(chain, type = "l", col = "#0D2D40", xlab = "Iteration",
     ylab = "Fuel Level (L)", main = "Metropolis Hastings Chain for Fuel Level",
     ylim = c(0, 120))
abline(h = 34, col = "#9AD5F8", lty = 5)
legend("topright", legend = c("Chain", "Sensor Reading at 34 L"),
       col = c("#0D2D40", "#9AD5F8"), lty = c(4, 4), bty = "n")

### PART 3 ###
KDE <- function(iteration, MCMC_sample){
  modes <- numeric(length(iteration))
  
  for (i in seq_along(iteration)) {
    index_i <- iteration[i]
    dens <- density(MCMC_sample[1:index_i])
    modes[i] <- dens$x[which.max(dens$y)]
  }
  return(modes)
}

iteration <- seq(50, length(chain), by = 50)
modes <- KDE(iteration, chain)

plot(iteration, modes, type = "l", col = "#0D2D40", lwd = 2,
     xlab = "Iteration", ylab = "Estimated Mode (L)",
     main = "Estimated Mode Convergence using Kernel Density Estimation")
abline(h = 34, col = "#9AD5F8", lty = 2, lwd = 2)
legend("topright", legend = c("Estimated Mode", "True Mode at 34 L"),
       col = c("#0D2D40", "#9AD5F8"), lty = c(1,2), lwd = c(2,2), bty = "n")

### PART 4 ###
mu <- 60
std <- 10
lower_bound <- 0
upper_bound <- 182
pc_sample <- numeric(0)
while(length(pc_sample) < n) {
  temp <- rnorm(10000, mean = mu, sd = std)
  valid <- temp[temp > lower_bound & temp < upper_bound]
  pc_sample <- c(pc_sample, valid)
}

pc_sample <- pc_sample[1:n] 
pc_modes <- KDE(iteration, pc_sample)

plot(iteration, pc_modes, type = "l", col = "#0D2D40", lwd = 2,
     xlab = "Iteration", ylab = "Estimated Mode (L)",
     main = "Positive Control on KDE Mode Estimation")
abline(h = mu, col = "#9AD5F8", lty = 2, lwd = 2)
legend("topright", legend = c("KDE Estimated Mode", "True Mode at 60 L"),
       col = c("#0D2D40", "#9AD5F8"), lty = c(1, 2), lwd = c(2, 2), bty = "n")

final_estimate <- pc_modes[length(pc_modes)]
error <- abs(final_estimate - mu)
print(paste("Final KDE mode estimate:", round(final_estimate, 2)))
print(paste("Absolute error:", round(error, 2)))

### PART 5 ###
library(ggplot2)
set.seed(107)
true_mode <- 60

prior_sample <- runif(n, min = 0, max = 182)
w <- dnorm(prior_sample, mean = true_mode, sd = 20)
w <- w / sum(w)
bmc_samples <- sample(prior_sample, size = n, replace = TRUE, prob = w)
bmc_modes <- KDE(iteration, bmc_samples)

epsilon <- 0.05
convergence <- function(modes, iteration, true_value, epsilon, points = 10) {
  conv_run <- NA
  for (j in 1:(length(modes) - points + 1)) {
    if (all(abs(modes[j:(j + points - 1)] - true_value) < epsilon)) {
      conv_run <- iteration[j]
      break
    }
    else
      conv_run <- "Did not converge"
  }
  return(conv_run)
}

conv_mh  <- convergence(pc_modes, iteration, true_mode, epsilon)
conv_bmc <- convergence(bmc_modes, iteration, true_mode, epsilon)

comparison_table <- data.frame(
  Method = c("Metropolis Hastings", "Bayesian Monte Carlo"),
  Convergence_Runs = c(conv_mh, conv_bmc)
)
kable(comparison_table)

df_compare <- data.frame(
  iter = rep(iteration, 2),
  mode_est = c(pc_modes, bmc_modes),
  method = rep(c("MH", "BMC"), each = length(iteration))
)

ggplot(df_compare, aes(x = iter, y = mode_est, color = method)) +
  geom_line(size = 1) +
  geom_hline(yintercept = true_mode, linetype = "dashed", color = "#9AD5F8", size = 2) +
  labs(title = "Estimated Mode Convergence Comparison: MH vs. BMC",
       x = "Number of Runs",
       y = "Estimated Mode (L)",
       color = "MCMC Method") +
  scale_color_manual(values = c("BMC" = "#5E74DD", "MH" = "#0D2D40")) +
  theme_minimal()