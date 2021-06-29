# Input data -------------------------------------------------------------------
make_strategies <- function() {
  x <- data.frame(
    strategy_id = 1:5
  )
  
  # Create indicator variables for each treatment strategy
  for (i in 2:nrow(x)) {
    x[, paste0("strategy", i)] <- ifelse(x$strategy_id == i, 1, 0)
  }
  x
}

make_patients <- function() {
  n <- 10
  data.frame(
    patient_id = 1:n,
    age = rnorm(n, mean = 65, sd = 10),
    female = rbinom(n, 1, .5)
  )
} 

make_input_data <- function() {
  # Expand so there is one row for each patient and strategy
  # Like the expand() function in hesim
  patients <- make_patients()
  strategies <- make_strategies()
  x <- merge(patients, strategies)
  
  # Return and set order
  id_cols <- c("strategy_id", "patient_id")
  nonid_cols <- colnames(x)[!colnames(x) %in% id_cols]
  x <- x[, c(id_cols, nonid_cols)]
  attr(x, "n_strategies") <- nrow(strategies)
  x
}

# Parameters -------------------------------------------------------------------
# Set underlying parameter values
set_sick_params <- function(gamma_mean, gamma_sd, 
                            delta_mean, delta_sd) {

  if (is.null(names(gamma_mean))) {
    names(gamma_mean) <- paste0("strategy", 2:(length(gamma_mean) + 1))
  }
  beta_mean <- c(gamma_mean, delta_mean)
  beta_sd <- c(gamma_sd, delta_sd)
  list(mean = beta_mean, sd = beta_sd)
}

set_params <- function(input_data) {
  n_strategies <- attr(input_data, "n_strategies")
  
  # (1) Transitions from sick state 
  # A multinomial logistic regression
  
  ## (1a) Sick -> Sicker
  sick_sicker <- set_sick_params(
    gamma_mean <- runif(n_strategies - 1, min = log(.6), max = log(1)),
    gamma_sd = rep(.05, n_strategies - 1),
    delta_mean <- c(intercept = log(.2), age = log(1.001), female = log(.9)),
    delta_sd <- c(intercept = .07, age = .001, female = .01)
  )
  
  # (1b) Sick -> Death
  sick_death <- set_sick_params(
    gamma_mean <- runif(n_strategies - 1, min = log(.8), max = log(1)),
    gamma_sd = rep(.05, n_strategies - 1),
    delta_mean <- c(intercept = log(.1), age = log(1.003), female = log(.9)),
    delta_sd <- c(intercept = .07, age = .002, female = .01)
  )
  
  # (2) Transitions from sicker state 
  # Modeling of the log(rate) parameter of an exponential distribution
  
  # (2a) Sicker -> Death
  sicker_death <- list(
    mean = c(intercept = log(.3), age = log(1.002), female = log(.9)),
    sd = c(intercept = .05, age = .001, female = .01)
  )
  
  # Return all
  list(
    sick_sicker = sick_sicker,
    sick_death = sick_death,
    sicker_death = sicker_death
  )
}


# Sample parameters
sample_coefs <- function(n = 100, params) {
  mean <- params$mean 
  sd <- params$sd
  k <- length(mean) # Number of total covariates
  sim <- rnorm(n * k, mean = mean, sd = sd)
  out <- t(matrix(sim, nrow = k))
  colnames(out) <- names(mean)
  out
}

sample_params <- function(n = 100, params) {
  lapply(params, function(z) sample_coefs(n = n, params = z))
}

# Simulate state probabilities (i.e., Markov trace) ----------------------------
make_x <- function(data, params) {
  data[, "intercept"] <- 1 # Add intercept
  as.matrix(data[, colnames(params)])
}

mlogit_probs <- function(x, beta) {
  prob <- rep(NA, length(beta) + 1)
  
  # Reference category
  exb <- sapply(beta, function(beta) exp(x %*% t(beta)))
  exb_sum <- sum(exb)
  prob[1] <- 1/(1 + exb_sum)
  
  # Remaining categories
  for (i in 1:length(beta)) {
    prob[i + 1] <- prob[1] * exp(x %*% t(beta[[i]]))
  }
  prob
}

tp_sick <- function(input_data, params) {
  beta <- params[c("sick_sicker", "sick_death")]
  x <- make_x(input_data, beta[[1]])
  mlogit_probs(x, beta)
}

tp_sicker <- function(input_data, params) {
  beta <- params[["sicker_death"]]
  x <- make_x(input_data, beta)
  rate <- exp(x %*% t(beta))
  prob_death <- 1 - exp(-rate)
  c(sick = 0, sicker = 1 - prob_death, death = prob_death)
}

tpmatrix <- function(input_data, params) {
  rbind(
    tp_sick(input_data, params),
    tp_sicker(input_data, params),
    c(0, 0, 1) # Death is an absorbing state
  )
}

sim_markov_chain <- function (x0, p, n_cycles) {
  x <- matrix(NA, ncol = length(x0), nrow = n_cycles)
  x <- rbind(x0, x)
  colnames(x) <- colnames(p)
  rownames(x) <- 0:n_cycles
  for (t in 1:n_cycles) {
      x[t + 1, ] <- x[t, ] %*% p
  }
  data.frame(cycle = 0:n_cycles, x)
}

sim_stateprobs1 <- function(input_data, params,
                            n_cycles = 20,
                            x0 = c(1, 0, 0)) {
  params_sample <- sample_params(n = 1, params) 
  tpmat <- tpmatrix(input_data, params_sample)
  sim_markov_chain(x0, p = tpmat, n_cycles = n_cycles)
}

sim_stateprobs <- function(input_data, params,
                           n_cycles = 20, n_samples = 100,
                           x0 = c(1, 0, 0)) {

  out <- vector(mode = "list", length = nrow(input_data) * n_samples)
  it <- 1 # Counter for loop
  for (s in 1:n_samples) { # Number of parameter samples for PSA
    for (i in 1:nrow(input_data)) { # Loop over strategies & patients
      
      # Simulate state probabilities
      out[[it]] <- sim_stateprobs1(input_data[i, ], params = params, 
                                   n_cycles = n_cycles, x0 = x0)
      
      # Store the ID variables
      out[[it]]$sample <- s
      out[[it]]$strategy_id <- input_data[i, "strategy_id"]
      out[[it]]$patient_id <- input_data[i, "patient_id"]
      
      # Iterate counter
      it <- it + 1
      
    } # End strategy and patient loop
  } # End parameter loop
  out <- do.call("rbind", out)
  out[, c("sample", "strategy_id", "patient_id",
          "cycle", "sick", "sicker", "death")]
}