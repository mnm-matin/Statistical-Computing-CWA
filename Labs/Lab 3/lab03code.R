#' ---
#' title: Statistical Computing, Lab 3, Supplementary code
#' author: Finn Lindgren
#' date: 28/1/2020
#' output: pdf_document
#' ---

#' See the Prediction and Proper Scoring Rules longform notes
#' for examples of how to use these functions.

# To access these functions inside other R files, use
#   source("lab3code.R")

#' Functions for linear expectation and log-variance models
#' ----

#' model_Z() constructs the Z-matrices for the example model from Lecture 3.
#' To use different Z-definitions, define a new model_Z function before
#' calling other functions that use model_Z(), such as model_predict()

# Input:
#   x: vector of covariate values
# Output:
#   A list with two named elements:
#     ZE: A matrix with rows of the form (1, x[i], 0, 0)
#     ZV: A matrix with rows of the form (0, 0, 1, x[i])
model_Z <- function(x) {
  Z0 <- model.matrix(~ 1 + x)
  list(ZE = cbind(Z0, Z0 * 0), ZV = cbind(Z0 * 0, Z0))
}

#' neg_log_lik() computes the negative log-likelihood for the general model class from Lecture 4.
# Input:
#   theta: parameter vector of length ncol(Z$ZE)
#   Z: a list of the type returned by model_Z()
#   y: data vector
# Output:
#   A negative-log-likelihood value
neg_log_lik <- function(theta, Z, y) {
  -sum(dnorm(y, mean = Z$ZE %*% theta, sd = exp(Z$ZV %*% theta)^0.5, log = TRUE))
}

#' model_predict() computes either point predictions and model component confidence
#' intervals, or data predictions.

# Input:
#   theta: parameter vector
#   data: covariate information, in whatever format the model_Z() wants it.
#   Sigma_theta: The parameter estimation error covariance matrix
#                If NULL, theta is treated as being a fixed, known parameter vector
#   type: What to compute?
#         'expectation': Point estimates, std.dev., confidence intervals for the expectation model
#         'log-variance': Point estimates, std.dev., confidence intervals for the log-variance model
#         'observation': Point predictions, std.dev., and prediction intervals
#   alpha: The target error probability for the intervals
#   df: Degrees of freedom for t-quantiles for the interval construction. Inf=Normal distribution.
#   nonlinear.correction: TRUE/FALSE. If FALSE, turns off the nonlinear correction term sigma^2/2
#                         in the exponential component of the prediction variance
# Output:
#   A data.frame with columns 'mu', 'sigma', 'lwr', 'upr'
model_predict <- function(theta, data, Sigma_theta = NULL,
                          type = c("expectation", "log-variance", "observation"),
                          alpha = 0.05, df = Inf,
                          nonlinear.correction = TRUE) {
  type <- match.arg(type)
  Z <- model_Z(data) ## Note: Will use model_Z() defined in the global workspace!
  fit_E <- Z$ZE %*% theta
  fit_V <- Z$ZV %*% theta
  if (is.null(Sigma_theta)) {
    ZE_var <- 0
    ZV_var <- 0
  } else {
    ZE_var <- rowSums(Z$ZE * (Z$ZE %*% Sigma_theta))
    ZV_var <- rowSums(Z$ZV * (Z$ZV %*% Sigma_theta))
  }
  if (type == "expectation") {
    fit <- fit_E
    sigma <- ZE_var^0.5
  } else if (type == "log-variance") {
    fit <- fit_V
    sigma <- ZV_var^0.5
  } else if (type == "observation") { ## observation predictions
    fit <- fit_E
    sigma <- (exp(fit_V + ZV_var / 2 * nonlinear.correction) + ZE_var)^0.5
  }
  q <- qt(1 - alpha / 2, df = df)
  lwr <- fit - q * sigma
  upr <- fit + q * sigma
  data.frame(mu = fit, sigma, lwr, upr)
}


#' General score functions
#' ----

# Input:
#   pred: data.frame with (at least) a column 'mu'
#   y: data vector
# Output:
#   A vector of Squared Error score values
score_se <- function(pred, y) {
  (y - pred$mu)^2
}

# Input:
#   pred: data.frame with (at least) columns 'mu' and 'sigma'
#   y: data vector
# Output:
#   A vector of Dawid-Sebastiani score values
score_ds <- function(pred, y) {
  ((y - pred$mu) / pred$sigma)^2 + 2 * log(pred$sigma)
}


#' Generate synthetic data for the 3D-printer problem
#' ----

generate_printer_data <- function(n) {
  data <-
    data.frame(actual = NA,
               cad = runif(n, 10, 300),
               colour = as.factor(sample(c('black', 'red', 'blue'),
                                         prob = c(0.5, 0.25, 0.25),
                                         size = n,
                                         replace = TRUE))
    )
  mean_values <-
    (data$colour == "black") * 2 +
    (data$colour == "red") * 1 +
    (data$colour == "blue") * 2 +
    data$cad * (
      (data$colour == "black") * (1.05) +
      (data$colour == "red") * (1.1) +
      (data$colour == "blue") * (0.9)
    )
  logv_values <-
    (data$colour == "black") * (2) +
    (data$colour == "red") * (1) +
    (data$colour == "blue") * (1) +
    data$cad * (
      (data$colour == "black") * (0.01) +
      (data$colour == "red") * (0.02) +
      (data$colour == "blue") * (0.02)
    )
  data$actual <- rnorm(n, mean = mean_values, sd = exp(logv_values / 2))

  data
}

set.seed(12345L)
printer_data <- generate_printer_data(100)

# # Simple data plotting.
# ggplot(printer_data) +
#   geom_point(aes(cad, actual, col = colour))
#
# # The default colours don't match the material colours, but we can override them:
# ggplot(printer_data) +
#   geom_point(aes(cad, actual, col = colour)) +
#   scale_color_manual(values = c(black = "black",
#                                 red = "red",
#                                 blue = "blue"))
