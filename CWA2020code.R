#' ---
#' title: Statistical Computing, Coursework A, Supplementary code
#' author: Finn Lindgren
#' date: 4/2/2020
#' output: pdf_document
#' ---

#' See the Prediction and Proper Scoring Rules longform notes
#' for examples of how to use these functions.

# To access these functions inside other R files, use
#   source("CWA2020code.R")

#' Functions for linear expectation and log-variance models
#' ----

#' model_Z() constructs the Z-matrices for models of the type used in lab 3,
#' extended to allow different formulas for expectation and variance parts
#' of the model.

# Input:
#   formulas: A list with two named elements:
#     E: The formula for the expectation part, e.g. ~ 1 + x
#     V: The formula for the log-variance part, e.g. ~ 1 + x
#   data: A data.frame of covariate values used by the formulas
# Output:
#   A list with two named elements:
#     ZE: The matrix defining the link between linear
#         predictor and parameters for the expectations
#     ZV: The matrix defining the link between linear
#         predictor and parameters for the log-variances
model_Z <- function(formulas, data) {
  ZE <- model.matrix(formulas$E, data = data)
  ZV <- model.matrix(formulas$V, data = data)
  list(ZE = cbind(ZE, ZV * 0), ZV = cbind(ZE * 0, ZV))
}

#' neg_log_lik() computes the negative log-likelihood for the general model
#' class from Lab 3.
# Input:
#   theta: parameter vector of length ncol(Z$ZE)
#   Z: a list of the type returned by model_Z()
#   y: the response data vector
# Output:
#   A negative-log-likelihood value
neg_log_lik <- function(theta, Z, y) {
  -sum(dnorm(y, mean = Z$ZE %*% theta, sd = exp(Z$ZV %*% theta)^0.5, log = TRUE))
}

#' model_predict() computes either point predictions and model component confidence
#' intervals, or data predictions.

# Input:
#   theta: parameter vector
#   formulas: a list of formulas, of the kind needed by model_Z()
#   Sigma_theta: The parameter estimation error covariance matrix
#                If NULL, theta is treated as being a fixed, known parameter vector
#   newdata: covariate information, in whatever format the model_Z() wants it.
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
model_predict <- function(theta, formulas, Sigma_theta = NULL,
                          newdata,
                          type = c("expectation", "log-variance", "observation"),
                          alpha = 0.05, df = Inf,
                          nonlinear.correction = TRUE) {
  type <- match.arg(type)
  Z <- model_Z(formulas, data = newdata)
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

# Input:
#   pred: data.frame with (at least) columns 'lwr' and 'upr'
#   y: data vector
#   alpha: The target probability mass to fall outside the intervals
# Output:
#   A vector of Interval score values
score_interval <- function(pred, y, alpha = 0.05) {
  L <- pred$lwr
  U <- pred$upr
  U - L + 2 / alpha * (pmax(0, L - y) + pmax(0, y - U))
}

