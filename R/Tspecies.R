#' @title Correcting Species Divergence Time via Ks-Distribution Variance
#'
#' @description Tspecies uses a pre-fitted GAM model to predict Ne based on
#'     a given Ks value and a neutral mutation rate miu. The Ks refers to a vector
#'     representing d/L, where d is the base difference between sequences, and L
#'     is the sequence length. We will corrected the Ks value automatically by using
#'     multiple substitution and 2Ne based on our data.
#'
#' @param Ks A vector value representing d/L.
#' @param miu An numeric value representing neutral mutation rate.
#' @return A numeric value representing the divergence time of species in generations.
#' @export
#'
#' @importFrom mgcv predict.gam gam
#' @importFrom stats optimize
#'
#' @examples
#' Tspecies(c(0.257, 0.202, 0.066, 0.197, 0.278, 0.089), 0.00000008)
#'
Tspecies <- function(Ks, miu) {
  corrected_Ks <- (-3/4)*log(1-((4*Ks)/3))
  Ks_variance <- var(corrected_Ks)
  Ks_mean <- mean(corrected_Ks)
  package_path <- system.file(package = "Tspecies")
  if (miu >= 0.0000001) {
    model_file <- system.file("extdata/gam_model1.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model1
  } else if (miu < 0.0000001 & miu >= 0.00000009505) {
    model_file <- system.file("extdata/gam_model2.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model2
  } else if (miu < 0.00000009505 & miu >= 0.0000000901) {
    model_file <- system.file("extdata/gam_model3.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model3
  } else if (miu < 0.0000000901 & miu >= 0.00000008515) {
    model_file <- system.file("extdata/gam_model4.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model4
  } else if (miu < 0.00000008515 & miu >= 0.0000000802) {
    model_file <- system.file("extdata/gam_model5.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model5
  } else if (miu < 0.0000000802 & miu >= 0.00000007525) {
    model_file <- system.file("extdata/gam_model6.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model6
  } else if (miu < 0.00000007525 & miu >= 0.0000000703) {
    model_file <- system.file("extdata/gam_model7.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model7
  } else if (miu < 0.0000000703 & miu >= 0.00000006535) {
    model_file <- system.file("extdata/gam_model8.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model8
  } else if (miu < 0.00000006535 & miu >= 0.0000000604) {
    model_file <- system.file("extdata/gam_model9.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model9
  } else if (miu < 0.0000000604 & miu >= 0.00000005545) {
    model_file <- system.file("extdata/gam_model10.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model10
  } else if (miu < 0.00000005545 & miu >= 0.0000000505) {
    model_file <- system.file("extdata/gam_model11.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model11
  } else if (miu < 0.0000000505 & miu >= 0.00000004555) {
    model_file <- system.file("extdata/gam_model12.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model12
  } else if (miu < 0.00000004555 & miu >= 0.0000000406) {
    model_file <- system.file("extdata/gam_model13.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model13
  } else if (miu < 0.0000000406 & miu >= 0.00000003565) {
    model_file <- system.file("extdata/gam_model14.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model14
  } else if (miu < 0.00000003565 & miu >= 0.0000000307) {
    model_file <- system.file("extdata/gam_model15.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model15
  } else if (miu < 0.0000000307 & miu >= 0.00000002575) {
    model_file <- system.file("extdata/gam_model16.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model16
  } else if (miu < 0.00000002575 & miu >= 0.0000000208) {
    model_file <- system.file("extdata/gam_model17.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model17
  } else if (miu < 0.0000000208 & miu >= 0.00000001585) {
    model_file <- system.file("extdata/gam_model18.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model18
  } else if (miu < 0.00000001585 & miu >= 0.0000000109) {
    model_file <- system.file("extdata/gam_model19.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model19
  } else if (miu < 0.0000000109 & miu >= 0.00000000595) {
    model_file <- system.file("extdata/gam_model20.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model20
  } else if (miu < 0.00000000595 & miu > 0.00000000055) {
    model_file <- system.file("extdata/gam_model21.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model21
  } else if (miu <= 0.00000000055 & miu > 0.0000000001) {
    model_file <- system.file("extdata/gam_model22.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model22
  } else if (miu <= 0.0000000001) {
    model_file <- system.file("extdata/gam_model23.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model23
  } else {
    stop("Invalid import.")
  }

  predict_Ne <- function(model, y_target, ne_range = c(1000, 1000000)) {
    objective <- function(ne) {
      pred <- predict.gam(model, newdata = data.frame(Ne = ne), type = "response")
      if (is.na(pred)) stop("The result of prediction is NA, please check the scope of input.")
      return((pred - y_target)^2)
    }
    result <- optimize(f = objective, interval = ne_range)
    return(result$minimum)
  }

  Ne <- predict_Ne(model, Ks_variance)
  cat("The predicted number of effective population size is around", Ne, "\n")
  tspc <- Ks_mean/(2*miu) - Ne
  if (tspc <= 0) {
    tspc <- 0
  }
  cat("Estimated divegence time of the species:", tspc, "(in generations)", "\n")
}
