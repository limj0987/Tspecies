#' @title Correcting Species Divergence Time via Ks-Distribution Variance
#'
#' @description Tspecies uses a pre-fitted GAM model to predict Ne based on
#'     a given Ks value and a neutral mutation rate miu. The Ks refers to a vector
#'     representing d/L, where d is the base difference between sequences, and L
#'     is the sequence length. We will corrected the Ks value automatically by using
#'     multiple substitution and 2Ne based on our data.
#'
#' @param Ks A vector value representing d/L.
#' @param miu An numeric value representing neutral mutation rate in generation.
#' @param g An numeric value of generation time in years, with a default value of 1.
#' @return A numeric value representing the divergence time of species in generations.
#' @export
#'
#' @importFrom mgcv
#' @importFrom stats
#'
#' @examples
#' Tspecies(c(0.257, 0.202, 0.066, 0.197, 0.278, 0.089), 0.00000008)
#'
Tspecies <- function(Ks, miu, g) {
  mul_corrected_Ks <- (-3/4)*log(1-((4*Ks)/3))
  Ks_variance <- var(mul_corrected_Ks)
  Ks_mean <- mean(mul_corrected_Ks)
  package_path <- system.file(package = "Tspecies")
  if (miu >= 0.0000001) {
    model_file <- system.file("extdata/gam_model1_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model1
  } else if (miu < 0.0000001 & miu >= 0.00000009505) {
    model_file <- system.file("extdata/gam_model2_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model2
  } else if (miu < 0.00000009505 & miu >= 0.0000000901) {
    model_file <- system.file("extdata/gam_model3_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model3
  } else if (miu < 0.0000000901 & miu >= 0.00000008515) {
    model_file <- system.file("extdata/gam_model4_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model4
  } else if (miu < 0.00000008515 & miu >= 0.0000000802) {
    model_file <- system.file("extdata/gam_model5_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model5
  } else if (miu < 0.0000000802 & miu >= 0.00000007525) {
    model_file <- system.file("extdata/gam_model6_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model6
  } else if (miu < 0.00000007525 & miu >= 0.0000000703) {
    model_file <- system.file("extdata/gam_model7_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model7
  } else if (miu < 0.0000000703 & miu >= 0.00000006535) {
    model_file <- system.file("extdata/gam_model8_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model8
  } else if (miu < 0.00000006535 & miu >= 0.0000000604) {
    model_file <- system.file("extdata/gam_model9_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model9
  } else if (miu < 0.0000000604 & miu >= 0.00000005545) {
    model_file <- system.file("extdata/gam_model10_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model10
  } else if (miu < 0.00000005545 & miu >= 0.0000000505) {
    model_file <- system.file("extdata/gam_model11_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model11
  } else if (miu < 0.0000000505 & miu >= 0.00000004555) {
    model_file <- system.file("extdata/gam_model12_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model12
  } else if (miu < 0.00000004555 & miu >= 0.0000000406) {
    model_file <- system.file("extdata/gam_model13_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model13
  } else if (miu < 0.0000000406 & miu >= 0.00000003565) {
    model_file <- system.file("extdata/gam_model14_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model14
  } else if (miu < 0.00000003565 & miu >= 0.0000000307) {
    model_file <- system.file("extdata/gam_model15_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model15
  } else if (miu < 0.0000000307 & miu >= 0.00000002575) {
    model_file <- system.file("extdata/gam_model16_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model16
  } else if (miu < 0.00000002575 & miu >= 0.0000000208) {
    model_file <- system.file("extdata/gam_model17_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model17
  } else if (miu < 0.0000000208 & miu >= 0.00000001585) {
    model_file <- system.file("extdata/gam_model18_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model18
  } else if (miu < 0.00000001585 & miu >= 0.0000000109) {
    model_file <- system.file("extdata/gam_model19_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model19
  } else if (miu < 0.0000000109 & miu >= 0.00000000595) {
    model_file <- system.file("extdata/gam_model20_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model20
  } else if (miu < 0.00000000595 & miu > 0.00000000055) {
    model_file <- system.file("extdata/gam_model21_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model21
  } else if (miu <= 0.00000000055 & miu > 0.0000000001) {
    model_file <- system.file("extdata/gam_model22_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model22
  } else if (miu <= 0.0000000001) {
    model_file <- system.file("extdata/gam_model23_1.0.RData", package = "Tspecies", mustWork = TRUE)
    load(model_file)
    model <- gam_model23
  } else {
    stop("Invalid miu import.")
  }

  Ne_formula_calculated <- sqrt(Ks_variance / (16 * (miu^2)))
  Ne <- Ne_formula_calculated

  predict_Ne <- function(model, y_target, ne_range = c(1000, 1000000)) {
    objective <- function(ne_candidate) {
      pred <- NA
      tryCatch({
        pred <- predict.gam(model, newdata = data.frame(Ne = ne_candidate), type = "response")
      }, error = function(e){
        warning(paste("predict.gam error for Ne =", ne_candidate, ":", e$message))
      })

      if (is.na(pred)) {
        return(Inf)
      }
      return((pred - y_target)^2)
    }

    opt_result <- NULL
    tryCatch({
      opt_result <- optimize(f = objective, interval = ne_range)
    }, error = function(e) {
      warning(paste("Optimization failed:", e$message))
    })

    if (is.null(opt_result) || !is.finite(opt_result$objective)) {
      return(NA_real_)
    }
    return(opt_result$minimum)
  }

  Ne_2 <- predict_Ne(model, Ks_variance)

  if (is.na(Ne_2) || Ne_2 >= 900000) {
    cat("Info: Ne will be calculated based on the preset formula.\n")
  } else {
    Ne <- Ne_2 / 2
  }

  cat("The mutation rate is", miu, "per site per generation","\n")
  cat("The predicted number of effective population size is around", Ne, "\n")


  if (missing(g)) {
    tspc <- Ks_mean/(2*miu) - Ne_2
    corrected_Ks <- 2 * tspc * miu
    if (tspc <= 0) {
      tspc <- 0
    }
    cat("Estimated divegence time of the species:", tspc, "(generations)", "\n")
  } else {
    cat("The generation time of this species is", g, "years", "\n")
    tspc <- (Ks_mean/(2*miu) - Ne_2)*g
    corrected_Ks <- 2 * tspc * miu / g
    if (tspc <= 0) {
      tspc <- 0 }
    cat("Estimated divegence time of the species:", tspc, "(years)", "\n")
  }
  cat("The corrected_Ks is", corrected_Ks, "\n")
}




#' @title Create a Plot of Ks Distribution
#' @description A quick and simple way to create a plot of Ks distribution.
#' @param Ks A vector value representing d/L.
#' @param bin An numeric value of grouping, with a default value of 0.006.
#' @return A plot generated by ggplot2.
#' @export
#'
#' @importFrom ggplot2
#'
#' @examples
#' TspeciesPlot(c(0.257, 0.202, 0.066, 0.197, 0.278, 0.089), 0.06)

TspeciesPlot <- function(Ks, bin = 0.006) {
  df <- data.frame(Ks = Ks)
  ggplot(df, aes(x = Ks)) +
    geom_histogram(binwidth = bin, fill = "lightblue") +
    labs(x = "Ks", y = "Frequency") +
    theme_classic()
}
