#' Get GAMLSS-based scores for new data, using gamlss2 fit
#'
#' Compute scores on new data based on a gamlss2/gamlss fit. Gets scores for
#' new data input as `newdata`. If specified, fits an adjustment model for
#' the new dataset (e.g. adjusting for batch effects). Adjustment model is fit
#' using reference data input as `refdata` and specified in `newformula`.
#'
#' @param object Object of class `gamlss2`, typically output of
#'   \link[gamlss2]{gamlss2}
#' @param newdata Data frame containing new observations
#' @param refdata Data frame containing reference observations, typically
#'   healthy individuals. Can overlap with `newdata`. If NULL, uses newdata
#' @param rm.term String indicating term to remove from original model fit,
#'   usually a batch variable (e.g. setting random effects to zero)
#' @param adjust Whether to fit adjustment parameters on the new data (e.g.
#'   getting batch adjustments for a new site)
#' @param newformula If `adjust=TRUE`, adjustment fit uses this formula. Need to
#'   include offset terms for desired parameter, e.g. y ~ offset(mu) uses fitted
#'   mu values from `object` fit, y ~ offset(mu) | offset(sigma) uses fitted
#'   mu and sigma values from `object`. See Details
#' @param which.params Parameters to adjust
#' @param type Type of score to compute. See Details
#'
#' @return Vector of scores of length equal to number of rows in `newdata`
#' @export
#'
#' @details
#' When `adjust=FALSE`, essentially the same as \link[gamlss2]{predict.gamlss2}.
#' When `adjust=TRUE`, fixes the predictions from `object` as offsets for an
#' adjustment model. These offsets need to be specified in `newformula`. Then,
#' the adjustment model is fit, adjustment parameters are combined with the
#' offsets for `which.params`, and scores are computed. See references for more
#' details. \link[gamlss]{gamlss} support is still limited to fixed effects,
#' but workarounds are available for \link[gamlss]{random} effects.
#'
#' For `type="cent"`, this is called out-of-sample centile scoring.
#' `type="resid"` returns the studentized residuals and `type="quantile"`
#' returns the quantile residuals (matching quantiles to an N(0,1) distribution).
#'
#' @rdname predict_score
#' @examples
#' ft <- gamlss2(Sepal.Length ~ Sepal.Width + Species | Species,
#'   family = BCCG(), data = iris[1:100,])
#' predict_score(ft, iris[-(1:100),], rm.term = "Species") |> hist()
#'
#' # example for gamlss
#' if (require("gamlss")) {
#'   ft <- gamlss(Sepal.Length ~ Sepal.Width + Species, ~ Species,
#'     family = BCCG(), data = iris[1:100,])
#'   predict_score(ft, iris[-(1:100),], rm.term = "Species") |> hist()
#' }
#'
#' @references Dinga, R., Fraza, C. J., Bayer, J. M. M., Kia, S. M., Beckmann, C. F., & Marquand, A. F. (2021). Normative modeling of neuroimaging data using generalized additive models of location scale and shape (p. 2021.06.14.448106). bioRxiv. https://doi.org/10.1101/2021.06.14.448106
#'
#' Bethlehem, R. A. I., Seidlitz, J., White, S. R., Vogel, J. W., Anderson, K. M., Adamson, C., Adler, S., Alexopoulos, G. S., Anagnostou, E., Areces-Gonzalez, A., Astle, D. E., Auyeung, B., Ayub, M., Bae, J., Ball, G., Baron-Cohen, S., Beare, R., Bedford, S. A., Benegal, V., … Alexander-Bloch, A. F. (2022). Brain charts for the human lifespan. Nature, 604(7906), Article 7906. https://doi.org/10.1038/s41586-022-04554-y
predict_score <- function(object, ...) UseMethod("predict_score")

#' @rdname predict_score
#' @export
predict_score.gamlss2 <-
  function(object, newdata, refdata = NULL,
           type = c("cent", "resid", "quantile"),
           adjust = TRUE, rm.term = NULL,
           newformula = y ~ offset(mu) | offset(sigma),
           which.params = c("mu", "sigma")) {
    type = match.arg(type)

    feat <- names(object$model)[1]
    mterms <- c("Intercept", setdiff(names(object$model), c(feat, rm.term)))
    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    if (is.null(refdata)) {
      refdata <- newdata
    }

    if (!is.null(rm.term)) {
      oldlevels <- object$xlevels[[rm.term]]
      newlevels <- setdiff(levels(refdata[[rm.term]]), oldlevels)
      object$xlevels[[rm.term]] <- c(oldlevels, newlevels)
    }

    if (adjust) {
      refdata$y <- refdata[[feat]]
      # get offsets for refdata and fit adjustment model
      pred2 <- predict(object, newdata = refdata, type = "link", terms = mterms)
      refdata <- cbind(refdata, pred2)
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = object$family)

      # get offsets for newdata
      newdata$y <- refdata[[feat]]
      pred <- predict(object, newdata = newdata, type = "link", terms = mterms)
      newdata <- cbind(newdata, pred)

      # apply fit2 estimates, which are shifts to the original parameters
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0
      params <- family(fit2)$map2par(pred + shift)
    } else {
      params <- predict(object, newdata = newdata, type = "parameter", terms = mterms)
    }

    switch(
      type,
      "cent" = object$family$cdf(q = newdata[,feat], par = params),
      "resid" = (newdata[,feat] - params[,1])/params[,2],
      "quantile" = object$family$rqres(newdata[,feat], par = params)
    )
  }

#' @rdname predict_score
#' @export
predict_score.gamlss <-
  function(object, newdata, refdata = NULL,
           type = c("cent", "resid"),
           adjust = TRUE, rm.term = NULL,
           newformula = y ~ offset(mu) | offset(sigma),
           which.params = c("mu", "sigma")) {
    type = match.arg(type)

    feat <- as.character(object$mu.formula[[2]])
    mterms <- c("Intercept", setdiff(names(object$model), c(feat, rm.term)))
    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    if (is.null(refdata)) {
      refdata <- newdata
    }

    if (!is.null(rm.term)) {
      oldlevels <- object$xlevels[[rm.term]]
      newlevels <- setdiff(levels(refdata[[rm.term]]), oldlevels)
      object$xlevels[[rm.term]] <- c(oldlevels, newlevels)
    }

    if (adjust) {
      # get offsets for refdata and fit adjustment model
      pred2 <- predictAll(object, newdata = refdata, type = "link", terms = mterms)[object$parameters]
      refdata <- cbind(refdata, pred2)
      refdata$y <- refdata[[feat]]
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = object$family[1])

      # get offsets for newdata
      pred <- predictAll(object, newdata = newdata, type = "link", terms = mterms)[object$parameters]
      newdata <- cbind(newdata, pred)
      newdata$y <- refdata[[feat]]

      # apply fit2 estimates, which are shifts to the original parameters
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0
      params <- family(fit2)$map2par(data.frame(pred) + shift)
    } else {
      links <- predictAll(object, newdata = newdata, type = "link", terms = mterms)[which.params]
      params <- lapply(setNames(nm = c("mu", "sigma", "nu", "tau")[which.params]), function(x) {
        link <- object[[paste0(x, ".link")]]
        switch(link,
               "identity" = links[[x]],
               "log" = exp(links[[x]]))
      })
    }

    switch(
      type,
      "cent" = do.call(get(paste0("p", object$family[1])), c(list(q = newdata[,feat]), params)),
      "resid" = (newdata[,feat] - params[,1])/params[,2]
    )
  }
