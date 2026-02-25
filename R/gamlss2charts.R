predict_score <- function(object, ...) UseMethod("predict_score")

#' Get GAMLSS-based scores for new data, using gamlss2 fit
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
#' @return
#' @export
#'
#' @examples
predict_score.gamlss2 <-
  function(object, newdata, refdata = NULL,
           type = c("cent", "resid", "quantile", "conformal", "psr"),
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
      "quantile" = object$family$rqres(newdata[,feat], par = params),
      "conformal" = {
        cfs <- object$family$rqres(refdata[,feat], par = params)
        cfs2 <- object$family$rqres(newdata[,feat], par = params)
        (1 + sapply(cfs2, function(x) sum(x > cfs)))/(length(cfs) + 1)
      },
      "psr" = {
        2*object$family$cdf(q = newdata[,feat], par = params) - 1
      }
    )
  }

#' Get GAMLSS-based scores for new data, using gamlss fit
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
#' @return
#' @export
#'
#' @examples
predict_score.gamlss <-
  function(object, newdata, refdata = NULL,
           type = c("cent", "resid", "quantile", "conformal", "psr"),
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
      "cent" = do.call(get(paste0("p", object$family[1])), c(list(q = newdata[,feat]), params))
    )
  }
