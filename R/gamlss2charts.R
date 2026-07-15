#' Get GAMLSS-based scores for new data, using gamlss2 fit
#'
#' Compute scores on new data based on a gamlss2/gamlss fit. Gets scores for
#' new data input as `newdata`. If specified, fits an adjustment model for
#' the new dataset (e.g. adjusting for batch effects). Adjustment model is fit
#' using reference data input as `refdata` and specified in `newformula`.
#'
#' @param object Object of class `gamlss2`, `gamlss`, or a `list` of params.
#'   See details
#' @param newdata Data frame containing new observations
#' @param refdata Data frame containing reference observations, typically
#'   healthy individuals. Can overlap with `newdata`. If NULL, uses newdata
#' @param rm.term String indicating factor to remove from original model fit,
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
#' but workarounds are available for \link[gamlss]{random} effects. Lists of
#' parameters can be provided as `object`. The list needs to have names `new`
#' and/or `ref`, where parameters correspond to `newdata` and `refdata`
#' respectively. `feat` and `family` need to be specified, and should match with
#' the `newdata` and distribution used to fit the original parameters.
#'
#' For `type="cent"`, this is called out-of-sample centile scoring.
#' `type="resid"` returns the residuals and `type="zscore"` returns z-scores.
#' `type="quantile"` returns the quantile residuals (matching centile scores to
#' an N(0,1) distribution). `type="parameter"` returns the fitted parameters.
#'
#' @rdname predict_score
#' @examples
#' train <- iris[1:100,]
#' train$Species <- droplevels(train$Species)
#' ft <- gamlss2(Sepal.Length ~ Sepal.Width + Species | Species,
#'   family = BCCG(), data = train)
#' predict_score(ft, iris[-(1:100),], rm.term = "Species") |> hist()
#'
#' # example for gamlss
#' if (require("gamlss")) {
#'   train <- iris[1:100,]
#'   train$Species <- droplevels(train$Species)
#'   ft <- gamlss(Sepal.Length ~ Sepal.Width + Species, ~ Species,
#'     family = BCCG(), data = train)
#'   predict_score(ft, iris[-(1:100),], rm.term = "Species") |> hist()
#' }
#'
#' @references Bethlehem, R. A. I., Seidlitz, J., White, S. R., Vogel, J. W., Anderson, K. M., Adamson, C., Adler, S., Alexopoulos, G. S., Anagnostou, E., Areces-Gonzalez, A., Astle, D. E., Auyeung, B., Ayub, M., Bae, J., Ball, G., Baron-Cohen, S., Beare, R., Bedford, S. A., Benegal, V., … Alexander-Bloch, A. F. (2022). Brain charts for the human lifespan. Nature, 604(7906), Article 7906. https://doi.org/10.1038/s41586-022-04554-y
#'
#' Dinga, R., Fraza, C. J., Bayer, J. M. M., Kia, S. M., Beckmann, C. F., & Marquand, A. F. (2021). Normative modeling of neuroimaging data using generalized additive models of location scale and shape (p. 2021.06.14.448106). bioRxiv. https://doi.org/10.1101/2021.06.14.448106
# S3 generic. Dispatches on the class of `object` to one of the three methods
# below (.gamlss2, .gamlss, or .list). `...` forwards all other arguments.
predict_score <- function(object, ...) UseMethod("predict_score")

#' @rdname predict_score
#' @export
# ---- Method for gamlss2 fits -------------------------------------------------
predict_score.gamlss2 <-
  function(object, newdata, refdata = NULL,
           # `type` default is the first element, "cent" (centile score).
           type = c("cent", "resid", "zscore", "quantile", "parameter"),
           adjust = TRUE, rm.term = NULL,
           # gamlss2 formula syntax: `|` separates the model for each distribution
           # parameter. offset(mu)/offset(sigma) freeze the original model's link-
           # scale predictions as fixed offsets in the adjustment model.
           newformula = y ~ offset(mu) | offset(sigma),
           which.params = c("mu", "sigma")) {
    # Validate `type` against the allowed set and reduce to a single value.
    type = match.arg(type)

    # A z-score (y - mu)/sigma is only a true N(0,1) deviate for the Normal
    # family ("NO"); warn if the user requests it for any other distribution.
    if(object$family$family != "NO" & type == "zscore") {
      warning("Z-scores may not be valid for families other than NO")
    }

    # Response variable name = first variable in the model formula.
    feat <- all.vars(formula(object))[1]
    # Terms to predict from: all formula variables except the response and the
    # term being removed. Excluding `rm.term` here is how its effect (e.g. a
    # batch / site random effect) is zeroed out of the predictions.
    mterms <- c("Intercept", setdiff(all.vars(formula(object)), c(feat, rm.term)))
    # Turn parameter names into named integer column indices (mu=1, sigma=2,
    # nu=3, tau=4); used later to pick which parameter columns get adjusted.
    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    # If no reference data is given, use newdata as its own reference.
    if (is.null(refdata)) {
      refdata <- newdata
    }

    # match rm.term factor levels to original fit
    # newdata/refdata may contain factor levels the model never saw (e.g. a new
    # site). predict() would error on unseen levels, so append the new levels to
    # the fit's stored xlevels and reorder refdata's levels to match. The effect
    # is still excluded from prediction via `mterms`. - UPDATED
    if (!is.null(rm.term)) {
      oldlevels <- object$xlevels[[rm.term]]
      newlevels <- setdiff(levels(refdata[[rm.term]]), oldlevels)
      object$xlevels[[rm.term]] <- c(oldlevels, newlevels)
    }

    if (adjust) {
      refdata$y <- refdata[[feat]]  # `newformula`'s LHS is `y`, so copy the response into `y`
      # get offsets for refdata and fit adjustment model
      # Original model's link-scale predictions on refdata (rm.term excluded).
      # These are the frozen offsets fed to the adjustment model. ####HOW DOES THIS WORK IF refdata <- newdata???
      pred2 <- predict(object, newdata = refdata, type = "link", terms = mterms)
      refdata <- cbind(refdata, pred2)  # add mu/sigma... columns for the offset() terms
      # Fit the adjustment model. Because the formula is all offsets, fit2 only
      # estimates the residual shift (typically a per-parameter intercept) on top
      # of the frozen offsets -- i.e. how this reference sample deviates from `object`.
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = object$family)

      # get offsets for newdata
      # Same offset extraction for the observations we actually want to score.
      newdata$y <- newdata[[feat]]
      pred <- predict(object, newdata = newdata, type = "link", terms = mterms)
      newdata <- cbind(newdata, pred)

      # apply fit2 estimates, which are shifts to the original parameters
      # Adjustment model's link-scale prediction on newdata.
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0  # zero out any parameters not requested for adjustment
      # Combine frozen offsets with estimated shifts on the link scale, then
      # map2par() inverts the links back to natural parameter scale.
      # NOTE (potential bug): fit2's linear predictor already includes the
      # offset() terms, so predict(fit2, type="link") returns offset + shift.
      # Adding `pred` again here may double-count the offset. Flagging only.
      params <- family(fit2)$map2par(pred + shift)
    } else {
      # No adjustment: ask the original fit directly for natural-scale parameters
      # (equivalent to predict.gamlss2).
      params <- predict(object, newdata = newdata, type = "parameter", terms = mterms)
    }

    # Compute the requested score using the family's own distribution functions.
    switch(
      type,
      "cent" = object$family$cdf(q = newdata[,feat], par = params),    # CDF at y -> centile in [0,1]
      "resid" = newdata[,feat] - params[,1],                           # raw residual y - mu
      "zscore" = (newdata[,feat] - params[,1])/params[,2],             # (y - mu)/sigma
      "quantile" = object$family$rqres(newdata[,feat], par = params),  # (randomized) quantile residuals on N(0,1) scale
      "parameter" = params                                             # the fitted parameters themselves
    )
  }

#' @rdname predict_score
#' @export
# ---- Method for gamlss fits (older gamlss package) ---------------------------
# Same overall logic as the gamlss2 method, but adapted to the gamlss object
# layout and API: predictAll() instead of predict(type="parameter"), a character
# `family` vector, per-parameter link functions, and distribution functions
# looked up by name (pNO, pBCCG, ...).
predict_score.gamlss <-
  function(object, newdata, refdata = NULL,
           type = c("cent", "resid", "zscore", "quantile", "parameter"),
           adjust = TRUE, rm.term = NULL,
           newformula = y ~ offset(mu) | offset(sigma),
           which.params = c("mu", "sigma")) {
    type = match.arg(type)

    # In gamlss, $family is a character vector; [1] is the family code (e.g. "NO").
    if(object$family[1] != "NO" & type == "zscore") {
      warning("Z-scores may not be valid for families other than NO")
    }

    # Response name: LHS of the mu formula. mu.formula[[2]] is the response symbol.
    feat <- as.character(object$mu.formula[[2]])
    # Terms to predict from: columns of the stored model frame, minus response and rm.term.
    
    ###EDIT: object$model returns NULL, borrowing code from gamlss2 method - may
    #need to update to list_predictors depending on robustness to smooths, models 
    #saved elsewhere, etc
    mterms <- c("Intercept", setdiff(all.vars(formula(object)), c(feat, rm.term)))
    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    if (is.null(refdata)) {
      refdata <- newdata
    }

    # Append unseen factor levels for rm.term.
    if (!is.null(rm.term)) {
      # rm.term may enter several distribution parameters ("moments"); gamlss
      # stores levels per moment as <param>.xlevels. Update the stored levels in
      # every moment where rm.term appears so predict() accepts the new site.
      updated <- FALSE
      for (p in object$parameters) {                      # e.g. c("mu","sigma")
        slot <- paste0(p, ".xlevels")
        if (rm.term %in% names(object[[slot]])) {
          oldlevels <- object[[slot]][[rm.term]]
          newlevels <- setdiff(levels(newdata[[rm.term]]), oldlevels)
          object[[slot]][[rm.term]] <- c(oldlevels, newlevels)
          updated <- TRUE
        }
      }
      stopifnot(updated)                                  # rm.term not found in any moment
    }

    if (adjust) {
      # get offsets for refdata and fit adjustment model
      # predictAll() predicts all distribution parameters at once; subset by
      # object$parameters (the fitted parameter names) to keep them in order.
      pred2 <- predictAll(object, newdata = refdata, type = "link", terms = mterms)[object$parameters]
      refdata <- cbind(refdata, pred2)     # add mu/sigma... columns for the offsets
      refdata$y <- refdata[[feat]]         # `newformula` LHS is `y`
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = object$family[1])

      # get offsets for newdata
      pred <- predictAll(object, newdata = newdata, type = "link", terms = mterms)[object$parameters]
      newdata <- cbind(newdata, pred)
      # POTENTIAL BUG: assigns the REFDATA response into newdata$y (should almost
      # certainly be newdata[[feat]]). It happens to be harmless in this method
      # because newdata$y is never read afterwards (predict(fit2) uses the offset
      # columns, and the final switch uses newdata[,feat]) -- but it is wrong, and
      # will error if refdata and newdata differ in number of rows.
      # newdata$y <- refdata[[feat]] #commented out for testing

      # apply fit2 estimates, which are shifts to the original parameters
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0
      # data.frame(pred): the offset columns as a data.frame so it adds elementwise
      # to `shift`. Same potential offset double-counting note as the gamlss2 method.
      params <- family(fit2)$map2par(data.frame(pred) + shift)
    } else {
      # No adjustment. predictAll(type="link") returns link-scale values, so invert
      # each parameter's link by hand using the stored mu.link/sigma.link/... .
      links <- predictAll(object, newdata = newdata, type = "link", terms = mterms)[which.params]
      params <- lapply(setNames(nm = c("mu", "sigma", "nu", "tau")[which.params]), function(x) {
        link <- object[[paste0(x, ".link")]]
        switch(
          link,
          "identity" = links[[x]],                       # inverse of identity link
          "log" = exp(links[[x]]),                        # inverse of log link
          "logit" = exp(links[[x]])/(1+exp(links[[x]]))   # inverse of logit link
        )
      })
      # NOTE: here `params` is a named LIST (one element per parameter), whereas in
      # the adjust=TRUE branch it is a matrix/data.frame from map2par(). That matters
      # for the switch() below -- see the "resid"/"zscore" cases.
    }

    # Score the observations. gamlss has no unified family object, so the CDF is
    # looked up by name: get(paste0("p", family)) e.g. pBCCG, called via do.call
    # with the parameter list.
    switch(
      type,
      "cent" = do.call(get(paste0("p", object$family[1])), c(list(q = newdata[,feat]), params)),
      # POTENTIAL BUG: params[,1] / params[,2] assume a matrix/data.frame. That
      # holds in the adjust=TRUE branch, but in the adjust=FALSE branch params is a
      # list, so params[,1] would error. "resid"/"zscore" with adjust=FALSE look broken.
      "resid" = newdata[,feat] - params[,1],                          # raw residual y - mu
      # quantile residuals: push the CDF value through qnorm() to get an N(0,1)-scale value.
      "quantile" = qnorm(do.call(get(paste0("p", object$family[1])), c(list(q = newdata[,feat]), params))),
      "zscore" = (newdata[,feat] - params[,1])/params[,2],            # (y - mu)/sigma
      "parameter" = params
    )
  }

#' @rdname predict_score
#' @export
# ---- Method for a plain list of pre-computed parameters -----------------------
# Use when you have no fitted model object, only parameters. Pass a list with
# $new (params for newdata) and optionally $ref (params for refdata), plus `feat`
# and `family` explicitly, since there is no model object to read them from.
predict_score.list <-
  function(object, newdata, refdata = NULL, feat, family,
           type = c("cent", "resid", "zscore", "quantile", "parameter"),
           adjust = TRUE, rm.term = NULL,
           newformula = y ~ offset(mu) | offset(sigma),
           which.params = c("mu", "sigma")) {
    type = match.arg(type)

    if(family != "NO" & type == "zscore") {
      warning("Z-scores may not be valid for families other than NO")
    }

    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    # If no refdata, reuse newdata and use the supplied $new params as the reference
    # offsets. POTENTIAL BUG: `pred2` is only assigned inside this branch. If refdata
    # IS supplied, pred2 is never set, so cbind(refdata, pred2) below errors. This
    # likely should fall back to object$ref when refdata is provided.
    if (is.null(refdata)) {
      refdata <- newdata
      pred2 <- object$new
    }

    newdata$y <- newdata[[feat]]  # `y` columns for the trivial fit / offsets
    refdata$y <- refdata[[feat]]
    pred <- object$new            # supplied link-scale params for newdata
    # Fit a trivial intercept-only model purely to obtain a family object exposing
    # the $map2par / $cdf / $rqres helpers used below.
    fit <- gamlss2::gamlss2(y ~ 1, data = newdata, family = family)

    if (adjust) {
      # get offsets for refdata and fit adjustment model
      refdata <- cbind(refdata, pred2)
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = family)

      # get offsets for newdata
      newdata <- cbind(newdata, pred)

      # apply fit2 estimates, which are shifts to the original parameters
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0
      # Same potential offset double-counting note as the other methods.
      params <- fit2$family$map2par(pred + shift)
    } else {
      # No adjustment: map the supplied link-scale params straight to natural scale.
      params <- fit$family$map2par(pred)
    }

    # Same scoring switch as the gamlss2 method, using the trivial fit's family object.
    switch(
      type,
      "cent" = fit$family$cdf(q = newdata[,feat], par = params),
      "resid" = newdata[,feat] - params[,1],
      "zscore" = (newdata[,feat] - params[,1])/params[,2],
      "quantile" = fit$family$rqres(newdata[,feat], par = params),
      "parameter" = params
    )
  }
