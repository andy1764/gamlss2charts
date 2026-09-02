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
#'   mu and sigma values from `object`. If `NULL` (the default) AND
#'   `which.params` is also `NULL`, it is derived automatically. See Details
#' @param which.params Parameters to adjust. If `NULL` (the default) AND
#'   `newformula` is also `NULL`, defaults to the moments that contain `rm.term`.
#' @param type Type of score to compute. See Details
#' @param traindata Data originally used to fit `object` (`gamlss` method only,
#'   see Details).
#'
#' @return Vector of scores of length equal to number of rows in `newdata`
#' @export
#'
#' @details
#' When `adjust=FALSE`, essentially the same as \link[gamlss2]{predict.gamlss2}.
#' When `adjust=TRUE`, fixes the predictions from `object` as offsets for an
#' adjustment model. These offsets need to be specified in `newformula`. Then,
#' the adjustment model is fit, adjustment parameters are combined with the
#' offsets for `which.params`, and scores are computed. By default, both the
#' batch offset formula (`newformula`) and parameters to adjust (`which.params`)
#' are derived automatically from model object. Supplying either `newformula` or
#' `which.params` disables auto-derivation and the other falls back to its legacy
#' mu/sigma default (NOTE: should be updated!). See references for more details.
#'
#' For \link[gamlss]{gamlss} fits with `adjust=TRUE`, the batch term
#' `rm.term` may be a fixed factor, a smooth, or a \link[gamlss]{random} effect:
#' offsets set `rm.term` to its baseline (the population mean for a centered
#' random effect, the reference level for a fixed factor) and its site-specific
#' shift is re-estimated by the adjustment model. Because the adjustment model's
#' intercept absorbs any constant baseline difference, the choice of baseline does
#' not affect `adjust=TRUE` scores. (With `adjust=FALSE`, random effects are not
#' supported, since predicting an unseen level returns `NA`.) If passing data with
#' levels already included in the training data, `rm.term` and `adjust` will
#' throw out the model's fitted batch estimate and recompute the offset.
#'
#' By default a \link[gamlss]{gamlss} fit is scored WITHOUT the original fitting
#' data whenever possible: batch-baseline offsets are rebuilt from the parametric coefficients,
#' the stored spline coefficients and interpolation functions, and the stored
#' per-level random-effect BLUPs (levels unseen in the fit default to the
#' population value 0). This is exact for `adjust=TRUE` (any fixed-factor baseline
#' difference is absorbed by the adjustment model). Some \link[gamlss]{gamlss} smoother types
#' (`cs()`, `ps()`, `ga()`/`s()`, ...) are not supported, so models with such smooths
#' outside `rm.term` (which is dropped) need `traindata`. (This concerns only the \link[gamlss]{gamlss}
#' method; \link[gamlss2]{gamlss2} fits predict their `s()`/`ga()` smooths without the
#' original data natively.)
#'
#' Lists of parameters can be provided as `object`. The list needs to have names `new`
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
predict_score.gamlss2 <-
  function(object, newdata, refdata = NULL,
           # `type` default is the first element, "cent" (centile score).
           type = c("cent", "resid", "zscore", "quantile", "parameter"),
           adjust = TRUE, rm.term = NULL,
           newformula = NULL,
           which.params = NULL) {
    type = match.arg(type)

    if (!is.null(rm.term) && !adjust) {
      warning("`rm.term` is set with `adjust = FALSE`: the '", rm.term, "' effect ",
              "is removed but not re-estimated/batch-adjusted for the new data. For a random ",
              "effect this sets it to the population mean (0); for a fixed factor ",
              "it sets it to the reference level. Use `adjust = TRUE` to estimate ",
              "and apply a batch adjustment.")
    }

    if(object$family$family != "NO" & type == "zscore") {
      warning("Z-scores may not be valid for families other than NO")
    }

    feat <- all.vars(formula(object))[1]
    # Terms to predict from: all variables except feat and rm.term
    mterms <- c("Intercept", setdiff(all.vars(formula(object)), c(feat, rm.term)))
    # Resolve the adjustment spec: honor explicit args, else auto-derive from the
    # moments containing rm.term (adjust those, freeze the rest).
    spec <- .resolve_adjust_spec(newformula, which.params,
                                 .moment_formulas_gamlss2(object), rm.term)
    newformula   <- spec$newformula
    which.params <- spec$which.params
    # Turn parameter names into named integer column indices (mu=1, sigma=2,
    # nu=3, tau=4); used later to pick which parameter columns get adjusted.
    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    # If no reference data is given, use newdata as its own reference.
    if (is.null(refdata)) {
      refdata <- newdata
    } else {
      #confirm all newdata levels are in refdata
      stopifnot("Not all levels of rm.term in newdata are present in refdata " = levels(newdata[[rm.term]]) %in% levels(refdata[[rm.term]]))
    }

    # match rm.term factor levels to original fit, appending unseen levels
    # to the fit's stored xlevels and reorder refdata's levels to match. The effect
    # is still excluded from prediction via `mterms`
    if (!is.null(rm.term)) {
      for (par in c("mu", "sigma", "nu", "tau")) {
        oldlevels <- object$xlevels[[par]][[rm.term]]
        newlevels <- setdiff(levels(refdata[[rm.term]]), oldlevels)
        object$xlevels[[par]][[rm.term]] <- c(oldlevels, newlevels)
      }
    }

    if (adjust) {
      # The adjustment model estimates ONE batch shift for all of `refdata`, so
      # `rm.term` must be a single batch level shared by newdata and refdata --
      # multiple levels would be silently pooled into one (averaged) adjustment.
      if (!is.null(rm.term)) {
        new_lv <- unique(as.character(newdata[[rm.term]]))
        ref_lv <- unique(as.character(refdata[[rm.term]]))
        if (length(new_lv) != 1L || length(ref_lv) != 1L || !identical(new_lv, ref_lv))
          stop("with adjust = TRUE, `rm.term` ('", rm.term, "') must be a single ",
               "batch level shared by newdata and refdata (the adjustment estimates ",
               "one shift for that batch). newdata level(s): ",
               paste(new_lv, collapse = ", "), "; refdata level(s): ",
               paste(ref_lv, collapse = ", "), ".")
      }
      refdata$y <- refdata[[feat]]  # `newformula`'s LHS is `y`, so copy the response into `y`
      # get offsets for refdata and fit adjustment model.
      pred2 <- predict(object, newdata = refdata, type = "link", terms = mterms)
      refdata <- cbind(refdata, pred2)  # add mu/sigma... columns for the offset() terms
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = object$family)

      # get offsets for newdata
      newdata$y <- newdata[[feat]]
      pred <- predict(object, newdata = newdata, type = "link", terms = mterms)
      newdata <- cbind(newdata, pred)

      # apply fit2 estimates (shifts to the original parameters) to newdata
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0  # zero out any parameters not requested for adjustment
      # Combine frozen offsets with estimated shifts on the link scale, then
      # map2par() inverts the links back to natural parameter scale.
      # NOTE (potential bug): fit2's linear predictor already includes the
      # offset() terms, so predict(fit2, type="link") returns offset + shift.
      # Adding `pred` again here may double-count the offset. Flagging only.
      params <- family(fit2)$map2par(pred + shift)
    } else {
      # No adjustment
      params <- predict(object, newdata = newdata, type = "parameter", terms = mterms)
    }

    # Compute the requested score using the family's own distribution functions.
    switch(
      type,
      "cent" = object$family$cdf(newdata[,feat], par = params),        # CDF at y -> centile in [0,1]
      "resid" = newdata[,feat] - params[,1],                           # raw residual y - mu
      "zscore" = (newdata[,feat] - params[,1])/params[,2],             # (y - mu)/sigma
      "quantile" = object$family$rqres(newdata[,feat], par = params),  # (randomized) quantile residuals on N(0,1) scale
      "parameter" = params                                             # the fitted parameters themselves
    )
  }

#' @rdname predict_score
#' @export
predict_score.gamlss <-
  function(object, newdata, refdata = NULL,
           type = c("cent", "resid", "zscore", "quantile", "parameter"),
           adjust = TRUE, rm.term = NULL,
           newformula = NULL,
           which.params = NULL,
           traindata = NULL) {
    type = match.arg(type)

    if (!is.null(rm.term) && !adjust) {
      warning("`rm.term` is set with `adjust = FALSE`: the '", rm.term, "' effect ",
              "is removed but not re-estimated/batch-adjusted. For a random ",
              "effect this sets it to the population mean (0); for a fixed factor ",
              "it sets it to the reference level. Use `adjust = TRUE` to estimate ",
              "and apply a batch adjustment.")
    }

    if(object$family[1] != "NO" & type == "zscore") {
      warning("Z-scores may not be valid for families other than NO")
    }

    feat <- as.character(object$formula[[2]])

    # Resolve the adjustment spec: honor explicit args, else auto-derive from the
    # moments containing rm.term (adjust those, freeze the rest).
    spec <- .resolve_adjust_spec(newformula, which.params,
                                 .moment_formulas_gamlss(object), rm.term)
    newformula   <- spec$newformula
    which.params <- spec$which.params
    which.params <- setNames(1:4, c("mu", "sigma", "nu", "tau"))[which.params]
    if (is.null(refdata)) {
      refdata <- newdata
    } else {
      #confirm all newdata levels are in refdata
      stopifnot("Not all levels of rm.term in newdata are present in refdata " = levels(newdata[[rm.term]]) %in% levels(refdata[[rm.term]]))
    }

    # Informative error for unseen levels of a NON-rm.term factor
    known_levels <- list()
    for (p in object$parameters) {
      xl <- object[[paste0(p, ".xlevels")]]
      for (fv in setdiff(names(xl), rm.term))
        known_levels[[fv]] <- union(known_levels[[fv]], xl[[fv]])
    }
    for (fv in names(known_levels)) {
      for (nm in c("newdata", "refdata")) {
        d <- if (nm == "newdata") newdata else refdata
        if (fv %in% names(d)) {
          unseen <- setdiff(na.omit(as.character(unique(d[[fv]]))), known_levels[[fv]])
          if (length(unseen))
            stop("factor `", fv, "` in ", nm, " has new level", paste(unseen, collapse = ", "),
                 ". Only `rm.term` may introduce unseen levels.")
        }
      }
    }

    ###NOTE: may make more sense to add option to override use_datafree by providing traindata -
    # basically what im implementing now in gamlssTools, so you can compare datafree vs traindata outputs
    use_datafree <- .datafree_eligible_gamlss(object, rm.term)
      if (!use_datafree && is.null(traindata)) {
      stop("This model has a kept smoother type (cs/ps/ga/s) that is not yet ",
           "reconstructed data-free. Supply `traindata` (the original fitting ",
           "data) to use the predict-based path.")
    }

    # Find rm.term in any moment and append unseen factor levels.
    # Not needed for the data-free path (it sets a fixed-factor rm.term to its
    # reference level directly), and appending levels there would misalign the
    # design with coef(); so only run this for the predict-based path.
    # A random effect (random(site)) never appears in *.xlevels, so the loop
    # simply matches nothing and we proceed.
    if (!is.null(rm.term) && !use_datafree) {
      for (p in object$parameters) {
        slot <- paste0(p, ".xlevels")
        if (rm.term %in% names(object[[slot]])) {
          oldlevels <- object[[slot]][[rm.term]]
          newlevels <- setdiff(levels(refdata[[rm.term]]), oldlevels)
          object[[slot]][[rm.term]] <- c(oldlevels, newlevels)
        }
      }
    }

    if (adjust) {
      # `rm.term` must be a single batch level shared by newdata and refdata for adjustment
      if (!is.null(rm.term)) {
        new_lv <- unique(as.character(newdata[[rm.term]]))
        ref_lv <- unique(as.character(refdata[[rm.term]]))
        if (length(new_lv) != 1L || length(ref_lv) != 1L || !identical(new_lv, ref_lv))
          stop("with adjust = TRUE, `rm.term` ('", rm.term, "') must be a single ",
               "batch level shared by newdata and refdata (the adjustment estimates ",
               "one shift for that batch). newdata level(s): ",
               paste(new_lv, collapse = ", "), "; refdata level(s): ",
               paste(ref_lv, collapse = ", "), ".")
      }
      # get offsets for refdata and fit adjustment model
      off_ref <- .pop_offset_gamlss(object, refdata, rm.term, object$parameters,
                                    traindata, datafree = use_datafree)
      refdata <- cbind(refdata, off_ref)   # add mu/sigma... columns for the offsets
      refdata$y <- refdata[[feat]]         # `newformula` LHS is `y`
      fit2 <- gamlss2::gamlss2(newformula, data = refdata, family = object$family[1])

      # get offsets for newdata (same computation on the observations to score)
      off_new <- .pop_offset_gamlss(object, newdata, rm.term, object$parameters,
                                    traindata, datafree = use_datafree)
      newdata <- cbind(newdata, off_new)

      # apply fit2 estimates, which are shifts to the original parameters
      shift <- predict(fit2, newdata = newdata, type = "link")
      shift[,-which.params] <- 0
      params <- family(fit2)$map2par(off_new + shift)
    } else {
      #no adjustment
      pnames <- c("mu", "sigma", "nu", "tau")[which.params]
      if (use_datafree) {
        links <- setNames(lapply(pnames, function(pp)
          .lp_nodata_gamlss(object, pp, newdata, drop.term = rm.term)), pnames)
      } else {
        mterms <- c("Intercept", setdiff(all.vars(formula(object)), c(feat, rm.term)))
        links <- predictAll(object, newdata = newdata, type = "link", terms = mterms)[which.params]
        names(links) <- pnames
      }
      params <- lapply(setNames(nm = pnames), function(x) {
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
