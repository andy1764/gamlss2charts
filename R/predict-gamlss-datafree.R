# Set of internal helper functions that allow the predict_score.gamlss method to
# get predictions without the original model-fitting data

# ---- internal: is a gamlss fit eligible for data-free prediction? ------------
# rebuilds parametric terms (coef + design), pb() smooths (stored knots + coefficients)
# and random() effects (stored per-level BLUPs). Other gamlss smoother types (
# cs, ps, ga/s, ...) are not currently implemented. TRUE = the model has NO kept
# (non-dropped) smoother of an unsupported type
.datafree_eligible_gamlss <- function(object, rm.term) {
  ok <- TRUE
  for (p in object$parameters) {
    sm <- colnames(object[[paste0(p, ".s")]])
    for (lab in sm) {
      supported <- grepl("^pb\\(", lab) || grepl("^random\\(", lab)
      dropped   <- !is.null(rm.term) && rm.term %in% all.vars(str2lang(lab))
      if (!supported && !dropped) ok <- FALSE
    }
  }
  ok
}

# ---- internal: data-free link-scale linear predictor for one parameter -------
# Rebuilds parameter `p`'s link-scale predictor on `newdata` WITHOUT the original
# fitting data, dropping `drop.term`. The parametric part is aligned to coef() BY NAME
# (coef also carries an entry per smoother, so positional alignment is unsafe).
# Each kept pb() term adds its linear coefficient * x plus the stored
# interpolation function getSmo(...)$fun(x); each kept random() effect adds the
# stored per-level BLUP getSmo(...)$coef[level] (unseen levels -> 0, the
# population value). Only valid when model is data-free eligible (see
# .datafree_eligible_gamlss())
.lp_nodata_gamlss <- function(object, p, newdata, drop.term = NULL) {
  cf   <- coef(object, p)
  tl   <- attr(terms(object[[paste0(p, ".formula")]]), "term.labels")
  sm   <- colnames(object[[paste0(p, ".s")]]); if (is.null(sm)) sm <- character(0)
  pb_lab     <- sm[grepl("^pb\\(", sm)]
  random_lab <- sm[grepl("^random\\(", sm)]
  param_lab  <- setdiff(tl, sm)                      # genuine parametric terms
  xlev <- object[[paste0(p, ".xlevels")]]

  # factor covariates -> training levels
  for (fv in names(xlev)) {
    if (!is.null(drop.term) && fv == drop.term) {
      # drop.term -> constant
      newdata[[fv]] <- factor(xlev[[fv]][1], levels = xlev[[fv]],
                              ordered = is.ordered(newdata[[fv]]))
    } else if (fv %in% names(newdata)) {
      # align levels with training data (factor() keeps an already-ordered class)
      newdata[[fv]] <- factor(newdata[[fv]], levels = xlev[[fv]],
                              ordered = is.ordered(newdata[[fv]]))
    }
  }

  # parametric part (intercept + genuine parametric terms), aligned by name
  pfo <- if (length(param_lab)) stats::reformulate(param_lab) else ~1
  mf  <- stats::model.frame(pfo, newdata, na.action = stats::na.pass)
  Xp  <- stats::model.matrix(pfo, mf)
  lp  <- as.numeric(Xp %*% cf[colnames(Xp)])

  # pb() smooths: linear part (coef * x) + stored nonlinear interpolation
  for (lab in pb_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # pb on the batch var -> dropped
    v  <- vars[1]
    lp <- lp + cf[[lab]] * newdata[[v]] +
      gamlss::getSmo(object, p, which = match(lab, sm))$fun(newdata[[v]])
  }

  # random() effects: add the stored per-level BLUP (unseen levels -> population 0)
  for (lab in random_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # dropped batch random effect
    v    <- vars[1]
    blup <- gamlss::getSmo(object, p, which = match(lab, sm))$coef
    b    <- as.numeric(blup[as.character(newdata[[v]])])
    if (anyNA(b)) {
      warning("random(", v, "): ", sum(is.na(b)),
              " level(s) not seen in the fit; their effect is set to 0 (population).")
      b[is.na(b)] <- 0
    }
    lp <- lp + b
  }
  lp
}

# ---- internal: batch-baseline offsets for the gamlss method ------------------
# Per-parameter link-scale predictions with the batch term set to its baseline
# (population mean 0 for a random effect, reference level for a fixed factor). When
# `datafree = TRUE` (default whenever every kept smoother is pb()/random(),
# including purely parametric models) the reconstruction uses stored coefficients
# + pb interpolation + random BLUPs and needs no original data; otherwise it goes
# through predict.gamlss (type = "terms"/"link"), which requires the training
# data in scope or via `traindata`.
.pop_offset_gamlss <- function(object, scoredata, rm.term, params,
                               traindata = NULL, datafree = FALSE) {
  # Keep only the variables the model actually uses
  model_vars <- unique(unlist(lapply(params, function(p)
    all.vars(object[[paste0(p, ".formula")]]))))
  scoredata <- scoredata[, intersect(model_vars, names(scoredata)), drop = FALSE]

  #loop over moments
  out <- lapply(params, function(p) {
    #if non training data required, use alternate function
    if (datafree) return(.lp_nodata_gamlss(object, p, scoredata, drop.term = rm.term))
    fo <- object[[paste0(p, ".formula")]]
    drop_term <- !is.null(rm.term) && !is.null(fo) && rm.term %in% all.vars(fo)
    #eval if rm.term is present
    if (drop_term) {
      #use training data to predict without rm.term
      args <- list(object, what = p, newdata = scoredata, type = "terms")
      if (!is.null(traindata)) args$data <- traindata
      tm <- do.call(predict, args)
      ic <- attr(tm, "constant"); if (is.null(ic)) ic <- 0
      # drop only columns whose variables include rm.term (matches "random(site)"
      # and "site", but not look-alikes such as "prestige_site")
      drop <- vapply(colnames(tm), function(cn) {
        v <- tryCatch(all.vars(stats::reformulate(cn)), error = function(e) character(0))
        rm.term %in% v
      }, logical(1))
      ic + rowSums(tm[, !drop, drop = FALSE])
    } else {
      #use training data to predict
      args <- list(object, what = p, newdata = scoredata, type = "link")
      if (!is.null(traindata)) args$data <- traindata
      as.numeric(do.call(predict, args))
    }
  })
  as.data.frame(setNames(out, params))
}
