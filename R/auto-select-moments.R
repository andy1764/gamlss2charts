
# Set of internal helper functions that auto-select the correct moments to have offsets
# calculated and adjusted

# ---- internal: auto-derive the adjustment spec from rm.term ------------------
# The moments whose linear predictor contains `rm.term`.
.params_with_term <- function(moment_formulas, rm.term) {
  if (is.null(rm.term)) return(character(0))
  keep <- vapply(moment_formulas, function(fo)
    !is.null(fo) && rm.term %in% all.vars(fo), logical(1))
  names(moment_formulas)[keep]
}

# Per-moment formulas for a gamlss2 fit: split the stored `Formula` by its `|`
# components, aligned to the family's parameter names. Parameters beyond the
# supplied components (modelled intercept-only) get NULL.
.moment_formulas_gamlss2 <- function(object) {
  mom <- object$family$names
  F   <- Formula::Formula(object$formula)
  n   <- length(attr(F, "rhs"))
  stats::setNames(lapply(seq_along(mom), function(i)
    if (i <= n) stats::formula(F, lhs = 0, rhs = i) else NULL), mom)
}

# Per-moment formulas for a gamlss fit: the stored <param>.formula slots.
.moment_formulas_gamlss <- function(object) {
  mom <- object$parameters
  stats::setNames(lapply(mom, function(p) object[[paste0(p, ".formula")]]), mom)
}

# Build the adjustment `newformula`. Each moment carrying `rm.term` gets an
# estimated batch shift via `offset(param)` (implicit intercept); every other
# moment is frozen at the original fit with `offset(param) - 1` (no free
# intercept)
.build_adjust_formula <- function(moments, with_term) {
  rhs <- ifelse(moments %in% with_term,
                sprintf("offset(%s)", moments),
                sprintf("offset(%s) - 1", moments))
  stats::as.formula(paste("y ~", paste(rhs, collapse = " | ")))
}

# Resolve `newformula` / `which.params`. If the caller supplied either, stay in
# manual mode for now (legacy default for other param). If both are NULL,
# auto-derive from the moments that contain `rm.term`: adjust exactly those and
# freeze the rest

# Function to create newformula from `which.params` and `rm.terms`
# If neither provided, defaults to adjusting all parameters. If `which.params`
# not provided, adjusts parameters that originally include `rm.term`. If
# `which.params` is provided, adjusts specified parameters.
.resolve_adjust_spec <- function(moment_formulas, which.params, rm.term) {
  if (is.null(which.params) & is.null(rm.term)) {
    # all parameters
    which_params <- names(moment_formulas)
  } else if (is.null(which.params)) {
    # parameters that originally included rm.term
    which_params <- .params_with_term(moment_formulas, rm.term)
    if (length(which_params) == 0L) {
      stop("`rm.term` ('", rm.term, "') was not found in any moment")
    }
  } else {
    which_params <- which.params
  }

  return(list(newformula = .build_adjust_formula(names(moment_formulas), which_params),
              which.params = which_params))
}
