# check_unseen_site.R
# -----------------------------------------------------------------------------
# Validation of the out-of-sample site adjustment in predict_score().
#
# Idea: a group of subjects belongs to "siteNEW". We compute their centiles two
# ways and check they agree:
#
#   (A) GOLD STANDARD - siteNEW is part of the ORIGINAL fitted data, so the model
#       estimates a site term for it directly. Centiles come straight from that
#       full model.
#
#   (B) UNSEEN SITE - the original model is fit WITHOUT siteNEW. siteNEW is then
#       treated as a new batch: predict_score(adjust = TRUE, rm.term = "site")
#       estimates the site shift on siteNEW reference data and applies it.
#
# If the adjustment machinery is correct, centiles from (A) and (B) should match
# up to sampling noise in the shared coefficients. A large, systematic gap would
# point at a bug (e.g. the offset being double-counted).
#
# Every condition calls record() to log the A-vs-B agreement (max |diff|,
# mean |diff|, correlation); the accumulated results table is printed at the end.
# -----------------------------------------------------------------------------

suppressMessages(library(gamlss2))
suppressMessages(library(gamlss))     # for the gamlss method + pNO()
source("~/Desktop/BGD_Repos/gamlss2charts/R/gamlss2charts.R") #source edited versions
if (!requireNamespace("gamlssTools", quietly = TRUE))
  devtools::install_github("BGDlab/gamlssTools", build_vignettes = FALSE)
library(gamlssTools)

## ---- results table ---------------------------------------------------------
## record() logs the agreement between gold-standard (A) and out-of-sample (B)
## centiles for one condition/engine. Printed as a table at the very end.
results <- data.frame(condition = character(), engine = character(),
                      max_abs = numeric(), mean_abs = numeric(),
                      cor = numeric(), n = integer())
record <- function(condition, engine, centA, centB) {
  results <<- rbind(results, data.frame(
    condition = condition, engine = engine,
    max_abs   = max(abs(centA - centB)),
    mean_abs  = mean(abs(centA - centB)),
    cor       = cor(centA, centB),
    n         = length(centA)))
}

## ---- simulate toy data ------------------------------------------------------
## Response depends on x with a common slope/scale across sites; each site only
## shifts the mean (a pure mu intercept shift = the "batch effect").
set.seed(1)
sites      <- c("site1", "site2", "site3", "siteNEW")
n_per      <- 500
site_shift <- c(site1 = 0, site2 = 0.8, site3 = -0.6, siteNEW = 1.5)

dat <- do.call(rbind, lapply(sites, function(s) {
  x  <- runif(n_per, 0, 10)
  mu <- 2 + 0.5 * x + site_shift[[s]]
  data.frame(x = x, y = rnorm(n_per, mu, sd = 1), site = s)
}))
dat$site <- factor(dat$site)

new_dat <- droplevels(subset(dat, site == "siteNEW"))   # subjects to score
train_dat <- droplevels(subset(dat, site != "siteNEW"))   # only the "known" sites

#check
unique(new_dat$site)
unique(train_dat$site)
unique(dat$site)

## ============================================================================
### BASIC EXAMPLE
## ============================================================================

## gamlss2
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA2  <- gamlss2(y ~ x + site | 1, family = NO, data = dat)
parA2  <- predict(fitA2, newdata = new_dat, type = "parameter")
centA2 <- fitA2$family$cdf(q = new_dat$y, par = parA2)
centA2b <- pred_og_centile(fitA2, new_dat) #gamlssTools method - confirm identical
#max diff
max(abs(centA2 - centA2b))
#cor
cor(centA2, centA2b)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB2  <- gamlss2(y ~ x + site | 1, family = NO, data = train_dat)
centB2 <- predict_score(fitB2, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Basic", "gamlss2", centA2, centB2)

cat("\n===================== gamlss2 =====================\n")
cat("head(A) :", round(head(centA2), 4), "\n")
cat("head(B) :", round(head(centB2), 4), "\n")
cat("summary(A - B):\n"); print(summary(centA2 - centB2))
cat(sprintf("max|A-B| = %.4g   cor = %.5f\n",
            max(abs(centA2 - centB2)), cor(centA2, centB2)))

plot(centA2, centB2, main = "gamlss2", xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = dat,
                 trace = FALSE)
pA1    <- predictAll(fitA1, newdata = new_dat)
centA1 <- pNO(new_dat$y, mu = pA1$mu, sigma = pA1$sigma)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Basic", "gamlss", centA1, centB1)

cat("\n===================== gamlss ======================\n")
cat("head(A) :", round(head(centA1), 4), "\n")
cat("head(B) :", round(head(centB1), 4), "\n")
cat("summary(A - B):\n"); print(summary(centA1 - centB1))
cat(sprintf("max|A-B| = %.4g   cor = %.5f\n",
            max(abs(centA1 - centB1)), cor(centA1, centB1)))
plot(centA1, centB1, main = "gamlss",  xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

#test adjust=FALSE
centB3 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust=FALSE, rm.term = "site", which.params = "mu")
max(abs(centA1 - centB3)) #does poorly as expected

## ============================================================================
### BATCH EFFECTS IN MU & SIGMA
## ============================================================================
## NOTE: site enters both mu and sigma. which.params = c("mu","sigma") re-estimates
## the siteNEW shift for BOTH parameters (offset(mu)|offset(sigma) supply the
## population baselines), so agreement should stay tight.

## gamlss2
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA2  <- gamlss2(y ~ x + site | x + site, family = NO, data = dat)
parA2  <- predict(fitA2, newdata = new_dat, type = "parameter")
centA2 <- fitA2$family$cdf(q = new_dat$y, par = parA2)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB2  <- gamlss2(y ~ x + site | x + site, family = NO, data = train_dat)
centB2 <- predict_score(fitB2, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu) | offset(sigma), which.params = c("mu", "sigma"))
record("Batch mu&sigma", "gamlss2", centA2, centB2)

#max diff
max(abs(centA2 - centB2)) #slightly less precise
#cor
cor(centA2, centB2)

plot(centA2, centB2, main = "gamlss2", xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ x + site, sigma.formula = ~ x + site, family = NO, data = dat,
                 trace = FALSE)
pA1    <- predictAll(fitA1, newdata = new_dat)
centA1 <- pNO(new_dat$y, mu = pA1$mu, sigma = pA1$sigma)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ x + site, sigma.formula = ~ x + site, family = NO, data = train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu) | offset(sigma), which.params = c("mu", "sigma"))
record("Batch mu&sigma", "gamlss", centA1, centB1)

#max diff
max(abs(centA1 - centB1))
#cor
cor(centA1, centB1)

plot(centA1, centB1, main = "gamlss",  xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## ============================================================================
### SMOOTHS
## ============================================================================

## gamlss2
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA2  <- gamlss2(y ~ s(x) + site | 1, family = NO, data = dat)
parA2  <- predict(fitA2, newdata = new_dat, type = "parameter")
centA2 <- fitA2$family$cdf(q = new_dat$y, par = parA2)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB2  <- gamlss2(y ~ s(x) + site | 1, family = NO, data = train_dat)
centB2 <- predict_score(fitB2, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Smooths", "gamlss2", centA2, centB2)

#max diff
max(abs(centA2 - centB2))
#cor
cor(centA2, centB2)

plot(centA2, centB2, main = "gamlss2", xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ pb(x) + site, sigma.formula = ~1, family = NO, data = dat,
                 trace = FALSE)
pA1    <- predictAll(fitA1, newdata = new_dat)
centA1 <- pNO(new_dat$y, mu = pA1$mu, sigma = pA1$sigma)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ pb(x) + site, sigma.formula = ~1, family = NO, data = train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Smooths", "gamlss", centA1, centB1)

#max diff
max(abs(centA1 - centB1))
#cor
cor(centA1, centB1)

plot(centA1, centB1, main = "gamlss",  xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)


## ============================================================================
### RANDOM BATCH EFFECTS
## ============================================================================

## gamlss2 - NOT SUPPORTED. predict_score.gamlss2 has no type="terms"
## population-offset path, so s(site, bs = "re") on an unseen level triggers a
## "contrasts ... 2 or more levels" error. Logged as NA below.
fitA2  <- gamlss2(y ~ x + s(site, bs = "re") | 1, family = NO, data = dat)
gamlss2_re <- tryCatch({
  fitB2 <- gamlss2(y ~ x + s(site, bs = "re") | 1, family = NO, data = train_dat)
  predict_score(fitB2, newdata = new_dat, refdata = new_dat, type = "cent",
                adjust = TRUE, rm.term = "site",
                newformula = y ~ offset(mu), which.params = "mu")
}, error = function(e) e)
if (inherits(gamlss2_re, "error")) {
  cat("gamlss2 random effect NOT supported:", conditionMessage(gamlss2_re), "\n")
  results <- rbind(results, data.frame(condition = "Random effect", engine = "gamlss2",
                                       max_abs = NA, mean_abs = NA, cor = NA, n = NA))
}

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ x + random(site), sigma.formula = ~1, family = NO, data = dat,
                 trace = FALSE)
centA1 <- pred_og_centile(fitA1, dat, new.data=new_dat)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ x + random(site), sigma.formula = ~1, family = NO, data = train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Random effect", "gamlss", centA1, centB1)

#max diff
max(abs(centA1 - centB1))
#cor
cor(centA1, centB1)

plot(centA1, centB1, main = "gamlss",  xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

#does passing traindata make a difference?
centB3 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu", traindata = train_dat)

max(abs(centA1 - centB3)) #not really

## ============================================================================
### CASES & CONTROLS USING REF DATA
## ============================================================================
## Patients (newdata) differ from the reference controls (refdata): fewer/more
## rows and added noise. The site adjustment is estimated on the controls and
## applied to score the patients.

set.seed(2)
n_pat <- 200
idx <- sample(nrow(new_dat), n_pat, replace = TRUE)   # replace=TRUE lets n_pat exceed nrow
pt_dat <- new_dat[idx, ]
pt_dat$y <- pt_dat$y + rnorm(n_pat, 0, 0.75) # measurement/patient noise
rownames(pt_dat) <- NULL

## gamlss2
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA2  <- gamlss2(y ~ x + site | 1, family = NO, data = dat)
parA2  <- predict(fitA2, newdata = pt_dat, type = "parameter")
centA2 <- fitA2$family$cdf(q = pt_dat$y, par = parA2)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB2  <- gamlss2(y ~ x + site | 1, family = NO, data = train_dat)
centB2 <- predict_score(fitB2, newdata = pt_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Cases & controls", "gamlss2", centA2, centB2)

#max diff
max(abs(centA2 - centB2))
#cor
cor(centA2, centB2)

plot(centA2, centB2, main = "gamlss2", xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ x + site, sigma.formula =  ~1, family = NO, data = dat,
                 trace = FALSE)
pA1    <- predictAll(fitA1, newdata = pt_dat)
centA1 <- pNO(pt_dat$y, mu = pA1$mu, sigma = pA1$sigma)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = pt_dat, refdata = new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Cases & controls", "gamlss", centA1, centB1)

#max diff
max(abs(centA1 - centB1))
#cor
cor(centA1, centB1)

plot(centA1, centB1, main = "gamlss",  xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)


## ============================================================================
### DIFFERENT FAMILY
## ============================================================================

#simulate GG data
set.seed(1)
sites      <- c("site1", "site2", "site3", "siteNEW")
n_per <- 500
site_shift <- c(site1=0, site2=0.3, site3=-0.2, siteNEW=0.5)   # log-mu scale now
gg_dat <- do.call(rbind, lapply(sites, function(s){
  x  <- runif(n_per, 0, 10)
  mu <- exp(0.5 + 0.1*x + site_shift[[s]])        # exp() -> mu > 0
  data.frame(x = x, y = rGG(n_per, mu = mu, sigma = 0.3, nu = 2), site = s)
}))
gg_dat$site <- factor(gg_dat$site)
gg_new_dat <- droplevels(subset(gg_dat, site == "siteNEW"))   # subjects to score
gg_train_dat <- droplevels(subset(gg_dat, site != "siteNEW"))   # only the "known" sites

## gamlss2
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA2  <- gamlss2(y ~ x + site, family = GG, data = gg_dat)
centA2 <- pred_og_centile(fitA2, gg_dat, new.data = gg_new_dat)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB2  <- gamlss2(y ~ x + site | 1, family = GG, data = gg_train_dat)
centB2 <- predict_score(fitB2, newdata = gg_new_dat, refdata = gg_new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Diff family (GG)", "gamlss2", centA2, centB2)

#max diff
max(abs(centA2 - centB2))
#cor
cor(centA2, centB2)

plot(centA2, centB2, main = "gamlss2", xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ x + site, family = GG, data = gg_dat,
                 trace = FALSE)
# og.data = full data (all sites), new.data = subjects to score
centA1 <- pred_og_centile(fitA1, gg_dat, new.data = gg_new_dat)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ x + site, sigma.formula = ~1, family = GG, data = gg_train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = gg_new_dat, refdata = gg_new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Diff family (GG)", "gamlss", centA1, centB1)

#max diff
max(abs(centA1 - centB1))
#cor
cor(centA1, centB1)

plot(centA1, centB1, main = "gamlss",  xlab = "A: siteNEW in model",
     ylab = "B: siteNEW unseen", pch = 16, col = "#00000055"); abline(0, 1, col = 2)

## ============================================================================
### SMALL SAMPLE
## ============================================================================
set.seed(1)
sites      <- c("site1", "site2", "site3", "siteNEW")
n_per      <- 50
site_shift <- c(site1 = 0, site2 = 0.8, site3 = -0.6, siteNEW = 1.5)

small_dat <- do.call(rbind, lapply(sites, function(s) {
  x  <- runif(n_per, 0, 10)
  mu <- 2 + 0.5 * x + site_shift[[s]]
  data.frame(x = x, y = rnorm(n_per, mu, sd = 1), site = s)
}))
small_dat$site <- factor(small_dat$site)

small_new_dat <- droplevels(subset(small_dat, site == "siteNEW"))   # subjects to score
small_train_dat <- droplevels(subset(small_dat, site != "siteNEW"))   # only the "known" sites

## gamlss2
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA2  <- gamlss2(y ~ x + site | 1, family = NO, data = small_dat)
parA2  <- predict(fitA2, newdata = small_new_dat, type = "parameter")
centA2 <- fitA2$family$cdf(q = small_new_dat$y, par = parA2)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB2  <- gamlss2(y ~ x + site | 1, family = NO, data = small_train_dat)
centB2 <- predict_score(fitB2, newdata = small_new_dat, refdata = small_new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Small sample", "gamlss2", centA2, centB2)

max(abs(centA2 - centB2)) #higher but not bad
cor(centA2, centB2)

## gamlss
## (A) original fit INCLUDES siteNEW -> gold-standard centiles
fitA1  <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = small_dat,
                 trace = FALSE)
pA1    <- predictAll(fitA1, newdata = small_new_dat)
centA1 <- pNO(small_new_dat$y, mu = pA1$mu, sigma = pA1$sigma)

## (B) original fit EXCLUDES siteNEW -> treat siteNEW as an unseen site
fitB1  <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = small_train_dat,
                 trace = FALSE)
centB1 <- predict_score(fitB1, newdata = small_new_dat, refdata = small_new_dat,
                        type = "cent", adjust = TRUE, rm.term = "site",
                        newformula = y ~ offset(mu), which.params = "mu")
record("Small sample", "gamlss", centA1, centB1)

max(abs(centA1 - centB1))
cor(centA1, centB1)

## ============================================================================
### REALISTIC CHALLENGE
## ============================================================================
## A BCTo model with smooths + a second covariate + a random site effect. We
## sweep the number of siteNEW subjects (the reference sample used to estimate
## the batch adjustment) to see how agreement degrades as that sample shrinks.
## The 6 training sites are generated ONCE and held fixed; only siteNEW varies.
set.seed(1)
train_sites <- c("site1", "site2", "site3", "site4", "site5", "site6")
train_n     <- c(site1 = 100, site2 = 400, site3 = 200, site4 = 50, site5 = 500, site6 = 300)
train_shift <- c(site1 = 0, site2 = 0.8, site3 = 1.6, site4 = 0.2, site5 = 0.5, site6 = 0.3)
siteNEW_shift <- 0.5

gen_site <- function(s, n, shift) {
  x <- runif(n, 0, 10); z <- runif(n, 0, 3)
  data.frame(x = x, y = rBCT(n, 2 + 0.5 * x + shift + z), z = z, site = s)
}

# fixed training data + the "known sites" fit reused for every siteNEW size
train_dat <- do.call(rbind, lapply(train_sites,
                                   function(s) gen_site(s, train_n[[s]], train_shift[[s]])))
train_dat$site <- factor(train_dat$site)
fitB1 <- gamlss(y ~ pb(x) + z + random(site), sigma.formula = ~pb(x) + z + random(site),
                family = BCTo, data = train_dat, trace = FALSE)

for (n_new in c(25, 50, 75, 100)) {
  new_raw  <- gen_site("siteNEW", n_new, siteNEW_shift)
  full_dat <- rbind(train_dat, new_raw)          # training + this siteNEW sample
  full_dat$site <- factor(full_dat$site)
  new_dat  <- droplevels(subset(full_dat, site == "siteNEW"))

  ## (A) gold standard: fit INCLUDES siteNEW
  fitA1  <- gamlss(y ~ pb(x) + z + random(site), sigma.formula = ~pb(x) + z + random(site),
                   family = BCTo, data = full_dat, trace = FALSE)
  centA1 <- pred_og_centile(fitA1, full_dat, new.data = new_dat)

  ## (B) unseen site: score siteNEW off the fixed training fit
  centB1 <- predict_score(fitB1, newdata = new_dat, refdata = new_dat,
                          type = "cent", adjust = TRUE, rm.term = "site",
                          newformula = y ~ offset(mu) | offset(sigma),
                          which.params = c("mu", "sigma"))
  record(sprintf("Realistic (BCTo) n=%d", n_new), "gamlss", centA1, centB1)
}

## ============================================================================
### SUMMARY TABLE
## ============================================================================
## A-vs-B agreement across every condition tested. High cor and small max/mean
## |diff| mean the unseen-site machinery reproduces the gold-standard centiles.

results_print <- results
results_print[, c("max_abs", "mean_abs", "cor")] <-
  round(results_print[, c("max_abs", "mean_abs", "cor")], 4)
cat("\n============================ SUMMARY ============================\n")
print(results_print, row.names = FALSE)
