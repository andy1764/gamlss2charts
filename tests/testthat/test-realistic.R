# A harder end-to-end case: BCTo family, smooth + a second covariate + a random
# site effect. Slow (fits a random-effect BCTo model), so skipped on CRAN.

test_that("Realistic: BCTo with pb(x) + z + random(site)", {
  skip_on_cran()
  skip_if_not_installed("gamlss")

  gen_site <- function(s, n, shift, seed) {
    set.seed(seed)
    x <- runif(n, 0, 10); z <- runif(n, 0, 3)
    data.frame(x = x, y = gamlss.dist::rBCTo(n, 2 + 0.5 * x + shift + z), z = z, site = s)
  }

  train_shift <- c(site1 = 0, site2 = 0.8, site3 = 1.6, site4 = 0.5)
  train <- do.call(rbind, Map(function(s, i) gen_site(s, 150, train_shift[[s]], i),
                              names(train_shift), seq_along(train_shift)))
  train$site <- factor(train$site)
  new <- gen_site("siteNEW", 100, 0.6, 99); new$site <- factor("siteNEW")
  full <- rbind(train, new); full$site <- factor(full$site)

  form  <- y ~ pb(x) + z + gamlss::random(site)
  fitA  <- gamlss(form, sigma.formula = ~ pb(x) + z + gamlss::random(site), family = BCTo,
                  data = full, trace = FALSE)
  pA    <- predictAll(fitA, newdata = new, data = full)
  centA <- gamlss.dist::pBCTo(new$y, mu = pA$mu, sigma = pA$sigma, nu = pA$nu, tau = pA$tau)

  fitB  <- gamlss(form, sigma.formula = ~ pb(x) + z + gamlss::random(site), family = BCTo,
                  data = train, trace = FALSE)
  b <- score_both(fitB, new, NULL, NULL)
  # looser tolerances: small siteNEW n + random-effect shrinkage, and the z-scale
  # is wider than [0,1]. Checks centiles AND z-scores (quantile residuals).
  expect_agree(centA, b$cent, cor_min = 0.98, max_diff = 0.12, mean_diff = 0.06,
               zA = qnorm(centA), zB = b$quant,
               z_cor_min = 0.98, z_max_diff = 0.8, z_mean_diff = 0.25)
})
