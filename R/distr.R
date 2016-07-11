## Statistical distributions

distr <-
function(distr = "bernoulli", n = 100,
unif_a = -1, unif_b = 1,
binom_p = 0.3, binom_size = 10,
poisson_lambda = 5,
normal_mu = 0, normal_var = 1,
beta_shape1 = 1, beta_shape2 = 1,
gamma_shape = 1, gamma_rate = 1,
seed = 0)
{
    distr <- match.arg(tolower(distr),
        c("bernoulli", "binomial", "poisson", "normal",
          "lognormal", "uniform", "beta", "gamma"))
    if (distr == "uniform" && unif_b < unif_a)
        stop("unif_b must be greater than unif_a")
    if (binom_p < 0 || binom_p > 1)
        stop("p must be in 0-1")
    if (n < 1 || n > 1000)
        stop("n must be on 1-1000")
    op <- par(las = 1)
    set.seed(seed)
    y <- switch(distr,
        "bernoulli" = rbinom(1000, 1, binom_p),
        "binomial" = rbinom(1000, binom_size, binom_p),
        "poisson" = rpois(1000, poisson_lambda),
        "normal" = rnorm(1000, normal_mu, sqrt(normal_var)),
        "lognormal" = rlnorm(1000, normal_mu, sqrt(normal_var)),
        "uniform" = runif(1000, unif_a, unif_b),
        "beta" = rbeta(1000, beta_shape1, beta_shape2),
        "gamma" = rgamma(1000, gamma_shape, gamma_rate))
    yy <- y[1:n]
    x <- switch(distr,
        "bernoulli" = c(0,1),
        "binomial" = seq(0, max(yy)+1, by = 1),
        "poisson" = seq(0, max(yy)+1, by = 1),
        "normal" = seq(min(yy)-1, max(yy)+1, length.out = 1000),
        "lognormal" = seq(0.0001, max(yy)+1, length.out = 1000),
        "uniform" = seq(unif_a+0.0001, unif_b-0.0001, length.out = 1000),
        "beta" = seq(0.0001, 0.9999, length.out = 1000),
        "gamma" = seq(0.0001, max(yy), length.out = 1000))
    d <- switch(distr,
        "bernoulli" = dbinom(x, 1, binom_p),
        "binomial" = dbinom(x, binom_size, binom_p),
        "poisson" = dpois(x, poisson_lambda),
        "normal" = dnorm(x, normal_mu, sqrt(normal_var)),
        "lognormal" = dlnorm(x, normal_mu, sqrt(normal_var)),
        "uniform" = dunif(x, unif_a, unif_b),
        "beta" = dbeta(x, beta_shape1, beta_shape2),
        "gamma" = dgamma(x, gamma_shape, gamma_rate))
    xlab <- "x"
    ylab <- "Density"
    main <- paste0(toupper(substring(distr, 1, 1)), substring(distr, 2),
        " distribution (n = ", n, ")")
    if (distr %in% c("bernoulli", "binomial", "poisson")) {
        tmp <- table(yy) / n
        plot(tmp, ylim=c(0, max(tmp, d)),
            ylab = ylab, xlab = xlab, main = main,
            col = "#cccccc", lwd = 10)
        points(x, d, pch = 21, col = "#c7254e", type = "b",
            lty = 2, cex = 2)
    } else {
        tmp <- hist(yy, plot = FALSE)
        hist(yy, freq = FALSE, ylim=c(0, max(tmp$density, d)),
            ylab = ylab, xlab = xlab, main = main,
            col = "#ecf0f1", border = "#cccccc")
        lines(x, d, lwd = 2, col = "#c7254e")
    }
    par(op)
    invisible(list(y = y, x = x, d = d))
}
