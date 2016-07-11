## Frequentist inference, Bernoulli model

mle <-
function(p = 0.3, n = 10, seed = 0)
{
    if (p < 0 || p > 1)
        stop("p must be in 0-1")
    if (n < 1 || n > 1000)
        stop("n must be on 1-1000")
    set.seed(seed)
    y <- rbinom(n = 1000, size = 1, p = p)
    pt <- seq(0, 1, by = 0.0005)
    L <- sapply(pt, function(z)
        prod(dbinom(y[1:n], size = 1, prob = z)))
    op <- par(las = 1)
    plot(pt, L, type = "l", col="#3498db",
        ylab = "Likelihood", xlab="p",
        sub=paste0("Mean = ", round(mean(y[1:n]), 2), " (",
            sum(1-y[1:n]), " 0s & ", sum(y[1:n]), " 1s)"),
        main = paste("Estimate =", round(pt[which.max(L)], 2)))
    abline(v = p, lwd = 2, col = "#c7254e")
    abline(v = pt[which.max(L)], lwd = 2, col = "#18bc9c")
    par(op)
    invisible(list(y = y, pt = pt, L = L))
}

