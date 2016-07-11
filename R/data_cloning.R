## Bernoulli model, data cloning

data_cloning <-
function(p = 0.3, n = 10, a = 1, b = 1,
K = 1, scale = c("prob", "logit"), seed = 0)
{
    scale <- match.arg(scale)
    if (p < 0 || p > 1)
        stop("p must be in 0-1")
    if (n < 1 || n > 1000)
        stop("n must be on 1-1000")
    set.seed(seed)
    y <- rbinom(n = n, size = 1, p = p)
    yk <- rep(y, K)
    BY <- 0.0005
    pval <- seq(0.001, 0.999, by = BY)
    fLik <- function(p, y)
        sum(dbinom(y, size = 1, prob = p, log=TRUE))
    Lik <- exp(sapply(pval, fLik, y=yk))
    if (all(Lik <= 0)) {
        est <- optimize(fLik, c(0.001, 0.999), y=yk, maximum=TRUE)$maximum
        Lik[which.min(abs(pval - est))] <- 1
    }
    if (scale == "prob") {
        p <- p
        fPri <- function(p, shape1=0.5, shape2=0.5)
            dbeta(p, shape1, shape2)
        Pri <- sapply(pval, fPri, a, b)
    } else {
        p <- qlogis(p)
        N <- 10^5
        x <- rbeta(N, a, b)
        br <- c(0.001, seq(0.001+BY/2, 0.999-BY/2, by = BY), 0.999)
        d <- as.numeric(table(cut(x, breaks=br))) / N
        pval <- qlogis(pval)
        g <- diff(qlogis(br))
        gy <-  d / g
        Pri <- smooth.spline(pval, gy)$y
    }
    Pos <- Lik * Pri
    M <- cbind(Pri=Pri/max(Pri),
        Lik=Lik/max(Lik),
        Pos=Pos/max(Pos))
    Col <- c("#cccccc", "#3498db", "#f39c12")
    op <- par(las = 1)
    matplot(pval, M, type = "l",
        col=Col, lwd=2, lty=1,
        ylab = "Density", xlab="p",
        sub=paste0("Mean = ", round(mean(y[1:n]), 2), " (",
            sum(1-y[1:n]), " 0s & ", sum(y[1:n]), " 1s)"),
        main = paste0("True value = ", round(p, 2),
            ", Posterior mode = ", round(pval[which.max(Pos)], 2)))
    abline(v = p, lwd = 2, col = "#c7254e")
    abline(v = pval[which.max(Pos)], lwd = 2, col = "#18bc9c")
    legend("topleft",lty=1, lwd=2, col=Col, bty="n",
           legend=c("Prior","Likelihood","Posterior"))
    par(op)
    if (scale == "prob") {
        M <- cbind(p=pval, M)
    } else {
        M <- cbind(logit_p=pval, M)
    }
    invisible(data.frame(M))
}
