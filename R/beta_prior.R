## Bernoulli model, Beta prior

beta_prior <-
function(p = 0.3, n = 10, a = 1, b = 1, scale = c("prob", "logit"), seed = 0)
{
    scale <- match.arg(scale)
    if (p < 0 || p > 1)
        stop("p must be in 0-1")
    if (n < 1 || n > 1000)
        stop("n must be on 1-1000")
    set.seed(seed)
    y <- rbinom(n = 1000, size = 1, p = p)
    BY <- 0.0005
    pval <- seq(0.001, 0.999, by = BY)
    fLik <- function(p, y)
        prod(dbinom(y, size = 1, prob = p))
    Lik <- sapply(pval, fLik, y=y[1:n])
    fPri <- function(p, shape1=0.5, shape2=0.5)
        dbeta(p, shape1, shape2)
    Pri <- sapply(pval, fPri, a, b)
    if (scale == "prob") {
        p <- p
    } else {
        p <- qlogis(p)
        br <- c(0.001, seq(0.001+BY/2, 0.999-BY/2, by = BY), 0.999)
        dx <- diff(pval)
        dx <- c(dx[1], dx)
        d <- Pri * dx / diff(qlogis(br))
        Pri <- smooth.spline(pval, d)$y
        pval <- qlogis(pval)
    }
    Pos <- Lik * Pri
    M <- cbind(Pri=Pri/max(Pri),
        Lik=Lik/max(Lik),
        Pos=Pos/max(Pos))
    Col <- c("#cccccc", "#3498db", "#f39c12")
    op <- par(las = 1)
    matplot(pval, M, type = "l",
        col=Col, lwd=2, lty=1,
        ylab = "Density",
        xlab=ifelse(scale == "logit", "logit(p)","p"),
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
