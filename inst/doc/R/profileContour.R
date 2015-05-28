n <- 5;
J <- 8;
N <- n * J;
g <- rep(1:J, rep(n, J));
mu <- 0;
sigma.the <- 1;
sigma.eps <- 1;

set.seed(0);

theta <- rnorm(J, 0, sigma.the * sigma.eps);
y <- rnorm(N, mu + theta[g], sigma.eps);

y.bar <- mean(y);
y.bar.j <- sapply(1:J, function(j) mean(y[g == j]));
S.w <- sum((y - y.bar.j[g])^2);
S.b <- sum((y.bar.j - y.bar)^2);

profiledLikelihood <- function(sigma.the, sigma.eps) {
  sigma.the.sq <- sigma.the^2;
  sigma.eps.sq <- sigma.eps^2;
  v <- sigma.the.sq + 1 / n;
  
  -0.5 * N * log(sigma.eps.sq) - 0.5 * J * log(v) - 0.5 *
    (S.w + S.b / v) / sigma.eps.sq;
}

sigma.eps.hat <- function(sigma.the) {
  v <- sigma.the^2 + 1 / n;
  sqrt((S.w + S.b / v) / N);
}

profiledProfiledLikelihood <- function(sigma.the) {
  v <- sigma.the^2 + 1 / n;
  -(N / 2) * log(S.w + S.b / v) - (J / 2) * log(v);
}

sigma.thes <- seq(0.25, 2.5, length.out = 101);
sigma.epss <- seq(0.6, 1.25, length.out = 101);

zValues <- sapply(sigma.epss, function(y) profiledLikelihood(sigma.thes, y));
zValues <- exp(zValues - median(zValues));
zValues <- zValues / sum(zValues * (sigma.thes[2] - sigma.thes[1]) * (sigma.epss[2] - sigma.epss[1]));

pdf(imgPath, height = defaultImgHeight, width = defaultImgWidth * 2);
par(mfrow = c(1, 2));
par(defaultPars);
contour(sigma.thes, sigma.epss, zValues,
        xlab = expression(sigma[b]), ylab=expression(sigma),
        main = "Likelihood", drawlabels = FALSE);

yValues <- sigma.eps.hat(sigma.thes);
lines(sigma.thes, yValues, col="gray");
text(sigma.thes[2], yValues[2], expression(hat(sigma)(sigma[b])),
     adj = c(-0.1, 0.5), cex = defaultPars$cex.lab * defaultPars$cex);

yValues <- profiledProfiledLikelihood(sigma.thes);
yValues <- exp(yValues - median(yValues));
plotPars <- defaultPars;
plotPars$mar[2] <- 0.1;
par(plotPars);
plot(sigma.thes, yValues, type = 'l', xlab=expression(sigma[b]),
    yaxt='n', ylab = "", bty = 'n', main="Profiled Likelihood", yaxs = 'i');
ignored <- dev.off();
