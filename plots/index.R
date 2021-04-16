#' ---
#' title: Improved plots for intactness and sector effects
#' author: Peter, Dave, Emily, and others
#' output: word_document
#' ---
#'
#' The goal of this script is to provide functions to create intactness and sector effect graphics replacing the violin plots used previously.
#'
#' ## Install and load dependencies
#'
#' We will need to beesplot package to be able to make the plot using glorious base R graphics.
#'
#'+ results='hide',warning=FALSE,message=FALSE
if (!require("beeswarm", quietly=TRUE))
  install.packages("beeswarm")
library(beeswarm)
#'
#' ## Simulated data
#'
#' - `group`: groups of species or sectors to be shown on the x axis
#' - `value`: values to show in the ye axis
set.seed(1)
n <- c(200, 400)
x <- data.frame(
  group=rep(c("Group 1", "Group 2"), n),
  value=c(rnorm(n[1], 0.1, 0.1),
          rnorm(n[2]/2, -0.2, 0.025), rnorm(n[2]/2, -0.05, 0.05)))
x$value <- plogis(x$value) * 100
#'
#' ## Univariate plots
#'
#' Unimodal for `Group 1`
hist(x$value[x$group=="Group 1"])
#' Bimodal for `Group 2`
hist(x$value[x$group=="Group 2"])
#'
#' ## Conditional plot
#'
#' We want to show the two conditional distribution (histograms from above) in the same graph. One approach is boxplot
boxplot(value ~ group, x, range=0)
#' This does not show the bimodality, can we show to dots?
stripchart(value ~ group, x, vertical = TRUE, method = 'jitter', pch=21)
#' We are getting somewhere, the plot is really busy, let's improve
stripchart(value ~ group, x, vertical = TRUE, method = 'jitter',
           pch=19, col="#00000044")
#' Wouldn't it be lovely of the dots would follow the histogram shapes?
#' We can use beeswarm plot for that
beeswarm(value ~ group, x,
         priority="random", corral="random")
#' Now combine beeswarm with boxplot and color the groups
beeswarm(value ~ group, x,
         priority="random", corral="random",
        corralWidth=0.25,
         cex=0.5, pch=19, col=1:2)
b <- boxplot(value ~ group, x, range=0, add=TRUE,
        col=c("#00000044", "#ff000044"), border=NA,
        pars = list(boxwex = 0.5))
boxplot(b$stats[3,] ~ as.factor(b$names), add=TRUE,
        border=c("#00000044", "#ff000044"),
        pars = list(boxwex = 0.5))
#'
#' ## Wrapping it up in a function
#'
beesbox <- function(x, y,
  bw=0.4, sw=0.3,
  xlab="", ylab="Value", main="", col=NULL, ...) {

  x <- data.frame(group=as.factor(x), value=y)
  m <- length(unique(x$group))
  if (is.null(col))
    col <- hcl.colors(m, "Pastel 1")
  col0 <- substr(col, 1, 7)
  col1 <- paste0(col0, "44")
  col2 <- paste0(col0, "88")
  b <- boxplot(value ~ group, x, range=0,
          col=col1, border=NA,
          pars=list(boxwex = 2*bw),
          xlab=xlab, ylab=ylab, main=main, axes=FALSE, ...)
  a <- beeswarm(value ~ group, x, add=TRUE,
           priority="random", corral="random",
           cex=0.5, pch=19, col=1:2, do.plot=FALSE)
  g <- as.factor(b$names)
  for (i in seq_along(g)) {
    lines(c(i-bw, i+bw), b$stats[3,c(i,i)], col=col2[i], lwd=2, lend=1)
    a$x[a$x.orig == g[i]] <- i+sw*(a$x[a$x.orig == g[i]]-i) /
      max(abs(a$x[a$x.orig == g[i]]-i))
  }
  a$col <- col0[match(a$x.orig, g)]
  points(a$x, a$y, cex=0.5, pch=19, col=a$col)
  axis(2)
  axis(1, seq_along(g), g, lwd=0)
  invisible(a)
}
#' Test of the pudding
beesbox(x$group, x$value, ylim=c(40, 60), ylab="Intactness",
        col=hcl.colors(length(unique(x$group)), "Dark 2"))
abline(h=50)

