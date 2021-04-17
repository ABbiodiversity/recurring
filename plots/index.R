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

if (!require("wordcloud", quietly=TRUE))
  install.packages("wordcloud")
library(wordcloud)

if (!require("intrval", quietly=TRUE))
  install.packages("intrval")
library(intrval)
#'
#' ## Simulated data
#'
#' - `group`: groups of species or sectors to be shown on the x axis
#' - `value`: values to show in the ye axis
set.seed(1)
n <- c(200, 400)
x <- data.frame(
  group=as.factor(rep(c("Group 1", "Group 2"), n)),
  value=c(rnorm(n[1], 0.1, 0.1),
          rnorm(n[2]/2, -0.2, 0.025), rnorm(n[2]/2, -0.05, 0.05)))
x$value <- plogis(x$value) * 100
rownames(x) <- paste("Point", 1:nrow(x))
head(x)
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
#' ## Use labels
#'
#' Specify which rows to highlight: we take the extremes by groups
lab <- c(1:3, nrow(x)-(2:0))
#' Not the best because labels overlap with the swarm.
#' Try workdcloud to make sure text does not overlap
a <- beeswarm(value ~ group, x,
         priority="random", corral="random",
        corralWidth=0.25,
         cex=0.5, pch=19, col=1:2)

l <- a[lab,]
l$g <- as.integer(x$group[lab])
l$x.rel <- abs(l$x - l$g) + l$g
wl <- as.data.frame(wordlayout(l$x.rel, l$y, rownames(l)))
wl$off <- 0.4
wl$x.orig <- a$x[lab]
wl$y.orig <- a$y[lab]
wl

coltxt=4
cextxt=0.8
text(wl$x+wl$off, wl$y, rownames(wl), cex=cextxt, col=coltxt)
segments(wl$x+wl$off-wl$width/2, wl$y, wl$x.orig, wl$y.orig,
         col=coltxt)
#'
#' ## Wrapping it up in a function
#'
#' - x: group
#' - y: values
#' - bw: width for the box
#' - bs: width for the swarm
#' - ba: alpha level for the box
#' - sa: alpha level for the points and the median line
#' - lab: row IDs
#' - collab, cexlab, offlab: color, cex, and x/y offset for labels
#' - others are graphical params
beesbox <- function(x, y,
  bw=0.4, sw=0.3,
  ba=0.33, sa=0.66,
  xlab="", ylab="Value", main="", col=NULL, ...) {

  x <- data.frame(group=as.factor(x), value=y)
  #x <- x[order(x$group, x$value),]
  m <- length(unique(x$group))
  if (is.null(col))
    col <- hcl.colors(m, "Pastel 1")
  a1 <- substr(hcl.colors(m, "Pastel 1", alpha=ba), 8, 9)
  a2 <- substr(hcl.colors(m, "Pastel 1", alpha=sa), 8, 9)
  col0 <- substr(col, 1, 7)
  col1 <- paste0(col0, a1)
  col2 <- paste0(col0, a2)
  b <- boxplot(value ~ group, x, range=0,
          col=col1, border=NA,
          pars=list(boxwex = 2*bw),
          xlab=xlab, ylab=ylab, main=main, axes=FALSE, ...)
  a <- beeswarm(value ~ group, x, add=TRUE,
           priority="random", corral="random",
           cex=0.5, pch=19, col=1:2, do.plot=FALSE)
  rownames(a) <- rownames(x)
  a$g <- as.integer(x$group)
  g <- as.factor(b$names)
  for (i in seq_along(g)) {
    lines(c(i-bw, i+bw), b$stats[3,c(i,i)], col=col2[i], lwd=2, lend=1)
    a$x[a$x.orig == g[i]] <- i+sw*(a$x[a$x.orig == g[i]]-i) /
      max(abs(a$x[a$x.orig == g[i]]-i))
  }
  a$col <- col2[match(a$x.orig, g)]
  points(a$x, a$y, cex=0.5, pch=19, col=a$col)
  axis(2)
  axis(1, seq_along(g), g, lwd=0)
  invisible(a)
}
addlabels <- function(a,
  lab=NULL, collab=4, cexlab=0.8, xofflab=0.4, yofflab=0, ...) {
  if (!is.null(lab)) {
    if (is.character(lab))
      lab <- which(rownames(a) %in% lab)
    l <- a[lab,]
    l$g <- as.integer(x$group[lab])
    l$x.rel <- abs(l$x - l$g) + l$g
    wl <- as.data.frame(wordlayout(l$x.rel, l$y, rownames(l)))
    wl$x.off <- xofflab
    wl$y.off <- yofflab
    wl$x.orig <- a$x[lab]
    wl$y.orig <- a$y[lab]
    text(wl$x+wl$x.off, wl$y+wl$y.off, rownames(wl),
         cex=cexlab, col=collab, ...)
    segments(wl$x+wl$x.off-wl$width/2, wl$y+wl$y.off,
             wl$x.orig, wl$y.orig,
             col=collab)
  }
  invisible(a)
}

#' Test of the pudding
beesbox(x$group, x$value, ylim=c(40, 60), ylab="Intactness",
        col=hcl.colors(length(unique(x$group)), "Dark 2"))
abline(h=50)
#' With labels
lab <- unlist(sapply(unique(x$group), function(i) {
  v <- x$value[x$group == i]
  q <- quantile(v, c(0.005, 0.995))
  which(x$group == i & x$value %)(% q)
}))
a <- beesbox(x$group, x$value, ylim=c(40, 60), ylab="Intactness",
        col=hcl.colors(length(unique(x$group)), "Dark 2"))
abline(h=50)
addlabels(a, lab=lab, yofflab=2)
