library(ape)
#Plotting phylogenetic time axis with dates instead of numerical scale.
#' Title
#'
#' @param side
#' @param root.time
#' @param backward
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
axisPhyloDate <- function(side = 1, root.time = NULL, backward = TRUE, ...)
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type

  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")

  if (is.null(root.time)) root.time <- lastPP$root.time

  if (type %in% c("phylogram", "cladogram")) {
    xscale <-
      if (lastPP$direction %in% c("rightwards", "leftwards")) range(lastPP$xx)
    else range(lastPP$yy)

    tmp <- lastPP$direction %in% c("leftwards", "downwards")

    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp)) tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward) tscale <- tscale - xscale[2]
    }

    ## the linear transformation between the x-scale and the time-scale:
    beta <- diff(xscale) / diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]

    lab <- pretty(tscale)
    x <- beta * lab + alpha
    axis(side = side, at = x, labels = FALSE,...)
    text(x,  par("usr")[3], labels =as.Date(date_decimal(lab)), srt=60,adj=1.1, xpd=TRUE)
  } else { # type == "fan"
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]; yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    firstandlast <- c(1, n)
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0); y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360*theta0/(2*pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0 # the length of tick marks
    r <- r0
    while (r > 1e-8) {
      x <- r * cos(theta0); y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2); thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa); ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5, 1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}



