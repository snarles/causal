####
##  Code used for triplots
####


subsample <- function(dat, p = 0.3) {
  dat[sample(dim(dat)[1], floor(dim(dat)[1]*p), replace = FALSE), ]
}

plot1 <- function(fmla, data, ...) plot(fmla, data = data,
                                        ylim = yl, xlim = xl, col = col0, ...)
plot2 <- function(fmla, data, fn = plot, ...) fn(fmla, data = data,
                                                 ylim = yl, xlim = xl, ...)

col0 <- rgb(0, 0, 0, 0.8)
col3 <- hsv(h = (1:3)/3, s = 1, v = 1, 0.8)

triplot <- function(fmla, dat) {
  par(bg = grey(0.3))
  plot2(fmla, dat[[1]], col = col3[1])
  plot2(fmla, dat[[2]], col = col3[2], fn = points)
  plot2(fmla, dat[[3]], col = col3[3], fn = points)
  plot2(fmla, subsample(dat[[2]], 0.4), col = col3[1], fn = points)
  plot2(fmla, subsample(dat[[2]], 0.3), col = col3[2], fn = points)
  plot2(fmla, subsample(dat[[2]], 0.2), col = col3[3], fn = points)
}

####
##  Find small windows with high concentration
##  Dedicated to XS--Xie Shuo!!
####

marginal_plot <- function() {
  smoothScatter(fm, pch = 19, cex = 0.3)
  points(fm, pch = ".", col = grey(0, 0.7))
  cc <- cor(model.frame(fm, data = sachs_all))[1, 2]
  cc <- substr(as.character(cc), 1, 5)
  title(paste("R=", cc))
}

getwind <- function() {
  plot(cf, data = sachs_all, pch = ".")
  xs <- locator(n = 2)
  xs$x <- sort(xs$x); xs$y <- sort(xs$y)
  polygon(c(xs$x, rev(xs$x)), rep(xs$y, each = 2), border = "red", lwd = 2, col = NA)
  lala <- model.frame(cf, data = sachs_all)[, c(2, 1)]
  filt <- (lala[, 1] > xs$x[1] & lala[, 1] < xs$x[2]) & 
    (lala[, 2] > xs$y[1] & lala[, 2] < xs$y[2])
  print(sum(filt))
  points(lala[filt, ], pch = ".", col = rgb(1, 0, 0, 0.3))
  xs
}

condition_w <- function(w, marginal = FALSE, ...) {
  lala <- model.frame(cf, data = sachs_all)[, c(2, 1)]
  filt <- (lala[, 1] > w$x[1] & lala[, 1] < w$x[2]) & 
    (lala[, 2] > w$y[1] & lala[, 2] < w$y[2])
  cc <- cor(model.frame(fm, data = sachs_all)[filt, ])[1, 2]
  cc <- substr(as.character(cc), 1, 5)
  if (marginal) {
    smoothScatter(fm, pch = ".")
    points(fm, data = sachs_all[filt, ], pch = 19, cex = 0.4)
  } else {
    plot(fm, data = sachs_all[filt, ], main = paste("R=", cc), ...)
  }
}

printw <- function(w) {
  cat(paste0("list(x=c(", w$x[1], ",", w$x[2], "),y=c(", w$y[1], ",", w$y[2], "))"))
}

####
##  Invariance stuff
####

quadd <- function(X) {
  cbind(X, x11 = X[,1]^2, x22 = X[, 2]^2, x12 = X[, 1] * X[, 2])
}
