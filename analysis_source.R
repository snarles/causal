####
##  Find small windows with high concentration
##  Dedicated to XS--Xie Shuo!!
####

getwind <- function() {
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
