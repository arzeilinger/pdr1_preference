#### Functions to plot rate parameters and equlilibrial probabilities using lattice

### Panel function to remove plot border
mpanel = function(...) {grid.segments(x0 = c(0,0), x1 = c(1,0), y0 = c(0,0),
                                      y1 = c(0,1))
  panel.xyplot(...)}

### Panel function to produce error bars
prepanel.ci <- function(x, y, ly, uy, subscripts, ...)
{
  y <- as.numeric(y)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  list(ylim = range(y, uy, ly, finite = TRUE))
}

panel.ci <- function(x, y, ly, uy, subscripts, pch = 16, ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  panel.arrows(x, ly, x, uy, col = 'black',
               length = 0.1, unit = "native",
               angle = 90, code = 3)
  panel.xyplot(x, y, pch = pch, ...)
}
