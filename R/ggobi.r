# Open ggobi with matching R plot
# Experimental interface to interact with scagnostics output with GGobi
#
# This is an experimental interface to explore the scagnostics output
# with GGobi, taking advantage of the fact that each point in the output
# of \code{\link{scagnostics}} corresponds to a plot.  To use it, switch into
# identify mode and hover over the plot (= point) that you want to see.
# 
# @keyword dynamic
ggobi.scagdf <- function(data, ...) {
  vars <- attr(data, "vars")
  orig <- attr(data, "data")
  data <- rbind(data, rep(0, ncol(data)))
  data <- rbind(data, rep(1, ncol(data)))
  
  g <- ggobi(data)
  
  shadowed(g[1]) <- c(rep(FALSE, nrow(data) - 2), TRUE, TRUE)
  glyph_type(g[1]) <- c(rep(6, nrow(data) - 2), 1, 1)

  gSignalConnect(g, "identify-point", function(gg, plot, id, d) { 
    if (id != -1) {
      v <- vars[id, , drop=TRUE]
      if (is.null(v)) return()
      plot(orig[, v$x], orig[, v$y], xlab=names(orig)[v$x], ylab=names(orig)[v$y])
    }  
  })

  rownames(vars) <- rownames(data)[1:(nrow(data) - 2)]
  g$vars <- vars

  invisible(g)
}