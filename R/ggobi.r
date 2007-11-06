# Open ggobi with matching R plot
#
# @keyword internal
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