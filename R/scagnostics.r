# Calculcate scagnostics for a scatterplot
# Scagnostics summarise potentially interesting patterns in 2d scatterplot
# 
# Current scagnostics are: 
# 
# \itemize{
#   \item Outlying
#   \item Skewed
#   \item Clumpy
#   \item Sparse
#   \item Striated
#   \item Convex
#   \item Skinny
#   \item Stringy
#   \item Monotonic
# }
# 
# These are described in more detail in: 
# Graph-Theoretic Scagnostics, Leland Wilkinson, Anushka 
# Anand, Robert Grossman. 
# \url{http://www.rgrossman.com/dl/proc-094.pdf}
# 
# You can call the function with two 1d vectors to get a single vector
# of scagnostics, or with a 2d structure (matrix or data frame) to get 
# scagnostics for every combination of the variables.
# 
# @arguments object to calculate scagnostics on: a vector, a matrix or a data.frame
# @arguments ...
# @keyword hplot
# @alias scagnostics.default
# @alias scagnostics.matrix
# @alias scagnostics.data.frame
# @alias scagnostics_2d
#X scagnostics(1:10, 1:10)
#X scagnostics(rnorm(100), rnorm(100))
#X scagnostics(mtcars)
#X scagnostics(as.matrix(mtcars))
#X 
#X if (require(rggobi)) ggobi(scagnostics(mtcars))
scagnostics <- function(x, ...) UseMethod("scagnostics", x)


scagnostics.default <- function(x, y, bins=50, ...) {
  stopifnot(length(x) == length(y))

  complete <- !is.na(x) & !is.na(y)  
  x <- x[complete]
  y <- y[complete]
  
  x <- (x - min(x)) / diff(range(x))
  y <- (y - min(y)) / diff(range(y))
  
  results <- rep(0, 9 + 3 * 1000)
  r <- .C("scagnostics",
    x = as.double(x),
    y = as.double(y),
    length = as.integer(length(x)),
    bins = as.integer(bins),
    results = as.double(results)
  )$results
  
  n <- r[10]
  s <- r[1:9]
  bins <- matrix(r[11:(10 + n * 3)], ncol=3)
  
  names(s) <- c("Outlying", "Skewed", "Clumpy", "Sparse", "Striated", "Convex", "Skinny", "Stringy", "Monotonic")

  list(s=s, bins=bins)
}

scagnostics.matrix <- function(x, ...) {
  scagnostics_2d(x, ...)
}
scagnostics.data.frame <- function(x, ...) {
  scagnostics_2d(x, ...)
}
scagnostics_2d <- function(x, ...) {
  vars <- expand.grid(x=1:ncol(x), y=1:ncol(x))
  vars <- vars[vars$x < vars$y, ]
  
  each <- lapply(seq_len(nrow(vars)), function(i) {
    scagnostics(x[,vars[i, 1]], x[,vars[i, 2]])$s
  })
  scag <- do.call("rbind", each)
  
  rownames(scag) <- apply(vars, 1, 
    function(v) paste(colnames(x)[v], collapse=" vs "))

  attr(scag, "vars") <- vars
  attr(scag, "data") <- x
  structure(scag, class = c("scagdf"))
}

# Print scagnostics data structure
# @keyword internal
print.scagdf <- function(x, ...) {
  attr(x, "vars") <- NULL
  attr(x, "data") <- NULL
  print.default(x, ...)
}

# Calculate scagnostics while tour is running
# @keyword internal
scagnostics.tour <- function(x, ...) { 
  if (!require("rggobi")) stop("rggobi required for the tour")
  
  g <- ggobi(x)
  d <- displays(g)[[1]]

  id <- gTimeoutAdd(1000, function(x) {print("."); TRUE}, data = NULL) 

  gSourceRemove(id)
}