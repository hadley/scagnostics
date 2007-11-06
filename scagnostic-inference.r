library(scagnostics)
set.seed(10)

n <- 250
ymany <- as.data.frame(t(replicate(100, scagnostics(rnorm(n), rnorm(n))$s)))


# qplot(value, data=melt(ymany), facets=. ~ variable)

qplot(ymany$Outlying , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Skewed   , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Clumpy   , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Sparse   , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Striated , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Convex   , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Skinny   , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Stringy  , geom="histogram", xlim=c(0, 1), binwidth=0.05)
qplot(ymany$Monotonic, geom="histogram", xlim=c(0, 1), binwidth=0.05)

