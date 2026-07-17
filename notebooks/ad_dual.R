## ---------------------------------------------------------------------------
## Forward-mode automatic differentiation by operator overloading.
## A dual number carries a value v and a vector d of partial derivatives with
## respect to 5 seed directions: (T, S, CT, AT, h).
## ---------------------------------------------------------------------------
ND <- 5L
# v is a numeric vector of length n; d is an n x ND matrix of partials.
dual    <- function(v, d) structure(list(v = v, d = d), class = "dual")
as.dual <- function(x, n) if (inherits(x, "dual")) x else
                          dual(rep(x, length.out = n), matrix(0, n, ND))
seed    <- function(v, i) { d <- matrix(0, length(v), ND); d[, i] <- 1; dual(v, d) }

Ops.dual <- function(e1, e2) {
  if (missing(e2)) {                                   # unary + and -
    if (.Generic == "-") return(dual(-e1$v, -e1$d))
    if (.Generic == "+") return(e1)
    stop("unsupported unary op: ", .Generic)
  }
  n <- max(if (inherits(e1,"dual")) NROW(e1$d) else 1L,
           if (inherits(e2,"dual")) NROW(e2$d) else 1L)
  a <- as.dual(e1, n); b <- as.dual(e2, n)
  switch(.Generic,
    "+" = dual(a$v + b$v, a$d + b$d),
    "-" = dual(a$v - b$v, a$d - b$d),
    "*" = dual(a$v * b$v, a$d * b$v + a$v * b$d),
    "/" = dual(a$v / b$v, (a$d * b$v - a$v * b$d) / b$v^2),
    "^" = { v <- a$v^b$v
            dual(v, v * (b$d * log(a$v) + b$v * a$d / a$v)) },
    stop("unsupported op: ", .Generic))
}

Math.dual <- function(x, ...) {
  switch(.Generic,
    "exp"  = dual(exp(x$v),  exp(x$v) * x$d),
    "log"  = dual(log(x$v),  x$d / x$v),
    "sqrt" = dual(sqrt(x$v), x$d / (2 * sqrt(x$v))),
    stop("unsupported math fn: ", .Generic))
}
val <- function(x) if (inherits(x, "dual")) x$v else x
