
n <- 50
p <- 10
x <- matrix(rnorm(n * p), n, p)
theta.true <- rnorm(p)
y <- x %*% theta.true + 0.1 * rnorm(n)

cvx.setup.dir <- "/usr/local/MATLAB/cvx/"

cvxcode <- paste("variables theta(p)",
                 "minimize(square_pos(norm(y - x * theta, 2)) / 2 + lam * norm(theta, 1))",
                 sep=";")
lasso <- CallCVX(cvxcode, const.vars=list(p=p, y=y, x=x, lam=2),
                 opt.var.names="theta", setup.dir=cvx.setup.dir)

