
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

#             res = solve.QP(R1, dvec, Amat, bvec, meq=0, factorized=T)     
#             u   = res$solution
#             u   = matrix(u, ncol = tau, byrow=T)
#             val = as.numeric(res$value + t(g) %*% Q %*% g)
#             if (is.null(.Object@SOLVER$cvx.setup.dir)) return(0)
#            
#             cvx.setup.dir <- "/usr/local/MATLAB/cvx/"
#             cvxcode <- paste("variables u(n)",
#                              "minimize(-dvec' * u + u' * H * u / 2 + norm(u,1))",
#                              "subject to\n u >= -1; u <= 1",
#                              sep=";")
#             res <- CallCVX(cvxcode, const.vars=list(dvec = dvec, H = as.matrix(H), n = npars),
#                              opt.var.names="u", setup.dir=cvx.setup.dir)
