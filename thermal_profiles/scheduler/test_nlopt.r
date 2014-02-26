eval_f0 <- function( x, a, b ){ 
  return( sqrt(x[2]) )
}

# constraint function
eval_g0 <- function( x, a, b ) {
  return( (a*x[1] + b)^3 - x[2] )
}

# gradient of objective function
eval_grad_f0 <- function( x ){ 
  return( c( 0, .5/sqrt(x[2]) ) )
}

# jacobian of constraint
eval_jac_g0 <- function( x, a, b ) {
  return( rbind( c( 3*a[1]*(a[1]*x[1] + b[1])^2, -1.0 ), 
                 c( 3*a[2]*(a[2]*x[1] + b[2])^2, -1.0 ) ) )
}


# functions with gradients in objective and constraint function
# this can be useful if the same calculations are needed for
# the function value and the gradient
eval_f1 <- function( x, a, b ){ 
  return( list("objective"=sqrt(x[2]), 
               "gradient"=c(0,.5/sqrt(x[2])) ) )
}

eval_g1 <- function( x, a, b ) {
  return( list( "constraints"=(a*x[1] + b)^3 - x[2],
                "jacobian"=rbind( c( 3*a[1]*(a[1]*x[1] + b[1])^2, -1.0 ), 
                                  c( 3*a[2]*(a[2]*x[1] + b[2])^2, -1.0 ) ) ) )
}


# define parameters
a <- c(2,-1)
b <- c(0, 1)

# Solve using NLOPT_LD_MMA with gradient information supplied in separate function
res0 <- nloptr( x0=c(1.234,5.678), 
                eval_f=eval_f0, 
                eval_grad_f=eval_grad_f0,
                lb = c(-Inf,0), 
                ub = c(Inf,Inf), 
                eval_g_ineq = eval_g0,
                eval_jac_g_ineq = eval_jac_g0,                
                opts = list("algorithm"="NLOPT_LD_MMA"),
                a = a, 
                b = b )
print( res0 )
