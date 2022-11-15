
#the arguments passed to the newt function
#tol <- #the convergence tolerance
#fscale <- #a rough estimate of the magnitude of func near the optimum - used in convergence testing.
#maxit <-#the maximum number of Newton iterations to try before giving up
#max.half <- #the maximum number of times a step should be halved before concluding that the step has failed toimprove the objective.
#eps <- #the finite difference intervals to use when a Hessian function is not provided.
#newt <- # should return a list containing:
#f <- #
#iter #the number of iterations taken to reach the minimum.
#theta <- vector() #  the value of the parameters at the minimum
#g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
#Hi the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).

# This will be used as the vector passed for the function
# ((t(A) + A) / 2)


# This function is implementing Newton's method for minimizing of functions.
# It will return a list
newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  
  # If the objective or derivatives are not finite at the initial theta
  
#  if (is.infinite(theta) ) # error using stop or warning
  
  # If the step fails to reduce the objective despite trying max.half step halvings # error using stop or warning
  
  # If maxit is reached without convergence; # error using stop or warning
  
  #If the Hessian is not positive definite at convergence# error using stop or warning
  
  
  
}


optimalvec <- vector()
func <- function(optimalvec, ...) {
  
  
  
}

# This is the gradient function which return the gradient vector of the object w.r.t
grad <- function(theta,...) {
  
}

# This is the Hessian matrix function which return the Hessian matrix of the objective w.r.t
hess <- function(theta,...) {
  
}

# test function
rb <- function(th,k=2) {
  k*(th[2]-th[1]ˆ2)ˆ2 + (1-th[1])ˆ2
}

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]ˆ2),k*2*(th[2]-th[1]ˆ2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]ˆ2) - 4*th[1]ˆ2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}












