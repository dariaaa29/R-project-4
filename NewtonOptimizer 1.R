## names: Yiqian Dai(s2329960); Maki Yoshida(S2346588)
## github repo address: https://github.com/dariaaa29/R-project-4.git
## team members contributions:

# Code to implemente Newton’s method for minimization of functions.
# As an iterative algorithm, Newton's method requires Hessian matrix
# and gradient values for the magnitude of each parameter update.
# In this task, the Hessian matrix may not be provided. At this time,
# the method calc_hess() to calculate the Hessian matrix based on the
# gradient value difference is implemented. Another big challenge is that
# inverse matrices cannot be solved directly using solve(), so Choleskey
# decomposition is used here to calculate the inverse of the matrix.
# Besides, in order to make the program more user-friendly, corresponding
# warnings and error prompts are added.


# test function
rb <- function(th, k = 2) {
  k * (th[2] - th[1] ^ 2) ^ 2 + (1 - th[1]) ^ 2
}

# Gradient function of the test function
gb <- function(th, k = 2) {
  c(-2 * (1 - th[1]) - k * 4 * th[1] * (th[2] - th[1] ^ 2), k * 2 * (th[2] -
                                                                       th[1] ^ 2))
}

# second derivative of test function
hb <- function(th, k = 2) {
  h <- matrix(0, 2, 2)
  h[1, 1] <- 2 - k * 2 * (2 * (th[2] - th[1] ^ 2) - 4 * th[1] ^ 2)
  h[2, 2] <- 2 * k
  h[1, 2] <- h[2, 1] <- -4 * k * th[1]
  h
}

# Find the Hessian matrix based on the difference of the gradient function
calc_hess = function(theta, grad, ..., eps = 1e-6) {
  # theta is the parameter value, grad is the gradient function
  # len gets the dimension of the parameter
  len = length(theta)
  # Substitute the current parameter into the gradient function to obtain the subtrahend
  minus_tmp = grad(theta, ...)
  Hi = NULL
  for (i in 1:len) {
    # Copy theta to avoid modification in each step
    theta_copy = theta
    # Calculate the partial derivative corresponding to each column of the Hessian matrix
    theta_copy[i] = theta_copy[i] + eps
    A = grad(theta_copy, ...) - minus_tmp
    # merge the value in the Hi
    Hi = cbind(Hi, A)
  }
  # The obtained Hi is finally uniformly divided by eps to obtain the entire Hessian matrix
  Hi / eps
}

newt = function(theta,
                func,
                grad,
                hess = NULL,
                ...,
                tol = 1e-8,
                fscale = 1,
                maxit = 100,
                max.half = 20,
                eps = 1e-6) {
  # Newton's method main function, optimize the objective function func and return the corresponding theta
  # First determine whether the objective function is defined at the initial point
  # if not, an error will be reported directly
  test_theta = tryCatch({
    obj_value = func(theta, ...)
    grad_value = grad(theta, ...)
  }, error = function(e) {
    stop("The objective or derivatives are not finite at the initial theta.")
    
  })
  # obj_values store all the vlues of objective function with each step
  obj_values = NULL
  # minus_flag is used to adjust the update direction
  minus_flag = 1
  # it is the number of iterations
  it = 0
  # theta_new is used to store the parameter theta after each update
  theta_new = theta
  while (it < maxit) {
    # Current objective function value
    obj_value = func(theta_new, ...)
    # Store all objective function values to determine whether the convergence is correct
    obj_values = cbind(obj_values, obj_value)
    # Calculate the gradient value (this is a column vector)
    grad_value = matrix(grad(theta_new, ...))
    # calculate the Hessian matrix
    if (is.null(hess)) {
      hess_value = calc_hess(theta_new, grad, ..., eps = eps)
    } else{
      hess_value = hess(theta_new, ...)
    }
    
    # Compute the inverse of a Hessian matrix using Choleskey decomposition
    test_theta = tryCatch({
      hess_inv = chol2inv(chol(hess_value))
    }, error = function(e) {
      # An error occurs indicating that the Hessian matrix does not satisfy positive definiteness
      # an error is reported
      stop("the Hessian is not positive definite.")
    })
    # Determine whether the convergence end condition is met
    if (all(grad_value < tol * (obj_value + fscale))) {
      break
    }
    # If the value of the objective function does not drop by half after the allowed number of steps
    # edit the minus_flag
    if (it > max.half) {
      if (obj_value > obj_values[it - max.half] / 2) {
        minus_flag = -1
      }
      # a warning is reported
      warning(
        "The step fails to reduce the objective despite during max.half step halvings.",
        immediate. = TRUE
      )
    }
    # update the theta
    theta_new = theta_new - minus_flag * hess_inv %*% grad_value
    it = it + 1
  }
  # If the maximum allowable steps of the loop are reached, a warning is displayed
  if (it == maxit) {
    warning("The maxit is reached without convergence.", immediate. = TRUE)
  }
  # return a list
  list(
    f = obj_value,
    theta = theta_new,
    iter = it,
    g = grad_value,
    Hi = hess_inv
  )
}

newt(c(2, 3), rb, gb, hb, 3)

