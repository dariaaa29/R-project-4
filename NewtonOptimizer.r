# names:
# github repo address
# team members contributions:



rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

calc_hess=function(theta,grad,...,eps=1e-6){
  # theta是参数值，grad是梯度函数
  len = length(theta)
  minus_tmp = grad(theta,...)
  Hi = NULL
  for(i in 1:len){
    theta_copy = theta
    theta_copy[i]=theta_copy[i]+eps
    A = grad(theta_copy,...)-minus_tmp
    Hi = cbind(Hi,A)
  }
  Hi/eps
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
  test_theta = tryCatch({
    obj_value = func(theta,...)
    grad_value = grad(theta,...)
  },error = function(e){
    stop("The objective or derivatives are not finite at the initial theta.");
  })
  
  obj_values = NULL
  minus_flag = 1
  it = 0
  theta_new = theta
  while (it < maxit) {
    # 目标函数值
    obj_value = func(theta_new, ...)
    # 将所有的目标函数值储存起来，用于判断是否收敛正确
    obj_values = cbind(obj_values,obj_value)
    # 梯度值，列向量
    grad_value = matrix(grad(theta_new, ...))
    # Hessian矩阵
    if (is.null(hess)) {
      hess_value = calc_hess(theta_new, grad, ..., eps = eps)
    } else{
      hess_value = hess(theta_new, ...)
    }
    
    # 对Hessian矩阵求逆
    # 可以用Choleskey分解来计算矩阵的逆，这时候可以用到函数chol2inv()
    hess_inv = chol2inv(chol(hess_value))
    # 判断是否结束
    if (all(grad_value < tol * (obj_value + fscale))) {
      break
    }
    if (it>max.half){
      if(obj_value>obj_values[it-max.half]/2){
        minus_flag = -1
      }
    }
    # 迭代更新参数
    theta_new = theta_new - minus_flag * hess_inv %*% grad_value
    it = it + 1
    #print(theta_new)
  }
  if (it==maxit){
    warning("maxit is reached without convergence.", immediate. = TRUE)
  }
  list(
    f = obj_value,
    theta = theta_new,
    iter = it,
    g = grad_value,
    Hi = hess_inv
  )
}

newt(c(2, 3), rb, gb, hb, 3)
