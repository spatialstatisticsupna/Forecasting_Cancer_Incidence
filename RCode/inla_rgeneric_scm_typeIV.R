inla.rgeneric.SCM = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL){
  
  envir = parent.env(environment())
  
  interpret.theta = function(){
    return(
      list(prec = exp(theta[1L]),
           a = exp(theta[2L])
      )
    )
    
  }
  
  graph = function(){
    return(Q())
  }
  
  Q = function(){
    param = interpret.theta()
    a.inv <- solve(param$a)[1]
    D <- as.vector(apply(W, 1, sum))
    C0 <- Diagonal(x = D) - W
    ##RW1
    D1 <- diff(diag(k), differences=1)
    R <- t(D1)%*%D1
    C1 <- kronecker(R,C0)
    
    C <- bdiag(param$prec * a.inv**2 * C1 + exp(15) * a.inv**4 * C1,
               exp(15) * C1)
    C[1:(k*nrow(W)),(k*nrow(W)+1):(nrow(W)*k*2)] <- - exp(15) * a.inv**2 * C1
    C[(k*nrow(W)+1):(nrow(W)*k*2), 1:(k*nrow(W))] <- - exp(15) * a.inv**2 * C1
    #t <- k
    Q <- inla.as.sparse(C)
    return(Q)
  }
  
  mu = function(){
    return(numeric(0))
  }
  
  log.norm.const = function(){
    return(numeric(0))
  }
  
  log.prior = function(){
    param = interpret.theta()
    
    res <- dgamma(param$prec, 1, 5e-05, log = TRUE) + theta[1L] +
      dgamma(param$a, 10, 10, log = TRUE) + theta[2L]
    return(res)
  }
  
  initial = function(){
    return(as.vector(initial.values))
  }
  
  quit = function(){
    return(invisible())
  }
  
  if (!length(theta)){
    theta <- initial()
  }
  val <- do.call(match.arg(cmd), args = list())
  return(val)
}