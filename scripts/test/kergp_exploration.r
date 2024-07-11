#
# kergp_exploration.r
# Exploring the `kergp` Gaussian process (GP) package. In particular, investigating 
# a method to include a fixed jitter in the covariance during maximum 
# likelihood estimation (MLE) and GP prediction. kergp includes the argument 
# `noise` which can be set during MLE to either estimate a nugget variance, 
# or exclude it entirely, but no way to fix it at a desired value. It includes 
# the argument `forceInterp` during prediction which incorporates a white 
# noise kernel with stationary variace equal to the nugget variance, which is 
# not the behavior I typically want when predicting with a nugget variance/jitter. 
#
# Andrew Roberts
#

library(kergp)

# ------------------------------------------------------------------------------
# Synthetic data for GP regressiont tests. 
# ------------------------------------------------------------------------------

# Settings.   
N_design <- 4
N_test <- 100
f <- function(x) 10 * 2*(x - 10)^2

# Training data. 
X <- matrix(seq(-1,1, length.out=N_design), ncol=1)
colnames(X) <- "x"
y <- f(X) # + 100*rnorm(n=N_design)
df <- data.frame(x=drop(X), y=drop(y))

# Test data. 
X_test <- matrix(seq(-1, 1, length.out=N_test), nrow=N_test)
y_test <- f(X_test)
colnames(X_test) <- "x"

# ------------------------------------------------------------------------------
# Custom kernel functions.  
# ------------------------------------------------------------------------------

# Gaussian with jitter added to diagonal. This is a slighly modified version 
# kergp's `kGauss` function found in kergp/R/kerenelNorm.R. The method of 
# adding a jitter here is a hack for convenience, and may not always be reliable. 
# In particular, the jitter is added whenever the two input matrices `x1` and 
# `x2` have the same number of rows. In reality, it should only be added 
# when the matrices actually contain the same exact observations. When using 
# this hack, one should therefore avoid predicting at a number of inputs 
# equal to the number of design points. Also, the jitter function in this 
# kernel is fixed - it cannot be learned. 

kGauss_jitter <- function(d = 1, jitter=sqrt(.Machine$double.eps), 
                          input_names=paste("theta", 1L:d, sep = "_")) {
  kernFun <- function(x1, x2, par) {
    K12 <- kergp:::kNormFun(x1, x2, par, k1FunGauss)
    
    # Hack: add jitter to diagonal if number of rows are equal. 
    if(nrow(x1) == nrow(x2)) K12 <- hetGP:::add_diag(K12, rep(jitter, nrow(x1)))
    return(K12)
  }
  
  k <- covMan( 
    kernel = kernFun,
    hasGrad = TRUE,
    acceptMatrix = TRUE,
    d = d,
    par = c(rep(1, d), 1),
    parLower = rep(1e-8, d + 1L),
    parUpper = rep(Inf, d + 1L),
    parNames = c(input_names, "sigma2"),
    label = "Gauss kernel"
  )
  
  inputNames(k) <- input_names
  return(k)
}

# ------------------------------------------------------------------------------
# Fit GP.   
# ------------------------------------------------------------------------------

ker_no_jitter <- kGauss()
inputNames(ker_no_jitter) <- colnames(X)
gp_no_jitter <- kergp::gp(y~1, data=df, inputs="x", cov=ker_no_jitter, estim=TRUE,
                          noise=FALSE, trace=TRUE, multistart=10)

ker_jitter <- kGauss_jitter(input_names=colnames(X))
gp_jitter <- kergp::gp(y~1, data=df, inputs="x", cov=ker_jitter, estim=TRUE, 
                       noise=FALSE, trace=TRUE, multistart=10)

# ------------------------------------------------------------------------------
# Predict GP.   
# ------------------------------------------------------------------------------

pred_no_jitter <- kergp:::predict.gp(gp_no_jitter, newdata=X_test, seCompute=TRUE, 
                                     covCompute=TRUE, forceInterp=FALSE)

pred_jitter <- kergp:::predict.gp(gp_jitter, newdata=X_test, seCompute=TRUE, 
                                  covCompute=TRUE, forceInterp=FALSE)

plot(X_test, y_test, col="red", type="l")
points(X, y)
lines(X_test, pred_no_jitter$mean, col="blue")
lines(X_test, pred_no_jitter$mean + 1.64*pred_no_jitter$sd, col="gray")
lines(X_test, pred_no_jitter$mean - 1.64*pred_no_jitter$sd, col="gray")

plot(X_test, y_test, col="red", type="l")
points(X, y)
lines(X_test, pred_jitter$mean, col="blue")
lines(X_test, pred_jitter$mean + 1.64*pred_no_jitter$sd, col="gray")
lines(X_test, pred_jitter$mean - 1.64*pred_no_jitter$sd, col="gray")


# ------------------------------------------------------------------------------
# Quadratic kernel.    
# ------------------------------------------------------------------------------

jitter <- sqrt(.Machine$double.eps)

quad_ker_func <- function(x1, x2, par) { 
  affine_comb <- sum(x1 * x2) + par[1]
  kern <- affine_comb^2
  attr(kern, "gradient") <- c(cst=2*affine_comb)
  
  # Hack: add jitter to diagonal if number of rows are equal. 
  if(nrow(x1) == nrow(x2)) kern <- hetGP:::add_diag(kern, rep(jitter, nrow(x1)))
  
  return(kern)
}

quad_ker <- covMan(kernel = quad_ker_func,
                   hasGrad = TRUE,
                   d = 1L,
                   parLower = c(cst=-Inf),
                   parUpper = c(cst=Inf),
                   parNames = c("cst"),
                   label = "quadratic kernel with jitter", 
                   inputs = "x", 
                   par = c(cst=1.0))








