rm(list=ls())
set.seed(33)

library(ggplot2)

## Overview
# Define kernel
# Define mean function
# Collocate X values
# Calculate covariance matrix
# Draw the jointly multivariate normal samples
# Plot that ish

# Parameters
t = c(-1,1) # plotting range
D = 50*floor(abs(t[1]-t[2]))
lkern = 0.1 # kernel length-scale hyperparameter
sigkern = 2 # kernel varaince hyperparameter

# Define kernel
kernel = function(x0, x1){
    # Using squared exponential loss
    return(sigkern^2 * exp(-1/2 * (x0-x1)^2 / lkern))
}

# Define mean function
mu = function(x){
    # Constant mean
    return(0)
}

# Box-Meuler
boxm = function(us){
    # takes 2 unif(0,1) and returns 2 indep std normal
    u0 = us[1]
    u1 = us[2]
    r = sqrt(-2*log(u0))
    thet = 2*pi*u1
    zs = r * c(cos(thet), sin(thet))
    return(zs)
}

# Collocate X values
X = seq(t[1], t[2], length.out=D)

# Calculate covariance matrix
S = outer(X, X, kernel)

## Draw joint multivariate gaussian samples
# Draw an even number of uniforms
d = D + (D %% 2)
u = matrix(runif(d, 0, 1), d/2, 2)
# Use them to draw std normals from the Box-Meuler transform
z = c(apply(u, 1, boxm))
