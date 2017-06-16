rm(list=ls())
set.seed(33)

library(ggplot2)
library(MASS)
library(reshape)

## Overview
# Define kernel
# Define mean function
# Collocate X values
# Calculate covariance matrix
# Draw the jointly multivariate normal samples
# Plot that ish

# Parameters
t = c(-1,1) # plotting range
D = 50*floor(abs(t[1]-t[2])) # number of evenlt spaced collocation points
K = 10 # number of curves
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
    return(0*x)
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
xs = seq(t[1], t[2], length.out=D)

# Calculate covariance matrix
S = outer(xs, xs, kernel)

## Draw joint multivariate gaussian samples
# # Draw an even number of uniforms
# d = D + (D %% 2)
# u = matrix(runif(d, 0, 1), d/2, 2)
# # Use them to draw std normals from the Box-Meuler transform
# z = c(apply(u, 1, boxm))
# # Calculate the joint multiv gaus
# A = chol(S)
# ys = mu(xs) + A %*% z
ys = t(mvrnorm(K, mu = mu(xs), Sigma=S))

# Plot
plt_df = data.frame(x=xs, y=ys)
plt_df_m = melt(plt_df, id=c("x"))
plt = ggplot(plt_df_m) +
    geom_line(aes(x=x, y=value, color=variable))
