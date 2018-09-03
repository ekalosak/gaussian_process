# Author: Eric Kalosa-Kenyon
# Supporting code for the 1D gaussian process blog post
#
# Methods: TODO:
# Plots: TODO:

## @knitr load_libraries
set.seed(33)

library(ggplot2)
library(MASS)
library(reshape)
library(latex2exp)

## @knitr set_params1
# Parameters
t = c(-pi,pi) # plotting range
D = 20*floor(abs(t[1]-t[2])) # number of evenlt spaced collocation points
K = 10 # number of curves
lkern = 0.2 # kernel length-scale hyperparameter
sigkern = 2.7 # kernel varaince hyperparameter

# Define kernel
kernel_sel = function(x0, x1){
    # Using squared exponential loss
    return(sigkern^2 * exp(-1/2 * (x0-x1)^2 / lkern))
}

# Define mean function
mu = function(x){
    # Constant mean
    return(0*x)
}

# Collocate X values
xs = seq(t[1], t[2], length.out=D)

# Calculate covariance matrix
S = outer(xs, xs, kernel_sel)

# Draw joint multivariate gaussian samples
ys = t(mvrnorm(K, mu = mu(xs), Sigma=S))

## @knitr plot_basic_gp
# Plot
plt_df = data.frame(x=xs, y=ys)
names(plt_df) = c("x", paste("GP", 1:K, sep=""))
plt_df_m = melt(plt_df, id=c("x"))
plt_basic = ggplot(plt_df_m) +
    geom_line(aes(x=x, y=value, color=variable)) +
    labs(x="X", y="Y",
         title=paste(K, "draws from a gaussian process"),
         color="Draw number") +
    guides(color=FALSE)

plt_basic

## @knitr alt_kernel
kernel_01 = function(x0, x1){
    # Using 0-1 loss
    return((x0 == x1) + 0)
}

# Calculate covariance matrix
S = outer(xs, xs, kernel_01)

# Draw joint multivariate gaussian samples
ys = t(mvrnorm(K, mu = mu(xs), Sigma=S))

## @knitr plot_01_kernel
# Plot
plt_df = data.frame(x=xs, y=ys)
names(plt_df) = c("x", paste("GP", 1:K, sep=""))
plt_df_m = melt(plt_df, id=c("x"))
plt_altk = ggplot(plt_df_m) +
    geom_line(aes(x=x, y=value, color=variable)) +
    labs(x="X", y="Y",
         title=paste(K, "draws from a gaussian process"),
         color="Draw number") +
    guides(color=FALSE) # remove legend for color variable

plt_altk

## @knitr plot_01_kernel_sideways
plt_altks = ggplot(plt_df_m) +
    geom_freqpoly(aes(value, ..density.., color=variable), bins=floor(D/5)) +
    stat_function(fun = dnorm, colour = "steelblue", size=1.1) +
    labs(title="Lateral density of the 01-kernel GP",
         x="Y", y="Density") +
    guides(color=FALSE)

plt_altks

## @knitr abs_kernel
kernel_abs = function(x0, x1){
    # Using exponential loss
    return(4^2 * exp(-1/2 * abs(x0-x1) / 0.5))
}

# Calculate covariance matrix
S = outer(xs, xs, kernel_abs)

# Draw joint multivariate gaussian samples
ys = t(mvrnorm(K, mu = mu(xs), Sigma=S))

## @knitr plot_abs_kernel
# Plot
plt_df = data.frame(x=xs, y=ys)
names(plt_df) = c("x", paste("GP", 1:K, sep=""))
plt_df_m = melt(plt_df, id=c("x"))
plt_absk = ggplot(plt_df_m) +
    geom_line(aes(x=x, y=value, color=variable)) +
    labs(x="X", y="Y",
         title=paste(K, "draws from a gaussian process"),
         color="Draw number") +
    guides(color=FALSE)

plt_absk

## @knitr gaus_spline
# Parameters
M = 7 # number of points to observe from the model
# NOTE: that if this gets too big, inverting the covariance matrix can bevome
# numerically unstable resulting in failures to invert (which is ok) or wildly
# divergent estimates (which could be catastrophic).

# Heterogeneous mean function
mu = function(x){
    return(x^2 * exp(-x^2) + sin(x))
}

ys = mu(xs)
# ixs = base::sample(1:D, M)
ixs = ceiling(seq(1, D, length.out=M))
xsobv = xs[ixs]
ysobv = ys[ixs] + rnorm(M, sd=0.4)
xu = t[2] + 1/2

## @knitr plot_regression_setup
# Plot the true function and the observed points
plt_reg = ggplot(data.frame(x=xs, y=ys)) +
    geom_line(aes(x=x, y=y), color="steelblue", size=1.1) +
    geom_point(data=data.frame(x=xsobv, y=ysobv), aes(x=x, y=y), color="coral") +
    geom_vline(xintercept=xu, color="darkviolet") +
    labs(
         title="Mean function and noisy observations",
         x="X", y="Y"
         )

plt_reg

## @knitr kreiging
gen_posterior = function(xu, xs, ys, ker=kernel_sel){
    ## Generates a posterior function

    # posterior probability f(yu | ys) = N(m, s)
    K0 = outer(xs, xs, ker) # kernel of the observed values
    K1 = outer(xu, xs, ker)
    K2 = outer(xu, xu, ker)
    K = rbind(cbind(K0, t(K1)), cbind(K1, K2))

    # mean and variance of posterior
    m = K1 %*% solve(K0) %*% ys
    s = K2 - K1 %*% solve(K0) %*% t(K1)

    return(
           function(n){
               mvrnorm(n, mu = m, Sigma = s)
           }
    )
}

posterior = gen_posterior(xu, xsobv, ysobv)

n = 100
post_obs = posterior(n)

## @knitr plot_posterior
plt_post = ggplot(data=data.frame(y=post_obs)) +
    geom_histogram(aes(x=y), bins=floor(n/4)) +
    labs(
         title="Draws from the posterior",
         x="Value", y="Count"
         )

plt_post

## @knitr credible_band
# Parameters
m = 50 # resolution of the posterior
n = 10 # number of draws from the posterior
d = t[2] - t[1] # length of the sampling support

# Collocated posterior points
xus = seq(t[1] - d/10, t[2] + d/10, length.out=m)

# Draw from the posterior for an illustrative plot
posterior_multi = gen_posterior(xus, xsobv, ysobv)
draws = posterior_multi(n)

# Draw for the emperical 90% band and calculate the band
n = 50
draws = posterior_multi(n)
qs = apply(t(draws), 1, (function(x) (return(quantile(x, c(0.05, 0.95))))))


## @knitr plot_credible_band
# Plot the aforementioned illustrative plot of draws from the posterior
df_post = data.frame(t(rbind(draws, xus)))
names(df_post) = c(paste("Draw", 1:n), "X")
df_post_m = melt(df_post, id="X")
plt_cb = ggplot(df_post_m) +
    geom_line(aes(x=X, y=value, color=variable)) +
    geom_point(data=data.frame(x=xsobv, y=ysobv),
               aes(x=x, y=y),
               color="coral") +
    labs(x="X", y="Y", title=paste(n, "draws from the posterior"),
         color="Draw number") +
    guides(color=FALSE)

# Plot the 90% band
df_post_band = data.frame(t(rbind(qs, xus)))
names(df_post_band) = c("q0", "q1", "X")
plt_pb = ggplot(df_post_band) +
    geom_ribbon(aes(ymin=q0, ymax=q1, x=X), alpha=0.6) +
    geom_point(data=data.frame(x=xsobv, y=ysobv),
               aes(x=x, y=y),
               color="coral") +
    geom_line(data=data.frame(x=xs, y=ys), aes(x=x, y=y), color="steelblue") +
    labs(x="X", y="Y", title=paste("90% credible band using", n, "draws"))

plt_cb
plt_pb


## @knitr tune_model
## Grid search to maximize the log likelihood
kern = function(x0, x1, sig, ell){
    # Using squared exponential loss
    return(sig^2 * exp(-1/2 * (x0-x1)^2 / ell))
}

loglik = function(sig, ell){
    kern_sl = function(x0, x1){
        return(kern(x0, x1, sig, ell))
    }
    x = xsobv
    y = ysobv
    n = length(x)
    K = outer(x, x, kern_sl)
    r = -1/2*t(y)%*%solve(K)%*%y - 1/2*log(det(K)) - n/2*log(pi*2)
    return(r)
}

## Parameterize and calculate log likelihoods
G = 50 # resolution of gridsearch
sigbs = c(0.5, 1.5) # bounds on parameters
ellbs = c(0.01, 1.5) # bounds on parameters
sigs = log(seq(exp(sigbs[1]), exp(sigbs[2]), length.out=G))
ells = log(seq(exp(ellbs[1]), exp(ellbs[2]), length.out=G))
params = expand.grid(sigs, ells) # generate all pairs of sigs and ells
names(params) = c("sig", "ell")

params$loglik = apply(params, 1, (function(r)(return(loglik(r[1], r[2])))))
params$loglik = log(-params$loglik)


## @knitr plot_parameters
# Plot the grid search for parameters
ix = which(params$loglik == min(params$loglik))[1]
sigopt = params$sig[ix]
ellopt = params$ell[ix]
plt_params = ggplot(data=params, aes(x=sig, y=ell, z=loglik, fill=loglik)) +
    geom_contour(aes(color=loglik), binwidth=0.05, color="steelblue") +
    geom_hline(aes(yintercept=ellopt), color="tomato") +
    geom_vline(aes(xintercept=sigopt), color="darkviolet") +
    labs(x=TeX("$\\sigma_k$"),
         y="l",
         title="Grid search for optimal kernel parameters")

plt_params


## @knitr draw_tuned_post
# Parameters
m = 50 # resolution of the posterior
n = 10 # number of draws from the posterior
d = t[2] - t[1] # length of the sampling support

# Collocated posterior points
xus = seq(t[1] - d/10, t[2] + d/10, length.out=m)

# Draw for the emperical 90% band and calculate the band
kern_sel_opt = function(x0, x1){
    return(sigopt^2 * exp(-1/2 * (x0-x1)^2 / ellopt))
}
posterior_opt = gen_posterior(xus, xsobv, ysobv, ker=kern_sel_opt)
draws = posterior_opt(n)
qs = apply(t(draws), 1, (function(x) (return(quantile(x, c(0.05, 0.95))))))

## @knitr plot_tuned_post
# Plot the aforementioned illustrative plot of draws from the posterior
df_post = data.frame(t(rbind(draws, xus)))
names(df_post) = c(paste("Draw", 1:n), "X")
df_post_m = melt(df_post, id="X")
plt_post_tuned = ggplot(df_post_m) +
    geom_line(aes(x=X, y=value, color=variable)) +
    geom_point(data=data.frame(x=xsobv, y=ysobv),
               aes(x=x, y=y),
               color="coral") +
    labs(x="X", y="Y", title=paste(n, "draws from the optimized posterior"),
         color="Draw number")


# Plot the 90% band
n = 50 # number of draws from the posterior
draws = posterior_opt(n)
qs = apply(t(draws), 1, (function(x) (return(quantile(x, c(0.05, 0.95))))))
df_post_band = data.frame(t(rbind(qs, xus)))
names(df_post_band) = c("q0", "q1", "X")
plt_post_band_tuned = ggplot(df_post_band) +
    geom_ribbon(aes(ymin=q0, ymax=q1, x=X), alpha=0.6) +
    geom_point(data=data.frame(x=xsobv, y=ysobv),
               aes(x=x, y=y),
               color="coral") +
    geom_line(data=data.frame(x=xs, y=ys), aes(x=x, y=y), color="steelblue") +
    labs(x="X", y="Y", title=paste("90% credible band using", n, "draws"))

plt_post_tuned
plt_post_band_tuned


## @knitr gp2d
t = c(-pi, pi)
D = 50 # D^2 number of grid points
ndraws = 6 # draw 6 times from the gaussian process
x1s = x2s = seq(t[1], t[2], length.out=D) # marginal collocation grid points
gd = expand.grid(x1s, x2s) # collocated grid
dm = dist(gd)
dmx = as.matrix(dm) # distance between each point on the collocated grid

mu2d = function(pts){
    # pts in R^(n x p) -> Ys in R^n
    # and p = 2 here i.e. we're working in 2D

    # ys = apply(pts, 1, prod) # arbitrary mean function
    ys = rep(0, dim(pts)[1]) # for simplicity use uniform 0 mean here
    return(ys)
}

mus = mu2d(gd)

kernel_on_dist = function(r){
    # r is the distance between points already computed by e.g. dist()
    ell = 1
    a = 1/sqrt(pi^2) * exp(-1/2*r/ell) # exponential kernel wrt Malhab dist
    # note that here the Malhabanobis distance is wrt the identity matrix
    # because dist() supports natively Euclidean distance and... well, a more
    # general distance function can be implemented later ;)
    return(a)
}

S = kernel_on_dist(dmx)

# NOTE: <S> is n by n, <gd> is n by p
stopifnot(dim(S)[1] == dim(S)[2])
stopifnot(dim(S)[1] == dim(gd)[1])

ys = t(mvrnorm(ndraws, mu = mus, Sigma=S))

## @knitr plot_gp2d
plot_2d_df = melt(
                  data.frame(
                             Y=ys,
                             X1=gd[,1],
                             X2=gd[,2]
                             ),
                  id=c("X1", "X2")
                  )

plt_2d = ggplot(plot_2d_df, aes(x=X1, y=X2)) +
    geom_raster(aes(fill=value), interpolate=TRUE) +
    facet_wrap(~variable) +
    scale_fill_gradientn(colours = heat.colors(20))

plt_2d
