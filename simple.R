## Preamble and imports
rm(list=ls())
seed = 33
set.seed(seed)
library(ggplot2)

## Parameters
t0 = -5
t1 = 5
N = 100 # plotting resolution
M = 10 # sampling density (M \le N)
sig = 0.2 # y = f(x) + e where e ~ N(0, sig)

## Define the model f(x)
mod = function(x){
    return(sin(x)*exp(-(x^2))*x^4)
}

## Generate some observations
xs_plot = seq(t0, t1, length.out=N)
ys_plot_no_noise = mod(xs_plot)
ys_plot = ys_plot_no_noise + rnorm(N, mean=0, sd=sig)
ixs = sample(1:N, M, replace=FALSE)
xs = xs_plot[ixs]
ys = ys_plot[ixs]

## Plot the true function and the noisy observations
plt_obs = ggplot(data.frame(x=xs_plot, yo=ys_plot, yt=ys_plot_no_noise)) +
    geom_line(aes(x=x, y=yt), color="steelblue", size=1.1) +
    geom_point(aes(x=x, y=yo), color="coral", shape=18, size=1) +
    labs(x="X", y="Y", title="Model and generated observations")
