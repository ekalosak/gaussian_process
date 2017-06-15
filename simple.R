## Preamble and imports
rm(list=ls())
seed = 33
set.seed(seed)
library(ggplot2)
library(MASS)
library(reshape)

## Parameters
t0 = -5
t1 = 5
N = 100 # plotting resolution
M = 30 # sampling density (M \le N)
sig = 0.15 # y = f(x) + e where e ~ N(0, sig)

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
    geom_point(aes(x=x, y=yo), color="coral", size=1) +
    labs(x="X", y="Y", title="Model and generated observations")


# src: https://www.cs.toronto.edu/~hinton/csc2515/notes/gp_slides_fall08.pdf
kernel = function(x1, x2){
    return(exp(-abs(x1-x2)))
}
S = outer(xs, xs, kernel) # covariance matrix computation
# means are the observed values

## Plot a handful of draws from this process
K = 17 # number of draws
gpdraws = mvrnorm(K, ys, S)
df_gp_ex = data.frame(t(gpdraws))
names(df_gp_ex) = paste(rep("GP", K), 1:K, sep="")
df_gp_ex$X = xs
df_gp_ex_m = melt(df_gp_ex, id=c("X"))
df_gp_ex_qs = data.frame(
    t(apply(gpdraws, 2, function(x){quantile(x, c(0.1, 0.9))}))
)
names(df_gp_ex_qs) = c("q0.1", "q0.9")
df_gp_ex_qs$x = xs
plt_ex_gp = ggplot(data=df_gp_ex) +
    geom_line(data=df_gp_ex_m, aes(x=X, y=value, color=variable)) +
    geom_point(data=data.frame(x=xs, y=ys),
               aes(x=x, y=y), color="coral", size=1.1) +
    geom_line(data=data.frame(x=xs_plot, y=ys_plot_no_noise),
              aes(x=x, y=y), color="steelblue", size=1.2) +
    geom_ribbon(data=df_gp_ex_qs,
                aes(ymin=q0.1, ymax=q0.9, x=x), alpha=0.4, color="grey")
