# gaussian_process

Some snippets for practicing gaussian process models

A Gaussian Process is a continuum of random variables $X_t$ for
$t\in[\tau_0, \tau_1]\subset \mathbf{R}$ where every finite subset of random variables
$\{X_1,\dots,X_k\} \subset X_t$ has a multivaraite gaussian joint distribution
with a correlation matrix $\Sigma$ typically generated by a distance-sensitive
kernel e.g. $\Sigma(x,x^\star) = x^\top x^\star$.
