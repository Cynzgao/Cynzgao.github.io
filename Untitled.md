
---
title: 'Blog Post number 1'
date: 2019-06-14
permalink: /posts/2019/06/blog-post-1/
tags:
  - cool posts
  - numerical methods
---

## Introduction
First developed by Harten, Engquist, Osher, and Chakravarthy in 1987, the Essentially Non-Oscillatory (ENO) methods are part of high-resolution schemes in finding the numerical solution of partial differential equations. Since hyperbolic PDEs may contain singularities in their solutions, the ENO methods are desirable because they are able to provide high-order accuracy when the function is smooth while avoiding Gibbs phenomenon at discontinuities. The ENO schemes could achieve such goals by reconstructing a piece-wise polynomial of the solution from its cell average. In other words, the procedure uses the smoothest stencil to approximate fluxes at cell boundaries to avoid oscillations near shocks.

As researchers later pointed out, the ENO methods, although uniformly high order accurate, can have certain drawbacks. Its reconstruction procedure from cell averages could be complicated and computationally expensive. In addition, its adaptive stencil could change from perturbations and might not be in the regions where solution is smooth. In 1994, a new version of ENO schemes was proposed by Liu, Osher, and Chan. This new modification is called weighted ENO (WENO). The main idea is that, unlike ENO that chooses the smoothest stencil for reconstruction, WENO will use a convex combination of all candidate stencils. A weight will be assigned to each candidate, which represents how much this stencil contributes to the final solution. The scheme will remain essentially non-oscillatory as ENO, while adding one more order of accuracy.

Two years later, an improved version of WENO method was proposed by Jiang and Shu. The previous WENO scheme requires weights on candidate stencils, which depends on their relative smoothness. Thus, following the previous WENO scheme, Jiang and Shu devised a new way of evaluating the smoothness of a stencil. This new smoothness indicator achieves the highest possible order of accuracy (fifth order) when $r = 3$. We will now present the derivations of two WENO schemes.

## WENO, version 1
### Problem set-up
Recall the hyperbolic conservation law $\bf{u}_t + div \bf{f}(\bf{u}) = 0$. In 1 dimension, the conservation law becomes
$$u_t + f(u)_x = 0.$$

For a uniform discretization in space $\cdots < x_0 < x_1 < x_2 < \cdots$, $\Delta x = x_{i+1} - x_i$. We use $I_i = [x_{i - 1/2}, x_{i+1/2}]$ to represent i-th cells. The cell average can be written as

$$\bar{u}_i =\frac{1}{\Delta x} \int_{I_i}u(x, t)dx
= \frac{1}{\Delta x} \int^{x_{i+1/2}}_{x_{i-1/2}}u(x, t)dx. $$

Integrating the conservation law over each cell, we obtain

$$\frac{d}{dt}\bar{u}_i(t) = - \frac{1}{\Delta x}(f(u(x_{i+1/2}, t)) - f(u(x_{i-1/2}, t))).$$

In order to evaluate $f(u(x, t))$ at each cell boundaries $x_{i+1/2}$, we can reconstruct a piecewise polynomial $R(x) = \{R_j(x)\}$ approximating $u(x, t)$ given the cell averages $\bar{u} = \{\bar{u}_j\}$. For a general flux, there are multiple methods used for flux-splitting, including Engquist-Osher, Godunov, Roe with entropy fix, and Lax-Friedrichs etc.

Putting the time discretization in a more compact form, we can write
$$
\begin{align*}
    \frac{d}{dt}\bar{u}_i(t) &\approx L_i(\bar{u})\\
    L_i(\bar{u}) &= -\frac{1}{\Delta x}[\hat{f}(R_i(x_{i+1/2}), R_{i+1}(x_{i+1/2}) - \hat{f}(R_{i-1}(x_{x_{i-1/2}}), R_i(x_{i-1/2}))],
\end{align*}
$$
where $\hat{f}(R_i(x_{i+1/2}), R_{i+1}(x_{i+1/2})$ approximates $f(u(x_{i+1/2}, t))$, $\hat{f}(R_{i-1}(x_{x_{i-1/2}}), R_i(x_{i-1/2}))$ approximates $f(u(x_{i-1/2}, t))$.


### Reconstruction
We denote the r candidate stencils by $S_k$, $k = 0, 1, \cdots, r-1$, where

$$S_k = \{I_{i+k-r+1}, I_{i+k-r+2}, \cdots, I_{i+k}\}$$

For example, when $r = 3$, the candidate stencils would be
$$
\begin{align*}
    S_0 &= \{I_{i - 2}, I_{i-1}, I_i\} = \{x_{i-5/2}, x_{i-3/2}, x_{i-1/2}, x_{i+1/2}\}\\
    S_1 &= \{I_{i - 1}, I_{i}, I_{i+1}\} = \{x_{i-3/2}, x_{i-1/2}, x_{i+1/2}, x_{i+3/2}\}\\
    S_2 &= \{I_{i}, I_{i+1}, I_{i+2}\} = \{x_{i-1/2}, x_{i+1/2}, x_{i+3/2}, x_{i+5/2}\}\\
\end{align*}
$$
The three candidate stencils are illustrated below. Each color represents a different stencil.


We would then obtain the interpolation polynomials on each of the stencils such that

$$
\bar{u}_i = \frac{1}{\Delta x_i} \int_{I_i} p_j (x)dx.
$$

Again, using $r = 3$ as an example, we have three corresponding quadratic polynomials,
$$
\begin{align*}
    p_i(x) &= \frac{\bar{u}_i - 2\bar{u}_{i-1} + \bar{u}_{i-2}}{2 \Delta x^2} (x - x_{i-1})^2 + \frac{\bar{u}_i - \bar{u}_{i-2}}{2 \Delta x}(x - x_{i-1}) + \bar{u}_{i-1} - \frac{\bar{u}_i - 2\bar{u}_{i-1} + \bar{u}_{i-2}}{24}\\
    p_{i+1}(x) &= \frac{\bar{u}_{i+1} - 2\bar{u}_{i} + \bar{u}_{i-1}}{2 \Delta x^2} (x - x_{i})^2 + \frac{\bar{u}_{i+1} - \bar{u}_{i-1}}{2 \Delta x}(x - x_{i}) + \bar{u}_{i} - \frac{\bar{u}_{i+1} - 2\bar{u}_{i} + \bar{u}_{i-1}}{24}\\
    p_{i+2}(x) &= \frac{\bar{u}_{i+2} - 2\bar{u}_{i+1} + \bar{u}_{i}}{2 \Delta x^2} (x - x_{i+1})^2 + \frac{\bar{u}_{i+2} - \bar{u}_{i}}{2 \Delta x}(x - x_{i+1}) + \bar{u}_{i+1} - \frac{\bar{u}_{i+2} - 2\bar{u}_{i+1} + \bar{u}_{i}}{24}.\\
\end{align*}
$$
Let the approximated value $u^{k}_{i+1/2} = p_k(x_{i+1/2})$, the expressions above can be simplified to
$$
\begin{align*}
    u^0_{i+1/2} &= \frac{1}{3} \bar{u}_{i-2} - \frac{7}{6}\bar{u}_{i-1} + \frac{11}{6}\bar{u}_{i}\\
    u^1_{i+1/2} &= -\frac{1}{6} \bar{u}_{i-1} + \frac{5}{6}\bar{u}_{i} + \frac{1}{3}\bar{u}_{i+1}\\
    u^3_{i+1/2} &= \frac{1}{3} \bar{u}_{i} + \frac{5}{6}\bar{u}_{i+1} - \frac{1}{6}\bar{u}_{i+2}\\
\end{align*}
$$
Note that this approximation is third order accurate, namely, $u^k_{i+1/2} - u(x_{i+1/2}) = \mathcal{O}(\Delta x^3)$. However, we are not satisfied with interpolation on a single stencil, we want to achieve more by using all of them together. In other words, we aim to find a convex combination of the polynomials that follows the conservation law and it is essentially non-oscillatory. Since the hyperbolic PDEs may contain shocks and discontinuities, we should use the smoother part of the region to form our solution. The smoother the stencil is, the heavier weight it will have on the combination. The difference between two versions of the WENO methods is largely based on how the smoothness indicator is defined.

### First Smoothness Measurement
The indicator of smoothness (IS, or $\beta$'s in some literatures) on $j$-th stencil proposed by Liu, Osher, and Chan is defined to be
$$
IS_j = \sum^{r-1}_{l=1} \left ( \sum^l_{k=1} (\triangle^{r-1}[u_{j-r+k}])^2 \right)/l,
$$

where $\triangle[u]$ comes from a table of differences of $\{\bar{u}_i\}$ on $S_i$ and satisfies $\triangle[u_l] = u_{l+1} - u_l, \triangle^k[u_l] = \triangle^{k-1}[u_{l+1}] - \triangle^{k-1}[u_l]$. This way, the smoothness indicator becomes the summation of all averages of squared same order differences. When $r = 3$,
$$IS_j= \frac{(\triangle[u_{j-2}])^2 + (\triangle[u_{j-1}])^2}{2} + (\triangle^2[u_{j-2}])^2.$$
When the true solution $u(x,t)$ is discontinuous on $S_j$, $IS_j \approx \mathcal{O}(1)$; if $u(x, t)$ is continuous on $S_j$, $IS_j \approx \mathcal{O}(\Delta x^2)$.

The convex combination $R_j(x)$ that comes from $r$ interpolating polynomials ($p_k$ obtained before) on $r$ stencils is of the form
$$
\begin{align*}
    R_j(x) &= \sum^{r-1}_{k=0} \omega_k p_{j+k}(x)\\
    \omega_k &= \frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l},
\end{align*}
$$
where $\sum^{r-1}_{k=0} \omega_k = 1, k = 0, 1, \cdots, r-1$ and $\alpha^j_k >0$. The ENO property is satisfied if the coefficients $\frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} = \mathcal{O}(1)$ if the stencil $S_{j+k}$ is in the smooth regions, and $\frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} \leq \mathcal{O}(\Delta x^r)$ if the stencil is in the discontinuous region of the solution.

Let us consider the coefficients defined as
$$
\alpha^j_k = \frac{C^j_k}{(\epsilon+ IS_{j+k})^r},
$$
where $\epsilon = (10^{-5} \sim 10^{-7})$ is a small positive number added here to ensure that the denominator will not become zero because $IS$ could potentially be zero. The coefficients $C_k^r$'s can be found by computing the weights of polynomials on each stencil to form the polynomial across all candidate stencils, ie.
$$
Q(x_{i+1/2}) = \sum^{r-1}_{k=0}C_k^r p_k(x_{i+1/2}),
$$
where $Q(x)$ is the interpolating polynomial on the largest stencil.

For example, when $r = 3$, we have 3 candidate stencils and the largest stencil that includes all contains 5 cells $\{I_{i-2}, I_{i-1}, I_{i}, I_{i+1}, I_{i+2}\}$. $Q(x)$, evaluated at $x_{i+1/2}$, has the expression
$$
Q(x_{i+1/2}) = \frac{1}{30} \bar{u}_{i-2} - \frac{13}{60}\bar{u}_{i-1} + \frac{47}{60} \bar{u}_i + \frac{9}{20} \bar{u}_{i+1} - \frac{1}{20}\bar{u}_{i+2}
$$

Combined with previous expressions for interpolation polynomials $p_k$'s we found before, it is easy to see that the coefficients are given by $C_0^3 = \frac{1}{10}, C_1^3 = \frac{6}{10}, C_2^3 = \frac{3}{10}$. When $r = 2$, $C_0^2 = \frac{1}{3}, C_1^2 = \frac{2}{3}$. The optimal weight for $C_k^r$ leads to one order of improvement in accuracy.  Now it is clear to see that
$$
\begin{align*}
    \frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} &= \mathcal{O}(1), \text{when stencil is smooth} \\
    \frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} &\leq \max(\mathcal{O}(\epsilon^r), \mathcal{O}((\Delta x)^r)), \text{when stencil meets discontinuities}.
\end{align*}
$$
Thus, these weights $\{\alpha^j_k\}^{r-1}_{k=0}$ satisfies the ENO properties.

Going back to our example when $r = 3$, the reconstructed solution is expressed as
$$
R_j (x) = \omega_0 p_j(x) + \omega_1 p_{j+1}(x) + \omega_2 p_{j+2}(x),
$$
where $p_k(x)$'s are calculated above and the weights are given by
$$
\omega_0 = \frac{\alpha_0^j}{\alpha_0^j + \alpha_1^j + \alpha_2^j}, \omega_1 = \frac{\alpha_1^j}{\alpha_0^j + \alpha_1^j + \alpha_2^j},
\omega_2 = \frac{\alpha_2^j}{\alpha_0^j + \alpha_1^j + \alpha_2^j}
$$
and
$$
\alpha_0^j = \frac{1}{10(\epsilon + IS_j)^3},
\alpha_1^j = \frac{6}{10(\epsilon + IS_{j+1})^3},
\alpha_2^j = \frac{3}{10(\epsilon + IS_{j+2})^3}.
$$

The interpolating method described above satisfies
$$
\begin{align*}
    &u(x, t) = R_j(x) + \mathcal{O}((\Delta x)^r)\\
    &u(x_j^*, t) = R_j(x_j^*) + \mathcal{O}((\Delta x)^{r+1}),
\end{align*}
$$
where $x_j^*$'s are chosen to be points at the cell boundaries. For a general upwind method, when $f'(R(x)) > 0$,
$$
L_j(\bar{u}) = -\frac{1}{\Delta x}[f(R_j(x_{j+1/2})) - f(R_{j-1}(x_{j-1/2}))].
$$
The points chosen are $x_j^* = x_{j+1/2}, x_{j-1}^* = x_{j - 1/2}$. Following the reconstruction scheme presented above, we obtain
$$
\begin{align*}
    &R_i(x_{j+1/2}) - u(x_{j+1/2}, t) = \mathcal{O}(\Delta x^{r+1})\\
    &R_{j-1}(x_{j-1/2}) - u(x_{j-1/2}, t) = \mathcal{O}(\Delta x^{r+1}).
\end{align*}
$$
Hence,
$$
L_j(\bar{u}) = -\frac{1}{\Delta x}[f(u(x_{j+1/2}, t)) - f(u(x_{j-1/2}, t))] + \mathcal{O}(\Delta x^{r+1}).
$$
Similarly, for $f'(R(x)) < 0$, we have
$$
L_j(\bar{u}) = -\frac{1}{\Delta x}[f(R_{j+1}(x_{j+1/2})) - f(R_{j}(x_{j-1/2}))].
$$
Choosing $x_j^* = x_{j-1/2}, x_{j+1}^* = x_{j+1/2}$, we achieve the same order of accuracy. Therefore, in teh smooth regions, the spatial discretization $L$ is able to approximate $\frac{d \bar{u}}{dt}$ to the $r+1$-th order.
