# LOO Cross-Validation and Variance Estimation for the LP Estimator

## Context

This document derives the leave-one-out (LOO) cross-validation score and the
variance estimator for the local polynomial density estimator implemented in
`densityLP`. It serves as the mathematical basis for the functions
`cv_density_lp` and the variance output of `density_lp_point`.

---

## 1. Additive structure of the estimator

Fix a target point $t \in \mathcal{D}$, a multi-index set $\gamma = (m, h)$,
and let $B_\gamma$ be the Gram matrix estimated by Monte Carlo quadrature.
Write $L_\gamma$ for its lower Cholesky factor ($B_\gamma = L_\gamma
L_\gamma^\top$) and set

$$H_\gamma(u) = B_\gamma^{-1/2}\,\Phi_\gamma(u/h), \qquad
H_0 = H_\gamma(0) = L_\gamma^{-\top} e_1,$$

where $e_1$ is the first standard basis vector (the constant monomial evaluates
to 1 at $u = 0$, all others to 0).

The local polynomial estimator (equation (4) of the article) is

$$\hat{f}_\gamma(t)
  = \frac{1}{n}\sum_{i=1}^n Z_i^\gamma(t),
  \qquad
  Z_i^\gamma(t) = h^{-d}\,H_0^\top H_\gamma(X_i - t)\,\mathbf{1}_{X_i \in V_t(h)},$$

where $V_t(h) = \mathcal{D} \cap (t + [-h,h]^d)$.  The terms
$Z_i^\gamma(t)$ are i.i.d. (the Gram matrix $B_\gamma$ is computed from
quadrature independent of the observations).

**In the code.** In `lp_estimator_cpp`, the quantity $H_0^\top H_\gamma(X_i-t)$
is the dot product of `H_0` with column $i$ of `H_V`. The final estimator is

```
(h^{-d} / n) * dot(H_0, rowSums(H_V))   =   (1/n) * sum_i Z_i
```

---

## 2. Variance estimator

Since the $Z_i$ are i.i.d.,

$$\operatorname{Var}(\hat f_\gamma(t))
  = \frac{1}{n}\operatorname{Var}(Z_1^\gamma(t))
  \approx \frac{1}{n}\,\mathbb{E}[(Z_1^\gamma(t))^2].$$

The natural empirical estimator (equation (5) of the article) is

$$\hat{v}_\gamma(t) = \frac{1}{n^2}\sum_{i=1}^n (Z_i^\gamma(t))^2.$$

**Implementation.** The individual $Z_i$ are available inside `lp_estimator_cpp`
as the column-wise dot products `H_0^T H_V[:, i]` (times $h^{-d}$). Computing
$\hat v_\gamma$ costs a single additional pass over these values — essentially
free.

---

## 3. LOO correction

**Key formula.** The observation $X_i$ contributes to $\hat f_\gamma(X_i)$ via
the self-influence term

$$Z_i^\gamma(X_i) = h^{-d}\,H_0^{(X_i)\top} H_\gamma^{(X_i)}(0)
  = h^{-d}\,e_1^\top B_{X_i}^{-1} e_1
  = h^{-d}\,(B_{X_i}^{-1})_{11},$$

where the superscript $(X_i)$ indicates that $B_\gamma$ is computed with $t =
X_i$. Since the Gram matrix is independent of the observations, removing $X_i$
from the training set only removes its contribution from the sum:

$$\hat f_{\gamma,-i}(X_i)
  = \frac{n\,\hat f_\gamma(X_i) - Z_i^\gamma(X_i)}{n - 1}.$$

**No additional Cholesky decomposition is needed.** The term
$(B_{X_i}^{-1})_{11}$ equals $\|L_{X_i}^{-\top} e_1\|^2 = \|H_0^{(X_i)}\|^2$,
which is the squared norm of `H_0` — already computed during the estimation at
$t = X_i$.

**Implementation plan.**

1. For each $i$, call `density_lp_point` with $t = X_i$, and request that it
   also returns $\|H_0\|^2 = (B_{X_i}^{-1})_{11}$.
2. The LOO estimate follows algebraically:

$$\hat f_{\gamma,-i}(X_i) = \frac{n\,\hat f_\gamma(X_i) - h^{-d}\|H_0^{(X_i)}\|^2}{n-1}.$$

**Cost.** $n$ evaluations of `density_lp_point` at the data points, each with
$N_\text{quad}$ quadrature draws. No extra Cholesky or triangular solve beyond
what the standard estimation already requires.

---

## 4. LOO CV score

For a pair $\gamma = (m, h)$, the LOO log-likelihood score is

$$\operatorname{CV}(\gamma)
  = -\frac{1}{n}\sum_{i=1}^n \log\hat f_{\gamma,-i}(X_i),$$

with the convention $\log(0) = -\infty$ (zero estimates are penalised
automatically). The selected pair is

$$\hat\gamma = \arg\min_{\gamma \in \Gamma} \operatorname{CV}(\gamma).$$

**Grid.** In practice, $\Gamma = \mathcal{M} \times \mathcal{H}$ where
$\mathcal{M} \subseteq \{0, 1, 2, 3\}$ and $\mathcal{H}$ is a user-supplied
bandwidth grid. The score is evaluated at each point of the grid; the
minimising pair $\hat\gamma$ is returned along with the full score matrix for
diagnostic purposes.

---

## 5. API sketch

```r
cv_density_lp(X, h_grid, m_grid = 0:3, domain, N_quad = 500L)
# Returns: list(gamma_hat, h_hat, m_hat, scores = matrix |m_grid| x |h_grid|)

cv_density_lp_ppp(pp, h_grid, m_grid = 0:3, N_quad = 500L)
# Returns: same, domain derived from Window(pp)
```

The score matrix allows the user to plot the CV surface and inspect the
selection.

---

## 6. Connection to the article's selection procedure

The article uses the Goldenshluger-Lepski (GL) selection rule (Section 2.3),
which minimises $\hat A_\gamma + \hat U_\gamma$ — a bias upper bound plus a
penalised standard deviation bound. The GL procedure is theoretically optimal
(adaptive minimax) but requires tuning the constant $\delta > 1$.

LOO CV is a practical alternative that:
- requires no tuning constant,
- targets predictive performance rather than pointwise risk,
- is consistent for bandwidth selection under mild conditions.

Both procedures rely on the same $\hat v_\gamma$ variance estimator. The GL
selection will be implemented separately as `select_lp` for theoretical
reproducibility of the article's simulations.
