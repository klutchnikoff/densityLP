# Design document — Package R `densityLP`

> Architecture pour l'implémentation de l'estimateur de densité par polynômes
> locaux sur domaines complexes (Bertin, Klutchnikoff, Ouimet — *Bernoulli*, à
> paraître).

---

## Table des matières

1. [Contexte statistique](#1-contexte-statistique)
2. [Architecture générale](#2-architecture-générale)
3. [Cœur algorithmique générique](#3-cœur-algorithmique-générique)
4. [Samplers — générateurs uniformes sur V(h)](#4-samplers--générateurs-uniformes-sur-vh)
5. [Wrappers géospatiaux (d = 2)](#5-wrappers-géospatiaux-d--2)
6. [Intégration Rcpp / RcppArmadillo](#6-intégration-rcpp--rcpparmadillo)
7. [Structure du package](#7-structure-du-package)
8. [Dépendances](#8-dépendances)
9. [Feuille de route](#9-feuille-de-route)

---

## 1. Contexte statistique

### L'estimateur de base

Soit $X_1, \ldots, X_n$ un échantillon i.i.d. de densité $f$ sur un domaine
$\mathcal{D} \subset \mathbb{R}^d$. Pour $\gamma = (m, h) \in \mathbb{N}_0
\times (0, \rho]$ et un point cible $t \in \mathcal{D}$, l'estimateur
polynomial local est :

$$\hat{f}_\gamma(t)
  = \frac{1}{n} \sum_{i=1}^n H_\gamma^\top(0)\, H_\gamma(X_i - t)\, w_h(X_i - t)$$

**Ingrédients :**

| Objet | Définition | Rôle |
|---|---|---|
| $\mathcal{V}(h)$ | $(\mathcal{D} - t) \cap h[-1,1]^d$ | Voisinage de $t$ |
| $K(u)$ | $\mathbf{1}_{[-1,1]^d}(u)$ | Noyau uniforme |
| $w_h(u)$ | $h^{-d} K(h^{-1}u)\,\mathbf{1}_\mathcal{D}(t+u)$ | Poids restreint au domaine |
| $\Phi_\gamma(u)$ | $\{\varphi_\alpha(u/h)\}_{|\alpha| \leq m}$ | Monomiales rescalées ($D_m \times 1$) |
| $\mathcal{B}_\gamma$ | $\int_{\mathcal{V}(h)} \Phi_\gamma \Phi_\gamma^\top w_h\,\mathrm{d}u$ | Matrice de Gram ($D_m \times D_m$) |
| $H_\gamma(u)$ | $\mathcal{B}_\gamma^{-1/2} \Phi_\gamma(u)$ | Base orthonormée (Cholesky) |

avec $D_m = \binom{m+d}{d}$ et les multi-indices $\alpha \in \mathbb{N}_0^d$
ordonnés par degré total puis lexicographiquement.

### Réécriture Monte Carlo

Si $u_1, \ldots, u_N \sim \mathrm{Unif}(\mathcal{V}(h))$, la matrice de Gram
s'estime par :

$$\hat{\mathcal{B}}_\gamma
  = \frac{\mathrm{vol}(\mathcal{V}(h))}{N \cdot h^d}
    \sum_{k=1}^N \Phi_\gamma(u_k)\,\Phi_\gamma^\top(u_k)$$

**Conséquence architecturale :** toute la complexité géométrique est absorbée
par un générateur uniforme sur $\mathcal{V}(h)$. Le volume
$\mathrm{vol}(\mathcal{V}(h))$ est obtenu gratuitement par ce même générateur.

---

## 2. Architecture générale

```
┌──────────────────────────────────────────────────────────────────┐
│                   Wrapper géospatial (d = 2)                     │
│                                                                  │
│              density_lp_ppp(pp, h, m, nx, ny, N_quad)           │
│              ppp/owin → grille automatique depuis Window(pp)     │
│                          │                                       │
│                 → "density_lp_ppp" (S3)                         │
└──────────────────────────┬───────────────────────────────────────┘
                           │
┌──────────────────────────▼───────────────────────────────────────┐
│               Interface utilisateur principale                   │
│                                                                  │
│   density_lp(X, t_grid, h, m, domain, N_quad)                   │
│   X      : matrice n × d   (données brutes, pas de conversion)   │
│   t_grid : matrice p × d   (points cibles)                       │
│   domain : objet "lp_domain" (constructeur dédié)               │
│   → objet S3 "density_lp"                                        │
│                                                                  │
│   Boucle sur t_grid → appels à density_lp_point() (interne)     │
└──────────────────────────┬───────────────────────────────────────┘
                           │
┌──────────────────────────▼───────────────────────────────────────┐
│                  Cœur (interne, non exporté)                     │
│                                                                  │
│   density_lp_point(X, t, h, m, domain, N_quad)  → scalaire      │
│                                                                  │
│   ┌─────────────┐   ┌────────────────┐   ┌──────────────────┐   │
│   │build_alphas │   │ gram_matrix()  │   │ lp_estimator()   │   │
│   │(multi-idx)  │   │ (MC interne)   │   │ (somme pondérée) │   │
│   └─────────────┘   └────────────────┘   └──────────────────┘   │
└──────────────────────────┬───────────────────────────────────────┘
                           │
┌──────────────────────────▼───────────────────────────────────────┐
│          Objets domaine "lp_domain" (constructeurs S3)           │
│                                                                  │
│  domain_Rd(d)          → R^d (pas de contrainte)                 │
│  domain_from_indicator(f, d)     → domaine via is_in_domain : n×d→logical  │
│  domain_sector(k)      → secteur D_k = {y ≤ x^k, x ∈ [0,1]}    │
│                                                                  │
│  Chaque objet porte : d, sampler interne, label                  │
│  Samplers : entièrement internes, jamais exposés                 │
└──────────────────────────────────────────────────────────────────┘
```

**Principe de séparation des responsabilités :**

- `density_lp` est l'unique fonction exportée du cœur. Elle travaille sur
  une **grille** de points cibles et retourne un objet S3.
- `density_lp_point` est **interne** : un point cible → un scalaire.
  Cela rend la parallélisation future triviale (remplacer la boucle par
  `future.apply::future_lapply` sans changer l'interface).
- Le **domaine** est un objet S3 `"lp_domain"` créé par un constructeur.
  Le sampler est encapsulé dedans et invisible pour l'utilisateur.
- `X` est toujours une **matrice numérique** : pas de conversion implicite
  depuis tibble ou data.frame. La dimension `ncol(X) == domain$d` est
  vérifiée explicitement avec un message d'erreur clair.
- Le wrapper `density_lp_ppp` est **séparé** : il construit la grille
  depuis `Window(pp)` et retourne un objet S3 `"density_lp_ppp"`.

---

## 3. Cœur algorithmique générique

### 3.1 Multi-indices

```r
# Génère tous les alpha ∈ N_0^d avec |alpha| ≤ m
# Ordonnés par degré total, puis lexicographiquement
build_alphas <- function(m, d) {
  grille  <- do.call(expand.grid, rep(list(0:m), d))
  valides <- rowSums(grille) <= m
  alphas  <- as.matrix(grille[valides, , drop = FALSE])
  # Tri : degré total croissant, puis lexicographique
  ord <- order(rowSums(alphas),
               apply(alphas, 1, paste, collapse = "-"))
  alphas[ord, , drop = FALSE]
}
# D_m = nrow(build_alphas(m, d)) = choose(m + d, d)
```

### 3.2 Monomiales rescalées

```r
# U    : matrice d × N  (points u_k)
# m, h : degré, bandwidth
# Retourne : matrice D_m × N
build_Phi <- function(U, m, h, alphas = NULL) {
  d      <- nrow(U)
  if (is.null(alphas)) alphas <- build_alphas(m, d)
  Dm     <- nrow(alphas)
  N      <- ncol(U)
  U_sc   <- U / h   # d × N  (u/h)

  # Évaluation vectorisée : φ_α(u/h) = ∏_j (u_j/h)^{α_j}
  Phi <- matrix(NA_real_, Dm, N)
  for (j in seq_len(Dm)) {
    Phi[j, ] <- apply(U_sc ^ alphas[j, ], 2, prod)
  }
  Phi
}
# Note : pour d=2, m petit, une version explicite sans boucle est préférable.
```

### 3.3 Matrice de Gram

```r
# s : sortie d'un sampler — list(points = d×N, vol = scalaire)
gram_matrix <- function(s, m, h, alphas = NULL) {
  d   <- nrow(s$points)
  N   <- ncol(s$points)
  if (is.null(alphas)) alphas <- build_alphas(m, d)

  Phi <- build_Phi(s$points, m, h, alphas)     # Dm × N
  w   <- s$vol / (N * h^d)                     # poids MC uniforme
  w * tcrossprod(Phi)                           # Dm × Dm
}
```

### 3.4 Estimateur principal

```r
lp_density_core <- function(X,          # n × d
                             t,          # d-vecteur (point cible)
                             h,          # bandwidth > 0
                             m,          # degré polynomial ≥ 0
                             sampler,    # function(N, t, h) → list(points, vol)
                             N_quad = 500) {
  d      <- ncol(X)
  n      <- nrow(X)
  alphas <- build_alphas(m, d)

  # ── 1. Matrice de Gram via sampler ────────────────────────────────────────
  s <- sampler(N_quad, t, h)
  if (ncol(s$points) < 2)
    stop("Voisinage V(h) trop petit ou vide pour h = ", h)

  B <- gram_matrix(s, m, h, alphas)

  # ── 2. Cholesky → base orthonormée ───────────────────────────────────────
  L     <- chol(B)                           # upper triangular : B = L^T L
  L_inv <- backsolve(L, diag(nrow(L)))       # B^{-1/2} = L^{-1}

  # ── 3. H_γ(0) ────────────────────────────────────────────────────────────
  # Φ_γ(0) = e_1 (seul le monôme constant vaut 1 en u=0)
  Phi_0 <- matrix(0, nrow(alphas), 1)
  Phi_0[1, 1] <- 1.0
  H_0 <- L_inv %*% Phi_0                    # Dm × 1

  # ── 4. Observations dans V(h) ────────────────────────────────────────────
  U    <- sweep(X, 2, t)                     # n × d  (X_i - t)
  in_V <- apply(abs(U) <= h, 1, all)         # noyau uniforme sur [-h,h]^d
  U_V  <- t(U[in_V, , drop = FALSE])         # d × N_V

  if (ncol(U_V) == 0) return(0)

  # ── 5. Estimateur ─────────────────────────────────────────────────────────
  Phi_V <- build_Phi(U_V, m, h, alphas)      # Dm × N_V
  H_V   <- L_inv %*% Phi_V                   # Dm × N_V
  (h^(-d) / n) * sum(crossprod(H_0, H_V))
}
```

---

## 4. Samplers — générateurs uniformes sur V(h)

### Interface contractuelle

Tout sampler est une **fonction-fabrique** retournant une **closure** :

```r
# Appel : sampler(N, t, h)
# Retour : list(
#   points = matrice d × N  (coordonnées u = x - t),
#   vol    = scalaire        (vol(V(h)) estimé ou exact)
# )
```

Le volume `vol` est fourni gratuitement par le sampler — c'est l'un des
avantages majeurs de cette architecture.

---

### 4.1 `sampler_rejection` — fallback universel

Valide pour tout domaine représentable par `is_in_domain : ℝ^d → {TRUE,FALSE}`.

```r
sampler_rejection <- function(is_in_domain, n_estimate = 2000) {
  function(N, t, h) {
    d <- length(t)

    # Estimation préliminaire du taux d'acceptation
    U_pre    <- matrix(runif(n_estimate * d, -h, h), n_estimate, d)
    in_V_pre <- is_in_domain(sweep(U_pre, 2, t, "+"))
    p_acc    <- mean(in_V_pre)
    vol      <- (2 * h)^d * p_acc

    # Rejection sampling
    accepted <- matrix(NA_real_, d, 0)
    while (ncol(accepted) < N) {
      need <- ceiling(1.5 * (N - ncol(accepted)) / max(p_acc, 0.01))
      U    <- matrix(runif(need * d, -h, h), need, d)
      in_V <- is_in_domain(sweep(U, 2, t, "+"))
      accepted <- cbind(accepted, t(U[in_V, , drop = FALSE]))
    }
    list(points = accepted[, seq_len(N), drop = FALSE], vol = vol)
  }
}
```

---

### 4.2 `sampler_spatstat` — exact pour d = 2

Exploite `spatstat.geom::runifpoint` pour un tirage exact sur le polygone
$\mathcal{V}(h) = \mathcal{D} \cap ([t_1-h, t_1+h] \times [t_2-h, t_2+h])$.

```r
sampler_spatstat <- function(win) {
  # win : objet owin représentant D (coordonnées métriques projetées)
  function(N, t, h) {
    box  <- spatstat.geom::owin(c(t[1]-h, t[1]+h),
                                c(t[2]-h, t[2]+h))
    V_h  <- spatstat.geom::intersect.owin(win, box)
    vol  <- spatstat.geom::area.owin(V_h)        # exact (polygone)
    pts  <- spatstat.geom::runifpoint(N, win = V_h)
    list(
      points = rbind(pts$x - t[1], pts$y - t[2]),  # recentré en t
      vol    = vol
    )
  }
}
```

---

### 4.3 `sampler_qmc` — quasi-Monte Carlo (Sobol)

Convergence en $O((\log N)^d / N)$ même pour des intégrandes discontinus.
Recommandé pour $d \geq 2$ avec `is_in_domain` quelconque.

```r
sampler_qmc <- function(is_in_domain, scrambling = 1) {
  function(N, t, h) {
    d       <- length(t)
    N_sobol <- ceiling(3 * N)
    sob     <- randtoolbox::sobol(N_sobol, dim = d, scrambling = scrambling)
    U       <- sweep(sob * 2*h, 2, h, "-")           # N_sobol × d

    in_V    <- is_in_domain(sweep(U, 2, t, "+"))
    vol     <- (2*h)^d * mean(in_V)
    pts     <- t(U[in_V, , drop = FALSE])             # d × N_V

    # Complétion par rejection si taux d'acceptation trop faible
    if (ncol(pts) < N) {
      warning("QMC : taux d'acceptation faible (",
              round(mean(in_V)*100, 1), "%), complétion par MC standard.")
      while (ncol(pts) < N) {
        U2   <- matrix(runif(N * d, -h, h), N, d)
        in_V2 <- is_in_domain(sweep(U2, 2, t, "+"))
        pts  <- cbind(pts, t(U2[in_V2, , drop = FALSE]))
      }
    }
    list(points = pts[, seq_len(N), drop = FALSE], vol = vol)
  }
}
```

---

### 4.4 `sampler_analytic` — secteur polynomial $\mathcal{D}_k$

Pour les domaines $\mathcal{D}_k = \{(x,y) : 0 \leq y \leq x^k,\, 0 \leq x
\leq 1\}$ utilisés dans les simulations de l'article.

```r
sampler_sector <- function(k) {
  function(N, t, h) {
    tx <- t[1]; ty <- t[2]
    x_lo <- max(0,  tx - h)
    x_hi <- min(1,  tx + h)

    # Volume exact par intégration 1D
    vol <- integrate(function(x) {
      pmax(0, pmin(x^k, ty + h) - pmax(0, ty - h))
    }, x_lo, x_hi)$value

    # Rejection sampling (efficace : domaine monotone en x)
    accepted <- matrix(NA_real_, 2, 0)
    while (ncol(accepted) < N) {
      x  <- runif(2*N, x_lo, x_hi)
      y  <- runif(2*N, max(0, ty - h), ty + h)
      ok <- (y >= 0) & (y <= x^k) & (abs(x - tx) <= h)
      accepted <- cbind(accepted, rbind(x[ok] - tx, y[ok] - ty))
    }
    list(points = accepted[, seq_len(N), drop = FALSE], vol = vol)
  }
}
```

---

### 4.5 Tableau de synthèse des samplers

| Sampler | Dimension | Volume | Précision | Dépendances |
|---|---|---|---|---|
| `sampler_spatstat` | d = 2 | exact | exact | `spatstat.geom` |
| `sampler_qmc` | d ≥ 2 | MC | quasi-exacte | `randtoolbox` |
| `sampler_rejection` | d ≥ 1 | MC | MC standard | — |
| `sampler_sector` | d = 2 | exact (1D) | rejection | — |
| `sampler_har` | d >> 2 | MC | MC (hit-and-run) | à implémenter |

---

## 5. Wrappers géospatiaux (d = 2)

### 5.1 Principes de projection

Travailler en **coordonnées métriques projetées** (Lambert-93 EPSG:2154, ou
UTM) est indispensable pour que le bandwidth $h$ ait une interprétation
physique en mètres. La reprojection se fait en amont via `sf`.

```r
# Systèmes de coordonnées recommandés :
# - France métropolitaine : EPSG:2154 (Lambert-93)
# - Usage général         : UTM (choisir le fuseau adapté)
# - Monde entier          : EPSG:3857 (Web Mercator, à éviter près des pôles)
```

### 5.2 Wrapper `lp_density_ppp`

```r
lp_density_ppp <- function(pp,           # objet ppp (coordonnées déjà projetées)
                            t,            # c(x, y) en unités de projection
                            h,            # bandwidth (mètres si projection métrique)
                            m,            # degré polynomial
                            sampler_type = c("spatstat", "qmc", "rejection"),
                            N_quad = 500,
                            ...) {
  sampler_type <- match.arg(sampler_type)
  win <- spatstat.geom::Window(pp)
  X   <- cbind(pp$x, pp$y)

  sampler <- switch(sampler_type,
    spatstat  = sampler_spatstat(win),
    qmc       = sampler_qmc(owin_indicator(win)),
    rejection = sampler_rejection(owin_indicator(win))
  )
  lp_density_core(X, t, h, m, sampler, N_quad, ...)
}

# Indicatrice depuis un owin — évaluation vectorisée
owin_indicator <- function(win) {
  function(X_mat) {
    spatstat.geom::inside.owin(x = X_mat[, 1],
                               y = X_mat[, 2],
                               w = win)
  }
}
```

### 5.3 Wrapper `lp_density_sf`

```r
lp_density_sf <- function(sf_points,     # objet sf de type POINT
                           sf_domain,    # objet sf de type POLYGON / MULTIPOLYGON
                           t_lonlat,     # c(lon, lat) — point cible (WGS84)
                           h,            # bandwidth en mètres
                           m,
                           crs_metric = 2154,   # EPSG Lambert-93 par défaut
                           sampler_type = "spatstat",
                           N_quad = 500, ...) {
  crs_m    <- sf::st_crs(crs_metric)

  pts_m    <- sf::st_transform(sf_points, crs_m)
  dom_m    <- sf::st_transform(sf_domain, crs_m)
  t_sf     <- sf::st_transform(
                 sf::st_sfc(sf::st_point(t_lonlat), crs = 4326),
                 crs_m)
  t_m      <- as.numeric(sf::st_coordinates(t_sf))

  X_mat    <- sf::st_coordinates(pts_m)[, 1:2, drop = FALSE]

  # Conversion sf → owin pour spatstat
  win      <- spatstat.geom::as.owin(sf::as_Spatial(dom_m))
  pp       <- spatstat.geom::as.ppp(X_mat, W = win)

  lp_density_ppp(pp, t_m, h, m, sampler_type, N_quad, ...)
}
```

### 5.4 Surface de densité sur une grille

```r
lp_density_surface <- function(pp, h, m,
                                nx = 100, ny = 100,
                                sampler_type = "spatstat",
                                N_quad = 300, ...) {
  win  <- spatstat.geom::Window(pp)
  bbox <- spatstat.geom::boundingbox(win)
  xs   <- seq(bbox$xrange[1], bbox$xrange[2], length.out = nx)
  ys   <- seq(bbox$yrange[1], bbox$yrange[2], length.out = ny)

  grid <- as.matrix(expand.grid(x = xs, y = ys))
  # Restreindre aux points dans D
  in_D <- spatstat.geom::inside.owin(grid[,1], grid[,2], win)
  grid_D <- grid[in_D, , drop = FALSE]

  # Calcul (parallélisable)
  est <- apply(grid_D, 1, function(t_vec) {
    tryCatch(
      lp_density_ppp(pp, t_vec, h, m, sampler_type, N_quad, ...),
      error = function(e) NA_real_
    )
  })

  # Recomposition en image spatstat
  vals <- rep(NA_real_, nrow(grid))
  vals[in_D] <- est
  spatstat.geom::im(matrix(vals, nx, ny), xcol = xs, yrow = ys)
}
```

---

## 6. Intégration Rcpp / RcppArmadillo

### 6.1 Stratégie

La séparation entre R et C++ suit le principe :

- **R** : appel à `is_in_domain` (peut être n'importe quelle fonction R), 
  génération des points du sampler, gestion des erreurs.
- **C++** : accumulation de la matrice de Gram, décomposition de Cholesky,
  calcul de la somme pondérée — toutes les opérations sur données numériques
  denses.

```
Côté R                          Côté C++ (RcppArmadillo)
────────────────────────────    ─────────────────────────────────────
sampler(N, t, h)           →    gram_matrix_cpp(U_V, alphas, h, w)
is_in_domain(X)            →    chol(), inv(), somme pondérée
filtrage in_V              →    lp_estimator_cpp(U_V, U_Q, alphas, h, vol, n)
```

### 6.2 Code C++ principal

Fichier : `src/lp_core.cpp`

```cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ── Évaluation de Φ_γ(u/h) pour un vecteur u ───────────────────────────────
// alphas : Dm × d (entiers)
// u      : d-vecteur
inline arma::vec phi_eval(const arma::vec& u, const arma::imat& alphas,
                           double h) {
  int Dm = alphas.n_rows;
  int d  = alphas.n_cols;
  arma::vec phi(Dm);
  arma::vec uh = u / h;
  for (int j = 0; j < Dm; ++j) {
    double val = 1.0;
    for (int l = 0; l < d; ++l)
      val *= std::pow(uh(l), alphas(j, l));
    phi(j) = val;
  }
  return phi;
}

// ── Matrice de Gram ─────────────────────────────────────────────────────────
// U_Q    : d × N_quad  (points de quadrature, recentrés en t)
// alphas : Dm × d
// h      : bandwidth
// weight : vol(V(h)) / (N_quad * h^d)
// [[Rcpp::export]]
arma::mat gram_matrix_cpp(arma::mat U_Q, arma::imat alphas,
                           double h, double weight) {
  int Dm = alphas.n_rows;
  int N  = U_Q.n_cols;
  arma::mat B(Dm, Dm, arma::fill::zeros);

  for (int k = 0; k < N; ++k) {
    arma::vec phi = phi_eval(U_Q.col(k), alphas, h);
    B += phi * phi.t();
  }
  return weight * B;
}

// ── Estimateur complet (une seule passe C++) ────────────────────────────────
// U_Q    : d × N_quad  (points de quadrature dans V(h), recentrés)
// U_V    : d × N_V     (X_i - t pour X_i dans V(h))
// alphas : Dm × d
// h, vol, n : paramètres
// [[Rcpp::export]]
double lp_estimator_cpp(arma::mat U_Q, arma::mat U_V, arma::imat alphas,
                         double h, double vol, int n) {
  int Dm   = alphas.n_rows;
  int N_Q  = U_Q.n_cols;
  int N_V  = U_V.n_cols;
  int d    = U_Q.n_rows;

  // 1. Gram matrix
  double weight = vol / (N_Q * std::pow(h, d));
  arma::mat B   = gram_matrix_cpp(U_Q, alphas, h, weight);

  // 2. Cholesky : B = L^T L (upper), puis L^{-1}
  arma::mat L     = arma::chol(B, "upper");
  arma::mat L_inv = arma::inv(arma::trimatu(L));

  // 3. H_γ(0) = L_inv * e_1  (monôme constant = 1 en u=0)
  arma::vec e1(Dm, arma::fill::zeros);
  e1(0) = 1.0;
  arma::vec H0 = L_inv * e1;

  // 4. Somme : Σ_i H_0^T H(u_i) w_h(u_i)
  double result = 0.0;
  double w_obs  = std::pow(h, -d) / n;

  for (int k = 0; k < N_V; ++k) {
    arma::vec phi = phi_eval(U_V.col(k), alphas, h);
    arma::vec Hk  = L_inv * phi;
    result       += arma::dot(H0, Hk);
  }
  return w_obs * result;
}
```

### 6.3 Wrapper R utilisant le backend C++

```r
lp_density_core_rcpp <- function(X, t, h, m, sampler, N_quad = 500) {
  d      <- ncol(X)
  n      <- nrow(X)
  alphas <- build_alphas(m, d)
  # Conversion entière pour Rcpp
  alphas_int <- structure(as.integer(alphas),
                          dim = dim(alphas))

  # Sampler → points de quadrature (côté R)
  s   <- sampler(N_quad, t, h)
  U_Q <- s$points          # d × N_quad
  vol <- s$vol

  # Filtrage observations (côté R — is_in_domain peut être owin etc.)
  U    <- sweep(X, 2, t)
  in_V <- apply(abs(U) <= h, 1, all)
  U_V  <- t(U[in_V, , drop = FALSE])   # d × N_V

  if (ncol(U_V) == 0) return(0)

  # Délégation à C++
  lp_estimator_cpp(U_Q, U_V, alphas_int, h, vol, n)
}
```

---

## 7. Structure du package

```
densityLP/
├── DESCRIPTION
├── NAMESPACE
├── README.md
│
├── R/
│   ├── alphas.R          # build_alphas(), build_Phi()           [interne]
│   ├── gram.R            # gram_matrix()                         [interne]
│   ├── samplers.R        # sampler_rejection(), sampler_qmc(),   [interne]
│   │                     # sampler_spatstat(), sampler_sector()
│   ├── domain.R          # domain_Rd(), domain_from_indicator(),           [exporté]
│   │                     # domain_sector() — constructeurs S3 "lp_domain"
│   ├── core.R            # density_lp_point()                    [interne]
│   ├── density_lp.R      # density_lp()                          [exporté]
│   ├── density_lp_ppp.R  # density_lp_ppp()                      [exporté]
│   ├── methods.R         # print/plot pour "density_lp",         [exporté]
│   │                     # "density_lp_ppp"
│   └── utils.R           # vérifications, messages d'erreur      [interne]
│
├── src/
│   ├── lp_core.cpp       # gram_matrix_cpp(), lp_estimator_cpp()
│   └── Makevars          # flags de compilation (OpenMP si parallèle)
│
├── tests/
│   └── testthat/
│       ├── test-alphas.R       # multi-indices et monomiales
│       ├── test-core.R         # estimateur sur domaine simple (R^d)
│       ├── test-samplers.R     # contrat vol et dimension
│       ├── test-gram.R         # symétrie, définie-positivité
│       ├── test-wrappers.R     # cohérence ppp ↔ sf
│       └── test-rcpp.R         # parité R / C++
│
├── vignettes/
│   ├── getting-started.Rmd    # exemple minimal 2D
│   ├── spatial-data.Rmd       # usage avec spatstat/sf, données GPS
│   └── simulations.Rmd        # reproduction des simulations de l'article
│
└── data/
    └── sector_domain.rda      # owin du secteur D_{2.1} pour les exemples
```

### DESCRIPTION (squelette)

```
Package: densityLP
Title: Local Polynomial Density Estimation on Complicated Domains
Version: 0.1.0
Authors@R: c(
    person("Karine",  "Bertin",      role = c("aut"), ...),
    person("Nicolas", "Klutchnikoff", role = c("aut", "cre"), ...),
    person("Frédéric","Ouimet",      role = c("aut"), ...))
Description: Nonparametric local polynomial estimator for multivariate
    density functions on known domains of arbitrary dimension, including
    complicated domains with sharp concavities and holes. Implements the
    procedure of Bertin, Klutchnikoff and Ouimet (Bernoulli, forthcoming).
License: GPL-3
Imports:
    spatstat.geom,
    sf,
    randtoolbox,
    cubature,
    Rcpp,
    RcppArmadillo
Suggests:
    testthat (>= 3.0.0),
    ggplot2,
    patchwork
LinkingTo: Rcpp, RcppArmadillo
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.x.x
```

---

## 8. Dépendances

| Package | Rôle | Obligatoire |
|---|---|---|
| `spatstat.geom` | `owin`, `ppp`, `runifpoint`, `intersect.owin` | Non (wrapper d=2) |
| `sf` | Reprojection CRS, conversion polygon | Non (wrapper GPS) |
| `randtoolbox` | Séquence de Sobol (QMC) | Non (sampler_qmc) |
| `cubature` | `hcubature` adaptatif | Non (alternative) |
| `Rcpp` | Interface C++ | Recommandé |
| `RcppArmadillo` | Algèbre linéaire dense (Cholesky, BLAS) | Recommandé |

Le cœur pur R (`lp_density_core`) n'a **aucune dépendance obligatoire**.

---

## 9. Feuille de route

### Phase 1 — Cœur R (priorité maximale)
- [ ] `build_alphas()` + tests unitaires (ordre, cardinalité $D_m$)
- [ ] `build_Phi()` + tests (valeurs exactes pour m=0,1,2 et d=1,2)
- [ ] `gram_matrix()` + tests (symétrie, définie-positivité)
- [ ] Constructeurs `domain_Rd()`, `domain_from_indicator()`, `domain_sector()`
- [ ] `density_lp_point()` avec `sampler_rejection` uniquement
- [ ] `density_lp()` (boucle sur grille → objet S3 `"density_lp"`)
- [ ] `print` / `plot` pour `"density_lp"`
- [ ] Validation sur $\mathcal{D} = [0,1]^2$ avec `domain_Rd()`

### Phase 2 — Samplers additionnels
- [ ] `sampler_spatstat()` + test de convergence du volume
- [ ] `sampler_qmc()` + comparaison variance vs rejection
- [ ] `sampler_sector()` pour $\mathcal{D}_1$ et $\mathcal{D}_{2.1}$
- [ ] Validation numérique de `domain_sector()`

### Phase 3 — Backend Rcpp
- [ ] `gram_matrix_cpp()` + parité numérique avec version R
- [ ] `lp_estimator_cpp()` (une passe complète)
- [ ] Benchmark R vs C++ (objectif : ×10 minimum)

### Phase 4 — Wrapper spatstat
- [ ] `density_lp_ppp()` : grille depuis `Window(pp)`, objet S3 `"density_lp_ppp"`
- [ ] `print` / `plot` pour `"density_lp_ppp"` (via `spatstat::im`)

### Phase 5 — Simulations de l'article
- [ ] Reproduction des 4 densités ($f_1, f_{2.1}, g_1, g_{2.1}$)
- [ ] Boucle sur $n \in \{200, 500, 1000, 2000\}$, $R = 5000$ réplications
- [ ] Comparaison oracle LP vs oracle `sparr`
- [ ] Parallélisation via `future.apply` ou `parallel`

### Phase 6 — Sélection adaptative (Goldenshluger-Lepski)
- [ ] Calcul de $\hat{v}_\gamma$, $\hat{A}_\gamma$, $\hat{\mathbb{U}}_\gamma$
- [ ] Grille $\Gamma_n$ sur échelle logarithmique
- [ ] `density_lp_adaptive()` — estimateur final sans paramètre de réglage

---

*Document généré lors de la conception du package, avril 2026.*
*Référence : Bertin K., Klutchnikoff N., Ouimet F. — "A new adaptive local*
*polynomial density estimation procedure on complicated domains", Bernoulli.*
