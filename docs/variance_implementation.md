# Estimation de la Variance dans densityLP

Ce document décrit l'implémentation actuelle de l'estimateur de variance ponctuelle dans le package `densityLP`, cohérente avec l'approche par polynômes locaux sur domaines complexes (Bertin, Klutchnikoff & Ouimet).

## Modèle et Noyau

`densityLP` utilise exclusivement un **noyau "box"** (uniforme) sur l'hypercube $[-h, h]^d$. Ce choix est motivé par :
1. La robustesse théorique sur des domaines $\mathcal{D}$ aux frontières irrégulières.
2. L'efficacité du calcul de la matrice de Gram via échantillonnage de Monte Carlo (rejection sampling).

## Formulation de l'Estimateur

L'estimateur polynomial local de degré $m$ au point $t$ s'écrit comme une moyenne de contributions individuelles des observations $X_i$ situées dans le voisinage $V(h) = \mathcal{D} \cap [t-h, t+h]^d$ :

$$\hat{f}(t) = \frac{1}{n} \sum_{i=1}^n Z_i$$

où la contribution $Z_i$ est calculée via le vecteur de base $\Phi$ et la matrice de Gram $B_\gamma$ :
- Si $X_i \notin V(h)$, $Z_i = 0$.
- Si $X_i \in V(h)$, $Z_i = h^{-d} e_1^T B_\gamma^{-1} \Phi((X_i - t)/h)$.

## Calcul de la Variance Ponctuelle

À bandwidth $h$ fixé, la variance de l'estimateur est estimée par la variance empirique des contributions $Z_i$ :

$$\widehat{\text{Var}}(\hat{f}(t)) = \frac{1}{n^2} \sum_{i=1}^n Z_i^2$$

*(Note : On utilise ici la somme des carrés sans soustraire la moyenne au carré pour rester conservateur et cohérent avec l'implémentation C++ optimisée).*

### Détails de l'implémentation C++

Dans le backend `src/gram_matrix.cpp`, le calcul est effectué efficacement lors de l'évaluation de la densité :
1. Calcul de $H_0 = L^{-1} e_1$ (via la décomposition de Cholesky de la matrice de Gram).
2. Calcul des scores $S = H_V^T H_0$ pour toutes les observations dans le voisinage.
3. $\text{variance} = \|S\|^2 / (n^2 h^{2d})$.

## Utilisation des résultats

L'objet retourné par `density_lp()` contient un vecteur `$variance` de même taille que `$estimate`. Ces valeurs peuvent être utilisées pour construire des bandes de variabilité stochastique :

$$\hat{f}(t) \pm \sqrt{\widehat{\text{Var}}(\hat{f}(t))}$$

> **Attention :** Cette variance ne prend pas en compte l'incertitude liée à la sélection automatique de $\hat{h}$ par validation croisée (voir `docs/research_h_variability.md` pour des pistes théoriques à ce sujet).
