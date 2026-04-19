# Estimation de la variance de $\hat{f}_{\hat{h}}(x)$ — Guide d'implémentation

## Contexte et objectif

On dispose d'un $n$-échantillon $X_1,\ldots,X_n \overset{\text{iid}}{\sim} f$ dans $\mathbb{R}^d$, où $f \in \mathcal{H}(s,L)$ est
potentiellement très irrégulière ($s$ faible). L'estimateur à noyau est

$$\hat{f}_h(x) = \frac{1}{n}\sum_{i=1}^n Z_i^h(x), \qquad Z_i^h(x) = \frac{1}{h^d}K\!\left(\frac{x - X_i}{h}\right)$$

avec $K$ noyau symétrique d'ordre 2 (typiquement gaussien), et $\hat{h}$ sélectionné par validation croisée
LOO.

**Objectif :** estimer ponctuellement

$$\sigma(x) = \sqrt{\text{Var}(\hat{f}_{\hat{h}}(x))}$$

et tracer la bande $\hat{f}_{\hat{h}}(x) \pm \hat{\sigma}(x)$ — visualisation de la variabilité stochastique,
sans prétention d'intervalle de confiance.

---

## Principe : estimateur par LGN

À $h$ **fixé**, les $Z_i^h(x)$ sont i.i.d. et

$$\text{Var}(\hat{f}_h(x)) = \frac{1}{n}\,\text{Var}(Z_1^h(x)) = \frac{1}{n}\!\left[\mathbb{E}\bigl[(Z_1^h(x))^2\bigr] - \mathbb{E}\bigl[Z_1^h(x)\bigr]^2\right]$$

Les deux termes sont des espérances estimables directement par LGN, **sans approximation sur la régularité
de $f$** — contrairement au plug-in $\hat{f}(x)\|K\|_2^2/(nh^d)$ qui suppose $f$ quasi-constante à l'échelle $h$.

### Problème avec $h = \hat{h}$

Avec $\hat{h} = \hat{h}(X_1,\ldots,X_n)$, les $Z_i^{\hat{h}}(x)$ **ne sont plus i.i.d.** : ils partagent
$\hat{h}$ qui dépend de l'échantillon entier. L'estimateur LGN naïf n'estime alors que la variance
conditionnelle $\text{Var}(\hat{f}_{\hat{h}}(x)\mid\hat{h})$, et manque le terme

$$\text{Var}\!\bigl(\mathbb{E}[\hat{f}_{\hat{h}}(x)\mid\hat{h}]\bigr) \asymp h_{\text{opt}}^{2(s-1)}\,\text{Var}(\hat{h})$$

qui est non négligeable quand $s < 1$ dans les zones irrégulières.

### Solution : pseudo-valeurs LOO

On pose $\hat{h}_{-i} = \hat{h}(X_1,\ldots,X_{i-1},X_{i+1},\ldots,X_n)$ et

$$\tilde{Z}_i(x) = \frac{1}{\hat{h}_{-i}^d}\,K\!\left(\frac{x - X_i}{\hat{h}_{-i}}\right)$$

Les $\tilde{Z}_i(x)$ sont **genuinement i.i.d.** : $\hat{h}_{-i}$ est calculé sans $X_i$, donc $\hat{h}_{-i}
\perp X_i$. L'estimateur

$$\boxed{\hat{\sigma}^2(x) = \frac{1}{n(n-1)}\sum_{i=1}^n \bigl(\tilde{Z}_i(x) - \bar{\tilde{Z}}(x)\bigr)^2,
\qquad \bar{\tilde{Z}}(x) = \frac{1}{n}\sum_{i=1}^n \tilde{Z}_i(x)}$$

est consistant pour $\text{Var}(\hat{f}_{\hat{h}}(x))$ en totalité (variance conditionnelle + variabilité de
la sélection), par LGN standard, sans aucune hypothèse supplémentaire sur $f$.

---

## Approximation de $\hat{h}_{-i}$ par la méthode des M-estimateurs

Calculer exactement les $n$ sélecteurs LOO $\hat{h}_{-i}$ coûterait $O(n^3)$. On exploite la structure
M-estimateur du sélecteur.

### Critère LOO et condition du premier ordre

Le critère LOO s'écrit

$$\text{CV}(h) = \int \hat{f}_h^2(x)\,dx - \frac{2}{n}\sum_{i=1}^n \hat{f}_{h,-i}(X_i)$$

$\hat{h}$ annule sa dérivée. En notant la contribution individuelle de $X_i$ :

$$\psi(h, X_i) = \frac{d}{dh}\bigl[\text{contribution de } X_i \text{ dans } \text{CV}(h)\bigr]$$

on a $\frac{1}{n}\sum_i \psi(\hat{h}, X_i) = 0$.

### Développement de Taylor au premier ordre

En écrivant la condition d'annulation sur l'échantillon sans $X_i$ et en développant autour de $\hat{h}$ :

$$0 = \frac{1}{n-1}\sum_{j\neq i}\psi(\hat{h}_{-i}, X_j) \approx \frac{1}{n-1}\sum_{j\neq i}\psi(\hat{h},
X_j) + (\hat{h}_{-i} - \hat{h})\,\widehat{H}$$

Puisque $\sum_j \psi(\hat{h}, X_j) = 0$, on a $\sum_{j\neq i}\psi(\hat{h}, X_j) = -\psi(\hat{h}, X_i)$.
On obtient :

$$\boxed{\hat{h}_{-i} \approx \hat{h} + \frac{\psi(\hat{h}, X_i)}{(n-1)\,\widehat{H}}}$$

où $\widehat{H} = \frac{d^2}{dh^2}\text{CV}(\hat{h})$ est la courbure du critère LOO au minimum.

### Précision de l'approximation

L'erreur de Taylor est $O\bigl((\hat{h}_{-i} - \hat{h})^2\bigr) = O(n^{-2})$. L'erreur induite sur
$\hat{\sigma}^2(x)$ est $O(n^{-2}h^{-d-1})$, négligeable devant la variance de l'estimateur lui-même
qui est $O(n^{-1})$.

---

## Algorithme complet

### Étape 1 — Sélection de bandwidth par LOO : $O(n^2)$

1. Évaluer $\text{CV}(h)$ sur une grille (ou par optimisation continue).
2. À $\hat{h}$ : calculer et **stocker** pour chaque $i$ la contribution individuelle
   $\psi(\hat{h}, X_i)$ à la dérivée du critère.
3. Calculer $\widehat{H} = \frac{d^2}{dh^2}\text{CV}(\hat{h})$ (courbure au minimum, disponible
   gratuitement si l'optimisation évalue aussi la dérivée seconde).

### Étape 2 — Approximation des $\hat{h}_{-i}$ : $O(n)$

Pour $i = 1,\ldots,n$ :

$$\hat{h}_{-i} \leftarrow \hat{h} + \frac{\psi(\hat{h}, X_i)}{(n-1)\,\widehat{H}}$$

> **Contrôle :** vérifier que $\hat{h}_{-i} > 0$ pour tout $i$. Si quelques $\hat{h}_{-i}$ sont négatifs
> ou aberrants (critère LOO mal conditionné), recalculer exactement ces quelques cas.

### Étape 3 — Calcul des $\tilde{Z}_i(x)$ sur grille : $O(nG)$

Sur une grille $(x_k)_{k=1}^G$ :

$$\tilde{Z}_i(x_k) = \frac{1}{\hat{h}_{-i}^d}\,K\!\left(\frac{x_k - X_i}{\hat{h}_{-i}}\right)$$

### Étape 4 — Estimateur de variance : $O(nG)$

$$\bar{\tilde{Z}}(x_k) = \frac{1}{n}\sum_{i=1}^n \tilde{Z}_i(x_k)$$

$$\hat{\sigma}^2(x_k) = \frac{1}{n(n-1)}\sum_{i=1}^n \bigl(\tilde{Z}_i(x_k) - \bar{\tilde{Z}}(x_k)\bigr)^2$$

La bande à tracer est $\hat{f}_{\hat{h}}(x_k) \pm \hat{\sigma}(x_k)$.

> **Remarque :** $\bar{\tilde{Z}}(x_k) \neq \hat{f}_{\hat{h}}(x_k)$ en général (les $\hat{h}_{-i}$
> diffèrent de $\hat{h}$), mais la différence est $O(n^{-1})$.

### Coût total

| Étape | Coût |
|---|---|
| Sélection LOO + stockage de $\psi_i$ | $O(n^2)$ |
| Approximation des $\hat{h}_{-i}$ | $O(n)$ |
| Calcul de $\tilde{Z}_i(x_k)$ | $O(nG)$ |
| Estimateur $\hat{\sigma}^2(x_k)$ | $O(nG)$ |
| **Total** | $O(n^2 + nG)$ |

Le coût est dominé par la sélection LOO, identique à ce qu'on ferait sans estimer la variance.

---

## Variante : sous-échantillonnage pour $n$ grand

Si $n$ est grand et le calcul de $\psi(\hat{h}, X_i)$ coûteux, on peut tirer un sous-échantillon
$\mathcal{I} \subset \{1,\ldots,n\}$ de taille $m \ll n$ et poser

$$\hat{\sigma}^2_m(x) = \frac{1}{m(m-1)}\sum_{i\in\mathcal{I}}\bigl(\tilde{Z}_i(x) -
\bar{\tilde{Z}}_{\mathcal{I}}(x)\bigr)^2$$

Cet estimateur est consistant dès que $m \to \infty$, avec variance d'ordre $O(m^{-1})$.

---

## Notes d'implémentation

**Noyau gaussien :** $K(u) = (2\pi)^{-d/2}\exp(-\|u\|^2/2)$, $\|K\|_2^2 = (4\pi)^{-d/2}$.

**Vectorisation :** les étapes 3 et 4 se vectorisent naturellement en représentant
$(\tilde{Z}_i(x_k))_{i,k}$ comme une matrice $n \times G$, dont on calcule la variance colonne par colonne.

**Forme analytique de $\psi(h, X_i)$ (noyau gaussien, $d=1$) :**

$$\psi(h, X_i) = \frac{d}{dh}\text{CV}(h) = -\frac{2}{n^2 h^2}\sum_{j\neq i}\left[\phi_{0,2h^2}(X_i - X_j)
+ \frac{(X_i-X_j)^2 - 2h^2}{2h^2}\phi_{0,2h^2}(X_i-X_j)\right] + \ldots$$

En pratique, on différencie numériquement le critère LOO analytique si la dérivée exacte est fastidieuse.

**LOO analytique :** pour le noyau gaussien, $\int \hat{f}_h^2 = \frac{1}{n^2}\sum_{i,j}\phi_{0,2h^2}(X_i-X_j)$
et $\hat{f}_{h,-i}(X_i) = \frac{1}{(n-1)h}\sum_{j\neq i}K\!\left(\frac{X_i-X_j}{h}\right)$
sont tous deux calculables en $O(n^2)$.
