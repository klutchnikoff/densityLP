# Design document — Package R `densityLP`

> Architecture pour l'implémentation de l'estimateur de densité par polynômes
> locaux sur domaines complexes (Bertin, Klutchnikoff, Ouimet — *Bernoulli*, 2026).

---

## État du Projet (Avril 2026)

Le package est fonctionnel et dispose d'un backend C++ performant (RcppArmadillo). La structure S3 est en place pour les domaines et les estimations. La transition Quarto → RMarkdown est terminée.

---

## 1. Architecture Générale (Réelle)

L'implémentation actuelle suit une séparation stricte entre la géométrie (R) et le calcul numérique (C++) :

- **R (Interface & Géométrie) :**
  - Objets `domain_lp` : Portent la logique de la frontière via des samplers.
  - `density_lp` : Orchestre le calcul sur une grille (boucle R).
  - `density_lp_ppp` : Intégration avec l'écosystème `spatstat`.
- **C++ (Calcul Intensif) :**
  - `lp_estimator_loo_cpp` : Calcule ponctuellement l'estimation, la variance et le terme d'influence LOO en une seule passe.
  - `cv_lp_fixed_h_m_cpp` : Version vectorisée pour la cross-validation (accélération ×50).

---

## 2. Feuille de route mise à jour

### Phase 1 — Socle (Terminée ✅)
- [x] `build_alphas()` & `build_Phi()`
- [x] Backend C++ pour la matrice de Gram et Cholesky
- [x] Estimation de la variance ponctuelle $\widehat{\text{Var}}(\hat{f}(t))$
- [x] Systèmes de domaines (`Rd`, `indicator`, `sector`, `owin`)
- [x] Méthodes S3 `print` et `plot` de base

### Phase 2 — Optimisation et Robustesse (En cours 🚧)
- [x] Migration de la boucle de Cross-Validation (LOO) vers le C++.
- [x] Migration de Quarto vers RMarkdown pour la compatibilité CRAN.
- [ ] **Batch Processing :** Passer du calcul point par point à un calcul par blocs en C++ pour réduire l'overhead `.Call`.
- [ ] **Parallélisation :** Intégration native de `future.apply` dans `density_lp` pour les grilles haute résolution.
- [ ] **API Cleanup :** Retourner des objets (listes nommées) plutôt que des vecteurs pour les résultats internes.

### Phase 3 — Sélection Adaptative (Prochaine étape 🎯)
- [ ] **Goldenshluger-Lepski :**
  - [ ] Calcul du terme de biais $\hat{A}_\gamma$ par comparaison de paires d'estimateurs.
  - [ ] Grille de bandwidths $\Gamma_n$ sur échelle logarithmique.
  - [ ] Implémentation de `density_lp_adaptive()`.
- [ ] **Validation :** Reproduire les simulations de l'article (Phase 5 du plan initial).

### Phase 4 — Maturité & Alignement `sparr` (Ambition 🚀)
- [ ] **Interface :** Support complet des objets `im` et `ppp` de `spatstat` (marques, fenêtres multiples).
- [ ] **Visualisation :** Développer des fonctions de plot avancées (intervalles de confiance, surfaces 3D).
- [ ] **Noyaux :** Ajouter le support optionnel de noyaux à support compact (Epanechnikov) pour un meilleur lissage visuel.
- [ ] **Documentation :** Créer une "Galerie des domaines complexes" montrant la force du package.

---

## 3. Notes Techniques sur la Variance

La variance est calculée à $h$ fixé par la variance empirique des contributions individuelles $Z_i$ au point $t$.
$$ \widehat{\text{Var}}(\hat{f}(t)) = \frac{1}{n^2} \sum_{i=1}^n Z_i^2 $$
Le backend C++ fournit cette valeur "gratuitement" lors du calcul de l'estimateur. Pour plus de détails, voir `docs/variance_implementation.md`.

---

## 4. Organisation du Code

- `R/` : Logique métier, interface S3, et wrappers.
- `src/` : Backend Armadillo (Gram matrix, LOO correction, variance).
- `vignettes/` : Guides d'utilisation au format RMarkdown.
- `docs/` : Notes théoriques et roadmap.

---
*Dernière mise à jour : 27 Avril 2026.*
