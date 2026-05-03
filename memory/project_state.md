---
name: Project state and roadmap
description: État du package après le refactoring unified-api et prochaines fonctionnalités prévues
type: project
---

État du package après le refactoring unified-api (mergé) : introduction de `pp_lp` comme objet d'entrée unifié, `density_lp()` et `cv_density_lp()` sont devenus des génériques S3. CI passe (0 erreur, 0 warning R CMD check).

Prochaine fonctionnalité prévue : sélecteur de bandwidth **Goldenshluger-Lepski**, nommé `gl_density_lp` ou similaire.

**Why:** Même design que `cv_density_lp` — renvoie un objet de sélection (bandwidth(s) optimal(aux) + diagnostics) sans estimer la densité finale. Cohérent avec l'approche de `sparr` (LSCV.density/LIK.density).

**How to apply:** Retourner un objet `gl_density_lp` (classe dédiée) avec `h_hat`, les scores/critères intermédiaires, et le `call`. L'utilisateur appelle ensuite `density_lp(x, h = gl$h_hat, ...)`.
