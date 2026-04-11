# Modèle ChloricemSCM — Pénétration des chlorures avec modèle de complexation de surface non-électrostatique

> **Fichier source :**
> `src/Models/ModelFiles/ChloricemSCM.cpp`
>
> **Exemple :**
> `base/ChloricemSCM/ChloricemSCM` · `base/ChloricemSCM/ChloricemSCM.msh`
>
> **Modèle parent :** Chloricem (isotherme de Langmuir)
>
> **Voir aussi :** ChloricemEla (modèle d'Elakneswaran)

---

## Table des matières

1. [Contexte et objectif](#1-contexte-et-objectif)
2. [Différence avec Chloricem](#2-différence-avec-chloricem)
3. [Modèle de complexation de surface NE-SCM](#3-modèle-de-complexation-de-surface-ne-scm)
   - 3.1 [Réactions de surface](#31-réactions-de-surface)
   - 3.2 [Formule d'adsorption](#32-formule-dadsorption)
   - 3.3 [Dépendance au pH](#33-dépendance-au-ph)
   - 3.4 [Lien avec l'isotherme de Langmuir](#34-lien-avec-lisotherme-de-langmuir)
4. [Paramètres du modèle SCM](#4-paramètres-du-modèle-scm)
5. [Fichier d'entrée](#5-fichier-dentrée)
6. [Calibration et valeurs indicatives](#6-calibration-et-valeurs-indicatives)
7. [Références bibliographiques](#7-références-bibliographiques)

---

## 1. Contexte et objectif

Le modèle **ChloricemSCM** est une variante du modèle Chloricem dans laquelle l'isotherme de Langmuir pour la fixation des chlorures par les C-S-H est remplacée par un **modèle de complexation de surface non-électrostatique (NE-SCM)**.

Toutes les autres équations — transport multi-ionique (Nernst-Planck), conservation des éléments (Ca, Si, Na, K, Cl), chimie de la solution (portlandite, C-S-H, calcite, sel de Friedel) — sont identiques à celles de Chloricem.

---

## 2. Différence avec Chloricem

Dans **Chloricem**, la quantité de chlorures adsorbés sur les C-S-H est donnée par une isotherme de Langmuir dont les paramètres $\alpha$ et $\beta$ dépendent du rapport Ca/Si ($x_\text{csh}$) :

$$\Gamma_\text{Cl}^\text{Langmuir}(c_\text{Cl}, x_\text{csh}) = \frac{\alpha(x_\text{csh})\,c_\text{Cl}}{1 + \beta(x_\text{csh})\,c_\text{Cl}}$$

Cette formule présente deux limites :

1. Elle ne tient **pas compte du pH** : la fixation des chlorures par les C-S-H est pourtant connue pour dépendre fortement du pH de la solution interstitielle.
2. Elle est **purement empirique** : les paramètres $\alpha$ et $\beta$ n'ont pas de sens mécaniste direct.

Dans **ChloricemSCM**, ces paramètres sont remplacés par des constantes d'équilibre de surface dont la signification physique est explicite.

---

## 3. Modèle de complexation de surface NE-SCM

### 3.1 Réactions de surface

Le modèle considère les **sites silanols** (≡SiOH) à la surface des C-S-H, avec trois réactions d'équilibre :

| Réaction | Constante |
|----------|-----------|
| $\equiv\text{SiOH} \rightleftharpoons \equiv\text{SiO}^- + \text{H}^+$ | $K_\text{acid}$ (mol/L) |
| $\equiv\text{SiOH} + \text{H}^+ \rightleftharpoons \equiv\text{SiOH}_2^+$ | $K_\text{bas}$ (L/mol) |
| $\equiv\text{SiOH}_2^+ + \text{Cl}^- \rightleftharpoons \equiv\text{SiOH}_2\text{Cl}$ | $K_\text{Cl}$ (L/mol) |

La première réaction (déprotonation) et la deuxième (protonation) décrivent l'amphotérisme des silanols. La troisième décrit la **complexation de surface de l'ion chlorure** sur le site protoné, le seul chargé positivement susceptible de fixer Cl⁻.

### 3.2 Formule d'adsorption

En écrivant le bilan de masse sur les sites de surface ($\Gamma_\text{max}$ = densité totale de sites, en mol par mol de C-S-H) :

$$\Gamma_\text{max} = [\equiv\text{SiOH}] + [\equiv\text{SiO}^-] + [\equiv\text{SiOH}_2^+] + [\equiv\text{SiOH}_2\text{Cl}]$$

et en exprimant chaque espèce par rapport à $[\equiv\text{SiOH}]$ via les lois d'action de masse, on obtient la quantité de chlorures adsorbés par mole de C-S-H :

$$\boxed{
\Gamma_\text{Cl}(c_\text{Cl}, c_\text{H}, x_\text{csh}) =
  \frac{\Gamma_\text{max}(x_\text{csh})\;K_\text{bas}\;K_\text{Cl}\;c_\text{H}\;c_\text{Cl}}
       {1 + \dfrac{K_\text{acid}}{c_\text{H}} + K_\text{bas}\,c_\text{H}\,(1 + K_\text{Cl}\,c_\text{Cl})}
}$$

où $c_\text{H} = [\text{H}^+]$ est calculé à chaque pas par le solveur `HardenedCementChemistry`.

### 3.3 Dépendance au pH

Le comportement limite selon le pH est le suivant :

- **pH très bas** ($c_\text{H} \gg K_\text{acid}/K_\text{bas}$, sites dominés par $\equiv\text{SiOH}_2^+$) :
  $$\Gamma_\text{Cl} \approx \frac{\Gamma_\text{max}\,K_\text{Cl}\,c_\text{Cl}}{1 + K_\text{Cl}\,c_\text{Cl}}$$
  isotherme de Langmuir pure avec capacité $\Gamma_\text{max}$.

- **pH intermédiaire** (autour du point de charge nulle, $c_\text{H} \approx \sqrt{K_\text{acid}/K_\text{bas}}$) :
  adsorption maximale.

- **pH très élevé** ($c_\text{H} \ll \sqrt{K_\text{acid}/K_\text{bas}}$, sites dominés par $\equiv\text{SiO}^-$) :
  $$\Gamma_\text{Cl} \approx \frac{\Gamma_\text{max}\,K_\text{bas}\,K_\text{Cl}\,c_\text{H}^2\,c_\text{Cl}}{K_\text{acid}}$$
  adsorption proportionnelle à $c_\text{H}^2$ — très sensible au pH (cohérent avec le fait que la surface très négative repousse Cl⁻).

### 3.4 Lien avec l'isotherme de Langmuir

À **pH fixé** (solution tampon), le NE-SCM se réduit à une isotherme de Langmuir effective :

$$\Gamma_\text{Cl} \approx \frac{\alpha_\text{eff}\,c_\text{Cl}}{1 + \beta_\text{eff}\,c_\text{Cl}}$$

avec :

$$\alpha_\text{eff}(c_\text{H}) = \frac{\Gamma_\text{max}\,K_\text{bas}\,K_\text{Cl}\,c_\text{H}}{1 + K_\text{acid}/c_\text{H} + K_\text{bas}\,c_\text{H}}$$

$$\beta_\text{eff}(c_\text{H}) = \frac{K_\text{bas}\,K_\text{Cl}\,c_\text{H}}{1 + K_\text{acid}/c_\text{H} + K_\text{bas}\,c_\text{H}}$$

Ainsi, les paramètres de Langmuir $\alpha$ et $\beta$ peuvent être interprétés comme des valeurs effectives à un pH de référence.

---

## 4. Paramètres du modèle SCM

Les quatre paramètres peuvent dépendre du rapport Ca/Si $x_\text{csh}$, et sont fournis sous forme de courbes dans le fichier d'entrée :

| Paramètre | Courbe Bil | Description | Unité |
|-----------|-----------|-------------|-------|
| $\Gamma_\text{max}$ | `gamma_max` | Densité totale de sites de surface | mol/mol$_\text{CSH}$ |
| $K_\text{acid}$ | `k_acid` | Constante d'acidité $K_{a1}$ (≡SiOH → ≡SiO⁻ + H⁺) | mol/L |
| $K_\text{bas}$ | `k_bas` | Constante de protonation $K_b$ (≡SiOH + H⁺ → ≡SiOH₂⁺) | L/mol |
| $K_\text{Cl}$ | `k_cl` | Constante de complexation du chlorure | L/mol |

---

## 5. Fichier d'entrée

La différence avec le fichier `Chloricem` tient uniquement à la ligne `Curves` d'adsorption des chlorures :

**Chloricem (Langmuir) :**
```
Curves = Adscl  x = Range{x1=0.5, x2=1.5, n=100}
         alpha = Expressions(1){alpha = 3.192;}
         beta  = Expressions(1){beta  = 26.6;}
```

**ChloricemSCM (NE-SCM) :**
```
Curves = AdSCM  x = Range{x1=0.5, x2=1.5, n=100}
         gamma_max = Expressions(1){gamma_max = 0.2 ;}
         k_acid    = Expressions(1){k_acid    = 1.e-14 ;}
         k_bas     = Expressions(1){k_bas     = 1.e12 ;}
         k_cl      = Expressions(1){k_cl      = 100. ;}
```

Tous les autres paramètres (maillage, matériau, conditions initiales et aux limites, pas de temps) sont identiques à ceux de l'exemple Chloricem.

---

## 6. Calibration et valeurs indicatives

Les valeurs proposées dans le cas test sont calées pour **reproduire approximativement** le comportement de l'isotherme de Langmuir de référence ($\alpha = 3.192$, $\beta = 26.6$) à un pH de référence de 12.5 ($c_\text{H} = 3.16 \times 10^{-13}$ mol/L).

| Paramètre | Valeur exemple | Commentaire |
|-----------|---------------|-------------|
| $\Gamma_\text{max}$ | 0.2 mol/mol$_\text{CSH}$ | Densité totale de sites, à calibrer sur BET + titrimetrie |
| $K_\text{acid}$ | $10^{-14}$ mol/L | p$K_{a1}$ = 14 — site très peu acide, favorable à la protonation en milieu alcalin |
| $K_\text{bas}$ | $10^{12}$ L/mol | p$K_b^{-1}$ = 12 — site très basique, $\equiv$SiOH$_2^+$ stablement formé à pH 12–13 |
| $K_\text{Cl}$ | 100 L/mol | $\log K_\text{Cl}$ = 2 — complexation modérée |

> **Attention :** ces valeurs n'ont pas de justification thermodynamique directe ; elles sont issues d'un recalage sur l'isotherme de Langmuir. Pour une approche rigoureuse, utiliser le modèle d'Elakneswaran (voir `ChloricemEla`), dont les constantes sont directement mesurables par électrophorèse et titrimetrie de surface.

---

## 7. Références bibliographiques

- Hao, I. (1986). *Ion binding on mineral surfaces*. Adv. Colloid Interface Sci.
- Stumm, W. (1992). *Chemistry of the Solid-Water Interface*. Wiley.
- Viallis-Terrisse, H., Nonat, A., Petit, J.-C. (2001). Zeta-potential study of calcium silicate hydrates interacting with alkaline cations. *J. Colloid Interface Sci.*, **244**, 58–65. doi:[10.1006/jcis.2001.7897](https://doi.org/10.1006/jcis.2001.7897)
- Elakneswaran, Y., Nawa, T., Kurumisawa, K. (2009). Electrokinetic potential of hydrated cement in relation to adsorption of chlorides. *Cement and Concrete Research*, **39**, 340–344. doi:[10.1016/j.cemconres.2009.01.006](https://doi.org/10.1016/j.cemconres.2009.01.006)
