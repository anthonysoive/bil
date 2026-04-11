# Modèle ChloricemEla — Pénétration des chlorures avec le modèle de complexation de surface d'Elakneswaran

> **Fichier source :**
> `src/Models/ModelFiles/ChloricemEla.cpp`
>
> **Exemple :**
> `base/ChloricemEla/ChloricemEla` · `base/ChloricemEla/ChloricemEla.msh`
>
> **Modèle parent :** Chloricem (isotherme de Langmuir)
>
> **Référence principale :** @elakneswaran2009
>
> **Voir aussi :** ChloricemSCM (NE-SCM générique sur silanols protonés)

---

## Table des matières

1. [Contexte et objectif](#1-contexte-et-objectif)
2. [Différence avec Chloricem et ChloricemSCM](#2-différence-avec-chloricem-et-chloricemscm)
3. [Modèle de complexation de surface d'Elakneswaran](#3-modèle-de-complexation-de-surface-delakneswaran)
   - 3.1 [Mécanisme : Ca²⁺ comme cation pontant](#31-mécanisme--ca²⁺-comme-cation-pontant)
   - 3.2 [Réactions de surface](#32-réactions-de-surface)
   - 3.3 [Dérivation de la formule d'adsorption](#33-dérivation-de-la-formule-dadsorption)
   - 3.4 [Dépendance au pH et à la concentration en Ca²⁺](#34-dépendance-au-ph-et-à-la-concentration-en-ca²⁺)
   - 3.5 [Lien avec l'isotherme de Langmuir](#35-lien-avec-lisotherme-de-langmuir)
4. [Paramètres du modèle](#4-paramètres-du-modèle)
5. [Fichier d'entrée](#5-fichier-dentrée)
6. [Valeurs des paramètres](#6-valeurs-des-paramètres)
7. [Comparaison des trois modèles d'adsorption](#7-comparaison-des-trois-modèles-dadsorption)
8. [Références bibliographiques](#8-références-bibliographiques)

---

## 1. Contexte et objectif

Le modèle **ChloricemEla** est une variante du modèle Chloricem dans laquelle l'isotherme de Langmuir pour la fixation des chlorures par les C-S-H est remplacée par le **modèle de complexation de surface non-électrostatique (NE-SCM) proposé par Elakneswaran et al. (2009)**.

La spécificité de ce modèle par rapport au modèle NE-SCM générique (ChloricemSCM) est le rôle central joué par le **cation calcium Ca²⁺**, dont la concentration en solution est calculée dynamiquement par le solveur `HardenedCementChemistry` à chaque pas de temps. La fixation de Cl⁻ dépend donc à la fois du pH et de la composition ionique de la solution interstitielle.

Toutes les autres équations — transport multi-ionique, conservation des éléments, chimie des phases solides — sont identiques à celles de Chloricem.

---

## 2. Différence avec Chloricem et ChloricemSCM

| Modèle | Formule d'adsorption | Variables | Paramètres |
|--------|----------------------|-----------|------------|
| **Chloricem** | Langmuir empirique | $c_\text{Cl}$ | $\alpha(x_\text{csh}),\;\beta(x_\text{csh})$ |
| **ChloricemSCM** | NE-SCM sur ≡SiOH₂⁺ | $c_\text{Cl},\;c_\text{H}$ | $\Gamma_\text{max}, K_\text{acid}, K_\text{bas}, K_\text{Cl}$ |
| **ChloricemEla** | NE-SCM d'Elakneswaran | $c_\text{Cl},\;c_\text{H},\;\mathbf{c_\text{Ca^{2+}}}$ | $\Gamma_\text{max}, K_{a1}, K_\text{Ca}, K_\text{Cl}$ |

La différence essentielle entre ChloricemSCM et ChloricemEla tient au **mécanisme de fixation** :

- **ChloricemSCM** : Cl⁻ se lie directement sur un site protoné ≡SiOH₂⁺ (favorable en conditions acides/neutres).
- **ChloricemEla** : Cl⁻ se lie via le Ca²⁺ adsorbé sur un site déprotoné ≡SiO⁻ (favorable en milieu **alcalin** — cohérent avec la solution interstitielle du ciment, pH 12–13).

---

## 3. Modèle de complexation de surface d'Elakneswaran

### 3.1 Mécanisme : Ca²⁺ comme cation pontant

À pH 12–13, la surface des C-S-H est **fortement chargée négativement** (les silanols sont déprotonés : ≡SiO⁻ domine). Cette charge négative attire les cations de la solution (Ca²⁺, Na⁺, K⁺). Elakneswaran et al. ont montré par mesures de potentiel ζ que le Ca²⁺ s'adsorbe préférentiellement et peut même **inverser le signe du potentiel de surface**, créant des sites ≡SiOCa⁺ chargés positivement capables de fixer l'anion Cl⁻ par complexe de sphère externe.

### 3.2 Réactions de surface

Trois réactions d'équilibre sur les sites silanols des C-S-H :

| N° | Réaction | Constante |
|----|----------|-----------|
| (1) | $\equiv\text{SiOH} \rightleftharpoons \equiv\text{SiO}^- + \text{H}^+$ | $K_{a1}$ (mol/L) |
| (2) | $\equiv\text{SiO}^- + \text{Ca}^{2+} \rightleftharpoons \equiv\text{SiOCa}^+$ | $K_\text{Ca}$ (L/mol) |
| (3) | $\equiv\text{SiOCa}^+ + \text{Cl}^- \rightleftharpoons \equiv\text{SiOCaCl}$ | $K_\text{Cl}$ (L/mol) |

La réaction (3) est un **complexe de sphère externe** (outer-sphere) : l'ion Cl⁻ reste en solution et interagit électrostatiquement avec le site ≡SiOCa⁺, sans liaison covalente.

### 3.3 Dérivation de la formule d'adsorption

Les concentrations de surface des quatre espèces s'expriment par rapport à $[\equiv\text{SiOH}]$ via les lois d'action de masse :

$$[\equiv\text{SiO}^-]     = \frac{K_{a1}}{c_\text{H}}\,[\equiv\text{SiOH}]$$

$$[\equiv\text{SiOCa}^+]   = K_\text{Ca}\,c_\text{Ca}\,[\equiv\text{SiO}^-]
                            = \frac{K_{a1}\,K_\text{Ca}\,c_\text{Ca}}{c_\text{H}}\,[\equiv\text{SiOH}]$$

$$[\equiv\text{SiOCaCl}]   = K_\text{Cl}\,c_\text{Cl}\,[\equiv\text{SiOCa}^+]
                            = \frac{K_{a1}\,K_\text{Ca}\,K_\text{Cl}\,c_\text{Ca}\,c_\text{Cl}}{c_\text{H}}\,[\equiv\text{SiOH}]$$

Le bilan de masse sur les sites totaux $\Gamma_\text{max}$ donne :

$$\Gamma_\text{max} = [\equiv\text{SiOH}]\left(1 + \frac{K_{a1}}{c_\text{H}}\Bigl[1 + K_\text{Ca}\,c_\text{Ca}\,(1 + K_\text{Cl}\,c_\text{Cl})\Bigr]\right)$$

La quantité de chlorures adsorbés par mole de C-S-H est alors :

$$\boxed{
\Gamma_\text{Cl}(c_\text{Cl}, c_\text{H}, c_\text{Ca}, x_\text{csh}) =
  \frac{\Gamma_\text{max}(x_\text{csh})\;\dfrac{K_{a1}}{c_\text{H}}\;K_\text{Ca}\;K_\text{Cl}\;c_\text{Ca}\;c_\text{Cl}}
       {1 + \dfrac{K_{a1}}{c_\text{H}}\,\Bigl[1 + K_\text{Ca}\,c_\text{Ca}\,(1 + K_\text{Cl}\,c_\text{Cl})\Bigr]}
}$$

où $c_\text{H} = [\text{H}^+]$ et $c_\text{Ca} = [\text{Ca}^{2+}]$ sont calculés à chaque pas par `HardenedCementChemistry`.

### 3.4 Dépendance au pH et à la concentration en Ca²⁺

**Influence du pH** (à $c_\text{Ca}$ fixé) :

- **pH élevé** ($c_\text{H} \ll K_{a1}$) : $K_{a1}/c_\text{H} \gg 1$ — les sites ≡SiO⁻ dominent. Si en plus $K_\text{Ca}\,c_\text{Ca} \gg 1$ :
  $$\Gamma_\text{Cl} \approx \frac{\Gamma_\text{max}\,K_\text{Cl}\,c_\text{Cl}}{1 + K_\text{Cl}\,c_\text{Cl}}$$
  Isotherme de Langmuir avec capacité $\Gamma_\text{max}$ — **adsorption maximale en milieu très alcalin et riche en Ca²⁺**.

- **pH faible** ($c_\text{H} \gg K_{a1}$) : $K_{a1}/c_\text{H} \ll 1$ — les sites ≡SiOH dominent, peu de ≡SiO⁻ disponibles pour adsorber Ca²⁺, donc peu de Cl⁻ adsorbé.
  $$\Gamma_\text{Cl} \approx \frac{\Gamma_\text{max}\,K_{a1}\,K_\text{Ca}\,K_\text{Cl}\,c_\text{Ca}\,c_\text{Cl}}{c_\text{H}}$$

**Influence de $c_\text{Ca}$** (à pH fixé) :

- **$c_\text{Ca} \to 0$** : aucun site ≡SiOCa⁺ disponible → $\Gamma_\text{Cl} \to 0$. C'est la différence fondamentale avec ChloricemSCM : **l'adsorption de Cl⁻ est nulle en l'absence de Ca²⁺**.
- **$c_\text{Ca}$ élevé** ($K_\text{Ca}\,c_\text{Ca} \gg 1$) : la formule se réduit à une isotherme de Langmuir en $c_\text{Cl}$ à capacité $\Gamma_\text{max}$.

Ce comportement est cohérent avec les observations expérimentales qui montrent que la substitution du Ca²⁺ par Na⁺ ou K⁺ dans la solution interstitielle réduit significativement la fixation des chlorures.

### 3.5 Lien avec l'isotherme de Langmuir

À pH et $c_\text{Ca}$ fixes, le modèle se réduit à une isotherme de Langmuir effective :

$$\Gamma_\text{Cl} = \frac{\alpha_\text{eff}\,c_\text{Cl}}{1 + \beta_\text{eff}\,c_\text{Cl}}$$

avec :

$$\alpha_\text{eff} = \frac{\Gamma_\text{max}\,(K_{a1}/c_\text{H})\,K_\text{Ca}\,c_\text{Ca}}
                         {1 + (K_{a1}/c_\text{H})\,(1 + K_\text{Ca}\,c_\text{Ca})}$$

$$\beta_\text{eff} = \frac{(K_{a1}/c_\text{H})\,K_\text{Ca}\,K_\text{Cl}\,c_\text{Ca}}
                        {1 + (K_{a1}/c_\text{H})\,(1 + K_\text{Ca}\,c_\text{Ca})}$$

Le rapport $\alpha_\text{eff}/\beta_\text{eff} = \Gamma_\text{max}$ est bien la capacité d'adsorption maximale, indépendante des conditions chimiques.

---

## 4. Paramètres du modèle

Les quatre paramètres peuvent dépendre du rapport Ca/Si $x_\text{csh}$ et sont fournis sous forme de courbes :

| Paramètre | Courbe Bil | Description | Unité |
|-----------|-----------|-------------|-------|
| $\Gamma_\text{max}$ | `gamma_max` | Densité totale de sites silanols | mol/mol$_\text{CSH}$ |
| $K_{a1}$ | `k_a1` | Constante de déprotonation de ≡SiOH | mol/L |
| $K_\text{Ca}$ | `k_ca` | Constante de complexation de Ca²⁺ sur ≡SiO⁻ | L/mol |
| $K_\text{Cl}$ | `k_cl` | Constante de complexation de Cl⁻ sur ≡SiOCa⁺ | L/mol |

La concentration en Ca²⁺ libre ($c_\text{Ca}$) n'est **pas un paramètre** : elle est calculée dynamiquement par `HardenedCementChemistry` en fonction de l'état des phases solides (dissolution de la portlandite, décalcification des C-S-H, précipitation de la calcite).

---

## 5. Fichier d'entrée

La différence avec le fichier `Chloricem` tient uniquement à la ligne `Curves` d'adsorption des chlorures et au nom du modèle :

**Chloricem (Langmuir) :**
```
Model = Chloricem
...
Curves = Adscl  x = Range{x1=0.5, x2=1.5, n=100}
         alpha = Expressions(1){alpha = 3.192;}
         beta  = Expressions(1){beta  = 26.6;}
```

**ChloricemEla (Elakneswaran) :**
```
Model = ChloricemEla
...
Curves = AdEla  x = Range{x1=0.5, x2=1.5, n=100}
         gamma_max = Expressions(1){gamma_max = 0.056 ;}
         k_a1      = Expressions(1){k_a1      = 5.01e-13 ;}
         k_ca      = Expressions(1){k_ca      = 100. ;}
         k_cl      = Expressions(1){k_cl      = 10.2 ;}
```

Tous les autres paramètres (maillage, matériau, conditions initiales et aux limites, pas de temps) sont identiques à ceux de l'exemple Chloricem.

---

## 6. Valeurs des paramètres

Les valeurs proposées dans le cas test sont directement issues des mesures d'Elakneswaran et al. (2009) sur des pâtes de ciment Portland :

| Paramètre | Valeur | Origine |
|-----------|--------|---------|
| $\Gamma_\text{max}$ | 0.056 mol/mol$_\text{CSH}$ | 2 sites/nm² × $S_\text{BET}$ = 150 m²/g × $M_\text{CSH}$ ≈ 170 g/mol |
| $K_{a1}$ | $5.01 \times 10^{-13}$ mol/L | p$K_{a1}$ = 12.3, mesuré par titrimetrie de surface |
| $K_\text{Ca}$ | 100 L/mol | $\log K_\text{Ca}$ = 2.0, ajusté sur potentiel ζ |
| $K_\text{Cl}$ | 10.2 L/mol | $\log K_\text{Cl}$ = 1.01, ajusté sur potentiel ζ |

**Calcul de $\Gamma_\text{max}$** :

$$\Gamma_\text{max} = \frac{N_s \times S_\text{BET} \times M_\text{CSH}}{N_\text{Av}}
= \frac{2 \times 10^{18}\,\text{sites/m}^2 \times 150\,\text{m}^2/\text{g} \times 170\,\text{g/mol}}{6.022 \times 10^{23}\,\text{mol}^{-1}}
\approx 0.085\;\text{mol/mol}_\text{CSH}$$

La valeur retenue (0.056 mol/mol$_\text{CSH}$) correspond à $S_\text{BET}$ = 100 m²/g, valeur conservatrice pour un C-S-H à Ca/Si = 0.8–1.2.

> **Remarque sur la surface spécifique :** $S_\text{BET}$ varie selon la méthode de mesure (N₂ ou eau), le degré d'hydratation et le rapport Ca/Si. Pour tenir compte de la variation avec $x_\text{csh}$, on peut définir `gamma_max` comme une courbe dépendant de $x_\text{csh}$ plutôt que comme une constante.

---

## 7. Comparaison des trois modèles d'adsorption

| Caractéristique | Langmuir | ChloricemSCM | ChloricemEla |
|----------------|----------|--------------|--------------|
| Dépendance en $c_\text{Cl}$ | ✓ | ✓ | ✓ |
| Dépendance en pH | ✗ | ✓ | ✓ |
| Dépendance en $c_\text{Ca^{2+}}$ | ✗ | ✗ | ✓ |
| Paramètres | 2 | 4 | 4 |
| Mesurabilité des paramètres | faible | modérée | élevée |
| Adsorption en l'absence de Ca²⁺ | non nulle | non nulle | **nulle** |
| Cohérence à pH 12–13 | ad hoc | discutable | bonne |
| Référence | empirique | générique | @elakneswaran2009 |

---

## 8. Références bibliographiques

- Elakneswaran, Y., Nawa, T., Kurumisawa, K. (2009). Electrokinetic potential of hydrated cement in relation to adsorption of chlorides. *Cement and Concrete Research*, **39**(4), 340–344. doi:[10.1016/j.cemconres.2009.01.006](https://doi.org/10.1016/j.cemconres.2009.01.006)
- Viallis-Terrisse, H., Nonat, A., Petit, J.-C. (2001). Zeta-potential study of calcium silicate hydrates interacting with alkaline cations. *Journal of Colloid and Interface Science*, **244**, 58–65. doi:[10.1006/jcis.2001.7897](https://doi.org/10.1006/jcis.2001.7897)
- Stumm, W. (1992). *Chemistry of the Solid-Water Interface*. Wiley-Interscience, New York.
