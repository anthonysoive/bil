# M2 — Écoulement diphasique liquide-gaz 1D : portage Bil → Julia/VoronoiFVM

**Fichiers sources :**
- Modèle C : `src/Models/ModelFiles/M2.c`
- Cas test Bil : `base/M2/m2`
- Courbes de rétention : `base/M2/sl.c` / `base/M2/sol`
- Port Julia : `docs/M2_VoronoiFVM.jl`

---

## 1. Physique et description du cas test

### 1.1 Problème physique

On simule le **drainage gravitaire d'une colonne de sol verticale** initialement saturée en eau. Une colonne de 1 m est soumise à la gravité : l'eau s'écoule vers le bas tandis que l'air atmosphérique pénètre par le haut. Le modèle couple la conservation de la masse du liquide et du gaz dans un milieu poreux déformable.

**Géométrie :** domaine 1D $z \in [0,\, L]$ avec $L = 1\,\text{m}$, axe vertical, $z = 0$ en bas.

**Conditions aux limites :**

| Frontière | Condition |
|---|---|
| $z = 0$ (bas) | Dirichlet : $p_l = p_\text{atm} = 10^5\,\text{Pa}$ |
| $z = L$ (haut) | Dirichlet : $p_g = p_\text{atm} = 10^5\,\text{Pa}$ |

**Condition initiale :** colonne saturée, $p_l = p_g = p_\text{atm}$ partout, soit $p_c = 0$ et $S_l = 1$.

### 1.2 Paramètres matériau

| Symbole | Valeur | Unité | Description |
|---|---|---|---|
| $g$ | $-9.81$ | m/s² | Accélération gravitationnelle (vers le bas) |
| $\phi$ | $0.3$ | — | Porosité |
| $\rho_l$ | $1000$ | kg/m³ | Masse volumique du liquide |
| $k_\text{int}$ | $4.4 \times 10^{-13}$ | m² | Perméabilité intrinsèque |
| $\mu_l$ | $10^{-3}$ | Pa·s | Viscosité dynamique liquide |
| $\mu_g$ | $1.8 \times 10^{-5}$ | Pa·s | Viscosité dynamique gaz |
| $M_g/RT$ | $28.8 \times 10^{-3} / 2436$ | kg/J | Constante gaz idéal (air, $T = 293\,\text{K}$) |

---

## 2. Équations mathématiques

### 2.1 Variables primaires

Deux inconnues en chaque point :

$$u = (p_l,\; p_g)$$

La pression capillaire est définie par :

$$p_c = p_g - p_l \;\geq 0$$

### 2.2 Équations de conservation

**Conservation de la masse du liquide :**

$$\frac{\partial}{\partial t}\!\left(\rho_l\,\phi\,S_l\right) + \nabla \cdot \mathbf{W}_l = 0$$

**Conservation de la masse du gaz :**

$$\frac{\partial}{\partial t}\!\left(\rho_g\,\phi\,S_g\right) + \nabla \cdot \mathbf{W}_g = 0$$

avec $S_g = 1 - S_l$.

### 2.3 Loi de Darcy généralisée (flux)

Les flux de masse sont donnés par :

$$\mathbf{W}_l = -K_L\,\nabla H_l, \qquad \mathbf{W}_g = -K_G\,\nabla H_g$$

où les **têtes hydrauliques** incluent la contribution gravitationnelle :

$$H_l = p_l - \rho_l\, g\, z, \qquad H_g = p_g - \rho_g\, g\, z$$

Les **conductivités hydrauliques** sont :

$$K_L = \frac{\rho_l\, k_\text{int}\, k_{rl}(p_c)}{\mu_l}, \qquad K_G = \frac{\rho_g\, k_\text{int}\, k_{rg}(p_c)}{\mu_g}$$

### 2.4 Loi des gaz parfaits

La masse volumique du gaz est liée à sa pression par :

$$\rho_g = p_g \cdot \frac{M_g}{RT}$$

ce qui rend le terme de stockage $\rho_g\,\phi\,S_g$ non linéaire en $p_g$.

### 2.5 Courbes de rétention (matériau "sol")

Les relations $S_l(p_c)$, $k_{rl}(p_c)$ et $k_{rg}(p_c)$ sont définies analytiquement dans `base/M2/sl.c` :

**Saturation :**

$$S_l^\text{brut}(p_c) = 1 - 0.096863\left(\frac{p_c}{9806.65}\right)^{2.4279495}$$

**Perméabilité relative liquide :**

$$k_{rl}(p_c) = 1 - 2.207\,(1 - S_l)^{1.0121}$$

**Saturation effective et perméabilité relative gaz :**

$$\tilde{S}_l = \frac{S_l - 0.2}{0.8}, \qquad k_{rg}(p_c) = \left(1 - \tilde{S}_l\right)^2 \left(1 - \tilde{S}_l^{5/3}\right)$$

### 2.6 Régularisation de la saturation au voisinage de $S_l = 1$

Au voisinage de la saturation totale ($p_c \to 0$), $\mathrm{d}S_l/\mathrm{d}p_c \to 0$ et la matrice de capacité devient singulière. Le modèle Bil applique un **raccordement exponentiel** pour $p_c < p_{c3}$ (avec $p_{c3} = 300\,\text{Pa}$) :

$$S_l(p_c) = \begin{cases}
S_l^\text{brut}(p_c) & \text{si } p_c \geq p_{c3} \\[4pt]
1 - \bigl(1 - S_l^\text{brut}(p_{c3})\bigr)\,\exp\!\left(\dfrac{p_c - p_{c3}}{p_{c3}}\right) & \text{si } p_c < p_{c3}
\end{cases}$$

La continuité en $p_c = p_{c3}$ est assurée par construction. La dérivée $\mathrm{d}S_l/\mathrm{d}p_c$ reste non nulle pour $p_c < p_{c3}$, ce qui rend la matrice de capacité définie positive.

---

## 3. Discrétisation — Volumes finis (FVM)

Bil utilise une discrétisation **volumes finis centrés sur les cellules** (méthode deux-points, `src/Models/Methods/FVM.h`). Pour chaque arête $(i,j)$ séparant les cellules $i$ et $j$ de distance $d_{ij}$, le flux discret est :

$$W_{l,ij} = -K_L^\text{moy} \cdot \frac{H_{l,j} - H_{l,i}}{d_{ij}}$$

Le **coefficient de transfert** $K_L^\text{moy}$ est calculé comme la moyenne arithmétique des valeurs aux deux nœuds (fonctions `ComputeTransferCoefficients` et `ComputeFluxes` dans `M2.c`).

Le terme de stockage est intégré par quadrature de point médian sur chaque cellule de contrôle.

---

## 4. Portage Bil → Julia/VoronoiFVM

### 4.1 Correspondance des composants

| Composant Bil (C) | Équivalent Julia (VoronoiFVM) |
|---|---|
| `pm()` + `dm2()` — paramètres matériau | `const` globales |
| `Curve_ComputeValue()` — courbes tabulées | Fonctions analytiques `Sl_raw`, `krl`, `krg` |
| `saturation()` — régularisation | Fonction `Sl(pc)` |
| `ComputeTransferCoefficients()` | Calcul en ligne dans `flux!` |
| `ComputeFluxes()` / `Fluxes()` | Callback `flux!(f, u, edge, _)` |
| `ComputeImplicitTerms()` | Callback `storage!(f, u, _, __)` |
| `ComputeLoads()` / conditions limites | Callback `bcondition!(f, u, node, _)` |
| `ComputeMatrix()` + `TangentCoefficients()` | Jacobien automatique via `ForwardDiff.jl` |
| `ComputeResidu()` | Résidu automatique dans `VoronoiFVM.solve` |
| `TEMP` + `ALGO Dtini/Dtmax/Tol` | `VoronoiFVM.SolverControl` |

### 4.2 Encodage des courbes : tabulé → analytique

Dans Bil, les courbes de rétention sont lues depuis un fichier texte (`courbes = sol`) et évaluées par interpolation linéaire via `Curve_ComputeValue`. Le fichier `base/M2/sol` est généré par le programme `base/M2/sl.c` qui contient les formules analytiques explicites.

**Choix retenu :** utiliser directement les formules analytiques de `sl.c` dans Julia, ce qui évite de lire un fichier externe et rend le script autonome.

### 4.3 Convention de flux VoronoiFVM

VoronoiFVM appelle `flux!(f, u, edge, _)` pour chaque arête et **divise ensuite** `f[s]` par la longueur de l'arête $(x_2 - x_1)$. Pour reproduire un flux physique $W = -K\,\mathrm{d}H/\mathrm{d}x$, il faut donc fournir :

$$f[s] = K\,(H_1 - H_2)$$

En développant la tête hydraulique $H = p - \rho\,g\,z$ :

$$f[s] = K\,(p_1 - p_2) \underbrace{-\, K\,\rho\,g\,(z_1 - z_2)}_{\text{terme gravitaire}}$$

Le terme gravitaire est obtenu naturellement sans traitement spécial, contrairement à certains solveurs où il doit être ajouté séparément comme terme source.

> Ce point est identique à l'implémentation dans `M1_VoronoiFVM.jl`.

### 4.4 Densité du gaz variable : traitement à l'arête

Dans `M2.c`, la densité du gaz $\rho_g = p_g \cdot M_g/RT$ est calculée à chaque nœud et moyennée à l'arête pour les coefficients de transfert (`ComputeTransferCoefficients`). En Julia, on applique le même traitement :

```julia
rho_g_m = (pg1 * MgsRT + pg2 * MgsRT) / 2
```

La tête hydraulique du gaz utilise cette densité moyenne, ce qui introduit une légère approximation (termes en $\rho_g \cdot g \cdot z \ll p_g$) tout à fait justifiée car $M_g g z / (RT) \approx 10^{-4}$ pour $z \leq 1\,\text{m}$.

### 4.5 Jacobien : différentiation automatique

Dans Bil, le Jacobien est calculé par **différences finies directionnelles** (`ComputeComponentDerivatives` avec incrément `dxi = 1e2`). VoronoiFVM utilise **`ForwardDiff.jl`** (différentiation automatique en mode direct) : les callbacks `flux!` et `storage!` sont évalués sur des `Dual` numbers, ce qui donne un Jacobien exact au niveau machine, sans réglage d'incrément.

Cela implique que toutes les expressions dans les callbacks doivent être différentiables. Les fonctions `clamp`, `max`, `min`, `exp`, `^` sont toutes compatibles `ForwardDiff`.

### 4.6 Paramètres du solveur temporel

| Paramètre Bil (`base/M2/m2`) | Paramètre `SolverControl` Julia | Valeur |
|---|---|---|
| `Dtini = 0.1` | `Δt` | $0.1\,\text{s}$ |
| `Dtmax = 3600` | `Δt_max` | $3600\,\text{s}$ |
| `Tol = 1e-10` | `reltol` | $10^{-10}$ |
| `OBJE p_l = 1e3, p_g = 1e3` | `Δu_opt` | $10^3\,\text{Pa}$ |
| `TEMP 0 600 3600` | `times = tsave` | $[0, 600, 3600]\,\text{s}$ |

### 4.7 Conditions aux limites

Dans Bil, les conditions sont déclarées via les blocs `COND` et `FONC` :
```
Reg = 1 Inc = p_l Champ = 2   # p_l = 1e5 Pa en x=0
Reg = 3 Inc = p_g Champ = 2   # p_g = 1e5 Pa en x=L
```

`simplexgrid` de `ExtendableGrids.jl` assigne automatiquement :
- `BFaceRegion = 1` à $x = 0$
- `BFaceRegion = 2` à $x = L$

La correspondance Bil → Julia est donc directe via `boundary_dirichlet!`.

---

## 5. Résultats et validation

### 5.1 Sorties à t = 0, 600, 3600 s (nœud haut, $z = L$)

| $t$ (s) | $p_l$ (Pa) | $p_g$ (Pa) | $p_c$ (Pa) | $S_l$ (-) |
|---|---|---|---|---|
| 0 | $1.0000 \times 10^5$ | $1.0000 \times 10^5$ | 0 | 0.999993 |
| 600 | $9.1623 \times 10^4$ | $1.0000 \times 10^5$ | 8377 | 0.9339 |
| 3600 | $9.0876 \times 10^4$ | $1.0000 \times 10^5$ | 9124 | 0.9187 |

### 5.2 Comparaison avec la référence Bil (t = 600 s, $z = L$)

| Grandeur | VoronoiFVM | Bil `m2.t1` | Écart relatif |
|---|---|---|---|
| $p_l$ (Pa) | $9.162 \times 10^4$ | $9.155 \times 10^4$ | < 0.1 % |
| $p_c$ (Pa) | 8 377 | 8 450 | < 1 % |
| $S_l$ (-) | 0.9339 | ~0.933 | < 0.1 % |

### 5.3 Cohérence physique

- **Front de drainage** : à $t = 3600\,\text{s}$, la désaturation est localisée entre $z = 0.9\,\text{m}$ et $z = 1.0\,\text{m}$ ; la colonne inférieure reste quasi-saturée ($S_l \approx 0.99999$).
- **Profil de $p_l$** dans la zone saturée : gradient quasi-hydrostatique $\approx -\rho_l\,g \approx -9810\,\text{Pa/m}$ (9.1 Pa/cm observé, valeur attendue 9.81 Pa/cm — légère différence due au flux non nul).
- **Conservation** : $p_g = p_\text{atm}$ strictement respecté au nœud haut (condition de Dirichlet).
- **Nombre de pas de temps** : 87 (adaptatif, partant de $\Delta t = 0.1\,\text{s}$).

---

## 6. Limitations et pistes d'amélioration

| Limitation | Description | Solution envisageable |
|---|---|---|
| Courbes non génériques | Les formules analytiques de `sl.c` sont codées en dur | Lire le fichier `sol` tabulé avec `Interpolations.jl` |
| Gaz à densité variable | Approximation $\rho_g^\text{moy}$ à l'arête pour la tête hydraulique | Formuler le flux en $p_g$ uniquement et traiter $\rho_g$ implicitement |
| 1D uniquement | Le script est limité à une colonne 1D | Généraliser à 2D/3D avec `simplexgrid` ou maillage `gmsh` |
| Pas de sortie VTK | Seul un affichage console est produit | Ajouter `VTKView.jl` ou `WriteVTK.jl` |
