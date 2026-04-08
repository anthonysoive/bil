# Architecture et fonctionnement de Bil

> **Sources :**
> `src/Models/Model.h` · `src/Models/Model.c` · `src/Models/Models.h`
> `src/Models/ListOfModels.h` · `src/Models/ModelFiles/Template.c`
> `src/Models/ModelFiles/TemplateFEM.c` · `src/Models/ModelFiles/TemplateFVM.c`
> `src/Modules/ModuleFiles/Monolithic.c`
> `src/Models/Methods/MaterialPointMethod.h` · `src/Models/Methods/CustomValues.h`

---

## Table des matières

1. [Vue d'ensemble](#1-vue-densemble)
2. [Principe fondamental : architecture plugin en C](#2-principe-fondamental--architecture-plugin-en-c)
3. [L'interface `Model_t` — la vtable en C](#3-linterface-model_t--la-vtable-en-c)
4. [Le registre de modèles — `ListOfModels.h`](#4-le-registre-de-modèles--listofmodelsh)
5. [Résolution au démarrage — `Model_Initialize`](#5-résolution-au-démarrage--model_initialize)
6. [Ce que doit implémenter un modèle physique](#6-ce-que-doit-implémenter-un-modèle-physique)
7. [La boucle de calcul — `Monolithic.c`](#7-la-boucle-de-calcul--monolithicc)
8. [Flux de données complet](#8-flux-de-données-complet)
9. [Ajouter un nouveau modèle](#9-ajouter-un-nouveau-modèle)
10. [Interface moderne `MaterialPointMethod.h` (C++)](#10-interface-moderne-materialpointmethodh-c)
11. [Bases de données physico-chimiques](#11-bases-de-données-physico-chimiques)
12. [Méthode FEM² — homogénéisation numérique](#12-méthode-fem--homogénéisation-numérique)
13. [Gestion de l'historique des solutions](#13-gestion-de-lhistorique-des-solutions)
14. [Comparaison avec d'autres patterns](#14-comparaison-avec-dautres-patterns)

---

## 1. Vue d'ensemble

Bil est organisé en **deux couches strictement séparées** :

| Couche | Fichiers | Rôle |
|--------|----------|------|
| **Framework générique** | `Monolithic.c`, `SNIA.c`, `Solver.c`, `FEM.c`, `Model.c`… | Boucle temporelle, assemblage de la matrice globale, résolution du système linéaire, I/O — **ne sait rien de la physique** |
| **Plugin physique** | `M7.c`, `BBM.c`, `Richards.c`, `Fick.c`… | Implémente les équations : termes implicites, matrice tangente, résidu, sorties — **ne sait rien du solveur** |

Le lien entre les deux est un **contrat d'interface** matérialisé par le struct `Model_t`, qui est une table de pointeurs de fonctions (équivalent C d'une *vtable* C++).

---

## 2. Principe fondamental : architecture plugin en C

La question naturelle est : « puisque chaque simulation correspond à un fichier `.c` différent, où est la généricité ? »

La réponse est que **la généricité est dans le framework**, pas dans les modèles. Le framework appelle toujours les mêmes signatures de fonctions, peu importe la physique :

```
framework → model->computematrix(el, t, dt, K)
framework → model->computeimplicitterms(el, t, dt)
framework → model->computeresidu(el, t, dt, r)
```

C'est exactement le même mécanisme qu'une interface en Java ou une classe abstraite en C++, mais réalisé en C pur avec des pointeurs de fonctions.

```
                 ┌─────────────────────────────────────┐
                 │         Framework générique          │
                 │  Monolithic.c  ·  SNIA.c  ·  FEM.c  │
                 │                                     │
                 │   Mesh_ComputeMatrix(mesh, ...)      │
                 │   Mesh_ComputeResidu(mesh, ...)      │
                 │   Solver_Solve(solver, ...)          │
                 └──────────────┬──────────────────────┘
                                │ appels via pointeurs de fonctions
                                │ (Model_t comme interface)
          ┌─────────────────────┼──────────────────────────┐
          │                     │                          │
   ┌──────▼──────┐      ┌───────▼──────┐      ┌───────────▼──────┐
   │    M7.c     │      │   BBM.c      │      │   Richards.c     │
   │ Poroélas-   │      │ Barcelone    │      │ Écoulement non   │
   │ ticité non  │      │ Basic Model  │      │ saturé (Darcy)   │
   │ saturée     │      │ (plasticité) │      │                  │
   └─────────────┘      └──────────────┘      └──────────────────┘
```

---

## 3. L'interface `Model_t` — la vtable en C

Défini dans `src/Models/Model.h`, le struct `Model_t` contient **exclusivement des pointeurs de fonctions** :

```c
/* src/Models/Model.h — signatures des méthodes de l'interface */

typedef int  (Model_SetModelProperties_t)      (Model_t*) ;
typedef int  (Model_ReadMaterialProperties_t)  (Material_t*, DataFile_t*) ;
typedef int  (Model_PrintModelProperties_t)    (Model_t*, FILE*) ;
typedef int  (Model_DefineElementProperties_t) (Element_t*, IntFcts_t*, ShapeFcts_t*) ;
typedef int  (Model_ComputeInitialState_t)     (Element_t*, double) ;
typedef int  (Model_ComputeExplicitTerms_t)    (Element_t*, double) ;
typedef int  (Model_ComputeImplicitTerms_t)    (Element_t*, double, double) ;
typedef int  (Model_ComputeMatrix_t)           (Element_t*, double, double, double*) ;
typedef int  (Model_ComputeResidu_t)           (Element_t*, double, double, double*) ;
typedef int  (Model_ComputeLoads_t)            (Element_t*, double, double, Load_t*, double*) ;
typedef int  (Model_ComputeOutputs_t)          (Element_t*, double, double*, Result_t*) ;
typedef int  (Model_ComputePropertyIndex_t)    (const char*) ;
typedef void (Model_ComputeMaterialProperties_t)(Element_t*, double) ;

struct Model_t {
  Model_SetModelProperties_t*        setmodelprop ;        /* initialise les autres pointeurs */
  Model_ReadMaterialProperties_t*    readmatprop ;         /* lit les paramètres du fichier d'entrée */
  Model_PrintModelProperties_t*      printmodelprop ;      /* affiche l'aide du modèle */
  Model_DefineElementProperties_t*   defineelementprop ;   /* alloue les tableaux par élément */
  Model_ComputeInitialState_t*       computeinitialstate ; /* calcule l'état initial */
  Model_ComputeExplicitTerms_t*      computeexplicitterms ;/* termes explicites (perméabilité…) */
  Model_ComputeImplicitTerms_t*      computeimplicitterms ;/* termes implicites (stockage, flux…) */
  Model_ComputeMatrix_t*             computematrix ;       /* matrice de rigidité/couplage */
  Model_ComputeResidu_t*             computeresidu ;       /* résidu de Newton-Raphson */
  Model_ComputeLoads_t*              computeloads ;        /* vecteur de chargement */
  Model_ComputeOutputs_t*            computeoutputs ;      /* champs de post-traitement */
  Model_ComputePropertyIndex_t*      computepropertyindex ;/* indice d'une propriété matériau */
  Model_ComputeMaterialProperties_t* ComputeMaterialProperties ; /* propriétés au point de Gauss */
  /* + métadonnées : codename, shorttitle, authors, nbofequations, nameofequations… */
} ;
```

Le framework ne manipule que des `Model_t*`. Il ignore totalement si le modèle derrière est M7, BBM ou un modèle personnalisé.

---

## 4. Le registre de modèles — `ListOfModels.h`

`src/Models/ListOfModels.h` est le **fichier de déclaration de tous les modèles compilés**. Il utilise la technique de la **X-macro** pour générer à la fois la liste des noms et la liste des pointeurs de fonctions à partir d'une seule source de vérité :

```c
/* src/Models/ListOfModels.h */

#define ListOfModels_Nb   39

#define ListOfModels_Names \
  "BBM","BBMgas","BExM","CHMBWP","Cryspom","DWS1","Elasd","Elast","Fick", \
  "Frostaco","Frostaco3d","Gascoal","HydrateThermoPoroplasticity",         \
  "M1","M10","M2","M7","NSFSHao","Plastold","Shen",                        \
  "Sulfaco","Sulfaco3d","SulfacoESA3d","Sulfacocl","Sulfuricem",            \
  "TVIThermoPoroplast","Thermoporoplast","Yuan1","hydrapel","usoil",        \
  "BBMGas","Chloricem","Duracem","Frostsoil","MechaMic",                    \
  "Plast","Poroplast","Richards","Sulfaconew"

#define ListOfModels_Methods(m) \
  BBM##m, BBMgas##m, BExM##m, CHMBWP##m, Cryspom##m, DWS1##m,             \
  Elasd##m, Elast##m, Fick##m, Frostaco##m, Frostaco3d##m, Gascoal##m,    \
  HydrateThermoPoroplasticity##m,                                           \
  M1##m, M10##m, M2##m, M7##m, NSFSHao##m, Plastold##m, Shen##m,          \
  Sulfaco##m, Sulfaco3d##m, SulfacoESA3d##m, Sulfacocl##m,                 \
  Sulfuricem##m, TVIThermoPoroplast##m, Thermoporoplast##m,                 \
  Yuan1##m, hydrapel##m, usoil##m, BBMGas##m, Chloricem##m, Duracem##m,    \
  Frostsoil##m, MechaMic##m, Plast##m, Poroplast##m, Richards##m,          \
  Sulfaconew##m
```

**Comment fonctionne la X-macro :**

L'appel `ListOfModels_Methods(_SetModelProp)` est développé par le préprocesseur en :

```c
BBM_SetModelProp, BBMgas_SetModelProp, ..., M7_SetModelProp, ..., Sulfaconew_SetModelProp
```

Ce qui permet de construire, **statiquement à la compilation**, deux tableaux parallèles indexés de 0 à 38 :

```
index │ nom (modelnames[i])             │ fonction (xModel_SetModelProperties[i])
──────┼─────────────────────────────────┼─────────────────────────────────────────
  0   │ "BBM"                           │ BBM_SetModelProp
  1   │ "BBMgas"                        │ BBMgas_SetModelProp
  …   │ …                               │ …
 16   │ "M7"                            │ M7_SetModelProp
  …   │ …                               │ …
 38   │ "Sulfaconew"                    │ Sulfaconew_SetModelProp
```

---

## 5. Résolution au démarrage — `Model_Initialize`

Lorsque Bil lit `Model = M7` dans le fichier d'entrée, il appelle `Model_Initialize` (`src/Models/Model.c`) :

```c
/* src/Models/Model.c — simplifié */

Model_t* Model_Initialize(Model_t* model, const char* codename, ...)
{
  int n_models = Models_NbOfModels ;                               // 38
  const char* modelnames[] = {Models_ListOfNames} ;               // {"BBM","BBMgas",...,"M7",...}
  Model_SetModelProperties_t* xModel_SetModelProperties[] =
                              {Models_ListOfSetModelProp} ;        // {BBM_SetModelProp,...,M7_SetModelProp,...}

  int i = 0 ;
  while(i < n_models && strcmp(modelnames[i], codename)) i++ ;   // cherche "M7" → i=15

  Model_GetSetModelProperties(model) = xModel_SetModelProperties[i] ; // assigne M7_SetModelProp
  Model_SetModelProperties(model) ;   // appelle M7_SetModelProp(model)
                                      // → peuple TOUS les autres pointeurs de Model_t

  return(model) ;
}
```

Après cet appel, `model->computematrix` pointe vers `k7`, `model->computeresidu` vers `c7`, etc. — toutes les fonctions de `M7.c`.

---

## 6. Ce que doit implémenter un modèle physique

Chaque fichier modèle doit fournir une fonction `SetModelProp` et plusieurs fonctions de calcul. Voici le squelette minimal tiré de `src/Models/ModelFiles/Template.c` :

```c
/* MonModele.c */
#include "CommonModel.h"
#include "FEM.h"                 /* ou "FVM.h" selon la discrétisation */

#define TITLE   "Nom du modèle"
#define AUTHORS "Auteur"

#include "PredefinedModelMethods.h"  /* macro qui déclare SetModelProp automatiquement */

#define NEQ  2   /* nombre d'équations = nombre d'inconnues nodales */
#define NVI  9   /* termes implicites par point de Gauss */
#define NVE  2   /* termes explicites par point de Gauss */

/* --- 1. Registre des paramètres matériaux --- */
int pm(const char* s)
{
  if(strcmp(s,"young")   == 0) return 0 ;
  if(strcmp(s,"poisson") == 0) return 1 ;
  return -1 ;  /* paramètre inconnu */
}

/* --- 2. Propriétés du modèle (peuple Model_t) --- */
int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  Model_CopyNameOfEquation(model, 0, "meca") ;
  Model_CopyNameOfEquation(model, 1, "mass") ;
  Model_CopyNameOfUnknown(model,  0, "u")   ;
  Model_CopyNameOfUnknown(model,  1, "p_l") ;
  return 0 ;
}

/* --- 3. Lecture des paramètres du fichier d'entrée --- */
int ReadMatProp(Material_t* mat, DataFile_t* datafile)
{
  Material_ScanProperties(mat, datafile, pm) ;
  return 2 ;  /* nb de paramètres scalaires */
}

/* --- 4. Termes explicites : calculés avec les valeurs au pas précédent --- */
int ComputeExplicitTerms(Element_t* el, double t) { ... }

/* --- 5. Termes implicites : stockage, flux, contraintes… --- */
int ComputeImplicitTerms(Element_t* el, double t, double dt) { ... }

/* --- 6. Matrice de rigidité tangente --- */
int ComputeMatrix(Element_t* el, double t, double dt, double* k) { ... }

/* --- 7. Résidu pour Newton-Raphson --- */
int ComputeResidu(Element_t* el, double t, double dt, double* r) { ... }

/* --- 8. Post-traitement : champs de sortie --- */
int ComputeOutputs(Element_t* el, double t, double* s, Result_t* r) { ... }
```

**Correspondance avec les noms historiques de M7.c** (interface ancienne via `OldMethods.h`) :

| Méthode interface moderne | Nom dans M7.c (ancienne API) | Rôle |
|---------------------------|------------------------------|------|
| `SetModelProp`            | `dm7` + initialisation       | Déclare les équations et inconnues |
| `ComputeImplicitTerms`    | `mxnd`                       | Calcule $M_l$, $W_l$, $\sigma$, $\phi$ aux points de Gauss |
| `ComputeMatrix`           | `k7`                         | Assemble la matrice de rigidité $K$ |
| `ComputeResidu`           | `c7`                         | Calcule le résidu $R = F_\text{ext} - F_\text{int}$ |
| `ComputeInitialState`     | `rsnd`                       | Initialise les champs au temps $t=0$ |
| `ComputeOutputs`          | `so7` (via `Views`)          | Produit les champs de post-traitement |

---

## 7. La boucle de calcul — `Monolithic.c`

Le module `src/Modules/ModuleFiles/Monolithic.c` est le chef d'orchestre. Sa boucle interne (simplifiée) est :

```
Pour chaque pas de temps [t_n → t_{n+1}] :
│
├── Mesh_ComputeExplicitTerms(mesh, t_n)
│     └── pour chaque élément : model->computeexplicitterms(el, t_n)
│         [calcule perméabilité, saturation… à partir du pas précédent]
│
├── Itérations de Newton-Raphson :
│   │
│   ├── Mesh_ComputeImplicitTerms(mesh, t_{n+1}, dt)
│   │     └── pour chaque élément : model->computeimplicitterms(el, t, dt)
│   │         [calcule flux, stockage, contraintes avec les inconnues courantes]
│   │
│   ├── Mesh_ComputeMatrix(mesh, a, t_{n+1}, dt)
│   │     └── pour chaque élément : model->computematrix(el, t, dt, K_e)
│   │         [assemble la matrice tangente élémentaire K_e dans K_global]
│   │
│   ├── Mesh_ComputeResidu(mesh, r, loads, t_{n+1}, dt)
│   │     └── pour chaque élément : model->computeresidu(el, t, dt, r_e)
│   │         [assemble le résidu r_global = F_ext - F_int]
│   │
│   ├── Solver_Solve(solver, K_global, r_global, du)
│   │     [résout K · Δu = r par SuperLU / PETSc / Crout…]
│   │
│   └── Mise à jour des inconnues : u ← u + Δu
│       Vérification de convergence (||Δu|| / ||u|| < tol)
│
└── Si convergence : avance t_n ← t_{n+1}, écrit les sorties si t ∈ Dates
```

À **aucun moment** `Monolithic.c` ne référence `M7`, `BBM` ou tout autre nom de modèle. Il ne passe que par les pointeurs de `Model_t`.

---

## 8. Flux de données complet

```
Fichier d'entrée (M7-1)
  "Model = M7"
        │
        ▼
  Parser (Flex/Bison)
        │ lit la chaîne "M7"
        ▼
  Model_Initialize("M7")                    [Model.c]
        │ recherche linéaire dans modelnames[]
        │ assigne xModel_SetModelProperties[15] = M7_SetModelProp
        │
        ▼
  M7_SetModelProp(model)                    [M7.c]
        │ model->computematrix        = k7
        │ model->computeresidu        = c7
        │ model->computeimplicitterms = mxnd
        │ model->computeinitialstate  = rsnd
        │ model->readmatprop          = dm7 (via pm)
        │ model->computeoutputs       = so7
        │
        ▼
  Monolithic.c — boucle temporelle
        │ appelle model->computeexplicitterms(el, t)  → physique M7
        │ appelle model->computeimplicitterms(el,t,dt)→ physique M7
        │ appelle model->computematrix(el, t, dt, K)  → physique M7
        │ appelle model->computeresidu(el, t, dt, r)  → physique M7
        │
        ▼
  Solver.c — résolution K·Δu = r           [SuperLU / PETSc / …]
        │
        ▼
  Fichiers de résultats (.t0, .t1, …)
```

---

## 9. Ajouter un nouveau modèle

La procédure est minimaliste par conception :

**Étape 1** — Créer `src/Models/ModelFiles/MonModele.c` (ou `.cpp`) en s'appuyant sur `Template.c`, `TemplateFEM.c` ou `TemplateFVM.c`.

**Étape 2** — Enregistrer le modèle dans `src/Models/ListOfModels.h` :

```c
/* Avant */
#define ListOfModels_Nb   39
#define ListOfModels_Names  ..., "Sulfaconew"
#define ListOfModels_Methods(m)  ..., Sulfaconew##m

/* Après */
#define ListOfModels_Nb   40
#define ListOfModels_Names  ..., "Sulfaconew", "MonModele"
#define ListOfModels_Methods(m)  ..., Sulfaconew##m, MonModele##m
```

**Étape 3** — Recompiler :

```bash
cd build && make
```

Le fichier d'entrée peut maintenant utiliser `Model = MonModele`. **Aucun autre fichier du framework n'est modifié.**

---

## 10. Interface moderne `MaterialPointMethod.h` (C++)

Pour les modèles C++, Bil fournit une interface de haut niveau qui simplifie considérablement l'implémentation. Au lieu d'écrire manuellement les 11 fonctions, on dérive une classe `MPM_t` depuis `MaterialPointMethod_t`.

### 10.1 Définir les variables internes avec `CustomValues.h`

```cpp
#include "CustomValues.h"

template<typename T> struct ImplicitValues_t;  /* varient à chaque itération */
template<typename T> struct ExplicitValues_t;  /* varient à chaque pas de temps */
template<typename T> struct ConstantValues_t;  /* constantes au cours du temps */

template<typename T>
using V = CustomValues_t<T, ImplicitValues_t, ExplicitValues_t, ConstantValues_t, ...>;
```

Les trois types de variables internes correspondent aux tableaux par point de Gauss :

| Type | Accès dans l'ancienne API | Quand recalculé |
|------|--------------------------|-----------------|
| Implicites (`NVI`) | `vi[k]` | À chaque itération de Newton |
| Explicites (`NVE`) | `ve[k]` | Une fois par pas de temps (avec les valeurs au pas précédent) |
| Constants (`NV0`) | `v0[k]` | Jamais après l'initialisation |

### 10.2 Dériver la classe `MPM_t`

```cpp
#include "MaterialPointMethod.h"

struct MPM_t : public MaterialPointMethod_t<V> {

  /* Initialise val avec les inconnues nodales et leurs gradients */
  V<double>* SetInputs(Element_t* el, double const& t, int const& p,
                       double const* const* u, V<double>& val);

  /* Initialise toutes les variables internes (état initial) */
  V<double>* Initialize(Element_t* el, double const& t, V<double>& val);

  /* Intègre la loi de comportement de t-dt à t */
  V<double>* Integrate(Element_t* el, double const& t, double const& dt,
                       V<double> const& val_n, V<double>& val);

  /* Remplit la k-ième colonne de la matrice tangente */
  int SetTangentMatrix(Element_t* el, double const& t, double const& dt,
                       int const& p, V<double> const& val, V<double> const& dval,
                       int const& k, double* c);

  /* Remplit la matrice de transfert (FVM) */
  int SetTransferMatrix(Element_t* el, double const& dt, int const& p,
                        V<double> const& val, double* c);

  /* Calcule les flux entre les nœuds i et j (FVM) */
  V<double>* SetFluxes(Element_t* el, double const& t, int const& i, int const& j,
                       V<double> const& grdval, V<double>* val);

  /* Donne l'indice dans V de chaque inconnue primaire */
  void SetIndexOfPrimaryVariables(Element_t* el, int* ind);

  /* Donne l'incrément Δu utilisé pour les dérivées numériques */
  void SetIncrementOfPrimaryVariables(Element_t* el, double* dui);
};
```

### 10.3 Différentiation automatique avec `autodiff.h`

La matrice tangente peut être obtenue par **différentiation automatique** plutôt que par différences finies :

```cpp
#define USE_AUTODIFF
#include "autodiff.h"
```

Quand `USE_AUTODIFF` est défini, l'opérateur `Differentiate` calcule analytiquement la dérivée de `Integrate` par rapport aux inconnues primaires, sans implémentation manuelle du Jacobien.

Sinon, la dérivée numérique est utilisée automatiquement via `SetIncrementOfPrimaryVariables` :

```cpp
void SetIncrementOfPrimaryVariables(Element_t* el, double* dui) {
  ObVal_t* obval = Element_GetObjectiveValue(el);
  dui[0] = 1.e-2 * ObVal_GetValue(obval + 0);  /* incrément = 1% de la variation objective */
}
```

---

## 11. Bases de données physico-chimiques

Bil fournit une bibliothèque de données physico-chimiques accessible depuis n'importe quel modèle via des en-têtes d'inclusion. Ces données sont stockées dans `src/Models/DataBases/`.

```c
/* Viscosité de l'eau en fonction de la température (en K) */
#include "WaterViscosity.h"
double mu_w = WaterViscosity(293);   /* Pa·s à 293 K */

/* Coefficient de diffusion d'un ion dans l'eau */
#include "DiffusionCoefficientOfMoleculeInWater.h"
double d_ca = DiffusionCoefficientOfMoleculeInWater(Ca, 293);  /* m²/s */

/* Coefficient de diffusion dans l'air */
#include "DiffusionCoefficientOfMoleculeInAir.h"
double d_co2_air = DiffusionCoefficientOfMoleculeInAir(CO2, 293);
```

Parmi les bases disponibles :

| En-tête | Contenu |
|---------|---------|
| `WaterViscosity.h` | Viscosité de l'eau en fonction de T |
| `AirViscosity.h` | Viscosité de l'air |
| `DiffusionCoefficientOfMoleculeInWater.h` | Coefficients de diffusion ionique en solution aqueuse |
| `DiffusionCoefficientOfMoleculeInAir.h` | Coefficients de diffusion en phase gazeuse |
| `HardenedCementChemistry.h` | Chimie des liants hydrauliques hydratés |
| `CementSolutionChemistry.h` | Chimie de la solution interstitielle cimentaire |
| `ElectricChargeOfIonInWater.h` | Charges électriques des ions |
| `EquilibriumConstantOfHomogeneousReactionInWater.h` | Constantes d'équilibre des réactions en solution |
| `IonizationConstantOfWater.h` | Constante d'ionisation de l'eau |
| `HenrysLawConstantForSolubilityOfGasInWater.h` | Constante de Henry (solubilité des gaz) |
| `AtmosphericPressure.h` | Pression atmosphérique de référence |

---

## 12. Méthode FEM² — homogénéisation numérique

Bil dispose d'une capacité de **couplage multi-échelle FEM²** (ou FEM-in-FEM) : à chaque point de Gauss du maillage macro, un problème micro est résolu sur une cellule représentative.

Ce module est utilisé notamment dans le modèle `MechaMic.c` pour l'homogénéisation mécanique de matériaux hétérogènes. Les conditions aux limites périodiques sur la cellule micro sont gérées via le bloc `PERIODICITIES` du fichier d'entrée :

```
PERIODICITIES
2
MasterRegion = 10  SlaveRegion = 11  PeriodVector = 2 0 0
MasterRegion = 20  SlaveRegion = 21  PeriodVector = 0 2 0
```

---

## 13. Gestion de l'historique des solutions

Le module de résolution maintient un **historique circulaire** des solutions aux derniers instants $t_n, t_{n-1}, t_{n-2}, \ldots$ sous forme d'une liste chaînée circulaire. Cette structure permet à `ComputeExplicitTerms` d'accéder aux valeurs au pas précédent sans copie supplémentaire.

```
Avant avance temporelle :    tn-3 ← tn-2 ← tn-1 ← tn ←•

Après avance tn → tn+1 :     tn-2 ← tn-1 ← tn ← tn+1 ←•
                              (tn-3 est réutilisé pour stocker tn+1)
```

Le nombre de solutions conservées en mémoire simultanément est contrôlé par l'option :

```bash
bil -with "Monolithic N" fichier   # N = nombre de solutions en mémoire (défaut 2)
```

---

## 14. Comparaison avec d'autres patterns

| Concept | En C (Bil) | En C++ | En Python |
|---------|-----------|--------|-----------|
| Interface / contrat | `struct Model_t` avec pointeurs de fonctions | Classe abstraite avec méthodes `virtual` | Classe de base avec méthodes à redéfinir |
| Implémentation | Fonctions C libres assignées dans `SetModelProp` | Sous-classe qui override les méthodes virtuelles | Sous-classe qui redéfinit les méthodes |
| Registre | `ListOfModels.h` (X-macro, statique, compile-time) | Factory pattern, map de `std::string` → constructeur | Dictionnaire nom → classe |
| Dispatch | Appel via pointeur de fonction `model->computematrix(...)` | Dispatch via vtable C++ (automatique) | Dispatch via MRO Python |
| Découverte | Déclaration manuelle dans `ListOfModels.h` | Idem ou RTTI | Idem ou `importlib` dynamique |

Le choix du C avec X-macros garantit **zéro overhead à l'exécution** (pas de lookup dynamique) et une compatibilité maximale avec les compilateurs C89/C++17 utilisés par le projet.
