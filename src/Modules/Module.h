#ifndef MODULE_H
#define MODULE_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Module_t; //typedef struct Module_t        Module_t ;
struct DataSet_t;
struct OutputFiles_t;
struct Solutions_t;
struct Solver_t;


/* Declaration of Macros, Methods and Structures */

extern Module_t*  (Module_New)(void) ;
extern void       (Module_Delete)(void*) ;
extern Module_t*  (Module_Initialize)(Module_t*,const char*) ;


#if 0
//extern int    Module_StoreCurrentSolution(Mesh_t*,double,char*) ;
//extern int    Module_LoadCurrentSolution(Mesh_t*,double*,char*) ;
//extern void   Module_InitializeCurrentPointers(Mesh_t*,Solution_t*) ;
//extern void   Module_SetCurrentUnknownsWithBoundaryConditions(Mesh_t*,BConds_t*,double) ;
//extern void   Module_UpdateCurrentUnknowns(Mesh_t*,Solver_t*) ;
#endif


#define Module_MaxLengthOfKeyWord        (30)
#define Module_MaxLengthOfFileName       (60)
#define Module_MaxLengthOfShortTitle     (80)
#define Module_MaxLengthOfAuthorNames    (80)


#define Module_GetCodeNameOfModule(MOD)   ((MOD)->codename)
#define Module_GetShortTitle(MOD)         ((MOD)->shorttitle)
#define Module_GetNameOfAuthors(MOD)      ((MOD)->authors)
#define Module_GetSetModuleProp(MOD)      ((MOD)->setmoduleprop)
#define Module_GetComputeProblem(MOD)     ((MOD)->computeproblem)
#define Module_GetSolveProblem(MOD)       ((MOD)->solveproblem)
#define Module_GetIncrement(MOD)          ((MOD)->increment)
//#define Module_GetStepForward(MOD)        ((MOD)->stepforward)
#define Module_GetInitializeProblem(MOD)  ((MOD)->initializeproblem)
#define Module_GetSequentialIndex(MOD)    ((MOD)->sequentialindex)
#define Module_GetNbOfSequences(MOD)      ((MOD)->nbofsequences)


#define Module_SetCodeNameOfModule(MOD,A)\
        do {\
          Module_GetCodeNameOfModule(MOD) = A;\
        } while(0)
        
#define Module_SetShortTitle(MOD,A)\
        do {\
          Module_GetShortTitle(MOD) = A;\
        } while(0)
        
#define Module_SetNameOfAuthors(MOD,A)\
        do {\
          Module_GetNameOfAuthors(MOD) = A;\
        } while(0)
        
#define Module_SetSetModuleProp(MOD,A)\
        do {\
          Module_GetSetModuleProp(MOD) = A;\
        } while(0)
        
#define Module_SetComputeProblem(MOD,A)\
        do {\
          Module_GetComputeProblem(MOD) = A;\
        } while(0)
        
#define Module_SetSolveProblem(MOD,A)\
        do {\
          Module_GetSolveProblem(MOD) = A;\
        } while(0)
        
#define Module_SetIncrement(MOD,A)\
        do {\
          Module_GetIncrement(MOD) = A;\
        } while(0)
        
#define Module_SetInitializeProblem(MOD,A)\
        do {\
          Module_GetInitializeProblem(MOD) = A;\
        } while(0)
        
#define Module_SetSequentialIndex(MOD,A)\
        do {\
          Module_GetSequentialIndex(MOD) = A;\
        } while(0)
        
#define Module_SetNbOfSequences(MOD,A)\
        do {\
          Module_GetNbOfSequences(MOD) = A;\
        } while(0)
        



/* Copy operations */
#define Module_CopyCodeNameOfModule(MOD,codename) \
        (strcpy(Module_GetCodeNameOfModule(MOD),codename))

#define Module_CopyShortTitle(MOD,title) \
        (strcpy(Module_GetShortTitle(MOD),title))

#define Module_CopyNameOfAuthors(MOD,authors) \
        (strcpy(Module_GetNameOfAuthors(MOD),authors))


/* Short hands */
#define Module_SetModuleProp(MOD) \
        (Module_GetSetModuleProp(MOD)) ? \
        Module_GetSetModuleProp(MOD)(MOD) : -1

#define Module_ComputeProblem(MOD,...) \
        (Module_GetComputeProblem(MOD)) ? \
        Module_GetComputeProblem(MOD)(__VA_ARGS__) : -1

#define Module_SolveProblem(MOD,...) \
        (Module_GetSolveProblem(MOD)) ? \
        Module_GetSolveProblem(MOD)(__VA_ARGS__) : -1

#define Module_Increment(MOD,...) \
        (Module_GetIncrement(MOD)) ? \
        Module_GetIncrement(MOD)(__VA_ARGS__) : -1

#define Module_InitializeProblem(MOD,...) \
        (Module_GetInitializeProblem(MOD)) ? \
        Module_GetInitializeProblem(MOD)(__VA_ARGS__) : -1


/*  Typedef names of Methods */
typedef int   Module_SetModuleProp_t(Module_t*) ;
typedef int   Module_ComputeProblem_t(DataSet_t*) ;
typedef int   Module_SolveProblem_t(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*) ;
typedef int   Module_StepForward_t(DataSet_t*,Solutions_t*,Solver_t*,double,double) ;
typedef int   Module_Increment_t(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*,double,double) ;
typedef int   Module_InitializeProblem_t(DataSet_t*,Solutions_t*) ;


struct Module_t {              /* module */
  Module_SetModuleProp_t*      setmoduleprop ;
  Module_ComputeProblem_t*     computeproblem ;
  Module_SolveProblem_t*       solveproblem ;
  Module_Increment_t*          increment ;
  //Module_StepForward_t*        stepforward ;
  Module_InitializeProblem_t*  initializeproblem ;
  int nbofsequences ;
  int sequentialindex ;
  
  char*  codename ;          /* Code name of the module */
  char*  authors ;           /* Authors of this module */
  char*  shorttitle ;        /* Short title of the module */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
