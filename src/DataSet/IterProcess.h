#ifndef ITERPROCESS_H
#define ITERPROCESS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct IterProcess_t; //typedef struct IterProcess_t  IterProcess_t ;
struct DataFile_t;
struct ObVals_t;
struct Node_t;
struct Nodes_t;
struct Solver_t;


extern IterProcess_t*  (IterProcess_New)(void) ;
extern IterProcess_t*  (IterProcess_Create)(DataFile_t*,ObVals_t*) ;
extern void            (IterProcess_Delete)(void*) ;
extern int             (IterProcess_SetCurrentError)(IterProcess_t*,Nodes_t*,Solver_t*) ;
extern void            (IterProcess_PrintCurrentError)(IterProcess_t*) ;


#define IterProcess_GetNbOfIterations(IPR)           ((IPR)->niter)
#define IterProcess_GetNbOfRepetitions(IPR)          ((IPR)->nrecom)
#define IterProcess_GetTolerance(IPR)                ((IPR)->tol)
#define IterProcess_GetRepetitionIndex(IPR)          ((IPR)->irecom)
#define IterProcess_GetIterationIndex(IPR)           ((IPR)->iter)
#define IterProcess_GetError(IPR)                    ((IPR)->error)
#define IterProcess_GetObValIndexOfCurrentError(IPR) ((IPR)->obvalindex)
#define IterProcess_GetNodeIndexOfCurrentError(IPR)  ((IPR)->nodeindex)
#define IterProcess_GetNodeOfCurrentError(IPR)       ((IPR)->node)
#define IterProcess_GetObVals(IPR)                   ((IPR)->obvals)


#define IterProcess_SetNbOfIterations(IPR,A) \
        do {\
          IterProcess_GetNbOfIterations(IPR) = A;\
        } while(0)
        
#define IterProcess_SetNbOfRepetitions(IPR,A) \
        do {\
          IterProcess_GetNbOfRepetitions(IPR) = A;\
        } while(0)
        
#define IterProcess_SetTolerance(IPR,A) \
        do {\
          IterProcess_GetTolerance(IPR) = A;\
        } while(0)
        
#define IterProcess_SetRepetitionIndex(IPR,A) \
        do {\
          IterProcess_GetRepetitionIndex(IPR) = A;\
        } while(0)
        
#define IterProcess_SetIterationIndex(IPR,A) \
        do {\
          IterProcess_GetIterationIndex(IPR) = A;\
        } while(0)
        
#define IterProcess_SetError(IPR,A) \
        do {\
          IterProcess_GetError(IPR) = A;\
        } while(0)
        
#define IterProcess_SetObValIndexOfCurrentError(IPR,A) \
        do {\
          IterProcess_GetObValIndexOfCurrentError(IPR) = A;\
        } while(0)
        
#define IterProcess_SetNodeIndexOfCurrentError(IPR,A) \
        do {\
          IterProcess_GetNodeIndexOfCurrentError(IPR) = A;\
        } while(0)
        
#define IterProcess_SetNodeOfCurrentError(IPR,A) \
        do {\
          IterProcess_GetNodeOfCurrentError(IPR) = A;\
        } while(0)
        
#define IterProcess_SetObVals(IPR,A) \
        do {\
          IterProcess_GetObVals(IPR) = A;\
        } while(0)
        


#define IterProcess_GetObVal(IPR) \
        ObVals_GetObVal(IterProcess_GetObVals(IPR))
        

/* Operations on iterations */
#define IterProcess_IncrementIterationIndex(IPR) \
        (IterProcess_GetIterationIndex(IPR)++)

#define IterProcess_LastIterationIsNotReached(IPR) \
        (IterProcess_GetIterationIndex(IPR) < IterProcess_GetNbOfIterations(IPR))

#define IterProcess_InitializeIterations(IPR) \
        (IterProcess_GetIterationIndex(IPR) = 0)


/* Operations on repetitions */
#define IterProcess_IncrementRepetitionIndex(IPR) \
       (IterProcess_GetRepetitionIndex(IPR)++)

#define IterProcess_LastRepetitionIsNotReached(IPR) \
        (IterProcess_GetRepetitionIndex(IPR) < IterProcess_GetNbOfRepetitions(IPR))

#define IterProcess_InitializeRepetitions(IPR) \
        (IterProcess_GetRepetitionIndex(IPR) = 0)


/* Operations on convergence */
#define IterProcess_ConvergenceIsAttained(IPR) \
        (IterProcess_GetError(IPR) < IterProcess_GetTolerance(IPR))

#define IterProcess_ConvergenceIsNotAttained(IPR) \
        (!IterProcess_ConvergenceIsAttained(IPR))


/* Error on which unknown? */
#define IterProcess_GetNameOfTheCurrentError(IPR) \
        (ObVal_GetNameOfUnknown(IterProcess_GetObVal(IPR) + IterProcess_GetObValIndexOfCurrentError(IPR)))


struct IterProcess_t {        /* Iterative process */
  int    niter ;              /* Max nb of iterations */
  int    iter ;               /* Current iteration index */
  int    nrecom ;             /* Max nb of repetitions */
  int    irecom ;             /* Current repetition index */
  double tol ;                /* Tolerance */
  double error ;              /* Current error */
  int    obvalindex ;         /* Objective value index pertaining to the greatest error */
  size_t nodeindex ;          /* Node index pertaining to the greatest error */
  Node_t* node ;              /* Node pertaining to the greatest error */
  ObVals_t* obvals ;          /* Objective variations */
} ;


#ifdef __CPLUSPLUS
}
#endif

#include "ObVals.h"
#include "ObVal.h"
#endif
