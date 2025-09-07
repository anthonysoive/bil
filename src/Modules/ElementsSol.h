#ifndef ELEMENTSSOL_H
#define ELEMENTSSOL_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* class-like structures "ElementsSol_t" and attributes */

/* Forward declarations */
struct ElementsSol_t; //typedef struct ElementsSol_t  ElementsSol_t ;
struct ElementSol_t;
struct Mesh_t;


extern ElementsSol_t*   (ElementsSol_Create)(Mesh_t*) ;
extern void             (ElementsSol_Delete)(void*) ;
extern void             (ElementsSol_AllocateMemoryForImplicitTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_AllocateMemoryForExplicitTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_AllocateMemoryForConstantTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_Copy)(ElementsSol_t*,ElementsSol_t*) ;


#define ElementsSol_GetNbOfElements(ESS)         ((ESS)->NbOfElements)
#define ElementsSol_GetElementSol(ESS)           ((ESS)->elementsol)

#define ElementsSol_SetNbOfElements(ESS,A) \
        do {\
          ElementsSol_GetNbOfElements(ESS) = A;\
        } while(0)
        
#define ElementsSol_SetElementSol(ESS,A) \
        do {\
          ElementsSol_GetElementSol(ESS) = A;\
        } while(0)




/* Delete data */
#define ElementsSol_DeleteExplicitGenericData(ESS) \
        do { \
          size_t NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          for(size_t ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_DeleteExplicitGenericData(elementsol + ElementsSol_i) ; \
          } \
        } while(0)
        
        
#define ElementsSol_DeleteConstantGenericData(ESS) \
        do { \
          size_t NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          for(size_t ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_DeleteConstantGenericData(elementsol + ElementsSol_i) ; \
          } \
        } while(0)
        
        

/* Share data */
#define ElementsSol_ShareExplicitGenericData(ESS_DEST,ESS_SRC) \
        do { \
          size_t NbOfElements = ElementsSol_GetNbOfElements(ESS_SRC) ; \
          ElementSol_t* elementsol_s = ElementsSol_GetElementSol(ESS_SRC) ; \
          ElementSol_t* elementsol_d = ElementsSol_GetElementSol(ESS_DEST) ; \
          for(size_t ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            GenericData_t* gdat = ElementSol_GetExplicitGenericData(elementsol_s + ElementsSol_i) ; \
            ElementSol_GetExplicitGenericData(elementsol_d + ElementsSol_i) = gdat ; \
          } \
        } while(0)


#define ElementsSol_ShareConstantGenericData(ESS_DEST,ESS_SRC) \
        do { \
          size_t NbOfElements = ElementsSol_GetNbOfElements(ESS_SRC) ; \
          ElementSol_t* elementsol_s = ElementsSol_GetElementSol(ESS_SRC) ; \
          ElementSol_t* elementsol_d = ElementsSol_GetElementSol(ESS_DEST) ; \
          for(size_t ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            GenericData_t* gdat = ElementSol_GetConstantGenericData(elementsol_s + ElementsSol_i) ; \
            ElementSol_GetConstantGenericData(elementsol_d + ElementsSol_i) = gdat ; \
          } \
        } while(0)



/* Set pointers to null */
#define ElementsSol_DiscardExplicitGenericData(ESS) \
        do { \
          size_t NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          for(size_t ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_GetExplicitGenericData(elementsol + ElementsSol_i) = NULL ; \
          } \
        } while(0)


#define ElementsSol_DiscardConstantGenericData(ESS) \
        do { \
          size_t NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          for(size_t ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_GetConstantGenericData(elementsol + ElementsSol_i) = NULL ; \
          } \
        } while(0)






struct ElementsSol_t {
  size_t NbOfElements ;
  ElementSol_t* elementsol ;
} ;


/* For the macros */
#include "ElementSol.h"
#include "GenericData.h"

#ifdef __CPLUSPLUS
}
#endif
#endif
