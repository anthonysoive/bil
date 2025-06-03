#ifndef NODESOL_H
#define NODESOL_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct NodeSol_t; //typedef struct NodeSol_t      NodeSol_t ;



extern NodeSol_t* (NodeSol_Create)(const int) ;
extern void       (NodeSol_Delete)(void*) ;
extern void       (NodeSol_Copy)(NodeSol_t*,NodeSol_t*) ;


#define NodeSol_GetNbOfUnknowns(NS)       ((NS)->nu)
#define NodeSol_GetUnknown(NS)            ((NS)->u)
//#define NodeSol_GetPreviousNodeSol(NS)    ((NS)->prev)
//#define NodeSol_GetNextNodeSol(NS)        ((NS)->next)


#define NodeSol_SetNbOfUnknowns(NS,A) \
        do {\
          NodeSol_GetNbOfUnknowns(NS) = A;\
        } while(0)
        
#define NodeSol_SetUnknown(NS,A) \
        do {\
          NodeSol_GetUnknown(NS) = A;\
        } while(0)



struct NodeSol_t {            /* Nodal Solutions */
  unsigned int nu ;     /* Nb of unknowns */
  double* u ;                 /* Nodal Unknowns */
  //NodeSol_t* prev ;           /* Previous Nodal Solutions */
  //NodeSol_t* next ;           /* Next Nodal Solutions */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
