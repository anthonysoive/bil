#ifndef NODESSOL_H
#define NODESSOL_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct NodesSol_t; //typedef struct NodesSol_t     NodesSol_t ;
struct Mesh_t;
struct Nodes_t;
struct NodeSol_t;


extern NodesSol_t*    (NodesSol_Create)(Mesh_t*) ;
extern void           (NodesSol_Delete)(void*) ;
extern void           (NodesSol_Copy)(NodesSol_t*,NodesSol_t*) ;
extern void           (NodesSol_CopySelectedSequentialUnknowns)(NodesSol_t*,NodesSol_t*,const int) ;
 
 

#define NodesSol_GetNodes(NSS)                 ((NSS)->nodes)
//#define NodesSol_GetNbOfDOF(NSS)               ((NSS)->NbOfDOF)
#define NodesSol_GetNbOfNodes(NSS)             ((NSS)->NbOfNodes)
#define NodesSol_GetNodeSol(NSS)               ((NSS)->nodesol)
//#define NodesSol_GetNodalValue(NSS)            ((NSS)->NodalValue)


#define NodesSol_SetNodes(NSS,A)\
        do {\
          NodesSol_GetNodes(NSS) = A;\
        } while(0)

#define NodesSol_SetNbOfNodes(NSS,A)\
        do {\
          NodesSol_GetNbOfNodes(NSS) = A;\
        } while(0)
        
#define NodesSol_SetNodeSol(NSS,A)\
        do {\
          NodesSol_GetNodeSol(NSS) = A;\
        } while(0)
        





/* complete the structure types by using the typedef */
struct NodesSol_t {           /* Nodal Solutions */
  Nodes_t* nodes ;
  size_t NbOfNodes ;
  NodeSol_t* nodesol ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
