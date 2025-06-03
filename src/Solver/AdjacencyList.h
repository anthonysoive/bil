#ifndef ADJACENCYLIST_H
#define ADJACENCYLIST_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct AdjacencyList_t; //typedef struct AdjacencyList_t AdjacencyList_t ;


extern AdjacencyList_t* (AdjacencyList_Create)(int,int*) ;
extern void             (AdjacencyList_Delete)(void*) ;


#define AdjacencyList_GetNbOfNeighbors(adj)       ((adj)->ndest)
#define AdjacencyList_GetNeighbor(adj)            ((adj)->dest)



struct AdjacencyList_t {      /* Format */
  unsigned int  ndest ;       /* Nb of neighbors */
  int* dest ;                 /* Neighbors */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
