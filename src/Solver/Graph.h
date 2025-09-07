#ifndef GRAPH_H
#define GRAPH_H

#ifdef __CPLUSPLUS
extern "C" {
#endif



/* Forward declarations */
struct Graph_t; //typedef struct Graph_t Graph_t ;
struct AdjacencyList_t;


extern Graph_t*  (Graph_Create)(size_t,unsigned short int*) ;
extern void      (Graph_Delete)(void*) ;


#define Graph_GetNbOfVertices(G)              ((G)->nvertices)
#define Graph_GetNbOfEdges(G)                 ((G)->nedges)
#define Graph_GetAdjacencyList(G)             ((G)->adj)


#define Graph_GetDegreeOfVertex(G,I) \
        AdjacencyList_GetNbOfNeighbors(Graph_GetAdjacencyList(G) + I)
        
#define Graph_GetMaxDegreeOfVertex(G,I) \
        AdjacencyList_GetMaxNbOfNeighbors(Graph_GetAdjacencyList(G) + I)
        
#define Graph_GetNeighborOfVertex(G,I) \
        AdjacencyList_GetNeighbor(Graph_GetAdjacencyList(G) + I)

#define Graph_MaxDegreeNotAttainedAtVertex(G,I) \
        AdjacencyList_MaxNbOfNeighborsIsNotAttained(Graph_GetAdjacencyList(G) + I)


#define Graph_AddEdge(G,I,J) \
        do {\
          AdjacencyList_AddNeighbor(Graph_GetAdjacencyList(G) + I,J);\
          AdjacencyList_AddNeighbor(Graph_GetAdjacencyList(G) + J,I);\
        } while(0)


#define Graph_UpdateTheNbOfEdges(G) \
        do { \
          size_t nvert = Graph_GetNbOfVertices(G) ; \
          size_t nedges = 0 ; \
          for(size_t i = 0 ; i < nvert ; i++) { \
            nedges += Graph_GetDegreeOfVertex(G,i) ; \
          } \
          Graph_GetNbOfEdges(G) = nedges/2 ; \
        } while(0)


struct Graph_t {              /* Graph */
  size_t  nvertices ;   /* Nb of vertices */
  size_t  nedges ;      /* Nb of edges */
  AdjacencyList_t* adj ;      /* Adjacency list */
} ;


#ifdef __CPLUSPLUS
}
#endif

#include "AdjacencyList.h"
#endif
