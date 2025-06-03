#ifndef BUFFERS_H
#define BUFFERS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Buffers_t; //typedef struct Buffers_t     Buffers_t ;
struct Buffer_t;


#include <stdio.h>
extern Buffers_t*  (Buffers_Create)(size_t) ;
extern void        (Buffers_Delete)(void*) ;


#define Buffers_GetNbOfBuffers(BFS)   ((BFS)->nbofbuffers)
#define Buffers_GetBuffer(BFS)        ((BFS)->buffer)


#define Buffers_MaxNbOfBuffers   SharedMS_MaxNbOfThreads


#define Buffers_GetBufferOfCurrentThread(BFS) \
        ((SharedMS_CurrentThreadId < Buffers_GetNbOfBuffers(BFS)) ? \
        (Buffers_GetBuffer(BFS) + SharedMS_CurrentThreadId) : \
        NULL)



struct Buffers_t {
  int nbofbuffers ;
  Buffer_t* buffer ;
} ;


#ifdef __CPLUSPLUS
}
#endif

/* For the macros */
#include "SharedMS.h"
#endif
