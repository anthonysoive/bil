#ifndef BUFFER_H
#define BUFFER_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Buffer_t; //typedef struct Buffer_t     Buffer_t ;

#include <stdio.h>
extern Buffer_t*  (Buffer_Create)(size_t) ;
extern void       (Buffer_Delete)(void*) ;
extern void*      (Buffer_Allocate)(Buffer_t*,size_t) ;
extern void       (Buffer_Free)(Buffer_t*) ;
extern void       (Buffer_FreeFrom)(Buffer_t*,char*) ;


#define Buffer_GetBeginOfBuffer(BUF)    ((BUF)->buffer_begin)
#define Buffer_GetEndOfBuffer(BUF)      ((BUF)->buffer_end)
#define Buffer_GetHeadOfBuffer(BUF)     ((BUF)->head)
#define Buffer_GetTailOfBuffer(BUF)     ((BUF)->tail)
#define Buffer_GetSize(BUF)             ((BUF)->sz)
#define Buffer_GetAvailableSize(BUF)    ((BUF)->available_sz)


struct Buffer_t {           /* Circular buffer */
  void* buffer_begin ;      /* Data buffer */
  void* buffer_end ;        /* End of data buffer */
  void* head ;              /* Pointer to head */
  void* tail ;              /* Pointer to tail */
  size_t sz ;               /* Size of the buffer */
  size_t available_sz ;     /* Available size in the buffer */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
