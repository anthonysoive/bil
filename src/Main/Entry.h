#ifndef ENTRY_H
#define ENTRY_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Entry_t; //typedef struct Entry_t  Entry_t ;
struct Context_t;


extern Entry_t*    (Entry_Create)   (int,char**) ;
extern void        (Entry_Delete)   (void*) ;
extern int         (Entry_Execute)  (Entry_t*) ;
extern int         (Entry_Main)     (int,char**) ;


#define Entry_GetContext(E)          ((E)->context)

struct Entry_t {
  Context_t*     context ;
} ;

#ifdef __CPLUSPLUS
}
#endif
#endif
