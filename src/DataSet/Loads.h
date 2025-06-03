#ifndef LOADS_H
#define LOADS_H


/* Forward declarations */
struct Loads_t; //typedef struct Loads_t        Loads_t ;
struct DataFile_t;
struct Fields_t;
struct Functions_t;
struct Load_t;


extern Loads_t* (Loads_New)(const int) ;
extern Loads_t* (Loads_Create)(DataFile_t*,Fields_t*,Functions_t*) ;
extern void     (Loads_Delete)(void*) ;

#define Loads_GetNbOfLoads(LOADS)        ((LOADS)->n_cg)
#define Loads_GetLoad(LOADS)             ((LOADS)->cg)


struct Loads_t {              /* loadings */
  unsigned int n_cg ;         /* nb */
  Load_t* cg ;                /* loading */
} ;


#endif
