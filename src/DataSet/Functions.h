#ifndef FUNCTIONS_H
#define FUNCTIONS_H


/* Forward declarations */
struct Functions_t; //typedef struct Functions_t    Functions_t ;
struct Function_t;
struct DataFile_t;


extern Functions_t* (Functions_New)     (const int) ;
extern Functions_t* (Functions_Create)  (DataFile_t*) ;
extern void         (Functions_Delete)  (void*) ;


#define Functions_GetNbOfFunctions(FCTS)      ((FCTS)->n_fn)
#define Functions_GetFunction(FCTS)           ((FCTS)->fn)




struct Functions_t {          /* time functions */
  int n_fn ;         /* nb of functions */
  Function_t* fn ;            /* function */
} ;

#endif
