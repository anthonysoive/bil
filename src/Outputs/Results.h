#ifndef RESULTS_H
#define RESULTS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Results_t; //typedef struct Results_t      Results_t ;
struct Result_t;

extern Results_t* (Results_Create)(int) ;
extern void       (Results_Delete)(void*) ;


#define Results_GetNbOfResults(results)      ((results)->nbofresults)
#define Results_GetResult(results)           ((results)->result)



struct Results_t {            /* Results */
  unsigned int nbofresults ;  /* nb of results */
  Result_t* result ;          /* result */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
