#ifndef SOLVERS_H
#define SOLVERS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Solvers_t; //typedef struct Solvers_t       Solvers_t ;
struct Mesh_t;
struct Options_t;
struct Solver_t;


extern Solvers_t*  (Solvers_Create)(Mesh_t*,Options_t*,const int) ;
extern void        (Solvers_Delete)(void*) ;


#define Solvers_GetNbOfSolvers(SV)        ((SV)->nbofsolvers)
#define Solvers_GetSolver(SV)             ((SV)->solver)


/* complete the structure types by using the typedef */
struct Solvers_t {
  int nbofsolvers ;
  Solver_t* solver ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
