#ifndef BCONDS_H
#define BCONDS_H


/* Forward declarations */
struct BConds_t; //typedef struct BConds_t       BConds_t ;
struct DataFile_t;
struct Functions_t;
struct Fields_t;
struct Mesh_t;
struct BCond_t;


extern BConds_t* (BConds_New)(const int) ;
extern BConds_t* (BConds_Create)(DataFile_t*,Fields_t*,Functions_t*) ;
extern void      (BConds_Delete)(void*) ;
extern void      (BConds_EliminateMatrixRowColumnIndexes)(BConds_t*,Mesh_t*) ;
extern void      (BConds_AssignBoundaryConditions)(BConds_t*,Mesh_t*,double) ;


#define BConds_GetNbOfBConds(BCS)        ((BCS)->n_cl)
#define BConds_GetBCond(BCS)             ((BCS)->cl)


struct BConds_t {             /* boundary conditions */
  int n_cl ;         /* nb */
  BCond_t* cl ;               /* boundary condition */
} ;

#endif
