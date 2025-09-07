#ifndef PETSCAIJFORMAT_H
#define PETSCAIJFORMAT_H


/* Forward declarations */
struct PetscAIJFormat_t; //typedef struct PetscAIJFormat_t PetscAIJFormat_t ;
struct Mesh_t;


extern PetscAIJFormat_t* (PetscAIJFormat_Create)(Mesh_t*,const int) ;
extern void              (PetscAIJFormat_Delete)(void*) ;
extern size_t            (PetscAIJFormat_AssembleElementMatrix)(PetscAIJFormat_t*,double*,int*,int*,int) ;
extern void              (PetscAIJFormat_PrintMatrix)(PetscAIJFormat_t*,const char*) ;



#if 1
/** The getters */
#define PetscAIJFormat_GetNbOfNonZeroValues(a)              ((a)->nnz)
#define PetscAIJFormat_GetNonZeroValue(a)                   ((a)->nzval)
#define PetscAIJFormat_GetIndexOfRow(a)                     ((a)->idwrow)
#define PetscAIJFormat_GetIndexOfColumn(a)                  ((a)->idxcol)
#define PetscAIJFormat_GetStorage(a)                        ((a)->store)
#endif



/* complete the structure types by using the typedef */
#if 1
struct PetscAIJFormat_t {
  size_t    nnz ;          /* nb of non zero values */
  double* nzval ;       /* Non zero values */
  int* idxrow ;         /* Indices of the rows */
  int* idxcol ;         /* Indices of the columns */
  void* store;          /* pointer to the actual storage of the matrix */
} ;
#endif



#endif
