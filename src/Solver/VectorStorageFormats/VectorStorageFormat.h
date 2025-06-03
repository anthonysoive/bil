#ifndef VECTORSTORAGEFORMAT_H
#define VECTORSTORAGEFORMAT_H

#ifdef __CPLUSPLUS
extern "C" {
#endif



enum VectorStorageFormat_e {        /* format of the vector to be stored */
  VectorStorageFormat_Array,        /* Array of doubles */
  VectorStorageFormat_PetscVec,     /* PETSc Vec format */
  VectorStorageFormat_NULL
} ;



/* Forward declarations */
typedef enum VectorStorageFormat_e VectorStorageFormat_e ;

struct VectorStorageFormat_t; //typedef struct VectorStorageFormat_t  VectorStorageFormat_t ;
struct Options_t;


extern VectorStorageFormat_t* (VectorStorageFormat_Create)(Options_t*) ;
extern void                   (VectorStorageFormat_Delete)(void*) ;


/** The getters */
#define VectorStorageFormat_GetType(MSF)              ((MSF)->type)
#define VectorStorageFormat_GetOptions(MSF)           ((MSF)->options)



#define VectorStorageFormat_Is(MSF,KEY) \
        (VectorStorageFormat_GetType(MSF) == Utils_CAT(VectorStorageFormat_,KEY))

#define VectorStorageFormat_Type(KEY) \
        ((VectorStorageFormat_e) Utils_CAT(VectorStorageFormat_,KEY))



/* complete the structure types by using the typedef */

struct VectorStorageFormat_t {
  VectorStorageFormat_e type ;
  Options_t* options ;
} ;


#ifdef __CPLUSPLUS
}
#endif

#include "Utils.h"
#endif
