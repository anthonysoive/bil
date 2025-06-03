#ifndef ICONDS_H
#define ICONDS_H


/* Forward declarations */
struct IConds_t; //typedef struct IConds_t       IConds_t ;
struct DataFile_t;
struct Functions_t;
struct Mesh_t;
struct Fields_t;
struct ICond_t;


//extern IConds_t* IConds_Create(DataFile_t*,Mesh_t*,Fields_t*) ;
extern IConds_t* (IConds_Create)(DataFile_t*,Fields_t*,Functions_t*) ;
extern void      (IConds_Delete)(void*) ;
extern void      (IConds_AssignInitialConditions)(IConds_t*,Mesh_t*,double) ;


#define IConds_MaxLengthOfKeyWord        (30)
#define IConds_MaxLengthOfFileName       (60)


#define IConds_GetNbOfIConds(ICS)              ((ICS)->n_ic)
#define IConds_GetICond(ICS)                   ((ICS)->ic)
#define IConds_GetFileNameOfNodalValues(ICS)   ((ICS)->file)


struct IConds_t {             /* Initial Conditions */
  char*  file ;               /* Name of file containing nodal values */
  unsigned int n_ic ;         /* Nb of IC */
  ICond_t* ic ;               /* Initial condition */
} ;


#endif
