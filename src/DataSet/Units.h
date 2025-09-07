#ifndef UNITS_H
#define UNITS_H


/* Forward declarations */
struct Units_t; //typedef struct Units_t        Units_t ;
struct DataFile_t;
struct Unit_t;
 
 
extern Units_t* (Units_New)     (void) ;
extern Units_t* (Units_Create)  (DataFile_t*) ;
extern void     (Units_Delete)  (void*) ;


#define Units_MaxNbOfUnits   (7)


#define Units_GetNbOfUnits(U)         ((U)->n_units)
#define Units_GetUnit(U)              ((U)->unit)


struct Units_t {              /* units */
  int n_units ;      /* nb of units */
  Unit_t*  unit ;             /* Unit */
} ;

#endif
