#ifndef UNIT_H
#define UNIT_H



/* Forward declarations */
struct Unit_t; //typedef struct Unit_t         Unit_t ;
struct DataFile_t;


extern Unit_t* (Unit_New)     (void) ;
extern void    (Unit_Delete)  (void*) ;
extern int     (Unit_Scan)    (Unit_t*,DataFile_t*) ;



#define Unit_MaxLengthOfKeyWord   (30)



#define Unit_GetName(U)         ((U)->name)
#define Unit_GetValue(U)        ((U)->value)


struct Unit_t {           /* unit */
  char*  name ;           /* Its name */
  double value ;          /* Its value */
} ;


#endif
