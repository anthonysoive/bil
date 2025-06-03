#ifndef FUNCTION_H
#define FUNCTION_H


/* Forward declarations */
struct Function_t; //typedef struct Function_t     Function_t ;
struct DataFile_t;


extern Function_t*  (Function_New)          (const int) ;
extern void         (Function_Delete)       (void*) ;
extern int          (Function_Scan)         (Function_t*,DataFile_t*) ;
extern double       (Function_ComputeValue) (Function_t*,double) ;


#define Function_MaxLengthOfFileName       (200)
#define Function_MaxLengthOfTextLine       (500)


#define Function_GetNbOfPoints(FCT)      ((FCT)->n)
#define Function_GetXValue(FCT)          ((FCT)->t)
#define Function_GetFValue(FCT)          ((FCT)->f)



struct Function_t {           /* fonction du temps */
  int    n ;                  /* nombre de points */
  double* t ;                 /* temps */
  double* f ;                 /* valeurs f(t) */
} ;


/* Old notations which I try to eliminate little by little */
#define fonc_t           Function_t
#define fonction(a,b)    Function_ComputeValue(&(b),(a))

#endif
