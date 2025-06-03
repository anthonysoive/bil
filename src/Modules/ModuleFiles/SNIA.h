#ifndef SNIA_H
#define SNIA_H


struct Solutions_t;
struct DataSet_t;
struct Solver_t;
struct OutputFiles_t;


extern int (SNIA_Initialize)(DataSet_t*,Solutions_t*) ;
extern int (SNIA_Increment)(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*,double,double) ;
extern int (SNIA_StepForward)(DataSet_t*,Solutions_t*,Solver_t*,double,double) ;
extern int (SNIA_Iterate)(DataSet_t*,Solutions_t*,Solver_t*) ;


#endif
