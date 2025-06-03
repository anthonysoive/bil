#ifndef MONOLITHIC_H
#define MONOLITHIC_H


struct Solutions_t;
struct DataSet_t;
struct Solver_t;
struct OutputFiles_t;


extern int (Monolithic_Initialize)(DataSet_t*,Solutions_t*) ;
extern int (Monolithic_Increment)(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*,double,double) ;
extern int (Monolithic_StepForward)(DataSet_t*,Solutions_t*,Solver_t*,double,double) ;
extern int (Monolithic_Iterate)(DataSet_t*,Solutions_t*,Solver_t*) ;

#endif
