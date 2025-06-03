#ifndef ICOND_H
#define ICOND_H


/* Forward declarations */
struct ICond_t; //typedef struct ICond_t        ICond_t ;
struct DataFile_t;
struct Functions_t;
struct Function_t;
struct Fields_t;
struct Field_t;


extern ICond_t* (ICond_New)    (void) ;
extern void     (ICond_Delete) (void*) ;
extern void     (ICond_Scan)   (ICond_t*,DataFile_t*) ;



#define ICond_MaxLengthOfRegionName      Region_MaxLengthOfRegionName
#define ICond_MaxLengthOfKeyWord        (30)
#define ICond_MaxLengthOfFileName       (60)


#define ICond_GetRegionTag(IC)               ((IC)->RegionTag)
#define ICond_GetRegionName(IC)              ((IC)->RegionName)
#define ICond_GetNameOfUnknown(IC)           ((IC)->NameOfUnknown)
#define ICond_GetFunction(IC)                ((IC)->Function)
#define ICond_GetField(IC)                   ((IC)->Field)
#define ICond_GetFileNameOfNodalValues(IC)   ((IC)->FileNameOfNodalValues)
#define ICond_GetEntryInDataFile(IC)         ((IC)->EntryInDataFile)
#define ICond_GetFunctionIndex(IC)           ((IC)->FunctionIndex)
#define ICond_GetFieldIndex(IC)              ((IC)->FieldIndex)
#define ICond_GetFunctions(IC)               ((IC)->Functions)
#define ICond_GetFields(IC)                  ((IC)->Fields)





struct ICond_t {              /* Initial condition */
  Function_t* Function ;
  Field_t* Field ;
  Functions_t* Functions ;
  Fields_t* Fields ;
  char*  EntryInDataFile ;       /* Entry point in the string of data */
  int    RegionTag ;
  char*  RegionName ;
  char*  NameOfUnknown ;
  int    FunctionIndex ;
  int    FieldIndex ;
  char*  FileNameOfNodalValues ;
} ;


#include "Region.h"
#endif
