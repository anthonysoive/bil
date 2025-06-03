#ifndef PERIODICITY_H
#define PERIODICITY_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Periodicity_t;         //typedef struct Periodicity_t        Periodicity_t ;
struct DataFile_t;

extern Periodicity_t* (Periodicity_New)(void) ;
extern void           (Periodicity_Delete)(void*) ;
extern void           (Periodicity_Scan)(Periodicity_t*,DataFile_t*) ;


#define Periodicity_MaxLengthOfRegionName      Region_MaxLengthOfRegionName

#define Periodicity_GetMasterRegion(P)         ((P)->MasterRegion)
#define Periodicity_GetSlaveRegion(P)          ((P)->SlaveRegion)
#define Periodicity_GetMasterRegionName(P)     ((P)->MasterRegionName)
#define Periodicity_GetSlaveRegionName(P)      ((P)->SlaveRegionName)
#define Periodicity_GetPeriodVector(P)         ((P)->PeriodVector)



struct Periodicity_t {
  int     MasterRegion ;               /* Master region index */
  int     SlaveRegion ;                /* Slave region index */
  char*   MasterRegionName ;
  char*   SlaveRegionName ;
  double* PeriodVector ;
} ;


#ifdef __CPLUSPLUS
}
#endif

/* Needs for the macros */
#include "Region.h"
#endif
