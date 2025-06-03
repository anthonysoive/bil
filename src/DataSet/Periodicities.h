#ifndef PERIODICITIES_H
#define PERIODICITIES_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* Forward declaration */
struct Periodicities_t;       //typedef struct Periodicities_t      Periodicities_t ;
struct DataFile_t;
struct Mesh_t;
struct Graph_t;
struct Periodicity_t;


extern Periodicities_t* (Periodicities_Create)(DataFile_t*) ;
extern Periodicities_t* (Periodicities_New)(const int) ;
extern void             (Periodicities_Delete)(void*) ;
extern void             (Periodicities_EliminateMatrixRowColumnIndexes)(Mesh_t*) ;
extern void             (Periodicities_UpdateMatrixRowColumnIndexes)(Mesh_t*) ;
extern void             (Periodicities_UpdateGraph)(Mesh_t*,Graph_t*) ;


#define Periodicities_GetNbOfPeriodicities(PS) ((PS)->nbperiod)
#define Periodicities_GetPeriodicity(PS)       ((PS)->periodicity)


struct Periodicities_t {            /* Periodicities */
  unsigned int   nbperiod ;         /* Nb of periodicities */
  Periodicity_t* periodicity ;      /* Periodicity */
} ;



#ifdef __CPLUSPLUS
}
#endif
#endif
