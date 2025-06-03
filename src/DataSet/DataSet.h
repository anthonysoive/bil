#ifndef DATASET_H
#define DATASET_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct DataSet_t; //typedef struct DataSet_t      DataSet_t ;
struct Options_t;
struct Units_t;
struct DataFile_t;
struct Geometry_t;
struct Mesh_t;
struct Materials_t;
struct Dates_t;
struct Points_t;
struct IConds_t;
struct BConds_t;
struct Loads_t;
struct Functions_t;
struct Fields_t;
struct ObVals_t;
struct TimeStep_t;
struct IterProcess_t;
struct Module_t;
struct Models_t;


extern DataSet_t*  (DataSet_Create)    (char*,Options_t*) ;
extern void        (DataSet_Delete)    (void*) ;
//extern DataSet_t*  (DataSet_Create1)   (char*,Options_t*) ;
extern void        (DataSet_PrintData) (DataSet_t*,char*) ;



#define DataSet_New() \
        (DataSet_t*) Mry_New(DataSet_t)




#define DataSet_GetUnits(DS)          ((DS)->units)
#define DataSet_GetDataFile(DS)       ((DS)->datafile)
#define DataSet_GetGeometry(DS)       ((DS)->geometry)
#define DataSet_GetMesh(DS)           ((DS)->mesh)
#define DataSet_GetMaterials(DS)      ((DS)->materials)
#define DataSet_GetDates(DS)          ((DS)->dates)
#define DataSet_GetPoints(DS)         ((DS)->points)
#define DataSet_GetIConds(DS)         ((DS)->iconds)
#define DataSet_GetBConds(DS)         ((DS)->bconds)
#define DataSet_GetLoads(DS)          ((DS)->loads)
#define DataSet_GetFunctions(DS)      ((DS)->functions)
#define DataSet_GetFields(DS)         ((DS)->fields)
/* #define DataSet_GetIntFcts(DS)        ((DS)->intfcts) */
#define DataSet_GetObVals(DS)         ((DS)->obvals)
#define DataSet_GetModels(DS)         ((DS)->models)
#define DataSet_GetTimeStep(DS)       ((DS)->timestep)
#define DataSet_GetIterProcess(DS)    ((DS)->iterprocess)
#define DataSet_GetOptions(DS)        ((DS)->options)
//#define DataSet_GetModules(DS)        ((DS)->modules)
#define DataSet_GetModule(DS)         ((DS)->module)





#define DataSet_GetSequentialIndex(DS) \
        Module_GetSequentialIndex(DataSet_GetModule(DS))
        
        
#define DataSet_GetNbOfSequences(DS) \
        Module_GetNbOfSequences(DataSet_GetModule(DS))



struct DataSet_t {               /* set of data for the problem to work out */
  Units_t*       units ;         /* Units */
  DataFile_t*    datafile ;      /* data file */
  Geometry_t*    geometry ;      /* Geometry */
  Mesh_t*        mesh ;          /* Mesh */
  Materials_t*   materials ;     /* Materials */
  Dates_t*       dates ;         /* Dates */
  Points_t*      points ;        /* Points */
  IConds_t*      iconds ;        /* Initial Conditions */
  BConds_t*      bconds ;        /* Boundary Conditions */
  Loads_t*       loads ;         /* Loadings */
  Functions_t*   functions ;     /* Functions */
  Fields_t*      fields ;        /* Fields */
  /* IntFcts_t*     intfcts ; */       /* Interpolation Functions */
  ObVals_t*      obvals ;        /* Objective Values */
  Models_t*      models ;        /* Models */
  TimeStep_t*    timestep ;      /* time step managing */
  IterProcess_t* iterprocess ;   /* iterative process */
  Options_t*     options ;       /* options */
  //Modules_t*     modules ;       /* modules */
  Module_t*      module ;        /* module */
} ;



#ifdef __CPLUSPLUS
}
#endif

/* For the macros */
#include "Mry.h"
#include "Module.h"
#endif
