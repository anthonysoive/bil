#ifndef POSFILESFORGMSH_H
#define POSFILESFORGMSH_H


/* Forward declarations */
struct PosFilesForGMSH_t; //typedef struct PosFilesForGMSH_t PosFilesForGMSH_t ;
struct DataSet_t;
struct OutputFiles_t;

extern PosFilesForGMSH_t*  (PosFilesForGMSH_Create)(DataSet_t*) ;
extern void                (PosFilesForGMSH_Delete)(void*) ;
extern void                (PosFilesForGMSH_ParsedFileFormat)(PosFilesForGMSH_t*) ;
extern void                (PosFilesForGMSH_ASCIIFileFormat)(PosFilesForGMSH_t*) ;



#define PosFilesForGMSH_GetDataSet(PFG)                ((PFG)->dataset)
#define PosFilesForGMSH_GetOutputFiles(PFG)            ((PFG)->outputfiles)
#define PosFilesForGMSH_GetBilVersion(PFG)             ((PFG)->version)


#define PosFilesForGMSH_GetBilMajorVersion(PFG) \
        atoi(PosFilesForGMSH_GetBilVersion(PFG))
        
#define PosFilesForGMSH_GetBilMinorVersion(PFG) \
        atoi(String_FindChar(PosFilesForGMSH_GetBilVersion(PFG),'.')+1)


/* complete the structure types by using the typedef */
struct PosFilesForGMSH_t {
  char version[10] ;
  DataSet_t* dataset ;
  OutputFiles_t* outputfiles ;  /* The output files */
} ;

#include "String_.h"


#endif
