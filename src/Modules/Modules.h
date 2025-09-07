#ifndef MODULES_H
#define MODULES_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Vacuous declarations and typedef names */

/* class-like structure "Modules_t" */
struct Modules_t       ; typedef struct Modules_t       Modules_t ;
struct Module_t;


extern void       (Modules_Delete)(void*) ;
extern void       (Modules_Print)(char*) ;
extern Module_t*  (Modules_FindModule)(Modules_t*,const char*) ;



#define Modules_NbOfModules               (ListOfModules_Nb)
#define Modules_ListOfNames               ListOfModules_Names
#define Modules_ListOfSetModuleProp       ListOfModules_Methods(_SetModuleProp)


#define Modules_GetNbOfModules(MODS)    ((MODS)->n_modules)
#define Modules_GetModule(MODS)         ((MODS)->module)


#define Modules_SetNbOfModules(MODS,A) \
        do {\
          Modules_GetNbOfModules(MODS) = A;\
        } while(0)
        
#define Modules_SetModule(MODS,A) \
        do {\
          Modules_GetModule(MODS) = A;\
        } while(0)


struct Modules_t {              /* modules */
  int n_modules ;      /* nb of modules */
  Module_t* module ;            /* module */
} ;


/* For the macros */
#include "ListOfModules.h"

#ifdef __CPLUSPLUS
}
#endif
#endif
