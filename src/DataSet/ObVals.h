#ifndef OBVALS_H
#define OBVALS_H


/* Forward declarations */
struct ObVals_t; //typedef struct ObVals_t       ObVals_t ;
struct DataFile_t;
struct Mesh_t;
struct Materials_t;
struct ObVal_t;


extern ObVals_t*  (ObVals_New)(const int) ;
extern ObVals_t*  (ObVals_Create)(DataFile_t*,Mesh_t*,Materials_t*) ;
extern void       (ObVals_Delete)(void*) ;
extern int        (ObVals_FindObValIndex)(ObVals_t*,char*) ;



#define ObVals_GetNbOfObVals(OVS)    ((OVS)->n_obj)
#define ObVals_GetObVal(OVS)         ((OVS)->obj)



struct ObVals_t {             /* objective variations */
  int n_obj ;        /* nb */
  ObVal_t* obj ;              /* objective variation */
} ;


#endif
