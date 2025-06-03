#ifndef FIELDS_H
#define FIELDS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Fields_t; //typedef struct Fields_t       Fields_t ;
struct Field_t;
struct DataFile_t;


extern Fields_t*   (Fields_New)     (const int) ;
extern Fields_t*   (Fields_Create)  (DataFile_t*) ;
extern void        (Fields_Delete)  (void*) ;


#define Fields_GetNbOfFields(FLDS)    ((FLDS)->n_ch)
#define Fields_GetField(FLDS)         ((FLDS)->ch)




struct Fields_t {             /* fields */
  unsigned int n_ch ;         /* nb of fields */
  Field_t* ch ;               /* field */
} ;


#ifdef __CPLUSPLUS
}
#endif

#endif
