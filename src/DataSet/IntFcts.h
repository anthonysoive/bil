#ifndef INTFCTS_H
#define INTFCTS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif



/* Forward declarations */
struct IntFcts_t; //typedef struct IntFcts_t      IntFcts_t ;
struct IntFct_t;

/* Declaration of Macros, Methods and Structures */


/* 1. IntFcts_t 
 * ------------*/
extern IntFcts_t*  (IntFcts_Create)(void) ;
extern void        (IntFcts_Delete)(void*) ;
extern int         (IntFcts_FindIntFct)(IntFcts_t*,unsigned short int,unsigned short int,const char*) ;
extern int         (IntFcts_AddIntFct)(IntFcts_t*,unsigned short int,unsigned short int,const char*) ;


#define IntFcts_MaxNbOfIntFcts             (4)

#define IntFcts_GetNbOfIntFcts(IFCTS)    ((IFCTS)->n_fi)
#define IntFcts_GetIntFct(IFCTS)         ((IFCTS)->fi)


struct IntFcts_t {            /* interpolations */
  int n_fi ;         /* nb of interpolation function */
  IntFct_t* fi ;              /* interpolation function */
} ;



#ifdef __CPLUSPLUS
}
#endif
#endif
