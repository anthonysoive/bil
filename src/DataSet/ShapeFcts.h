#ifndef SHAPEFCTS_H
#define SHAPEFCTS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif



/* Forward declarations */
struct ShapeFcts_t; //typedef struct ShapeFcts_t      ShapeFcts_t ;
struct ShapeFct_t;


extern ShapeFcts_t*  (ShapeFcts_Create)(void) ;
extern void          (ShapeFcts_Delete)(void*) ;
extern int           (ShapeFcts_FindShapeFct)(ShapeFcts_t*,int,int) ;
extern int           (ShapeFcts_AddShapeFct)(ShapeFcts_t*,int,int) ;


#define ShapeFcts_MaxNbOfShapeFcts             (10)

#define ShapeFcts_GetNbOfShapeFcts(SFS)    ((SFS)->n_sh)
#define ShapeFcts_GetShapeFct(SFS)         ((SFS)->sh)


struct ShapeFcts_t {          /* Shape functions */
  unsigned int n_sh ;         /* Number of shape functions */
  ShapeFct_t*  sh ;           /* Shape function */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
