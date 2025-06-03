#ifndef POINT_H
#define POINT_H


/* Forward declarations */
struct Point_t; //typedef struct Point_t        Point_t ;
struct Mesh_t;
struct Element_t;


extern Point_t*  (Point_New)                  (void) ;
extern void      (Point_Delete)               (void*) ;
extern void      (Point_SetEnclosingElement)  (Point_t*,Mesh_t*) ;
extern void      (Point_Scan)                 (Point_t*,char*) ;



#define Point_MaxLengthOfRegionName      Region_MaxLengthOfRegionName


#define Point_GetCoordinate(PT)                    ((PT)->Coordinate)
#define Point_GetEnclosingElement(PT)              ((PT)->EnclosingElement)
#define Point_GetRegionTag(PT)                     ((PT)->RegionTag)
#define Point_GetRegionName(PT)                    ((PT)->RegionName)


struct Point_t {
  double* Coordinate ;
  int   RegionTag ;
  char* RegionName ;
  Element_t* EnclosingElement ;  /* Element inside which the point lies */
} ;


#include "Region.h"
#endif
