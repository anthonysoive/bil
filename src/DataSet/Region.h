#ifndef REGION_H
#define REGION_H

#ifdef __CPLUSPLUS
extern "C" {
#endif



/* Forward declarations */
struct Region_t; //typedef struct Region_t         Region_t ;


extern Region_t*  (Region_New)   (void) ;
extern void       (Region_Delete)(void*) ;


        
#define Region_MaxLengthOfRegionName  (50)


/* Accessors */
#define Region_GetRegionName(R)           ((R)->RegionName)


struct Region_t {
  char* RegionName ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
