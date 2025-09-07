#ifndef REGIONS_H
#define REGIONS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Regions_t; //typedef struct Regions_t        Regions_t ;
struct Region_t;


extern Regions_t*  (Regions_New)            (void) ;
extern void        (Regions_Delete)         (void*) ;
extern int         (Regions_FindRegionIndex)(Regions_t*,const char*) ;
extern Region_t*   (Regions_FindRegion)     (Regions_t*,const char*) ;


        
#define Regions_MaxNbOfRegions  (20)


#define Regions_GetNbOfRegions(RS)          ((RS)->NbOfRegions)
#define Regions_GetRegion(RS)               ((RS)->Region)



struct Regions_t {
  Region_t*    Region ;
  int NbOfRegions ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
