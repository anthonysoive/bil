#ifndef DATES_H
#define DATES_H


/* Forward declarations */
struct Dates_t; //typedef struct Dates_t        Dates_t ;
struct DataFile_t;
struct Date_t;


extern Dates_t*  (Dates_New)    (const int) ;
extern Dates_t*  (Dates_Create) (DataFile_t*) ;
extern void      (Dates_Delete) (void*) ;



#define Dates_GetNbOfDates(DATES)       ((DATES)->n_dates)
#define Dates_GetDate(DATES)            ((DATES)->date)



struct Dates_t {
  int n_dates ;               /* nb of dates */
  Date_t* date ;              /* date */
} ;

#endif
