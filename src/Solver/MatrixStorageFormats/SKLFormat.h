#ifndef SKLFORMAT_H
#define SKLFORMAT_H


/* Forward declarations */
struct SKLFormat_t; //typedef struct SKLFormat_t   SKLFormat_t ;



/** The getters */
#define SKLFormat_GetNbOfNonZeroValues(a)                    ((a)->nnz)
#define SKLFormat_GetNonZeroValue(a)                         ((a)->nzval)
#define SKLFormat_GetPointerToLowerRow(a)                    ((a)->l)
#define SKLFormat_GetPointerToUpperColumn(a)                 ((a)->u)
#define SKLFormat_GetFirstNonZeroValueIndexOfLowerRow(a)     ((a)->rowptr)
#define SKLFormat_GetFirstNonZeroValueIndexOfUpperColumn(a)  ((a)->colptr)



/* complete the structure types by using the typedef */

/* Skyline format */
struct SKLFormat_t {          /* Skyline storage format */
  size_t    nnz ;       /* Nb of non zero values */
  double* nzval ;             /* Non zero values */
  double* l ;                 /* Strictly lower triangular matrix values */
  double* u ;                 /* Upper triangular matrix values including diagonal */
  size_t* colptr ;      /* Index of element in u which starts a column */
  size_t* rowptr ;      /* Index of element in l which starts a row */
} ;

#endif
