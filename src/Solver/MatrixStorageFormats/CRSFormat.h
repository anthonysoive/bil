#ifndef CRSFORMAT_H
#define CRSFORMAT_H



/* Forward declarations */
struct CRSFormat_t; //typedef struct CRSFormat_t      CRSFormat_t ;


/** The getters */
#define CRSFormat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define CRSFormat_GetNonZeroValue(a)                    ((a)->nzval)
#define CRSFormat_GetColumnIndexOfNonZeroValue(a)       ((a)->colind)
#define CRSFormat_GetFirstNonZeroValueIndexOfRow(a)     ((a)->rowptr)
                                                      


/* complete the structure types by using the typedef */

/* Compressed row storage format
 * If a_ij = nzval[k] then colind[k] = j and rowptr[i] <= k < rowptr[i + 1] */
struct CRSFormat_t {
  size_t    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  size_t* colind ;      /* Column indices of the non zeros */
  size_t* rowptr ;      /* Index of element in nzval which starts a row */
} ;

#endif
