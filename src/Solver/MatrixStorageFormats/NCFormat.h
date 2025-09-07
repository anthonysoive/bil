#ifndef NCFORMAT_H
#define NCFORMAT_H



/* Forward declarations */
struct NCFormat_t; //typedef struct NCFormat_t      NCFormat_t ;
//#define NCFormat_t        NCformat     /* Mimic the NCformat struct defined in SuperLu */
struct Mesh_t;


//extern NCFormat_t* (NCFormat_Create)(Mesh_t*) ;
extern NCFormat_t* (NCFormat_Create)(Mesh_t*,const int) ;
extern void        (NCFormat_Delete)(void*) ;
extern size_t      (NCFormat_AssembleElementMatrix)(NCFormat_t*,double*,int*,int*,int,int*,int) ;
extern void        (NCFormat_PrintMatrix)(NCFormat_t*,size_t,const char*) ;




/** The getters */
#define NCFormat_GetNbOfNonZeroValues(F)               ((F)->nnz)
#define NCFormat_GetNonZeroValue(F)                    ((F)->nzval)
#define NCFormat_GetRowIndexOfNonZeroValue(F)          ((F)->rowind)
#define NCFormat_GetFirstNonZeroValueIndexOfColumn(F)  ((F)->colptr)

/*
#define NCFormat_NbOfNonZeroValuesInColumn(F,j)        ((F)->colptr[j + 1] - (F)->colptr[j])
#define NCFormat_RowIndexStartingColumn(F,j)           ((F)->rowind[(F)->colptr[j]])
#define NCFormat_UpperColumn(F,j)                      ((F)->nzval + (F)->colptr[j + 1] - (j))
*/


/* complete the structure types by using the typedef */

/* Same as CCSFormat
 * Compressed column storage format (known as Harwell-Boeing sparse matrix format)
 * If a_ij = nzval[k] then rowind[k] = i and colptr[j] <= k < colptr[j + 1] */
struct NCFormat_t {
  size_t    nnz ;             /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  size_t* rowind ;            /* Row indices of the non zeros */
  size_t* colptr ;            /* Index of element in nzval which starts a column */
} ;


#endif
