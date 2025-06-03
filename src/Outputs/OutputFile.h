#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct OutputFile_t; //typedef struct OutputFile_t  OutputFile_t ;
struct TextFile_t;



extern char    OutputFile_TypeOfCurrentFile ;

extern OutputFile_t*  (OutputFile_Create)(char*) ;
extern void           (OutputFile_Delete)(void*) ;


#define OutputFile_GetTextFile(OF)            ((OF)->textfile)


#define OutputFile_IsPointType      (OutputFile_TypeOfCurrentFile == 'p')
#define OutputFile_IsTimeType       (OutputFile_TypeOfCurrentFile == 't')




struct OutputFile_t {             /* Output file */
  TextFile_t* textfile ;          /* The file */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
