#ifndef COMMANDLINE_H
#define COMMANDLINE_H


/* Forward declarations */
struct CommandLine_t; //typedef struct CommandLine_t  CommandLine_t ;


extern CommandLine_t*    (CommandLine_Create)(int,char**) ;
extern void              (CommandLine_Delete)(void*) ;


#define CommandLine_GetNbOfArg(cmd)            ((cmd)->argc)
#define CommandLine_GetArg(cmd)                ((cmd)->argv)



struct CommandLine_t {        /* Command line */
  int    argc ;               /* Nb of command line arguments */
  char** argv ;               /* Command line arguments */
} ;


#endif
