#ifndef STRING_H
#define STRING_H

#ifdef __CPLUSPLUS
extern "C" {
#endif



/* Forward declarations */
//struct String_t; //typedef struct String_t       String_t ;



static int   String_bytes;
static char* String_pchar;
static char const* String_pcchar;
static void* String_void;


//extern char*       (String_Create)         (const char*) ;
//extern void        (String_Delete)         (void*) ;
extern char* (String_FindToken)        (char*,const char*) ;
extern char const* (String_FindToken)  (char const*,const char*) ;
extern char* (String_FindToken)        (char*,const char*,const char*) ;
extern char const* (String_FindToken)  (char const*,const char*,const char*) ;
extern char* (String_FindAndSkipToken) (char*,const char*) ;
extern char const* (String_FindAndSkipToken) (char const*,const char*) ;
extern char* (String_FindAndSkipToken) (char*,const char*,const char*) ;
extern char const* (String_FindAndSkipToken) (char const*,const char*,const char*) ;
extern char* (String_FindNthToken)     (char*,const char*,const int) ;
extern char* (String_FindNthToken)     (char*,const char*,const char*,const int) ;
extern int   (String_CountTokens)         (char*,const char*) ;
extern int   (String_CountTokens)         (char*,const char*,const char*) ;
extern int   (String_CountTokensAloneInOneLine)(char*,const char*) ;
extern int   (String_CountTokensAloneInOneLine)(char*,const char*,const char*) ;
extern char* (String_CopyLine)             (const char*) ;
//extern const char* (String_SkipRemainingComments)(const char*) ;
//extern int         (String_NbOfUncommentedLines) (const char*,const char*) ;
extern int   (String_FindPositionIndex)    (const char*,const char* const*,const int) ;
extern char* (String_RemoveComments)       (char const*,char*) ;



#define String_MaxLengthOfKeyWord     (30)

#define String_MaxNbOfKeyWords        (10)
#define String_MaxLengthOfKeyWords    (String_MaxNbOfKeyWords*String_MaxLengthOfKeyWord)
#define String_MaxLengthOfLine        (500)






/* Scan string
 * ----------- */
/** Scan a string with a given format and return the nb of characters read. */
#define String_Scan(STR,...) \
        Logic_IF(Logic_GE(Arg_NARG(__VA_ARGS__),2))\
        (String_ScanN,String_Scan2)(STR,__VA_ARGS__)


#define String_ScanStringUntil(STR,KEY,END) \
        String_Scan(STR,"%*[ ]%[^" END "]",KEY)


#define String_ScanAffectedKeyphrase(STR,KEY) \
        String_Scan(STR,"%*[ ]%[^=]",KEY)


#define String_FindAndScanExp(STR,TOK,DEL, ...) \
        ((String_pchar  = String_FindAndSkipToken(STR,TOK,DEL)) ? \
        ((String_pchar  = String_FindToken(String_pchar," "))   ? \
        ((String_Scan(String_pchar,__VA_ARGS__))) : 0) : 0)


/** Reads N data from the string STR with the format "FMT" and store into V . */
#define String_ScanArray(STR,N,FMT,V) \
        do { \
          char* String_c = STR ; \
          for(size_t String_i = 0 ; String_i < (N) ; String_i++) { \
            String_c += String_Scan(String_c,FMT,(V) + String_i) ; \
          } \
        } while(0)
        
/** Reads N data for each entries from the string STR with the format "FMT". */
#define String_ScanArrays(STR,N,FMT,...) \
        do { \
          char* String_c = STR ; \
          for(size_t String_i = 0 ; String_i < (N) ; String_i++) { \
            String_c += String_Scan(String_c,FMT,String_IncrementAll(__VA_ARGS__)) ; \
          } \
        } while(0)
        
#define String_ScanArrays0(STR,N,FMT,...) \
        do { \
          char* String_c = STR ; \
          for(size_t String_i = 0 ; String_i < (N) ; String_i++) { \
            String_c += String_Scan(String_c,FMT,__VA_ARGS__) ; \
            Algos_SEPWITH(String_IncrementAll0(__VA_ARGS__),;) ; \
          } \
        } while(0)


/** Gets the advanced position in the string. */
#define String_GetAdvancedPosition \
        (String_pchar)


/* Implementation */
#define String_Scan2(STR,FMT) \
        (sscanf(STR,String_Fmt(FMT),&String_bytes) , \
         String_pchar = STR + String_bytes , \
         String_bytes)
        
#define String_ScanN(STR,FMT, ...) \
        (sscanf(STR,String_Fmt(FMT),__VA_ARGS__,&String_bytes) , \
         String_pchar = STR + String_bytes , \
         String_bytes)

#define String_Fmt(FMT)  FMT"%n"

#define String_Increment(a) \
        ((a) + String_i)
        
#define String_IncrementAll(...) \
        Tuple_SEQ(Algos_MAP(Tuple_TUPLE(__VA_ARGS__),String_Increment))

#define String_Increment0(a) \
        (a++)
        
#define String_IncrementAll0(...) \
        Algos_MAP(Tuple_TUPLE(__VA_ARGS__),String_Increment0)

//        std::strcat(strcpy(String_Save,FMT),"%n")




/* Find characters
 * --------------- */
#define String_FindChar(STR,C) \
        ((STR) ? std::strchr(STR,C) : NULL)
        

#define String_FindAnyChar(STR,Cs) \
        ((STR) ? std::strpbrk(STR,Cs) : NULL)
        

#define String_FindEndOfLine(STR) \
        String_FindChar(STR,'\n')
        

#define String_FindEndOfString(STR) \
        String_FindChar(STR,'\0')



/* Skip tokens/characters
 * ---------------------- */
#define String_SpaceChars \
        " \f\n\r\t\v"

#define String_BlankChars \
        " \t\r"


#define String_SkipAnyChars(STR,Cs) \
        ((STR) ? (STR) + std::strspn(STR,Cs) : NULL)
        

#define String_SkipAnyOtherChars(STR,Cs) \
        ((STR) ? (STR) + std::strcspn(STR,Cs) : NULL)
        

#define String_SkipBlankChars(STR) \
        String_SkipAnyChars(STR,String_BlankChars)
        

#define String_SkipNonBlankChars(STR) \
        String_SkipAnyOtherChars(STR,String_BlankChars)
        

#define String_SkipSpaceChars(STR) \
        String_SkipAnyChars(STR,String_SpaceChars)
        

#define String_SkipNonSpaceChars(STR) \
        String_SkipAnyOtherChars(STR,String_SpaceChars)
        

#define String_SkipLine(STR) \
        ((String_pchar = String_FindEndOfLine(STR)),String_SkipSpaceChars(String_pchar))


#define String_SkipNextToken(STR) \
        String_SkipNonBlankChars(String_SkipBlankChars(STR))
        



/* Compare with characters
 * ----------------------- */
#define String_Is(...) \
        Utils_CAT_NARG(String_Is,__VA_ARGS__)(__VA_ARGS__)
        
#define String_IsNot(...) \
        !String_Is(__VA_ARGS__)

#define String_CaseIgnoredIs(...) \
        Utils_CAT_NARG(String_CaseIgnoredIs,__VA_ARGS__)(__VA_ARGS__)
        
#define String_BeginsWithAnyChar(STR,Cs) \
        ((STR) ? std::strspn(STR,Cs) : 0)

#define String_BeginsWithSingleLineComment(STR) \
        (String_Is(STR,"#",1) || String_Is(STR,"//",2))

#define String_BeginsWithMultiLineComment(STR) \
        String_Is(STR,"/*",2)

#define String_SkipMultiLineComment(STR) \
        String_FindAndSkipToken(STR,"*/")


/* Implementation */
#define String_Is2(STR,...) \
        ((STR) ? (!std::strcmp(STR,__VA_ARGS__)) : 0)
        
#define String_Is3(STR,...) \
        ((STR) ? (!std::strncmp(STR,__VA_ARGS__)) : 0)
        
#define String_CaseIgnoredIs2(STR,...) \
        ((STR) ? (!strcasecmp(STR,__VA_ARGS__)) : 0)
        
#define String_CaseIgnoredIs3(STR,...) \
        ((STR) ? (!strncasecmp(STR,__VA_ARGS__)) : 0)



#if 0
#define String_GetStringLength(STR)                 ((STR)->length)
#define String_GetStringContent(STR)                ((STR)->head)
#define String_GetCurrentPositionInString(STR)      ((STR)->current)




struct String_t {
  char*     head ;         /* String content */
  int       length ;        /* String length */
  char*     current ;       /* Current position in the string */
} ;
#endif



#ifdef __CPLUSPLUS
}
#endif

#include <cstring>
//#include <strings.h>
#include <stdarg.h>
#include <stdio.h>
#include "Arg.h"
#include "Tuple.h"
#include "Algos.h"
#include "Logic.h"
#include "Utils.h"
#endif
