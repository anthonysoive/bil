#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
//#include <float.h>
//#include <ctype.h>
#include <string.h>
//#include <strings.h>
//#include <assert.h>
#include "Message.h"
#include "Expression.h"

/*
 * Evaluate expressions in a string
 * C code created by a parser generator from AnaGram
 * http://www.parsifalsoft.com/
 */

#ifndef strdup
#define strdup strdup0

static char* (strdup0)(const char* s)
{
  size_t slen = strlen(s);
  char* result = (char*) malloc(slen + 1);
  if(result == NULL)
  {
    return NULL;
  }

  memcpy(result, s, slen+1);
  return result;
}
#endif



#include "Libraries/Expression/ForTheRecord/evaluateExpression/evalwrap.c"
#include "Libraries/Expression/ForTheRecord/evaluateExpression/evalkern.c"




double (Expression_Evaluate)(char* variablename,char* expressionstrings)
{
  double val ;
  int errorFlag = evaluateExpression(expressionstrings) ;
  
  if(errorFlag) {
    Message_RuntimeError("Expression_Evaluate: syntax error\n\
                          %s at line %d, column %d\n",\
                        errorRecord.message,\
                        errorRecord.line,\
                        errorRecord.column) ;
  }
  
  {
    int i = 0 ;
    
    while(strncmp(variablename,variable[i].name,strlen(variablename))) i++ ;
    
    val = variable[i].value ;
  }
  
  return(val) ;
}
