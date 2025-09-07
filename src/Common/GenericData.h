#ifndef GENERICDATA_H
#define GENERICDATA_H


/* Forward declarations */
struct GenericData_t; //typedef struct GenericData_t     GenericData_t ;
struct TypeId_t;






#define GenericData_MaxLengthOfKeyWord           (50)



#define GenericData_GetTypeId(GD)                ((GD)->GetTypeId())
#define GenericData_GetName(GD)                  ((GD)->GetName())
#define GenericData_GetNbOfData(GD)              ((GD)->GetNbOfData())
//#define GenericData_GetData(GD)                  ((GD)->GetData())

#define GenericData_GetNextGenericData(GD)       ((GD)->GetNextGenericData())
#define GenericData_GetPreviousGenericData(GD)   ((GD)->GetPreviousGenericData())
#define GenericData_GetDelete(GD)                ((GD)->GetDelete())

#define GenericData_SetTypeId(GD,X)              ((GD)->SetTypeId(X))
#define GenericData_SetName(GD,X)                ((GD)->SetName(X))
#define GenericData_SetNbOfData(GD,X)            ((GD)->SetNbOfData(X))
//#define GenericData_SetData(GD,X)                ((GD)->SetData(X))

#define GenericData_SetNextGenericData(GD,X)     ((GD)->SetNextGenericData(X))
#define GenericData_SetPreviousGenericData(GD,X) ((GD)->SetPreviousGenericData(X))
#define GenericData_SetDelete(GD,X)              ((GD)->SetDelete(X))


#define GenericData_GetData(GD) \
        TypeId_GetData(GenericData_GetTypeId(GD))
        
#define GenericData_Set(GD,X) \
        TypeId_Set(GenericData_GetTypeId(GD),X)


#define GenericData_GetSize(GD) \
        TypeId_GetSize(GenericData_GetTypeId(GD))

#define GenericData_GetIdNumber(GD) \
        TypeId_GetIdNumber(GenericData_GetTypeId(GD))


        
        
#define GenericData_Is(GD,NAME) \
        (!strncmp(GenericData_GetName(GD),NAME,GenericData_MaxLengthOfKeyWord))
        
        
#define GenericData_FindData(GD,NAME) \
        GenericData_GetData_(GenericData_Find(GD,NAME))
        
        
#define GenericData_FindNbOfData(GD,NAME) \
        GenericData_GetNbOfData_(GenericData_Find(GD,NAME))
        
        
#define GenericData_Merge(A,B) \
        GenericData_Append(A,GenericData_First(B))





/* Implementation */
#define GenericData_GetNbOfData_(GD) \
        ((GD) ? GenericData_GetNbOfData(GD) : 0)

#define GenericData_GetData_(GD) \
        ((GD) ? GenericData_GetData(GD) : NULL)



/* Typedef names of methods */
//typedef void* GenericData_Allocate(int,TypedId_t) ;



#include "GenericObject.h"
#include <stdio.h>

/* Generic data */
struct GenericData_t {
  private:
  TypeId_t* _typ;                /* The type id of data */
  char* _name;                   /* Name of the data */
  size_t _n;                     /* Nb of data */
  //void* _data; /* The data */
  GenericData_t* _prev;          /* Previous generic data */
  GenericData_t* _next;          /* Next generic data */
  GenericObject_Delete_t* _delete;
  
  public:
  /* The getters */
  TypeId_t* GetTypeId(void){return(_typ);}
  char* GetName(void){return(_name);}
  size_t GetNbOfData(void){return(_n);}
  //auto GetData(void){return(_data);}
  GenericData_t* GetNextGenericData(void){return(_next);}
  GenericData_t* GetPreviousGenericData(void){return(_prev);}
  GenericObject_Delete_t* GetDelete(void){return(_delete);}
  
  /* The setters */
  void SetTypeId(TypeId_t* x){_typ = x;}
  void SetName(char* x){_name = x;}
  void SetNbOfData(size_t x){_n = x;}
  //template<typename T> void SetData(T x){_data = x;}
  void SetNextGenericData(GenericData_t* x){_next = x;}
  void SetPreviousGenericData(GenericData_t* x){_prev = x;}
  void SetDelete(void (*x)(void*)){_delete = x;}
} ;





inline GenericData_t* (GenericData_New)       (void) ;
template<typename T>
inline GenericData_t* (GenericData_Create)    (size_t,T*,const char*) ;
template<typename T>
inline void           (GenericData_Initialize)(GenericData_t*,size_t,T*,const char*) ;
inline GenericData_t* (GenericData_Append)    (GenericData_t*,GenericData_t*) ;
inline GenericData_t* (GenericData_First)     (GenericData_t*) ;
inline GenericData_t* (GenericData_Last)      (GenericData_t*) ;
inline GenericData_t* (GenericData_Find)      (GenericData_t*,const char*) ;

//inline void           (GenericData_DeleteData)   (TypeId_t*,void*) ;
inline void           (GenericData_Delete)       (void*) ;


/* For the macros */
#include "TypeId.h"

#include "GenericData.h.in"

#endif
