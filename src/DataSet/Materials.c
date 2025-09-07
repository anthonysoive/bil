#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <assert.h>
#include "Message.h"
#include "DataFile.h"
#include "Materials.h"
#include "Material.h"
#include "Models.h"
#include "Model.h"
#include "Curves.h"
#include "Mry.h"
#include "Geometry.h"
#include "Fields.h"
#include "Functions.h"


/* Extern functions */

Materials_t* (Materials_New)(const int n_mats,Models_t* models)
{
  Materials_t* materials   = (Materials_t*) Mry_New(Materials_t) ;


  Materials_GetNbOfMaterials(materials) = n_mats ;

  /* Allocate the materials */
  {
    Material_t* material   = (Material_t*) Mry_New(Material_t,n_mats) ;
    int    i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat   = Material_New() ;
      
      material[i] = mat[0] ;
      free(mat) ;
    }
    
    Materials_GetMaterial(materials) = material ;
  }
  
  
  /* Allocate the space for the models used by the materials */
  {
    /* We create the space for n_mats models max */
    int n_models = n_mats ;
    Models_t* usedmodels = (models) ? models : Models_New(n_models) ;
    
    Materials_GetUsedModels(materials) = usedmodels ;
  }
  
  
  /* All materials share the same pointer to usedmodels */
  {
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    int    i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat   = Materials_GetMaterial(materials) + i ;
      
      Material_GetUsedModels(mat) = usedmodels ;
    }
  }
  
  
  return(materials) ;
}



Materials_t* (Materials_Create)(DataFile_t* datafile,Geometry_t* geom,Fields_t* fields,Functions_t* functions,Models_t* models)
{
  int n_mats = DataFile_CountTokens(datafile,"MATE,Material",",") ;
  Materials_t* materials = Materials_New(n_mats,models) ;
  
  
  Message_Direct("Enter in %s","Materials") ;
  Message_Direct("\n") ;
  

  /* Scan the datafile */
  {
    int i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
      char* c = DataFile_FindNthToken(datafile,"MATE,Material",",",i + 1) ;
      
      c = String_SkipLine(c) ;
      
      DataFile_SetCurrentPositionInFileContent(datafile,c) ;
  
      Message_Direct("Enter in %s %d","Material",i+1) ;
      Message_Direct("\n") ;
      
      Material_GetFields(mat) = fields ;
      Material_GetFunctions(mat) = functions ;
      
      Material_Scan(mat,datafile,geom) ;
    }
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(materials) ;
}



void (Materials_Delete)(void* self)
{
  Materials_t* materials = (Materials_t*) self ;
  
  {
    int n_mats = Materials_GetNbOfMaterials(materials) ;
    Material_t* material = Materials_GetMaterial(materials) ;
    
    if(material) {
      int i ;
      
      for(i = 0 ; i < n_mats ; i++) {
        Material_Delete(material+i) ;
      }
      //Mry_Delete(material,n_mats,Material_Delete) ;
      free(material) ;
      Materials_GetMaterial(materials) = NULL ;
    }
  }
  
  {
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    
    if(usedmodels) {
      Models_Delete(usedmodels) ;
      free(usedmodels) ;
      Materials_GetUsedModels(materials) = NULL ;
    }
  }
}
