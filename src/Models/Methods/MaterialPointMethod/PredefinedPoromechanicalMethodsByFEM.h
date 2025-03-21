#ifndef PREDEFINEDPOROMECHANICALMETHODSBYFEM_H
#define PREDEFINEDPOROMECHANICALMETHODSBYFEM_H

#include "PredefinedModelMethods.h"
#include "BaseName.h"
#include "CustomValues.h"
#include "MaterialPointMethod.h"


#define ImplicitValues_t BaseName(_ImplicitValues_t)
#define ExplicitValues_t BaseName(_ExplicitValues_t)
#define ConstantValues_t BaseName(_ConstantValues_t)


template<typename T>
struct ImplicitValues_t ;

template<typename T>
struct ExplicitValues_t;

template<typename T>
struct ConstantValues_t;


#define Values_t    BaseName(_Values_t)

template<typename T>
using Values_t = CustomValues_t<T,ImplicitValues_t,ExplicitValues_t,ConstantValues_t> ;


#define MPM_t      BaseName(_MPM_t)

struct MPM_t: public MaterialPointMethod_t<Values_t> {
  MaterialPointMethod_SetInputs_t<Values_t> SetInputs;
  template<typename T>
  MaterialPointMethod_Integrate_t<Values_t,T> Integrate;
  Values_t<double>* Integrate(Element_t* el,double const& t,double const& dt,Values_t<double> const& val_n,Values_t<double>& val) {return(Integrate<double>(el,t,dt,val_n,val));}
  #ifdef USE_AUTODIFF
  Values_t<real>* Integrate(Element_t* el,double const& t,double const& dt,Values_t<double> const& val_n,Values_t<real>& val) {return(Integrate<real>(el,t,dt,val_n,val));}
  #endif
  MaterialPointMethod_Initialize_t<Values_t>  Initialize;
  MaterialPointMethod_SetTangentMatrix_t<Values_t> SetTangentMatrix;
  MaterialPointMethod_SetTransferMatrix_t<Values_t> SetTransferMatrix;
  MaterialPointMethod_SetIndexes_t SetIndexes;
  MaterialPointMethod_SetIncrements_t SetIncrements;
} ;



#define Values_d    BaseName(_Values_d)
using Values_d = Values_t<double> ;
#define Values_Index(V)  CustomValues_Index(Values_d,V,double)

#endif
