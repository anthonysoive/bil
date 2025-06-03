#ifndef FVM_H
#define FVM_H



/* Forward eclarations */
struct FVM_t;  //typedef struct FVM_t     FVM_t ;
struct Buffers_t;
struct Element_t;

/*  Typedef names of Methods */
//typedef void     FVM_ComputeFluxes_t(FVM_t*,double*,double*,int,int) ;


#define FVM_MaxNbOfNodes                         (Element_MaxNbOfNodes)
#define FVM_MaxNbOfDOF                           (Element_MaxNbOfDOF)
#define FVM_MaxShift                             (Model_MaxNbOfEquations*Model_MaxNbOfEquations)
#define FVM_MaxNbOfMatrices                      (4)

#define FVM_MaxSizeOfMatrix                      (FVM_MaxNbOfDOF*FVM_MaxNbOfDOF*sizeof(double))
#define FVM_MaxSizeOfOutput                      (FVM_MaxNbOfDOF*FVM_MaxNbOfDOF*sizeof(double))
#define FVM_MaxSizeOfInput                       (FVM_MaxNbOfNodes*FVM_MaxNbOfNodes*FVM_MaxShift*sizeof(double))
#define FVM_SizeOfBuffer                         (FVM_MaxNbOfMatrices*FVM_MaxSizeOfMatrix)





#define FVM_GetElement(fvm)                      ((fvm)->GetElement())
#define FVM_GetInput(fvm)                        ((fvm)->GetInput())
#define FVM_GetOutput(fvm)                       ((fvm)->GetOutput())
#define FVM_GetShiftOfInput(fvm)                 ((fvm)->GetShiftOfInput())
#define FVM_GetBuffers(fvm)                      ((fvm)->GetBuffers())
#define FVM_GetCellVolumes(fvm)                  ((fvm)->GetCellVolumes())
#define FVM_GetCellSurfaceAreas(fvm)             ((fvm)->GetCellSurfaceAreas())
#define FVM_GetIntercellDistances(fvm)           ((fvm)->GetIntercellDistances())


#define FVM_SetElement(fvm,A)                    ((fvm)->SetElement(A))
#define FVM_SetInput(fvm,A)                      ((fvm)->SetInput(A))
#define FVM_SetOutput(fvm,A)                     ((fvm)->SetOutput(A))
#define FVM_SetShiftOfInput(fvm,A)               ((fvm)->SetShiftOfInput(A))
#define FVM_SetBuffers(fvm,A)                    ((fvm)->SetBuffers(A))
#define FVM_SetCellVolumes(fvm,A)                ((fvm)->SetCellVolumes(A))
#define FVM_SetCellSurfaceAreas(fvm,A)           ((fvm)->SetCellSurfaceAreas(A))
#define FVM_SetIntercellDistances(fvm,A)         ((fvm)->SetIntercellDistances(A))




/* Buffer */
#define FVM_GetBuffer(fvm) \
        Buffers_GetBufferOfCurrentThread(FVM_GetBuffers(fvm))



#define FVM_AllocateInBuffer(fvm,sz) \
        (Buffer_Allocate(FVM_GetBuffer(fvm),(sz)))
        
#define FVM_FreeBuffer(fvm) \
        (Buffer_Free(FVM_GetBuffer(fvm)))
        
#define FVM_FreeBufferFrom(fvm,p) \
        (Buffer_FreeFrom(FVM_GetBuffer(fvm),(char*) (p)))



#include "GenericObject.h"

struct FVM_t {
  public:
  
  void SetElement(Element_t* el) {_el = el;}
  void SetInput(void* input) {_input = input;}
  void SetOutput(void* output) {_output = output;}
  void SetShiftOfInput(int shift) {_shift = shift;}
  void SetBuffers(Buffers_t* buffers) {_buffers = buffers;}
  void SetCellVolumes(double* cellvolumes) {_cellvolumes = cellvolumes;}
  void SetCellSurfaceAreas(double* cellsurfaceareas) {_cellsurfaceareas = cellsurfaceareas;}
  void SetIntercellDistances(double* celldistances) {_celldistances = celldistances;}
  
  /* Accessors */
  Element_t* GetElement() {return(_el);}
  void*      GetInput() {return(_input);}
  void*      GetOutput() {return(_output);}
  int        GetShiftOfInput() {return(_shift);}
  Buffers_t* GetBuffers() {return(_buffers);}
  double*    GetCellVolumes() {return(_cellvolumes);}
  double*    GetCellSurfaceAreas() {return(_cellsurfaceareas);}
  double*    GetIntercellDistances() {return(_celldistances);}
  
  //GenericObject_Delete_t* Delete ;
  
  private:
  Element_t* _el ;             /* Element */
  void*      _input ;          /* Input */
  void*      _output ;         /* Output*/
  int        _shift ;          /* Shift of input */
  Buffers_t* _buffers ;        /* Buffer */
  double*    _cellvolumes ;
  double*    _cellsurfaceareas ;
  double*    _celldistances ;
} ;



struct Load_t;
struct Mesh_t;
struct IntFct_t;

extern FVM_t*     (FVM_GetInstance)(Element_t*) ;
extern void       (FVM_Delete)(void*) ;
extern double*    (FVM_ComputeMassMatrix)(FVM_t*,double*,int) ;
extern double*    (FVM_ComputeIsotropicConductionMatrix)(FVM_t*,double*,int) ;
extern double*    (FVM_ComputeMassAndIsotropicConductionMatrix)(FVM_t*,double*,int) ;
extern double*    (FVM_ComputeSurfaceLoadResidu)(FVM_t*,Load_t*,double,double) ;
extern double*    (FVM_ComputeBodyForceResidu)(FVM_t*,double const*,int const) ;
extern double*    (FVM_ComputeFluxResidu)(FVM_t*,double const*,int const) ;
extern double*    (FVM_ComputeMassAndFluxResidu)(FVM_t*,double const*,int const) ;
extern double*    (FVM_ComputeMassBalanceEquationResidu)(FVM_t*,double const*,double const*,double const) ;
extern double*    (FVM_ComputeCellVolumes)(FVM_t*) ;
extern double*    (FVM_ComputeCellSurfaceAreas)(FVM_t*) ;
extern double*    (FVM_ComputeCellVolumesAndSurfaceAreas)(FVM_t*) ;
extern short int  (FVM_FindLocalCellIndex)(FVM_t*,double*) ;
extern double*    (FVM_ComputeIntercellDistances)(FVM_t*) ;
extern double*    (FVM_ComputeTheNodalFluxVector)(FVM_t*,double*) ;
extern double     (FVM_AverageCurrentImplicitTerm)(Mesh_t*,const char*,const int,const int) ;
extern double     (FVM_AveragePreviousImplicitTerm)(Mesh_t*,const char*,const int,const int) ;
extern double*    (FVM_ComputeGradient)(FVM_t*,double*,IntFct_t*,int,int);


/* For the macros */
#include "Element.h"
#include "Model.h"
#include "Buffers.h"
#include "Buffer.h"


//#include "FVM.h.in"

#endif
