#ifndef MATH_H
#define MATH_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Forward declarations */
struct Math_t; //typedef struct Math_t   Math_t;
struct Buffers_t;


/* Constants */
#define Math_Pi           (3.14159265358979323846)
#define Math_Ln10         (2.3025850929940456840179914546843642076011)



#define Math_SolveByGaussEliminationWithoutPartialPivoting(...) \
        Math_SolveByGaussEliminationKIJ(__VA_ARGS__,NULL)

#define Math_SolveByGaussElimination(...) \
        Utils_CAT_NARG(Math_SolveByGaussElimination,__VA_ARGS__)(__VA_ARGS__)

/* Implementation */
#define Math_SolveByGaussElimination3 \
        Math_SolveByGaussEliminationWithPartialPivoting
        
#define Math_SolveByGaussElimination4 \
        Math_SolveByGaussEliminationKIJ


#include "BilConfig.h"

#ifdef HAVE_LAPACK
#if defined(__cplusplus)
  extern "C" {
#endif

extern int dgeev_(const char*,const char*,int*,double*,int*,double*, double*,double*,int*,double*,int*,double*,int*,int*);

#if defined(__cplusplus)
  }
#endif
#endif


#define Math_ComputeSecondDeviatoricStrainInvariant \
        Math_ComputeSecondDeviatoricStressInvariant
        
#define Math_ComputeDeviatoricStrain \
        Math_ComputeDeviatoricStress
        
#define Math_ComputeMatrixInverse(a,n) \
        Math_SolveByGaussJordanElimination(a,NULL,n,0)
        
#define Math_InvertMatrix(A,N)  \
        Math_SolveByGaussJordanElimination(A,NULL,N,0)

#define Math_Max(a,b)          (((a) > (b)) ? (a) : (b))
#define Math_Min(a,b)          (((a) < (b)) ? (a) : (b))

#define Math_Swap(a,b,Type_t) \
        do { \
          Type_t tmp = (a) ; \
          (a) = (b) ; \
          (b) = tmp ; \
        } while(0)
        
#define Math_SwapDouble(a,b)   Math_Swap(a,b,double)
#define Math_SwapInt(a,b)      Math_Swap(a,b,int)

#define Math_Sign(a)           ((a < 0) ? -1 : (a > 0))


#ifndef M_PI
#define M_PI       Math_Pi
#endif

#ifndef MAX
#define MAX        Math_Max
#endif

#ifndef MIN
#define MIN        Math_Min
#endif



#define Math_MaxNbOfMatrices                      (4)
#define Math_MaxSizeOfMatrix                      (9*sizeof(double))
#define Math_SizeOfBuffer                         (Math_MaxNbOfMatrices*Math_MaxSizeOfMatrix)


#define Math_GetOutputs(M)             ((M)->outputs)
#define Math_GetBuffers(M)             ((M)->buffers)
#define Math_GetDelete(M)              ((M)->Delete)



#define Math_GetBuffer(M) \
        Buffers_GetBufferOfCurrentThread(Math_GetBuffers(M))

/* Operations on buffer */
#define Math_AllocateInBuffer(M,sz) \
        (Buffer_Allocate(Math_GetBuffer(M),(sz)))
        
#define Math_FreeBuffer(M) \
        (Buffer_Free(Math_GetBuffer(M)))
        
#define Math_FreeBufferFrom(M,p) \
        (Buffer_FreeFrom(Math_GetBuffer(M),(char*) (p)))



#include "GenericObject.h"

struct Math_t {          /* some math methods */
  void* outputs ;
  Buffers_t*  buffers ;         /* Buffer */
  GenericObject_Delete_t* Delete ;
} ;



/* Old notations that I try to eliminate little by little */
#define j2       Math_ComputeSecondDeviatoricStressInvariant
//#define gausse   Math_SolveByGaussElimination


inline Math_t*  (Math_GetInstance)(void) ;
inline Math_t*  (Math_Create)(void) ;
inline int      (Math_ComputePolynomialEquationRootsOfDegree4)(double*) ;
inline int      (Math_ComputePolynomialEquationRootsOfDegree3)(double*) ;
inline int      (Math_ComputePolynomialEquationRootsOfDegree2)(double*) ;
inline void     (Math_Delete)(void*) ;
inline double   (Math_ComputeFirstStressInvariant)(const double*) ;
inline double   (Math_ComputeSecondStressInvariant)(const double*) ;
inline double   (Math_ComputeSecondDeviatoricStressInvariant)(double const*);
inline double*  (Math_SolveByGaussEliminationWithPartialPivoting)(double*,double*,int) ;
inline double*  (Math_SolveByGaussJordanElimination)(double*,double*,int,int) ;
inline int      (Math_ComputePolynomialEquationRoots)(double*,int const) ;
inline int      (Math_PolishPolynomialEquationRoot)(double*,int,double*,double,int) ;
inline double   (Math_Compute3x3MatrixDeterminant)(const double*) ;
inline double*  (Math_Inverse3x3Matrix)(const double*) ;
inline double*  (Math_ComputePrincipalStresses)(const double*) ;
inline double*  (Math_ComputeRealEigenvaluesAndEigenvectorsOf3x3Matrix)(double*,const char) ;
inline void     (Math_PrintStiffnessTensor)(const double*) ;
inline void     (Math_PrintStressTensor)(const double*) ;
inline void     (Math_PrintMatrix)(const double*,const int) ;
inline void     (Math_PrintVector)(const double*,const int) ;
inline double*  (Math_ComputeDeviatoricStress)(const double*) ;
inline double*  (Math_SolveByGaussEliminationKIJ)(double*,double*,int,int*) ;
inline double*  (Math_SolveByGaussEliminationJIK)(double*,double*,int,int*) ;


/* For the macros */
#include "Utils.h"
#include "Buffers.h"
#include "Buffer.h"

#include "Math_.h.in"




template<typename T>
inline T        (Math_ComputeSecondDeviatoricStressInvariant)(T const*);
template<typename T>
inline int      (Math_ComputePolynomialEquationRoots)(T*,int const) ;

template<typename T>
inline T (Math_ComputeSecondDeviatoricStressInvariant)(T const* sig)
/** Second invariant of the deviatoric part of a stress tensor:
    J2 = 1/2 tr(dev.dev)  (dev = sig - 1/3 tr(sig) Id) */
{
#define SIG(i,j) (sig[3*(i)+(j)])
  T j2a = (SIG(0,0) - SIG(1,1))*(SIG(0,0) - SIG(1,1))
             + (SIG(1,1) - SIG(2,2))*(SIG(1,1) - SIG(2,2))
             + (SIG(2,2) - SIG(0,0))*(SIG(2,2) - SIG(0,0)) ;
  T j2b = SIG(0,1)*SIG(1,0) + SIG(1,2)*SIG(2,1) + SIG(2,0)*SIG(0,2) ;
  return(j2a/6. + j2b) ;
#undef SIG
}





/*
   Polynomial roots
*/
template<typename T>
inline int (Math_ComputePolynomialEquationRootsOfDegree2)(T* x)
{
  T b = x[1]/x[0] ;
  T c = x[2]/x[0] ;
  T delta = b*b - 4*c ;
  
  if(c == 0) {
    if(b == 0) {
      x[0] = 0 ;
      return(1) ;
    } else {
      x[0] = MAX(-b,0) ;
      x[1] = MIN(-b,0) ;
      return(2) ;
    }
  } else if(delta == 0) {
    T temp1 = - b*0.5 ;
    x[0] = temp1 ;
    return(1) ;
  } else if(delta > 0) {
    T temp1 = - b*0.5;
    T temp2 = sqrt(delta)*0.5 ;
    if(b < 0) {
      x[0] = temp1 + temp2 ;
      x[1] = c/x[0] ;
    } else {
      x[1] = temp1 - temp2 ;
      x[0] = c/x[1] ;
    }
    return(2) ;
  }
  
  return(0) ;
}



template<typename T>
inline int (Math_ComputePolynomialEquationRootsOfDegree3)(T* x)
{
  T b = x[1]/x[0] ;
  T c = x[2]/x[0] ;
  T d = x[3]/x[0] ;
  /* Depressed cubic equation : x3 + px + q = 0 */
  T p  = c - b*b/3 ;
  T q  = d + 2*b*b*b/27 - b*c/3 ;
  T p3 = p*p*p ;
  T q2 = q*q ;
  T r2 = p3/27 + q2/4 ;
  int n = 0 ;
  
  if(q == 0) {
    if(p >= 0) {
      x[0] = 0 ;
      n = 1 ;
    } else {
      T r = sqrt(-p) ;
      x[0] = r ;
      x[1] = 0 ;
      x[2] = - r ;
      n = 3 ;
    }
    
  } else if(r2 == 0) { /* q != 0 and p3/27 = -q2/4 (< 0) */
    T x1 = 3*q/p ;
    T x2 = - 0.5*x1 ; /* double root */
    x[0] = MAX(x1,x2) ;
    x[1] = MIN(x1,x2) ;
    n = 2 ;
    
  } else if(r2 > 0) {
    T r  = sqrt(r2) ;
    T r1 = (q < 0) ? r : -r ;
    T u3 = - 0.5*q + r1 ; /* so u3 is != 0 even for p = 0 */
    T u  = (u3 > 0) ? pow(u3,1/3.) : - pow(-u3,1/3.) ;
    x[0] = u - p/(3*u) ;
    n = 1 ;
    
  /* Trigonometric solution */
  } else { /* q != 0 and p3/27 < -q2/4 (< 0) */
    T srp3 = sqrt(-p/3) ;
    T k    = acos(1.5*q/(p*srp3))/3 ; /* 0 < k < pi/3 */
    T ck   = srp3*cos(k) ;
    T sk   = srp3*sin(k) ;
    T sr3  = sqrt(3.) ;
    T x1   = 2*ck ;
    T x2   = - ck + sr3*sk ;
    T x3   = - ck - sr3*sk ;
    x[0] = x1 ;
    x[1] = x2 ;
    x[2] = x3 ;
    n = 3 ;
  }
  
  {
    int i ;
    for(i = 0 ; i < n ; i++) x[i] -= b/3 ;
  }
  
  return(n) ;
}



template<typename T>
inline int (Math_ComputePolynomialEquationRootsOfDegree4)(T* x)
{
  T b = x[1]/x[0] ;
  T c = x[2]/x[0] ;
  T d = x[3]/x[0] ;
  T e = x[4]/x[0] ;
  /* Depressed quartic equation: x4 + px2 + qx + r = 0 */
  T b2 = b*b ;
  T b3 = b2*b ;
  T b4 = b2*b2 ;
  T p  = - 0.375*b2 + c ;
  T q  =   0.125*b3 - 0.5*b*c + d ;
  T r  = - 0.01171875*b4 + 0.0625*b2*c - 0.25*b*d + e ;
  int n = 0 ;
  
  if(e == 0) {
    int i ;
    n = Math_ComputePolynomialEquationRootsOfDegree3(x) ;
    x[n] = 0 ;
    for(i = n ; i > 0 ; i--) {
      if(x[i - 1] >= 0) break ;
      x[i] = x[i - 1] ;
      x[i - 1] = 0 ;
    }
    return(n + 1) ;
    
  /* Biquadractic equations */
  } else if(q == 0) { 
    T y[3] = {1,0,0} ;
    int i,m ;
    y[1] = p ;
    y[2] = r ;
    m = Math_ComputePolynomialEquationRootsOfDegree2(y) ;
    /* First keep only positive roots */
    for(i = 0 ; i < m ; i++) {
      if(y[i] >= 0) {
        x[i] = sqrt(y[i]) ;
        n++ ;
      }
    }
    /* Second put the negative roots behind */
    for(i = 0 ; i < n ; i++) {
        x[2*n - 1 - i] =  - x[i] ;
    }
    n *= 2 ;
    
  /* Ferrari's solution */
  } else { 
    /* The nested depressed cubic equation: x3 + p1x + q1 */
    T p1 = - p*p/12 - r ;
    T q1 = - p*p*p/108 + p*r/3 - 0.125*q*q ;
    T y[4] = {1,0,0,0} ;
    y[2] = p1 ;
    y[3] = q1 ;
    /* Solving for the nested depressed cubic equation */
    Math_ComputePolynomialEquationRootsOfDegree3(y) ;
    {
      T s = y[0] - 5./6*p ;
      T w2 = p + 2*s ;
      if(w2 <= 0) { /* Theoretically it's never met */
        arret("Math_ComputePolynomialEquationRootsOfDegree4: never met?") ;
        return(0) ;
      } else {
        T w = sqrt(w2) ;
        T h1 = - (3*p + 2*s + 2*q/w) ;
        T h2 = - (3*p + 2*s - 2*q/w) ;
        if(h1 == 0) {
          x[0] = 0.5*w ;
          n = 1 ;
        } else if(h1 > 0) {
          x[0] = 0.5*(w + sqrt(h1)) ;
          x[1] = 0.5*(w - sqrt(h1)) ;
          n = 2 ;
        }
        if(h2 == 0) {
          x[n++] = - 0.5*w ;
        } else if(h2 > 0) {
          x[n++] = 0.5*(- w + sqrt(h2)) ;
          x[n++] = 0.5*(- w - sqrt(h2)) ;
        }
      }
    }
  }
  
  {
    int i ;
    for(i = 0 ; i < n ; i++) x[i] -= 0.25*b ;
  }
  
  return(n) ;
}




template<typename T>
inline int (Math_ComputePolynomialEquationRoots)(T* x,int const n)
/** Real Roots of Polynomial Equation of Degree n
 *  Return the nb of real solutions, sorted 
 *  from max to min */
{
#define MAX_DEGREE  (4)
  //double y[MAX_DEGREE + 1] ;
  int m = 0 ;
  
  if(x[0] != 0) {
    if(n == 0) {
      arret("Math_ComputePolynomialEquationRoots: degree 0") ;
    } else if(n == 1) {
      x[0] = - x[1]/x[0] ;
      return(1) ;
    } else if(n == 2) {
      m = Math_ComputePolynomialEquationRootsOfDegree2(x) ;
      return(m) ;
    } else if(n == 3) {
      m = Math_ComputePolynomialEquationRootsOfDegree3(x) ;
    } else if(n == 4) {
      m = Math_ComputePolynomialEquationRootsOfDegree4(x) ;
    } else {
      arret("Math_ComputePolynomialEquationRoots: degree too big") ;
    }
  } else {
    int i ;
    
    m = Math_ComputePolynomialEquationRoots(x + 1,n - 1) ;
    for(i = 0 ; i < m ; i++) x[i] = x[i + 1] ;
    return(m) ;
  }
  
  return(m) ;
#undef MAX_DEGREE
}




template<typename T>
inline int (Math_PolishPolynomialEquationRoot)(T* x,int n,T* proot,T tol,int iterations)
{
  T root = proot[0] ;
  
  for(int it = 0 ; it < iterations ; it++) {
    T error = x[0] ;
    T derivative = 0 ;
    T droot ;
      
    for(int i = 0 ; i < n ; i++) {
      T a = error ;
      T b = x[i + 1] ;
      
      error = a*root + b ;
    }
      
    for(int i = 0 ; i < n ; i++) {
      T a = derivative ;
      T b = (n - i)*x[i] ;
      
      derivative = a*root + b ;
    }
      
    if(derivative == 0) {
      proot[0] = root ;
      return(0) ;
    }
      
    droot = - error / derivative ;
    root += droot ;
      
    if(fabs(droot) < tol) {
      proot[0] = root ;
      return(0) ;
    }
  }
  
  /* Raise an interrupt signal instead of exit */
  Message_Warning("Math_PolishPolynomialEquationRoot: no convergence") ;
  {    
    Message_Direct("\nthe %dth order polynomial equation:\n",n) ;
    
    for(int i = 0 ; i <= n ; i++) {
      Message_Direct("%d order coefficient: %lf\n",n-i,x[i]) ;
    }
    
    Message_Direct("\nhas not converged for the root %lf\n",root) ;
  }
  //Exception_Interrupt ;
  /*
  arret("Math_PolishPolynomialEquationRoot: not converged ") ;
  */
  return(-1) ;
}


#ifdef __CPLUSPLUS
}
#endif
#endif
