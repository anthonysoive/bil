# OS
#----
UNAME = Linux


# Compilers (C++11 for variadic macros, C++17 for autodiff)
#----------
CC     := gcc -x c   -std=c11
CPP    := gcc -x c++ -std=gnu++17
# Force c++ compiler
CC     := ${CPP}
F77    := gfortran -cpp
LEX    := flex #-d
YACC   := bison #-t


# Compiler and preprocessor flags (recommended by GSL)
# see www.gnu.org/software/gsl/manual/html_node/GCC-warning-options-for-numerical-programs.html)


# Warning flags
#--------------
# Warning flags common to C/C++
WFLAGSmini := -Wpedantic -Wall -Wextra

WFLAGS := ${WFLAGSmini} -Wformat=2 -Wconversion -Wtrampolines -Wimplicit-fallthrough -Wcast-align -Wcast-qual -Wpointer-arith -Wshadow -Wsign-conversion

#WFLAGS := -Wall -Wcast-align -Wcast-qual -Wnested-externs -Wpointer-arith -Wmissing-prototypes -Wstrict-prototypes -Wshadow -fno-common -Dinline=  -Wvariadic-macros #-W #-Wtraditional #-Wconversion #-Wwrite-strings #-Werror #-fshort-enums


# Warning flags for C
WCFLAGS := ${WFLAGS} -Wnested-externs -Wmissing-prototypes -Wstrict-prototypes

# Warning flags for C++
WCXXFLAGS := ${WFLAGS} -Wnon-virtual-dtor -Wzero-as-null-pointer-constant -Woverloaded-virtual -Wfloat-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wno-unused


# Directory flags (-I replaced by -iquote)
#----------------
# DFLAGS = -I${BIL_PATH}
DFLAGS = -iquote${BIL_PATH} -D '__BASEFILE__=${basename ${notdir $@}}'


# Preprocessor flags
#-------------------
IFLAGS-SRC = -iprefix ${BIL_PATH}/src/ -iwithprefix . -iwithprefix Common -iwithprefix DataSet -iwithprefix Help -iwithprefix Libraries -iwithprefix Macros -iwithprefix Main -iwithprefix Models -iwithprefix Modules -iwithprefix Outputs -iwithprefix Solver -iwithprefix Tools -iwithprefix Parser -iwithprefix Parallelization

# Additional flags for solver
IFLAGS-SOLVER = -iprefix ${BIL_PATH}/src/Solver/ -iwithprefix . -iwithprefix MatrixStorageFormats -iwithprefix ResolutionMethods -iwithprefix VectorStorageFormats

# Additional flags for models
IFLAGS-MODELS = -iprefix ${BIL_PATH}/src/Models/ -iwithprefix .  -iwithprefix DataBases  -iwithprefix Libraries  -iprefix ${BIL_PATH}/src/Models/Methods/ -iwithprefix .  -iwithprefix MaterialPointMethod -iwithprefix FiniteElementMethod -iwithprefix FiniteVolumeMethod -iprefix ${BIL_PATH}/src/Models/ConstitutiveLaws/ -iwithprefix . -iwithprefix PlasticityModels -iwithprefix DamageModels

# Additional flags for modules
IFLAGS-MODULES =  -iprefix ${BIL_PATH}/src/Modules/ -iwithprefix .

# Additional flags for tools
IFLAGS-TOOLS = -iprefix ${BIL_PATH}/src/Tools/ -iwithprefix . -iwithprefix Libraries/Expression


IFLAGS = ${IFLAGS-SRC} ${IFLAGS-SOLVER} ${IFLAGS-MODELS} ${IFLAGS-MODULES} ${IFLAGS-TOOLS}


# Option flags
#-------------

sinclude ${BIL_PATH}/make.config
OFLAGS :=
#${info ENABLE_OPENMP = ${ENABLE_OPENMP}}
ifdef ENABLE_OPENMP
  ifeq (${ENABLE_OPENMP},ON)
    OFLAGS += -fopenmp
    #${info OFLAGS = ${OFLAGS}}
  endif
endif


# C/CPP compiler flags
#---------------------

# C compiler flags
CFLAGS = -pedantic ${DFLAGS} ${WCFLAGS} ${IFLAGS} -fPIC ${OFLAGS} -fpermissive

# C++ compiler flags
# fpermissive is not satisfactory but at the present time it is needed
# because most of the implementations adopt a C-style which are not 
# allowed in C++.
CXXFLAGS = -pedantic ${DFLAGS} ${WCXXFLAGS} ${IFLAGS} -fPIC ${OFLAGS} -fpermissive
CFLAGS = ${CXXFLAGS}

# C/CPP compiler flags for models
CFLAGS-MODELS = ${CFLAGS} #${IFLAGS-MODELS}
CXXFLAGS-MODELS = ${CXXFLAGS} #${IFLAGS-MODELS}


# C/CPP compiler flags for modules
CFLAGS-MODULES = ${CFLAGS} #${IFLAGS-MODULES}
CXXFLAGS-MODULES = ${CXXFLAGS} #${IFLAGS-MODULES}


# Fortran compiler flags
#-----------------------
FFLAGS = -fPIC


# Optimization flags
#------------------
OPTIM  = -gdwarf-2 -g3     # Debug
#OPTIM  = -O3 -DNDEBUG      # Release
#OPTIM  = -O3 -g -gdwarf-3  # RelWithDebInfo




# Linker
#-------
LINKER = g++ ${CXXFLAGS}



# Linker flags
#-------------
LFLAGS = -lgfortran ${OFLAGS}




# Libraries
#----------

# Library extensions. Depends on which kind of library you want to install
# Choose either "so" (shared object) or "a" (static)
LIBBILEXT   = so

# How you create a static library on this machine
AR     = /usr/bin/ar ruvs
RANLIB = ranlib


# To create and install a shared library
CCSO = gcc -shared ${CXXFLAGS}
LDCONFIG  = /sbin/ldconfig 

