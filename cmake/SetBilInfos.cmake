# Bil version
# -----------
file(STRINGS ${BIL_PATH}/VERSION BIL_VERSION)
#message("BIL_VERSION = ${BIL_VERSION}")


# General informations
# --------------------
set(BIL_SHORT_LICENSE "GNU General Public License")
string(TIMESTAMP BIL_DATE)
string(TIMESTAMP BIL_YEAR %Y)
cmake_host_system_information(RESULT BIL_HOST QUERY HOSTNAME)
#set(BIL_HOST    ${shell hostname}: ${shell hostname -I})
#cmake_host_system_information(RESULT BIL_DISTRIB QUERY DISTRIB_INFO)
set(BIL_URL       "http://bil.ifsttar.fr")
set(BIL_EMAIL     "patrick.dangla@univ-eiffel.fr")
set(BIL_COPYRIGHT "Copyright \(C\) 2002 Patrick Dangla")
set(BIL_PROGNAME  "Bil, a modeling platform based on FEM/FVM")
execute_process(COMMAND date "+%Y%m%d" OUTPUT_VARIABLE BIL_DATE 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND hostname OUTPUT_VARIABLE BIL_HOSTNAME 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE BIL_PACKAGER 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

       
# The OS
# ------
if(APPLE)
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
    set(BIL_OS "MacOSARM")
  else()
    set(BIL_OS "MacOSX")
  endif()
elseif(CYGWIN OR MSYS)
  # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
  # as a pure Windows build and make cmake find pure Windows import libraries
  # (.lib)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
     CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                    OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CXX_COMPILER_MACHINE MATCHES "mingw")
      set(BIL_OS "Windows")
      set(WIN32 1)
      add_definitions(-DWIN32)
      #set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")
      #set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".lib" ".LIB" ".dll" ".DLL" ".dll.a")
    endif()
  endif()
endif()
if(NOT BIL_OS)
  #set(BIL_OS "${CMAKE_SYSTEM_NAME}")
  set(BIL_OS ${CMAKE_HOST_SYSTEM_NAME})
endif()


               

configure_file(${BIL_PATH}/BilInfo.h.in 
               ${BIL_PATH}/BilInfo.h)
configure_file(${BIL_PATH}/BilVersion.h.in
               ${BIL_PATH}/BilVersion.h)
configure_file(${BIL_PATH}/BilPath.h.in
               ${BIL_PATH}/BilPath.h)
