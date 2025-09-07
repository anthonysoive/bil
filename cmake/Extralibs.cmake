
# Re-build the library file
#[[
message("Updating the library file: ${BIL_PATH}/BilExtraLibs.h")
file(REMOVE ${BIL_PATH}/BilExtraLibs.h)
file(TOUCH  ${BIL_PATH}/BilExtraLibs.h)

file(APPEND ${BIL_PATH}/BilExtraLibs.h "#ifndef BILEXTRALIBS_H\n")
file(APPEND ${BIL_PATH}/BilExtraLibs.h "#define BILEXTRALIBS_H\n")

file(APPEND ${BIL_PATH}/BilExtraLibs.h "#include \"BilConfig.h\"\n")
#]]

set(BILEXTRALIBS_INSTALL_RPATH)
set(BIL_EXTRALIBS)

file(STRINGS ${BIL_PATH}/EXTRALIBS libs REGEX "^[^# ]")

foreach(lib IN ITEMS ${libs})
  #message("lib = ${lib}")
  string(REGEX MATCH "^[A-Z,0-9]+" lib_name "${lib}")
  string(REGEX REPLACE "^[ ]*[A-Z,0-9]+[ ]+=[ ]+|[ ]+$" "" lib_path "${lib}")

  #message("lib_name = ${lib_name}")
  #message("lib_path = ${lib_path}")
  
  # This is an alternative to option, only if ENABLE_${lib_name} is false.
  if(EXISTS "${lib_path}" AND NOT ENABLE_${lib_name})
  
    get_filename_component(lib_full_path ${lib_path} ABSOLUTE)
    get_filename_component(${lib_name}_DIR ${lib_full_path} DIRECTORY)
    
    list(APPEND BILEXTRALIBS_INSTALL_RPATH ${${lib_name}_DIR})
    list(APPEND BIL_EXTRALIBS ${lib_full_path})
    
    #file(APPEND ${BIL_PATH}/BilExtraLibs.h "#define ${lib_name}LIB  ${${lib_name}_DIR}\n")
    
    message("${lib_name} library found in " ${${lib_name}_DIR})
    #message("is added to the cache")

    #file(APPEND ${BIL_PATH}/BilExtraLibs.h "#define HAVE_${lib_name}\n")
    set(HAVE_${lib_name} TRUE)
  endif()
endforeach()

list(REMOVE_DUPLICATES BILEXTRALIBS_INSTALL_RPATH)
list(REMOVE_DUPLICATES BIL_EXTRALIBS)

#file(APPEND ${BIL_PATH}/BilExtraLibs.h "#endif")


if(NOT "${BIL_EXTRALIBS}" STREQUAL "")
  message("Some extra-libraries are used.")
  message("Full paths to these extra-libraries:")
  message("${BIL_EXTRALIBS}")
else()
  message("No extra-libraries are used")
endif()


#message(BILEXTRALIBS_INSTALL_RPATH = "${BILEXTRALIBS_INSTALL_RPATH}")

#[[
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
]]

#[[
#message(CMAKE_INSTALL_RPATH_USE_LINK_PATH = ${CMAKE_INSTALL_RPATH_USE_LINK_PATH})
message(CMAKE_INSTALL_RPATH = "${CMAKE_INSTALL_RPATH}")
]]


#welcome:
#       read -p "Do you use BLAS library (Y/N)?:" BLAS_USE
