macro(set_extralibs extralibs extralibs_install_rpath)

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
    
    list(APPEND ${extralibs_install_rpath} ${${lib_name}_DIR})
    list(APPEND ${extralibs} ${lib_full_path})
    
    message("${lib_name} library found in " ${${lib_name}_DIR})

    set(HAVE_${lib_name} TRUE)
  endif()
endforeach()

list(REMOVE_DUPLICATES ${extralibs_install_rpath})
list(REMOVE_DUPLICATES ${extralibs})

if(NOT "${${extralibs}}" STREQUAL "")
  message("Some extra-libraries are used.")
  message("Full paths to these extra-libraries:")
  message("${${extralibs}}")
else()
  message("No extra-libraries are used")
endif()

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
endmacro()
