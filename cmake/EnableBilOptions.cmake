file(STRINGS ${BIL_PATH}/OPTIONS opts REGEX "^[^# ]")

foreach(opt IN ITEMS ${opts})
  #message("opt = ${opt}")
  string(REGEX MATCH "^[^ ]*" opt_enablename "${opt}")
  string(REGEX REPLACE "^[ ]*[^ ]+[ ]+=[ ]+|[ ]+$" "" opt_bool "${opt}")

  #message("opt_enablename = ${opt_enablename}")
  #message("opt_bool = ${opt_bool}")
  
  if("${opt_bool}")
    string(REGEX REPLACE "ENABLE_" "" opt_name "${opt_enablename}")
    
    if(NOT "${opt_name}" STREQUAL "")
      set(ENABLE_${opt_name} ON CACHE BOOL "" FORCE)
    endif()
  endif()
endforeach()
