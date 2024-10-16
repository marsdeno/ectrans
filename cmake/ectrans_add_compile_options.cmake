# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.



if( NOT DEFINED ECTRANS_HAVE_CONTIGUOUS_ISSUE )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel"  )
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS_EQUAL 19)
      set( ECTRANS_HAVE_CONTIGUOUS_ISSUE True )
    endif()
  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU"  )
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL "9.2"
     OR CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL "12.2.0" )
      set( ECTRANS_HAVE_CONTIGUOUS_ISSUE True )
    endif()
  endif()
endif()

macro( ectrans_add_compile_options )
  set( options NOFAIL )
  set( single_value_args FLAGS BEHAVIOUR)
  set( multi_value_args SOURCES )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to ectrans_add_compile_flags(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()
  if(NOT _PAR_SOURCES)
    ecbuild_critical("SOURCES keyword missing to ectrans_add_compile_flags()")
  endif()
  if(NOT _PAR_FLAGS)
    ecbuild_critical("FLAGS keyword missing to ectrans_add_compile_flags()")
  endif()
  if(NOT _PAR_BEHAVIOUR)
    ecbuild_info("BEHAVIOUR keyword not provided to ectrans_add_compile_flags(), defaulting to APPEND")
    set( BEHAVIOUR "APPEND" )
  elseif( NOT ${_PAR_BEHAVIOUR} MATCHES "APPEND|REPLACE")
    ecbuild_critical("BEHAVIOUR keyword for ectrans_add_compile_flags() must be APPEND or REPLACE")
  else()
    set( BEHAVIOUR ${_PAR_BEHAVIOUR} )
  endif()

  if ( ${BEHAVIOUR} MATCHES "APPEND" )
    set( BASE_FLAGS ${ECTRANS_Fortran_FLAGS_${CMAKE_BUILD_TYPE_CAPS}} )
  else()
    set( BASE_FLAGS "" )
  endif()

  foreach( _file ${_PAR_SOURCES} )
    ecbuild_warn("Adding custom compile flags for file ${_file} : [${_PAR_FLAGS}]")
    set_source_files_properties( ${_file} PROPERTIES COMPILE_FLAGS "${BASE_FLAGS} ${_PAR_FLAGS}" )
  endforeach()
endmacro()

