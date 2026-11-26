# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# fft992 library
# Legacy FFT implementations (FFT992 and Bluestein algorithms)

function( generate_file )
  set( options )
  set( oneValueArgs INPUT OUTPUT BACKEND )
  set( multiValueArgs )
  cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set( output  ${_PAR_OUTPUT} )
  set( input   ${_PAR_INPUT} )
  set( backend ${_PAR_BACKEND} )
  set( sed_rules ${PROJECT_SOURCE_DIR}/src/sedrenames.txt )

  set( JPRB_dp JPRD )
  set( jprb_dp jprd )
  set( JPRB_sp JPRM )
  set( jprb_sp jprm )
  set( JPRB_gpu_dp JPRD )
  set( jprb_gpu_dp jprd )
  set( JPRB_gpu_sp JPRM )
  set( jprb_gpu_sp jprm )

  add_custom_command(
    OUTPUT  ${output}
    COMMAND  cat ${sed_rules} |
             sed -e "s/VARIANTDESIGNATOR/${backend}/g" |
             sed -e "s/TYPEDESIGNATOR_UPPER/${JPRB_${backend}}/g" |
             sed -e "s/TYPEDESIGNATOR_LOWER/${jprb_${backend}}/g" |
             sed -rf - ${input} > ${output}
    DEPENDS  ${input} ${sed_rules}
    COMMENT "Generating ${output}"
    VERBATIM
  )
endfunction( generate_file )

set(fft992_sources
  fft992.F90
  fft992_cc.F90
  set99.F90
  set99b.F90
  bluestein_mod.F90
)

foreach( prec dp sp )
  if( HAVE_${prec} )
    set( target fft992_${prec} )

    set( destination ${CMAKE_CURRENT_BINARY_DIR}/generated/fft992_${prec} )
    file( MAKE_DIRECTORY ${destination} )
    set( generated )
    foreach( file_i ${fft992_sources} )
      get_filename_component( outfile_name    ${file_i} NAME )
      set( outfile ${destination}/${outfile_name} )
      ecbuild_debug( "Generate ${outfile}" )
      generate_file( BACKEND ${prec} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile} )
      list( APPEND generated ${outfile} )
    endforeach( file_i )

    ecbuild_add_library(
      TARGET ${target}
      SOURCES ${generated}
      PUBLIC_LIBS parkind_${prec} fiat
      INSTALL_HEADERS LISTED
      HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/fft992
    )

    ecbuild_target_fortran_module_directory(
      TARGET ${target}
      MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${target}
      INSTALL_MODULE_DIRECTORY module/${target}
    )

    target_include_directories(${target} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/ectrans>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/fiat>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/${target}>
      $<INSTALL_INTERFACE:module/ectrans>
      $<INSTALL_INTERFACE:module/fiat>
      $<INSTALL_INTERFACE:module/${target}>
    )
  endif()
endforeach()
