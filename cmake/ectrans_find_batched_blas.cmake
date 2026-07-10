# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ectrans_find_batched_blas )
  if( DEFINED HAVE_BATCHED_BLAS AND NOT HAVE_BATCHED_BLAS )
    message( STATUS "Batched BLAS disabled (ENABLE_BATCHED_BLAS=OFF or default) - using sequential GEMM loop" )
  else()
    set( HAVE_BATCHED_BLAS OFF )

    if( HAVE_MKL )
      # MKL always supports batched BLAS
      set( HAVE_BATCHED_BLAS ON )
    else()
      # Try to detect batched BLAS (cblas_dgemm_batch) in LAPACK_LIBRARIES
      if( LAPACK_LIBRARIES )
        set( CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} )
        include( CheckSymbolExists )
        check_symbol_exists( cblas_dgemm_batch "cblas.h" _HAVE_CBLAS_DGEMM_BATCH )
        if( _HAVE_CBLAS_DGEMM_BATCH )
          set( HAVE_BATCHED_BLAS ON )
        endif()
      endif()
    endif()

    if( HAVE_BATCHED_BLAS )
      message( "Batched BLAS (cblas_dgemm_batch) detected" )
    else()
      message( "Batched BLAS NOT detected - using sequential GEMM loop" )
    endif()
  endif()

  ecbuild_debug_var( HAVE_BATCHED_BLAS )
endmacro()
