# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


if( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )

  set( ECTRANS_Fortran_FLAGS "-fp-model precise -fp-speculation=safe -fast-transcendentals -assume norealloc_lhs")
  set( ECTRANS_Fortran_FLAGS_BIT "-march=core-avx2 -no-fma -O2 -align array32byte")
  set( ECTRANS_Fortran_FLAGS_RELEASE "-march=core-avx2 -Ofast")
  set( ECTRANS_Fortran_FLAGS_RELWITHDEBINFO "-g1 ${FOO_Fortran_FLAGS_RELEASE}")
  set( ECTRANS_Fortran_FLAGS_DEBUG "-g -traceback -warn all")


endif()

