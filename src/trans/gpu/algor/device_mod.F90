! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DEVICE_MOD

#ifdef CUDAGPU
#define hipDeviceSynchronize cudaDeviceSynchronize
#define hipStreamSynchronize cudaStreamSynchronize
#define hipStreamDestroy cudaStreamDestroy
#define hipSetDevice cudaSetDevice
#define hipGetDevice cudaGetDevice
#define hipGetDeviceCount cudaGetDeviceCount
#endif

INTERFACE DEVICE_SYNC

INTEGER FUNCTION DEVICE_SYNCHRONIZE() BIND(C, NAME='hipDeviceSynchronize')
END FUNCTION DEVICE_SYNCHRONIZE

END INTERFACE DEVICE_SYNC

INTERFACE DEVICESTREAMSYNC

INTEGER FUNCTION DEVICE_STREAM_SYNCHRONIZE(STREAM) BIND(C, NAME='hipStreamSynchronize')
USE ISO_C_BINDING, ONLY: C_PTR
TYPE(C_PTR) :: STREAM
END FUNCTION DEVICE_STREAM_SYNCHRONIZE

END INTERFACE DEVICESTREAMSYNC

INTERFACE DEVICESTREAMDESTROY

INTEGER FUNCTION DEVICE_STREAM_DESTROY(STREAM) BIND(C, NAME='hipStreamDestroy')
USE ISO_C_BINDING, ONLY: C_PTR
TYPE(C_PTR) :: STREAM
END FUNCTION DEVICE_STREAM_DESTROY

END INTERFACE DEVICESTREAMDESTROY

INTERFACE DEVICESETDEVICE

INTEGER FUNCTION DEVICE_SETDEVICE(DEVNUM) BIND(C, NAME='hipSetDevice')
USE ISO_C_BINDING, ONLY: C_INT
INTEGER(C_INT), VALUE :: DEVNUM
END FUNCTION DEVICE_SETDEVICE

END INTERFACE DEVICESETDEVICE

INTERFACE DEVICEGETDEVICE

INTEGER FUNCTION DEVICE_GETDEVICE(DEVNUM) BIND(C, NAME='hipGetDevice')
USE ISO_C_BINDING, ONLY: C_INT
INTEGER(C_INT) :: DEVNUM
END FUNCTION DEVICE_GETDEVICE

END INTERFACE DEVICEGETDEVICE

INTERFACE DEVICEGETDEVICECOUNT

INTEGER FUNCTION DEVICE_GETDEVICECOUNT(DEVNUM) BIND(C, NAME='hipGetDeviceCount')
USE ISO_C_BINDING, ONLY: C_INT
INTEGER(C_INT) :: DEVNUM
END FUNCTION DEVICE_GETDEVICECOUNT

END INTERFACE DEVICEGETDEVICECOUNT

INTERFACE DEVICEGETMEMINFO

INTEGER FUNCTION DEVICE_MEMGETINFO(MEMFREE_MB, MEMTOTAL_MB) BIND(C, NAME='c_hipmemgetinfo')
USE ISO_C_BINDING, ONLY: C_INT
INTEGER(C_INT) :: MEMFREE_MB, MEMTOTAL_MB
END FUNCTION DEVICE_MEMGETINFO

END INTERFACE DEVICEGETMEMINFO

END MODULE DEVICE_MOD
