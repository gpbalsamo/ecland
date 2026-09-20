MODULE YOMFORC1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE
SAVE

!     -----------------------------------------------

!*    ATMOSPHERIC FORCING DATA FOR 1D SURFACE SCHEME, PILPS VERSION

!     -----------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPTYFC=12

REAL(KIND=JPRB),ALLOCATABLE,TARGET:: GFOR(:,:,:)
REAL(KIND=JPRB),POINTER :: UFI(:,:)
REAL(KIND=JPRB),POINTER :: VFI(:,:)
REAL(KIND=JPRB),POINTER :: TFI(:,:)
REAL(KIND=JPRB),POINTER :: QFI(:,:)
REAL(KIND=JPRB),POINTER :: CO2FI(:,:)
REAL(KIND=JPRB),POINTER :: PSFI(:,:)
REAL(KIND=JPRB),POINTER :: SRFFI(:,:)
REAL(KIND=JPRB),POINTER :: TRFFI(:,:)
REAL(KIND=JPRB),POINTER :: R30FI(:,:)
REAL(KIND=JPRB),POINTER :: S30FI(:,:)
REAL(KIND=JPRB),POINTER :: R30FI_C(:,:)
REAL(KIND=JPRB),POINTER :: S30FI_C(:,:)
REAL(KIND=JPRD) :: DTIMFC                ! Forcing frequency 
REAL(KIND=JPRD) :: RTSTFC                ! first step of loaded forcing ref. time 
INTEGER(KIND=JPIM) :: NSTPFC
INTEGER(KIND=JPIM) :: JPSTPFC           ! Dimension of the forcing set in the namelist
INTEGER(KIND=JPIM) :: DIMFORC            ! forcing dimension 1/2

! Windowed (streaming) forcing reads -- see NAMDIM's own comment in
! namdim1s.h, rdfvar.F90 and reload_forc1s.F90. NFORCWINDOW=0 (the default) is
! the original, unwindowed behaviour: JPSTPFC=NDFORC, GFOR holds the whole
! run's forcing, RDFVAR is called exactly once (from SUFCDF), and DTFORC's
! refill guard never triggers -- byte-for-byte identical to before this
! parameter existed. NFORCWINDOW>0 makes JPSTPFC the (smaller) buffer size;
! DTFORC calls RELOAD_FORC1S to keep it topped up as the run advances.
INTEGER(KIND=JPIM) :: NFORCWINDOW
! Set once the loaded window stops advancing (i.e. the forcing file has no
! further records): without it, a short final window would re-trigger
! DTFORC's refill guard on every remaining timestep -- still correct, but
! it would re-read the same tail repeatedly. Only ever consulted when
! NFORCWINDOW>0.
LOGICAL :: LFORCEOF

! Module-persistent mirrors of NAMFORC's per-variable file names (SUFCDF's
! own locals, declared via #include "namforc1s.h", are not visible outside
! that subroutine). Set once by SUFCDF right after it reads NAMFORC; read
! by RELOAD_FORC1S (reload_forc1s.F90) on every windowed refill, so the actual
! file-name namelist read itself only ever happens once, as before.
CHARACTER(LEN=255) :: CFFORCU,CFFORCV,CFFORCT,CFFORCQ,CFFORCC,CFFORCP,&
                      &CFFORCRAIN,CFFORCSNOW,CFFORCSW,CFFORCLW
! REAL(KIND=JPRB) :: PHISTA

!------------------------------------------------------------     
END MODULE YOMFORC1S
