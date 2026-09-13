! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!> Live fuel moisture content (LEFIRE): purely diagnostic, recomputed every
!> timestep from the current soil moisture and LAI (LFMC_MOD), with no memory
!> of its own -- unlike the dead-fuel-moisture and fuel-load reservoirs, it is
!> never read back, so it does not need the full prognostic-state/restart
!> machinery those use (see yos_soil.F90's LEFIRE comment). This module is
!> just a plain shared buffer from where CALLPAR1S computes it (LFMC_MOD) to
!> where the diagnostic-output code (UPDDIAG) reads it for o_fire.nc.

MODULE YOS_FIRE_DIAG_MOD

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE
SAVE

REAL(KIND=JPRB), ALLOCATABLE :: ZLFMCL(:,:) ! live fuel moisture content, low veg (%), (NPROMA,NBLOCKS)
REAL(KIND=JPRB), ALLOCATABLE :: ZLFMCH(:,:) ! live fuel moisture content, high veg (%), (NPROMA,NBLOCKS)

CONTAINS

SUBROUTINE FIRE_DIAG_ENSURE_ALLOC(KPROMA,KBLOCKS)
INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN) :: KBLOCKS
IF (.NOT. ALLOCATED(ZLFMCL)) THEN
  ALLOCATE(ZLFMCL(KPROMA,KBLOCKS))
  ALLOCATE(ZLFMCH(KPROMA,KBLOCKS))
  ZLFMCL(:,:) = 0.0_JPRB
  ZLFMCH(:,:) = 0.0_JPRB
ENDIF
END SUBROUTINE FIRE_DIAG_ENSURE_ALLOC

END MODULE YOS_FIRE_DIAG_MOD
