! (C) Copyright 2024- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE LFMC_MOD
CONTAINS
SUBROUTINE LFMC(KIDIA,KFDIA,LDLAND,KTVL,KTVH,PLAIL,PLAIH,PSSM,PLFMC_L,PLFMC_H)

!     PURPOSE
!     -------
!     THIS ROUTINE CALCULATES LIVE FUEL MOISTURE CONTENT FOR BOTH HIGH AND LOW VEGETATION.
!     PURELY DIAGNOSTIC: RECOMPUTED EACH CALL FROM THE CURRENT SOIL MOISTURE AND LAI, NO
!     PERSISTENT STATE. GATED BY LEFIRE (NAMPARVEG).

!     INTERFACE.
!     ----------
!     CALLED FROM *callpar1s*

!     METHOD.
!     LFMC is calculated based on vegetation type following parameterisation from McNorton et al., 2024

!     REFERENCE.
!     Original    J. McNorton      MAR 2024 (arpifs/sparky/lfmc_mod.F90)
!     Ported to ecLand offline     2026

!==============================================================================

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
LOGICAL,            INTENT(IN)    :: LDLAND(:)
INTEGER(KIND=JPIM), INTENT(IN)    :: KTVL(:)
INTEGER(KIND=JPIM), INTENT(IN)    :: KTVH(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PLAIL(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PLAIH(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PSSM(:,:)          ! soil moisture, layers 1-4 (m3 m-3)
REAL(KIND=JPRB),    INTENT(INOUT) :: PLFMC_L(:)         ! LIVE FUEL MOISTURE CONTENT LOW VEG (%)
REAL(KIND=JPRB),    INTENT(INOUT) :: PLFMC_H(:)         ! LIVE FUEL MOISTURE CONTENT HIGH VEG (%)

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZACO(20),ZBCO(20),ZCCO(20),ZDCO(20),ZECO(20)
REAL(KIND=JPRB) :: ZF1(20),ZF2(20),ZF3(20),ZF4(20)
INTEGER(KIND=JPIM) :: ITYL, ITYH
REAL(KIND=JPRB)    :: ZSML, ZSMH

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LFMC_MOD:LFMC',0,ZHOOK_HANDLE)

! Coefficients per IFS vegetation type (McNorton et al., 2024)
ZACO=(/110.76_JPRB, 197.4184_JPRB, 132.5899_JPRB, 132.5899_JPRB, 131.5411_JPRB, &
     & 131.5411_JPRB, 165.762_JPRB, 242.3179_JPRB,242.3179_JPRB, 110.76_JPRB, &
     & 242.3179_JPRB, 242.3179_JPRB, 242.3179_JPRB, 242.3179_JPRB, 242.3179_JPRB, &
     & 158.773_JPRB, 158.773_JPRB,242.3179_JPRB, 242.3179_JPRB, 242.3179_JPRB/)

ZBCO=(/1500.93_JPRB, 119.7407_JPRB, 65.24821_JPRB, 65.24821_JPRB, 78.71811_JPRB, &
     & 78.71811_JPRB, 148.2193_JPRB, 237.92_JPRB,237.92_JPRB, 1500.93_JPRB, &
     & 237.92_JPRB, 237.92_JPRB, 237.92_JPRB, 237.92_JPRB, 237.92_JPRB, 131.0207_JPRB, &
     & 131.0207_JPRB,237.92_JPRB, 237.92_JPRB, 237.92_JPRB/)

ZCCO=(/27.15421_JPRB, 0.2343546_JPRB, 1.902277_JPRB, 1.902277_JPRB, 1.459748e-14_JPRB, &
     & 1.459748e-14_JPRB, 4.305513_JPRB,1.774242_JPRB, 1.774242_JPRB, 27.15421_JPRB, &
     & 1.774242_JPRB, 1.774242_JPRB, 1.774242_JPRB, 1.774242_JPRB, 1.774242_JPRB, &
     & 3.50411_JPRB, 3.50411_JPRB, 1.774242_JPRB, 1.774242_JPRB, 1.774242_JPRB/)

ZDCO=(/4.476948e-24_JPRB, 0.2264072_JPRB, 3.360863e-11_JPRB, 3.360863e-11_JPRB, &
     & 1.213111e-20_JPRB, 1.213111e-20_JPRB,3.222054e-12_JPRB, 0.0326482_JPRB, &
     & 0.0326482_JPRB, 4.476948e-24_JPRB, 0.0326482_JPRB, 0.0326482_JPRB, &
     & 0.0326482_JPRB,0.0326482_JPRB, 0.0326482_JPRB, 1.234682e-10_JPRB, &
     & 1.234682e-10_JPRB, 0.0326482_JPRB, 0.0326482_JPRB, 0.0326482_JPRB/)

ZECO=(/2.207181e-15_JPRB, 0.5289065_JPRB, 1.237002_JPRB, 1.237002_JPRB, 3.827632_JPRB, &
     & 3.827632_JPRB, 0.5709138_JPRB,6.274949e-12_JPRB, 6.274949e-12_JPRB, &
     & 2.207181e-15_JPRB, 6.274949e-12_JPRB, 6.274949e-12_JPRB, 6.274949e-12_JPRB, &
     & 6.274949e-12_JPRB,6.274949e-12_JPRB, 1.722875_JPRB, 1.722875_JPRB, &
     & 6.274949e-12_JPRB, 6.274949e-12_JPRB, 6.274949e-12_JPRB/)

! Root-zone soil-moisture weights per layer (McNorton et al., 2024), calibrated
! for a 4-layer soil column. Not re-derived for NCSS=9/14 configurations.
ZF1 = (/24.0_JPRB, 35.0_JPRB, 26.0_JPRB, 26.0_JPRB, 24.0_JPRB, 25.0_JPRB, 27.0_JPRB, &
      & 100.0_JPRB, 47.0_JPRB, 24.0_JPRB,17.0_JPRB,100.0_JPRB, 25.0_JPRB, 100.0_JPRB, &
      & 100.0_JPRB, 23.0_JPRB, 23.0_JPRB, 19.0_JPRB, 19.0_JPRB, 100.0_JPRB/)
ZF2 = (/41.0_JPRB, 38.0_JPRB, 39.0_JPRB, 38.0_JPRB, 38.0_JPRB, 34.0_JPRB, 27.0_JPRB, &
      & 0.0_JPRB, 45.0_JPRB, 41.0_JPRB,31.0_JPRB,0.0_JPRB, 34.0_JPRB, 0.0_JPRB, &
      & 0.0_JPRB, 36.0_JPRB, 36.0_JPRB, 35.0_JPRB, 35.0_JPRB, 0.0_JPRB/)
ZF3 = (/31.0_JPRB, 23.0_JPRB, 29.0_JPRB, 29.0_JPRB, 31.0_JPRB, 27.0_JPRB, 27.0_JPRB, &
      & 0.0_JPRB, 8.0_JPRB, 31.0_JPRB, 33.0_JPRB,0.0_JPRB, 27.0_JPRB, 0.0_JPRB, &
      & 0.0_JPRB, 30.0_JPRB, 30.0_JPRB, 36.0_JPRB, 36.0_JPRB, 0.0_JPRB/)
ZF4 = (/4.0_JPRB, 4.0_JPRB, 6.0_JPRB, 7.0_JPRB, 7.0_JPRB, 14.0_JPRB, 19.0_JPRB, &
      & 0.0_JPRB, 0.0_JPRB, 4.0_JPRB, 19.0_JPRB, 0.0_JPRB, 14.0_JPRB, 0.0_JPRB, &
      & 0.0_JPRB, 11.0_JPRB, 11.0_JPRB, 10.0_JPRB, 10.0_JPRB, 0.0_JPRB/)

DO JL=KIDIA,KFDIA
 IF (.NOT.LDLAND(JL)) THEN
  PLFMC_L(JL) = 0.0_JPRB
  PLFMC_H(JL) = 0.0_JPRB
 ELSE
  ITYL = KTVL(JL)
  ITYH = KTVH(JL)

! KTVL/KTVH==0 means "no vegetation of this type" (e.g. bare-soil/desert
! sites, common in PLUMBER2's arid sites -- AU-ASM, US-SRM/SRG/Wkg/Whs,
! US-Myb, US-FPe) -- a legitimate value, not a missing-data sentinel, but
! ZF1..ZF4/ZACO..ZECO are 1-indexed per-type tables with no slot 0, so
! indexing them with ITYL/ITYH==0 is an out-of-bounds access (crashes with
! "Subscript #1 of the array ZF1 has value 0"). Guard the same way
! FUEL_MOD's KTVL(JL) > 0 check already does for the sibling LEFIRE fuel
! calculation -- no live fuel of a vegetation type that isn't present.
  IF((PSSM(JL,1)+PSSM(JL,2)+PSSM(JL,3)+PSSM(JL,4)) > 1E-7_JPRB) THEN
   IF (ITYL > 0) THEN
    ZSML = MAX((ZF1(ITYL)*PSSM(JL,1)+ZF2(ITYL)*PSSM(JL,2)+ZF3(ITYL)*PSSM(JL,3)+ZF4(ITYL)*PSSM(JL,4))/100.0_JPRB,1E-7_JPRB)
    PLFMC_L(JL) = MAX(ZACO(ITYL) - ZBCO(ITYL)*EXP(-(ZCCO(ITYL)*ZSML+ZDCO(ITYL)*PLAIL(JL)+ZECO(ITYL)*PLAIL(JL)*ZSML)),30.0_JPRB)
   ELSE
    PLFMC_L(JL) = 0.0_JPRB
   ENDIF
   IF (ITYH > 0) THEN
    ZSMH = MAX((ZF1(ITYH)*PSSM(JL,1)+ZF2(ITYH)*PSSM(JL,2)+ZF3(ITYH)*PSSM(JL,3)+ZF4(ITYH)*PSSM(JL,4))/100.0_JPRB,1E-7_JPRB)
    PLFMC_H(JL) = MAX(ZACO(ITYH) - ZBCO(ITYH)*EXP(-(ZCCO(ITYH)*ZSMH+ZDCO(ITYH)*PLAIH(JL)+ZECO(ITYH)*PLAIH(JL)*ZSMH)),30.0_JPRB)
   ELSE
    PLFMC_H(JL) = 0.0_JPRB
   ENDIF
  ELSE
!  Must be a real branch, not MERGE: MERGE is a function, so Fortran does not
!  guarantee short-circuiting and compilers routinely evaluate BOTH arguments
!  before applying the mask -- MERGE(ZACO(ITYL),...,ITYL > 0) therefore still
!  evaluates ZACO(0) at a point with no vegetation of that type, which is the
!  very out-of-bounds access the ITYL/ITYH guards exist to prevent. Reached
!  when SUM(PSSM(JL,1:4)) <= 1E-7, i.e. essentially-dry soil, which is
!  correlated with (not independent of) KTV==0 at bare-soil/desert points.
   IF (ITYL > 0) THEN
    PLFMC_L(JL) = ZACO(ITYL)
   ELSE
    PLFMC_L(JL) = 0.0_JPRB
   ENDIF
   IF (ITYH > 0) THEN
    PLFMC_H(JL) = ZACO(ITYH)
   ELSE
    PLFMC_H(JL) = 0.0_JPRB
   ENDIF
  ENDIF
 ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('LFMC_MOD:LFMC',1,ZHOOK_HANDLE)

END SUBROUTINE LFMC
END MODULE LFMC_MOD
