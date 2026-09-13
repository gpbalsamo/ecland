! (C) Copyright 2024- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE DFMC_MOD
CONTAINS
SUBROUTINE DFMC(KIDIA,KFDIA,PTSTEP,LDLAND,PT2M,PD2M,PPRE,PTSK,PPSU,PPPR, &
    & PDFMC_1,PDFMC_10,PDFMC_100,PDFMC_1000)

!     PURPOSE
!     -------
!     THIS ROUTINE CALCULATES DEAD FUEL MOISTURE CONTENT FOR BOTH FOLIAGE AND WOOD.
!     PROGNOSTIC: PDFMC_* AND PPPR PERSIST ACROSS TIMESTEPS (DFMC_1/10/100/1000,
!     PPR IN THE PROGNOSTIC STATE, GATED BY LEFIRE / NAMPARVEG).

!     INTERFACE.
!     ----------
!     CALLED FROM *callpar1s*

!     METHOD.
!     DFMC is calculated based on vegetation type following parameterisation from McNorton et al., 2024

!     REFERENCE.
!     Original    J. McNorton      MAR 2024 (arpifs/sparky/dfmc_mod.F90)
!     Ported to ecLand offline     2026 -- dropped the sparky LSP/CP-vs-previous-
!     accumulated-total differencing (PSIPR): ecLand's own Rainf is already a
!     per-timestep rate, unlike IFS's own accumulated-since-forecast-start
!     LSP/CP diagnostics, so no previous-total bookkeeping is needed here.

!==============================================================================

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
REAL(KIND=JPRB),    INTENT(IN)    :: PTSTEP
LOGICAL,            INTENT(IN)    :: LDLAND(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PT2M(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PD2M(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PPRE(:)            ! RAINFALL RATE (kg m-2 s-1)
REAL(KIND=JPRB),    INTENT(IN)    :: PTSK(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PPSU(:)
REAL(KIND=JPRB),    INTENT(INOUT) :: PPPR(:)            ! PREVIOUS-RAINFALL FLAG
REAL(KIND=JPRB),    INTENT(INOUT) :: PDFMC_1(:)         ! DEAD FUEL MOISTURE CONTENT FOLIAGE (1H)
REAL(KIND=JPRB),    INTENT(INOUT) :: PDFMC_10(:)        ! DEAD FUEL MOISTURE CONTENT FOLIAGE (10H)
REAL(KIND=JPRB),    INTENT(INOUT) :: PDFMC_100(:)       ! DEAD FUEL MOISTURE CONTENT WOOD (100H)
REAL(KIND=JPRB),    INTENT(INOUT) :: PDFMC_1000(:)      ! DEAD FUEL MOISTURE CONTENT WOOD (1000H)

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZPHA, ZPPRE, ZRH
REAL(KIND=JPRB) :: ZTa, ZTd, ZTs, ZPa, ZPd, ZPs
REAL(KIND=JPRB) :: ZA, ZHC, ZPR, ZSC, ZTOP, ZRHO, ZCA, ZA1, ZA10, ZA100, ZA1000, ZB1, ZB10, ZB100, ZB1000
REAL(KIND=JPRB) :: ZD, ZD1, ZD10, ZD100, ZD1000
REAL(KIND=JPRB) :: ZE1,ZE10,ZE100,ZE1000

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DFMC_MOD:DFMC',0,ZHOOK_HANDLE)

! CONSTANTS
ZA     = 0.000772_JPRB
ZHC    = 1.59E4_JPRB / 3600.0_JPRB ! convective heat transfer coefficient (J m-2 s-1 K-1)
ZPR    = 0.7_JPRB ! Prandtl number
ZSC    = 0.58_JPRB ! Schmidt number
ZTOP   = 0.662_JPRB * ZHC * (ZPR / ZSC) ** 0.667_JPRB
ZRHO   = 400.0_JPRB ! stick mass density (kg m-3) THIS VARIABLE IS FUEL-TYPE DEPENDENT
ZCA    = 1005.0_JPRB ! constant-pressure specific heat of air (J kg-1 K-1)
ZA1    = 0.0020_JPRB ! radius of stick (m) DRYING TIME (CARLSON 2007)
ZA10   = 0.0064_JPRB ! radius of stick (m) DRYING TIME (CARLSON 2007)
ZA100  = 0.02_JPRB   ! radius of stick (m) DRYING TIME (CARLSON 2007)
ZA1000 = 0.064_JPRB  ! radius of stick (m) DRYING TIME (CARLSON 2007)

ZB1    = ZRHO * ZCA * ZA1
ZB10   = ZRHO * ZCA * ZA10
ZB100  = ZRHO * ZCA * ZA100
ZB1000 = ZRHO * ZCA * ZA1000

DO JL=KIDIA,KFDIA
 IF (.NOT.LDLAND(JL)) THEN
  PDFMC_1(JL)   = 0.0_JPRB
  PDFMC_10(JL)  = 0.0_JPRB
  PDFMC_100(JL) = 0.0_JPRB
  PDFMC_1000(JL)= 0.0_JPRB
 ELSE

  ZPHA = PPSU(JL)
  ZPPRE = PPRE(JL)*3600.0_JPRB ! mm per hour (PPRE already a per-second rate)
  ZTa  = PT2M(JL)-273.15_JPRB
  ZTd  = PD2M(JL)-273.15_JPRB
  ZTs  = PTSK(JL)-273.15_JPRB
  ZRH = (EXP(17.625_JPRB * ZTd/(243.04_JPRB + ZTd)) / EXP(17.625_JPRB * ZTa/(243.04_JPRB + ZTa)))

  ZPa = (0.61094_JPRB*EXP((17.625_JPRB*ZTa)/(ZTa+243.04_JPRB)))/1000.0_JPRB ! August-Roche Magnus (kPa/1000 = Pa)
  ZPd = (0.61094_JPRB*EXP((17.625_JPRB*ZTd)/(ZTd+243.04_JPRB)))/1000.0_JPRB ! August-Roche Magnus (kPa/1000 = Pa)

  ZPs = ZPa*ZRH + ZPHA*ZA*(ZTa-ZTs)

  ZD = 0.32_JPRB * (1.0_JPRB - EXP(-1.885_JPRB * ZPPRE))
  ! Rainfall factors based on Carlson 2007 -- smaller moisture changes when
  ! rain follows rain, tracked via the previous-rainfall flag PPPR.
  IF (PPPR(JL) /= 0.0_JPRB) THEN
   ZD1    = MIN(0.3_JPRB, MAX(0.0_JPRB, 5.0_JPRB*ZD))
   ZD10   = MIN(0.3_JPRB, MAX(0.0_JPRB, 0.15_JPRB*ZD))
   ZD100  = MIN(0.2_JPRB, MAX(0.0_JPRB, 0.25_JPRB*ZD))
   ZD1000 = MIN(0.2_JPRB, MAX(0.0_JPRB, 0.25_JPRB*ZD))
  ELSE
   ZD1    = MIN(0.3_JPRB, MAX(0.0_JPRB, 5.0_JPRB*ZD))
   ZD10   = MIN(0.3_JPRB, MAX(0.0_JPRB, 0.55_JPRB*ZD))
   ZD100  = MIN(0.2_JPRB, MAX(0.0_JPRB, 0.5_JPRB*ZD))
   ZD1000 = MIN(0.2_JPRB, MAX(0.0_JPRB, 0.5_JPRB*ZD))
  ENDIF

  IF (ZPPRE /= 0.0_JPRB) THEN
   PPPR(JL)=1.0_JPRB
  ELSE
   PPPR(JL)=0.0_JPRB
  ENDIF

  IF (ZPHA /= 0.0_JPRB) THEN
    ZE1 = PTSTEP * (ZTOP / (ZPHA * ZB1)) * (ZPs - ZPd)
    ZE10 = PTSTEP * (ZTOP / (ZPHA * ZB10)) * (ZPs - ZPd)
    ZE100 = PTSTEP * (ZTOP / (ZPHA * ZB100)) * (ZPs - ZPd)
    ZE1000 = PTSTEP * (ZTOP / (ZPHA * ZB1000)) * (ZPs - ZPd)
  ELSE
    ZE1 = 0.0_JPRB
    ZE10 = 0.0_JPRB
    ZE100 = 0.0_JPRB
    ZE1000 = 0.0_JPRB
  ENDIF

  PDFMC_1(JL)   = MAX(1.E-7_JPRB,MIN(0.30_JPRB,PDFMC_1(JL)    +ZD1    +ZE1))
  PDFMC_10(JL)  = MAX(1.E-7_JPRB,MIN(0.30_JPRB,PDFMC_10(JL)   +ZD10   +ZE10))
  PDFMC_100(JL) = MAX(1.E-7_JPRB,MIN(0.20_JPRB,PDFMC_100(JL)  +ZD100  +ZE100))
  PDFMC_1000(JL)= MAX(1.E-7_JPRB,MIN(0.20_JPRB,PDFMC_1000(JL) +ZD1000 +ZE1000))

 ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('DFMC_MOD:DFMC',1,ZHOOK_HANDLE)

END SUBROUTINE DFMC
END MODULE DFMC_MOD
