! (C) Copyright 2024- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE FUEL_MOD
CONTAINS
SUBROUTINE FUEL(KIDIA,KFDIA,PTSTEP,KTVL,KTVH,PLAIL,PLAIH,PCVL,PCVH,PLLFL,PLWFL,PDFFL,PDWFL,PNEE)

!     PURPOSE
!     -------
!     THIS ROUTINE CALCULATES LIVE/DEAD LEAF/WOOD FUEL LOAD.
!     PROGNOSTIC: PLLFL/PLWFL/PDFFL/PDWFL PERSIST ACROSS TIMESTEPS (LEFIRE /
!     NAMPARVEG).

!     INTERFACE.
!     ----------
!     CALLED FROM *callpar1s*

!     METHOD.
!     1. Fuel total is updated using ecLand's own CTESSEL net CO2 flux
!     2. Then split between low and high vegetation according to PFT mass density and cover
!     3. Then split between leaf/wood live/dead based on literature and LAI

!     REFERENCE.
!     Original    J. McNorton      JUL 2024 (arpifs/sparky/fuel_mod.F90)
!     Ported to ecLand offline     2026 -- PNEE here is ecLand's own PCO2FLUX
!     (srfcotwo_mod.F90), positive = net release to the atmosphere (ALMA/
!     FLUXNET convention, kg CO2 m-2 s-1). sparky's original formula added its
!     PSNEE (a CAMS reanalysis product using the opposite, uptake-positive,
!     convention) directly to the fuel total; ported here as a subtraction to
!     keep ecLand's own sign convention: net release (PNEE>0) depletes the
!     fuel pool, net uptake (PNEE<0) grows it.

!==============================================================================

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
REAL(KIND=JPRB),    INTENT(IN)    :: PTSTEP
INTEGER(KIND=JPIM), INTENT(IN)    :: KTVL(:)
INTEGER(KIND=JPIM), INTENT(IN)    :: KTVH(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PLAIL(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PLAIH(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PCVL(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PCVH(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PNEE(:)            ! ecLand PCO2FLUX, kg CO2 m-2 s-1, +ve = to atmosphere

REAL(KIND=JPRB), INTENT(INOUT) :: PLLFL(:)     ! LIVE LEAF FUEL LOAD (KG M-2)
REAL(KIND=JPRB), INTENT(INOUT) :: PLWFL(:)     ! LIVE WOOD FUEL LOAD (KG M-2)
REAL(KIND=JPRB), INTENT(INOUT) :: PDFFL(:)     ! DEAD FOLIAGE FUEL LOAD (KG M-2)
REAL(KIND=JPRB), INTENT(INOUT) :: PDWFL(:)     ! DEAD WOOD FUEL LOAD (KG M-2)

REAL(KIND=JPRB) :: ZFD(19) ! VEGETATION SPECIFIC MASS DENSITY
REAL(KIND=JPRB) :: ZLMA(19)! VEGETATION SPECIFIC LEAF MASS DENSITY
REAL(KIND=JPRB) :: ZLAC(19)! VEGETATION SPECIFIC ALPHA COEFFICIENT (Harper et al., 2018 - JULES parameterisation)
REAL(KIND=JPRB) :: ZMSC(19)! VEGETATION SPECIFIC STATIC WOOD COEFFICIENT
REAL(KIND=JPRB) :: ZDFF(19)! VEGETATION SPECIFIC FRACTION OF DEAD FOLIAGE (REMAINING IS WOOD)

INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: ITYL, ITYH
REAL(KIND=JPRB)    :: ZTFL,ZLFL, ZHFL
REAL(KIND=JPRB)    :: ZLLL,ZLWL,ZHLL,ZHWL,ZLDF,ZLDW,ZHDF,ZHDW
REAL(KIND=JPRB)    :: ZNC
REAL(KIND=JPRB)    :: ZCVL, ZCVH ! cover of each vegetation type, zero where the type is absent

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FUEL_MOD:FUEL',0,ZHOOK_HANDLE)

! 1. Crops (L); 2. Short Grass (L); 3. Ev. Needleleaf (H); 4. De. Needleleaf (H);
! 5. De. Broadleaf (H); 6. Ev. Broadleaf (H); 7. Tall Grass/Mixed Crop (L);
! 8. Desert (L); 9. Tundra (L); 10. Irrigated Crops (L); 11. Sparse Veg (L);
! 12. Ice Caps (-); 13. Bogs/Marshes (L); 14. Inland Water (-); 15. Ocean (-);
! 16. Ev. Shrubs (L); 17. De. Shrubs (L); 18. Broadleaf Savannah (H);
! 19. Interrupted Forest (H)

ZFD=(/4.42_JPRB,1.55_JPRB,101.03_JPRB,76.98_JPRB,146.38_JPRB, &
    & 300.01_JPRB,7.54_JPRB,0.01_JPRB,9.56_JPRB,2.47_JPRB,1.18_JPRB, &
    & 0.01_JPRB,3.60_JPRB,0.01_JPRB,0.01_JPRB,63.87_JPRB,63.87_JPRB, &
    & 122.08_JPRB,122.08_JPRB/)

ZLMA=(/0.1370_JPRB,0.0495_JPRB,0.2263_JPRB,0.1006_JPRB,0.0823_JPRB, &
     & 0.1039_JPRB,0.1370_JPRB,0.1370_JPRB,0.0495_JPRB,0.1370_JPRB, &
     & 0.1370_JPRB,0.1370_JPRB,0.1370_JPRB,0.1370_JPRB,0.1370_JPRB, &
     & 0.1515_JPRB,0.0709_JPRB,0.0823_JPRB,0.0823_JPRB/)

ZLAC=(/0.005_JPRB,0.005_JPRB,0.65_JPRB,0.80_JPRB,0.78_JPRB,0.845_JPRB, &
     & 0.005_JPRB,0.005_JPRB,0.005_JPRB,0.005_JPRB,0.005_JPRB, &
     & 0.005_JPRB,0.005_JPRB,0.005_JPRB,0.005_JPRB,0.13_JPRB,0.13_JPRB, &
     & 0.78_JPRB,0.78_JPRB/)

ZMSC=(/0.62833981_JPRB, 0.58869273_JPRB, 0.46880262_JPRB, 0.33528824_JPRB, &
    &  0.47032268_JPRB, 0.41006731_JPRB, 0.48550652_JPRB,0.0_JPRB, 0.01_JPRB, &
    &  0.6960738_JPRB, 0.5855936_JPRB, 0.01_JPRB, 0.27677251_JPRB, 0.01_JPRB, &
    & 0.01_JPRB, 0.35474094_JPRB, 0.5985808_JPRB,0.43414621_JPRB, 0.01_JPRB/)

ZDFF=(/1.0_JPRB,1.0_JPRB,0.45_JPRB,0.6_JPRB,0.9_JPRB,0.1_JPRB,1.0_JPRB, &
     & 0.5_JPRB,0.5_JPRB,1.0_JPRB,0.5_JPRB,0.5_JPRB,0.5_JPRB,0.5_JPRB, &
     & 0.5_JPRB,0.5_JPRB,0.9_JPRB,0.9_JPRB,0.9_JPRB/)

DO JL=KIDIA,KFDIA
! KTVL/KTVH==0 means "no vegetation of this type" (bare-soil, or a point with
! only low or only high vegetation). ZFD..ZDFF are 1-indexed per-type tables
! with no slot 0, so indexing them with a zero type is an out-of-bounds read
! (silent in a release build, "Subscript #1 of the array ZFD has value 0"
! under -check bounds). Treat the absent type as zero cover and clamp its
! (then irrelevant) table index to 1: every weight and cap below multiplies
! by cover, so the absent type contributes exactly zero, and a point with
! only one vegetation type gets the whole fuel pool on that type. Points
! with both types present are unchanged.
  ZCVL = 0.0_JPRB
  ZCVH = 0.0_JPRB
  IF (KTVL(JL) > 0) ZCVL = PCVL(JL)
  IF (KTVH(JL) > 0) ZCVH = PCVH(JL)

  IF ((ZCVL + ZCVH) > 0.001_JPRB) THEN

   ITYL = MAX(1,KTVL(JL))
   ITYH = MAX(1,KTVH(JL))
   ZNC  = 1.0_JPRB/(ZCVH+ZCVL)

   ZTFL = PLLFL(JL) + PLWFL(JL) + PDFFL(JL) + PDWFL(JL)
   ZTFL = MAX(0.0_JPRB,ZTFL - (PNEE(JL)*PTSTEP/2.0_JPRB)) ! PNEE in kg CO2 m-2 s-1, assume total flux 50% carbon

   ZLFL = ZTFL * ( (ZFD(ITYL)*ZCVL*ZNC) / (ZFD(ITYL)*ZCVL*ZNC + ZFD(ITYH)*ZCVH*ZNC) )
   ZHFL = ZTFL * ( (ZFD(ITYH)*ZCVH*ZNC) / (ZFD(ITYL)*ZCVL*ZNC + ZFD(ITYH)*ZCVH*ZNC) )

   ZLLL = MIN(ZLFL, ZLMA(ITYL)*ZCVL*ZNC*PLAIL(JL)) ! Live Low Leaf Load
   ZHLL = MIN(ZHFL, ZLMA(ITYH)*ZCVH*ZNC*PLAIH(JL)) ! Live High Leaf Load

   ZLWL = MIN(ZLFL - ZLLL, (ZLAC(ITYL)*ZCVL*ZNC*PLAIL(JL)**1.667_JPRB + ZMSC(ITYL)*ZLFL)) ! Live Low Wood Load
   ZHWL = MIN(ZHFL - ZHLL, (ZLAC(ITYH)*ZCVH*ZNC*PLAIH(JL)**1.667_JPRB + ZMSC(ITYH)*ZHFL)) ! Live High Wood Load

   ZLDF = MAX(0.0_JPRB,ZDFF(ITYL)*(ZLFL-(ZLLL + ZLWL))) ! Dead Low Foliage Load
   ZHDF = MAX(0.0_JPRB,ZDFF(ITYH)*(ZHFL-(ZHLL + ZHWL))) ! Dead High Foliage Load

   ZLDW = MAX(0.0_JPRB,ZLFL - (ZLLL + ZLWL + ZLDF)) ! Dead Low Wood Load
   ZHDW = MAX(0.0_JPRB,ZHFL - (ZHLL + ZHWL + ZHDF)) ! Dead High Wood Load

   PLLFL(JL) = ZLLL+ZHLL
   PLWFL(JL) = ZLWL+ZHWL
   PDFFL(JL) = ZLDF+ZHDF
   PDWFL(JL) = ZLDW+ZHDW

  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('FUEL_MOD:FUEL',1,ZHOOK_HANDLE)
END SUBROUTINE FUEL
END MODULE FUEL_MOD
