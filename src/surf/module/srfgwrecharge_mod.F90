MODULE SRFGWRECHARGE_MOD
CONTAINS
SUBROUTINE SRFGWRECHARGE(KIDIA,KFDIA,KLON,KLEVS,&
 & PTSPHY,&
 & PSSDP3,&
 & YDSOIL,&
 & LDLAND,&
 & PWTDM1M,&
 & PGWDRAINFLUX,&
 & PSDOR,&
 & PRHSW,&
 & PWTDE1,&
 & PGWDRAINUSED,&
 & PGWCAPFLUX)

USE PARKIND1        , ONLY : JPIM, JPRB
USE YOMHOOK         , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_SOIL        , ONLY : TSOIL
USE YOMSURF_SSDP_MOD, ONLY : SSDP3D_ID
USE YOS_THF         , ONLY : RHOH2O

! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFGWRECHARGE* - BIDIRECTIONAL WATER-TABLE/SOIL EXCHANGE: CAPILLARY-
!                        RISE EXTRACTION FROM A FINITE PROGNOSTIC AQUIFER,
!                        PLUS FREE-DRAINAGE RECHARGE BACK INTO IT

!     PURPOSE.
!     --------
!     THIRD OF THE "DEPTH TRILOGY" PROTOTYPES ALONGSIDE LEUNIFORMROOT
!     (ROOTING DEPTH) AND LEBEDROCKLIM (BEDROCK-LIMITED DRAINAGE). THE
!     PROGNOSTIC WATER TABLE DEPTH (*PWTDM1M*/*PWTDE1*) IS A GENUINE,
!     RESTART-CHECKPOINTED ECLAND STATE VARIABLE (SEE YOMGP1S0/1/A'S
!     WTDNU0/1/A AND RDSUPR'S 'WTD' CVARS2D ENTRY), NOT A MODULE-LEVEL
!     SAVE. IT EVOLVES IN BOTH DIRECTIONS:

!     (1) EXTRACTION (DEEPENS THE WATER TABLE): A GARDNER (1958)-TYPE
!         EXPONENTIAL CAPILLARY-RISE FLUX FROM THE WATER TABLE INTO EVERY
!         LAYER ABOVE IT, DECAYING WITH HEIGHT ABOVE THE WATER TABLE AT
!         THE VAN GENUCHTEN ALPHA RATE -- 1/ALPHA IS THE CHARACTERISTIC
!         CAPILLARY LENGTH FOR THAT SOIL TEXTURE. USES ONLY PARAMETERS
!         ECLAND ALREADY COMPUTES (RWCONSM3D, RMVGALPHA3D).

!     (2) RECHARGE (SHALLOWS THE WATER TABLE): A STATE-DEPENDENT FRACTION
!         OF THE ACTUAL GRAVITY-DRAINAGE FLUX LEAVING THE BOTTOM SOIL LAYER,
!         AS DIAGNOSED BY *SRFWEXC_VG* (PSAWGFL/ZSAWGFL AT JK=KLEVS_WB --
!         ALREADY POST-LEBEDROCKLIM TAPER AND POST-LEFRZFLOOR FLOOR),
!         PASSED IN HERE AS *PGWDRAINFLUX*. NOTE *PROFD* -- WHAT THIS FLUX
!         WOULD OTHERWISE ENTIRELY BECOME -- IS QSB IN ECLAND'S OWN
!         ACCOUNTING, AND IS PURELY THIS GRAVITY-DRAINAGE TERM: THE
!         LESSRO/VIC SATURATION-EXCESS MECHANISM FEEDS *PROFS* (QS)
!         SEPARATELY, NOT QSB, SO THE TWO ARE NOT ENTANGLED HERE.
!
!         DIVERTING THE FULL FLUX TO THE AQUIFER (THE FIRST VERSION OF
!         THIS FIX) IS ITS OWN OVER-SIMPLIFICATION: PHYSICALLY, DRAINAGE
!         REACHING A WATER TABLE SPLITS BETWEEN RECHARGING LOCAL STORAGE
!         AND CONTINUING DOWN-GRADIENT AS LATERAL BASEFLOW, NOT ONE OR THE
!         OTHER. A FLAT, TERRAIN-ONLY PARTITION (THE SECOND VERSION) WAS
!         ALSO INCOMPLETE: A STATIC FRACTION HAS NO WAY TO STOP RECHARGING
!         AN AQUIFER THAT IS ALREADY NEAR ITS OWN "FULL" STATE, SO IT EITHER
!         COLLAPSES SITES WHERE 100% DIVERSION WOULD HAVE HELD A DEEPER,
!         STABLE WATER TABLE, OR IT STAYS TOO GENEROUS TO RECOVER MUCH
!         LATERAL FLOW AT ALL -- IT CANNOT DO BOTH AT ONCE. 42-SITE PLUMBER2
!         VALIDATION OF THAT VERSION CONFIRMED EXACTLY THIS: IT RECOVERED
!         ONLY ~28% OF THE OFF-VS-100%-DIVERSION QSB GAP, AND FLIPPED SEVERAL
!         PREVIOUSLY-HELD SITES (DE-THA, US-HA1, FI-HYY, ...) TO COLLAPSE.
!
!         THIS VERSION MAKES THE PARTITION STATE-DEPENDENT INSTEAD, REUSING
!         THE VIC (VARIABLE INFILTRATION CAPACITY) ALGEBRA *SRFWEXC_VG*'S
!         OWN LESSRO BLOCK USES FOR SURFACE SATURATION-EXCESS RUNOFF --
!         SAME SHAPE, APPLIED ONE LEVEL DOWN. LESSRO ASKS "GIVEN HOW WET THE
!         TOP *RSRDEP* METRES OF SOIL ALREADY ARE, HOW MUCH OF THIS RAINFALL
!         CAN STILL INFILTRATE VS. RUN OFF?"; HERE THE QUESTION IS "GIVEN
!         HOW FULL THE AQUIFER ALREADY IS RELATIVE TO ITS OWN EFFECTIVE
!         DEPTH *RGWLATDEPTH*, HOW MUCH OF THIS DRAINAGE CAN STILL RECHARGE
!         IT VS. CONTINUE AS LATERAL BASEFLOW?". CURRENT AQUIFER STORAGE IS
!         (RGWLATDEPTH-WTD)*RGWSPECYIELD (ZERO ONCE THE WATER TABLE IS AT
!         OR BELOW RGWLATDEPTH, MAXIMAL AS IT APPROACHES THE SURFACE) --
!         THE DIRECT ANALOGUE OF LESSRO'S SOIL-MOISTURE STATE, PLAYING THE
!         SAME ROLE THE USER'S OWN QUESTION IDENTIFIED ("THE SOIL MOISTURE
!         BEING ALSO PART OF THE EQUATION"). THE VIC SHAPE PARAMETER (B) IS
!         STILL SET BY *PSDOR* VIA THE SAME ZROEFF FORM LESSRO USES, BUT
!         THROUGH THE DEDICATED *RGWLATSIGMIN*/*RGWLATSIGMAX* RANGE BELOW
!         RATHER THAN LESSRO'S OWN *RSIGORMIN*/*RSIGORMAX* (THOSE, REUSED
!         AS-IS, WERE ALREADY SHOWN NOT TO DIFFERENTIATE REAL PLUMBER2
!         TERRAIN FOR THIS PURPOSE). *RGWLATDEPTH* IS THE NEW FREE
!         PARAMETER THIS ADDS, ANALOGOUS TO LESSRO'S OWN *RSRDEP* -- IT SETS
!         HOW QUICKLY (IN TERMS OF WATER-TABLE RISE) THE AQUIFER'S
!         "SATURATED FRACTION" APPROACHES 1 AND STARTS REJECTING RECHARGE
!         TO LATERAL FLOW. ALL THREE (RGWLATSIGMIN, RGWLATSIGMAX,
!         RGWLATDEPTH) ARE DELIBERATELY LEFT AS HARDCODED PARAMETERS FOR
!         MANUAL CALIBRATION, THE SAME WAY THE EARLIER TWO WERE.
!
!         THE RESULT IS SELF-LIMITING IN BOTH DIRECTIONS WITHOUT NEEDING THE
!         SEPARATE HARD TENDENCY CLAMP TO DO ALL THE WORK: AN EMPTY AQUIFER
!         (DEEP WTD) TAKES NEARLY ALL OF THE INCOMING FLUX AS RECHARGE, A
!         NEAR-FULL ONE (SHALLOW WTD) PASSES NEARLY ALL OF IT ON AS LATERAL
!         BASEFLOW, AND THE VIC ALGEBRA INTERPOLATES SMOOTHLY (AND, VIA THE
!         B-EXPONENT, HETEROGENEOUSLY BY TERRAIN) IN BETWEEN.
!
!         *PGWDRAINUSED* REPORTS THE PORTION ACTUALLY USED FOR RECHARGE SO
!         THE CALLER CAN SUBTRACT IT FROM WHAT REACHES *PROFD*/QSB --
!         OTHERWISE THE SAME WATER WOULD BE COUNTED TWICE.

!     BOTH TERMS CONVERT A WATER-EQUIVALENT DEPTH (M) TO A WATER-TABLE-
!     DEPTH CHANGE VIA THE SPECIFIC YIELD (DRAINABLE POROSITY, RGWSPECYIELD
!     -- SIZE OF THE UNCONFINED AQUIFER), THE SAME STORAGE-TO-DEPTH
!     RELATIONSHIP USED BY NIU ET AL. (2007)'S SIMPLE GROUNDWATER SCHEME
!     FOR CLM. THE RESULT IS EXPRESSED AS A TENDENCY (*PWTDE1*, M/S), NOT
!     AN IN-PLACE UPDATE: ECLAND INTEGRATES GP1=GP0+TDT*ZGPE_STACK FOR
!     EVERY PROGNOSTIC FIELD GENERICALLY (SEE CPG1S.F90), SO THE WATER
!     TABLE DEPTH IS CLAMPED HERE BY LIMITING THE TENDENCY ITSELF (NEVER
!     LET IT DEEPEN PAST RDBEDROCK -- SOLID ROCK -- OR SHALLOW PAST
!     RGWTD_MIN), NOT BY REASSIGNING STATE.

!     IT IS PURELY ADDITIVE, LIKE SRFWVAPOR: THE LIQUID DIFFUSIVITY
!     MATRIX (*PCFW*) BUILT BY *SRFWEXC_VG* (OR *SRFWEXC*) IS NOT
!     TOUCHED, AND *SRFWDIF*/*SRFWINC* DOWNSTREAM ARE UNCHANGED.

!     THE "NO DATA AVAILABLE" INITIAL/DEFAULT WATER TABLE DEPTH (100M, A
!     SINGLE GLOBAL VALUE FOR NOW, NOT YET PER-GRIDPOINT) LIVES WHERE THE
!     FIELD ITSELF IS INITIALISED (SUGP1S.F90 COLD START, RDSUPR.F90'S
!     "MISSING FROM RESTART FILE" FALLBACK), NOT IN YOS_SOIL.F90 -- SEE
!     LEGWRECHARGE'S COMMENT THERE. RDBEDROCK CAPS HOW DEEP THE WATER TABLE
!     CAN GO -- THE AQUIFER CANNOT DRAIN BELOW SOLID ROCK.

!**   INTERFACE.
!     ----------
!          *SRFGWRECHARGE* IS CALLED FROM *SURFTSTP_CTL*, ALONGSIDE
!     *SRFWVAPOR* -- AFTER *SRFWEXC_VG* or *SRFWEXC* (SO PRHSW ALREADY
!     HOLDS THE LIQUID BUDGET), BUT BEFORE *SRFWDIF*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF ACTIVE SOIL-WATER LAYERS (KLEVS_WB AT THE
!                 CALL SITE, MATCHING THE SUBSEQUENT *SRFWDIF* CALL)

!     INPUT PARAMETERS (REAL):
!    *PTSPHY*     TIME STEP                                       S
!    *PSSDP3*     SPATIALLY DISTRIBUTED SOIL PARAMETERS (VAN
!                 GENUCHTEN ALPHA, SATURATED/RESIDUAL MOISTURE,
!                 REFERENCE CONDUCTIVITY, ...)
!    *PWTDM1M*    WATER-TABLE DEPTH AT THE START OF THIS TIMESTEP    M
!    *PGWDRAINFLUX* GRAVITY-DRAINAGE FLUX OUT OF THE BOTTOM LAYER,
!                 AS DIAGNOSED BY SRFWEXC_VG (PSAWGFL AT KLEVS_WB,
!                 POST-LEBEDROCKLIM TAPER)                    KG/M**2/S
!    *PSDOR*      SUBGRID STANDARD DEVIATION OF OROGRAPHY -- SETS THE VIC
!                 SHAPE PARAMETER FOR THE LATERAL/VERTICAL PARTITION OF
!                 PGWDRAINFLUX, SAME PARAMETER SRFWEXC_VG'S OWN LESSRO/VIC
!                 BLOCK USES (BUT SCALED BY THE DEDICATED RGWLATSIGMIN/
!                 RGWLATSIGMAX BELOW, NOT LESSRO'S OWN RSIGORMIN/MAX)      M

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS (DERIVED TYPE):
!    *YDSOIL*     SOIL PARAMETERS (RDAW LAYER THICKNESSES, RDBEDROCK LOWER
!                 BOUND, RGWSPECYIELD AQUIFER SPECIFIC YIELD)

!     UPDATED PARAMETERS (REAL):
!    *PRHSW*      RIGHT-HAND SIDE OF SOIL MOISTURE EQUATIONS,
!                 AS BUILT BY *SRFWEXC_VG* or *SRFWEXC*; THE EXTRACTION
!                 AND RECHARGE FLUXES ARE ADDED TO IT IN PLACE. M**3/M**3
!    *PWTDE1*     WATER-TABLE DEPTH TENDENCY, ACCUMULATED IN PLACE   M/S

!     OUTPUT PARAMETERS (REAL):
!    *PGWCAPFLUX*  CAPILLARY-RISE FLUX DRAWN FROM THE AQUIFER INTO THE
!                 SOIL COLUMN THIS STEP (>=0)                     KG/M**2/S
!    *PGWDRAINUSED* PORTION OF *PGWDRAINFLUX* DIVERTED TO THE AQUIFER --
!                 SUBTRACT FROM THE SAME FLUX BEFORE IT REACHES PROFD/QSB
!                 TO AVOID DOUBLE-COUNTING                    KG/M**2/S

!     REFERENCE.
!     ----------
!     Gardner, W.R. (1958): Some steady-state solutions of the unsaturated
!     moisture flow equation with application to evaporation from a water
!     table. Soil Science, 85(4), 228-232.
!     Niu, G.-Y., Yang, Z.-L., Dickinson, R.E., Gulden, L.E., Su, H.
!     (2007): Development of a simple groundwater model for use in
!     climate models and evaluation with Gravity Recovery and Climate
!     Experiment data. JGR Atmospheres, 112, D07103.

!     Original
!     Prototype for the ecLand "depth trilogy" (rooting depth / bedrock /
!     water table), third of three. 2026-08-31
!     Made genuinely prognostic and bidirectional. 2026-08-31
!     Recharge now uses SRFWEXC_VG's own diagnosed drainage flux instead
!     of a separately-calibrated saturation threshold. 2026-09-01
!     Recharge now takes only an orography-set fraction of that flux,
!     leaving the rest as lateral baseflow (PROFD/Qsb) -- diverting all
!     of it was itself an over-simplification. 2026-09-02
!     Static orography fraction replaced by a state-dependent VIC-style
!     partition on aquifer fullness (new RGWLATDEPTH parameter), since the
!     static version could not both recover lateral flow and preserve
!     water-table stability at the same sites. 2026-09-02
!     RGWLATSIGMIN/RGWLATSIGMAX corrected from an ad hoc 10/150 (tuned only
!     against this branch's own 42-site PLUMBER2 PSDOR distribution) to
!     57/525, matching RSIGORMIN/RSIGORMAX now used for LESSRO itself at
!     this branch's actual TCo1279 (~9km) clim resolution (see
!     namelist_ecland_50R1_runoff_fix in plumber2-ecland) -- the old 100/
!     1000 Fortran default several sites had silently inherited (via the
!     _ctl namelist, which never set these) was for a much coarser
!     resolution than these clim files are. RGWLATDEPTH swept 5-80m on a
!     5-site subset at the corrected sigma values: unlike the static
!     partition, larger RGWLATDEPTH improved BOTH water-table realism and
!     Qsb preservation together (no forced tradeoff), still climbing at
!     80m, not yet saturated. 2026-09-02
!
!     KNOWN OPEN ISSUE, PARKED 2026-09-02: at that point RGWLATDEPTH=80 did
!     NOT hold across soil vertical discretizations. Same 5-site subset,
!     NCSS=4 vs 9 (identical 2.89m total depth, finer layers only) vs 14
!     (12m total depth): AU-Tum's final WTD came out 6.41m / 25.01m /
!     31.77m respectively -- NCSS=4->9 alone (pure refinement, same
!     physical depth) already differs by ~4x, so this is not simply
!     "deeper profile reaches a different equilibrium" but a genuine
!     discretization sensitivity. Suspected cause: PGWDRAINFLUX itself
!     (SRFWEXC_VG's diagnosed bottom-layer gravity-drainage flux, the
!     input to this whole partition) depends on the bottom layer's own
!     thickness/depth, which differs substantially across NCSS=4/9/14 --
!     so RGWLATDEPTH calibrated against one discretization's typical flux
!     magnitude does not transfer to another's. Not yet root-caused or
!     fixed; treat every calibration below as NCSS=4-specific until this
!     is resolved. (Confirmed independently, same date: NCSS=4->9 alone
!     also shifts total column water -5.5% network-mean, total ET -6.9%
!     (ESoil -12.0%, TVeg -5.7%), and Qsb +15% -- all from the near-
!     surface/ET side of the discretization change, NOT from RGWLATDEPTH:
!     identical at RGWLATDEPTH=80 and =300, so that part of the
!     sensitivity is baked into the core model, independent of this
!     partition's own calibration.)
!
!     42-SITE CALIBRATION AGAINST A REVISED TARGET, 2026-09-04: requirement
!     changed to "hold WTD within 20% of Fan et al. (2017) or deeper" (NOT
!     shallower), scored only against sites with Fan WTD >=10m (shallower
!     sites are typically lateral-inflow-maintained -- e.g. riparian/
!     valley-bottom settings -- which a 1D column cannot represent even in
!     principle) and excluding FR-Pue/CN-Din (Fan WTD 197m/257m, beyond
!     RDBEDROCK=100m -- structurally unreachable regardless of this
!     partition's tuning). 29 evaluable sites. Swept RGWLATSIGMIN/MAX and
!     RGWLATDEPTH across four regimes spanning the plausible range and
!     beyond -- (57,525,80), (57,525,300), (5,150,100), and an extreme
!     (0,10,100) pinning the VIC exponent near its 0.5 ceiling for nearly
!     every site -- and ALL FOUR landed on exactly 16/29 passing. The 13
!     that fail in every regime (DE-Gri, US-UMB, US-PFa, US-Blo, IT-Lav,
!     FI-Hyy, DE-Hai, US-MMS, AU-Ctr, AU-Tum, US-Me2, US-Ha1, DE-Tha) move
!     a little in the right direction under more aggressive settings (e.g.
!     US-UMB 4.35m->7.06m, DE-Tha 28.7m->33.0m going from (5,150,100) to
!     the extreme) but never cross their threshold. This is a genuine
!     ceiling, not an undertuned parameter: the VIC shape exponent only
!     materially changes the recharge/lateral split very close to full
!     saturation (the storage-deficit term dominates everywhere else,
!     regardless of exponent), and RGWLATDEPTH's own effect saturated
!     between 80 and 300 -- so these three parameters, however combined,
!     cap out at ~55% of evaluable sites for this target. Left at
!     (57,525,100) -- physically the most defensible of the four tested
!     (57/525 matches LESSRO's own resolution-correct calibration; 100
!     coincides with RDBEDROCK) -- since the more extreme settings bought
!     no additional passes. Reaching the remaining 13 sites, if possible
!     at all in a 1D column, likely needs a different lever (RGWSPECYIELD,
!     the extraction side of the balance, or accepting them as a
!     documented limitation alongside the shallow-WTD exclusion above).
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS

REAL(KIND=JPRB),    INTENT(IN)   :: PTSPHY
REAL(KIND=JPRB),    INTENT(IN)   :: PSSDP3(:,:,:)
TYPE(TSOIL),         INTENT(IN)  :: YDSOIL
LOGICAL,             INTENT(IN)  :: LDLAND(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWTDM1M(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PGWDRAINFLUX(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSDOR(:)

REAL(KIND=JPRB),    INTENT(INOUT):: PRHSW(:,:)
REAL(KIND=JPRB),    INTENT(INOUT):: PWTDE1(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PGWDRAINUSED(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PGWCAPFLUX(:)

!*         0.2    DECLARATION OF LOCAL VARIABLES.
!                 ----------- -- ----- ----------

! Floor on water-table depth: keeps the aquifer strictly below the active
! soil-water profile and away from a nonphysical depth-zero singularity in
! the tendency-clamp arithmetic below.
REAL(KIND=JPRB), PARAMETER :: RGWTD_MIN=0.5_JPRB

! Separate PARAMETERs from YDSOIL%RSIGORMIN/RSIGORMAX (not tied to the
! namelist) so this scheme's own calibration can move independently of
! LESSRO's, even though they start from the same numbers: PSDOR's own
! statistics are resolution-dependent, and RSIGORMIN=100/RSIGORMAX=1000
! (the Fortran default) is the WRONG scale for the TCo1279 (~9km) clim
! files this branch's PLUMBER2 sites are built from -- the resolution-
! correct operational values are RSIGORMIN=57/RSIGORMAX=525 (see
! namelist_ecland_50R1_runoff_fix). These two start at that same 57/525,
! since PSDOR's heterogeneity meaning at this resolution is exactly
! LESSRO's own calibration question -- free parameters, recalibrate here
! if the two schemes' required sensitivity to PSDOR turns out to differ.
REAL(KIND=JPRB), PARAMETER :: RGWLATSIGMIN=57.0_JPRB
REAL(KIND=JPRB), PARAMETER :: RGWLATSIGMAX=525.0_JPRB

! Effective depth (M) of the "fast" aquifer zone the VIC saturation-excess
! algebra below operates over -- the direct analogue of LESSRO's own
! RSRDEP (0.5m by default) for the top soil layer, just at the much larger
! scale appropriate to a water table that moves over months to years
! rather than a storm event. Aquifer storage is (RGWLATDEPTH-WTD)*
! RGWSPECYIELD: zero once WTD reaches or exceeds RGWLATDEPTH, maximal as
! WTD approaches the surface. Placeholder value below -- free parameter,
! calibrate alongside RGWLATSIGMIN/RGWLATSIGMAX.
REAL(KIND=JPRB), PARAMETER :: RGWLATDEPTH=100.0_JPRB

REAL(KIND=JPRB) :: ZDEPTH_UPPER, ZDEPTH_MID, ZDIST, ZCAPFLUX, ZEXTRACT
REAL(KIND=JPRB) :: ZRECHARGE, ZDWTDDT, ZROEFF
REAL(KIND=JPRB) :: ZBWS, ZB1, ZBM, ZWSTATE, ZWMAX, ZCONW1
REAL(KIND=JPRB) :: ZLYEPS, ZLIMRS, ZLYSIC, ZVOL, ZROS, ZDRAINDEPTH

INTEGER(KIND=JPIM) :: JK, JL, JKBOT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRFGWRECHARGE_MOD:SRFGWRECHARGE',0,ZHOOK_HANDLE)

JKBOT=KLEVS

ASSOCIATE(RDAW=>YDSOIL%RDAW, &
 & RWCONSM3D=>PSSDP3(:,:,SSDP3D_ID%NRWCONSM3D), &
 & RMVGALPHA3D=>PSSDP3(:,:,SSDP3D_ID%NRMVGALPHA3D) )

DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN

!*         1.    CAPILLARY-RISE EXTRACTION PER LAYER (ZERO FOR ANY LAYER
!                AT OR BELOW THIS POINT'S CURRENT WATER TABLE DEPTH).

    ZEXTRACT = 0.0_JPRB
    ZDEPTH_UPPER = 0.0_JPRB
    DO JK=1,KLEVS
      ZDEPTH_MID = ZDEPTH_UPPER + 0.5_JPRB*RDAW(JK)
      ZDIST = PWTDM1M(JL) - ZDEPTH_MID
      IF (ZDIST > 0.0_JPRB) THEN
        ZCAPFLUX = RWCONSM3D(JL,JK)*EXP(-RMVGALPHA3D(JL,JK)*ZDIST)
        PRHSW(JL,JK) = PRHSW(JL,JK) + PTSPHY*ZCAPFLUX/RDAW(JK)
        ZEXTRACT = ZEXTRACT + PTSPHY*ZCAPFLUX ! water-equivalent depth (m)
      ENDIF
      ZDEPTH_UPPER = ZDEPTH_UPPER + RDAW(JK)
    ENDDO

!*         2.    RECHARGE = A STATE-DEPENDENT (VIC-STYLE) FRACTION OF THE
!                ACTUAL GRAVITY-DRAINAGE FLUX SRFWEXC_VG DIAGNOSED LEAVING
!                THE BOTTOM LAYER -- THE REST STAYS LATERAL (PROFD/QSB).
!                SAME ALGEBRA AS SRFWEXC_VG'S OWN LESSRO BLOCK, APPLIED TO
!                AQUIFER FULLNESS (OVER RGWLATDEPTH) INSTEAD OF TOP-LAYER
!                SOIL MOISTURE (OVER RSRDEP): AN EMPTY AQUIFER ABSORBS
!                NEARLY ALL OF PGWDRAINFLUX AS RECHARGE; A FULL ONE PASSES
!                NEARLY ALL OF IT ON AS LATERAL FLOW; THE VIC SHAPE
!                PARAMETER (SET BY PSDOR VIA RGWLATSIGMIN/RGWLATSIGMAX)
!                CONTROLS HOW SHARP THAT TRANSITION IS, EXACTLY AS IT
!                CONTROLS LESSRO'S OWN SATURATION-EXCESS TRANSITION.

    ZDRAINDEPTH = PGWDRAINFLUX(JL)*PTSPHY/RHOH2O ! kg/m2/s -> m (this step)

    IF (ZDRAINDEPTH > 0.0_JPRB) THEN

      ZWSTATE = MAX(0.0_JPRB,RGWLATDEPTH-PWTDM1M(JL))*YDSOIL%RGWSPECYIELD
      ZWMAX   = RGWLATDEPTH*YDSOIL%RGWSPECYIELD

      ZROEFF = MAX(0.0_JPRB,(PSDOR(JL)-RGWLATSIGMIN))/(PSDOR(JL)+RGWLATSIGMAX)
      ZBWS = MAX(MIN(ZROEFF,0.5_JPRB),0.01_JPRB)
      ZB1 = 1.0_JPRB+ZBWS
      ZBM = 1.0_JPRB/ZB1

      ZCONW1 = ZWMAX*ZB1
      ZLYEPS = MAX(0.0_JPRB,ZWSTATE-ZWMAX)
      ZLIMRS = -1.0_JPRB*ZLYEPS
      IF (ZLYEPS > 0.1_JPRB*ZWMAX) THEN
        ZLIMRS = 0.0_JPRB
      ENDIF

!            VIC saturated-fraction of the aquifer's effective depth
      ZLYSIC = MIN(1.0_JPRB,MAX(0.0_JPRB,(ZWSTATE-ZLYEPS)/ZWMAX))
      ZVOL = MAX(0.0_JPRB,(1.0_JPRB-ZLYSIC)**ZBM-ZDRAINDEPTH/ZCONW1)

!            Portion of the incoming drainage that bypasses recharge and
!            stays lateral flow (PROFD/Qsb), by the same VIC algebra
!            LESSRO uses for saturation-excess surface runoff
      ZROS = ZDRAINDEPTH-MAX(ZWMAX-ZWSTATE,ZLIMRS)
      IF (ZVOL > 0.0_JPRB) THEN
        ZROS = ZROS+ZWMAX*ZVOL**ZB1
      ENDIF
      ZROS = MAX(ZROS,0.0_JPRB)

!            RWTDRECHARGE throttles what the VIC partition would otherwise
!            hand the aquifer. The partition alone reduces to
!            (WTD/RGWLATDEPTH)**(ZBWS/(1+ZBWS)) and ZBWS is clamped to
!            [0.01,0.5], so the exponent never exceeds 1/3 and the recharge
!            fraction stays between 0.65 and 0.99 across the whole legal range
!            of the shape parameter at a 27 m water table -- which is why
!            sweeping RGWLATSIGMIN/RGWLATSIGMAX/RGWLATDEPTH could not move the
!            42-site score off 16/29 (179176e). Whatever is withheld here stays
!            in the bottom-layer drainage and reaches PROFD/Qsb through
!            PGWDRAINUSED, so the balance closes with no extra bookkeeping.
!            Defaults to 1.0 = the pre-existing behaviour, bit for bit.
      ZRECHARGE = YDSOIL%RWTDRECHARGE*MAX(0.0_JPRB,ZDRAINDEPTH-ZROS) ! water-equivalent depth (m)
    ELSE
      ZRECHARGE = 0.0_JPRB
    ENDIF

!            Recharge is a DIVERSION of the bottom-layer drainage, not an extra
!            withdrawal from the soil. That water has already left the column
!            through the free-drainage bottom boundary -- it is what SRFWINC
!            diagnoses as PSAWGFL and what PGWDRAINFLUX above is built from --
!            so debiting PRHSW for it here as well removed it twice while it
!            arrived once. The budget
!              Rainf+Snowf+Evap+Qs+Qsb-DelIntercept-DelSoilMoist-DelSWE
!            then came out at exactly 2x the recharge instead of 1x (measured
!            2.000/2.000/2.001/2.002 over a tenfold RWTDRECHARGE range), against
!            0.016 mm/yr with LEGWRECHARGE off. Reporting the diversion through
!            PGWDRAINUSED, which SURFTSTP_CTL subtracts from PROFD, is the whole
!            of the bookkeeping needed.
    PGWDRAINUSED(JL) = ZRECHARGE*RHOH2O/PTSPHY ! m -> kg/m2/s
    PGWCAPFLUX(JL)   = ZEXTRACT*RHOH2O/PTSPHY  ! m -> kg/m2/s

!*         3.    NET WATER-TABLE-DEPTH TENDENCY, CLAMPED SO THE IMPLIED
!                END-OF-STEP DEPTH STAYS WITHIN [RGWTD_MIN, RDBEDROCK].

    ZDWTDDT = (ZEXTRACT-ZRECHARGE)/(YDSOIL%RGWSPECYIELD*PTSPHY)
    ZDWTDDT = MAX(ZDWTDDT,(RGWTD_MIN-PWTDM1M(JL))/PTSPHY)
    ZDWTDDT = MIN(ZDWTDDT,(YDSOIL%RDBEDROCK-PWTDM1M(JL))/PTSPHY)

    PWTDE1(JL) = PWTDE1(JL) + ZDWTDDT

  ELSE
    PGWDRAINUSED(JL) = 0.0_JPRB
    PGWCAPFLUX(JL)   = 0.0_JPRB
  ENDIF
ENDDO

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SRFGWRECHARGE_MOD:SRFGWRECHARGE',1,ZHOOK_HANDLE)

END SUBROUTINE SRFGWRECHARGE
END MODULE SRFGWRECHARGE_MOD
