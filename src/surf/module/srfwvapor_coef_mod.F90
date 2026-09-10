MODULE SRFWVAPOR_COEF_MOD
CONTAINS
SUBROUTINE SRFWVAPOR_COEF(PTHETA_A,PTHETA_B,PTEMP_A,PTEMP_B,PDZ,&
 & PWSAT_A,PWREST_A,PALPHA_A,PN_A,PM_A,&
 & PWSAT_B,PWREST_B,PALPHA_B,PN_B,PM_B,&
 & YDCST,&
 & PQVAP)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOS_CST   , ONLY : TCST
USE YOS_THF   , ONLY : R2ES, R3LES, R4LES, R3IES, R4IES, R5LES, R5IES, RHOH2O

! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFWVAPOR_COEF* -  POINT-WISE PHILIP & DE VRIES (1957) SOIL
!                         WATER-VAPOUR FLUX BETWEEN TWO ADJACENT LAYERS

!     PURPOSE.
!     --------
!     COMPUTES THE LIQUID-WATER-EQUIVALENT VOLUMETRIC FLUX OF WATER
!     VAPOUR DIFFUSING BETWEEN TWO SOIL LAYERS "A" (UPPER) AND "B"
!     (LOWER), DRIVEN JOINTLY BY THE MATRIC-POTENTIAL (ISOTHERMAL) AND
!     TEMPERATURE (THERMAL) GRADIENTS ACROSS THE INTERFACE. THIS IS
!     THE CLASSICAL TWO-TERM PHILIP & DE VRIES FORMULATION AND DOES
!     NOT INTRODUCE A NEW PROGNOSTIC STATE (NO EXPLICIT DRY-AIR
!     PRESSURE, UNLIKE THE FULLER 3-PHASE FORMULATIONS USED BY SOME
!     RESEARCH SOIL-VEGETATION-ATMOSPHERE MODELS, E.G. STEMMUS).

!     THE RESULT IS MEANT TO BE ADDED, AS AN EXPLICIT SOURCE TERM, TO
!     THE RIGHT-HAND SIDE OF THE EXISTING (LIQUID-ONLY) SOIL MOISTURE
!     TRIDIAGONAL SYSTEM -- SEE *SRFWVAPOR* (THE VECTORISED DRIVER
!     THAT CALLS THIS ROUTINE ONCE PER INTERFACE) FOR HOW IT IS USED.
!     THE LIQUID DIFFUSIVITY/CONDUCTIVITY MATRIX BUILT BY *SRFWEXC_VG*
!     IS LEFT ENTIRELY UNTOUCHED.

!     METHOD.
!     -------
!     Vapour density in the pore space:      rho_v = rho_vsat(T) * Hr
!     Kelvin relative humidity:               Hr    = exp(h*g/(Rv*T))
!     Effective vapour diffusivity in soil:   Dveff = tau*theta_a*Datm(T)
!       tau   = Millington & Quirk (1961) tortuosity = theta_a**(7/3)/theta_s**2
!       Datm  = binary diffusivity of vapour in air, Campbell (1985):
!               Datm(T) = 2.12E-5 * (T/273.16)**2                 [m2/s]
!     Flux (liquid-equivalent, positive DOWNWARD, m/s):
!       q_vap = -Kvt*dT/dz - Kvh*dh/dz
!       Kvt   = (Dveff/rho_liquid) * Hr * d(rho_vsat)/dT
!       Kvh   = (Dveff/rho_liquid) * rho_vsat * Hr*g/(Rv*T)
!     Matric head h(theta) from the van Genuchten retention curve
!     already used elsewhere in the soil column (*SRFWEXC_VG*),
!     evaluated independently on each side of the interface using
!     the hydraulic parameters of that layer.

!     KNOWN SIMPLIFICATIONS (v1, offline/research prototype -- see the
!     STEMMUS-SCOPE -> ecLand feasibility note this implements):
!     - No Cass, Campbell & Jones (1984) thermal-vapour "enhancement
!       factor" (typically 3-10x): would need a continuous clay
!       fraction not currently held per soil texture class in
!       *YOS_SOIL*. Omitting it makes this a conservative estimate.
!     - This kernel is INTERNAL soil-column redistribution only; the
!       surface (soil-atmosphere) evaporative flux, computed upstream
!       of the soil water budget, is not touched.

!     REFERENCE.
!     ----------
!     Philip, J.R. & de Vries, D.A. (1957), Trans. AGU 38(2).
!     Cass, A., Campbell, G.S. & Jones, T.L. (1984), Soil Sci. Soc. Am. J.
!     Campbell, G.S. (1985), Soil Physics with BASIC.
!     Millington, R.J. & Quirk, J.P. (1961), Trans. Faraday Soc. 57.

!     Original
!     Prototype for the STEMMUS-SCOPE -> ecLand dryland vapour-flux
!     feasibility study.  2026-08-28
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB),    INTENT(IN)   :: PTHETA_A
REAL(KIND=JPRB),    INTENT(IN)   :: PTHETA_B
REAL(KIND=JPRB),    INTENT(IN)   :: PTEMP_A
REAL(KIND=JPRB),    INTENT(IN)   :: PTEMP_B
REAL(KIND=JPRB),    INTENT(IN)   :: PDZ
REAL(KIND=JPRB),    INTENT(IN)   :: PWSAT_A
REAL(KIND=JPRB),    INTENT(IN)   :: PWREST_A
REAL(KIND=JPRB),    INTENT(IN)   :: PALPHA_A
REAL(KIND=JPRB),    INTENT(IN)   :: PN_A
REAL(KIND=JPRB),    INTENT(IN)   :: PM_A
REAL(KIND=JPRB),    INTENT(IN)   :: PWSAT_B
REAL(KIND=JPRB),    INTENT(IN)   :: PWREST_B
REAL(KIND=JPRB),    INTENT(IN)   :: PALPHA_B
REAL(KIND=JPRB),    INTENT(IN)   :: PN_B
REAL(KIND=JPRB),    INTENT(IN)   :: PM_B
TYPE(TCST),         INTENT(IN)   :: YDCST

REAL(KIND=JPRB),    INTENT(OUT)  :: PQVAP

!*         0.2    DECLARATION OF LOCAL VARIABLES.
!                 ----------- -- ----- ----------

REAL(KIND=JPRB) :: ZHA, ZHB, ZHINT, ZTINT, ZTHETAINT, ZWSATINT, ZTHETAAIR
REAL(KIND=JPRB) :: ZTAU, ZDATM, ZDVEFF
REAL(KIND=JPRB) :: ZESAT, ZDESDT, ZRHOSAT, ZDRHOSATDT
REAL(KIND=JPRB) :: ZHR, ZDHRDH, ZKVT, ZKVH
REAL(KIND=JPRB) :: ZDTDZ, ZDHDZ
REAL(KIND=JPRB), PARAMETER :: ZEPS = 1.E-8_JPRB

#include "fcsttre.h"

!     ------------------------------------------------------------------

ASSOCIATE(RTT=>YDCST%RTT, RG=>YDCST%RG, RV=>YDCST%RV)

!*         1.    MATRIC HEAD ON EACH SIDE OF THE INTERFACE (VAN GENUCHTEN).

ZHA = VGHEAD(PTHETA_A,PWREST_A,PWSAT_A,PALPHA_A,PN_A,PM_A)
ZHB = VGHEAD(PTHETA_B,PWREST_B,PWSAT_B,PALPHA_B,PN_B,PM_B)

!*         2.    INTERFACE STATE (SIMPLE AVERAGE OF THE TWO LAYERS).

ZTINT     = 0.5_JPRB*(PTEMP_A+PTEMP_B)
ZHINT     = 0.5_JPRB*(ZHA+ZHB)
ZWSATINT  = 0.5_JPRB*(PWSAT_A+PWSAT_B)
ZTHETAINT = 0.5_JPRB*(PTHETA_A+PTHETA_B)
ZTHETAAIR = MAX(ZWSATINT-ZTHETAINT,ZEPS)

!*         3.    EFFECTIVE VAPOUR DIFFUSIVITY IN THE SOIL AIR SPACE.

ZTAU  = ZTHETAAIR**(7.0_JPRB/3.0_JPRB)/MAX(ZWSATINT**2,ZEPS)
ZDATM = 2.12E-5_JPRB*(ZTINT/RTT)**2
ZDVEFF= ZTAU*ZTHETAAIR*ZDATM

!*         4.    SATURATION VAPOUR DENSITY AND ITS TEMPERATURE DERIVATIVE
!                (REUSING THE SAME TETENS ESAT FUNCTION USED ELSEWHERE IN
!                THE MODEL, FOR CONSISTENCY).

ZESAT      = FOEEW(ZTINT)
ZDESDT     = ZESAT*FOEDESU(ZTINT)
ZRHOSAT    = ZESAT/(RV*ZTINT)
ZDRHOSATDT = ZRHOSAT*(FOEDESU(ZTINT)-1.0_JPRB/ZTINT)

!*         5.    KELVIN RELATIVE-HUMIDITY FACTOR AND ITS HEAD DERIVATIVE.

ZHR    = EXP(ZHINT*RG/(RV*ZTINT))
ZDHRDH = ZHR*RG/(RV*ZTINT)

!*         6.    THERMAL AND ISOTHERMAL VAPOUR "CONDUCTIVITIES".

ZKVT = (ZDVEFF/RHOH2O)*ZHR*ZDRHOSATDT
ZKVH = (ZDVEFF/RHOH2O)*ZRHOSAT*ZDHRDH

!*         7.    NET LIQUID-EQUIVALENT VAPOUR FLUX (POSITIVE DOWNWARD).

ZDTDZ = (PTEMP_B-PTEMP_A)/PDZ
ZDHDZ = (ZHB-ZHA)/PDZ

PQVAP = -(ZKVT*ZDTDZ+ZKVH*ZDHDZ)

END ASSOCIATE

!     ------------------------------------------------------------------

CONTAINS

REAL(KIND=JPRB) FUNCTION VGHEAD(PTHETA,PWREST,PWSAT,PALPHA,PN,PM) RESULT(ZH)

! Van Genuchten (1980) matric head h(theta), inverted from the same
! retention-curve parameters (alpha, n, m) already used for the liquid
! diffusivity/conductivity in *SRFWEXC_VG*. Se is clipped away from 0
! and 1 to keep h finite at/near saturation and residual moisture.

REAL(KIND=JPRB), INTENT(IN) :: PTHETA, PWREST, PWSAT, PALPHA, PN, PM
REAL(KIND=JPRB) :: ZSE
REAL(KIND=JPRB), PARAMETER :: ZSEEPS = 1.E-4_JPRB

ZSE = (PTHETA-PWREST)/MAX(PWSAT-PWREST,1.E-6_JPRB)
ZSE = MIN(MAX(ZSE,ZSEEPS),1.0_JPRB-ZSEEPS)
ZH  = -(1.0_JPRB/PALPHA)*((ZSE**(-1.0_JPRB/PM)-1.0_JPRB)**(1.0_JPRB/PN))

END FUNCTION VGHEAD

END SUBROUTINE SRFWVAPOR_COEF
END MODULE SRFWVAPOR_COEF_MOD
