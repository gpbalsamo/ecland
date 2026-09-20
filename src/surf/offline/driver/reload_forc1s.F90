SUBROUTINE RELOAD_FORC1S

USE YOMFORC1S, ONLY : UFI      ,VFI      ,TFI      ,QFI      ,&
     &            PSFI     ,SRFFI    ,TRFFI    ,R30FI    ,S30FI    ,&
     &            R30FI_C  ,S30FI_C  ,CO2FI    ,NSTPFC   ,&
     &            CFFORCU  ,CFFORCV  ,CFFORCT  ,CFFORCQ  ,CFFORCC  ,&
     &            CFFORCP  ,CFFORCRAIN,CFFORCSNOW,CFFORCSW,CFFORCLW
USE YOMCST   , ONLY : RTT ,RMD, RMCO2
USE YOMLUN1S , ONLY : NULOUT
USE YOMDPHY  , ONLY : NPOI
USE YOEPHY   , ONLY : LEAIRCO2COUP
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : DR_HOOK, JPHOOK, LHOOK

#ifdef DOC
! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *RELOAD_FORC1S * - READ ONE WINDOW OF ATMOSPHERIC FORCING DATA
!
!     PURPOSE.
!     --------
!        Read (or re-read) the currently-needed window of every forcing
!        variable via RDFVAR, and apply the same variable-derivation logic
!        SUFCDF has always applied (2D-vs-scalar wind, snowfall deduced
!        from Tair when absent, convective-precipitation-fraction split,
!        CO2 unit conversion). Factored out of SUFCDF so the exact same
!        sequence can be called again mid-run by DTFORC's windowed-refill
!        guard (NFORCWINDOW>0) without the two ever drifting out of sync.
!
!        Deliberately NOT included here: the NAMFORC namelist read and the
!        LOADIAB adiabatic-height correction, both of which stay in SUFCDF
!        only -- LOADIAB is explicitly incompatible with windowed reads
!        (SUFCDF aborts up front if both are requested together; see its
!        own header comment).
!
!     AUTHOR.
!     -------
!        Added 2026-09-20 alongside NFORCWINDOW (see yomforc1s.F90,
!        rdfvar.F90, dtforc.F90) -- see wfde5-ecland's PLAN.md Milestone
!        4/5 for the real-world need (long/high-resolution runs otherwise
!        holding the entire declared forcing period in memory at once).
#endif
IMPLICIT NONE
#include "rdfvar.intfb.h"
#include "netcdf.inc"

CHARACTER*100 CNAME
INTEGER NCID,IERR,NVARS,IDUM1,IDUM2,IDUM3,IDUMAR(20),IVAR
INTEGER(KIND=JPIM) :: JL,JT
REAL(KIND=JPRB) :: ZFRAC
LOGICAL LWIND2D,LSNOWF,LCTPF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RELOAD_FORC1S',0,ZHOOK_HANDLE)

!* 2D (Wind_E/Wind_N) vs scalar (Wind) wind forcing
LWIND2D=.FALSE.
NCID = NCOPN(CFFORCU, NCNOWRIT, IERR)
IF( IERR == 0 ) THEN
  CALL NCINQ(NCID, IDUM1, NVARS, IDUM2, IDUM3, IERR)
  DO IVAR=1,NVARS
    CALL NCVINQ(NCID,IVAR,CNAME,IDUM1,IDUM2,IDUMAR,IDUM3,IERR)
    IF(CNAME(1:6).EQ.'Wind_E') LWIND2D=.TRUE.
  ENDDO
ENDIF
CALL NCCLOS(NCID,IERR)

IF(LWIND2D)THEN
  CNAME='Wind_E'
  CALL RDFVAR(CFFORCU,CNAME,UFI)
  CNAME='Wind_N'
  CALL RDFVAR(CFFORCV,CNAME,VFI)
ELSE
  CNAME='Wind'
  CALL RDFVAR(CFFORCU,CNAME,UFI)
  VFI=0.0_JPRB
ENDIF

!* temperature
CNAME='Tair'
CALL RDFVAR(CFFORCT,CNAME,TFI)

!* specific humidity
CNAME='Qair'
CALL RDFVAR(CFFORCQ,CNAME,QFI)

IF (LEAIRCO2COUP) THEN
  CNAME='CO2air'
  CALL RDFVAR(CFFORCC,CNAME,CO2FI)
ENDIF

!* surface pressure
CNAME='PSurf'
CALL RDFVAR(CFFORCP,CNAME,PSFI)

!* shortwave / longwave radiation
CNAME='SWdown'
CALL RDFVAR(CFFORCSW,CNAME,SRFFI)
CNAME='LWdown'
CALL RDFVAR(CFFORCLW,CNAME,TRFFI)

!* rainfall (all large scale)
CNAME='Rainf'
CALL RDFVAR(CFFORCRAIN,CNAME,R30FI)

!* snowfall: read if present, else deduce from Tair/Rainf
LSNOWF=.FALSE.
NCID = NCOPN(CFFORCSNOW, NCNOWRIT, IERR)
IF( IERR == 0 ) THEN
  CALL NCINQ(NCID, IDUM1, NVARS, IDUM2, IDUM3, IERR)
  DO IVAR=1,NVARS
    CALL NCVINQ(NCID,IVAR,CNAME,IDUM1,IDUM2,IDUMAR,IDUM3,IERR)
    IF(CNAME(1:5).EQ.'Snowf') LSNOWF=.TRUE.
  ENDDO
ENDIF
CALL NCCLOS(NCID,IERR)

IF (LSNOWF) THEN
  CNAME='Snowf'
  CALL RDFVAR(CFFORCSNOW,CNAME,S30FI)
ELSE
  WRITE(NULOUT,*) " SNOWF DEDUCED FROM TAIR AND RAINF "
  S30FI=0.0_JPRB
  DO JL=1,NPOI
  DO JT=1,NSTPFC
    IF(TFI(JL,JT).LT.RTT)THEN
      S30FI(JL,JT)=R30FI(JL,JT)
      R30FI(JL,JT)=0.0_JPRB
    ENDIF
  ENDDO
  ENDDO
ENDIF

!* convective precipitation fraction: split Rainf/Snowf if present, else zero
LCTPF=.FALSE.
NCID = NCOPN(CFFORCRAIN, NCNOWRIT, IERR)
IF( IERR == 0 ) THEN
  CALL NCINQ(NCID, IDUM1, NVARS, IDUM2, IDUM3, IERR)
  DO IVAR=1,NVARS
    CALL NCVINQ(NCID,IVAR,CNAME,IDUM1,IDUM2,IDUMAR,IDUM3,IERR)
    IF(CNAME(1:4).EQ.'Ctpf') LCTPF=.TRUE.
  ENDDO
ENDIF
CALL NCCLOS(NCID,IERR)

IF (LCTPF) THEN
  CNAME='Ctpf'
  ! R30FI_C/S30FI_C double as the read buffer here, exactly as SUFCDF did
  ! (ALLOCATE(ZREAL3D(NPOI,JPSTPFC)) there was only ever used transiently
  ! for this one read) -- overwritten immediately below with the true
  ! convective-fraction split, so no separate scratch array is needed.
  CALL RDFVAR(CFFORCRAIN,CNAME,R30FI_C)
  IF (MINVAL(R30FI_C(:,1:NSTPFC)).LT.0.0_JPRB) THEN
    WRITE(NULOUT,*) " CONVECTIVE FRACTION < 0 "
  ENDIF
  IF (MAXVAL(R30FI_C(:,1:NSTPFC)).GT.1.0_JPRB) THEN
    WRITE(NULOUT,*) " CONVECTIVE FRACTION > 1 "
  ENDIF
  WRITE(NULOUT,*) "Correcting Rainf/Snoww using Ctpf "
  DO JL=1,NPOI
    DO JT = 1, NSTPFC
      ZFRAC=MIN(1.0_JPRB,MAX(0._JPRB,R30FI_C(JL,JT)))
      S30FI_C(JL,JT)=ZFRAC*S30FI(JL,JT)
      R30FI_C(JL,JT)=ZFRAC*R30FI(JL,JT)
      R30FI(JL,JT)=(1.0_JPRB - ZFRAC)*R30FI(JL,JT)
      S30FI(JL,JT)=(1.0_JPRB - ZFRAC)*S30FI(JL,JT)
    ENDDO
  ENDDO
ELSE
  WRITE(NULOUT,*) " CONVECTIVE TP NOT PRESENT, SET TO 0. "
  R30FI_C(:,1:NSTPFC)=0.0_JPRB
  S30FI_C(:,1:NSTPFC)=0.0_JPRB
ENDIF

!* convert atmospheric CO2 from ppm to kg/kg
IF (LEAIRCO2COUP) THEN
  DO JL=1,NPOI
  DO JT=1,NSTPFC
    CO2FI(JL,JT)=CO2FI(JL,JT)*RMCO2/(RMD*1000000._JPRB)
  ENDDO
  ENDDO
ENDIF

WRITE(NULOUT,*) " FORCING DATA (RE)LOADED FOR ",NSTPFC," FORCING STEPS"

IF (LHOOK) CALL DR_HOOK('RELOAD_FORC1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE RELOAD_FORC1S
