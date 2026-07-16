SUBROUTINE SURFDEALLO(YDSURF)

! (C) Copyright 2001- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURFDEALLO * - Routine to deallocate space for global variables
!                     from surface subroutine
!     Purpose.
!     --------
!           Deallocate space for the global fields.

!**   Interface.
!     ----------
!        *CALL* *SURFDEALLO*

!     Explicit arguments :
!     --------------------
!        YDSURF : surface metadata object (unused; kept for API consistency)

!     Implicit arguments :
!     --------------------
!        Pointers of comdecks

!     Method.
!     -------
!        Extraction from DEALLO

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        J.F. Estrade *ECMWF* 03-10-01
!     Modifications.
!     --------------
!        Original : 03-10-01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB
USE YOS_SURF, ONLY : TSURF

!ifndef INTERFACE

USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
!ifndef INTERFACE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SURFDEALLO',0,ZHOOK_HANDLE)
IF (LHOOK) CALL DR_HOOK('SURFDEALLO',1,ZHOOK_HANDLE)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE SURFDEALLO
