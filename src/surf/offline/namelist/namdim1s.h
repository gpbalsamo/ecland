! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!*    -----------------------------------------------------------------
! NFORCWINDOW: number of forcing records to keep resident at once, 0 (the
! default) means "load the whole NDFORC-length series", exactly as before
! this parameter existed. See yomforc1s.F90, rdfvar.F90, reload_forc1s.F90
! and dtforc.F90 for the windowed-read design this enables.
NAMELIST/NAMDIM/ NLON,NLAT,NDFORC,NCOOR,NPROMA,NCSNEC,NCSS,NFORCWINDOW
!     -----------------------------------------------------------------
