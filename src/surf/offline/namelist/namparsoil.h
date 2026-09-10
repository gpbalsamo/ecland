! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! * SOIL PARAMETERS
!     ------------------------------------------------------------------
NAMELIST/NAMPARSOIL/&
! * soil related constants
& RTF1, RTF2, RLAMBDAICE, RLAMBDAWAT, RKERST1, RKERST2&
& , RKERST3, RSRDEP, RSIGORMIN, RSIGORMAX, RWLMAX, RPSFR, RBARPWP, RVGBARCAP&
& , RBARCAP, RCLU, RRSF1A, RRSF1B, RRSF1C, RRHOSM, RLAMBDAQ, RLAMBDAO&
& , RLAMBDAMIN, RLAMBDAMAX, RLAMBDRYMA, RLAMBDRYMB, RLAMBDRYMC&
& , RTFREEZSICECEL, RTMELTSICECEL&
! * "depth trilogy" prototype switches -- see yos_soil.F90 for what each does.
!   All default to .FALSE. below (RDNML_SOIL), i.e. unset in a namelist means
!   plain baseline ecLand -- see the block comment ahead of the WRITE lines in
!   RDNML_SOIL for why running all six .TRUE. together is still expected to
!   work.
& , LEWVFLUX, LEWPFLOOR, LEFRZFLOOR, LEUNIFORMROOT, LEBEDROCKLIM, LEGWRECHARGE&
! * LEGWRECHARGE recharge throttle -- see yos_soil.F90. Defaults to 1.0 in
!   RDNML_SOIL, i.e. unset means the pre-existing behaviour bit for bit.
& , RWTDRECHARGE&
! * frozen-soil macropore-permeability cap -- see yos_soil.F90. Defaults to
!   .FALSE. below (RDNML_SOIL), i.e. unset means plain baseline ecLand.
& , LEFROZENSOILFIX
!     ------------------------------------------------------------------
