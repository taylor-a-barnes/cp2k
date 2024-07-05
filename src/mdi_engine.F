!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief support for CP2K to run as an MDI engine
!> \par History
!>      - created (2024-07, Taylor Barnes)
!> \author Taylor Barnes
! **************************************************************************************************
MODULE mdi_engine

   USE kinds,                           ONLY: default_path_length
   USE machine,                         ONLY: default_output_unit
   USE mdi,                             ONLY: MDI_Init

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: init_mdi

   LOGICAL, PUBLIC :: use_mdi
   CHARACTER(LEN=default_path_length), PUBLIC :: mdi_options

CONTAINS

! **************************************************************************************************
!> \brief Perform initialization to run as an MDI engine
!> \param ...
! **************************************************************************************************
   SUBROUTINE init_mdi()
      INTEGER :: ierr

      WRITE (default_output_unit, *) "AAA INITIALIZING MDI"      

      CALL MDI_Init(mdi_options, ierr)
      IF ( ierr .NE. 0 ) THEN
         WRITE (default_output_unit, "(/,T2,A)") &
            "ERROR: MDI_Init returned a non-zero error code"
      END IF

      WRITE (default_output_unit, *) "BBB INITIALIZING MDI"    

   END SUBROUTINE init_mdi

END MODULE mdi_engine
