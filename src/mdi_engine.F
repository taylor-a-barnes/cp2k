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

   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length
   USE machine,                         ONLY: default_output_unit
   USE string_utilities,                ONLY: integer_to_string

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: init_mdi

   LOGICAL, PUBLIC :: use_mdi

CONTAINS

! **************************************************************************************************
!> \brief Perform initialization to run as an MDI engine
!> \param ...
! **************************************************************************************************
   SUBROUTINE init_mdi()

      WRITE (default_output_unit, *) "AAA INITIALIZING MDI"      

   END SUBROUTINE init_mdi

END MODULE mdi_engine
