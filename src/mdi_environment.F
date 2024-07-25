!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods and functions on the MDI environment
!> \par History
!>      - created (2024-07, Taylor Barnes)
!> \author Taylor Barnes
! **************************************************************************************************
MODULE mdi_environment

   USE kinds,                           ONLY: default_path_length
   USE machine,                         ONLY: default_output_unit
   USE mdi,                             ONLY: MDI_Init,&
                                              MDI_Accept_communicator,&
                                              MDI_MPI_get_world_comm,&
                                              MDI_Register_node,&
                                              MDI_Register_command

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: initialize_mdi

   LOGICAL, PUBLIC :: use_mdi, mdi_forces_current
   CHARACTER(LEN=default_path_length), PUBLIC :: mdi_options
   INTEGER, PUBLIC :: mdi_comm, mdi_env_id

CONTAINS

! **************************************************************************************************
!> \brief Perform initialization to run as an MDI engine
!> \param ...
! **************************************************************************************************
   SUBROUTINE initialize_mdi()
      INTEGER :: ierr

      CALL MDI_Init(mdi_options, ierr)
      IF ( ierr .NE. 0 ) THEN
         WRITE (default_output_unit, "(/,T2,A)") &
            "ERROR: MDI_Init returned a non-zero error code"
         STOP 1
      END IF

      CALL MDI_Register_node("@DEFAULT", ierr);
      CALL MDI_Register_command("@DEFAULT", "EXIT", ierr);
      CALL MDI_Register_command("@DEFAULT", "<ENERGY", ierr);
      CALL MDI_Register_command("@DEFAULT", "<PE", ierr);
      CALL MDI_Register_command("@DEFAULT", "<NATOMS", ierr);
      CALL MDI_Register_command("@DEFAULT", "<NPARTICLES", ierr);
      CALL MDI_Register_command("@DEFAULT", "<FORCES", ierr);
      CALL MDI_Register_command("@DEFAULT", "<CELL", ierr);
      CALL MDI_Register_command("@DEFAULT", ">CELL", ierr);
      CALL MDI_Register_command("@DEFAULT", "<COORDS", ierr);
      CALL MDI_Register_command("@DEFAULT", ">COORDS", ierr);
      CALL MDI_Register_command("@DEFAULT", ">VELOCITIES", ierr);
      CALL MDI_Register_command("@DEFAULT", "<STRESS", ierr);

      CALL MDI_Accept_communicator(mdi_comm, ierr)
      IF ( ierr .NE. 0 ) THEN
         WRITE (default_output_unit, "(/,T2,A)") &
            "ERROR: MDI_Accept_communicator returned a non-zero error code"
         STOP 1
      END IF

      mdi_forces_current = .FALSE.

   END SUBROUTINE initialize_mdi

END MODULE mdi_environment
