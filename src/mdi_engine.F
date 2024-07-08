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
   USE mdi,                             ONLY: MDI_Init,&
                                              MDI_Accept_communicator,&
                                              MDI_Recv,&
                                              MDI_Send,&
                                              MDI_Recv_command,&
                                              MDI_INT,&
                                              MDI_DOUBLE,&
                                              MDI_COMMAND_LENGTH
   USE ISO_C_binding

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: init_mdi, mdi_listen

   LOGICAL, PUBLIC :: use_mdi
   CHARACTER(LEN=default_path_length), PUBLIC :: mdi_options

   LOGICAL, PRIVATE :: terminate_flag
   INTEGER, PRIVATE :: mdi_comm

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
         STOP 1
      END IF

      WRITE (default_output_unit, *) "BBB INITIALIZING MDI"

      CALL MDI_Accept_communicator(mdi_comm, ierr)
      IF ( ierr .NE. 0 ) THEN
         WRITE (default_output_unit, "(/,T2,A)") &
            "ERROR: MDI_Accept_communicator returned a non-zero error code"
         STOP 1
      END IF

   END SUBROUTINE init_mdi

   SUBROUTINE mdi_listen()
      CHARACTER(len=:), ALLOCATABLE :: command
      INTEGER                       :: ierr

      TYPE(C_PTR)                   :: class_obj

      ALLOCATE( character(MDI_COMMAND_LENGTH) :: command )

      ! Respond to the driver's commands
      response_loop: DO

         ! Receive a command from the driver and broadcast it to all ranks
         CALL MDI_Recv_command(command, mdi_comm, ierr)
         !CALL MPI_Bcast(command, MDI_COMMAND_LENGTH, MPI_CHAR, 0, world_comm, ierr)

         ierr = execute_command(command, mdi_comm, class_obj)

         IF ( terminate_flag ) EXIT

      END DO response_loop

      DEALLOCATE( command )

   END SUBROUTINE mdi_listen

   FUNCTION execute_command(command, comm, class_obj)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: command
      INTEGER, INTENT(IN)           :: comm
      TYPE(C_PTR), VALUE            :: class_obj
      INTEGER(KIND=C_INT)           :: execute_command
      !INTEGER                       :: ierr

      WRITE (default_output_unit, "(/,T2,A)") &
        "MDI Command received: ",TRIM(command)

      SELECT CASE( TRIM(command) )
      CASE( "EXIT" )
         terminate_flag = .true.
      !CASE( "<NATOMS" )
      !   CALL MDI_Send(natoms, 1, MDI_INT, comm, ierr)
      !CASE( "<COORDS" )
      !   CALL MDI_Send(coords, 3 * natoms, MDI_DOUBLE, comm, ierr)
      !CASE( ">COORDS" )
      !   CALL MDI_Recv(coords, 3 * natoms, MDI_DOUBLE, comm, ierr)
      !CASE( "<FORCES" )
      !   CALL MDI_Send(forces, 3 * natoms, MDI_DOUBLE, comm, ierr)
      CASE DEFAULT
         WRITE (default_output_unit, "(/,T2,A)") &
            "ERROR: Command not recognized."
         STOP 1
      END SELECT

      execute_command = 0
   END FUNCTION execute_command


END MODULE mdi_engine
