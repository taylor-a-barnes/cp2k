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
                                              dp
   USE machine,                         ONLY: default_output_unit
   USE mdi,                             ONLY: MDI_Init,&
                                              MDI_Accept_communicator,&
                                              MDI_Recv,&
                                              MDI_Send,&
                                              MDI_Recv_command,&
                                              MDI_INT,&
                                              MDI_DOUBLE,&
                                              MDI_COMMAND_LENGTH
   USE f77_interface,                   ONLY: calc_energy_force,&
                                              get_cell,&
                                              get_energy,&
                                              get_natom,&
                                              get_nparticle,&
                                              get_force,&
                                              get_pos,&
                                              get_stress_tensor,&
                                              set_cell,&
                                              set_pos,&
                                              set_vel
   USE ISO_C_binding

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: init_mdi, mdi_listen

   LOGICAL, PUBLIC :: use_mdi
   CHARACTER(LEN=default_path_length), PUBLIC :: mdi_options

   LOGICAL, PRIVATE :: terminate_flag, mdi_forces_current
   INTEGER, PRIVATE :: mdi_comm, mdi_env_id
   REAL(KIND=dp), ALLOCATABLE :: mdi_coords_buffer(:)

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

      mdi_forces_current = .FALSE.

   END SUBROUTINE init_mdi

! **************************************************************************************************
!> \brief Begin listening for commands from the driver
!> \param env_id id of the force_env
! **************************************************************************************************
   SUBROUTINE mdi_listen(env_id)
      CHARACTER(len=:), ALLOCATABLE :: command
      INTEGER, INTENT(IN)           :: env_id
      INTEGER                       :: ierr

      TYPE(C_PTR)                   :: class_obj

      ALLOCATE( character(MDI_COMMAND_LENGTH) :: command )
      mdi_env_id = env_id

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

! **************************************************************************************************
!> \brief Ensure that the energy and forces are up to date
!> \param ...
! **************************************************************************************************
   SUBROUTINE mdi_update_forces()
      IMPLICIT NONE

      INTEGER                       :: ierr

      IF ( .NOT. mdi_forces_current ) THEN
         CALL calc_energy_force(mdi_env_id, .TRUE., ierr)
      END IF
      mdi_forces_current = .TRUE.

   END SUBROUTINE mdi_update_forces

! **************************************************************************************************
!> \brief Execute an MDI command received from the driver
!> \param command command received from the driver
!> \param comm MDI communicator used to communicate with the driver
!> \param class_obj dummy argument used as part of the MDI plugin system
! **************************************************************************************************
   FUNCTION execute_command(command, comm, class_obj)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)   :: command
      INTEGER, INTENT(IN)            :: comm
      TYPE(C_PTR), VALUE             :: class_obj
      INTEGER(KIND=C_INT)            :: execute_command
      INTEGER                        :: ierr

      REAL(kind=dp)                  :: energy
      INTEGER                        :: natoms, nparticles
      REAL(KIND=dp), DIMENSION(3, 3) :: cell_tensor
      REAL(KIND=dp)                  :: cell_vector(9)
      INTEGER, DIMENSION(3)          :: periodicity

      WRITE (default_output_unit, "(/,T2,A)") &
        "MDI Command received: ",TRIM(command)

      IF ( .NOT. ALLOCATED(mdi_coords_buffer) ) THEN
         CALL get_natom(mdi_env_id, natoms, ierr)
         CALL get_nparticle(mdi_env_id, nparticles, ierr)
         IF ( natoms .NE. nparticles ) THEN
            WRITE (default_output_unit, "(/,T2,A)") &
               "ERROR: The MDI interface requires that the number of atoms equal the number of particles"
            STOP 1
         END IF
         ALLOCATE( mdi_coords_buffer(3*nparticles) )
      END IF

      SELECT CASE( TRIM(command) )
      CASE( "EXIT" )
         terminate_flag = .true.
      CASE( "<ENERGY" )
         CALL mdi_update_forces()
         CALL get_energy(mdi_env_id, energy, ierr)
         CALL MDI_Send(energy, 1, MDI_DOUBLE, comm, ierr)
      CASE( "<PE" )
         CALL mdi_update_forces()
         CALL get_energy(mdi_env_id, energy, ierr)
         CALL MDI_Send(energy, 1, MDI_DOUBLE, comm, ierr)
      CASE( "<NATOMS" )
         CALL get_natom(mdi_env_id, natoms, ierr)
         CALL MDI_Send(natoms, 1, MDI_INT, comm, ierr)
      CASE( "<NPARTICLES" )
         CALL get_nparticle(mdi_env_id, nparticles, ierr)
         CALL MDI_Send(nparticles, 1, MDI_INT, comm, ierr)
      CASE( "<FORCES" )
         CALL mdi_update_forces()
         CALL get_nparticle(mdi_env_id, nparticles, ierr)
         CALL get_force(mdi_env_id, mdi_coords_buffer, 3 * nparticles, ierr)
         CALL MDI_Send(mdi_coords_buffer, 3 * nparticles, MDI_DOUBLE, comm, ierr)
      CASE( "<CELL" )
         CALL get_cell(mdi_env_id, cell_tensor, periodicity, ierr)
         cell_vector(1) = cell_tensor(1,1)
         cell_vector(2) = cell_tensor(1,2)
         cell_vector(3) = cell_tensor(1,3)
         cell_vector(4) = cell_tensor(2,1)
         cell_vector(5) = cell_tensor(2,2)
         cell_vector(6) = cell_tensor(2,3)
         cell_vector(7) = cell_tensor(3,1)
         cell_vector(8) = cell_tensor(3,2)
         cell_vector(9) = cell_tensor(3,3)
         CALL MDI_Send(cell_vector, 9, MDI_DOUBLE, comm, ierr)
      CASE( ">CELL" )
         mdi_forces_current = .FALSE.
         CALL MDI_Recv(cell_vector, 9, MDI_DOUBLE, comm, ierr)
         cell_tensor(1,1) = cell_vector(1)
         cell_tensor(1,2) = cell_vector(2)
         cell_tensor(1,3) = cell_vector(3)
         cell_tensor(2,1) = cell_vector(4)
         cell_tensor(2,2) = cell_vector(5)
         cell_tensor(2,3) = cell_vector(6)
         cell_tensor(3,1) = cell_vector(7)
         cell_tensor(3,2) = cell_vector(8)
         cell_tensor(3,3) = cell_vector(9)
         CALL set_cell(mdi_env_id, cell_tensor, ierr)
      CASE( "<COORDS" )
         CALL get_nparticle(mdi_env_id, nparticles, ierr)
         CALL get_pos(mdi_env_id, mdi_coords_buffer, 3 * nparticles, ierr)
         CALL MDI_Send(mdi_coords_buffer, 3 * nparticles, MDI_DOUBLE, comm, ierr)
      CASE( ">COORDS" )
         mdi_forces_current = .FALSE.
         CALL get_nparticle(mdi_env_id, nparticles, ierr)
         CALL MDI_Recv(mdi_coords_buffer, 3 * nparticles, MDI_DOUBLE, comm, ierr)
         CALL set_pos(mdi_env_id, mdi_coords_buffer, 3 * nparticles, ierr)
      CASE( ">VELOCITIES" )
         CALL get_nparticle(mdi_env_id, nparticles, ierr)
         CALL MDI_Recv(mdi_coords_buffer, 3 * nparticles, MDI_DOUBLE, comm, ierr)
         CALL set_vel(mdi_env_id, mdi_coords_buffer, 3 * nparticles, ierr)
      CASE( "<STRESS" )
         CALL mdi_update_forces()
         CALL get_stress_tensor(mdi_env_id, cell_tensor, ierr)
         cell_vector(1) = cell_tensor(1,1)
         cell_vector(2) = cell_tensor(1,2)
         cell_vector(3) = cell_tensor(1,3)
         cell_vector(4) = cell_tensor(2,1)
         cell_vector(5) = cell_tensor(2,2)
         cell_vector(6) = cell_tensor(2,3)
         cell_vector(7) = cell_tensor(3,1)
         cell_vector(8) = cell_tensor(3,2)
         cell_vector(9) = cell_tensor(3,3)
         CALL MDI_Send(cell_vector, 9, MDI_DOUBLE, comm, ierr)
      CASE DEFAULT
         WRITE (default_output_unit, "(/,T2,A)") &
            "ERROR: Command not recognized."
         STOP 1
      END SELECT

      execute_command = 0
   END FUNCTION execute_command

END MODULE mdi_engine
