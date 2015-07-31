!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Driver for the PdV update.
!>  @author Wayne Gaudin
!>  @details Invokes the user specified kernel for the PdV update.

MODULE PdV_module

CONTAINS

SUBROUTINE PdV(predict)

  USE clover_module
  USE report_module
  USE PdV_kernel_module
  USE revert_module
  USE update_halo_module
  USE ideal_gas_module

  IMPLICIT NONE

  LOGICAL :: predict

  INTEGER :: prdct

  INTEGER :: t
  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: kernel_time,timer

  error_condition=0 ! Not used yet due to issue with OpenA reduction

  IF(profiler_on) kernel_time=timer()

  IF(use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL PdV_kernel(predict,                  &
                    chunk%tiles(t)%field%x_min,      &
                    chunk%tiles(t)%field%x_max,      &
                    chunk%tiles(t)%field%y_min,      &
                    chunk%tiles(t)%field%y_max,      &
                    dt,                         &
                    chunk%tiles(t)%field%xarea,      &
                    chunk%tiles(t)%field%yarea,      &
                    chunk%tiles(t)%field%volume ,    &
                    chunk%tiles(t)%field%density0,   &
                    chunk%tiles(t)%field%density1,   &
                    chunk%tiles(t)%field%energy0,    &
                    chunk%tiles(t)%field%energy1,    &
                    chunk%tiles(t)%field%pressure,   &
                    chunk%tiles(t)%field%viscosity,  &
                    chunk%tiles(t)%field%xvel0,      &
                    chunk%tiles(t)%field%xvel1,      &
                    chunk%tiles(t)%field%yvel0,      &
                    chunk%tiles(t)%field%yvel1,      &
                    chunk%tiles(t)%field%work_array1 )
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

  CALL clover_check_error(error_condition)
  IF(profiler_on) profiler%PdV=profiler%PdV+(timer()-kernel_time)

  IF(error_condition.EQ.1) THEN
    CALL report_error('PdV','error in PdV')
  ENDIF

  IF(predict)THEN
    CALL ideal_gas(.TRUE.)
    fields=0
    fields(FIELD_PRESSURE)=1
    CALL update_halo(fields,1)
  ENDIF

  IF ( predict ) THEN
    CALL revert()
  ENDIF

END SUBROUTINE PdV

END MODULE PdV_module
