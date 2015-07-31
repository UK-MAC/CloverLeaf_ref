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

!>  @brief Ideal gas kernel driver
!>  @author Wayne Gaudin
!>  @details Invokes the user specified kernel for the ideal gas equation of
!>  state using the specified time level data.

MODULE ideal_gas_module

CONTAINS

SUBROUTINE ideal_gas(predict)

  USE clover_module
  USE ideal_gas_kernel_module

  IMPLICIT NONE

  INTEGER :: t

  LOGICAL :: predict
  REAL(KIND=8) :: kernel_time, timer

  IF(profiler_on) kernel_time=timer()

  IF (.NOT. predict) THEN
    IF(use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
      DO t=1,tiles_per_task
        CALL ideal_gas_kernel(chunk%tiles(t)%field%x_min,    &
                          chunk%tiles(t)%field%x_max,      &
                          chunk%tiles(t)%field%y_min,      &
                          chunk%tiles(t)%field%y_max,      &
                          chunk%tiles(t)%field%density0,   &
                          chunk%tiles(t)%field%energy0,    &
                          chunk%tiles(t)%field%pressure,   &
                          chunk%tiles(t)%field%soundspeed  )
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF
  ELSE
    IF(use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
      DO t=1,tiles_per_task
        CALL ideal_gas_kernel(chunk%tiles(t)%field%x_min,    &
                          chunk%tiles(t)%field%x_max,      &
                          chunk%tiles(t)%field%y_min,      &
                          chunk%tiles(t)%field%y_max,      &
                          chunk%tiles(t)%field%density1,   &
                          chunk%tiles(t)%field%energy1,    &
                          chunk%tiles(t)%field%pressure,   &
                          chunk%tiles(t)%field%soundspeed  )
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF
  ENDIF

  IF(profiler_on) profiler%ideal_gas=profiler%ideal_gas+(timer()-kernel_time)

END SUBROUTINE ideal_gas

END MODULE ideal_gas_module

