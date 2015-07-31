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

!>  @brief Driver for the viscosity kernels
!>  @author Wayne Gaudin
!>  @details Selects the user specified kernel to caluclate the artificial 
!>  viscosity.

MODULE viscosity_module

CONTAINS

SUBROUTINE viscosity()

  USE clover_module
  USE viscosity_kernel_module
  
  IMPLICIT NONE

  INTEGER :: t

  IF(use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL viscosity_kernel(chunk%tiles(t)%field%x_min,                   &
                          chunk%tiles(t)%field%x_max,                     &
                          chunk%tiles(t)%field%y_min,                     &
                          chunk%tiles(t)%field%y_max,                     &
                          chunk%tiles(t)%field%celldx,                    &
                          chunk%tiles(t)%field%celldy,                    &
                          chunk%tiles(t)%field%density0,                  &
                          chunk%tiles(t)%field%pressure,                  &
                          chunk%tiles(t)%field%viscosity,                 &
                          chunk%tiles(t)%field%xvel0,                     &
                          chunk%tiles(t)%field%yvel0                      )
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE viscosity

END MODULE viscosity_module
