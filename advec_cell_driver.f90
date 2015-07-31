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

!>  @brief Cell centred advection driver.
!>  @author Wayne Gaudin
!>  @details Invokes the user selected advection kernel.

MODULE  advec_cell_driver_module

CONTAINS

SUBROUTINE advec_cell_driver(sweep_number,dir)

  USE clover_module
  USE advec_cell_kernel_module

  IMPLICIT NONE

  INTEGER :: t,sweep_number,dir

  IF(use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
    DO t=1,tiles_per_task
      CALL advec_cell_kernel(chunk%tiles(t)%field%x_min,               &
                           chunk%tiles(t)%field%x_max,                 &
                           chunk%tiles(t)%field%y_min,                 &
                           chunk%tiles(t)%field%y_max,                 &
                           dir,                                       &
                           sweep_number,                              &
                           chunk%tiles(t)%field%vertexdx,              &
                           chunk%tiles(t)%field%vertexdy,              &
                           chunk%tiles(t)%field%volume,                &
                           chunk%tiles(t)%field%density1,              &
                           chunk%tiles(t)%field%energy1,               &
                           chunk%tiles(t)%field%mass_flux_x,           &
                           chunk%tiles(t)%field%vol_flux_x,            &
                           chunk%tiles(t)%field%mass_flux_y,           &
                           chunk%tiles(t)%field%vol_flux_y,            &
                           chunk%tiles(t)%field%work_array1,           &
                           chunk%tiles(t)%field%work_array2,           &
                           chunk%tiles(t)%field%work_array3,           &
                           chunk%tiles(t)%field%work_array4,           &
                           chunk%tiles(t)%field%work_array5,           &
                           chunk%tiles(t)%field%work_array6,           &
                           chunk%tiles(t)%field%work_array7            )
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE advec_cell_driver

END MODULE  advec_cell_driver_module

