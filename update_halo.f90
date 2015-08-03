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

!>  @brief Driver for the halo updates
!>  @author Wayne Gaudin
!>  @details Invokes the kernels for the internal and external halo cells for
!>  the fields specified.

MODULE update_halo_module

  USE definitions_module
  USE clover_module
  USE update_halo_kernel_module
  USE update_internal_halo_kernel_module

CONTAINS

SUBROUTINE update_halo(fields,depth)

  IMPLICIT NONE

  INTEGER :: fields(NUM_FIELDS),depth
  REAL(KIND=8) :: timer,halo_time

  IF (profiler_on) halo_time=timer()
  CALL clover_exchange(fields,depth)
  IF (profiler_on) profiler%halo_exchange = profiler%halo_exchange + (timer() - halo_time)

  CALL update_boundary(fields, depth)

  CALL update_tile_boundary(fields, depth)

END SUBROUTINE update_halo

SUBROUTINE update_boundary(fields,depth)

  IMPLICIT NONE

  INTEGER :: t,fields(NUM_FIELDS),depth
  REAL(KIND=8) :: timer,halo_time

  IF (profiler_on) halo_time=timer()

  IF (ANY(chunk%chunk_neighbours .EQ. EXTERNAL_FACE)) THEN
    IF (use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
      DO t=1,tiles_per_task
        CALL update_halo_kernel(chunk%tiles(t)%field%x_min,          &
                                chunk%tiles(t)%field%x_max,          &
                                chunk%tiles(t)%field%y_min,          &
                                chunk%tiles(t)%field%y_max,          &
                                chunk%chunk_neighbours,     &
                                chunk%tiles(t)%tile_neighbours,     &
                  chunk%tiles(t)%field%density0,                      &
                  chunk%tiles(t)%field%energy0,                       &
                  chunk%tiles(t)%field%pressure,                      &
                  chunk%tiles(t)%field%viscosity,                     &
                  chunk%tiles(t)%field%soundspeed,                    &
                  chunk%tiles(t)%field%density1,                      &
                  chunk%tiles(t)%field%energy1,                       &
                  chunk%tiles(t)%field%xvel0,                         &
                  chunk%tiles(t)%field%yvel0,                         &
                  chunk%tiles(t)%field%xvel1,                         &
                  chunk%tiles(t)%field%yvel1,                         &
                  chunk%tiles(t)%field%vol_flux_x,                    &
                  chunk%tiles(t)%field%vol_flux_y,                    &
                  chunk%tiles(t)%field%mass_flux_x,                   &
                  chunk%tiles(t)%field%mass_flux_y,                   &
                                fields,                         &
                                depth                           )
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF
  ENDIF

  IF (profiler_on) profiler%halo_update = profiler%halo_update + (timer() - halo_time)

END SUBROUTINE update_boundary

SUBROUTINE update_tile_boundary(fields, depth)

  IMPLICIT NONE

  INTEGER :: t,fields(NUM_FIELDS),depth, right_idx, up_idx
  REAL(KIND=8) :: timer,halo_time

  IF (profiler_on) halo_time=timer()

  IF (tiles_per_task .GT. 1) THEN
    IF (use_fortran_kernels)THEN
!$OMP PARALLEL PRIVATE(right_idx, up_idx)
!$OMP DO
      DO t=1,tiles_per_task
        right_idx = chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT)

        IF (right_idx .NE. EXTERNAL_FACE) THEN
          CALL update_internal_halo_kernel(                &
                                  chunk%tiles(t)%field%x_min,          &
                                  chunk%tiles(t)%field%x_max,          &
                                  chunk%tiles(t)%field%y_min,          &
                                  chunk%tiles(t)%field%y_max,          &
                  chunk%tiles(t)%field%density0,                      &
                  chunk%tiles(t)%field%energy0,                       &
                  chunk%tiles(t)%field%pressure,                      &
                  chunk%tiles(t)%field%viscosity,                     &
                  chunk%tiles(t)%field%soundspeed,                    &
                  chunk%tiles(t)%field%density1,                      &
                  chunk%tiles(t)%field%energy1,                       &
                  chunk%tiles(t)%field%xvel0,                         &
                  chunk%tiles(t)%field%yvel0,                         &
                  chunk%tiles(t)%field%xvel1,                         &
                  chunk%tiles(t)%field%yvel1,                         &
                  chunk%tiles(t)%field%vol_flux_x,                    &
                  chunk%tiles(t)%field%vol_flux_y,                    &
                  chunk%tiles(t)%field%mass_flux_x,                   &
                  chunk%tiles(t)%field%mass_flux_y,                   &
                                  chunk%tiles(right_idx)%field%x_min,          &
                                  chunk%tiles(right_idx)%field%x_max,          &
                                  chunk%tiles(right_idx)%field%y_min,          &
                                  chunk%tiles(right_idx)%field%y_max,          &
                  chunk%tiles(right_idx)%field%density0,                      &
                  chunk%tiles(right_idx)%field%energy0,                       &
                  chunk%tiles(right_idx)%field%pressure,                      &
                  chunk%tiles(right_idx)%field%viscosity,                     &
                  chunk%tiles(right_idx)%field%soundspeed,                    &
                  chunk%tiles(right_idx)%field%density1,                      &
                  chunk%tiles(right_idx)%field%energy1,                       &
                  chunk%tiles(right_idx)%field%xvel0,                         &
                  chunk%tiles(right_idx)%field%yvel0,                         &
                  chunk%tiles(right_idx)%field%xvel1,                         &
                  chunk%tiles(right_idx)%field%yvel1,                         &
                  chunk%tiles(right_idx)%field%vol_flux_x,                    &
                  chunk%tiles(right_idx)%field%vol_flux_y,                    &
                  chunk%tiles(right_idx)%field%mass_flux_x,                   &
                  chunk%tiles(right_idx)%field%mass_flux_y,                   &
                                  fields,                         &
                                  depth,                            &
                                  CHUNK_LEFT                           )
        ENDIF
      ENDDO
!$OMP END DO NOWAIT

!$OMP BARRIER

!$OMP DO
      DO t=1,tiles_per_task
        up_idx = chunk%tiles(t)%tile_neighbours(CHUNK_TOP)

        IF (up_idx .NE. EXTERNAL_FACE) THEN
          CALL update_internal_halo_kernel(                &
                                  chunk%tiles(t)%field%x_min,          &
                                  chunk%tiles(t)%field%x_max,          &
                                  chunk%tiles(t)%field%y_min,          &
                                  chunk%tiles(t)%field%y_max,          &
                  chunk%tiles(t)%field%density0,                      &
                  chunk%tiles(t)%field%energy0,                       &
                  chunk%tiles(t)%field%pressure,                      &
                  chunk%tiles(t)%field%viscosity,                     &
                  chunk%tiles(t)%field%soundspeed,                    &
                  chunk%tiles(t)%field%density1,                      &
                  chunk%tiles(t)%field%energy1,                       &
                  chunk%tiles(t)%field%xvel0,                         &
                  chunk%tiles(t)%field%yvel0,                         &
                  chunk%tiles(t)%field%xvel1,                         &
                  chunk%tiles(t)%field%yvel1,                         &
                  chunk%tiles(t)%field%vol_flux_x,                    &
                  chunk%tiles(t)%field%vol_flux_y,                    &
                  chunk%tiles(t)%field%mass_flux_x,                   &
                  chunk%tiles(t)%field%mass_flux_y,                   &
                                  chunk%tiles(up_idx)%field%x_min,          &
                                  chunk%tiles(up_idx)%field%x_max,          &
                                  chunk%tiles(up_idx)%field%y_min,          &
                                  chunk%tiles(up_idx)%field%y_max,          &
                  chunk%tiles(up_idx)%field%density0,                      &
                  chunk%tiles(up_idx)%field%energy0,                       &
                  chunk%tiles(up_idx)%field%pressure,                      &
                  chunk%tiles(up_idx)%field%viscosity,                     &
                  chunk%tiles(up_idx)%field%soundspeed,                    &
                  chunk%tiles(up_idx)%field%density1,                      &
                  chunk%tiles(up_idx)%field%energy1,                       &
                  chunk%tiles(up_idx)%field%xvel0,                         &
                  chunk%tiles(up_idx)%field%yvel0,                         &
                  chunk%tiles(up_idx)%field%xvel1,                         &
                  chunk%tiles(up_idx)%field%yvel1,                         &
                  chunk%tiles(up_idx)%field%vol_flux_x,                    &
                  chunk%tiles(up_idx)%field%vol_flux_y,                    &
                  chunk%tiles(up_idx)%field%mass_flux_x,                   &
                  chunk%tiles(up_idx)%field%mass_flux_y,                   &
                                  fields,                         &
                                  depth,                            &
                                  CHUNK_BOTTOM)
        ENDIF
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF
  ENDIF

  IF (profiler_on) profiler%internal_halo_update = profiler%internal_halo_update + (timer() - halo_time)

END SUBROUTINE update_tile_boundary

END MODULE update_halo_module

