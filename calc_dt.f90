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

!>  @brief Driver for the timestep kernels
!>  @author Wayne Gaudin
!>  @details Invokes the user specified timestep kernel.

MODULE calc_dt_module

CONTAINS

SUBROUTINE calc_dt(local_dt,local_control,x_pos,y_pos)

  USE definitions_module
  USE calc_dt_kernel_module

  IMPLICIT NONE

  INTEGER          :: t
  REAL(KIND=8)     :: local_dt,x_pos,y_pos
  CHARACTER(LEN=8) :: local_control

  INTEGER          :: jldt,kldt
  REAL(KIND=8)     :: xl_pos,yl_pos,dtlp

  INTEGER          :: l_control
  INTEGER          :: small

  local_dt=g_big

  small = 0

  IF(use_fortran_kernels)THEN
!$OMP PARALLEL PRIVATE(dtlp)
!$OMP DO
    DO t=1,tiles_per_task
      CALL calc_dt_kernel(chunk%tiles(t)%field%x_min,     &
                          chunk%tiles(t)%field%x_max,     &
                          chunk%tiles(t)%field%y_min,     &
                          chunk%tiles(t)%field%y_max,     &
                          g_small,                       &
                          g_big,                         &
                          dtmin,                         &
                          dtc_safe,                      &
                          dtu_safe,                      &
                          dtv_safe,                      &
                          dtdiv_safe,                    &
                          chunk%tiles(t)%field%xarea,     &
                          chunk%tiles(t)%field%yarea,     &
                          chunk%tiles(t)%field%cellx,     &
                          chunk%tiles(t)%field%celly,     &
                          chunk%tiles(t)%field%celldx,    &
                          chunk%tiles(t)%field%celldy,    &
                          chunk%tiles(t)%field%volume,    &
                          chunk%tiles(t)%field%density0,  &
                          chunk%tiles(t)%field%energy0,   &
                          chunk%tiles(t)%field%pressure,  &
                          chunk%tiles(t)%field%viscosity, &
                          chunk%tiles(t)%field%soundspeed,&
                          chunk%tiles(t)%field%xvel0,     &
                          chunk%tiles(t)%field%yvel0,     &
                          chunk%tiles(t)%field%work_array1,&
                          local_dt,                      &
                          l_control,                     &
                          xl_pos,                        &
                          yl_pos,                        &
                          jldt,                          &
                          kldt,                          &
                          small                          )

!$OMP CRITICAL
    IF(dtlp.LE.dt) THEN
      dt=dtlp
      x_pos=xl_pos
      y_pos=yl_pos
      jdt=jldt
      kdt=kldt
    ENDIF
!$OMP END CRITICAL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF 

  IF(l_control.EQ.1) local_control='sound'
  IF(l_control.EQ.2) local_control='xvel'
  IF(l_control.EQ.3) local_control='yvel'
  IF(l_control.EQ.4) local_control='div'

END SUBROUTINE calc_dt

END MODULE calc_dt_module
