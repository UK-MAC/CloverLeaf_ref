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

!>  @brief Calculate the minimum timestep for all mesh chunks.
!>  @author Wayne Gaudin
!>  @details Invokes the kernels needed to calculate the timestep and finds
!>  the minimum across all chunks. Checks if the timestep falls below the
!>  user specified limitand outputs the timestep information.

MODULE timestep_module

CONTAINS

SUBROUTINE timestep()

  USE clover_module
  USE report_module
  USE update_halo_module
  USE viscosity_module
  USE calc_dt_module
  USE ideal_gas_module
  USE definitions_module

  IMPLICIT NONE

  INTEGER :: t

  REAL(KIND=8)    :: x_pos,y_pos

  REAL(KIND=8)    :: kernel_time,timer

  CHARACTER(LEN=8) :: dt_control

  INTEGER :: small

  INTEGER :: fields(NUM_FIELDS)

!$ INTEGER :: OMP_GET_THREAD_NUM

  dt    = g_big
  small=0

  CALL ideal_gas(.FALSE.)

  fields=0
  fields(FIELD_PRESSURE)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_DENSITY0)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  CALL update_halo(fields,1)

  CALL viscosity()

  fields=0
  fields(FIELD_VISCOSITY)=1
  CALL update_halo(fields,1)

  IF(profiler_on) kernel_time=timer()
  CALL calc_dt(dt_control,x_pos,y_pos)

  dt = MIN(dt, (dtold * dtrise), dtmax)

  CALL clover_min(dt)
  IF(profiler_on) profiler%timestep=profiler%timestep+(timer()-kernel_time)

  IF(dt.LT.dtmin) small=1

  IF (parallel%boss) THEN
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,"(' Step ', i7,' time ', f11.7,' control ',a11,' timestep  ',1pe9.2,i8,',',i8,' x ',1pe9.2,' y ',1pe9.2)") &
                      step,time,dt_control,dt,jdt,kdt,x_pos,y_pos
      WRITE(0,"(' Step ', i7,' time ', f11.7,' control ',a11,' timestep  ',1pe9.2,i8,',',i8,' x ',1pe9.2,' y ',1pe9.2)") &
                      step,time,dt_control,dt,jdt,kdt,x_pos,y_pos
!$  ENDIF
  ENDIF

  IF(small.EQ.1) THEN
    CALL report_error('timestep','small timestep')
  ENDIF

  dtold = dt

END SUBROUTINE timestep

END MODULE timestep_module
