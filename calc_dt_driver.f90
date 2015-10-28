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

!>  @brief standalone driver for the calc_dt kernels
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode


PROGRAM calc_dt_driver

  USE set_data_module
  USE calc_dt_kernel_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  INTEGER :: numargs,iargc,i
  CHARACTER (LEN=20)  :: command_line,temp

  INTEGER :: x_size,y_size

  REAL(KIND=8) :: kernel_time,timer,calc_dt_time

  LOGICAL :: use_fortran_kernels,use_C_kernels
  INTEGER :: x_min,x_max,y_min,y_max,its,iteration
  REAL(KIND=8) :: dt,g_small,g_big,dtmin,dtc_safe,dtu_safe,dtv_safe,dtdiv_safe,xl_pos,yl_pos,local_dt
  INTEGER :: l_control,jldt,kldt,small
  REAL(KIND=8),ALLOCATABLE :: xarea(:,:),yarea(:,:),volume(:,:)
  REAL(KIND=8),ALLOCATABLE :: cellx(:),celly(:)
  REAL(KIND=8),ALLOCATABLE :: celldx(:),celldy(:)
  REAL(KIND=8),ALLOCATABLE :: density0(:,:),energy0(:,:)
  REAL(KIND=8),ALLOCATABLE :: pressure(:,:),soundspeed(:,:),viscosity(:,:)
  REAL(KIND=8),ALLOCATABLE :: xvel0(:,:),yvel0(:,:),work_array1(:,:)

!$OMP PARALLEL
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
!$    WRITE(*,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
!$  ENDIF
!$OMP END PARALLEL

  x_size=100
  y_size=100
  its=1
  use_fortran_kernels=.TRUE.
  use_C_kernels=.FALSE.
  g_small    =1.0e-16
  g_big      =1.0e+21
  dtmin      =0.0000001_8
  dtc_safe   =0.7_8
  dtu_safe   =0.5_8
  dtv_safe   =0.5_8
  dtdiv_safe =0.7_8

  numargs = iargc()

  DO i=1,numargs,2
    CALL GETARG(i,command_line)
    SELECT CASE (command_line)
      CASE("-help")
        WRITE(*,*) "Usage -nx 100 -ny 100 -its 10 -kernel fortran|c"
        stop
      CASE("-nx")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") x_size
      CASE("-ny")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") y_size
      CASE("-its")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") its
      CASE("-kernel")
        CALL GETARG(i+1,temp)
        IF(temp.EQ."fortran") THEN
          use_fortran_kernels=.TRUE.
          use_C_kernels=.FALSE.
        ENDIF
        IF(temp.EQ."c") THEN
          use_fortran_kernels=.FALSE.
          use_C_kernels=.TRUE.
        ENDIF
    END SELECT
  ENDDO

  x_min=1
  y_min=1
  x_max=x_size
  y_max=y_size

  WRITE(*,*) "Calc dt Kernel"
  WRITE(*,*) "Mesh size ",x_size,y_size
  WRITE(*,*) "Iterations ",its

  kernel_time=timer()

  CALL set_data(x_min,x_max,y_min,y_max, &
                xarea=xarea,             &
                yarea=yarea,             &
                cellx=cellx,             &
                celly=celly,             &
                celldx=celldx,           &
                celldy=celldy,           &
                volume=volume,           &
                density0=density0,       &
                energy0=energy0,         &
                pressure=pressure,       &
                soundspeed=soundspeed,   &
                viscosity=viscosity,     &
                xvel0=xvel0,             &
                yvel0=yvel0,             &
                work_array1=work_array1, &
                dt=dt                    )

  WRITE(*,*) "Setup time ",timer()-kernel_time

  WRITE(*,*) "Data initialised"

  IF(use_fortran_kernels) THEN
    WRITE(*,*) "Running Fortran kernel"
  ENDIF

  IF(use_C_kernels) THEN
    WRITE(*,*) "Running C kernel"
  ENDIF

  kernel_time=timer()

  IF(use_fortran_kernels)THEN
    DO iteration=1,its
      CALL calc_dt_kernel(x_min,      &
                          x_max,      &
                          y_min,      &
                          y_max,      &
                          g_small,    &
                          g_big,      &
                          dtmin,      &
                          dtc_safe,   &
                          dtu_safe,   &
                          dtv_safe,   &
                          dtdiv_safe, &
                          xarea,      &
                          yarea,      &
                          cellx,      &
                          celly,      &
                          celldx,     &
                          celldy,     &
                          volume,     &
                          density0,   &
                          energy0,    &
                          pressure,   &
                          viscosity,  &
                          soundspeed, &
                          xvel0,      &
                          yvel0,      &
                          work_array1,&
                          local_dt,   &
                          l_control,  &
                          xl_pos,     &
                          yl_pos,     &
                          jldt,       &
                          kldt,       &
                          small       )
    ENDDO
  ELSEIF(use_C_kernels)THEN
    DO iteration=1,its
      CALL calc_dt_kernel_c(x_min,    &
                          x_max,      &
                          y_min,      &
                          y_max,      &
                          g_small,    &
                          g_big,      &
                          dtmin,      &
                          dtc_safe,   &
                          dtu_safe,   &
                          dtv_safe,   &
                          dtdiv_safe, &
                          xarea,      &
                          yarea,      &
                          cellx,      &
                          celly,      &
                          celldx,     &
                          celldy,     &
                          volume,     &
                          density0,   &
                          energy0,    &
                          pressure,   &
                          viscosity,  &
                          soundspeed, &
                          xvel0,      &
                          yvel0,      &
                          work_array1,&
                          local_dt,   &
                          l_control,  &
                          xl_pos,     &
                          yl_pos,     &
                          jldt,       &
                          kldt,       &
                          small       )
    ENDDO
  ENDIF

  calc_dt_time=(timer()-kernel_time)

  WRITE(*,*) "Calc dt time ",calc_dt_time 
  WRITE(*,*) "local_dt ",local_dt

  DEALLOCATE(xarea)
  DEALLOCATE(yarea)
  DEALLOCATE(cellx)
  DEALLOCATE(celly)
  DEALLOCATE(celldx)
  DEALLOCATE(celldy)
  DEALLOCATE(volume)
  DEALLOCATE(density0)
  DEALLOCATE(energy0)
  DEALLOCATE(pressure)
  DEALLOCATE(soundspeed)
  DEALLOCATE(viscosity)
  DEALLOCATE(xvel0)
  DEALLOCATE(yvel0)
  DEALLOCATE(work_array1)

END PROGRAM calc_dt_driver

