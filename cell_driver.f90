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

!>  @brief standalone driver for the cell advection kernels
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode


PROGRAM cell_driver

  USE set_data_module
  USE advec_cell_kernel_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  INTEGER :: numargs,iargc,i
  CHARACTER (LEN=20)  :: command_line,temp

  INTEGER :: x_size,y_size

  REAL(KIND=8) :: kernel_time,timer,cell_time

  LOGICAL :: use_fortran_kernels,use_C_kernels,reset_data
  INTEGER :: x_min,x_max,y_min,y_max,its,iteration,direction,sweep_number
  REAL(KIND=8) :: dt
  REAL(KIND=8),ALLOCATABLE :: celldx(:),celldy(:)
  REAL(KIND=8),ALLOCATABLE :: vertexdx(:),vertexdy(:)
  REAL(KIND=8),ALLOCATABLE :: xarea(:,:),yarea(:,:)

  REAL(KIND=8),ALLOCATABLE :: volume(:,:)
  REAL(KIND=8),ALLOCATABLE :: density1(:,:),energy1(:,:), density0(:,:), pressure(:,:), soundspeed(:,:), energy0(:,:)

  REAL(KIND=8),ALLOCATABLE :: xvel0(:,:),yvel0(:,:),xvel1(:,:),yvel1(:,:)
  REAL(KIND=8),ALLOCATABLE :: vol_flux_x(:,:),vol_flux_y(:,:),mass_flux_x(:,:),mass_flux_y(:,:)

  REAL(KIND=8),ALLOCATABLE :: post_ener(:,:),post_mass(:,:),pre_mass(:,:)
  REAL(KIND=8),ALLOCATABLE :: advec_vol(:,:),ener_flux(:,:),pre_vol(:,:),post_vol(:,:)
  REAL(KIND=8),ALLOCATABLE :: density_orig(:,:),energy_orig(:,:)

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

  numargs = iargc()

  DO i=1,numargs,2
    CALL GETARG(i,command_line)
    SELECT CASE (command_line)
      CASE("-help")
        WRITE(*,*) "Usage -nx 100 -ny 100 -its 10 -kernel fortran|c -reset off|on"
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
      CASE("-reset")
        CALL GETARG(i+1,temp)
        IF(temp.EQ."on") THEN
          reset_data=.TRUE.
        ENDIF
        IF(temp.EQ."off") THEN
          reset_data=.FALSE.
        ENDIF
    END SELECT
  ENDDO

  x_min=1
  y_min=1
  x_max=x_size
  y_max=y_size

  WRITE(*,*) "Advec Cell Kernel"
  WRITE(*,*) "Mesh size ",x_size,y_size
  WRITE(*,*) "Iterations ",its

  kernel_time=timer()

  CALL set_data(x_min,x_max,y_min,y_max,   &
                vertexdx=vertexdx,         &
                vertexdy=vertexdy,         &
                volume=volume,             &
                density1=density1,         &
                energy1=energy1,           &
                vol_flux_x=vol_flux_x,     &
                vol_flux_y=vol_flux_y,     &
                mass_flux_x=mass_flux_x,   &
                mass_flux_y=mass_flux_y,   &
                work_array1=pre_vol,       &
                work_array2=post_vol,       &
                work_array3=pre_mass,       &
                work_array4=post_mass,      &
                work_array5=advec_vol,      &
                work_array6=post_ener,      &
                work_array7=ener_flux,      &
                xarea=xarea,                &
                yarea=yarea,                &
                xvel0=xvel0,                &
                xvel1=xvel1,                &
                yvel0=yvel0,                &
                yvel1=yvel1,                &
                dt=dt,                      &
                soundspeed=soundspeed,      &
                density0=density0,          &
                pressure=pressure,          &
                energy0=energy0,            &
                celldx=celldx,              &
                celldy=celldy             )

  IF(reset_data) THEN
    ALLOCATE(density_orig(x_min-2:x_max+2,y_min-2:y_max+2))
    ALLOCATE(energy_orig(x_min-2:x_max+2,y_min-2:y_max+2))
    density_orig=density1
    energy_orig=energy1
  ENDIF

  WRITE(*,*) "Setup time ",timer()-kernel_time

  WRITE(*,*) "Data initialised"

  IF(use_fortran_kernels) THEN
    WRITE(*,*) "Running Fortran kernel"
  ENDIF

  IF(use_C_kernels) THEN
    WRITE(*,*) "Running C kernel"
  ENDIF

  IF(reset_data) THEN
    WRITE(*,*) "Resetting data for each iteration"
  ELSE
    WRITE(*,*) "Not resetting data for each iteration"
  ENDIF

  WRITE(*,*)
  WRITE(*,*) "Density Before: ",SUM(density1)
  WRITE(*,*) "Energy Before: ",SUM(energy1)
  WRITE(*,*) 

  cell_time=0.0

  IF(use_fortran_kernels)THEN
    DO iteration=1,its
      kernel_time=timer()
      direction=1
      sweep_number=1
      CALL advec_cell_kernel(x_min,x_max,y_min,y_max,   &
                             direction,                       &
                             sweep_number,              &
                             vertexdx,                  &
                             vertexdy,                  &
                             volume,                    &
                             density1,                  &
                             energy1,                   &
                             mass_flux_x,               &
                             vol_flux_x,                &
                             mass_flux_y,               &
                             vol_flux_y,                &
                             pre_vol,                   &
                             post_vol,                  &
                             pre_mass,                  &
                             post_mass,                 &
                             advec_vol,                 &
                             post_ener,                 &
                             ener_flux                  )
      direction=2
      sweep_number=2
      CALL advec_cell_kernel(x_min,x_max,y_min,y_max,   &
                             direction,                       &
                             sweep_number,              &
                             vertexdx,                  &
                             vertexdy,                  &
                             volume,                    &
                             density1,                  &
                             energy1,                   &
                             mass_flux_x,               &
                             vol_flux_x,                &
                             mass_flux_y,               &
                             vol_flux_y,                &
                             pre_vol,                   &
                             post_vol,                  &
                             pre_mass,                  &
                             post_mass,                 &
                             advec_vol,                 &
                             post_ener,                 &
                             ener_flux                  )
      cell_time=cell_time+(timer()-kernel_time)
      IF(reset_data.AND.iteration.NE.its) THEN
        density1=density_orig
        energy1=energy_orig
      ENDIF
    ENDDO
  ELSEIF(use_C_kernels)THEN
    DO iteration=1,its
      kernel_time=timer()
      direction=1
      sweep_number=1
      CALL advec_cell_kernel_c(x_min,x_max,y_min,y_max, &
                               direction,               &
                               sweep_number,            &
                               vertexdx,                &
                               vertexdy,                &
                               volume,                  &
                               density1,                &
                               energy1,                 &
                               mass_flux_x,             &
                               vol_flux_x,              &
                               mass_flux_y,             &
                               vol_flux_y,              &
                               pre_vol,                   &
                               post_vol,                  &
                               pre_mass,                  &
                               post_mass,                 &
                               advec_vol,                 &
                               post_ener,                 &
                               ener_flux                  )
      direction=2
      sweep_number=2
      CALL advec_cell_kernel_c(x_min,x_max,y_min,y_max, &
                               direction,                     &
                               sweep_number,            &
                               vertexdx,                &
                               vertexdy,                &
                               volume,                  &
                               density1,                &
                               energy1,                 &
                               mass_flux_x,             &
                               vol_flux_x,              &
                               mass_flux_y,             &
                               vol_flux_y,              &
                               pre_vol,                   &
                               post_vol,                  &
                               pre_mass,                  &
                               post_mass,                 &
                               advec_vol,                 &
                               post_ener,                 &
                               ener_flux                  )
      cell_time=cell_time+(timer()-kernel_time)
      IF(reset_data.AND.iteration.NE.its) THEN
        density1=density_orig
        energy1=energy_orig
      ENDIF
    ENDDO
  ENDIF

  WRITE(*,*) "Advec cell time ",cell_time 
  WRITE(*,*)
  WRITE(*,*) "Density ",SUM(density1)
  WRITE(*,*) "Energy ",SUM(energy1)

  DEALLOCATE(vertexdx)
  DEALLOCATE(vertexdy)
  DEALLOCATE(volume)
  DEALLOCATE(density1)
  DEALLOCATE(energy1)
  DEALLOCATE(vol_flux_x)
  DEALLOCATE(vol_flux_y)
  DEALLOCATE(mass_flux_x)
  DEALLOCATE(mass_flux_y)
  DEALLOCATE(pre_vol)
  DEALLOCATE(post_vol)
  DEALLOCATE(pre_mass)
  DEALLOCATE(post_mass)
  DEALLOCATE(advec_vol)
  DEALLOCATE(post_ener)
  DEALLOCATE(ener_flux)
  DEALLOCATE(xarea)
  DEALLOCATE(yarea)
  DEALLOCATE(xvel0)
  DEALLOCATE(xvel1)
  DEALLOCATE(yvel0)
  DEALLOCATE(yvel1)
  DEALLOCATE(soundspeed)
  DEALLOCATE(density0)
  DEALLOCATE(pressure)
  DEALLOCATE(energy0)
  DEALLOCATE(celldx)
  DEALLOCATE(celldy)
  IF(reset_data) THEN
    DEALLOCATE(density_orig)
    DEALLOCATE(energy_orig)
  ENDIF

END PROGRAM cell_driver

