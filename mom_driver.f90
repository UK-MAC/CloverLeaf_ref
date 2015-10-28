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

!>  @brief standalone driver for the momentum kernels
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode


PROGRAM mom_driver

  USE set_data_module
  USE advec_mom_kernel_mod

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  INTEGER :: numargs,iargc,i
  CHARACTER (LEN=20)  :: command_line,temp

  INTEGER :: x_size,y_size

  REAL(KIND=8) :: kernel_time,timer,mom_time

  LOGICAL :: use_fortran_kernels,use_C_kernels
  INTEGER :: x_min,x_max,y_min,y_max,its,iteration,direction,sweep_number,which_vel
  REAL(KIND=8) :: dt
  REAL(KIND=8),ALLOCATABLE :: celldx(:),celldy(:)
  REAL(KIND=8),ALLOCATABLE :: xarea(:,:),yarea(:,:)
  REAL(KIND=8),ALLOCATABLE :: volume(:,:)
  REAL(KIND=8),ALLOCATABLE :: density0(:,:),energy0(:,:),pressure(:,:),soundspeed(:,:),density1(:,:)
  REAL(KIND=8),ALLOCATABLE :: xvel0(:,:),yvel0(:,:),xvel1(:,:),yvel1(:,:)
  REAL(KIND=8),ALLOCATABLE :: vol_flux_x(:,:),vol_flux_y(:,:),mass_flux_x(:,:),mass_flux_y(:,:)
  REAL(KIND=8),ALLOCATABLE :: node_flux(:,:),node_mass_post(:,:),node_mass_pre(:,:)
  REAL(KIND=8),ALLOCATABLE :: mom_flux(:,:),pre_vol(:,:),post_vol(:,:)

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

  WRITE(*,*) "Advec Mom  Kernel"
  WRITE(*,*) "Mesh size ",x_size,y_size
  WRITE(*,*) "Iterations ",its

  kernel_time=timer()

  CALL set_data(x_min,x_max,y_min,y_max,   &
                celldx=celldx,             &
                celldy=celldy,             &
                xarea=xarea,               &
                yarea=yarea,               &
                volume=volume,             &
                density0=density0,         &
                energy0=energy0,           &
                density1=density1,         &
                pressure=pressure,         &
                soundspeed=soundspeed,     &
                xvel0=xvel0,               &
                xvel1=xvel1,               &
                yvel0=yvel0,               &
                yvel1=yvel1,               &
                vol_flux_x=vol_flux_x,     &
                vol_flux_y=vol_flux_y,     &
                mass_flux_x=mass_flux_x,   &
                mass_flux_y=mass_flux_y,   &
                work_array1=node_flux,     &
                work_array2=node_mass_post,&
                work_array3=node_mass_pre, &
                work_array4=mom_flux,      &
                work_array5=pre_vol,       &
                work_array6=post_vol,      &
                dt=dt                      )

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

      direction=1
      which_vel=1
      sweep_number=1
      CALL advec_mom_kernel(x_min,x_max,y_min,y_max,   &
                            xvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
      direction=1
      which_vel=2
      sweep_number=1
      CALL advec_mom_kernel(x_min,x_max,y_min,y_max,   &
                            yvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
      direction=2
      which_vel=1
      sweep_number=2
      CALL advec_mom_kernel(x_min,x_max,y_min,y_max,   &
                            xvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
      direction=2
      which_vel=2
      sweep_number=2
      CALL advec_mom_kernel(x_min,x_max,y_min,y_max,   &
                            yvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
    ENDDO
  ELSEIF(use_C_kernels)THEN
    DO iteration=1,its
      direction=1
      which_vel=1
      sweep_number=1
      CALL advec_mom_kernel_c(x_min,x_max,y_min,y_max, &
                            xvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
      direction=1
      which_vel=2
      sweep_number=1
      CALL advec_mom_kernel_c(x_min,x_max,y_min,y_max, &
                            yvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
      direction=2
      which_vel=1
      sweep_number=2
      CALL advec_mom_kernel_c(x_min,x_max,y_min,y_max, &
                            xvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
      direction=2
      which_vel=2
      sweep_number=2
      CALL advec_mom_kernel_c(x_min,x_max,y_min,y_max, &
                            yvel1,                     &
                            mass_flux_x,               &
                            vol_flux_x,                &
                            mass_flux_y,               &
                            vol_flux_y,                &
                            volume,                    &
                            density1,                  &
                            node_flux,                 &
                            node_mass_post,            &
                            node_mass_pre,             &
                            mom_flux,                  &
                            pre_vol,                   &
                            post_vol,                  &
                            celldx,                    &
                            celldy,                    &
                            which_vel,                 &
                            sweep_number,              &
                            direction                  )
    ENDDO
  ENDIF

  mom_time=(timer()-kernel_time)

  WRITE(*,*) "Advec mom time ",mom_time 
  WRITE(*,*) "X velocity ",SUM(xvel1)
  WRITE(*,*) "Y velocity ",SUM(yvel1)

  DEALLOCATE(celldx)
  DEALLOCATE(celldy)
  DEALLOCATE(xarea)
  DEALLOCATE(yarea)
  DEALLOCATE(volume)
  DEALLOCATE(density0)
  DEALLOCATE(energy0)
  DEALLOCATE(soundspeed)
  DEALLOCATE(density1)
  DEALLOCATE(xvel0)
  DEALLOCATE(xvel1)
  DEALLOCATE(yvel0)
  DEALLOCATE(yvel1)
  DEALLOCATE(vol_flux_x)
  DEALLOCATE(vol_flux_y)
  DEALLOCATE(mass_flux_x)
  DEALLOCATE(mass_flux_y)
  DEALLOCATE(node_flux)
  DEALLOCATE(node_mass_post)
  DEALLOCATE(node_mass_pre)
  DEALLOCATE(mom_flux)
  DEALLOCATE(pre_vol)
  DEALLOCATE(post_vol)

END PROGRAM mom_driver

