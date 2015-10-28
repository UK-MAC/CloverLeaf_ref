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

!>  @brief standalone driver for the ideal gas kernels
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode


PROGRAM update_halo_driver

  USE set_data_module
  USE update_halo_kernel_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  INTEGER :: numargs,iargc,i
  CHARACTER (LEN=20)  :: command_line,temp

  INTEGER :: x_size,y_size

  REAL(KIND=8) :: kernel_time,timer,update_halo_time

  LOGICAL :: use_fortran_kernels,use_C_kernels
  INTEGER :: x_min,x_max,y_min,y_max,its,iteration
  REAL(KIND=8),ALLOCATABLE :: density0(:,:),  &
                              density1(:,:), &
                              energy0(:,:), &
                              energy1(:,:), &
                              pressure(:,:),  &
                              viscosity(:,:), &
                              soundspeed(:,:), &
                              xvel0(:,:), &
                              xvel1(:,:), &
                              yvel0(:,:), &
                              yvel1(:,:), &
                              vol_flux_x(:,:), &
                              mass_flux_x(:,:), &
                              vol_flux_y(:,:), &
                              mass_flux_y(:,:), &
                              xarea(:,:), &
                              yarea(:,:), &
                              celldy(:), celldx(:)
    REAL(KIND=8) :: dt

 INTEGER :: fields(15),depth
 INTEGER, DIMENSION(4) :: chunk_neighbours
 INTEGER, DIMENSION(4) :: tile_neighbours

    REAL(KIND=8) :: sum_density0, sum_density1, sum_energy0, sum_energy1, sum_pressure, sum_viscosity
    REAL(KIND=8) :: sum_soundspeed, sum_xvel0, sum_xvel1, sum_yvel0, sum_yvel1, sum_vol_flux_x, sum_mass_flux_x
    REAL(KIND=8) :: sum_vol_flux_y, sum_mass_flux_y

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

  WRITE(*,*) "Update Halo Kernel"
  WRITE(*,*) "Mesh size ",x_size,y_size
  WRITE(*,*) "Iterations ",its

  kernel_time=timer()

  CALL set_data(x_min,x_max,y_min,y_max, &
                density0=density0, &
                density1=density1,  &
                energy0=energy0,  &
                energy1=energy1, &
                viscosity=viscosity, &
                pressure=pressure, &
                soundspeed=soundspeed, &
                xvel0=xvel0, &
                xvel1=xvel1, &
                yvel0=yvel0, &
                yvel1=yvel1, &
                vol_flux_x=vol_flux_x, &
                vol_flux_y=vol_flux_y, &
                mass_flux_x=mass_flux_x, &
                mass_flux_y=mass_flux_y,    &
                xarea=xarea, &
                yarea=yarea, &
                celldx=celldx, &
                celldy=celldy, &
                dt=dt)

  WRITE(*,*) "Setup time(s) ",timer()-kernel_time

  WRITE(*,*) "Data initialised"
  WRITE(*,*) "Density0 Before: ",SUM(density0)
  WRITE(*,*) "Density1 Before: ",SUM(density1)
  WRITE(*,*) "Energy0 Before: ",SUM(energy0)
  WRITE(*,*) "Energy1 Before: ",SUM(energy1)
  WRITE(*,*) "Pressure Before: ",SUM(pressure)
  WRITE(*,*) "Viscosity Before: ",SUM(viscosity)
  WRITE(*,*) "Soundspeed Before: ",SUM(soundspeed)
  WRITE(*,*) "Xvel0 Before: ",SUM(xvel0)
  WRITE(*,*) "Xvel1 Before: ",SUM(xvel1)
  WRITE(*,*) "Yvel0 Before: ",SUM(yvel0)
  WRITE(*,*) "Yvel1 Before: ",SUM(yvel1)
  WRITE(*,*) "Vol_Flux_X Before: ",SUM(vol_flux_x)
  WRITE(*,*) "Mass_Flux_X Before: ",SUM(mass_flux_x)
  WRITE(*,*) "Vol_Flux_Y Before: ",SUM(vol_flux_y)
  WRITE(*,*) "Mass_Flux_Y Before: ",SUM(mass_flux_y)

    depth = 2
    fields = 1
    chunk_neighbours=-1
    tile_neighbours=-1

  IF(use_fortran_kernels) THEN
    WRITE(*,*) "Running Fortran kernel"
  ENDIF

  IF(use_C_kernels) THEN
    WRITE(*,*) "Running C kernel"
  ENDIF

  kernel_time=timer()

  IF(use_fortran_kernels)THEN
    DO iteration=1,its

        CALL update_halo_kernel(x_min,x_max,y_min,y_max,                            &
                              chunk_neighbours,                                           &
                              tile_neighbours,                                            &
                              density0,                                                   &
                              energy0,                                                    &
                              pressure,                                                   &
                              viscosity,                                                  &
                              soundspeed,                                                 &
                              density1,                                                   &
                              energy1,                                                    &
                              xvel0,                                                      &
                              yvel0,                                                      &
                              xvel1,                                                      &
                              yvel1,                                                      &
                              vol_flux_x,                                                 &
                              vol_flux_y,                                                 &
                              mass_flux_x,                                                &
                              mass_flux_y,                                                &
                              fields,                                                     &
                              depth                                                       )
    ENDDO
  ELSEIF(use_C_kernels)THEN
    DO iteration=1,its
        CALL update_halo_kernel_c(x_min,x_max,y_min,y_max,                            &
                              chunk_neighbours,                                           &
                              tile_neighbours,                                            &
                              density0,                                                   &
                              energy0,                                                    &
                              pressure,                                                   &
                              viscosity,                                                  &
                              soundspeed,                                                 &
                              density1,                                                   &
                              energy1,                                                    &
                              xvel0,                                                      &
                              yvel0,                                                      &
                              xvel1,                                                      &
                              yvel1,                                                      &
                              vol_flux_x,                                                 &
                              vol_flux_y,                                                 &
                              mass_flux_x,                                                &
                              mass_flux_y,                                                &
                              fields,                                                     &
                              depth                                                       )
    ENDDO
  ENDIF

  update_halo_time=(timer()-kernel_time)

    WRITE(*,*)
  WRITE(*,*) "Update Halo time(s) ",update_halo_time
    WRITE(*,*)

    sum_density0 = SUM(density0)
    sum_density1 = SUM(density1)
    sum_energy0 = SUM(energy0)
    sum_energy1 = SUM(energy1)
    sum_pressure = SUM(pressure)
    sum_viscosity = SUM(viscosity)
    sum_soundspeed = SUM(soundspeed)
    sum_xvel0 = SUM(xvel0)
    sum_xvel1 = SUM(xvel1)
    sum_yvel0 = SUM(yvel0)
    sum_yvel1 = SUM(yvel1)
    sum_vol_flux_x = SUM(vol_flux_x)
    sum_mass_flux_x = SUM(mass_flux_x)
    sum_vol_flux_y = SUM(vol_flux_y)
    sum_mass_flux_y = SUM(mass_flux_y)

    WRITE(*,*) "Density0 After: ", sum_density0
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Density0 Answer: 17088391.0284234"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Density1 After: ", sum_density1
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Density1 Answer: 17088391.0284234"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Energy0 After: ", sum_energy0
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Energy0 Answer: 17088391.0284234"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Energy1 After: ", sum_energy1
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Energy1 Answer: 17088391.0284234"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Pressure After: ", sum_pressure
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Pressure Answer: 8308062.15615501"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Viscosity After: ", sum_viscosity
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Viscosity Answer: 664.169227378228"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Soundspeed After: ", sum_soundspeed
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Soundspeed Answer: 11825376.6271950"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Xvel0 After: ", sum_xvel0
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Xvel0 Answer: 6537.05920469748"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Xvel1 After: ", sum_xvel1
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Xvel1 Answer: 6537.05920469748"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Yvel0 After: ", sum_yvel0
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Yvel0 Answer: -17066170.1154708"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Yvel1 After: ", sum_yvel1
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "Yvel1 Answer: -17066170.1154708"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Vol_Flux_X After: ",  sum_vol_flux_x
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "VolFluxX Answer: 2.932333442174252E-002"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Mass_Flux_X After: ", sum_mass_flux_x
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "MassFluxX Answer: 3.629576549212290E-002"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Vol_Flux_Y After: ",  sum_vol_flux_y
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "VolFluxY Answer: -76.5365042859362"
        WRITE(*,*)
    ENDIF

    WRITE(*,*) "Mass_Flux_Y After: ", sum_mass_flux_y
    IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840)) THEN
        WRITE(*,*) "MassFluxY Answer: 5.220883950012912E-002"
        WRITE(*,*)
    ENDIF


    !IF ((x_size .EQ. 3840) .AND. (y_size.EQ.3840))  THEN
    !    IF ((sum_density0.EQ.17088391.0284234) &
    !        .AND. (sum_density1.EQ.17088391.0284234) &
    !        .AND. (sum_energy0.EQ.17088391.0284234) &
    !        .AND. (sum_energy1.EQ.17088391.0284234) &
    !        .AND. (sum_pressure.EQ.8308062.15615501) &
    !        .AND. (sum_viscosity.EQ.664.169227378228) &
    !        .AND. (sum_soundspeed.EQ.11825376.6271950) &
    !        .AND. (sum_xvel0.EQ.6537.05920469748) &
    !        .AND. (sum_xvel1.EQ.6537.05920469748) &
    !        .AND. (sum_yvel0.EQ.-17066170.1154708) &
    !        .AND. (sum_yvel1.EQ.-17066170.1154708) &
    !        .AND. (sum_vol_flux_x.EQ.2.932333442174252E-002) &
    !        .AND. (sum_mass_flux_x.EQ.3.629576549212290E-002) &
    !        .AND. (sum_vol_flux_y.EQ.-76.5365042859362) &
    !        .AND. (sum_mass_flux_y.EQ.5.220883950012912E-002) &
    !       ) THEN
    !        WRITE(*,*) "TEST PASSED"
    !    ELSE
    !        WRITE(*,*) "TEST FAILED"
    !    ENDIF
    !ENDIF

    DEALLOCATE(density0)
    DEALLOCATE(energy0)
    DEALLOCATE(pressure)
    DEALLOCATE(viscosity)
    DEALLOCATE(soundspeed)
    DEALLOCATE(density1)
    DEALLOCATE(energy1)
    DEALLOCATE(xvel0)
    DEALLOCATE(yvel0)
    DEALLOCATE(xvel1)
    DEALLOCATE(yvel1)
    DEALLOCATE(vol_flux_x)
    DEALLOCATE(vol_flux_y)
    DEALLOCATE(mass_flux_x)
    DEALLOCATE(mass_flux_y)

END PROGRAM update_halo_driver

