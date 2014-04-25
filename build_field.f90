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

!>  @brief  Allocates the data for each mesh chunk
!>  @author Wayne Gaudin
!>  @details The data fields for the mesh chunk are allocated based on the mesh
!>  size.

SUBROUTINE build_field(chunk,x_cells,y_cells)

   USE clover_module

   IMPLICIT NONE

   INTEGER :: chunk,x_cells,y_cells, j, k

   chunks(chunk)%field%x_min=1
   chunks(chunk)%field%y_min=1

   chunks(chunk)%field%x_max=x_cells
   chunks(chunk)%field%y_max=y_cells

   ALLOCATE(chunks(chunk)%field%density0  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%density1  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%energy0   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%energy1   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%pressure  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%viscosity (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%soundspeed(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))

   ALLOCATE(chunks(chunk)%field%xvel0(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%xvel1(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%yvel0(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%yvel1(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%vol_flux_x (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%mass_flux_x(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%vol_flux_y (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%mass_flux_y(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%work_array1(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array2(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array3(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array4(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array5(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array6(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array7(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%cellx   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2))
   ALLOCATE(chunks(chunk)%field%celly   (chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%vertexx (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3))
   ALLOCATE(chunks(chunk)%field%vertexy (chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%celldx  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2))
   ALLOCATE(chunks(chunk)%field%celldy  (chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%vertexdx(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3))
   ALLOCATE(chunks(chunk)%field%vertexdy(chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%volume  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                         chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%xarea   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                         chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%yarea   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                         chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ! Zeroing isn't strictly neccessary but it ensures physical pages
   ! are allocated. This prevents first touch overheads in the main code
   ! cycle which can skew timings in the first step

!$OMP PARALLEL
!$OMP DO 
    DO k=chunks(chunk)%field%y_min-2,chunks(chunk)%field%y_max+3
        DO j=chunks(chunk)%field%x_min-2,chunks(chunk)%field%x_max+3
            chunks(chunk)%field%work_array1(j,k)=0.0
            chunks(chunk)%field%work_array2(j,k)=0.0
            chunks(chunk)%field%work_array3(j,k)=0.0
            chunks(chunk)%field%work_array4(j,k)=0.0
            chunks(chunk)%field%work_array5(j,k)=0.0
            chunks(chunk)%field%work_array6(j,k)=0.0
            chunks(chunk)%field%work_array7(j,k)=0.0

            chunks(chunk)%field%xvel0(j,k)=0.0
            chunks(chunk)%field%xvel1(j,k)=0.0
            chunks(chunk)%field%yvel0(j,k)=0.0
            chunks(chunk)%field%yvel1(j,k)=0.0
        ENDDO
    ENDDO
!$OMP END DO

!$OMP DO 
    DO k=chunks(chunk)%field%y_min-2,chunks(chunk)%field%y_max+2
        DO j=chunks(chunk)%field%x_min-2,chunks(chunk)%field%x_max+2
            chunks(chunk)%field%density0(j,k)=0.0
            chunks(chunk)%field%density1(j,k)=0.0
            chunks(chunk)%field%energy0(j,k)=0.0
            chunks(chunk)%field%energy1(j,k)=0.0
            chunks(chunk)%field%pressure(j,k)=0.0
            chunks(chunk)%field%viscosity(j,k)=0.0
            chunks(chunk)%field%soundspeed(j,k)=0.0
            chunks(chunk)%field%volume(j,k)=0.0
        ENDDO
    ENDDO
!$OMP END DO

!$OMP DO 
    DO k=chunks(chunk)%field%y_min-2,chunks(chunk)%field%y_max+2  
        DO j=chunks(chunk)%field%x_min-2,chunks(chunk)%field%x_max+3  
            chunks(chunk)%field%vol_flux_x(j,k)=0.0
            chunks(chunk)%field%mass_flux_x(j,k)=0.0
            chunks(chunk)%field%xarea(j,k)=0.0
        ENDDO
    ENDDO
!$OMP END DO
!$OMP DO 
    DO k=chunks(chunk)%field%y_min-2,chunks(chunk)%field%y_max+3
        DO j=chunks(chunk)%field%x_min-2,chunks(chunk)%field%x_max+2
            chunks(chunk)%field%vol_flux_y(j,k)=0.0
            chunks(chunk)%field%mass_flux_y(j,k)=0.0
            chunks(chunk)%field%yarea(j,k)=0.0
        ENDDO
    ENDDO
!$OMP END DO

!$OMP DO 
    DO j=chunks(chunk)%field%x_min-2,chunks(chunk)%field%x_max+2
        chunks(chunk)%field%cellx(j)=0.0
        chunks(chunk)%field%celldx(j)=0.0
    ENDDO
!$OMP END DO
!$OMP DO 
    DO k=chunks(chunk)%field%y_min-2,chunks(chunk)%field%y_max+2
        chunks(chunk)%field%celly(k)=0.0
        chunks(chunk)%field%celldy(k)=0.0
    ENDDO
!$OMP END DO

!$OMP DO 
    DO j=chunks(chunk)%field%x_min-2,chunks(chunk)%field%x_max+3
        chunks(chunk)%field%vertexx(j)=0.0
        chunks(chunk)%field%vertexdx(j)=0.0
    ENDDO
!$OMP END DO
!$OMP DO 
    DO k=chunks(chunk)%field%y_min-2,chunks(chunk)%field%y_max+3
        chunks(chunk)%field%vertexy(k)=0.0
        chunks(chunk)%field%vertexdy(k)=0.0
    ENDDO
!$OMP END DO

!$OMP END PARALLEL
  
END SUBROUTINE build_field
