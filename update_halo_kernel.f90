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

!>  @brief Fortran kernel to update the external halo cells in a chunk.
!>  @author Wayne Gaudin
!>  @details Updates halo cells for the required fields at the required depth
!>  for any halo cells that lie on an external boundary. The location and type
!>  of data governs how this is carried out. External boundaries are always
!>  reflective.

MODULE update_halo_kernel_module


  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,PRIVATE,PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,EXTERNAL_FACE=-1

  INTEGER,PRIVATE,PARAMETER :: FIELD_DENSITY0   = 1         &
                            ,FIELD_DENSITY1   = 2         &
                            ,FIELD_ENERGY0    = 3         &
                            ,FIELD_ENERGY1    = 4         &
                            ,FIELD_PRESSURE   = 5         &
                            ,FIELD_VISCOSITY  = 6         &
                            ,FIELD_SOUNDSPEED = 7         &
                            ,FIELD_XVEL0      = 8         &
                            ,FIELD_XVEL1      = 9         &
                            ,FIELD_YVEL0      =10         &
                            ,FIELD_YVEL1      =11         &
                            ,FIELD_VOL_FLUX_X =12         &
                            ,FIELD_VOL_FLUX_Y =13         &
                            ,FIELD_MASS_FLUX_X=14         &
                            ,FIELD_MASS_FLUX_Y=15       &
                            ,NUM_FIELDS       =15
CONTAINS

  SUBROUTINE update_halo_kernel(x_min,x_max,y_min,y_max,                            &
                        chunk_neighbours,                                           &
                        tile_neighbours,                                           &
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
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER, DIMENSION(4) :: chunk_neighbours, tile_neighbours

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,density1,energy0,energy1,pressure,soundspeed,viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x, mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y, mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0, xvel1, yvel0, yvel1

  INTEGER :: fields(NUM_FIELDS),depth

!$OMP PARALLEL

  ! Update values in external halo cells based on depth and fields requested
#define UPDATE_CELL_CHECK(array) IF (fields(FIELD_##array).EQ.1) CALL update_halo_cell(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, array, depth)
  
  UPDATE_CELL_CHECK(density0)
  UPDATE_CELL_CHECK(density1)
  UPDATE_CELL_CHECK(energy0)
  UPDATE_CELL_CHECK(energy1)
  UPDATE_CELL_CHECK(pressure)
  UPDATE_CELL_CHECK(soundspeed)
  UPDATE_CELL_CHECK(viscosity)

!$OMP END PARALLEL

END SUBROUTINE update_halo_kernel

SUBROUTINE update_halo_cell(x_min,x_max,y_min,y_max, &
                        chunk_neighbours,               &
                        tile_neighbours,               &
                        cell_mesh,                           &
                        depth                           )
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER, DIMENSION(4) :: chunk_neighbours, tile_neighbours
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: cell_mesh

  INTEGER :: depth

  INTEGER :: j,k

  IF (chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=y_min-depth,y_max+depth
      DO j=1,depth
        cell_mesh(1-j,k)=cell_mesh(0+j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF (chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=y_min-depth,y_max+depth
      DO j=1,depth
        cell_mesh(x_max+j,k)=cell_mesh(x_max+1-j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

!$OMP BARRIER

  IF (chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=1,depth
      DO j=x_min-depth,x_max+depth
        cell_mesh(j,1-k)=cell_mesh(j,0+k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF (chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=1,depth
      DO j=x_min-depth,x_max+depth
        cell_mesh(j,y_max+k)=cell_mesh(j,y_max+1-k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

END SUBROUTINE update_halo_cell

END MODULE update_halo_kernel_module
