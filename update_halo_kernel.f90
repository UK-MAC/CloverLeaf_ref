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

  IF (fields(FIELD_DENSITY0   ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, density0,    depth, 0, 0,  1.0_8,  1.0_8)
  IF (fields(FIELD_DENSITY1   ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, energy0,     depth, 0, 0,  1.0_8,  1.0_8)
  IF (fields(FIELD_ENERGY0    ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, pressure,    depth, 0, 0,  1.0_8,  1.0_8)
  IF (fields(FIELD_ENERGY1    ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, viscosity,   depth, 0, 0,  1.0_8,  1.0_8)
  IF (fields(FIELD_PRESSURE   ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, soundspeed,  depth, 0, 0,  1.0_8,  1.0_8)
  IF (fields(FIELD_VISCOSITY  ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, density1,    depth, 0, 0,  1.0_8,  1.0_8)
  IF (fields(FIELD_SOUNDSPEED ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, energy1,     depth, 0, 0,  1.0_8,  1.0_8)

  IF (fields(FIELD_XVEL0      ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, xvel0,       depth, 1, 1, -1.0_8,  1.0_8)
  IF (fields(FIELD_XVEL1      ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, yvel0,       depth, 1, 1, -1.0_8,  1.0_8)

  IF (fields(FIELD_YVEL0      ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, xvel1,       depth, 1, 1,  1.0_8, -1.0_8)
  IF (fields(FIELD_YVEL1      ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, yvel1,       depth, 1, 1,  1.0_8, -1.0_8)

  IF (fields(FIELD_VOL_FLUX_X ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, vol_flux_x,  depth, 1, 0, -1.0_8,  1.0_8)
  IF (fields(FIELD_MASS_FLUX_X).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, mass_flux_x, depth, 1, 0, -1.0_8,  1.0_8)

  IF (fields(FIELD_VOL_FLUX_Y ).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, vol_flux_y,  depth, 0, 1,  1.0_8, -1.0_8)
  IF (fields(FIELD_MASS_FLUX_Y).EQ.1) CALL update_halo_inner(x_min, x_max, y_min, y_max, chunk_neighbours, tile_neighbours, mass_flux_y, depth, 0, 1,  1.0_8, -1.0_8)

!$OMP END PARALLEL

END SUBROUTINE update_halo_kernel

SUBROUTINE update_halo_inner(x_min,x_max,y_min,y_max, &
                        chunk_neighbours,               &
                        tile_neighbours,               &
                        mesh,                           &
                        depth,                      &
                        x_extra,y_extra,            &
                        x_invert,y_invert)
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER :: x_extra, y_extra
  INTEGER, DIMENSION(4) :: chunk_neighbours, tile_neighbours
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: mesh
  REAL(KIND=8) :: x_invert, y_invert

  INTEGER :: depth

  INTEGER :: j,k

  IF (chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=1,depth
      DO j=x_min-depth,x_max+depth
        mesh(j,1-k)=y_invert*mesh(j,y_extra+k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF (chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=1,depth
      DO j=x_min-depth,x_max+depth
        mesh(j,y_max+y_extra+k)=mesh(j,y_max+1-k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

!$OMP BARRIER

  IF (chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=y_min-depth,y_max+depth
      DO j=1,depth
        mesh(1-j,k)=x_invert*mesh(x_extra+j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF (chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=y_min-depth,y_max+depth
      DO j=1,depth
        mesh(x_max+x_extra+j,k)=mesh(x_max+1-j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

END SUBROUTINE update_halo_inner

END MODULE update_halo_kernel_module

