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
                            ,FIELD_MASS_FLUX_Y=15
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
  INTEGER, DIMENSION(4) :: chunk_neighbours
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure,viscosity,soundspeed
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x,mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y,mass_flux_y
  INTEGER :: fields(:),depth

  INTEGER :: fields(NUM_FIELDS),depth

!$OMP PARALLEL

  ! Update values in external halo cells based on depth and fields requested
  IF (fields(FIELD_DENSITY).EQ.1) THEN
    CALL update_halo_cell(x_min, x_max, y_min, y_max, halo_exchange_depth,  &
      chunk_neighbours, tile_neighbours, density, depth)
  ENDIF

  IF (fields(FIELD_ENERGY0).EQ.1) THEN
    CALL update_halo_cell(x_min, x_max, y_min, y_max, halo_exchange_depth,  &
      chunk_neighbours, tile_neighbours, energy0, depth)
  ENDIF

  IF (fields(FIELD_ENERGY1).EQ.1) THEN
    CALL update_halo_cell(x_min, x_max, y_min, y_max, halo_exchange_depth,  &
      chunk_neighbours, tile_neighbours, energy1, depth)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          energy1(j,1-k)=energy1(j,0+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          energy1(j,y_max+k)=energy1(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          energy1(1-j,k)=energy1(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          energy1(x_max+j,k)=energy1(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          pressure(j,1-k)=pressure(j,0+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          pressure(j,y_max+k)=pressure(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          pressure(1-j,k)=pressure(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          pressure(x_max+j,k)=pressure(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          viscosity(j,1-k)=viscosity(j,0+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          viscosity(j,y_max+k)=viscosity(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          viscosity(1-j,k)=viscosity(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          viscosity(x_max+j,k)=viscosity(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          soundspeed(j,1-k)=soundspeed(j,0+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+depth
        DO k=1,depth
          soundspeed(j,y_max+k)=soundspeed(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          soundspeed(1-j,k)=soundspeed(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          soundspeed(x_max+j,k)=soundspeed(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          xvel0(j,1-k)=xvel0(j,1+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          xvel0(j,y_max+1+k)=xvel0(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel0(1-j,k)=-xvel0(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel0(x_max+1+j,k)=-xvel0(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          xvel1(j,1-k)=xvel1(j,1+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          xvel1(j,y_max+1+k)=xvel1(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel1(1-j,k)=-xvel1(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel1(x_max+1+j,k)=-xvel1(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          yvel0(j,1-k)=-yvel0(j,1+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          yvel0(j,y_max+1+k)=-yvel0(j,y_max+1-k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          yvel0(1-j,k)=yvel0(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          yvel0(x_max+1+j,k)=yvel0(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          yvel1(j,1-k)=-yvel1(j,1+k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO j=x_min-depth,x_max+1+depth
        DO k=1,depth
          yvel1(j,y_max+1+k)=-yvel1(j,y_max+1-k)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF (chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=y_min-depth,y_max+depth
      DO j=1,depth
        mesh(x_max+j,k)=mesh(x_max+1-j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

  ! Don't need barrier if depth is only 1
!$  IF (depth .gt. 1) then
!$OMP BARRIER
!$  ENDIF

  IF (chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=1,depth
      DO j=x_min-depth,x_max+depth
        mesh(j,1-k)=mesh(j,0+k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF (chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE .AND. tile_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO k=1,depth
      DO j=x_min-depth,x_max+depth
        mesh(j,y_max+k)=mesh(j,y_max+1-k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

END SUBROUTINE update_halo_cell

END MODULE update_halo_kernel_module
