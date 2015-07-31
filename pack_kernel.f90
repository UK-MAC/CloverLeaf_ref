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

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

   ! These two need to be kept consistent with update_halo
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
                                ,FIELD_MASS_FLUX_Y=15         &
                                ,NUM_FIELDS       =15

   INTEGER,PRIVATE,PARAMETER :: CELL_DATA     = 1,        &
                                 VERTEX_DATA   = 2,        &
                                 X_FACE_DATA   = 3,        &
                                 y_FACE_DATA   = 4

CONTAINS

FUNCTION yincs(field_type) RESULT(y_inc)

  integer :: field_type, y_inc

  y_inc = 0

  IF (field_type.EQ.CELL_DATA) THEN
    y_inc=0
  ELSEIF (field_type.EQ.VERTEX_DATA) THEN
    y_inc=1
  ELSEIF (field_type.EQ.X_FACE_DATA) THEN
    y_inc=0
  ELSEIF (field_type.EQ.Y_FACE_DATA) THEN
    y_inc=1
  ENDIF

END FUNCTION

FUNCTION xincs(field_type) RESULT(x_inc)

  integer :: field_type, x_inc

  x_inc = 0

  IF (field_type.EQ.CELL_DATA) THEN
    x_inc=0
  ELSEIF (field_type.EQ.VERTEX_DATA) THEN
    x_inc=1
  ELSEIF (field_type.EQ.X_FACE_DATA) THEN
    x_inc=1
  ELSEIF (field_type.EQ.Y_FACE_DATA) THEN
    x_inc=0
  ENDIF

END FUNCTION

SUBROUTINE pack_all(x_min, x_max, y_min, y_max, &
    tile_neighbours, &
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
    fields, depth, face, packing, mpi_buffer, offsets, tile_offset)

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE pack_or_unpack(x_min,x_max,y_min,y_max,    &
                              field, mpi_buffer,          &
                              depth, x_inc, y_inc,        &
                              buffer_offset, edge_minus, edge_plus)

      IMPLICIT NONE

      INTEGER      :: depth,x_min,x_max,y_min,y_max,buffer_offset, x_inc, y_inc, edge_minus, edge_plus
      REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
      REAL(KIND=8) :: mpi_buffer(:)
    END SUBROUTINE
  END INTERFACE

  INTEGER      :: fields(:)
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER      :: face,tile_offset
  LOGICAL      :: packing
  INTEGER      :: depth,x_min,x_max,y_min,y_max, edge_minus, edge_plus
  INTEGER, DIMENSION(4) :: tile_neighbours

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,density1,energy0,energy1,pressure,soundspeed,viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x, mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y, mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0, xvel1, yvel0, yvel1

  PROCEDURE(pack_or_unpack), POINTER :: pack_func => NULL()

  SELECT CASE (face)
  CASE (CHUNK_LEFT, CHUNK_RIGHT)
    IF (tile_neighbours(CHUNK_BOTTOM) .EQ. EXTERNAL_FACE) THEN
      edge_minus = depth
    ELSE
      edge_minus = 0
    ENDIF

    IF (tile_neighbours(CHUNK_TOP) .EQ. EXTERNAL_FACE) THEN
      edge_plus = depth
    ELSE
      edge_plus = 0
    ENDIF
  CASE (CHUNK_BOTTOM, CHUNK_TOP)
    IF (tile_neighbours(CHUNK_LEFT) .EQ. EXTERNAL_FACE) THEN
      edge_minus = depth
    ELSE
      edge_minus = 0
    ENDIF

    IF (tile_neighbours(CHUNK_RIGHT) .EQ. EXTERNAL_FACE) THEN
      edge_plus = depth
    ELSE
      edge_plus = 0
    ENDIF
  END SELECT

  IF (packing .EQV. .TRUE.) THEN
    SELECT CASE (face)
    CASE (CHUNK_LEFT)
      pack_func => clover_pack_message_left
    CASE (CHUNK_RIGHT)
      pack_func => clover_pack_message_right
    CASE (CHUNK_BOTTOM)
      pack_func => clover_pack_message_bottom
    CASE (CHUNK_TOP)
      pack_func => clover_pack_message_top
    END SELECT
  ELSE
    SELECT CASE (face)
    CASE (CHUNK_LEFT)
      pack_func => clover_unpack_message_left
    CASE (CHUNK_RIGHT)
      pack_func => clover_unpack_message_right
    CASE (CHUNK_BOTTOM)
      pack_func => clover_unpack_message_bottom
    CASE (CHUNK_TOP)
      pack_func => clover_unpack_message_top
    END SELECT
  ENDIF

!$OMP PARALLEL
  IF (fields(FIELD_DENSITY0).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     density0,                 &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_DENSITY0),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_ENERGY0).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     energy0,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_ENERGY0),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_ENERGY1).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     energy1,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_ENERGY1),   &
                     edge_minus, edge_plus)
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE clover_pack_message_left(x_min,x_max,y_min,y_max,field,                &
                                 left_snd_buffer,                              &
                                 depth,x_inc, y_inc,                             &
                                 buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: left_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index=buffer_offset + j+(k+depth-1)*depth
      left_snd_buffer(index)=field(x_min+x_inc-1+j,k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_pack_message_left

SUBROUTINE clover_unpack_message_left(x_min,x_max,y_min,y_max,field,                &
                                   left_rcv_buffer,                              &
                                   depth,x_inc, y_inc,                             &
                                   buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: left_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index= buffer_offset + j+(k+depth-1)*depth
      field(x_min-j,k)=left_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_unpack_message_left

SUBROUTINE clover_pack_message_right(x_min,x_max,y_min,y_max,field,                &
                                  right_snd_buffer,                             &
                                   depth,x_inc, y_inc,                             &
                                  buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: right_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index= buffer_offset + j+(k+depth-1)*depth
      right_snd_buffer(index)=field(x_max+1-j,k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_pack_message_right

SUBROUTINE clover_unpack_message_right(x_min,x_max,y_min,y_max,field,                &
                                    right_rcv_buffer,                             &
                                   depth,x_inc, y_inc,                             &
                                    buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: right_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index= buffer_offset + j+(k+depth-1)*depth
      field(x_max+x_inc+j,k)=right_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_unpack_message_right

SUBROUTINE clover_pack_message_top(x_min,x_max,y_min,y_max,field,                &
                                top_snd_buffer,                               &
                                   depth,x_inc, y_inc,                             &
                                buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: top_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j+depth+(k-1)*(x_max+x_inc+(2*depth))
      top_snd_buffer(index)=field(j,y_max+1-k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_pack_message_top

SUBROUTINE clover_unpack_message_top(x_min,x_max,y_min,y_max,field,                &
                                  top_rcv_buffer,                               &
                                   depth,x_inc, y_inc,                             &
                                  buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: top_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j + depth+(k-1)*(x_max+x_inc+(2*depth))
      field(j,y_max+y_inc+k)=top_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_unpack_message_top

SUBROUTINE clover_pack_message_bottom(x_min,x_max,y_min,y_max,field,                &
                                   bottom_snd_buffer,                            &
                                   depth,x_inc, y_inc,                             &
                                   buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: bottom_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j+depth+(k-1)*(x_max+x_inc+(2*depth))
      bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_pack_message_bottom

SUBROUTINE clover_unpack_message_bottom(x_min,x_max,y_min,y_max,field,                &
                                     bottom_rcv_buffer,                            &
                                   depth,x_inc, y_inc,                             &
                                     buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-2:x_max+2+x_inc,y_min-2:y_max+2+y_inc)
  REAL(KIND=8) :: bottom_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j+depth+(k-1)*(x_max+x_inc+(2*depth))
      field(j,y_min-k)=bottom_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE clover_unpack_message_bottom

END MODULE pack_kernel_module
