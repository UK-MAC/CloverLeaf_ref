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

!>  @brief Main set up routine
!>  @author Wayne Gaudin
!>  @details Invokes the mesh decomposer and sets up chunk connectivity. It then
!>  allocates the communication buffers and call the chunk initialisation and
!>  generation routines. It calls the equation of state to calculate initial
!>  pressure before priming the halo cells and writing an initial field summary.

SUBROUTINE start

    USE clover_module
    USE parse_module
    USE update_halo_module
    USE ideal_gas_module

    IMPLICIT NONE

    INTEGER :: c, tile

    INTEGER :: x_cells,y_cells
    INTEGER:: right,left,top,bottom

    INTEGER :: fields(NUM_FIELDS) !, chunk_task_responsible_for 

    LOGICAL :: profiler_off

    IF(parallel%boss)THEN
       WRITE(g_out,*) 'Setting up initial geometry'
       WRITE(g_out,*)
    ENDIF

    time  = 0.0
    step  = 0
    dtold = dtinit
    dt    = dtinit

    CALL clover_barrier

    CALL clover_get_num_chunks(number_of_chunks)

    !ALLOCATE(chunk)

    !ALLOCATE(left(1:chunks_per_task))
    !ALLOCATE(right(1:chunks_per_task))
    !ALLOCATE(bottom(1:chunks_per_task))
    !ALLOCATE(top(1:chunks_per_task))

    CALL clover_decompose(grid%x_cells,grid%y_cells,left,right,bottom,top)

    !create the chunks 
      
    ! Needs changing so there can be more than 1 chunk per task
    chunk%task = parallel%task

    !chunk_task_responsible_for = parallel%task+1

    x_cells = right -left  +1
    y_cells = top   -bottom+1
      
    chunk%left    = left
    chunk%bottom  = bottom
    chunk%right   = right
    chunk%top     = top
    chunk%left_boundary   = 1
    chunk%bottom_boundary = 1
    chunk%right_boundary  = grid%x_cells
    chunk%top_boundary    = grid%y_cells
    chunk%x_min = 1
    chunk%y_min = 1
    chunk%x_max = x_cells
    chunk%y_max = y_cells
    
    
    

    !DEALLOCATE(left,right,bottom,top)


    ! create the tiles
    ALLOCATE( chunk%tiles(1:tiles_per_chunk) )

    CALL clover_tile_decompose(x_cells, y_cells)
    
    do tile=1, tiles_per_chunk
      !WRITE(*,*) parallel%task, tile, chunk%tiles(tile)%external_tile_mask(TILE_LEFT), chunk%tiles(tile)%external_tile_mask(TILE_RIGHT),chunk%tiles(tile)%external_tile_mask(TILE_BOTTOM), chunk%tiles(tile)%external_tile_mask(TILE_TOP)
      !WRITE(*,*) parallel%task, tile, chunk%tiles(tile)%t_xmin, chunk%tiles(tile)%t_xmax, chunk%tiles(tile)%t_ymin, chunk%tiles(tile)%t_ymax
      !WRITE(*,*) parallel%task, tile, chunk%tiles(tile)%t_left, chunk%tiles(tile)%t_right, chunk%tiles(tile)%t_bottom, chunk%tiles(tile)%t_top
    enddo

! set up the tiles within each chunk
!$OMP PARALLEL
!$OMP DO
    DO tile=1,tiles_per_chunk
        CALL build_field(tile)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL



  CALL clover_barrier

  CALL clover_allocate_buffers()

  IF(parallel%boss)THEN
     WRITE(g_out,*) 'Generating chunks'
  ENDIF

!$OMP PARALLEL
!$OMP DO
  DO tile=1,tiles_per_chunk
      CALL initialise_chunk(tile)
      CALL generate_chunk(tile)
  ENDDO
!$OMP END DO
!$OMP END PARALLEL



  advect_x=.TRUE.

  CALL clover_barrier

  ! Do no profile the start up costs otherwise the total times will not add up
  ! at the end
  profiler_off=profiler_on
  profiler_on=.FALSE.

!$OMP PARALLEL
!$OMP DO
  DO tile = 1, tiles_per_chunk
    CALL ideal_gas(tile,.FALSE.)
  END DO
!$OMP END DO
!$OMP END PARALLEL

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY0)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_PRESSURE)=1
  fields(FIELD_VISCOSITY)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1

  CALL update_halo(fields,2)

  IF(parallel%boss)THEN
     WRITE(g_out,*)
     WRITE(g_out,*) 'Problem initialised and generated'
  ENDIF

  CALL field_summary()

  IF(visit_frequency.NE.0) CALL visit()

  CALL clover_barrier

  profiler_on=profiler_off

END SUBROUTINE start
