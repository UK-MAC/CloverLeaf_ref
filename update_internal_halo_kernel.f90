
MODULE update_internal_halo_kernel_module

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

  SUBROUTINE update_internal_halo_left_right_kernel(                                &
                          x_min,x_max,y_min,y_max,                                    &
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
                          x_min_right,x_max_right,y_min_right,y_max_right,            &
                        density0_right,                                                   &
                        energy0_right,                                                    &
                        pressure_right,                                                   &
                        viscosity_right,                                                  &
                        soundspeed_right,                                                 &
                        density1_right,                                                   &
                        energy1_right,                                                    &
                        xvel0_right,                                                      &
                        yvel0_right,                                                      &
                        xvel1_right,                                                      &
                        yvel1_right,                                                      &
                        vol_flux_x_right,                                                 &
                        vol_flux_y_right,                                                 &
                        mass_flux_x_right,                                                &
                        mass_flux_y_right,                                                &
                          fields,                                                     &
                          depth                                                       )
    IMPLICIT NONE

    INTEGER :: fields(NUM_FIELDS),depth

    INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,density1,energy0,energy1,pressure,soundspeed,viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x, mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y, mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0, xvel1, yvel0, yvel1

    INTEGER :: x_min_right,x_max_right,y_min_right,y_max_right
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0_right,density1_right,energy0_right,energy1_right,pressure_right,soundspeed_right,viscosity_right
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x_right, mass_flux_x_right
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y_right, mass_flux_y_right
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0_right, xvel1_right, yvel0_right, yvel1_right

!$OMP PARALLEL
    ! FIXME
    IF (fields(FIELD_DENSITY0).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, density0, &
        x_min_right, x_max_right, y_min_right, y_max_right, density0_right, &
        depth)
    ENDIF
!$OMP END PARALLEL

  END SUBROUTINE

  SUBROUTINE update_internal_halo_cell_left_right(x_min_left,x_max_left,y_min_left,y_max_left,  &
                          mesh_left,   &
                          x_min_right,x_max_right,y_min_right,y_max_right,            &
                          mesh_right,                           &
                          depth                           )
    IMPLICIT NONE

    INTEGER :: depth

    INTEGER :: x_min_left,x_max_left,y_min_left,y_max_left
    REAL(KIND=8), DIMENSION(x_min_left-2:x_max_left+2,y_min_left-2:y_max_left+2) :: mesh_left

    INTEGER :: x_min_right,x_max_right,y_min_right,y_max_right
    REAL(KIND=8), DIMENSION(x_min_right-2:x_max_right+2,y_min_right-2:y_max_right+2) :: mesh_right

    INTEGER :: j,k

!$OMP DO
    DO k=y_min_left-depth,y_max_left+depth
!DIR$ IVDEP
      DO j=1,depth
        mesh_right(1-j,k)=mesh_left(x_max_left+1-j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP DO
    DO k=y_min_left-depth,y_max_left+depth
!DIR$ IVDEP
      DO j=1,depth
        mesh_left(x_max_left+j,k)=mesh_right(0+j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  END SUBROUTINE update_internal_halo_cell_left_right

  SUBROUTINE update_internal_halo_bottom_top_kernel(                                &
                          x_min,x_max,y_min,y_max,                                    &
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
                          x_min_top,x_max_top,y_min_top,y_max_top,            &
                        density0_top,                                                   &
                        energy0_top,                                                    &
                        pressure_top,                                                   &
                        viscosity_top,                                                  &
                        soundspeed_top,                                                 &
                        density1_top,                                                   &
                        energy1_top,                                                    &
                        xvel0_top,                                                      &
                        yvel0_top,                                                      &
                        xvel1_top,                                                      &
                        yvel1_top,                                                      &
                        vol_flux_x_top,                                                 &
                        vol_flux_y_top,                                                 &
                        mass_flux_x_top,                                                &
                        mass_flux_y_top,                                                &
                          fields,                                                     &
                          depth                                                       )
    IMPLICIT NONE

    INTEGER :: fields(NUM_FIELDS),depth

    INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,density1,energy0,energy1,pressure,soundspeed,viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x, mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y, mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0, xvel1, yvel0, yvel1

    INTEGER :: x_min_top,x_max_top,y_min_top,y_max_top
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0_top,density1_top,energy0_top,energy1_top,pressure_top,soundspeed_top,viscosity_top
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x_top, mass_flux_x_top
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y_top, mass_flux_y_top
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0_top, xvel1_top, yvel0_top, yvel1_top

!$OMP PARALLEL
    IF (fields(FIELD_DENSITY0).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, density0, &
        x_min_top, x_max_top, y_min_top, y_max_top, density0_top, &
        depth)
    ENDIF
!$OMP END PARALLEL

  END SUBROUTINE update_internal_halo_bottom_top_kernel

  SUBROUTINE update_internal_halo_cell_bottom_top(x_min_bottom,x_max_bottom,y_min_bottom,y_max_bottom,  &
                          mesh_bottom,   &
                          x_min_top,x_max_top,y_min_top,y_max_top,            &
                          mesh_top,                           &
                          depth                           )
    IMPLICIT NONE

    INTEGER :: depth

    INTEGER :: x_min_bottom,x_max_bottom,y_min_bottom,y_max_bottom
    REAL(KIND=8), DIMENSION(x_min_bottom-2:x_max_bottom+2,y_min_bottom-2:y_max_bottom+2) :: mesh_bottom

    INTEGER :: x_min_top,x_max_top,y_min_top,y_max_top
    REAL(KIND=8), DIMENSION(x_min_top-2:x_max_top+2,y_min_top-2:y_max_top+2) :: mesh_top

    INTEGER :: j,k

!$OMP DO
    DO k=1,depth
!DIR$ IVDEP
      DO j=x_min_bottom-depth,x_max_bottom+depth
        mesh_top(j,1-k)=mesh_bottom(j,y_max_bottom+1-k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP DO
    DO k=1,depth
!DIR$ IVDEP
      DO j=x_min_bottom-depth,x_max_bottom+depth
        mesh_bottom(j,y_max_bottom+k)=mesh_top(j,0+k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  END SUBROUTINE update_internal_halo_cell_bottom_top

END MODULE update_internal_halo_kernel_module


