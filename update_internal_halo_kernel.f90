
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

  SUBROUTINE update_internal_halo_kernel(                                &
                          x_min_1,x_max_1,y_min_1,y_max_1,               &
                        density0_1,                                                   &
                        energy0_1,                                                    &
                        pressure_1,                                                   &
                        viscosity_1,                                                  &
                        soundspeed_1,                                                 &
                        density1_1,                                                   &
                        energy1_1,                                                    &
                        xvel0_1,                                                      &
                        yvel0_1,                                                      &
                        xvel1_1,                                                      &
                        yvel1_1,                                                      &
                        vol_flux_x_1,                                                 &
                        vol_flux_y_1,                                                 &
                        mass_flux_x_1,                                                &
                        mass_flux_y_1,                                                &
                          x_min_2,x_max_2,y_min_2,y_max_2,            &
                        density0_2,                                                   &
                        energy0_2,                                                    &
                        pressure_2,                                                   &
                        viscosity_2,                                                  &
                        soundspeed_2,                                                 &
                        density1_2,                                                   &
                        energy1_2,                                                    &
                        xvel0_2,                                                      &
                        yvel0_2,                                                      &
                        xvel1_2,                                                      &
                        yvel1_2,                                                      &
                        vol_flux_x_2,                                                 &
                        vol_flux_y_2,                                                 &
                        mass_flux_x_2,                                                &
                        mass_flux_y_2,                                                &
                          fields,                                                     &
                          depth,                                                     &
                          face                                                           )

    IMPLICIT NONE

    INTEGER :: fields(NUM_FIELDS), depth, face

    INTEGER :: x_min_1,x_max_1,y_min_1,y_max_1
    REAL(KIND=8), DIMENSION(x_min_1-2:x_max_1+2,y_min_1-2:y_max_1+2) :: density0_1,density1_1,energy0_1,energy1_1,pressure_1,soundspeed_1,viscosity_1
    REAL(KIND=8), DIMENSION(x_min_1-2:x_max_1+3,y_min_1-2:y_max_1+2) :: vol_flux_x_1, mass_flux_x_1
    REAL(KIND=8), DIMENSION(x_min_1-2:x_max_1+2,y_min_1-2:y_max_1+3) :: vol_flux_y_1, mass_flux_y_1
    REAL(KIND=8), DIMENSION(x_min_1-2:x_max_1+3,y_min_1-2:y_max_1+3) :: xvel0_1, xvel1_1, yvel0_1, yvel1_1

    INTEGER :: x_min_2,x_max_2,y_min_2,y_max_2
    REAL(KIND=8), DIMENSION(x_min_2-2:x_max_2+2,y_min_2-2:y_max_2+2) :: density0_2,density1_2,energy0_2,energy1_2,pressure_2,soundspeed_2,viscosity_2
    REAL(KIND=8), DIMENSION(x_min_2-2:x_max_2+3,y_min_2-2:y_max_2+2) :: vol_flux_x_2, mass_flux_x_2
    REAL(KIND=8), DIMENSION(x_min_2-2:x_max_2+2,y_min_2-2:y_max_2+3) :: vol_flux_y_2, mass_flux_y_2
    REAL(KIND=8), DIMENSION(x_min_2-2:x_max_2+3,y_min_2-2:y_max_2+3) :: xvel0_2, xvel1_2, yvel0_2, yvel1_2

!$OMP PARALLEL

    IF (fields(FIELD_DENSITY0   ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, density0_1   , x_min_2, x_max_2, y_min_2, y_max_2, density0_2   , 0, 0, depth, face)
    IF (fields(FIELD_DENSITY1   ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, density1_1   , x_min_2, x_max_2, y_min_2, y_max_2, density1_2   , 0, 0, depth, face)
    IF (fields(FIELD_ENERGY0    ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, energy0_1    , x_min_2, x_max_2, y_min_2, y_max_2, energy0_2    , 0, 0, depth, face)
    IF (fields(FIELD_ENERGY1    ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, energy1_1    , x_min_2, x_max_2, y_min_2, y_max_2, energy1_2    , 0, 0, depth, face)
    IF (fields(FIELD_PRESSURE   ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, pressure_1   , x_min_2, x_max_2, y_min_2, y_max_2, pressure_2   , 0, 0, depth, face)
    IF (fields(FIELD_VISCOSITY  ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, viscosity_1  , x_min_2, x_max_2, y_min_2, y_max_2, viscosity_2  , 0, 0, depth, face)
    IF (fields(FIELD_SOUNDSPEED ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, soundspeed_1 , x_min_2, x_max_2, y_min_2, y_max_2, soundspeed_2 , 0, 0, depth, face)

    IF (fields(FIELD_XVEL0      ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, xvel0_1      , x_min_2, x_max_2, y_min_2, y_max_2, xvel0_2      , 1, 1, depth, face)
    IF (fields(FIELD_XVEL1      ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, xvel1_1      , x_min_2, x_max_2, y_min_2, y_max_2, xvel1_2      , 1, 1, depth, face)

    IF (fields(FIELD_YVEL0      ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, yvel0_1      , x_min_2, x_max_2, y_min_2, y_max_2, yvel0_2      , 1, 1, depth, face)
    IF (fields(FIELD_YVEL1      ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, yvel1_1      , x_min_2, x_max_2, y_min_2, y_max_2, yvel1_2      , 1, 1, depth, face)

    IF (fields(FIELD_VOL_FLUX_X ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, vol_flux_x_1 , x_min_2, x_max_2, y_min_2, y_max_2, vol_flux_x_2 , 1, 0, depth, face)
    IF (fields(FIELD_MASS_FLUX_X).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, mass_flux_x_1, x_min_2, x_max_2, y_min_2, y_max_2, mass_flux_x_2, 1, 0, depth, face)

    IF (fields(FIELD_VOL_FLUX_Y ).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, vol_flux_y_1 , x_min_2, x_max_2, y_min_2, y_max_2, vol_flux_y_2 , 0, 1, depth, face)
    IF (fields(FIELD_MASS_FLUX_Y).EQ.1) CALL update_internal_halo_inner(x_min_1, x_max_1, y_min_1, y_max_1, mass_flux_y_1, x_min_2, x_max_2, y_min_2, y_max_2, mass_flux_y_2, 0, 1, depth, face)

!$OMP END PARALLEL

  END SUBROUTINE

  SUBROUTINE update_internal_halo_inner(x_min_1,x_max_1,       &
                                        y_min_1,y_max_1,       &
                                        mesh_1,                                &
                                        x_min_2,x_max_2,    &
                                        y_min_2,y_max_2,    &
                                        mesh_2,                               &
                                        x_extra, y_extra,   &
                                        depth, face                     )

    IMPLICIT NONE

    INTEGER :: depth, face, x_extra, y_extra

    INTEGER :: x_min_1,x_max_1,y_min_1,y_max_1
    REAL(KIND=8), DIMENSION(x_min_1-2:x_max_1+2+x_extra,y_min_1-2:y_max_1+2+y_extra) :: mesh_1

    INTEGER :: x_min_2,x_max_2,y_min_2,y_max_2
    REAL(KIND=8), DIMENSION(x_min_2-2:x_max_2+2+x_extra,y_min_2-2:y_max_2+2+y_extra) :: mesh_2

    INTEGER :: j,k

    IF ((face .EQ. CHUNK_BOTTOM) .OR. (face .EQ. CHUNK_TOP)) THEN

!$OMP DO
      DO k=1,depth
        DO j=x_min_1-depth,x_max_1+depth+x_extra
          mesh_2(j,1-k)=mesh_1(j,y_max_1+1-k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

!$OMP DO
      DO k=1,depth
        DO j=x_min_1-depth,x_max_1+depth+x_extra
          mesh_1(j,y_max_1+y_extra+k)=mesh_2(j,y_extra+k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

    ELSEIF ((face .EQ. CHUNK_LEFT) .OR. (face .EQ. CHUNK_RIGHT)) THEN

!$OMP DO
      DO k=y_min_1-depth,y_max_1+depth+y_extra
        DO j=1,depth
          mesh_2(1-j,k)=mesh_1(x_max_1+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

!$OMP DO
      DO k=y_min_1-depth,y_max_1+depth+y_extra
        DO j=1,depth
          mesh_1(x_max_1+x_extra+j,k)=mesh_2(x_extra+j,k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

    ENDIF

  END SUBROUTINE update_internal_halo_inner

END MODULE update_internal_halo_kernel_module


