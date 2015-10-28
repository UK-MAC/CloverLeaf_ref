
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

!>  @brief allocate and sets data for standalone mode
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode

MODULE set_data_module

CONTAINS

SUBROUTINE set_data(x_min,x_max,y_min,y_max,     &
                    cellx,                       &
                    celly,                       &
                    vertexdx,                    &
                    vertexdy,                    &
                    celldx,                      &
                    celldy,                      &
                    xarea,                       &
                    yarea,                       &
                    volume,                      &
                    density0,                    &
                    density1,                    &
                    energy0,                     &
                    energy1,                     &
                    viscosity,                   &
                    pressure,                    &
                    soundspeed,                  &
                    xvel0,                       &
                    xvel1,                       &
                    yvel0,                       &
                    yvel1,                       &
                    vol_flux_x,                  &
                    vol_flux_y,                  &
                    mass_flux_x,                 &
                    mass_flux_y,                 &
                    work_array1,                 &
                    work_array2,                 &
                    work_array3,                 &
                    work_array4,                 &
                    work_array5,                 &
                    work_array6,                 &
                    work_array7,                 &
                    dt                           )

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8),OPTIONAL :: dt
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: cellx(:),celly(:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vertexdx(:),vertexdy(:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: celldx(:),celldy(:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: xarea(:,:),yarea(:,:),volume(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: density0(:,:),density1(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: energy0(:,:),energy1(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: pressure(:,:),viscosity(:,:),soundspeed(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: xvel0(:,:),yvel0(:,:),xvel1(:,:),yvel1(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_x(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_y(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_x(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_y(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: work_array1(:,:),work_array2(:,:),work_array3(:,:),work_array4(:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: work_array5(:,:),work_array6(:,:),work_array7(:,:)

  INTEGER :: j,k
  REAL(KIND=8),ALLOCATABLE :: vertexx(:),vertexy(:)
  REAL(KIND=8) :: radius,theta,x,y,mult,dx,dy,sound_speed_squared,v,pressurebyenergy,pressurebyvolume,width
  REAL(KIND=8) :: ugrad,vgrad,div,strain2,pgradx,pgrady,pgradx2,pgrady2,limiter,pgrad,xgrad,ygrad,grad,grad2


  ! Set the initial data

  dx=(10.0_8)/float(x_max-x_min+1)
  dy=(10.0_8)/float(y_max-y_min+1)

  ALLOCATE(vertexx(x_min-2:x_max+3))
  ALLOCATE(vertexy(y_min-2:y_max+3))
  IF(PRESENT(cellx)) THEN
    ALLOCATE(cellx(x_min-2:x_max+2))
  ENDIF
  IF(PRESENT(celly)) THEN
    ALLOCATE(celly(y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(vertexdx)) THEN
    ALLOCATE(vertexdx(x_min-2:x_max+3))
  ENDIF
  IF(PRESENT(vertexdy)) THEN
    ALLOCATE(vertexdy(y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(celldx)) THEN
    ALLOCATE(celldx(x_min-2:x_max+2))
  ENDIF
  IF(PRESENT(celldy)) THEN
    ALLOCATE(celldy(y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(xarea)) THEN
    ALLOCATE(xarea(x_min-2:x_max+3 ,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(yarea)) THEN
    ALLOCATE(yarea(x_min-2:x_max+2 ,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(volume)) THEN
    ALLOCATE(volume(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(density0)) THEN
    ALLOCATE(density0(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(density1)) THEN
    ALLOCATE(density1(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(energy0)) THEN
    ALLOCATE(energy0(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(energy1)) THEN
    ALLOCATE(energy1(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(pressure)) THEN
    ALLOCATE(pressure(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(soundspeed)) THEN
    ALLOCATE(soundspeed(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(viscosity)) THEN
    ALLOCATE(viscosity(x_min-2:x_max+2,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(xvel0)) THEN
    ALLOCATE(xvel0(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(yvel0)) THEN
    ALLOCATE(yvel0(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(xvel1)) THEN
    ALLOCATE(xvel1(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(yvel1)) THEN
    ALLOCATE(yvel1(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array1)) THEN
    ALLOCATE(work_array1(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array2)) THEN
    ALLOCATE(work_array2(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array3)) THEN
    ALLOCATE(work_array3(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array4)) THEN
    ALLOCATE(work_array4(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array5)) THEN
    ALLOCATE(work_array5(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array6)) THEN
    ALLOCATE(work_array6(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(work_array7)) THEN
    ALLOCATE(work_array7(x_min-2:x_max+3,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(vol_flux_x)) THEN
    ALLOCATE(vol_flux_x(x_min-2:x_max+3,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(vol_flux_y)) THEN
    ALLOCATE(vol_flux_y(x_min-2:x_max+2,y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(mass_flux_x)) THEN
    ALLOCATE(mass_flux_x(x_min-2:x_max+3,y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(mass_flux_y)) THEN
    ALLOCATE(mass_flux_y(x_min-2:x_max+2,y_min-2:y_max+3))
  ENDIF

!$OMP PARALLEL

!$OMP DO
  DO j=x_min-2,x_max+3
     vertexx(j)=0.0_8+dx*float(j-x_min)
  ENDDO
!$OMP ENDDO

!$OMP DO
  DO k=y_min-2,y_max+3
     vertexy(k)=0.0_8+dy*float(k-y_min)
  ENDDO
!$OMP ENDDO

  IF(PRESENT(cellx)) THEN
!$OMP DO
    DO j=x_min-2,x_max+2
      cellx(j)=0.5*(vertexx(j)+vertexx(j+1))
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(celly)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      celly(k)=0.5*(vertexy(k)+vertexy(k+1))
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vertexdx)) THEN
!$OMP DO
    DO j=x_min-2,x_max+3
      vertexdx(j)=dx
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vertexdy)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      vertexdy(k)=dy
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(celldx)) THEN
!$OMP DO
    DO j=x_min-2,x_max+2
      celldx(j)=dx
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(celldy)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      celldy(k)=dy
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(xarea)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        xarea(j,k)=celldy(k)
       ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(yarea)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        yarea(j,k)=celldx(j)
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(volume)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        volume(j,k)=dx*dy
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(density0)) THEN
!$OMP DO PRIVATE(radius)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        radius=sqrt((float(j)*dx-5.0_8)**2.0_8+(float(k)*dy-5.0_8)**2.0_8)
        IF(radius.LE.2.5_8) THEN
          density0(j,k)=2.0_8-(radius*2.0_8/10.0_8)
        ELSE
          density0(j,k)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(density1)) THEN
!$OMP DO PRIVATE(radius)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        radius=sqrt((float(j)*dx-5.0_8)**2.0_8+(float(k)*dy-5.0_8)**2.0_8)
        IF(radius.LE.2.5_8) THEN
          density1(j,k)=2.0_8-(radius*2.0_8/10.0_8)
        ELSE
          density1(j,k)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF
 
  IF(PRESENT(energy0)) THEN
!$OMP DO PRIVATE(radius)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        radius=sqrt((float(j)*dx-5.0_8)**2.0_8+(float(k)*dy-5.0_8)**2.0_8)
        IF(radius.LE.2.5_8) THEN
          energy0(j,k)=2.0_8-(radius*2.0_8/10.0_8)
        ELSE
          energy0(j,k)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(energy1)) THEN
!$OMP DO PRIVATE(radius)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        radius=sqrt((float(j)*dx-5.0_8)**2.0_8+(float(k)*dy-5.0_8)**2.0_8)
        IF(radius.LE.2.5_8) THEN
          energy1(j,k)=2.0_8-(radius*2.0_8/10.0_8)
        ELSE
          energy1(j,k)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(pressure)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        pressure(j,k)=(1.4_8-1.0_8)*density0(j,k)*energy0(j,k)
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(soundspeed)) THEN
!$OMP DO PRIVATE(v,pressurebyenergy,pressurebyvolume,sound_speed_squared)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        v=1.0_8/density0(j,k)
        pressurebyenergy=(1.4_8-1.0_8)*density0(j,k)
        pressurebyvolume=-density0(j,k)*pressure(j,k)
        sound_speed_squared=v*v*(pressure(j,k)*pressurebyenergy-pressurebyvolume)
        soundspeed(j,k)=SQRT(sound_speed_squared)
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

! OMP THIS WITH A REDUCTION
  IF(PRESENT(dt)) THEN
    dt=0.0_8
    width=MIN(dx,dy)
!$OMP DO REDUCTION(MAX : dt)
    DO k=y_min,y_max
      DO j=x_min,x_max
        IF(soundspeed(j,k).GT.dt) dt=soundspeed(j,k)
      ENDDO
    ENDDO
!$OMP ENDDO

!$OMP MASTER
    dt=width*0.7_8/dt
!$OMP END MASTER 
!$OMP BARRIER
  ENDIF

  IF(PRESENT(xvel0)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,mult)
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        x=(float(j)*dx-5.0_8)
        y=(float(k)*dy-5.0_8)
        radius=sqrt(x**2.0_8+y**2.0_8)
        IF(x.LE.0.0_8.AND.y.LE.0.0_8) mult=-1.0_8
        IF(x.LE.0.0_8.AND.y.GT.0.0_8) mult=-1.0_8
        IF(x.GT.0.0_8.AND.y.LE.0.0_8) mult=1.0_8
        IF(x.GT.0.0_8.AND.y.GT.0.0_8) mult=1.0_8
        IF(x.NE.0.0_8) THEN
          theta=atan(y/x)
        ELSE
          theta=atan(y/-0.000000001_8)
        ENDIF
        IF(radius.LE.2.5_8) THEN
          xvel0(j,k)=mult*(2.0_8-(radius*2.0_8/10.0_8))*sin(theta)
        ELSE
          xvel0(j,k)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*sin(theta)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(yvel0)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,mult)
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        x=(float(j)*dx-5.0_8)
        y=(float(k)*dy-5.0_8)
        radius=sqrt(x**2.0_8+y**2.0_8)
        IF(x.LE.0.0_8.AND.y.LE.0.0_8) mult=-1.0_8
        IF(x.LE.0.0_8.AND.y.GT.0.0_8) mult=-1.0_8
        IF(x.GT.0.0_8.AND.y.LE.0.0_8) mult=1.0_8
        IF(x.GT.0.0_8.AND.y.GT.0.0_8) mult=1.0_8
        IF(x.NE.0.0_8) THEN
          theta=atan(y/x)
        ELSE
          theta=atan(y/-0.000000001_8)
        ENDIF
        radius=sqrt((float(j)*dx-5.0_8)**2.0_8+(float(k)*dy-5.0_8)**2.0_8)
        IF(radius.LE.2.5_8) THEN
          yvel0(j,k)=mult*(2.0_8-(radius*2.0_8/10.0_8))*cos(theta)
        ELSE
          yvel0(j,k)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*cos(theta)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(xvel1)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,mult)
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        x=(float(j)*dx-5.0_8)
        y=(float(k)*dy-5.0_8)
        radius=sqrt(x**2.0_8+y**2.0_8)
        IF(x.LE.0.0_8.AND.y.LE.0.0_8) mult=-1.0_8
        IF(x.LE.0.0_8.AND.y.GT.0.0_8) mult=-1.0_8
        IF(x.GT.0.0_8.AND.y.LE.0.0_8) mult=1.0_8
        IF(x.GT.0.0_8.AND.y.GT.0.0_8) mult=1.0_8
        IF(x.NE.0.0_8) THEN
          theta=atan(y/x)
        ELSE
          theta=atan(y/-0.000000001_8)
        ENDIF
        IF(radius.LE.2.5_8) THEN
          xvel1(j,k)=mult*(2.0_8-(radius*2.0_8/10.0_8))*sin(theta)
        ELSE
          xvel1(j,k)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*sin(theta)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(yvel1)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,mult)
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        x=(float(j)*dx-5.0_8)
        y=(float(k)*dy-5.0_8)
        radius=sqrt(x**2.0_8+y**2.0_8)
        IF(x.LE.0.0_8.AND.y.LE.0.0_8) mult=-1.0_8
        IF(x.LE.0.0_8.AND.y.GT.0.0_8) mult=-1.0_8
        IF(x.GT.0.0_8.AND.y.LE.0.0_8) mult=1.0_8
        IF(x.GT.0.0_8.AND.y.GT.0.0_8) mult=1.0_8
        IF(x.NE.0.0_8) THEN
          theta=atan(y/x)
        ELSE
          theta=atan(y/-0.000000001_8)
        ENDIF
        radius=sqrt((float(j)*dx-5.0_8)**2.0_8+(float(k)*dy-5.0_8)**2.0_8)
        IF(radius.LE.2.5_8) THEN
          yvel1(j,k)=mult*(2.0_8-(radius*2.0_8/10.0_8))*cos(theta)
        ELSE
          yvel1(j,k)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*cos(theta)
        ENDIF
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(viscosity)) THEN
!$OMP DO PRIVATE(ugrad,vgrad,div,strain2,pgradx,pgrady,pgradx2,pgrady2,limiter,pgrad,xgrad,ygrad,grad,grad2)
    DO k=y_min,y_max
      DO j=x_min,x_max
        ugrad=(xvel0(j+1,k  )+xvel0(j+1,k+1))-(xvel0(j  ,k  )+xvel0(j  ,k+1))

        vgrad=(yvel0(j  ,k+1)+yvel0(j+1,k+1))-(yvel0(j  ,k  )+yvel0(j+1,k  ))

        div = (celldx(j)*(ugrad)+  celldy(k)*(vgrad))

        strain2 = 0.5_8*(xvel0(j,  k+1) + xvel0(j+1,k+1)-xvel0(j  ,k  )-xvel0(j+1,k  ))/celldy(k) &
                + 0.5_8*(yvel0(j+1,k  ) + yvel0(j+1,k+1)-yvel0(j  ,k  )-yvel0(j  ,k+1))/celldx(j)

        pgradx=(pressure(j+1,k)-pressure(j-1,k))/(celldx(j)+celldx(j+1))
        pgrady=(pressure(j,k+1)-pressure(j,k-1))/(celldy(k)+celldy(k+1))

        pgradx2 = pgradx*pgradx
        pgrady2 = pgrady*pgrady

        limiter = ((0.5_8*(ugrad)/celldx(j))*pgradx2+(0.5_8*(vgrad)/celldy(k))*pgrady2+strain2*pgradx*pgrady)  &
                /MAX(pgradx2+pgrady2,1.0e-16_8)

        IF ((limiter.GT.0.0).OR.(div.GE.0.0))THEN
          viscosity(j,k) = 0.0
        ELSE
          pgradx = SIGN(MAX(1.0e-16_8,ABS(pgradx)),pgradx)
          pgrady = SIGN(MAX(1.0e-16_8,ABS(pgrady)),pgrady)
          pgrad = SQRT(pgradx*pgradx+pgrady*pgrady)
          xgrad = ABS(celldx(j)*pgrad/pgradx)
          ygrad = ABS(celldy(k)*pgrad/pgrady)
          grad  = MIN(xgrad,ygrad)
          grad2 = grad*grad

          viscosity(j,k)=2.0_8*density0(j,k)*grad2*limiter*limiter
        ENDIF

    ENDDO
  ENDDO
!$OMP END DO
  ENDIF

  IF(PRESENT(work_array1)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array1(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array2)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array2(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array3)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array3(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array4)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array4(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array5)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array5(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array6)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array6(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array7)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      DO j=x_min-2,x_max+3
        work_array7(j,k)=0.0_8
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vol_flux_x)) THEN
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max+1 
        vol_flux_x(j,k)=0.25_8*dt*xarea(j,k)                  &
                       *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1))
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vol_flux_y)) THEN
!$OMP DO
    DO k=y_min,y_max+1
      DO j=x_min,x_max
        vol_flux_y(j,k)=0.25_8*dt*yarea(j,k)                  &
                       *(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k))
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(mass_flux_x)) THEN
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max+1 
        ! j-1 could be j, depending on the flow
        mass_flux_x(j,k)=vol_flux_x(j,k)*density1(j-1,k)
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(mass_flux_y)) THEN
!$OMP DO
    DO k=y_min,y_max+1
      DO j=x_min,x_max
        ! k-1 could be k, depending on the flow
        mass_flux_y(j,k)=vol_flux_x(j,k)*density1(j,k-1)
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

!$OMP END PARALLEL

END SUBROUTINE set_data

END MODULE set_data_module

