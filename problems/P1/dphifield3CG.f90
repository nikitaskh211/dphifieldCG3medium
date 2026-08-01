
! Log of changes:
! (APR-23-2026) - Commented unused variables in "efield" program.
! (APR-23-2026) - Clarified comments in "efield" program.

! P1.1. A subroutine for assigning an analyticaly pre-determined solution vector and equation's right-hand-side (case A).
SUBROUTINE testset(bmid, rhs)

  USE numberformat
  USE indata
  USE constants

  IMPLICIT none

  COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(OUT) :: bmid, rhs ! Output arrays.

  REAL(rk)    :: kx, ky, kz ! ku = pi/Lu,
  REAL(rk)    :: k0, k2     ! k0 = pi/L0, k2 = kx.kx + ky.ky + kz.kz.

  COMPLEX(rk) :: jx, jy, jz    ! Jext = (Jx, Jy, Jz),
  COMPLEX(rk) :: jk            ! Jext.k = Jx.kx + Jy.ky + Jz.kz.

  COMPLEX(rk) :: dx0, dy0, dz0 ! D0 = (Dx0, Dy0, Dz0),
  COMPLEX(rk) :: fi0           ! Phi0 = fi0,
  COMPLEX(rk) :: hx0, hy0, hz0 ! H0 = (Hx0, Hy0, Hz0),
  COMPLEX(rk) :: ax0, ay0, az0 ! A0 = (Ax0, Ay0, Az0).

  INTEGER     :: i, j, k ! Iterator variables for x, y, and z,
 !INTEGER     :: l       ! N/A.

  REAL(rk)    :: argx, argy, argz ! Arguments of the trigonometric functions,
  REAL(rk)    :: sinx, siny, sinz ! Values of SIN(Arguments),
  REAL(rk)    :: cosx, cosy, cosz ! Values of SIN(Arguments).

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Wavevector.
  kx = pi/sizex
  ky = pi/sizey
  kz = pi/sizez

  k0 = pi/size0
  k2 = kx**2 + ky**2 + kz**2

  ! Electric current density.
  jx = 1.0D0
  jy = 0.0D0
  jz = 0.0D0

  jk = kx*jx + ky*jy + kz*jz

  ! Maxwell's fields constants.
  ax0 = (0.0D0, 1.0D0)*k0*(jx - kx*jk/k2)/(k2 - k0**2)
  ay0 = (0.0D0, 1.0D0)*k0*(jy - ky*jk/k2)/(k2 - k0**2)
  az0 = (0.0D0, 1.0D0)*k0*(jz - kz*jk/k2)/(k2 - k0**2)

  fi0 = - jk/k2

  hx0 = (ky*az0 - kz*ay0)/(0.0D0, 1.0D0)*k0
  hy0 = (kz*ax0 - kx*az0)/(0.0D0, 1.0D0)*k0
  hz0 = (kx*ay0 - ky*ax0)/(0.0D0, 1.0D0)*k0

  dx0 = ax0 - kx*fi0/(0.0D0, 1.0D0)*k0
  dy0 = ay0 - ky*fi0/(0.0D0, 1.0D0)*k0
  dz0 = az0 - kz*fi0/(0.0D0, 1.0D0)*k0

  ! Array initialization.
  rhs  = (0.0D0, 0.0D0)
  bmid = (0.0D0, 0.0D0)

  !===================!
  ! Array assignment. !
  !===================!

  ! Interior nodes.
  DO i = 2, nptx - 1

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
    sinx = sin(argx)

    DO j = 2, npty - 1

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)

      DO k = 2, nptz - 1

        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(1 , i, j, k) = (jx + kx*jk/k0**2)*(0.0D0, 1.0D0)*k0*cosx*siny*sinz
        rhs(2 , i, j, k) = (jy + ky*jk/k0**2)*(0.0D0, 1.0D0)*k0*sinx*cosy*sinz
        rhs(3 , i, j, k) = (jz + kz*jk/k0**2)*(0.0D0, 1.0D0)*k0*sinx*siny*cosz

        ! rhs(4) = 0
        rhs(4 , i, j, k) = 0.0D0*sinx*siny*sinz

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(5 , i, j, k) = - (ky*jz - kz*jy)*sinx*cosy*cosz
        rhs(6 , i, j, k) = - (kz*jx - kx*jz)*cosx*siny*cosz
        rhs(7 , i, j, k) = - (kx*jy - ky*jx)*cosx*cosy*sinz

        ! rhs(8-10) = 0
        rhs(8 , i, j, k) = 0.0D0*cosx*siny*sinz
        rhs(9 , i, j, k) = 0.0D0*sinx*cosy*sinz
        rhs(10, i, j, k) = 0.0D0*sinx*siny*cosz

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(1 , i, j, k) = dx0*cosx*siny*sinz
        bmid(2 , i, j, k) = dy0*sinx*cosy*sinz
        bmid(3 , i, j, k) = dz0*sinx*siny*cosz

        ! x(4) = Phi
        bmid(4 , i, j, k) = fi0*sinx*siny*sinz

        ! x(5-7) = H
        bmid(5 , i, j, k) = hx0*sinx*cosy*cosz
        bmid(6 , i, j, k) = hy0*cosx*siny*cosz
        bmid(7 , i, j, k) = hz0*cosx*cosy*sinz

        ! x(8-10) = A
        bmid(8 , i, j, k) = ax0*cosx*siny*sinz
        bmid(9 , i, j, k) = ay0*sinx*cosy*sinz
        bmid(10, i, j, k) = az0*sinx*siny*cosz
      
      END DO
    END DO
  END DO

  ! Boubdary nodes of x = 0 and x = Lx.
  DO i = 1, nptx, nptx - 1

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
   !sinx = sin(argx)
    sinx = 0.0D0

    DO j = 1, npty

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)
      IF(j == 1 .OR. j == npty) siny = 0.0D0

      DO k = 1, nptz

        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)
        IF(k == 1 .OR. k == nptz) sinz = 0.0D0

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(1, i, j, k) = (jx + kx*jk/k0**2)*(0.0D0, 1.0D0)*k0*cosx*siny*sinz/2.0D0

        ! rhs(4) = 0

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(6, i, j, k) = - (kz*jx - kx*jz)*cosx*siny*cosz/2.0D0
        rhs(7, i, j, k) = - (kx*jy - ky*jx)*cosx*cosy*sinz/2.0D0

        ! rhs(8-10) = 0

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(1,i,j,k)=dx0*cosx*siny*sinz

        ! x(4) = Phi

        ! x(5-7) = H
        bmid(6,i,j,k)=hy0*cosx*siny*cosz
        bmid(7,i,j,k)=hz0*cosx*cosy*sinz

        ! x(8-10) = A
        bmid(8,i,j,k)=ax0*cosx*siny*sinz
    
      END DO
    END DO
  END DO

  ! Boubdary nodes of y = 0 and y = Ly.
  DO i = 1, nptx

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
    sinx = sin(argx)
    IF(i == 1 .OR. i == nptx) sinx = 0.0D0

    DO j = 1, npty, npty-1

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)
      IF(j == 1 .OR. j == npty) siny = 0.0D0

      DO k = 1, nptz

        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)
        IF(k == 1 .OR. k == nptz) sinz = 0.0D0

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(2, i, j, k) = (jy + ky*jk/k0**2)*(0.0D0, 1.0D0)*k0*sinx*cosy*sinz/2.0D0

        ! rhs(4) = 0

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(5, i, j, k) = - (ky*jz - kz*jy)*sinx*cosy*cosz/2.0D0
        rhs(7, i, j, k) = - (kx*jy - ky*jx)*cosx*cosy*sinz/2.0D0

        ! rhs(8-10) = 0

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(2, i, j, k) = dy0*sinx*cosy*sinz

        ! x(4) = Phi
        bmid(5, i, j, k) = hx0*sinx*cosy*cosz
        bmid(7, i, j, k) = hz0*cosx*cosy*sinz

        ! x(5-7) = H

        ! x(8-10) = A
        bmid(9, i, j, k) = ay0*sinx*cosy*sinz

      END DO
    END DO
  END DO

  ! Boubdary nodes of z = 0 and z = Lz.
  DO i = 1, nptx

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
    sinx = sin(argx)
    IF(i == 1 .OR. i == nptx) sinx = 0.0D0

    DO j = 1, npty

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)
      IF(j == 1 .OR. j == npty) siny = 0.0D0

      DO k = 1, nptz, nptz - 1
        
        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)
        IF(k == 1 .OR. k == nptz) sinz = 0.0D0

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(3, i, j, k) = (jz + kz*jk/k0**2)*(0.0D0, 1.0D0)*k0*sinx*siny*cosz/2.0D0

        ! rhs(4) = 0

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(5, i, j, k) = - (ky*jz - kz*jy)*sinx*cosy*cosz/2.0D0
        rhs(6, i, j, k) = - (kz*jx - kx*jz)*cosx*siny*cosz/2.0D0

        ! rhs(8-10) = 0

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(3, i, j, k) = dz0*sinx*siny*cosz

        ! x(4) = Phi

        ! x(5-7) = H
        bmid(5, i, j, k) = hx0*sinx*cosy*cosz
        bmid(6, i, j, k) = hy0*cosx*siny*cosz

        ! x(8-10) = A
        bmid(10, i, j, k) = az0*sinx*siny*cosz
  
      END DO
    END DO
  END DO

  !===============================!
  ! Unused zeroing of the arrays. !
  !===============================!
 !rhs(1, :, :, :) = (0.0D0,0.0D0)
 !rhs(2, :, :, :) = (0.0D0,0.0D0)
 !rhs(3, :, :, :) = (0.0D0,0.0D0)
 !rhs(4, :, :, :) = (0.0D0,0.0D0)
 !rhs(5, :, :, :) = (0.0D0,0.0D0)
 !rhs(6, :, :, :) = (0.0D0,0.0D0)
 !rhs(7, :, :, :) = (0.0D0,0.0D0)

 !bmid(2  , :, :, :) = (0.0D0,0.0D0)
 !bmid(4:7, :, :, :) = (0.0D0,0.0D0)

  !=====================================!
  ! Adjustments for the volume element. !
  !=====================================!

  ! Edge: Hx at y = {0, Ly} and z = {0, Lz}
  rhs(5, :, 1, 1)       = rhs(5, :, 1, 1)*0.5D0
  rhs(5, :, npty, 1)    = rhs(5, :, npty, 1)*0.5D0
  rhs(5, :, 1, nptz)    = rhs(5, :, 1, nptz)*0.5D0
  rhs(5, :, npty, nptz) = rhs(5, :, npty, nptz)*0.5D0

  ! Edge: Hy at x = {0, Lx} and z = {0, Lz}
  rhs(6, 1, :, 1)       = rhs(6, 1, :, 1)*0.5D0
  rhs(6, nptx, :, 1)    = rhs(6, nptx, :, 1)*0.5D0
  rhs(6, 1, :, nptz)    = rhs(6, 1, :, nptz)*0.5D0
  rhs(6, nptx, :, nptz) = rhs(6, nptx, :, nptz)*0.5D0

  ! Edge: Hz at x = {0, Lx} and y = {0, Ly}
  rhs(7, 1, 1, :)       = rhs(7, 1, 1, :)*0.5D0
  rhs(7, nptx, 1, :)    = rhs(7, nptx, 1, :)*0.5D0
  rhs(7, 1, npty, :)    = rhs(7, 1, npty, :)*0.5D0
  rhs(7, nptx, npty, :) = rhs(7, nptx, npty, :)*0.5D0

END SUBROUTINE testset

! P2.2. A subroutine for assigning an analyticaly pre-determined solution vector and equation's right-hand-side (case B).
SUBROUTINE testset1(bmid, rhs)

  USE numberformat
  USE indata
  USE constants

  IMPLICIT none

  COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(OUT) :: bmid, rhs ! Output arrays.

  REAL(rk)    :: kx, ky, kz ! ku = pi/Lu,
  REAL(rk)    :: k0         ! k0 = pi/L0.

  COMPLEX(rk) :: dx0, dy0, dz0 ! D0 = (Dx0, Dy0, Dz0),
  COMPLEX(rk) :: fi0           ! Phi0 = fi0,
  COMPLEX(rk) :: hx0, hy0, hz0 ! H0 = (Hx0, Hy0, Hz0),
  COMPLEX(rk) :: ax0, ay0, az0 ! A0 = (Ax0, Ay0, Az0).

  INTEGER     :: i, j, k ! Iterator variables for x, y, and z,
 !INTEGER     :: l       ! N/A.

  REAL(rk)    :: argx, argy, argz ! Arguments of the trigonometric functions,
  REAL(rk)    :: sinx, siny, sinz ! Values of SIN(Arguments),
  REAL(rk)    :: cosx, cosy, cosz ! Values of SIN(Arguments).

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Wavevector.
  kx = pi/sizex
  ky = pi/sizey
  kz = pi/sizez

  k0 = pi/size0

  ! Maxwell's fields constants.
  dx0 = (2.0D0, 0.0D0)
  dy0 = (0.5D0, 0.0D0)
  dz0 = (1.0D0, 0.0D0)

  fi0 = (0.0D0, 1.0D0)

  hx0 = (0.0D0, 2.0D0)
  hy0 = (0.0D0, 1.0D0)
  hz0 =  - (kx*hx0 + ky*hy0)/kz

  ax0 = (1.0D0, 0.0D0)
  ay0 = (2.0D0, 0.0D0)
  az0 = - (kx*ax0 + ky*ay0)/kz

  ! Array initialization.
  rhs  = (0.0D0, 0.0D0)
  bmid = (0.0D0, 0.0D0)

  !===================!
  ! Array assignment. !
  !===================!

  ! Interior nodes.
  DO i = 2, nptx - 1

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
    sinx = sin(argx)

    DO j = 2, npty - 1

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)

      DO k = 2, nptz - 1

        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(1, i, j, k) = - (kx**2*dx0 + kx*ky*dy0 + kx*kz*dz0 + 2.0D0*k0**2*dx0 - k0**2*ax0 - (0.0D0, 1.0D0)*k0*kx*fi0 + (0.0D0, 0.0D0)*k0*(ky*hz0 - kz*hy0))*cosx*siny*sinz
        rhs(2, i, j, k) = - (kx*ky*dx0 + ky**2*dy0 + ky*kz*dz0 + 2.0D0*k0**2*dy0 - k0**2*ay0 - (0.0D0, 1.0D0)*k0*ky*fi0 + (0.0D0, 0.0D0)*k0*(kz*hx0 - kx*hz0))*sinx*cosy*sinz
        rhs(3, i, j, k) = - (kx*kz*dx0 + ky*kz*dy0 + kz**2*dz0 + 2.0D0*k0**2*dz0 - k0**2*az0 - (0.0D0, 1.0D0)*k0*kz*fi0 + (0.0D0, 0.0D0)*k0*(kx*hy0 - ky*hx0))*sinx*siny*cosz
   
        ! rhs(4) = 0
        rhs(4, i, j, k) = - ((0.0D0, 1.0D0)*k0*(kx*dx0 + ky*dy0 + kz*dz0) - (0.0D0, 1.0D0)*k0*(kx*ax0 + ky*ay0 + kz*az0) + (kx**2 + ky**2 + kz**2)*fi0)*sinx*siny*sinz

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(5, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hx0 + (0.0D0, 1.0D0)*k0*(ky*dz0 - kz*dy0) + (0.0D0, 1.0D0)*k0*(ky*az0 - kz*ay0))*sinx*cosy*cosz
        rhs(6, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hy0 + (0.0D0, 1.0D0)*k0*(kz*dx0 - kx*dz0) + (0.0D0, 1.0D0)*k0*(kz*ax0 - kx*az0))*cosx*siny*cosz
        rhs(7, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hz0 + (0.0D0, 1.0D0)*k0*(kx*dy0 - ky*dx0) + (0.0D0, 1.0D0)*k0*(kx*ay0 - ky*ax0))*cosx*cosy*sinz
   
        ! rhs(8-10) = 0
        rhs(8 , i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*ax0 + (0.0D0, 1.0D0)*k0*(ky*hz0 - kz*hy0) - k0**2*dx0 + (0.0D0, 1.0D0)*k0*kx*fi0)*cosx*siny*sinz
        rhs(9 , i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*ay0 + (0.0D0, 1.0D0)*k0*(kz*hx0 - kx*hz0) - k0**2*dy0 + (0.0D0, 1.0D0)*k0*ky*fi0)*sinx*cosy*sinz
        rhs(10, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*az0 + (0.0D0, 1.0D0)*k0*(kx*hy0 - ky*hx0) - k0**2*dz0 + (0.0D0, 1.0D0)*k0*kz*fi0)*sinx*siny*cosz

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(1 , i, j, k) = dx0*cosx*siny*sinz
        bmid(2 , i, j, k) = dy0*sinx*cosy*sinz
        bmid(3 , i, j, k) = dz0*sinx*siny*cosz
   
        ! x(4) = Phi
        bmid(4 , i, j, k) = fi0*sinx*siny*sinz

        ! x(5-7) = H
        bmid(5 , i, j, k) = hx0*sinx*cosy*cosz
        bmid(6 , i, j, k) = hy0*cosx*siny*cosz
        bmid(7 , i, j, k) = hz0*cosx*cosy*sinz

        ! x(8-10) = A
        bmid(8 , i, j, k) = ax0*cosx*siny*sinz
        bmid(9 , i, j, k) = ay0*sinx*cosy*sinz
        bmid(10, i, j, k) = az0*sinx*siny*cosz
      
      END DO
    END DO
  END DO

  ! Boubdary nodes of x = 0 and x = Lx.
  DO i = 1, nptx, nptx - 1

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
   !sinx = sin(argx)
    sinx = 0.0D0

    DO j = 1, npty

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)
      IF(j == 1 .OR. j == npty) siny = 0.0D0

      DO k = 1, nptz

        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)
        IF(k == 1 .OR. k == nptz) sinz = 0.0D0

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(1, i, j, k) = - (kx**2*dx0 + kx*ky*dy0 + kx*kz*dz0 + 2.0D0*k0**2*dx0 - k0**2*ax0 - (0.0D0, 1.0D0)*k0*kx*fi0 + (0.0D0, 1.0D0)*k0*(ky*hz0 - kz*hy0))*cosx*siny*sinz/2.0D0

        ! rhs(4) = 0

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(6, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hy0 + (0.0D0, 1.0D0)*k0*(kz*dx0 - kx*dz0) + (0.0D0, 1.0D0)*k0*(kz*ax0 - kx*az0))*cosx*siny*cosz/2.0D0
        rhs(7, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hz0 + (0.0D0, 1.0D0)*k0*(kx*dy0 - ky*dx0) + (0.0D0, 1.0D0)*k0*(kx*ay0 - ky*ax0))*cosx*cosy*sinz/2.0D0

        ! rhs(8-10) = 0
        rhs(8, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*ax0 + (0.0D0, 1.0D0)*k0*(ky*hz0 - kz*hy0) - k0**2*dx0 + (0.0D0, 1.0D0)*k0*kx*fi0)*cosx*siny*sinz/2.0D0

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(1, i, j, k) = dx0*cosx*siny*sinz

        ! x(4) = Phi

        ! x(5-7) = H
        bmid(6, i, j, k) = hy0*cosx*siny*cosz
        bmid(7, i, j, k) = hz0*cosx*cosy*sinz

        ! x(8-10) = A
        bmid(8, i, j, k) = ax0*cosx*siny*sinz
    
      END DO
    END DO
  END DO

  ! Boubdary nodes of y = 0 and y = Ly.
  DO i = 1, nptx

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
    sinx = sin(argx)
    IF(i == 1 .OR. i == nptx) sinx = 0.0D0

    DO j = 1, npty, npty - 1

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)
      IF(j == 1 .OR. j == npty) siny = 0.0D0

      DO k = 1, nptz

        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)
        IF(k == 1 .OR. k == nptz) sinz = 0.0D0

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(2, i, j, k) = - (kx*ky*dx0 + ky**2*dy0 + ky*kz*dz0 + 2.0D0*k0**2*dy0 - k0**2*ay0 - (0.0D0, 1.0D0)*k0*ky*fi0 + (0.0D0, 1.0D0)*k0*(kz*hx0 - kx*hz0))*sinx*cosy*sinz/2.0D0

        ! rhs(4) = 0

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(5, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hx0 + (0.0D0, 1.0D0)*k0*(ky*dz0 - kz*dy0) + (0.0D0, 1.0D0)*k0*(ky*az0 - kz*ay0))*sinx*cosy*cosz/2.0D0
        rhs(7, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hz0 + (0.0D0, 1.0D0)*k0*(kx*dy0 - ky*dx0) + (0.0D0, 1.0D0)*k0*(kx*ay0 - ky*ax0))*cosx*cosy*sinz/2.0D0

        ! rhs(8-10) = 0
        rhs(9, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*ay0 + (0.0D0, 1.0D0)*k0*(kz*hx0 - kx*hz0) - k0**2*dy0 + (0.0D0, 1.0D0)*k0*ky*fi0)*sinx*cosy*sinz/2.0D0

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(2, i, j, k) = dy0*sinx*cosy*sinz

        ! x(4) = Phi

        ! x(5-7) = H
        bmid(5, i, j, k) = hx0*sinx*cosy*cosz
        bmid(7, i, j, k) = hz0*cosx*cosy*sinz

        ! x(8-10) = A
        bmid(9, i, j, k) = ay0*sinx*cosy*sinz
       
      END DO
    END DO
  END DO

  ! Boubdary nodes of z = 0 and z = Lz.
  DO i = 1, nptx

    ! Calculation of cos(x) and sin(x).
    argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
    cosx = cos(argx)
    sinx = sin(argx)
    IF(i == 1 .OR. i == nptx) sinx = 0.0D0

    DO j = 1, npty

      ! Calculation of cos(y) and sin(y).
      argy = pi*(j - 1.0D0)/(npty - 1.0D0)
      cosy = cos(argy)
      siny = sin(argy)
      IF(j == 1 .OR. j == npty) siny = 0.0D0

      DO k = 1, nptz, nptz - 1
        
        ! Calculation of cos(z) and sin(z).
        argz = pi*(k - 1.0D0)/(nptz - 1.0D0)
        cosz = cos(argz)
        sinz = sin(argz)
        IF(k == 1 .OR. k == nptz) sinz = 0.0D0

        !=============================!
        ! Right-hand-side assignment. !
        !=============================!

        ! rhs(1-3) = - i.k0.4.pi.Jext/c - 4.pi.grad(rhoext)
        rhs(3 , i, j, k) = - (kx*kz*dx0 + ky*kz*dy0 + kz**2*dz0 + 2.0D0*k0**2*dz0 - k0**2*az0 - (0.0D0, 1.0D0)*k0*kz*fi0 + (0.0D0, 1.0D0)*k0*(kx*hy0 - ky*hx0))*sinx*siny*cosz/2.0D0

        ! rhs(4) = 0

        ! rhs(5-7) = 4.pi.curl(Jext)/c
        rhs(5 , i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hx0 + (0.0D0, 1.0D0)*k0*(ky*dz0 - kz*dy0) + (0.0D0, 1.0D0)*k0*(ky*az0 - kz*ay0))*sinx*cosy*cosz/2.0D0
        rhs(6 , i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*hy0 + (0.0D0, 1.0D0)*k0*(kz*dx0 - kx*dz0) + (0.0D0, 1.0D0)*k0*(kz*ax0 - kx*az0))*cosx*siny*cosz/2.0D0

        ! rhs(8-10) = 0
        rhs(10, i, j, k) = - ((kx**2 + ky**2 + kz**2 + k0**2)*az0 + (0.0D0, 1.0D0)*k0*(kx*hy0 - ky*hx0) - k0**2*dz0 + (0.0D0, 1.0D0)*k0*kz*fi0)*sinx*siny*cosz/2.0D0

        !====================!
        ! Fields assignment. !
        !====================!

        ! x(1-3) = D
        bmid(3 , i, j, k) = dz0*sinx*siny*cosz

        ! x(4) = Phi

        ! x(5-7) = H
        bmid(5 , i, j, k) = hx0*sinx*cosy*cosz
        bmid(6 , i, j, k) = hy0*cosx*siny*cosz

        ! x(8-10) = A
        bmid(10, i, j, k) = az0*sinx*siny*cosz
  
      END DO
    END DO
  END DO

  !===============================!
  ! Unused zeroing of the arrays. !
  !===============================!
 !rhs(1,:,:,:)=(0.0D0,0.0D0)
 !rhs(2,:,:,:)=(0.0D0,0.0D0)
 !rhs(3,:,:,:)=(0.0D0,0.0D0)
 !rhs(4,:,:,:)=(0.0D0,0.0D0)
 !rhs(5,:,:,:)=(0.0D0,0.0D0)
 !rhs(6,:,:,:)=(0.0D0,0.0D0)
 !rhs(7,:,:,:)=(0.0D0,0.0D0)

 !rhs(2,:,:,:)=(0.0D0,0.0D0)
 !rhs(4,:,:,:)=(0.0D0,0.0D0)
 !bmid(2,:,:,:)=(0.0D0,0.0D0)
 !bmid(4:7,:,:,:)=(0.0D0,0.0D0)

END SUBROUTINE testset1

! P3.3. A program that calculates a solution to the test example in Moiseenko and Agren article.
PROGRAM efield

  USE numberformat
  USE indata

  IMPLICIT NONE

  !================================================!
  ! Discretized analytical solution dynamic array. !
  !================================================!
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: bmid ! Is to contain the analytical solution from the numeric example.
  
  !====================!
  ! CG dynamic arrays. !
  !====================!
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: rhs   ! Is to contain the right-hand-side of the positive definite form of equations,
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: xx    ! Is to contain the approximately computed solution,
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: rr    ! Is to contain the residual,
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: pp    ! Is to contain the search direction,
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: ap    ! Is to contain U(v) operator action for any arbitrary 10D vector v,
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: tempr ! Temporary 3D array for computation.

  !===============!
  ! CG iterators. !
  !===============!
  INTEGER, PARAMETER :: numiter = 20000 ! Maximum amount of iterations for a given problem,
  INTEGER            :: ii              ! Iterator variable for CG.

  !===============!
  ! CG variables. !
  !===============!
  COMPLEX(rk) :: alpha, beta                   ! The step size and conjugacy coefficient,
  REAL(rk)    :: rhsnorm, solnorm, dprr, dprrn ! |rhs|^2, |bmid|^2, |r0|^2, |r1|^2,
  REAL(rk)    :: err, errsol                   ! The normalized error and residual.

  !=====================!
  ! External functions. !
  !=====================!
  COMPLEX(rk), EXTERNAL :: scprod ! This is a function that computes the discrete inner product,
  REAL(rk), EXTERNAL    :: solerr ! This is a function that computes the normalized error.

  !===================!
  ! Unused variables. !
  !===================!
 !REAL(rk), EXTERNAL::errcalc
 !INTEGER :: imx,jmx,kmx,lmx
 !INTEGER :: i,j,k,l
 !REAL(rk)::rmx,rmxcur

  !================================!
  ! Reading the input from a file. !
  !================================!

  CALL datain

  !====================================!
  ! Allocating memory for computation. !
  !====================================!

  ALLOCATE(bmid(10, nptx, npty, nptz))

  ALLOCATE(xx(10, nptx, npty, nptz))
  ALLOCATE(rr(10, nptx, npty, nptz))
  ALLOCATE(pp(10, nptx, npty, nptz))
  ALLOCATE(rhs(10, nptx, npty, nptz))
  ALLOCATE(ap(10, nptx, npty, nptz))
  ALLOCATE(tempr(3, nptx, npty, nptz))

  !======================================!
  ! Assigning a pre-defined expressions. !
  !======================================!

  CALL testset(bmid, rhs)

 !rhs=-rhs
 !CALL mulmat(bmid,xx)
 !stop

  !==================!
  ! Computing norms. !
  !==================!

  rhsnorm = REAL(scprod(rhs,rhs))
  solnorm = REAL(scprod(bmid,bmid))

  !===================================!
  ! Doing an initial iteration of CG. !
  !===================================!

  CALL mulmat(bmid,xx,tempr)

 !rhs=-rhs

  errsol = SQRT(solerr(xx, rhs))/SQRT(solnorm)
  print *, errsol

 !OPEN(22, file='test.dat')
 !!
 !rmx=0.0D0
 !
 !DO i=1,nptx
 !
 ! DO j=1,npty
 !
 !  DO k=1,nptz
 !   !DO l=1,nrotrot
 !   DO l=1,10
 !    WRITE(22,'(4i5,4e15.5)'),i,j,k,l,xx(l,i,j,k),rhs(l,i,j,k)
 !
 !    rmxcur=ABS(xx(l,i,j,k)-rhs(l,i,j,k))
 !
 !    IF(rmxcur>rmx) THEN
 !     rmx=rmxcur
 !     imx=i
 !     jmx=j
 !     kmx=k
 !     lmx=l
 !    END IF
 !
 !   END DO
 !   
 !  END DO
 ! END DO
 !END DO
 !
 !PRINT *,imx,jmx,kmx,lmx,xx(lmx,imx,jmx,kmx),rhs(lmx,imx,jmx,kmx)
 !
 !CLOSE(22)
 !STOP
 !
 !err=scprod(rhs-epre,rhs-epre)/rhsnorm
 !CALL maxsearch(epre,rhs)
 !print*, err

 !xx=bmid
  xx = (0.0D0, 0.0D0)
  CALL mulmat(xx, rr, tempr)
  pp = - rhs - rr
  rr = pp

  dprr = REAL(scprod(rr,rr))

  !=========================!
  ! Doing the main CG loop. !
  !=========================!

  DO ii=1,numiter

    CALL mulmat(pp, ap, tempr)
    alpha = dprr/scprod(pp, ap)
    xx = xx + alpha*pp
    rr = rr - alpha*ap
    dprrn = REAL(scprod(rr, rr))
    beta = dprrn/dprr
    pp = rr + beta*pp
    dprr = dprrn

    OPEN (33, file='errs.dat')

    IF ((ii/1)*1==ii) THEN 
    err = SQRT(dprr)/SQRT(rhsnorm)
    errsol = SQRT(solerr(xx, bmid))/SQRT(solnorm)
    print *, ii, err, errsol
    WRITE (33,*) ii, err, errsol
  
    IF(err<accur) EXIT
    END IF
   !print *, ii

  END DO

  CLOSE(33)

 !OPEN(22, file='test.dat')
 !
 !DO i=1,nptx
 !
 ! DO j=1,npty
 !
 !  DO k=1,nptz
 !   !DO l=1,nrotrot
 !   DO l=1,10
 !    WRITE(22,'(4i5,4e15.5)'),i,j,k,l,xx(l,i,j,k),bmid(l,i,j,k)
 !    rmxcur=ABS(xx(l,i,j,k)-bmid(l,i,j,k))
 !
 !    IF(rmxcur>rmx) THEN
 !     rmx=rmxcur
 !     imx=i
 !     jmx=j
 !     kmx=k
 !     lmx=l
 !    END IF
 !
 !   END DO
 !   
 !  END DO
 ! END DO
 !END DO
 !
 !CLOSE(22)
 !
 !PRINT *,imx,jmx,kmx,lmx,xx(lmx,imx,jmx,kmx),bmid(lmx,imx,jmx,kmx)

END PROGRAM efield