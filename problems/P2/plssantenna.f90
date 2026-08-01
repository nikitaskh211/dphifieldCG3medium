! Log of changes:
! (APR-24-2026) - Introduced a module that contains all the necessary information about a single-strap antenna with a subroutine to get input values.
! (APR-24-2026) - Introduced a module that contains all the necessary information about a periodically inhomogeneous plasma with a subroutine to get input values.
! (APR-24-2026) - Implemented the CG algorithm in a program.
! (APR-24-2026) - Implemented functions for calculating injected and absorbed power : "injpow" and "abspow" respectfully.
! (APR-24-2026) - Implemented a function for calculating flux in power : "flxpow".
! (APR-26-2026) - Increased the upper limit for iterations from 1000 to 10000.
! (APR-27-2026) - CG stops when secondary residuals converge to a single value.

! P2.1.1. Module that defines constant parameters of a single-strap antenna.
MODULE ssantenna

   USE numberformat

   IMPLICIT NONE

   REAL(rk) :: ia       ! The total antenna's current,
   REAL(rk) :: xa, za   ! Central coordinates of the antenna on x- and z- coordinate,
   REAL(rk) :: dax, daz ! Half-widths of the antenna in x- and z- coordinate direction.

END MODULE ssantenna

! P2.1.2. A subroutine that reads antanna's parameters from a file.
SUBROUTINE antennain

  USE ssantenna

  IMPLICIT NONE

  OPEN (12, FILE = 'ssantenna_in.dat', STATUS = 'old')
  OPEN (13, FILE = 'ssantenna_out.dat')

  READ(12,*)  ia
  WRITE(13,*) ia

  READ(12,*)  xa, za
  WRITE(13,*) xa, za

  READ(12,*)  dax, daz
  WRITE(13,*) dax, daz

  CLOSE(12)
  CLOSE(13)

END SUBROUTINE antennain

! P2.1.3. A subroutine that sets the "jext" array.
SUBROUTINE setantenna(jext)

  USE numberformat
  USE constants
  USE indata
  USE ssantenna

  IMPLICIT NONE

  COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: jext ! Output array.

  INTEGER :: i, k ! Iterator variables for x and z.

  REAL(rk) :: hx, hz ! Spacings for x and z.

  REAL(rk) :: argx, argz ! Actual values of x and z.
  REAL(rk) :: cosx, cosz ! Actual values of functions of x and z.

  !===========================!
  ! Calculation of constants. !
  !===========================!

  hx = sizex/(nptx - 1.0D0)
  hz = sizez/(nptz - 1.0D0)

  !===================!
  ! Array assignment. !
  !===================!

  ! Double looped summation.
  DO i = 1, nptx

    ! Calculating argx and cosx.
    argx = (i - 1.0D0)*hx
    cosx = 0.0D0
    IF(ABS(argx - xa) < dax) cosx = COS((pi*(argx - xa))/(2.0D0*dax))**2

    DO k = 1, nptz

      ! Calculating argz and cosz.
      argz = (k - 1.0D0)*hz
      cosz = 0.0D0
      IF(ABS(argz - za) < daz) cosz = COS((pi*(argz - za))/(2.0D0*daz))**2

      ! Jextx = 0
      jext(1, i, :, k) = (0.0D0, 0.0D0)

      ! Jexty = Ia.cos(pi.(x - xa)/(2.dax)).cos(pi.(z - za)/(2.daz))/(dax*daz)
      jext(2, i, :, k) = (ia*cosx*cosz)/(dax*daz)

      ! Jextz = 0
      jext(3, i, :, k) = (0.0D0, 0.0D0)

    END DO
  END DO

END SUBROUTINE setantenna

! P2.1.4. A subroutine for calculating the right-hand-side of positive definite equations using analytical expressions for curl(Jext),
! bypassing all finite differences. Applies the same boundary conditions for Jext as setrhs.
! Note: div(Jext) = 0 analytically since Jext = Jy(x,z)*ey only, so grad(div[Jext]) vanishes exactly.
SUBROUTINE setrhs_analytical(rhs, jext)

  USE numberformat
  USE constants
  USE indata
  USE ssantenna

  IMPLICIT NONE

  COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(OUT) :: rhs  ! Output array,
  COMPLEX(rk), DIMENSION(3, nptx, npty, nptz),  INTENT(OUT) :: jext ! Output array.

  COMPLEX(rk) :: ik0 ! Wavenumber of the free space multiplied by an imaginary unit.
  REAL(rk)    :: cns ! A constant factor: cns = 4.pi/c

  INTEGER :: i, j, k ! Iterator variables for x, y and z.

  REAL(rk) :: hx, hz ! Spacings for x and z.

  REAL(rk) :: argx, argz ! Actual values of x and z.

  REAL(rk) :: fx,  fz  ! Actual values of cos^2 functions of x and z.
  REAL(rk) :: dfx, dfz ! Actual values of analytical derivatives of cos^2 functions of x and z.

  REAL(rk) :: prefac ! A prefactor: prefac = Ia/(dax*daz).

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Wavevector.
  ik0 = (0.0D0, 1.0D0)*(pi/size0)

  ! Constant factor.
  cns = 4.0D0*pi/lightspeed

  ! Prefactor.
  prefac = ia/(dax*daz)

  hx = sizex/(nptx - 1.0D0)
  hz = sizez/(nptz - 1.0D0)

  !===================!
  ! Array assignment. !
  !===================!

  jext = (0.0D0, 0.0D0)
  rhs  = (0.0D0, 0.0D0)

  ! Triple looped assignment.
  DO k = 1, nptz

    ! Calculating argz, fz and dfz.
    argz = (k - 1.0D0)*hz
    fz  = 0.0D0
    dfz = 0.0D0
    IF(ABS(argz - za) < daz) THEN
      fz  =   COS((pi*(argz - za))/(2.0D0*daz))**2
      dfz = - SIN((pi*(argz - za))/daz)*(pi/(2.0D0*daz))
    END IF

    DO i = 1, nptx

      ! Calculating argx, fx and dfx.
      argx = (i - 1.0D0)*hx
      fx  = 0.0D0
      dfx = 0.0D0
      IF(ABS(argx - xa) < dax) THEN
        fx  =   COS((pi*(argx - xa))/(2.0D0*dax))**2
        dfx = - SIN((pi*(argx - xa))/dax)*(pi/(2.0D0*dax))
      END IF

      DO j = 1, npty

        !===========================================================================!
        ! A set of assignments to calculate jext and rhs for the antenna's current. !
        !===========================================================================!

        ! Jextx = 0, Jextz = 0.
        jext(1, i, j, k) = (0.0D0, 0.0D0)
        jext(3, i, j, k) = (0.0D0, 0.0D0)

        ! Jexty = Ia.cos^2(pi.(x - xa)/(2.dax)).cos^2(pi.(z - za)/(2.daz))/(dax*daz).
        jext(2, i, j, k) = prefac*fx*fz

        !========================================================================!
        ! A set of assignments to calculate rhs to the governing equation for D. !
        !========================================================================!

        ! R(Jext){1:3} = - 4.pi.(i.k0.Jext + grad[div{Jext}]/[i.k0])/c = - 4.pi.i.k0.Jext/c,
        ! since div(Jext) = d/dy(Jy) = 0 exactly (Jext has no y-dependence).
        rhs(1, i, j, k) = - cns*ik0*jext(1, i, j, k) ! Dx-component: Jextx = 0.
        rhs(2, i, j, k) = - cns*ik0*jext(2, i, j, k) ! Dy-component: Jexty = Ia.fx.fz/(dax*daz).
        rhs(3, i, j, k) = - cns*ik0*jext(3, i, j, k) ! Dz-component: Jextz = 0.

        !==========================================================================!
        ! A set of assignments to calculate rhs to the governing equation for Phi. !
        !==========================================================================!

        rhs(4, i, j, k) = (0.0D0, 0.0D0) ! R(Jext){4} = 0.

        !========================================================================!
        ! A set of assignments to calculate rhs to the governing equation for H. !
        !========================================================================!

        ! R(Jext){5:7} = 4.pi.curl(Jext)/c, evaluated analytically.
        ! curl(Jext)_x = - d/dz(Jy) = - prefac.fx.dfz
        ! curl(Jext)_y =   0 (d/dx[Jx] - d/dz[Jz] = 0)
        ! curl(Jext)_z = + d/dx(Jy) = + prefac.dfx.fz
        rhs(5, i, j, k) = cns*(- prefac*fx*dfz)  ! Hx-component.
        rhs(6, i, j, k) = (0.0D0, 0.0D0)         ! Hy-component: curl(Jext)_y = 0.
        rhs(7, i, j, k) = cns*(+ prefac*dfx*fz)  ! Hz-component.

        !========================================================================!
        ! A set of assignments to calculate rhs to the governing equation for A. !
        !========================================================================!

        rhs(8,  i, j, k) = (0.0D0, 0.0D0)
        rhs(9,  i, j, k) = (0.0D0, 0.0D0)
        rhs(10, i, j, k) = (0.0D0, 0.0D0)

      END DO
    END DO
  END DO

  ! Neumann's BC's.
  jext(1, 1, :, :) = jext(1, 2, :, :)           ! On x = 0 jext,x's derivative = 0
  jext(1, nptx, :, :) = jext(1, nptx - 1, :, :) ! On x = Lx jext,x's derivative = 0

  jext(2, :, 1, :) = jext(2, :, 2, :)           ! On y = 0 jext,y's derivative = 0
  jext(2, :, npty, :) = jext(2, :, npty - 1, :) ! On y = Ly jext,y's derivative = 0

  jext(3, :, :, 1) = jext(3, :, :, 2)           ! On z = 0 jext,z's derivative = 0
  jext(3, :, :, nptz) = jext(3, :, :, nptz - 1) ! On z = Lz jext,z's derivative = 0

  rhs(1, 1, :, :) = rhs(1, 2, :, :)           ! On x = 0 jext,x's derivative = 0
  rhs(1, nptx, :, :) = rhs(1, nptx - 1, :, :) ! On x = Lx jext,x's derivative = 0

  rhs(2, :, 1, :) = rhs(2, :, 2, :)           ! On y = 0 jext,y's derivative = 0
  rhs(2, :, npty, :) = rhs(2, :, npty - 1, :) ! On y = Ly jext,y's derivative = 0

  rhs(3, :, :, 1) = rhs(3, :, :, 2)           ! On z = 0 jext,z's derivative = 0
  rhs(3, :, :, nptz) = rhs(3, :, :, nptz - 1) ! On z = Lz jext,z's derivative = 0

  ! Dirichlet BC's.
  jext(1, :, 1, :) = (0.0D0, 0.0D0)
  jext(1, :, :, 1) = (0.0D0, 0.0D0)
  jext(1, :, npty, :) = (0.0D0, 0.0D0)
  jext(1, :, :, nptz) = (0.0D0, 0.0D0)
  jext(2, 1, :, :) = (0.0D0, 0.0D0)
  jext(2, :, :, 1) = (0.0D0, 0.0D0)
  jext(2, nptx, :, :) = (0.0D0, 0.0D0)
  jext(2, :, :, nptz) = (0.0D0, 0.0D0)
  jext(3, 1, :, :) = (0.0D0, 0.0D0)
  jext(3, :, 1, :) = (0.0D0, 0.0D0)
  jext(3, nptx, :, :) = (0.0D0, 0.0D0)
  jext(3, :, npty, :) = (0.0D0, 0.0D0)
  rhs(1, :, 1, :) = (0.0D0, 0.0D0)
  rhs(1, :, :, 1) = (0.0D0, 0.0D0)
  rhs(1, :, npty, :) = (0.0D0, 0.0D0)
  rhs(1, :, :, nptz) = (0.0D0, 0.0D0)
  rhs(2, 1, :, :) = (0.0D0, 0.0D0)
  rhs(2, :, :, 1) = (0.0D0, 0.0D0)
  rhs(2, nptx, :, :) = (0.0D0, 0.0D0)
  rhs(2, :, :, nptz) = (0.0D0, 0.0D0)
  rhs(3, 1, :, :) = (0.0D0, 0.0D0)
  rhs(3, :, 1, :) = (0.0D0, 0.0D0)
  rhs(3, nptx, :, :) = (0.0D0, 0.0D0)
  rhs(3, :, npty, :) = (0.0D0, 0.0D0)

END SUBROUTINE setrhs_analytical

! P2.2.1. Module that defines constant parameter for periodically inhomogeneous plasma.
MODULE piplasma

   USE numberformat

   IMPLICIT NONE

   COMPLEX(rk) :: cnseps0 ! A constant parameter (namely epsilon0),
   REAL(rk)    :: intz    ! The location of the interface on z coordinate.

END MODULE piplasma

! P2.2.2. A subroutine that reads plasma's parameters from a file.
SUBROUTINE piplasmain

  USE numberformat
  USE piplasma

  IMPLICIT NONE

  REAL(rk) :: re_cnseps0, im_cnseps0

  OPEN (14, FILE = 'piplasma_in.dat', STATUS = 'old')
  OPEN (15, FILE = 'piplasma_out.dat')

  READ(14,*)  re_cnseps0, im_cnseps0
  WRITE(15,*) re_cnseps0, im_cnseps0

  cnseps0 = re_cnseps0 + (0.0D0, 1.0D0)*im_cnseps0

  READ(14,*)  intz
  WRITE(15,*) intz

  CLOSE(14)
  CLOSE(15)

END SUBROUTINE piplasmain

! P2.2.3. A subroutine that sets the "zeta" array.
SUBROUTINE setzeta(zeta)

  USE numberformat
  USE constants
  USE indata
  USE piplasma

  IMPLICIT NONE

  COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(OUT) :: zeta ! Output array.

  INTEGER :: i, j, k ! Iterator variables for x, y, and z.

  REAL(rk) :: hx, hy, hz ! Spacings for x, y and z.

  REAL(rk) :: kx, ky, kz ! Wavenumbers.

  REAL(rk) :: argx, argy, argz ! Actual values of x, y, and z.
  REAL(rk) :: sinx, siny, sinz ! Actual values of functions of x, y, and z.

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Calculation of node spacings.
  hx = sizex/(nptx - 1.0D0)
  hy = sizey/(npty - 1.0D0)
  hz = sizez/(nptz - 1.0D0)

  ! Wavevector.
  kx = pi/sizex
  ky = pi/sizey
  kz = pi/(sizez - intz)

  !===================!
  ! Array assignment. !
  !===================!

  DO i = 1, nptx

    ! Calculating argx and cosx.
    argx = (i - 1.0D0)*hx
    sinx = SIN(kx*argx)

    DO j = 1, npty

      ! Calculating argy and cosy.
      argy = (j - 1.0D0)*hy
      siny = SIN(ky*argy)

      DO k = 1, nptz

        ! Calculating argz and cosz.
        argz = (k - 1.0D0)*hz
        sinz = 0.0D0
        IF (argz > intz) sinz = SIN(kz*(argz - intz))
    
        ! ZetaXX = 1/(1 - cnseps0.sin[kx.x].sin[ky.y].sin[kz.{z - intz}])
        zeta(1, 1, i, j, k) = 1.0D0/(1.0D0 - cnseps0*sinx*siny*sinz)
        
        ! ZetaXY = 0
        zeta(1, 2, i, j, k) = (0.0D0, 0.0D0)

        ! ZetaXZ = 0
        zeta(1, 3, i, j, k) = (0.0D0, 0.0D0)
    
        ! ZetaYX = 0
        zeta(2, 1, i, j, k) = (0.0D0, 0.0D0)
        
        ! ZetaYY = 1/(1 - cnseps0.sin[kx.x].sin[ky.y].sin[kz.{z - intz}])
        zeta(2, 2, i, j, k) = 1.0D0/(1.0D0 - cnseps0*sinx*siny*sinz)

        ! ZetaYZ = 0
        zeta(2, 3, i, j, k) = (0.0D0, 0.0D0)
    
        ! ZetaZX = 0
        zeta(3, 1, i, j, k) = (0.0D0, 0.0D0)
        
        ! ZetaZY = 0
        zeta(3, 2, i, j, k) = (0.0D0, 0.0D0)

        ! ZetaZZ = 1/(1 - cnseps0.sin[kx.x].sin[ky.y].sin[kz.{z - intz}])
        zeta(3, 3, i, j, k) = 1.0D0/(1.0D0 - cnseps0*sinx*siny*sinz)

      END DO
    END DO
  END DO

END SUBROUTINE setzeta

! P2.3.1. A function that calculates the injected power.
FUNCTION injpow(zeta, xx, jext)

  USE numberformat
  USE indata

  IMPLICIT NONE

  COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN) :: zeta ! Input array,
  COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN)   :: xx   ! Input array,
  COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)    :: jext ! Input array.

  REAL(rk) :: injpow

  REAL(rk) :: hx, hy, hz ! Spacings for x, y and z.
  REAL(rk) :: wx, wy, wz ! Quadrature.

  INTEGER :: i, j, k ! Iterator variables for x, y, and z.

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Calculation of node spacings.
  hx = sizex/(nptx - 1.0D0)
  hy = sizey/(npty - 1.0D0)
  hz = sizez/(nptz - 1.0D0)

  !=========================!
  ! Accumulative summation. !
  !=========================!

  injpow = 0.0D0

  DO i = 1, nptx

    wx = 1.0D0
    IF (i == 1 .OR. i == nptx) wx = 0.5D0

    DO j = 1, npty

      wy = 1.0D0
      IF (j == 1 .OR. j == npty) wy = 0.5D0

      DO k = 1, nptz

        wz = 1.0D0
        IF (k == 1 .OR. k == nptz) wz = 0.5D0

        injpow = injpow + wx*wy*wz*REAL(DOT_PRODUCT(MATMUL(zeta(:, :, i, j, k), xx(1:3, i, j, k)), jext(1:3, i, j, k)))

      END DO
    END DO
  END DO

  injpow = injpow*hx*hy*hz

END FUNCTION injpow

! P2.3.2. A function that calculates the absorbed power.
FUNCTION abspow(zeta, xx)

  USE numberformat
  USE constants
  USE indata

  IMPLICIT NONE

  COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN) :: zeta ! Input array,
  COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN)   :: xx   ! Input array,

  REAL(rk) :: abspow

  REAL(rk) :: hx, hy, hz ! Spacings for x, y and z.
  REAL(rk) :: cns        ! Constant factor : cns = w/(4.pi)
  REAL(rk) :: wx, wy, wz ! Quadrature.

  INTEGER :: i, j, k ! Iterator variables for x, y, and z.

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Calculation of node spacings.
  hx = sizex/(nptx - 1.0D0)
  hy = sizey/(npty - 1.0D0)
  hz = sizez/(nptz - 1.0D0)

  ! Constant factor.
  cns = 0.25D0*lightspeed/size0

  !=========================!
  ! Accumulative summation. !
  !=========================!

  abspow = 0.0D0

  DO i = 1, nptx

    wx = 1.0D0
    IF (i == 1 .OR. i == nptx) wx = 0.5D0

    DO j = 1, npty

      wy = 1.0D0
      IF (j == 1 .OR. j == npty) wy = 0.5D0

      DO k = 1, nptz

        wz = 1.0D0
        IF (k == 1 .OR. k == nptz) wz = 0.5D0

        abspow = abspow - wx*wy*wz*AIMAG(DOT_PRODUCT(MATMUL(zeta(:, :, i, j, k), xx(1:3, i, j, k)), xx(1:3, i, j, k)))

      END DO
    END DO
  END DO

  abspow = cns*abspow*hx*hy*hz

END FUNCTION abspow

! P2.3.3. A function that calculates the absorbed power.
FUNCTION flxpow(zeta, xx)

  USE numberformat
  USE constants
  USE indata
  USE piplasma

  IMPLICIT NONE

  COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN) :: zeta ! Input array,
  COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN)   :: xx   ! Input array,

  REAL(rk) :: flxpow

  REAL(rk) :: hx, hy, hz ! Spacings for x, y and z.
  REAL(rk) :: wx, wy     ! Quadrature.
  REAL(rk) :: cns        ! Constant factor : cns = c/(4.pi)

  COMPLEX(rk), DIMENSION(3) :: evecc ! Temporary value of the conjugated electric field strength at node i, j, k.

  INTEGER :: i, j, k ! Iterator variables for x, y, and z.

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Calculation of node spacings.
  hx = sizex/(nptx - 1.0D0)
  hy = sizey/(npty - 1.0D0)
  hz = sizez/(nptz - 1.0D0)

  ! Constant factor.
  cns = lightspeed/(4.0D0*pi)

  ! Interface node.
  k = NINT(intz/hz)
  k = MAX(1, MIN(nptz, k))

  !=========================!
  ! Accumulative summation. !
  !=========================!

  flxpow = 0.0D0

  DO i = 1, nptx

    wx = 1.0D0
    IF (i == 1 .OR. i == nptx) wx = 0.5D0

    DO j = 1, npty

      wy = 1.0D0
      IF (j == 1 .OR. j == npty) wy = 0.5D0

      evecc = CONJG(MATMUL(zeta(:, :, i, j, k), xx(1:3, i, j, k)))

      flxpow = flxpow - wx*wy*REAL(evecc(1)*xx(6, i, j, k) - evecc(2)*xx(5, i, j, k))

    END DO
  END DO

  flxpow = cns*flxpow*hx*hy

END FUNCTION flxpow

! P2.4. A program that solves for fields in the current problem.
PROGRAM dfield_antenna

  USE numberformat
  USE constants
  USE indata
  USE ssantenna
  USE piplasma

  IMPLICIT NONE

  !======================!
  ! Declaring variables. !
  !======================!

  ! Physical variables.
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE    :: jext ! The electric current density.
  COMPLEX(rk), DIMENSION(:, :, :, :, :), ALLOCATABLE :: zeta ! The inverse of dielectric medium permittivity tensor.

  ! External functions.
  COMPLEX(rk), EXTERNAL :: scprod
  REAL(rk), EXTERNAL    :: resdivd, resdivh, resdiva, rescurlsrc
  REAL(rk), EXTERNAL    :: injpow, abspow, flxpow

  ! CG algorithm variables.
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: rhs, solution, residual, searchdr, temporary, Ux, Cx, Lx
  REAL(rk)                                        :: rhsnorm, oldresidualnorm, newresidualnorm, conjugacy
  COMPLEX(rk)                                     :: step
  INTEGER                                         :: iter, maxiter

  ! Physical stagnation.
  REAL(rk) :: currresdivd, currresdivh, currresdiva, currrescurlsrc
  REAL(rk) :: prevresdivd, prevresdivh, prevresdiva, prevrescurlsrc
  REAL(rk) :: physfloor, hx, hy, hz
  INTEGER  :: checkinterval

  !=========================!
  ! Setting up the program. !
  !=========================!

  ! Fetching pre-requisites.
  CALL datain
  CALL antennain
  CALL piplasmain

  ! Setting the maximum amount of iterations.
  maxiter = 10*nptx*npty*nptz ! Analytically this should be the hard limit.

  ! Allocating memory for physical variables.
  ALLOCATE(jext(3, nptx, npty, nptz))
  ALLOCATE(zeta(3, 3, nptx, npty, nptz))

  ! Initialization of jext and zeta.
  CALL setantenna(jext)
  CALL setzeta(zeta)

  ! Printing physical jext and zeta into a file.
  CALL write3dvector(jext, "jext_antenna.csv")
  CALL write3x3dtensor(zeta, "zeta_plasma.csv")

  !==================================!
  ! Stagnation of the physical laws. !
  !==================================!

  ! Discretization floor : O(k0*h)^2, take worst direction.
  hx = sizex/(nptx - 1.0D0)
  hy = sizey/(npty - 1.0D0)
  hz = sizez/(nptz - 1.0D0)
  physfloor = (pi*MAX(hx, hy, hz)/size0)**2

  ! Automatic check interval based on condition number estimate.
  checkinterval = MAX(100, INT(SQRT(1.0D0/physfloor) * &
                  SQRT(MAXVAL(ABS(zeta))/MINVAL(ABS(zeta)))))

  PRINT *, "[INF] Physical residual floor     = ", physfloor
  PRINT *, "[INF] Physical stagnation interval = ", checkinterval

  ! Initialization of previous residuals.
  prevresdivd    = HUGE(1.0D0)
  prevresdivh    = HUGE(1.0D0)
  prevresdiva    = HUGE(1.0D0)
  prevrescurlsrc = HUGE(1.0D0)

  !===============!
  ! CG algorithm. !
  !===============!

  ! Allocating memory for CG variables.
  ALLOCATE(rhs(10, nptx, npty, nptz))
  ALLOCATE(solution(10, nptx, npty, nptz))
  ALLOCATE(residual(10, nptx, npty, nptz))
  ALLOCATE(searchdr(10, nptx, npty, nptz))
  ALLOCATE(temporary(3, nptx, npty, nptz))
  ALLOCATE(Ux(10, nptx, npty, nptz))
  ALLOCATE(Cx(10, nptx, npty, nptz))
  ALLOCATE(Lx(10, nptx, npty, nptz))

  ! Setting the right-hand-side for CG computation.
  CALL setrhs(rhs, jext, temporary)
  !CALL setrhs_analytical(rhs, jext)
  rhsnorm = SQRT(REAL(scprod(rhs, rhs)))

  ! Initial guess and step of CG.
  solution = (0.0D0, 0.0D0)
  CALL mulmat(solution, Ux, temporary)
  CALL muladd(zeta, solution, Cx, temporary)
  Lx = Ux + Cx
  residual = rhs - Lx
  oldresidualnorm = REAL(scprod(residual, residual))
  IF (SQRT(oldresidualnorm) < accur) PRINT *, "[INF] Normalized residual norm is sufficiently small."
  searchdr = residual

  ! Main loop.
  PRINT *, "[INF] Information in format : (iteration), (oldresidualnorm/rhsnorm), (resdivd), (resdivh), (resdiva), (rescurlsrc)"
  OPEN(19, FILE = "errors.csv", STATUS = 'replace')
  WRITE(19, *) "Iter, NormResNorm, ResDivD, ResDivH, ResDivA, ResCurlSrc"
  DO iter = 1, maxiter

    ! Main steps of CG algorithm.
    CALL mulmat(searchdr, Ux, temporary)
    CALL muladd(zeta, searchdr, Cx, temporary)
    Lx = Ux + Cx
    step = oldresidualnorm/scprod(searchdr, Lx)
    solution = solution + step*searchdr
    residual = residual - step*Lx
    newresidualnorm = REAL(scprod(residual, residual))
    conjugacy = newresidualnorm/oldresidualnorm
    searchdr = residual + conjugacy*searchdr
    oldresidualnorm = newresidualnorm

    ! Calculating residuals of the physical laws.
    currresdivd    = resdivd(solution, jext)
    currresdivh    = resdivh(solution)
    currresdiva    = resdiva(solution)
    currrescurlsrc = rescurlsrc(solution, jext)

    ! Printing information about the current situation.
    PRINT *, "[INF]", iter, SQRT(oldresidualnorm)/rhsnorm, currresdivd, currresdivh, currresdiva, currrescurlsrc
    WRITE(19, "(I0, 5(ES24.16))") iter, SQRT(oldresidualnorm)/rhsnorm, currresdivd, currresdivh, currresdiva, currrescurlsrc

    ! CG convergence check.
    IF (SQRT(oldresidualnorm)/rhsnorm < accur) THEN
      PRINT *, "[INF] CG residual converged at iteration = ", iter
      EXIT
    END IF

    ! Physical stagnation check every checkinterval iterations.
    IF (MOD(iter, checkinterval) == 0) THEN

      ! Check if residuals have reached the discretization floor.
      IF (currresdivd    < physfloor .AND. &
          currresdivh    < physfloor .AND. &
          currresdiva    < physfloor .AND. &
          currrescurlsrc < physfloor) THEN
        PRINT *, "[INF] Physical residuals reached discretization floor at iteration = ", iter
        PRINT *, "[INF] Floor    = ", physfloor
        PRINT *, "[INF] Achieved = ", currresdivd, currresdivh, currresdiva, currrescurlsrc
        EXIT
      END IF

      ! Check if residuals have stagnated above the floor.
      IF (ABS(currresdivd    - prevresdivd)    / (prevresdivd    + TINY(1.0D0)) < accur .AND. &
          ABS(currresdivh    - prevresdivh)    / (prevresdivh    + TINY(1.0D0)) < accur .AND. &
          ABS(currresdiva    - prevresdiva)    / (prevresdiva    + TINY(1.0D0)) < accur .AND. &
          ABS(currrescurlsrc - prevrescurlsrc) / (prevrescurlsrc + TINY(1.0D0)) < accur) THEN
        PRINT *, "[INF] Physical residuals stagnated above floor at iteration = ", iter
        PRINT *, "[INF] Achieved = ", currresdivd, currresdivh, currresdiva, currrescurlsrc
        PRINT *, "[INF] Floor    = ", physfloor
        EXIT
      END IF

      ! Update previous values only at checkpoints.
      prevresdivd    = currresdivd
      prevresdivh    = currresdivh
      prevresdiva    = currresdiva
      prevrescurlsrc = currrescurlsrc

    END IF

  END DO
  CLOSE(19)

  ! Calculating powers.
  PRINT *, "[INF] Injected power = ", injpow(zeta, solution, jext)
  PRINT *, "[INF] Absorbed power = ", abspow(zeta, solution)
  PRINT *, "[INF] Outflux  power = ", flxpow(zeta, solution)

  ! Writing physical solution of the system to a file.
  CALL write10dvector(solution, "solution_antenna_plasma.csv")

  !=======================================!
  ! Memory cleanup for the whole program. !
  !=======================================!
  DEALLOCATE(rhs, solution, residual, searchdr, temporary, Ux, Cx, Lx)

END PROGRAM dfield_antenna