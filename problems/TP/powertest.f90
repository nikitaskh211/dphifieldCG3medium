! Log of changes:
! (MAY-04-2026) - Coppied base structure from "plssantenna.f90" file.
! (MAY-05-2026) - Implemented the power test in "powertest" program.

! TP.1.1. Module that defines constant parameters of a single-strap antenna.
MODULE ssantenna

   USE numberformat

   IMPLICIT NONE

   REAL(rk) :: ia       ! The total antenna's current,
   REAL(rk) :: xa, za   ! Central coordinates of the antenna on x- and z- coordinate,
   REAL(rk) :: dax, daz ! Half-widths of the antenna in x- and z- coordinate direction.

END MODULE ssantenna

! TP.1.2. A subroutine that reads antanna's parameters from a file.
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

! TP.1.3. A subroutine that sets the "jext" array.
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

! TP.2.1. Module that defines constant parameter for periodically inhomogeneous plasma.
MODULE piplasma

   USE numberformat

   IMPLICIT NONE

   COMPLEX(rk) :: cnseps0 ! A constant parameter (namely epsilon0),
   REAL(rk)    :: intz    ! The location of the interface on z coordinate.

END MODULE piplasma

! TP.2.2. A subroutine that reads plasma's parameters from a file.
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

! TP.2.3. A subroutine that sets the "zeta" array.
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

! TP.3.1. A function that calculates the injected power.
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

! TP.3.2. A function that calculates the absorbed power.
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

! TP.3.3. A function that calculates the absorbed power.
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

! TP.4. A program that solves for fields in the current problem.
PROGRAM powertest

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

  ! Field variables.
  COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: fields

  ! Pre-calculated power.
  REAL(rk) :: injpow0, abspow0, flxpow0

  !=========================!
  ! Setting up the program. !
  !=========================!

  ! Fetching pre-requisites.
  CALL datain
  CALL antennain
  CALL piplasmain

  ! Allocating memory for physical variables.
  ALLOCATE(jext(3, nptx, npty, nptz))
  ALLOCATE(zeta(3, 3, nptx, npty, nptz))

  ! Initialization of jext and zeta.
  CALL setantenna(jext)
  CALL setzeta(zeta)

  ! Printing physical jext and zeta into a file.
  CALL write3dvector(jext, "jext_antenna.csv")
  CALL write3x3dtensor(zeta, "zeta_plasma.csv")

  ! Allocating memory for field variables.
  ALLOCATE(fields(10, nptx, npty, nptz))

  ! Assigning constant values: D = (0.5 - 0.5i, - 0.4, 0.0), H = (0.0, 0.3, 1.0i)
  fields = (0.0D0, 0.0D0)
  fields(1, :, :, :) = (0.5D0, - 0.5D0)
  fields(2, :, :, :) = (- 0.4D0, 0.0D0)
  fields(6, :, :, :) = (0.3D0, 0.0D0)
  fields(7, :, :, :) = (0.0D0, 1.0D0)

  ! Calculating power analytically.
  injpow0 = REAL(fields(2, 1, 1, 1))*ia*sizey
  abspow0 = 2.0D0*lightspeed*sizex*sizey*(sizez - intz)*SUM(ABS(Fields(1:2, 1, 1, 1))**2)*1.22233991D-4/(size0*SQRT(pi**3))
  flxpow0 = - sizex*sizey*lightspeed*REAL(CONJG(fields(1, 1, 1, 1))*fields(6, 1, 1, 1))/(4.0D0*pi)

  ! Calculating powers.
  PRINT *, "[INF] Difference in injected power = ", injpow(zeta, fields, jext) - injpow0
  PRINT *, "[INF] Difference in absorbed power = ", abspow(zeta, fields) - abspow0
  PRINT *, "[INF] Difference in outflux  power = ", flxpow(zeta, fields) - flxpow0

  !=======================================!
  ! Memory cleanup for the whole program. !
  !=======================================!
  DEALLOCATE(jext, zeta, fields)

END PROGRAM powertest