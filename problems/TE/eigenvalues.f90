! Log of changes:
! (APR-26-2026) - Copied base structure from "plssantenna.f90" file.
! (APR-26-2026) - Removed all code related to antenna.
! (APR-26-2026) - Explicitly set Jext = 0, and solution = (0.0, 0.1)

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
    DO j = 1, npty
      DO k = 1, nptz

        injpow = injpow + REAL(DOT_PRODUCT(MATMUL(zeta(:, :, i, j, k), xx(1:3, i, j, k)), jext(1:3, i, j, k)))

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
    DO j = 1, npty
      DO k = 1, nptz

        abspow = abspow - AIMAG(DOT_PRODUCT(MATMUL(zeta(:, :, i, j, k), xx(1:3, i, j, k)), xx(1:3, i, j, k)))

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
    DO j = 1, npty

      evecc = CONJG(MATMUL(zeta(:, :, i, j, k), xx(1:3, i, j, k)))

      flxpow = flxpow - REAL(evecc(1)*xx(6, i, j, k) - evecc(2)*xx(5, i, j, k))

    END DO
  END DO

  flxpow = cns*flxpow*hx*hy

END FUNCTION flxpow

! P2.4. A program that solves for fields in the current problem.
PROGRAM dfield_antenna

  USE numberformat
  USE indata
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
  REAL(rk)                                        :: oldresidualnorm, newresidualnorm, conjugacy
  COMPLEX(rk)                                     :: step
  INTEGER                                         :: iter, maxiter

  !===========================!
  ! Calculation of constants. !
  !===========================!

  ! Setting the maximum amount of iterations.
  maxiter = 1000

  !=========================!
  ! Setting up the program. !
  !=========================!

  ! Fetching pre-requisites.
  CALL datain
  CALL piplasmain

  ! Allocating memory for physical variables.
  ALLOCATE(jext(3, nptx, npty, nptz))
  ALLOCATE(zeta(3, 3, nptx, npty, nptz))

  ! Initialization of jext and zeta.
 !CALL setantenna(jext)
  jext = (0.0D0)
  CALL setzeta(zeta)

  ! Printing physical jext and zeta into a file.
  CALL write3dvector(jext, "jext_antenna.csv")
  CALL write3x3dtensor(zeta, "zeta_plasma.csv")

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
 !rhs = (0.0D0)

  ! Initial guess and step of CG.
  solution = (1.0D0, 0.0D0)
  CALL mulmat(solution, Ux, temporary)
  CALL muladd(zeta, solution, Cx, temporary)
  Lx = Ux + Cx
  residual = rhs - Lx
  oldresidualnorm = REAL(scprod(residual, residual))
  IF (oldresidualnorm < accur) PRINT *, "[INF] Normalized residual norm is sufficiently small."
  searchdr = residual

  ! Main loop.
  PRINT *, "[INF] Information in format : (iteration), (oldresidualnorm), (resdivd), (resdivh), (resdiva), (rescurlsrc)"
  OPEN(19, FILE = "errors.csv", STATUS = 'replace')
  WRITE(19, *) "Iter, ResNorm, ResDivD, ResDivH, ResDivA, ResCurlSrc"
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
    IF (oldresidualnorm < accur) THEN
      PRINT *, "[INF] Normalized residual norm is sufficiently small. Solution converged to an accuracy = ", accur, ", at iteration = ", iter, "."
      EXIT
    END IF

    ! Printing informaiton about the current situation.
    PRINT *, "[INF]", iter, oldresidualnorm, resdivd(solution, jext), resdivh(solution), resdiva(solution), rescurlsrc(solution, jext)
    WRITE(19, "(I0, 4(ES24.16), ES24.16)") iter, oldresidualnorm, resdivd(solution, jext), resdivh(solution), resdiva(solution), rescurlsrc(solution, jext)

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