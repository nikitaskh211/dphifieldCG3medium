! Log of changes:
! (MAY-11-2026) - Removed "T0" folder as it is deprecated to this project.
! (MAY-11-2026) - Copied plasma modules and subroutines from "P2" folder.
! (MAY-11-2026) - Added a program that performs the test for

! TL.1. Module that defines constant parameter for periodically inhomogeneous plasma.
MODULE piplasma

   USE numberformat

   IMPLICIT NONE

   COMPLEX(rk) :: cnseps0 ! A constant parameter (namely epsilon0),
   REAL(rk)    :: intz    ! The location of the interface on z coordinate.

END MODULE piplasma

! TL.2. A subroutine that reads plasma's parameters from a file.
SUBROUTINE piplasmain

   USE numberformat
   USE piplasma

   IMPLICIT NONE

   REAL(rk) :: re_cnseps0, im_cnseps0

   OPEN (14, FILE='piplasma_in.dat', STATUS='old')
   OPEN (15, FILE='piplasma_out.dat')

   READ (14, *) re_cnseps0, im_cnseps0
   WRITE (15, *) re_cnseps0, im_cnseps0

   cnseps0 = re_cnseps0 + (0.0D0, 1.0D0)*im_cnseps0

   READ (14, *) intz
   WRITE (15, *) intz

   CLOSE (14)
   CLOSE (15)

END SUBROUTINE piplasmain

! TL.3. A subroutine that sets the "zeta" array.
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

! TL.4. A subroutine that calculates global, local interior and outer residual and prints it in a certain format.
SUBROUTINE respreset(vec1, vec2)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: vec1, vec2

   REAL(rk) :: num, den, candidate

   REAL(rk) :: respresetg, maxrespreseti, maxrespreseto

   REAL(rk) :: wx, wy, wz ! Quadrature.

   INTEGER :: i, j, k
   INTEGER :: maxii = -1, maxij = -1, maxik = -1
   INTEGER :: maxoi = -1, maxoj = -1, maxok = -1

   num = 0.0D0
   den = 0.0D0
   candidate = 0.0D0

   respresetg = 0.0D0
   maxrespreseti = 0.0D0
   maxrespreseto = 0.0D0

   DO i = 1, nptx

      wx = 1.0D0
      IF (i == 1 .OR. i == nptx) wx = 0.5D0

      DO j = 1, npty

         wy = 1.0D0
         IF (j == 1 .OR. j == npty) wy = 0.5D0

         DO k = 1, nptz

            wz = 1.0D0
            IF (k == 1 .OR. k == nptz) wz = 0.5D0

            num = num + SUM(ABS(vec1(:, i, j, k) - vec2(:, i, j, k))**2)*wx*wy*wz
            den = den + SUM(ABS(vec1(:, i, j, k))**2)*wx*wy*wz

            candidate = SUM(ABS(vec1(:, i, j, k) - vec2(:, i, j, k))**2)*wx*wy*wz
            IF (i > 1 .AND. i < nptx .AND. j > 1 .AND. j < npty .AND. k > 1 .AND. k < nptz) THEN

               IF (maxrespreseti < candidate) THEN

                  maxrespreseti = candidate
                  maxii = i
                  maxij = j
                  maxik = k

               END IF

            ELSE

               IF (maxrespreseto < candidate) THEN

                  maxrespreseto = candidate
                  maxoi = i
                  maxoj = j
                  maxok = k

               END IF

            END IF

         END DO
      END DO
   END DO

   IF (den > 0.0D0) THEN
      respresetg = SQRT(num/den)
      maxrespreseti = SQRT(nptx*npty*nptz*maxrespreseti/den)
      maxrespreseto = SQRT(nptx*npty*nptz*maxrespreseto/den)
   ELSE
      respresetg = -1.0D0
      maxrespreseti = -1.0D0
      maxrespreseto = -1.0D0
   END IF

   PRINT *, "[DBG] Global residual = ", respresetg
   PRINT *, "[DBG] Interior local residual = ", maxrespreseti, ", at node i, j, k = ", maxii, maxij, maxik
   PRINT *, "[DBG] Outer local residual = ", maxrespreseto, ", at node i, j, k = ", maxoi, maxoj, maxok

END SUBROUTINE respreset

! TL.5. A subroutine that sets the gradient of "zeta".
SUBROUTINE setdivzeta0(res, dvec)

   USE numberformat
   USE constants
   USE indata
   USE piplasma

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(nptx, npty, nptz), INTENT(OUT)   :: res  ! Output array.
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN) :: dvec ! Input array.

   INTEGER :: i, j, k

   REAL(rk) :: argx, argy, argz

   DO i = 1, nptx
      argx = pi*(i - 1.0D0)/(nptx - 1.0D0)
      DO j = 1, npty
         argy = pi*(j - 1.0D0)/(npty - 1.0D0)
         DO k = 1, nptz
            argz = pi*(k - 1.0D0)/(nptz - 1.0D0)

            res(i, j, k) = cnseps0*pi*( &
                           dvec(1, i, j, k)*COS(argx)*SIN(argy)*SIN(argz)/sizex + &
                           dvec(2, i, j, k)*SIN(argx)*COS(argy)*SIN(argz)/sizey + &
                           dvec(3, i, j, k)*SIN(argx)*SIN(argy)*COS(argz)/sizez &
                           )/((1.0D0 - cnseps0*SIN(argx)*SIN(argy)*SIN(argz))**2)

         END DO
      END DO
   END DO

END SUBROUTINE setdivzeta0

! TL.6. A set of subroutines that sets x, y, or z linearly.
SUBROUTINE setlinxyz(res)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res

   INTEGER  :: i, j, k
   REAL(rk) :: hx, hy, hz
   REAL(rk) :: argx, argy, argz

   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   DO i = 1, nptx
      argx = hx*(i - 1.0D0)
      DO j = 1, npty
         argy = hy*(j - 1.0D0)
         DO k = 1, nptz
            argz = hz*(k - 1.0D0)

            res(1, i, j, k) = argx
            res(2, i, j, k) = argy
            res(3, i, j, k) = argz

         END DO
      END DO
   END DO

END SUBROUTINE

! TL.7. A subroutine that sets x^2, y^2, or z^2 quadratically.
SUBROUTINE setquadxyz(res)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res

   INTEGER  :: i, j, k
   REAL(rk) :: hx, hy, hz
   REAL(rk) :: argx, argy, argz

   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   DO i = 1, nptx
      argx = hx*(i - 1.0D0)
      DO j = 1, npty
         argy = hy*(j - 1.0D0)
         DO k = 1, nptz
            argz = hz*(k - 1.0D0)

            res(1, i, j, k) = argx**2
            res(2, i, j, k) = argy**2
            res(3, i, j, k) = argz**2

         END DO
      END DO
   END DO

END SUBROUTINE

PROGRAM test_l_operator

   USE numberformat
   USE constants
   USE indata
   USE piplasma

   IMPLICIT NONE

   !===========================!
   ! Declaration of variables. !
   !===========================!

   ! Zeta array.
   COMPLEX(rk), DIMENSION(:, :, :, :, :), ALLOCATABLE :: zeta

   ! Vector of unknowns and analytical Lx = Ux + Cx
   COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: xa, Ux, Cx, Lxc, Lxa, tmp

   ! Testcase numbers.
   INTEGER            :: testcase
   INTEGER, PARAMETER :: testcases = 34

   ! Wavenumber of free space.
   REAL(rk)    :: k02
   COMPLEX(rk) :: ik0

   ! Constant parameters.
   COMPLEX(rk), PARAMETER :: D0 = (0.64D0, -0.36D0)
   COMPLEX(rk), PARAMETER :: Phi0 = (-0.25D0, 0.75D0)
   COMPLEX(rk), PARAMETER :: H0 = (0.0D0, 1.0D0)
   COMPLEX(rk), PARAMETER :: A0 = (25.0D0, 100.0D0)

   !========================!
   ! Reading configuration. !
   !========================!
   CALL datain
   CALL piplasmain

   !====================!
   ! Allocating memory. !
   !====================!
   ALLOCATE (zeta(3, 3, nptx, npty, nptz))
   ALLOCATE (xa(10, nptx, npty, nptz))
   ALLOCATE (Ux(10, nptx, npty, nptz))
   ALLOCATE (Cx(10, nptx, npty, nptz))
   ALLOCATE (Lxc(10, nptx, npty, nptz))
   ALLOCATE (Lxa(10, nptx, npty, nptz))
   ALLOCATE (tmp(10, nptx, npty, nptz))

   !===============!
   ! Zeta setting. !
   !===============!
   CALL setzeta(zeta)

   !============================!
   ! Setting runtime constants. !
   !============================!
   ik0 = (0.0D0, 1.0D0)*pi/size0
   k02 = ABS(ik0)**2

   !===================================!
   ! Assigning analytical expressions. !
   !===================================!
   DO testcase = 1, testcases

      !=====================!
      ! Initial assignment. !
      !=====================!
      xa = (0.0D0, 0.0D0)
      Ux = (0.0D0, 0.0D0)
      Cx = (0.0D0, 0.0D0)
      Lxc = (0.0D0, 0.0D0)
      Lxa = (0.0D0, 0.0D0)
      tmp = (0.0D0, 0.0D0)

      SELECT CASE (testcase)

         !===================!
         ! Constant D tests. !
         !===================!

      CASE (1) ! Corresponds to Dx = D0
         xa(1, :, :, :) = D0
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                   ! Lxa(1)    = k0^2.D0
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))    ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)  ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))        ! Lxc(1)    = div(Zeta.D0)
         Lxa(4, :, :, :) = -ik0*Lxc(7, :, :, :)                    ! Lxa(4)    = - i.k0.div(Zeta.D0)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)               ! Lxa(8:10) = - k0^2.Zeta.D0

      CASE (2) ! Corresponds to Dy = D0
         xa(2, :, :, :) = D0
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                   ! Lxa(1)    = k0^2.D0
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))    ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)  ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))        ! Lxc(1)    = div(Zeta.D0)
         Lxa(4, :, :, :) = -ik0*Lxc(7, :, :, :)                    ! Lxa(4)    = - i.k0.div(Zeta.D0)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)               ! Lxa(8:10) = - k0^2.Zeta.D0

      CASE (3) ! Corresponds to Dz = D0
         xa(3, :, :, :) = D0
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                   ! Lxa(1)    = k0^2.D0
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))    ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)  ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))        ! Lxc(1)    = div(Zeta.D0)
         Lxa(4, :, :, :) = -ik0*Lxc(7, :, :, :)                    ! Lxa(4)    = - i.k0.div(Zeta.D0)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)               ! Lxa(8:10) = - k0^2.Zeta.D0

         !====================!
         ! Constant Phi test. !
         !====================!

      CASE (4) ! Correspondent to Phi = Phi0
         xa(4, :, :, :) = Phi0

         !===================!
         ! Constant H tests. !
         !===================!

      CASE (5) ! Correspondent to Hx = H0
         xa(5, :, :, :) = H0
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)

      CASE (6) ! Correspondent to Hy = H0
         xa(6, :, :, :) = H0
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)

      CASE (7) ! Correspondent to Hz = H0
         xa(7, :, :, :) = H0
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)

         !===================!
         ! Constant A tests. !
         !===================!

      CASE (8) ! Correspondent to Ax = A0
         xa(8, :, :, :) = A0
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))

      CASE (9) ! Correspondent to Ay = A0
         xa(9, :, :, :) = A0
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))

      CASE (10) ! Correspondent to Az = A0
         xa(10, :, :, :) = A0
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))

         !=================!
         ! Linear D tests. !
         !=================!

      CASE (11) ! Corresponds to Dx = D0.x
         CALL setlinxyz(tmp)
         xa(1, :, :, :) = D0*tmp(1, :, :, :)
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                           ! Lxa(1)    = k0^2.D0.x
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))            ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :))         ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)          ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))                ! Lxc(1)    = div(Zeta.D0)
         Lxa(4, :, :, :) = -ik0*(Lxc(7, :, :, :) + zeta(1, 1, :, :, :)*D0) ! Lxa(4)    = - i.k0.D0.x.grad(Zeta) - i.k0.Zeta.div(D0.x)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)                       ! Lxa(8:10) = - k0^2.Zeta.D0.x

      CASE (12) ! Corresponds to Dy = D0.y
         CALL setlinxyz(tmp)
         xa(2, :, :, :) = D0*tmp(2, :, :, :)
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                           ! Lxa(1)    = k0^2.D0.y
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))            ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :))         ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)          ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))                ! Lxc(1)    = div(Zeta.D0.y)
         Lxa(4, :, :, :) = -ik0*(Lxc(7, :, :, :) + zeta(2, 2, :, :, :)*D0) ! Lxa(4)    = - i.k0.D0.y.grad(Zeta) - i.k0.Zeta.div(D0.y)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)                       ! Lxa(8:10) = - k0^2.Zeta.D0.y

      CASE (13) ! Corresponds to Dz = D0.z
         CALL setlinxyz(tmp)
         xa(3, :, :, :) = D0*tmp(3, :, :, :)
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                           ! Lxa(1)    = k0^2.D0.z
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))            ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :))         ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)          ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))                ! Lxc(1)    = div(Zeta.D0.z)
         Lxa(4, :, :, :) = -ik0*(Lxc(7, :, :, :) + zeta(3, 3, :, :, :)*D0) ! Lxa(4)    = - i.k0.D0.z.grad(Zeta) - i.k0.Zeta.div(D0.z)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)                       ! Lxa(8:10) = - k0^2.Zeta.D0.z

         !===================!
         ! Linear Phi tests. !
         !===================!

      CASE (14) ! Corresponds to Phi = Phi0.x
         CALL setlinxyz(tmp)
         xa(4, :, :, :) = Phi0*tmp(1, :, :, :)
         Lxa(1, :, :, :) = ik0*Phi0                                          ! Lxa(1)    = i.k0.d(Phi0.x)/dx
         Lxa(8, :, :, :) = ik0*Phi0                                          ! Lxa(8)    = i.k0.d(Phi0.x)/dx

      CASE (15) ! Corresponds to Phi = Phi0.y
         CALL setlinxyz(tmp)
         xa(4, :, :, :) = Phi0*tmp(2, :, :, :)
         Lxa(2, :, :, :) = ik0*Phi0                                          ! Lxa(2)    = i.k0.d(Phi0.y)/dy
         Lxa(9, :, :, :) = ik0*Phi0                                          ! Lxa(9)    = i.k0.d(Phi0.y)/dy

      CASE (16) ! Corresponds to Phi = Phi0.z
         CALL setlinxyz(tmp)
         xa(4, :, :, :) = Phi0*tmp(3, :, :, :)
         Lxa(3, :, :, :) = ik0*Phi0                                          ! Lxa(3)    = i.k0.d(Phi0.z)/dz
         Lxa(10, :, :, :) = ik0*Phi0                                         ! Lxa(10)   = i.k0.d(Phi0.z)/dz

         !=================!
         ! Linear H tests. !
         !=================!

      CASE (17) ! Corresponds to Hx = H0.x
         CALL setlinxyz(tmp)
         xa(5, :, :, :) = H0*tmp(1, :, :, :)
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)                           ! Lxa(5:7)  = k0^2.H0.x
         Lxa(3, :, :, :) = -ik0*H0                                          ! Lxa(3)    = i.k0.(curl A)_z = i.k0.dHz/dy - i.k0.dHy/dz = 0, but curl(H0.x ex)_y = -dHz/dx (no), actually curl(Hx ex)_z = dHx/dy = 0, curl(Hx ex)_y = -dHx/dz = 0, curl(Hx ex)_x = 0; only D eqs pick up: (curl H)_z for Dz eq: d(H0.x)/dy = 0; for Dx eq (eq 6): -ik0(curl H)_x = 0; for Dy eq: -ik0(curl H)_y = -ik0(-dHx/dz) = 0
         Lxa(3, :, :, :) = (0.0D0, 0.0D0)                                    ! Lxa(3)    = 0 (curl of H0.x.ex has no z component from x variation only)
         Lxa(7, :, :, :) = -ik0*H0                                          ! Lxa(7)    = -i.k0.(curl A)_z: for H eq (eq 7) +ik0(curl A)_z; here A=0 so 0; but D correction: ik0(curl D)_z = 0
         Lxa(7, :, :, :) = (0.0D0, 0.0D0)
         ! For Hx = H0.x: curl(H)_y = dHx/dz - dHz/dx = 0, curl(H)_z = dHy/dx - dHx/dy = 0
         ! So only k0^2.H contribution survives in eq (7), and eqs (6,8,9) get no curl H terms.
         Lxa(5, :, :, :) = k02*xa(5, :, :, :)                                ! Lxa(5)    = k0^2.H0.x
         Lxa(6, :, :, :) = (0.0D0, 0.0D0)
         Lxa(7, :, :, :) = (0.0D0, 0.0D0)

      CASE (18) ! Corresponds to Hy = H0.y
         CALL setlinxyz(tmp)
         xa(6, :, :, :) = H0*tmp(2, :, :, :)
         ! curl(H0.y.ey): curl_x = dHz/dy - dHy/dz = 0, curl_y = dHx/dz - dHz/dx = 0, curl_z = dHy/dx - dHx/dy = 0
         ! Only k0^2.H contribution survives.
         Lxa(5, :, :, :) = (0.0D0, 0.0D0)
         Lxa(6, :, :, :) = k02*xa(6, :, :, :)                                ! Lxa(6)    = k0^2.H0.y
         Lxa(7, :, :, :) = (0.0D0, 0.0D0)

      CASE (19) ! Corresponds to Hz = H0.z
         CALL setlinxyz(tmp)
         xa(7, :, :, :) = H0*tmp(3, :, :, :)
         ! curl(H0.z.ez): curl_x = dHz/dy - dHy/dz = 0, curl_y = dHx/dz - dHz/dx = 0, curl_z = dHy/dx - dHx/dy = 0
         ! Only k0^2.H contribution survives.
         Lxa(5, :, :, :) = (0.0D0, 0.0D0)
         Lxa(6, :, :, :) = (0.0D0, 0.0D0)
         Lxa(7, :, :, :) = k02*xa(7, :, :, :)                                ! Lxa(7)    = k0^2.H0.z

         !=================!
         ! Linear A tests. !
         !=================!

      CASE (20) ! Corresponds to Ax = A0.x
         CALL setlinxyz(tmp)
         xa(8, :, :, :) = A0*tmp(1, :, :, :)
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)                          ! Lxa(8:10) = k0^2.A0.x
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))       ! Lxa(1:3)  = -dagger(Zeta).k0^2.A0.x
         Lxa(4, :, :, :) = ik0*A0                                             ! Lxa(4)    = i.k0.div(A0.x ex) = i.k0.A0

      CASE (21) ! Corresponds to Ay = A0.y
         CALL setlinxyz(tmp)
         xa(9, :, :, :) = A0*tmp(2, :, :, :)
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)                          ! Lxa(8:10) = k0^2.A0.y
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))       ! Lxa(1:3)  = -dagger(Zeta).k0^2.A0.y
         Lxa(4, :, :, :) = ik0*A0                                             ! Lxa(4)    = i.k0.div(A0.y ey) = i.k0.A0

      CASE (22) ! Corresponds to Az = A0.z
         CALL setlinxyz(tmp)
         xa(10, :, :, :) = A0*tmp(3, :, :, :)
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)                          ! Lxa(8:10) = k0^2.A0.z
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))       ! Lxa(1:3)  = -dagger(Zeta).k0^2.A0.z
         Lxa(4, :, :, :) = ik0*A0                                             ! Lxa(4)    = i.k0.div(A0.z ez) = i.k0.A0

         !===================!
         ! Quadratic D tests. !
         !===================!

      CASE (23) ! Corresponds to Dx = D0.x^2
         CALL setquadxyz(tmp)
         xa(1, :, :, :) = D0*tmp(1, :, :, :)
         ! -nabla^2(D0.x^2.ex) = -D0.d^2(x^2)/dx^2.ex = -2.D0.ex
         Lxa(1, :, :, :) = -2.0D0*D0 + k02*xa(1, :, :, :)                  ! Lxa(1)    = (-nabla^2 + k0^2).D0.x^2
         Lxa(2, :, :, :) = k02*xa(2, :, :, :)
         Lxa(3, :, :, :) = k02*xa(3, :, :, :)
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))             ! Lxc(1:3)  = Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :))          ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxa(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)           ! Lxa(1:3)  = (-nabla^2 + k0^2 + dagger(Zeta).Zeta.(-nabla^2 + k0^2)).Dx
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))
         ! div(D0.x^2.ex) = d(D0.x^2)/dx = 2.D0.x; so Phi eq: -ik0.(Zeta.grad(D) + D.div(Zeta)) = -ik0.(zeta_xx.2.D0.x.d/dx...)
         ! Using setdivzeta0 for the Zeta.D term, plus the plain divergence contribution through zeta_xx
         CALL setlinxyz(Lxc(1:3, :, :, :))
         Lxa(4, :, :, :) = -ik0*(Lxc(7, :, :, :) + 2.0D0*D0*zeta(1, 1, :, :, :)*Lxc(1, :, :, :)) ! Lxa(4) = -i.k0.(div(Zeta.D) + Zeta_xx.div(D))
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)                        ! Lxa(8:10) = - k0^2.Zeta.D0.x^2

      CASE (24) ! Corresponds to Dy = D0.y^2
         CALL setquadxyz(tmp)
         xa(2, :, :, :) = D0*tmp(2, :, :, :)
         ! -nabla^2(D0.y^2.ey) = -2.D0.ey
         Lxa(1, :, :, :) = k02*xa(1, :, :, :)
         Lxa(2, :, :, :) = -2.0D0*D0 + k02*xa(2, :, :, :)                  ! Lxa(2)    = (-nabla^2 + k0^2).D0.y^2
         Lxa(3, :, :, :) = k02*xa(3, :, :, :)
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :))
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))
         CALL setlinxyz(Lxc(1:3, :, :, :))
         Lxa(4, :, :, :) = -ik0*(Lxc(7, :, :, :) + 2.0D0*D0*zeta(2, 2, :, :, :)*Lxc(2, :, :, :)) ! Lxa(4) = -i.k0.(div(Zeta.D) + Zeta_yy.div(D))
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)                        ! Lxa(8:10) = - k0^2.Zeta.D0.y^2

      CASE (25) ! Corresponds to Dz = D0.z^2
         CALL setquadxyz(tmp)
         xa(3, :, :, :) = D0*tmp(3, :, :, :)
         ! -nabla^2(D0.z^2.ez) = -2.D0.ez
         Lxa(1, :, :, :) = k02*xa(1, :, :, :)
         Lxa(2, :, :, :) = k02*xa(2, :, :, :)
         Lxa(3, :, :, :) = -2.0D0*D0 + k02*xa(3, :, :, :)                  ! Lxa(3)    = (-nabla^2 + k0^2).D0.z^2
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :))
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)
         CALL setdivzeta0(Lxc(7, :, :, :), xa(1:3, :, :, :))
         CALL setlinxyz(Lxc(1:3, :, :, :))
         Lxa(4, :, :, :) = -ik0*(Lxc(7, :, :, :) + 2.0D0*D0*zeta(3, 3, :, :, :)*Lxc(3, :, :, :)) ! Lxa(4) = -i.k0.(div(Zeta.D) + Zeta_zz.div(D))
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)                        ! Lxa(8:10) = - k0^2.Zeta.D0.z^2

         !======================!
         ! Quadratic Phi tests. !
         !======================!

      CASE (26) ! Corresponds to Phi = Phi0.x^2
         CALL setquadxyz(tmp)
         xa(4, :, :, :) = Phi0*tmp(1, :, :, :)
         CALL setlinxyz(Lxc(1:3, :, :, :))
         ! grad(Phi0.x^2) = 2.Phi0.x.ex; contributes to A and D eqs via i.k0.grad(Phi)
         Lxa(1, :, :, :) = ik0*2.0D0*Phi0*Lxc(1, :, :, :)                   ! Lxa(1)    = i.k0.d(Phi0.x^2)/dx = i.k0.2.Phi0.x
         Lxa(8, :, :, :) = ik0*2.0D0*Phi0*Lxc(1, :, :, :)                   ! Lxa(8)    = i.k0.d(Phi0.x^2)/dx
         ! -nabla^2(Phi0.x^2) = -2.Phi0; contributes to Phi eq: i.k0.div(A) - i.k0.div(D) = 0 since A=D=0

      CASE (27) ! Corresponds to Phi = Phi0.y^2
         CALL setquadxyz(tmp)
         xa(4, :, :, :) = Phi0*tmp(2, :, :, :)
         CALL setlinxyz(Lxc(1:3, :, :, :))
         ! grad(Phi0.y^2) = 2.Phi0.y.ey
         Lxa(2, :, :, :) = ik0*2.0D0*Phi0*Lxc(2, :, :, :)                   ! Lxa(2)    = i.k0.d(Phi0.y^2)/dy = i.k0.2.Phi0.y
         Lxa(9, :, :, :) = ik0*2.0D0*Phi0*Lxc(2, :, :, :)                   ! Lxa(9)    = i.k0.d(Phi0.y^2)/dy

      CASE (28) ! Corresponds to Phi = Phi0.z^2
         CALL setquadxyz(tmp)
         xa(4, :, :, :) = Phi0*tmp(3, :, :, :)
         CALL setlinxyz(Lxc(1:3, :, :, :))
         ! grad(Phi0.z^2) = 2.Phi0.z.ez
         Lxa(3, :, :, :) = ik0*2.0D0*Phi0*Lxc(3, :, :, :)                   ! Lxa(3)    = i.k0.d(Phi0.z^2)/dz = i.k0.2.Phi0.z
         Lxa(10, :, :, :) = ik0*2.0D0*Phi0*Lxc(3, :, :, :)                  ! Lxa(10)   = i.k0.d(Phi0.z^2)/dz

         !===================!
         ! Quadratic H tests. !
         !===================!

      CASE (29) ! Corresponds to Hx = H0.x^2
         CALL setquadxyz(tmp)
         xa(5, :, :, :) = H0*tmp(1, :, :, :)
         ! -nabla^2(H0.x^2.ex) + k0^2.H0.x^2.ex = -2.H0.ex + k0^2.H0.x^2.ex
         Lxa(5, :, :, :) = -2.0D0*H0 + k02*xa(5, :, :, :)                  ! Lxa(5)    = (-nabla^2 + k0^2).H0.x^2
         Lxa(6, :, :, :) = k02*xa(6, :, :, :)
         Lxa(7, :, :, :) = k02*xa(7, :, :, :)

      CASE (30) ! Corresponds to Hy = H0.y^2
         CALL setquadxyz(tmp)
         xa(6, :, :, :) = H0*tmp(2, :, :, :)
         ! -nabla^2(H0.y^2.ey) + k0^2.H0.y^2.ey = -2.H0.ey + k0^2.H0.y^2.ey
         Lxa(5, :, :, :) = k02*xa(5, :, :, :)
         Lxa(6, :, :, :) = -2.0D0*H0 + k02*xa(6, :, :, :)                  ! Lxa(6)    = (-nabla^2 + k0^2).H0.y^2
         Lxa(7, :, :, :) = k02*xa(7, :, :, :)

      CASE (31) ! Corresponds to Hz = H0.z^2
         CALL setquadxyz(tmp)
         xa(7, :, :, :) = H0*tmp(3, :, :, :)
         ! -nabla^2(H0.z^2.ez) + k0^2.H0.z^2.ez = -2.H0.ez + k0^2.H0.z^2.ez
         Lxa(5, :, :, :) = k02*xa(5, :, :, :)
         Lxa(6, :, :, :) = k02*xa(6, :, :, :)
         Lxa(7, :, :, :) = -2.0D0*H0 + k02*xa(7, :, :, :)                  ! Lxa(7)    = (-nabla^2 + k0^2).H0.z^2

         !===================!
         ! Quadratic A tests. !
         !===================!

      CASE (32) ! Corresponds to Ax = A0.x^2
         CALL setquadxyz(tmp)
         xa(8, :, :, :) = A0*tmp(1, :, :, :)
         ! -nabla^2(A0.x^2.ex) + k0^2.A0.x^2.ex = -2.A0.ex + k0^2.A0.x^2.ex
         Lxa(8, :, :, :) = -2.0D0*A0 + k02*xa(8, :, :, :)                  ! Lxa(8)    = (-nabla^2 + k0^2).A0.x^2
         Lxa(9, :, :, :) = k02*xa(9, :, :, :)
         Lxa(10, :, :, :) = k02*xa(10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))       ! Lxa(1:3)  = -dagger(Zeta).(-nabla^2 + k0^2).A0.x^2
         CALL setlinxyz(Lxc(1:3, :, :, :))
         Lxa(4, :, :, :) = ik0*2.0D0*A0*Lxc(1, :, :, :)                     ! Lxa(4)    = i.k0.div(A0.x^2.ex) = i.k0.2.A0.x

      CASE (33) ! Corresponds to Ay = A0.y^2
         CALL setquadxyz(tmp)
         xa(9, :, :, :) = A0*tmp(2, :, :, :)
         ! -nabla^2(A0.y^2.ey) + k0^2.A0.y^2.ey = -2.A0.ey + k0^2.A0.y^2.ey
         Lxa(8, :, :, :) = k02*xa(8, :, :, :)
         Lxa(9, :, :, :) = -2.0D0*A0 + k02*xa(9, :, :, :)                  ! Lxa(9)    = (-nabla^2 + k0^2).A0.y^2
         Lxa(10, :, :, :) = k02*xa(10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))       ! Lxa(1:3)  = -dagger(Zeta).(-nabla^2 + k0^2).A0.y^2
         CALL setlinxyz(Lxc(1:3, :, :, :))
         Lxa(4, :, :, :) = ik0*2.0D0*A0*Lxc(2, :, :, :)                     ! Lxa(4)    = i.k0.div(A0.y^2.ey) = i.k0.2.A0.y

      CASE (34) ! Corresponds to Az = A0.z^2
         CALL setquadxyz(tmp)
         xa(10, :, :, :) = A0*tmp(3, :, :, :)
         ! -nabla^2(A0.z^2.ez) + k0^2.A0.z^2.ez = -2.A0.ez + k0^2.A0.z^2.ez
         Lxa(8, :, :, :) = k02*xa(8, :, :, :)
         Lxa(9, :, :, :) = k02*xa(9, :, :, :)
         Lxa(10, :, :, :) = -2.0D0*A0 + k02*xa(10, :, :, :)                ! Lxa(10)   = (-nabla^2 + k0^2).A0.z^2
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))       ! Lxa(1:3)  = -dagger(Zeta).(-nabla^2 + k0^2).A0.z^2
         CALL setlinxyz(Lxc(1:3, :, :, :))
         Lxa(4, :, :, :) = ik0*2.0D0*A0*Lxc(3, :, :, :)                     ! Lxa(4)    = i.k0.div(A0.z^2.ez) = i.k0.2.A0.z

      END SELECT

      !=======================================================!
      ! Calculating operator action using finite differences. !
      !=======================================================!
      Lxc = (0.0D0, 0.0D0)
      CALL mulmat(xa, Ux, tmp)
      CALL muladd(zeta, xa, Cx, tmp)
      Lxc = Ux + Cx

      !====================================================!
      ! Calculating and printing the normalized residuals. !
      !====================================================!
      PRINT *, "[INF] Test = ", testcase
      CALL respreset(Lxa, Lxc)

   END DO

   !======================!
   ! Deallocating memory. !
   !======================!
   DEALLOCATE (zeta, xa, Ux, Cx, Lxc, Lxa, tmp)

END PROGRAM test_l_operator
