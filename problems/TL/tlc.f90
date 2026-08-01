! Log of changes:
! (2026-06-26) Initial version
! (2026-07-10) Added debugs for zeta tensor multiplication and removed excessive k02 for D.
! (2026-07-11) Added a global residual without boundary calculation.
! (2026-07-11) Added a global residual without boundary calculation.

! TLC.1. Module declaration for constant plasma.
MODULE const_plasma

   USE numberformat

   IMPLICIT NONE
   COMPLEX(rk) :: E0 = 1.0_rk + (1.0_rk, 2.0_rk)/3.0d0

END MODULE const_plasma

SUBROUTINE setzeta(zeta)

   USE numberformat
   USE indata
   USE const_plasma

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(OUT) :: zeta

   zeta(1, 1, :, :, :) = (1.0D0, 0.0D0)/E0
   zeta(1, 2, :, :, :) = (0.0D0, 0.0D0)
   zeta(1, 3, :, :, :) = (0.0D0, 0.0D0)
   zeta(2, 1, :, :, :) = (0.0D0, 0.0D0)
   zeta(2, 2, :, :, :) = (1.0D0, 0.0D0)/E0
   zeta(2, 3, :, :, :) = (0.0D0, 0.0D0)
   zeta(3, 1, :, :, :) = (0.0D0, 0.0D0)
   zeta(3, 2, :, :, :) = (0.0D0, 0.0D0)
   zeta(3, 3, :, :, :) = (1.0D0, 0.0D0)/E0

END SUBROUTINE

! TL.4. A subroutine that calculates global, local interior and outer residual and prints it in a certain format.
SUBROUTINE respreset(vec1, vec2)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: vec1, vec2

   REAL(rk) :: num1, num2, den1, den2, candidate

   REAL(rk) :: respreset_with_boundary, respreset_without_boundary, maxrespreseti, maxrespreseto

   REAL(rk) :: wx, wy, wz ! Quadrature.

   INTEGER :: i, j, k
   INTEGER :: maxii = -1, maxij = -1, maxik = -1
   INTEGER :: maxoi = -1, maxoj = -1, maxok = -1

   num1 = 0.0D0
   den1 = 0.0D0
   num2 = 0.0D0
   den2 = 0.0D0
   candidate = 0.0D0

   respreset_with_boundary = 0.0D0
   respreset_without_boundary = 0.0D0
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

            num1 = num1 + SUM(ABS(vec1(:, i, j, k) - vec2(:, i, j, k))**2)*wx*wy*wz
            den1 = den1 + SUM(ABS(vec1(:, i, j, k))**2)*wx*wy*wz

            candidate = SUM(ABS(vec1(:, i, j, k) - vec2(:, i, j, k))**2)*wx*wy*wz
            IF (i > 1 .AND. i < nptx .AND. j > 1 .AND. j < npty .AND. k > 1 .AND. k < nptz) THEN

               num2 = num2 + SUM(ABS(vec1(:, i, j, k) - vec2(:, i, j, k))**2)
               den2 = den2 + SUM(ABS(vec1(:, i, j, k))**2)

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

   IF (den1 > 0.0D0 .AND. den2 > 0.0D0) THEN
      respreset_with_boundary = SQRT(num1/den1)
      respreset_without_boundary = SQRT(num2/den2)
      maxrespreseti = SQRT(nptx*npty*nptz*maxrespreseti/den1)
      maxrespreseto = SQRT(nptx*npty*nptz*maxrespreseto/den1)
   ELSE
      respreset_with_boundary = -1.0D0
      respreset_without_boundary = -1.0D0
      maxrespreseti = -1.0D0
      maxrespreseto = -1.0D0
   END IF

   PRINT *, "[DBG] Global residual = ", respreset_with_boundary
   PRINT *, "[DBG] Global residual (excluding boundary) = ", respreset_without_boundary
   PRINT *, "[DBG] Interior local residual = ", maxrespreseti, ", at node i, j, k = ", maxii, maxij, maxik
   PRINT *, "[DBG] Outer local residual = ", maxrespreseto, ", at node i, j, k = ", maxoi, maxoj, maxok

END SUBROUTINE respreset

PROGRAM test_lc_operator

   USE numberformat
   USE constants
   USE indata

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
   INTEGER, PARAMETER :: testcases = 10 ! 10

   ! Wavenumber of free space.
   REAL(rk)    :: k02
   COMPLEX(rk) :: ik0

   ! Constant parameters.
   COMPLEX(rk), PARAMETER :: D0 = (0.64D0, -0.36D0)
   COMPLEX(rk), PARAMETER :: Phi0 = (-0.25D0, 0.75D0)
   COMPLEX(rk), PARAMETER :: H0 = (0.4D0, 1.0D0)
   COMPLEX(rk), PARAMETER :: A0 = (25.0D0, 100.0D0)

   !========================!
   ! Reading configuration. !
   !========================!
   CALL datain

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

      PRINT *, "[INF] Test = ", testcase

      SELECT CASE (testcase)

         !===================!
         ! Constant D tests. !
         !===================!

      CASE (1) ! Corresponds to Dx = D0
         PRINT *, "Dx = ", D0
         xa(1, :, :, :) = D0
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                   ! Lxa(1)    = k0^2.D0
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))    ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)  ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         Lxa(4, :, :, :) = -ik0*Lxc(7, :, :, :)                    ! Lxa(4)    = - i.k0.div(Zeta.D0)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)               ! Lxa(8:10) = - k0^2.Zeta.D0

      CASE (2) ! Corresponds to Dy = D0
         PRINT *, "Dy = ", D0
         xa(2, :, :, :) = D0
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                   ! Lxa(1)    = k0^2.D0
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))    ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)  ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         Lxa(4, :, :, :) = -ik0*Lxc(7, :, :, :)                    ! Lxa(4)    = - i.k0.div(Zeta.D0)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)               ! Lxa(8:10) = - k0^2.Zeta.D0

      CASE (3) ! Corresponds to Dz = D0
         PRINT *, "Dz = ", D0
         xa(3, :, :, :) = D0
         Lxa(1:3, :, :, :) = k02*xa(1:3, :, :, :)                   ! Lxa(1)    = k0^2.D0
         CALL multen(Lxc(1:3, :, :, :), zeta, Lxa(1:3, :, :, :))    ! Lxc(1:3)  = k0^2.Zeta.Lxa(1:3)
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! Lxc(4:6)  = dagger(Zeta).Zeta.Lxc(1:3)
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + Lxc(4:6, :, :, :)  ! Lxa(1:3)  = k0^2.(1 + dagger(Zeta).Zeta).D0
         Lxa(4, :, :, :) = -ik0*Lxc(7, :, :, :)                    ! Lxa(4)    = - i.k0.div(Zeta.D0)
         Lxa(8:10, :, :, :) = -Lxc(1:3, :, :, :)               ! Lxa(8:10) = - k0^2.Zeta.D0

         !====================!
         ! Constant Phi test. !
         !====================!

      CASE (4) ! Correspondent to Phi = Phi0
         PRINT *, "Phi = ", Phi0
         xa(4, :, :, :) = Phi0

         !===================!
         ! Constant H tests. !
         !===================!

      CASE (5) ! Correspondent to Hx = H0
         PRINT *, "Hx = ", H0
         xa(5, :, :, :) = H0
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)

      CASE (6) ! Correspondent to Hy = H0
         PRINT *, "Hy = ", H0
         xa(6, :, :, :) = H0
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)

      CASE (7) ! Correspondent to Hz = H0
         PRINT *, "Hz = ", H0
         xa(7, :, :, :) = H0
         Lxa(5:7, :, :, :) = k02*xa(5:7, :, :, :)

         !===================!
         ! Constant A tests. !
         !===================!

      CASE (8) ! Correspondent to Ax = A0
         PRINT *, "Ax = ", A0
         xa(8, :, :, :) = A0
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))

      CASE (9) ! Correspondent to Ay = A0
         PRINT *, "Ay = ", A0
         xa(9, :, :, :) = A0
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))

      CASE (10) ! Correspondent to Az = A0
         PRINT *, "Az = ", A0
         xa(10, :, :, :) = A0
         Lxa(8:10, :, :, :) = k02*xa(8:10, :, :, :)
         CALL muldagten(Lxa(1:3, :, :, :), zeta, -Lxa(8:10, :, :, :))

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
      CALL respreset(Lxa, Lxc)

   END DO

   !======================!
   ! Deallocating memory. !
   !======================!
   DEALLOCATE (zeta, xa, Ux, Cx, Lxc, Lxa, tmp)

END PROGRAM test_lc_operator
