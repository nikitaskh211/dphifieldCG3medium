! This module contains the definition of plasma.
MODULE const_plasma

   USE numberformat

   IMPLICIT NONE

   COMPLEX(rk) :: E0

END MODULE const_plasma

! TL.2. A subroutine that reads plasma's parameters from a file.
SUBROUTINE const_plasma_in

   USE numberformat
   USE const_plasma

   IMPLICIT NONE

   REAL(rk) :: re_E0, im_E0

   OPEN (14, FILE='plasma_in.dat', STATUS='old')
   OPEN (15, FILE='plasma_out.dat')

   READ (14, *) re_E0, im_E0
   WRITE (15, *) re_E0, im_E0

   E0 = 1.0D0 + re_E0 + (0.0D0, 1.0D0)*im_E0

   CLOSE (14)
   CLOSE (15)

END SUBROUTINE const_plasma_in

! This subroutine setts the zeta array.
SUBROUTINE set_zeta(zeta_arr)

   USE numberformat
   USE indata
   USE plasma

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(OUT) :: zeta_arr

   zeta_arr(1, 1, :, :, :) = 1.0D0/E0
   zeta_arr(1, 2, :, :, :) = 0.0D0
   zeta_arr(1, 3, :, :, :) = 0.0D0
   zeta_arr(2, 1, :, :, :) = 0.0D0
   zeta_arr(2, 2, :, :, :) = 1.0D0/E0
   zeta_arr(2, 3, :, :, :) = 0.0D0
   zeta_arr(3, 1, :, :, :) = 0.0D0
   zeta_arr(3, 2, :, :, :) = 0.0D0
   zeta_arr(3, 3, :, :, :) = 1.0D0/E0

END SUBROUTINE set_zeta

! This subroutine computes residuals of Mx_analytical and Mx_computed.
SUBROUTINE res_preset(vec1, vec2)

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

END SUBROUTINE res_preset

PROGRAM test_ll1_operator

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   ! Zeta array.
   COMPLEX(rk), DIMENSION(:, :, :, :, :), ALLOCATABLE :: zeta

   ! Vector of unknowns and analytical Lx = Ux + Cx
   ! Layout (10 components): 1:3 = D(x,y,z), 4 = Phi, 5:7 = H(x,y,z), 8:10 = A(x,y,z)
   COMPLEX(rk), DIMENSION(:, :, :, :), ALLOCATABLE :: xa, Ux, Cx, Lxc, Lxa, tmp

   ! Iterator numbers.
   INTEGER :: iterator

   ! Testcase numbers.
   INTEGER            :: testcase
   INTEGER, PARAMETER :: testcases = 60 ! 10 components of the vector: Dx,Dy,Dz,Phi,Hx,Hy,Hz,Ax,Ay,Az
                                        ! x 6 spatial profiles (x, y, z, i*x, i*y, i*z) = 60 unique test cases

   ! Wavenumber of free space.
   REAL(rk)    :: k02
   COMPLEX(rk) :: ik0

   ! Reading configuration.
   CALL datain
   CALL const_plasma_in

   ! Allocating memory.
   ALLOCATE(zeta(3, 3, nptx, npty, nptz))
   ALLOCATE(xa(10, nptx, npty, nptz))
   ALLOCATE(Ux(10, nptx, npty, nptz))
   ALLOCATE(Cx(10, nptx, npty, nptz))
   ALLOCATE(Lxc(10, nptx, npty, nptz))
   ALLOCATE(Lxa(10, nptx, npty, nptz))
   ALLOCATE(tmp(3, nptx, npty, nptz))

   ! Setting constants.
   ik0 = (0.0D0, 1.0D0)*pi/size0
   k02 = ABS(ik0)**2

   ! Setting zeta.
   CALL set_zeta(zeta)

   DO testcase = 1, testcases

      ! Initial assignment.
      xa = (0.0D0, 0.0D0)
      Ux = (0.0D0, 0.0D0)
      Cx = (0.0D0, 0.0D0)
      Lxc = (0.0D0, 0.0D0)
      Lxa = (0.0D0, 0.0D0)
      tmp = (0.0D0, 0.0D0)

      PRINT *, "[INF] Test = ", testcase

      SELECT CASE(testcase)

      CASE (1) ! Dx = x

         PRINT *, "Dx = x"
         DO iterator = 1, nptx

            xa(1, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = -ik0 * (1.0D0, 0.0D0) * zeta(1, 1, :, :, :)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (2) ! Dx = y

         PRINT *, "Dx = y"
         DO iterator = 1, npty

            xa(1, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = -ik0 ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (3) ! Dx = z

         PRINT *, "Dx = z"
         DO iterator = 1, nptz

            xa(1, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = ik0 ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (4) ! Dx = i*x

         PRINT *, "Dx = i*x"
         DO iterator = 1, nptx

            xa(1, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = -ik0 * (0.0D0, 1.0D0) * zeta(1, 1, :, :, :)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (5) ! Dx = i*y

         PRINT *, "Dx = i*y"
         DO iterator = 1, npty

            xa(1, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, -1.0D0)*ik0 ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (6) ! Dx = i*z

         PRINT *, "Dx = i*z"
         DO iterator = 1, nptz

            xa(1, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (7) ! Dy = x

         PRINT *, "Dy = x"
         DO iterator = 1, nptx

            xa(2, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = ik0 ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (8) ! Dy = y

         PRINT *, "Dy = y"
         DO iterator = 1, npty

            xa(2, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = -ik0 * (1.0D0, 0.0D0) * zeta(2, 2, :, :, :)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (9) ! Dy = z

         PRINT *, "Dy = z"
         DO iterator = 1, nptz

            xa(2, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = -ik0 ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (10) ! Dy = i*x

         PRINT *, "Dy = i*x"
         DO iterator = 1, nptx

            xa(2, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (11) ! Dy = i*y

         PRINT *, "Dy = i*y"
         DO iterator = 1, npty

            xa(2, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = -ik0 * (0.0D0, 1.0D0) * zeta(2, 2, :, :, :)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (12) ! Dy = i*z

         PRINT *, "Dy = i*z"
         DO iterator = 1, nptz

            xa(2, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, -1.0D0)*ik0 ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (13) ! Dz = x

         PRINT *, "Dz = x"
         DO iterator = 1, nptx

            xa(3, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = -ik0 ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (14) ! Dz = y

         PRINT *, "Dz = y"
         DO iterator = 1, npty

            xa(3, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = ik0 ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (15) ! Dz = z

         PRINT *, "Dz = z"
         DO iterator = 1, nptz

            xa(3, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = -ik0 * (1.0D0, 0.0D0) * zeta(3, 3, :, :, :)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (16) ! Dz = i*x

         PRINT *, "Dz = i*x"
         DO iterator = 1, nptx

            xa(3, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, -1.0D0)*ik0 ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (17) ! Dz = i*y

         PRINT *, "Dz = i*y"
         DO iterator = 1, npty

            xa(3, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (18) ! Dz = i*z

         PRINT *, "Dz = i*z"
         DO iterator = 1, nptz

            xa(3, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- D-eq (18): k0^2 D + k0^2 zeta^T.zeta.D  (grad-div term = 0 for linear input) ---
         Lxa(1:3, :, :, :) = k02 * xa(1:3, :, :, :) ! k0^2 * D
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D
         CALL muldagten(Lxc(4:6, :, :, :), zeta, Lxc(1:3, :, :, :)) ! zeta^T . zeta . D
         Lxa(1:3, :, :, :) = Lxa(1:3, :, :, :) + k02 * Lxc(4:6, :, :, :) ! + k0^2 * zeta^T.zeta.D

         ! --- Phi-eq (19): -i k0 * div(zeta . D) ---
         Lxa(4, :, :, :) = -ik0 * (0.0D0, 1.0D0) * zeta(3, 3, :, :, :)

         ! --- H-eq (17): i k0 * curl(D)  (constant, since D is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(D)[3]

         ! --- A-eq (16): -k0^2 * zeta . D ---
         CALL multen(Lxc(1:3, :, :, :), zeta, xa(1:3, :, :, :)) ! zeta . D (recompute, xa unchanged)
         Lxa(8:10, :, :, :) = -k02 * Lxc(1:3, :, :, :)

      CASE (19) ! Phi = x

         PRINT *, "Phi = x"
         DO iterator = 1, nptx

            xa(4, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): i k0 * grad(Phi)  (constant, since Phi is linear) ---
         Lxa(8, :, :, :) = ik0 ! i*k0*grad(Phi)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[3]

         ! --- H-eq (17): Phi does not enter ---
         Lxa(5:7, :, :, :) = (0.0D0, 0.0D0)

         ! --- D-eq (18): -i k0 * zeta^T . grad(Phi) ---
         Lxa(1:3, :, :, :) = (0.0D0, 0.0D0)
         Lxa(1, :, :, :) = -ik0 * (1.0D0, 0.0D0) * zeta(1, 1, :, :, :)

         ! --- Phi-eq (19): self-Laplacian = 0 for linear input, no other Phi terms ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (20) ! Phi = y

         PRINT *, "Phi = y"
         DO iterator = 1, npty

            xa(4, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): i k0 * grad(Phi)  (constant, since Phi is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[1]
         Lxa(9, :, :, :) = ik0 ! i*k0*grad(Phi)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[3]

         ! --- H-eq (17): Phi does not enter ---
         Lxa(5:7, :, :, :) = (0.0D0, 0.0D0)

         ! --- D-eq (18): -i k0 * zeta^T . grad(Phi) ---
         Lxa(1:3, :, :, :) = (0.0D0, 0.0D0)
         Lxa(2, :, :, :) = -ik0 * (1.0D0, 0.0D0) * zeta(2, 2, :, :, :)

         ! --- Phi-eq (19): self-Laplacian = 0 for linear input, no other Phi terms ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (21) ! Phi = z

         PRINT *, "Phi = z"
         DO iterator = 1, nptz

            xa(4, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): i k0 * grad(Phi)  (constant, since Phi is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[2]
         Lxa(10, :, :, :) = ik0 ! i*k0*grad(Phi)[3]

         ! --- H-eq (17): Phi does not enter ---
         Lxa(5:7, :, :, :) = (0.0D0, 0.0D0)

         ! --- D-eq (18): -i k0 * zeta^T . grad(Phi) ---
         Lxa(1:3, :, :, :) = (0.0D0, 0.0D0)
         Lxa(3, :, :, :) = -ik0 * (1.0D0, 0.0D0) * zeta(3, 3, :, :, :)

         ! --- Phi-eq (19): self-Laplacian = 0 for linear input, no other Phi terms ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (22) ! Phi = i*x

         PRINT *, "Phi = i*x"
         DO iterator = 1, nptx

            xa(4, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): i k0 * grad(Phi)  (constant, since Phi is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*grad(Phi)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[3]

         ! --- H-eq (17): Phi does not enter ---
         Lxa(5:7, :, :, :) = (0.0D0, 0.0D0)

         ! --- D-eq (18): -i k0 * zeta^T . grad(Phi) ---
         Lxa(1:3, :, :, :) = (0.0D0, 0.0D0)
         Lxa(1, :, :, :) = -ik0 * (0.0D0, 1.0D0) * zeta(1, 1, :, :, :)

         ! --- Phi-eq (19): self-Laplacian = 0 for linear input, no other Phi terms ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (23) ! Phi = i*y

         PRINT *, "Phi = i*y"
         DO iterator = 1, npty

            xa(4, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): i k0 * grad(Phi)  (constant, since Phi is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[1]
         Lxa(9, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*grad(Phi)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[3]

         ! --- H-eq (17): Phi does not enter ---
         Lxa(5:7, :, :, :) = (0.0D0, 0.0D0)

         ! --- D-eq (18): -i k0 * zeta^T . grad(Phi) ---
         Lxa(1:3, :, :, :) = (0.0D0, 0.0D0)
         Lxa(2, :, :, :) = -ik0 * (0.0D0, 1.0D0) * zeta(2, 2, :, :, :)

         ! --- Phi-eq (19): self-Laplacian = 0 for linear input, no other Phi terms ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (24) ! Phi = i*z

         PRINT *, "Phi = i*z"
         DO iterator = 1, nptz

            xa(4, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): i k0 * grad(Phi)  (constant, since Phi is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! i*k0*grad(Phi)[2]
         Lxa(10, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*grad(Phi)[3]

         ! --- H-eq (17): Phi does not enter ---
         Lxa(5:7, :, :, :) = (0.0D0, 0.0D0)

         ! --- D-eq (18): -i k0 * zeta^T . grad(Phi) ---
         Lxa(1:3, :, :, :) = (0.0D0, 0.0D0)
         Lxa(3, :, :, :) = -ik0 * (0.0D0, 1.0D0) * zeta(3, 3, :, :, :)

         ! --- Phi-eq (19): self-Laplacian = 0 for linear input, no other Phi terms ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (25) ! Hx = x

         PRINT *, "Hx = x"
         DO iterator = 1, nptx

            xa(5, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (26) ! Hx = y

         PRINT *, "Hx = y"
         DO iterator = 1, npty

            xa(5, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = ik0 ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = ik0 ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (27) ! Hx = z

         PRINT *, "Hx = z"
         DO iterator = 1, nptz

            xa(5, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = -ik0 ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = -ik0 ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (28) ! Hx = i*x

         PRINT *, "Hx = i*x"
         DO iterator = 1, nptx

            xa(5, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (29) ! Hx = i*y

         PRINT *, "Hx = i*y"
         DO iterator = 1, npty

            xa(5, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 1.0D0)*ik0 ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 1.0D0)*ik0 ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (30) ! Hx = i*z

         PRINT *, "Hx = i*z"
         DO iterator = 1, nptz

            xa(5, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, -1.0D0)*ik0 ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, -1.0D0)*ik0 ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (31) ! Hy = x

         PRINT *, "Hy = x"
         DO iterator = 1, nptx

            xa(6, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = -ik0 ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = -ik0 ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (32) ! Hy = y

         PRINT *, "Hy = y"
         DO iterator = 1, npty

            xa(6, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (33) ! Hy = z

         PRINT *, "Hy = z"
         DO iterator = 1, nptz

            xa(6, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = ik0 ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = ik0 ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (34) ! Hy = i*x

         PRINT *, "Hy = i*x"
         DO iterator = 1, nptx

            xa(6, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, -1.0D0)*ik0 ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, -1.0D0)*ik0 ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (35) ! Hy = i*y

         PRINT *, "Hy = i*y"
         DO iterator = 1, npty

            xa(6, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (36) ! Hy = i*z

         PRINT *, "Hy = i*z"
         DO iterator = 1, nptz

            xa(6, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 1.0D0)*ik0 ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 1.0D0)*ik0 ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (37) ! Hz = x

         PRINT *, "Hz = x"
         DO iterator = 1, nptx

            xa(7, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = ik0 ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = ik0 ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (38) ! Hz = y

         PRINT *, "Hz = y"
         DO iterator = 1, npty

            xa(7, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = -ik0 ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = -ik0 ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (39) ! Hz = z

         PRINT *, "Hz = z"
         DO iterator = 1, nptz

            xa(7, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (40) ! Hz = i*x

         PRINT *, "Hz = i*x"
         DO iterator = 1, nptx

            xa(7, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 1.0D0)*ik0 ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 1.0D0)*ik0 ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (41) ! Hz = i*y

         PRINT *, "Hz = i*y"
         DO iterator = 1, npty

            xa(7, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, -1.0D0)*ik0 ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, -1.0D0)*ik0 ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (42) ! Hz = i*z

         PRINT *, "Hz = i*z"
         DO iterator = 1, nptz

            xa(7, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): -i k0 * curl(H)  (constant, since H is linear) ---
         Lxa(8, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(9, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(10, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- H-eq (17): k0^2 * H (self term) ---
         Lxa(5:7, :, :, :) = k02 * xa(5:7, :, :, :)

         ! --- D-eq (18): -i k0 * curl(H)  (same constant as in A-eq) ---
         Lxa(1, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[1]
         Lxa(2, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[2]
         Lxa(3, :, :, :) = (0.0D0, 0.0D0) ! -i*k0*curl(H)[3]

         ! --- Phi-eq (19): H does not enter ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (43) ! Ax = x

         PRINT *, "Ax = x"
         DO iterator = 1, nptx

            xa(8, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = ik0

      CASE (44) ! Ax = y

         PRINT *, "Ax = y"
         DO iterator = 1, npty

            xa(8, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = -ik0 ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (45) ! Ax = z

         PRINT *, "Ax = z"
         DO iterator = 1, nptz

            xa(8, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = ik0 ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (46) ! Ax = i*x

         PRINT *, "Ax = i*x"
         DO iterator = 1, nptx

            xa(8, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 1.0D0)*ik0

      CASE (47) ! Ax = i*y

         PRINT *, "Ax = i*y"
         DO iterator = 1, npty

            xa(8, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, -1.0D0)*ik0 ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (48) ! Ax = i*z

         PRINT *, "Ax = i*z"
         DO iterator = 1, nptz

            xa(8, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (49) ! Ay = x

         PRINT *, "Ay = x"
         DO iterator = 1, nptx

            xa(9, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = ik0 ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (50) ! Ay = y

         PRINT *, "Ay = y"
         DO iterator = 1, npty

            xa(9, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = ik0

      CASE (51) ! Ay = z

         PRINT *, "Ay = z"
         DO iterator = 1, nptz

            xa(9, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = -ik0 ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (52) ! Ay = i*x

         PRINT *, "Ay = i*x"
         DO iterator = 1, nptx

            xa(9, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (53) ! Ay = i*y

         PRINT *, "Ay = i*y"
         DO iterator = 1, npty

            xa(9, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 1.0D0)*ik0

      CASE (54) ! Ay = i*z

         PRINT *, "Ay = i*z"
         DO iterator = 1, nptz

            xa(9, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, -1.0D0)*ik0 ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (55) ! Az = x

         PRINT *, "Az = x"
         DO iterator = 1, nptx

            xa(10, iterator, :, :) = (iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = -ik0 ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (56) ! Az = y

         PRINT *, "Az = y"
         DO iterator = 1, npty

            xa(10, :, iterator, :) = (iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = ik0 ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (57) ! Az = z

         PRINT *, "Az = z"
         DO iterator = 1, nptz

            xa(10, :, :, iterator) = (iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = ik0

      CASE (58) ! Az = i*x

         PRINT *, "Az = i*x"
         DO iterator = 1, nptx

            xa(10, iterator, :, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizex/(nptx - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, -1.0D0)*ik0 ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (59) ! Az = i*y

         PRINT *, "Az = i*y"
         DO iterator = 1, npty

            xa(10, :, iterator, :) = (0.0D0, 1.0D0)*(iterator - 1)*sizey/(npty - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 1.0D0)*ik0 ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 0.0D0)

      CASE (60) ! Az = i*z

         PRINT *, "Az = i*z"
         DO iterator = 1, nptz

            xa(10, :, :, iterator) = (0.0D0, 1.0D0)*(iterator - 1)*sizez/(nptz - 1)

         END DO

         ! --- A-eq (16): k0^2 * A (self term) ---
         Lxa(8:10, :, :, :) = k02 * xa(8:10, :, :, :)

         ! --- H-eq (17): i k0 * curl(A)  (constant, since A is linear) ---
         Lxa(5, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[1]
         Lxa(6, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[2]
         Lxa(7, :, :, :) = (0.0D0, 0.0D0) ! i*k0*curl(A)[3]

         ! --- D-eq (18): -k0^2 * zeta^T . A ---
         CALL muldagten(Lxc(1:3, :, :, :), zeta, xa(8:10, :, :, :)) ! zeta^T . A
         Lxa(1:3, :, :, :) = -k02 * Lxc(1:3, :, :, :)

         ! --- Phi-eq (19): i k0 * div(A) ---
         Lxa(4, :, :, :) = (0.0D0, 1.0D0)*ik0

      END SELECT

      !=======================================================!
      ! Calculating operator action using finite differences. !
      !=======================================================!
      Lxc = (0.0D0, 0.0D0)
     !CALL mulmat(xa, Ux, tmp)
     !CALL muladd(zeta, xa, Cx, tmp)
      CALL mul_l_vacuum(Ux, tmp, xa)
      CALL mul_l_correction(Cx, tmp, xa, zeta)
      Lxc = Ux + Cx

      !====================================================!
      ! Calculating and printing the normalized residuals. !
      !====================================================!
      CALL res_preset(Lxa, Lxc)

   END DO

   !======================!
   ! Deallocating memory. !
   !======================!
   DEALLOCATE(zeta, xa, Ux, Cx, Lxc, Lxa, tmp)

END PROGRAM