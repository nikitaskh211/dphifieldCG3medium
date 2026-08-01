! Log of changes:
! 1.

! TLC.1. Module declaration for constant plasma.
MODULE const_plasma

   USE numberformat

   IMPLICIT NONE
   COMPLEX(rk) :: E0 = 1.0_rk + (2.0_rk, -3.0_rk)

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

program test_l_constants

   use numberformat
   use constants
   use indata

   implicit none

end program test_l_constants
