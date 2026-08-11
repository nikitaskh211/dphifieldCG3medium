! Log of changes:
! (APR-21-2026) - Added a log of changes (changes started before this log was introduced hence why all the changes are "made" starting today).
! (APR-21-2026) - Added indentation for easier reading.
! (APR-21-2026) - Re-did all comments of the main code branch.
! (APR-21-2026) - Changed all constants from being set as "_rk" to an explicit [M.M]D[E] format to avoid conversion errors.
! (APR-21-2026) - Commented all the unused parameters in "constants" module (lines 30-39).
! (APR-21-2026) - Changed the speed of light value (line 28) for greater physical accuracy.
! (APR-21-2026) - Commented an unused "eps0" parameter in "indata" module (line 53).
! (APR-21-2026) - Moved "datain" subroutine right under "indata" module. Commented unused "chdum" and "j" parameters (lines 66-67) and obsolete instructions (lines 78-82).
! (APR-21-2026) - Changed unit number from "20" to "11" in "datain" subroutine.
! (APR-21-2026) - Changed all "CMPLX(0.0D0,k0)" to "(0.0D0, 1.0D0)*k0" to eliminate all GNU Fortran compiler warnings related to this.
! (APR-22-2026) - Commented unused variables in "graddiv" subroutine (lines 798-807).
! (APR-22-2026) - Commented an unused "l" variable in "laplas" subroutine (line 957).
! (APR-22-2026) - Identified all possible bugs found in "laplash" subroutine (lines 1168-1322). Explanation: there is an explicit "- k0^2.H" for each edge. If it is a bug variables and instructions for its instructions will get commented.
! (APR-22-2026) - Changed "central finite difference", "forward finite difference" and "backward finite difference" comments to "CFD", "FFD", "BFD" respectfully to ease reading and reduce clutter.
! (APR-22-2026) - Commented an unused "k02" variable in laplasa subroutine (line 1344) and obsolete instruction (line 1364).
! (APR-22-2026) - Changed all assignments of complex numbers from "0.0D0" to a more explicit "(0.0D0, 0.0D0)" to eliminate all GNU Fortran compiler warnings related to this.
! (APR-22-2026) - Commented an unused "l" variable in "nablaphi", "roth", "rotd", and "divd" subroutines (lines 1488, 1623, 1745, and 2006 respectfully).
! (APR-22-2026) - Added an explicit "INTENT(INOUT)" for "tempr" variable in "mulmat" subroutine (line 2057).
! (APR-22-2026) - Added two "ik0" and "k02" variables in "mulmat" subroutine (lines 2060-2061).
! (APR-22-2026) - Identified all possible bugs found in "mulmat" subroutine (lines 2145-2161). Explanation: normally just mass terms like D, H, or A must be halved in order to align with the CFD, however these half entire equations at edges which is not logical. Makes sense if all the equations get the same threatment.
! (APR-22-2026) - Identified and fixed a bug in "scprod" subroutine (line 2205). Explanation: Unintentional conversion from COMPLEX(rk) to REAL(rk) datatype.
! (APR-22-2026) - Identified an inefficiency in "scprod" subroutine. Explanation: Slow memory access (see detailed explanation at lines 2197-2199.)
! (APR-22-2026) - Only now i realized that with each comment i make in the header the line-numbers shift downwards. Make use of the description of changes.
! (APR-22-2026) - Added an explicit "INTENT(IN)" for "left" and "right" variables in "scprod" subroutine (line 2208).
! (APR-22-2026) - Identified and fixed a bug in "solerr" (line 2287). Explanation: Misplaced square.
! (APR-23-2026) - This file now contains exclusively essential definitions of modules and subroutines for problem-solving. All problems now must be defined in separate files. To compile use CLI: "gfortran fd_pde_pde_0bc.f90 [problem].f90 -o [executable].exe"
! (APR-23-2026) - Legacy "testset", "testset1" subroutines and "efield" program are now in a separate file "dphifieldCG3.f90".
! (APR-23-2026) - Introduced new code: "multen", "muldagten", "muladd" and "setrhs" subroutines with "resdivd", "resdivh", "resdiva", and "rescurlsrc" functions.
! (APR-24-2026) - To compile each problem for calculation use dedicated "compile.bat" file.
! (APR-26-2026) - Commented all supposed bugs in "laplash" and "mulmat" subroutines.
! (APR-26-2026) - Fixed boundary halving for A-like mass term in "muladd" subroutine.
! (APR-30-2026) - Added calculations for a pointwise maximum residual to secondary residual subroutines: "resdivd", "resdivh", "resdiva", "rescurlsrc".
! (MAY-08-2026) - Introduced correct boundary corrections to H governing equation in mulmat.

! Legacy code:
! Note: The following modules and subroutines are responsible for computing positive operator action of U for a single iteration.

! 1. Module that defines 64 bit (8 byte) floating point numbers.
MODULE numberformat

   IMPLICIT NONE

   INTEGER, PARAMETER :: rk = kind(1.0D0)

END MODULE numberformat

! 2. Module that defines physical constants in Gaussian system of units.
MODULE constants

   USE numberformat

   IMPLICIT NONE

   REAL(rk), PARAMETER :: lightspeed = 29979245800D0     ! Value of the speed of light,
   REAL(rk), PARAMETER :: pi = 3.141592654D0     ! Value of pi,
   !REAL(rk), PARAMETER :: mu0        = 4.0D0*pi*1D-9    ! Value of magnetic medium permeability of vacuum,
   !REAL(rk), PARAMETER :: epsylon0   = 8.85D-14          ! Value of electric medium permittivity of vacuum,
   !REAL(rk), PARAMETER :: ec         = 4.8088D-10        ! Value of the electron's charge,
   !REAL(rk), PARAMETER :: em         = 9.1094D-28        ! Value of the electron's mass,
   !REAL(rk), PARAMETER :: pm         = 1.6726D-24        ! Value of the positron's mas,
   !REAL(rk), PARAMETER :: pm         = 2.0D0*1.6726D-24 ! Doubled value of the positron's mass,
   !REAL(rk), PARAMETER :: CL         = 13                ! N/A,
   !REAL(rk), PARAMETER :: CONBOL     = 1.602D-12         ! N/A,
   !REAL(rk), PARAMETER :: SVEN       = 1D-7              ! N/A,
   !REAL(rk), PARAMETER :: SVIN       = 1D-8              ! N/A.

END MODULE constants

! 3.1. Module that defines collocated Cartesian 3D grid parameters, CG tolerance parameter, L0 parameter for k0.
MODULE indata

   USE numberformat

   IMPLICIT none

   INTEGER  :: nptx, npty, nptz    ! Amounts of nodes in each coordinate direction,

   REAL(rk) :: sizex, sizey, sizez ! Length dimensions of the cavity [0, Lx]x[0, Ly]x[0, Lz],
   !REAL(rk) :: eps0                ! N/A,
   REAL(rk) :: accur, size0        ! CG tolerance and L0.

END MODULE indata

! 3.2. A subroutine for reading values from a file "in3.dat" and writing values into a file "out3.dat" both for CG and debugging.
SUBROUTINE datain

   USE numberformat
   USE indata

   IMPLICIT NONE

   !CHARACTER(80)::chdum
   !INTEGER:: j

   OPEN (10, FILE='in3.dat', STATUS='old')
   OPEN (11, FILE='out3.dat')

   READ (10, *) nptx, npty, nptz
   WRITE (11, *) nptx, npty, nptz

   READ (10, *) sizex, sizey, sizez
   WRITE (11, *) sizex, sizey, sizez

   !READ(10, *) accur, eps0
   !WRITE(11, *) accur, eps0

   !READ(10, *) size0
   !WRITE(11, *) size0

   READ (10, *) accur, size0
   WRITE (11, *) accur, size0

   CLOSE (10)
   CLOSE (11)

END SUBROUTINE datain

! 3.3. A subroutine for calculating a gradient of a divergence from "dvec" array into "res" array.
! Note: Applies boundary conditions for D: Normal components have Neumann boundary conditions, tangent components have Dirichlet boundary conditions.
SUBROUTINE graddiv(res, dvec)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT none

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: dvec ! Input array.

   !REAL(rk) :: argx, argy, argz ! Arguments of the trigonometric functions,
   !REAL(rk) :: sinx, siny, sinz ! Values of SIN(Arguments),
   !REAL(rk) :: cosx, cosy, cosz ! Values of COS(Arguments).

   REAL(rk) :: hx, hy, hz ! Node spacings,
   REAL(rk) :: x2, y2, z2 ! Squared recipocals of node spacings,
   REAL(rk) :: xy, xz, yz ! Mixed recipocals of node spacings.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   !INTEGER :: l       ! N/A.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of squared recipocals of node spacings.
   x2 = 1.0D0/hx**2
   y2 = 1.0D0/hy**2
   z2 = 1.0D0/hz**2

   ! Calculation of mixed recipocals of node spacings.
   xy = 0.25D0/(hx*hy)
   xz = 0.25D0/(hx*hz)
   yz = 0.25D0/(hy*hz)

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1

            ! Discrete d2Dx/dx2 + d2Dy/dydx + d2Dz/dzdx all using CFDs.
            res(1, i, j, k) = (dvec(1, i + 1, j, k) - 2.0D0*dvec(1, i, j, k) + dvec(1, i - 1, j, k))*x2 &
                  + (dvec(2, i + 1, j + 1, k) + dvec(2, i - 1, j - 1, k) - dvec(2, i + 1, j - 1, k) - dvec(2, i - 1, j + 1, k))*xy &
                    + (dvec(3, i + 1, j, k + 1) + dvec(3, i - 1, j, k - 1) - dvec(3, i + 1, j, k - 1) - dvec(3, i - 1, j, k + 1))*xz

            ! Discrete d2Dx/dxdy + d2Dy/dy2 + d2Dz/dzdy all using CFDs.
  res(2, i, j, k) = (dvec(1, i + 1, j + 1, k) + dvec(1, i - 1, j - 1, k) - dvec(1, i + 1, j - 1, k) - dvec(1, i - 1, j + 1, k))*xy &
                              + (dvec(2, i, j + 1, k) - 2.0D0*dvec(2, i, j, k) + dvec(2, i, j - 1, k))*y2 &
                    + (dvec(3, i, j + 1, k + 1) + dvec(3, i, j - 1, k - 1) - dvec(3, i, j + 1, k - 1) - dvec(3, i, j - 1, k + 1))*yz

            ! Discrete d2Dx/dxdz + d2Dy/dydz + d2Dz/dz2 all using CFDs.
  res(3, i, j, k) = (dvec(1, i + 1, j, k + 1) + dvec(1, i - 1, j, k - 1) - dvec(1, i + 1, j, k - 1) - dvec(1, i - 1, j, k + 1))*xz &
                  + (dvec(2, i, j + 1, k + 1) + dvec(2, i, j - 1, k - 1) - dvec(2, i, j + 1, k - 1) - dvec(2, i, j - 1, k + 1))*yz &
                              + (dvec(3, i, j, k + 1) - 2.0D0*dvec(3, i, j, k) + dvec(3, i, j, k - 1))*z2

         END DO
      END DO
   END DO

   ! Boundary nodes for x = 0.
   i = 1
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete d2Dx/dx2 + d2Dy/dydx + d2Dz/dzdx using FFD for x and CFDs for y and z coordinate.
         res(1, i, j, k) = (dvec(1, i + 1, j, k) - dvec(1, i, j, k))*x2 &
                           + (dvec(2, i + 1, j + 1, k) - dvec(2, i + 1, j - 1, k))*xy &
                           + (dvec(3, i + 1, j, k + 1) - dvec(3, i + 1, j, k - 1))*xz

      END DO
   END DO

   ! Boundary nodes for x = Lx.
   i = nptx
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete d2Dx/dx2 + d2Dy/dydx + d2Dz/dzdx using BFD for x and CFDs for y and z coordinate.
         res(1, i, j, k) = (dvec(1, i - 1, j, k) - dvec(1, i, j, k))*x2 &
                           + (dvec(2, i - 1, j - 1, k) - dvec(2, i - 1, j + 1, k))*xy &
                           + (dvec(3, i - 1, j, k - 1) - dvec(3, i - 1, j, k + 1))*xz

      END DO
   END DO

   ! Boundary nodes for y = 0.
   j = 1
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete d2Dx/dxdy + d2Dy/dy2 + d2Dz/dzdy using FFD for y and CFDs for x and z coordinate.
         res(2, i, j, k) = (dvec(1, i + 1, j + 1, k) - dvec(1, i - 1, j + 1, k))*xy &
                           + (dvec(2, i, j + 1, k) - dvec(2, i, j, k))*y2 &
                           + (dvec(3, i, j + 1, k + 1) - dvec(3, i, j + 1, k - 1))*yz

      END DO
   END DO

   ! Boundary nodes for y = Ly.
   j = npty
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete d2Dx/dxdy + d2Dy/dy2 + d2Dz/dzdy using BFD for y and CFDs for x and z coordinate.
         res(2, i, j, k) = (dvec(1, i - 1, j - 1, k) - dvec(1, i + 1, j - 1, k))*xy &
                           + (dvec(2, i, j - 1, k) - dvec(2, i, j, k))*y2 &
                           + (dvec(3, i, j - 1, k - 1) - dvec(3, i, j - 1, k + 1))*yz

      END DO
   END DO

   ! Boundary nodes for z = 0.
   k = 1
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete d2Dx/dxdz + d2Dy/dydz + d2Dz/dz2 using FFD for z and CFDs for x and y coordinate.
         res(3, i, j, k) = (dvec(1, i + 1, j, k + 1) - dvec(1, i - 1, j, k + 1))*xz &
                           + (dvec(2, i, j + 1, k + 1) - dvec(2, i, j - 1, k + 1))*yz &
                           + (dvec(3, i, j, k + 1) - dvec(3, i, j, k))*z2

      END DO
   END DO

   ! Boundary nodes for z = Lz.
   k = nptz
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete d2Dx/dxdz + d2Dy/dydz + d2Dz/dz2 using BFD for z and CFDs for x and y coordinate.
         res(3, i, j, k) = (dvec(1, i - 1, j, k - 1) - dvec(1, i + 1, j, k - 1))*xz &
                           + (dvec(2, i, j - 1, k - 1) - dvec(2, i, j + 1, k - 1))*yz &
                           + (dvec(3, i, j, k - 1) - dvec(3, i, j, k))*z2

      END DO
   END DO

END SUBROUTINE graddiv

! 3.4. A subroutine for calculating a laplacian from "phiv" array into "res" array.
! Note: Applies boundary conditions for Phi: Dirichlet boundary conditions.
SUBROUTINE laplas(res, phiv)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(nptx, npty, nptz), INTENT(IN)  :: phiv ! Input array.

   REAL(rk) :: hx, hy, hz ! Node spacings,
   REAL(rk) :: x2, y2, z2 ! Squared recipocals of node spacings.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   !INTEGER :: l       ! N/A.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of squared recipocals of node spacings.
   x2 = 1.0D0/hx**2
   y2 = 1.0D0/hy**2
   z2 = 1.0D0/hz**2

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1

            ! Discrete d2Phi/dx2 + d2Phi/dy2 + d2Phi/dz2 all using CFDs.
            res(i, j, k) = (phiv(i + 1, j, k) - 2.0D0*phiv(i, j, k) + phiv(i - 1, j, k))*x2 &
                           + (phiv(i, j + 1, k) - 2.0D0*phiv(i, j, k) + phiv(i, j - 1, k))*y2 &
                           + (phiv(i, j, k + 1) - 2.0D0*phiv(i, j, k) + phiv(i, j, k - 1))*z2

         END DO
      END DO
   END DO

END SUBROUTINE laplas

! 3.5. A subroutine for calculating a laplacian from "mfld" array into "res" array.
! Note: Applies boundary conditions for H: Normal components have Dirichlet boundary conditions, tangent components have Neumann boundary conditions.
SUBROUTINE laplash(res, mfld)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT none

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: mfld ! Input array.

   REAL(rk) :: hx, hy, hz ! Node spacings,
   REAL(rk) :: x2, y2, z2 ! Squared recipocals of node spacings.

   REAL(rk) :: k02 ! Squared wavenumber of free space.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   INTEGER :: l       ! Iterator variable for a field component.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of squared recipocals of node spacings.
   x2 = 1.0D0/hx**2
   y2 = 1.0D0/hy**2
   z2 = 1.0D0/hz**2

   ! Wavevector.
   k02 = (pi/size0)**2/2.0D0

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1
            DO l = 1, 3

               ! Discrete (d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2, d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2, d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2) all using CFDs.
               res(l, i, j, k) = (mfld(l, i + 1, j, k) - 2.0D0*mfld(l, i, j, k) + mfld(l, i - 1, j, k))*x2 &
                                 + (mfld(l, i, j + 1, k) - 2.0D0*mfld(l, i, j, k) + mfld(l, i, j - 1, k))*y2 &
                                 + (mfld(l, i, j, k + 1) - 2.0D0*mfld(l, i, j, k) + mfld(l, i, j, k - 1))*z2

            END DO
         END DO
      END DO
   END DO

   ! Boundary nodes for x = 0.
   i = 1
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using FFD for x coordinate and CFDs for y and z coordinates.
         res(2, i, j, k) = (mfld(2, i + 1, j, k) - mfld(2, i, j, k))*x2 &
                           + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(2, i, j, k + 1) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j, k - 1))*z2/2.0D0

         ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using FFD for x coordinate and CFDs for y and z coordinates.
         res(3, i, j, k) = (mfld(3, i + 1, j, k) - mfld(3, i, j, k))*x2 &
                           + (mfld(3, i, j + 1, k) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for x = Lx.
   i = nptx
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using BFD for x coordinate and CFDs for y and z coordinates.
         res(2, i, j, k) = (mfld(2, i - 1, j, k) - mfld(2, i, j, k))*x2 &
                           + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(2, i, j, k + 1) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j, k - 1))*z2/2.0D0

         ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using BFD for x coordinate and CFDs for y and z coordinates.
         res(3, i, j, k) = (mfld(3, i - 1, j, k) - mfld(3, i, j, k))*x2 &
                           + (mfld(3, i, j + 1, k) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = 0.
   j = 1
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using FFD for y coordinate and CFDs for x and z coordinates.
         res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(1, i, j + 1, k) - mfld(1, i, j, k))*y2 &
                           + (mfld(1, i, j, k + 1) - 2.0D0*mfld(1, i, j, k) + mfld(1, i, j, k - 1))*z2/2.0D0

         ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using FFD for y coordinate and CFDs for x and z coordinates.
         res(3, i, j, k) = (mfld(3, i + 1, j, k) - 2.0D0*mfld(3, i, j, k) + mfld(3, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(3, i, j + 1, k) - mfld(3, i, j, k))*y2 &
                           + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = Ly.
   j = npty
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using BFD for y coordinate and CFDs for x and z coordinates.
         res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(1, i, j - 1, k) - mfld(1, i, j, k))*y2 &
                           + (mfld(1, i, j, k + 1) - 2.0D0*mfld(1, i, j, k) + mfld(1, i, j, k - 1))*z2/2.0D0

         ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using BFD for y coordinate and CFDs for x and z coordinates.
         res(3, i, j, k) = (mfld(3, i + 1, j, k) - 2.0D0*mfld(3, i, j, k) + mfld(3, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(3, i, j - 1, k) - mfld(3, i, j, k))*y2 &
                           + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for z = 0.
   k = 1
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using FFD for z coordinate and CFDs for x and y coordinate.
         res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(1, i, j + 1, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(1, i, j, k + 1) - mfld(1, i, j, k))*z2

         ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using FFD for z coordinate and CFDs for x and y coordinate.
         res(2, i, j, k) = (mfld(2, i + 1, j, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(2, i, j, k + 1) - mfld(2, i, j, k))*z2

      END DO
   END DO

   ! Boundary nodes for z = Lz.
   k = nptz
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using BFD for z coordinate and CFDs for x and y coordinate.
         res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(1, i, j + 1, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(1, i, j, k - 1) - mfld(1, i, j, k))*z2

         ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using BFD for z coordinate and CFDs for x and y coordinate.
         res(2, i, j, k) = (mfld(2, i + 1, j, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i - 1, j, k))*x2/2.0D0 &
                           + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                           + (mfld(2, i, j, k - 1) - mfld(2, i, j, k))*z2

      END DO
   END DO

   ! Boundary nodes for x = 0, y = 0.
   i = 1
   j = 1
   DO k = 2, nptz - 1

      ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using FFDs for x and y coordinates and CFD for z coordinate.
      res(3, i, j, k) = (mfld(3, i + 1, j, k) - mfld(3, i, j, k))*x2 &
                      + (mfld(3, i, j + 1, k) - mfld(3, i, j, k))*y2 &
                      + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0!&
      !- k02*mfld(3, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = 0, y = Ly.
   i = 1
   j = npty
   DO k = 2, nptz - 1

      ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using FFD for x coordinate, BFD for y coordinate and CFD for z coordinate.
      res(3, i, j, k) = (mfld(3, i + 1, j, k) - mfld(3, i, j, k))*x2 &
                      + (mfld(3, i, j - 1, k) - mfld(3, i, j, k))*y2 &
                      + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0!&
      !- k02*mfld(3, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = Lx, y = 0.
   i = nptx
   j = 1
   DO k = 2, nptz - 1

      ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using BFD for x coordinate, FFD for y coordinate and CFD for z coordinate.
      res(3, i, j, k) = (mfld(3, i - 1, j, k) - mfld(3, i, j, k))*x2 &
                      + (mfld(3, i, j + 1, k) - mfld(3, i, j, k))*y2 &
                      + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0!&
      !- k02*mfld(3, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = Lx, y = Ly.
   i = nptx
   j = npty
   DO k = 2, nptz - 1

      ! Discrete d2Hz/dx2 + d2Hz/dy2 + d2Hz/dz2 using BFDs for x and y coordinates and CFD for z coordinate.
      res(3, i, j, k) = (mfld(3, i - 1, j, k) - mfld(3, i, j, k))*x2 &
                      + (mfld(3, i, j - 1, k) - mfld(3, i, j, k))*y2 &
                      + (mfld(3, i, j, k + 1) - 2.0D0*mfld(3, i, j, k) + mfld(3, i, j, k - 1))*z2/2.0D0!&
      !- k02*mfld(3, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = 0, z = 0.
   i = 1
   k = 1
   DO j = 2, npty - 1

      ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using FFDs for x and z coordinates and CFD for y coordinate.
      res(2, i, j, k ) = (mfld(2, i + 1, j, k) - mfld(2, i, j, k))*x2 &
                       + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                       + (mfld(2, i, j, k + 1) - mfld(2, i, j, k))*z2                                   !&
      !- k02*mfld(2, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = 0, z = Lz.
   i = 1
   k = nptz
   DO j = 2, npty - 1

      ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using FFD for x coordinate, BFD for z coordinate and CFD for y coordinate.
      res(2, i, j, k) = (mfld(2, i + 1, j, k) - mfld(2, i, j, k))*x2 &
                      + (mfld(2, i, j +1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                      + (mfld(2, i, j, k - 1) - mfld(2, i, j, k))*z2                                  !&
      !- k02*mfld(2, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = Lx, z = 0.
   i = nptx
   k = 1
   DO j = 2, npty - 1

      ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using BFD for x coordinate, FFD for z coordinate and CFD for y coordinate.
      res(2, i, j, k) = (mfld(2, i - 1, j, k) - mfld(2, i, j, k))*x2 &
                      + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                      + (mfld(2, i, j, k + 1) - mfld(2, i, j, k))*z2                                   !&
      !- k02*mfld(2, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for x = Lx, z = Lz.
   i = nptx
   k = nptz
   DO j = 2, npty - 1

      ! Discrete d2Hy/dx2 + d2Hy/dy2 + d2Hy/dz2 using BFDs for x and z coordinates and CFD for y coordinate.
      res(2, i, j, k) = (mfld(2, i - 1, j, k) - mfld(2, i, j, k))*x2 &
                      + (mfld(2, i, j + 1, k) - 2.0D0*mfld(2, i, j, k) + mfld(2, i, j - 1, k))*y2/2.0D0 &
                      + (mfld(2, i, j, k - 1) - mfld(2, i, j, k))*z2                                   !&
      !- k02*mfld(2, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for y = 0, z = 0.
   j = 1
   k = 1
   DO i = 2, nptx - 1

      ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using FFDs for y and z coordinates and CFD for x coordinate.
      res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                      + (mfld(1, i, j + 1, k) - mfld(1, i, j, k))*y2 &
                      + (mfld(1, i, j, k + 1) - mfld(1, i, j, k))*z2                                   !&
      !- k02*mfld(1, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for y = 0, z = Lz.
   j = 1
   k = nptz
   DO i = 2, nptx - 1

      ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using FFD for y coordinate, BFD for z coordinate and CFD for x coordinate.
      res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                      + (mfld(1, i, j + 1, k) - mfld(1, i, j, k))*y2 &
                      + (mfld(1, i, j, k - 1) - mfld(1, i, j, k))*z2                                   !&
      !- k02*mfld(1, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for y = Ly, z = 0.
   j = npty
   k = 1
   DO i = 2, nptx - 1

      ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using BFD for y coordinate, FFD for z coordinate and CFD for x coordinate.
      res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                      + (mfld(1, i, j - 1, k) - mfld(1, i, j, k))*y2 &
                      + (mfld(1, i, j, k + 1) - mfld(1, i, j, k))*z2                                   !&
      !- k02*mfld(1, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

   ! Boundary nodes for y = Ly, z = Lz.
   j = npty
   k = nptz
   DO i = 2, nptx - 1

      ! Discrete d2Hx/dx2 + d2Hx/dy2 + d2Hx/dz2 using BFDs for y and z coordinates and CFD for x coordinate.
      res(1, i, j, k) = (mfld(1, i + 1, j, k) - 2.0D0*mfld(1, i, j, k) + mfld(1, i - 1, j, k))*x2/2.0D0 &
                      + (mfld(1, i, j - 1, k) - mfld(1, i, j, k))*y2 &
                      + (mfld(1, i, j, k - 1) - mfld(1, i, j, k))*z2                                   !&
      !- k02*mfld(1, i, j, k) ! <--- BUG: Likely a remnant of a subroutine that computes D'Alambertanian of H

   END DO

END SUBROUTINE laplash

! 3.6. A subroutine for calculating a laplacian from "afld" array into "res" array.
! Note: Applies boundary conditions for A: Normal components have Neumann boundary conditions, tangent components have Dirichlet boundary conditions.
SUBROUTINE laplasa(res, afld)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: afld ! Input array.

   REAL(rk) :: hx, hy, hz ! Node spacings,
   REAL(rk) :: x2, y2, z2 ! Squared recipocals of node spacings.

   !REAL(rk) :: k02 ! Squared wavenumber of free space.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   INTEGER :: l       ! Iterator variable for a field component.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of squared recipocals of node spacings.
   x2 = 1.0D0/hx**2
   y2 = 1.0D0/hy**2
   z2 = 1.0D0/hz**2

   ! Wavevector.
   !k02 = (pi/size0)**2/2.0D0

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1
            DO l = 1, 3

               ! Discrete (d2Ax/dx2 + d2Ax/dy2 + d2Ax/dz2, d2Ay/dx2 + d2Ay/dy2 + d2Ay/dz2, d2Az/dx2 + d2Az/dy2 + d2Az/dz2) all using CFDs.
               res(l, i, j, k) = (afld(l, i + 1, j, k) - 2.0D0*afld(l, i, j, k) + afld(l, i - 1, j, k))*x2 &
                                 + (afld(l, i, j + 1, k) - 2.0D0*afld(l, i, j, k) + afld(l, i, j - 1, k))*y2 &
                                 + (afld(l, i, j, k + 1) - 2.0D0*afld(l, i, j, k) + afld(l, i, j, k - 1))*z2

            END DO
         END DO
      END DO
   END DO

   ! Boundary nodes for x = 0.
   i = 1
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete d2Ax/dx2 + d2Ax/dy2 + d2Ax/dz2 using FFD for x coordinate and CFDs for y and z coordinates.
         res(1, i, j, k) = (afld(1, i + 1, j, k) - afld(1, i, j, k))*x2 &
                           + (afld(1, i, j + 1, k) - 2.0D0*afld(1, i, j, k) + afld(1, i, j - 1, k))*y2/2.0D0 &
                           + (afld(1, i, j, k + 1) - 2.0D0*afld(1, i, j, k) + afld(1, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for x = Lx.
   i = nptx
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete d2Ax/dx2 + d2Ax/dy2 + d2Ax/dz2 using BFD for x coordinate and CFDs for y and z coordinates.
         res(1, i, j, k) = (afld(1, i - 1, j, k) - afld(1, i, j, k))*x2 &
                           + (afld(1, i, j + 1, k) - 2.0D0*afld(1, i, j, k) + afld(1, i, j - 1, k))*y2/2.0D0 &
                           + (afld(1, i, j, k + 1) - 2.0D0*afld(1, i, j, k) + afld(1, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = 0.
   j = 1
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete d2Ay/dx2 + d2Ay/dy2 + d2Ay/dz2 using FFD for y coordinate and CFDs for x and z coordinates.
         res(2, i, j, k) = (afld(2, i + 1, j, k) - 2.0D0*afld(2, i, j, k) + afld(2, i - 1, j, k))*x2/2.0D0 &
                           + (afld(2, i, j + 1, k) - afld(2, i, j, k))*y2 &
                           + (afld(2, i, j, k + 1) - 2.0D0*afld(2, i, j, k) + afld(2, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = Ly.
   j = npty
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete d2Ay/dx2 + d2Ay/dy2 + d2Ay/dz2 using BFD for y coordinate and CFDs for x and z coordinates.
         res(2, i, j, k) = (afld(2, i + 1, j, k) - 2.0D0*afld(2, i, j, k) + afld(2, i - 1, j, k))*x2/2.0D0 &
                           + (afld(2, i, j - 1, k) - afld(2, i, j, k))*y2 &
                           + (afld(2, i, j, k + 1) - 2.0D0*afld(2, i, j, k) + afld(2, i, j, k - 1))*z2/2.0D0

      END DO
   END DO

   ! Boundary nodes for z = 0.
   k = 1
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete d2Az/dx2 + d2Az/dy2 + d2Az/dz2 using FFD for z coordinate and CFDs for x and y coordinates.
         res(3, i, j, k) = (afld(3, i + 1, j, k) - 2.0D0*afld(3, i, j, k) + afld(3, i - 1, j, k))*x2/2.0D0 &
                           + (afld(3, i, j + 1, k) - 2.0D0*afld(3, i, j, k) + afld(3, i, j - 1, k))*y2/2.0D0 &
                           + (afld(3, i, j, k + 1) - afld(3, i, j, k))*z2

      END DO
   END DO

   ! Boundary nodes for z = Lz.
   k = nptz
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete d2Az/dx2 + d2Az/dy2 + d2Az/dz2 using BFD for z coordinate and CFDs for x and y coordinates.
         res(3, i, j, k) = (afld(3, i + 1, j, k) - 2.0D0*afld(3, i, j, k) + afld(3, i - 1, j, k))*x2/2.0D0 &
                           + (afld(3, i, j + 1, k) - 2.0D0*afld(3, i, j, k) + afld(3, i, j - 1, k))*y2/2.0D0 &
                           + (afld(3, i, j, k - 1) - afld(3, i, j, k))*z2

      END DO
   END DO

END SUBROUTINE laplasa

! 3.7. A subroutine for calculating a gradient from "phiv" array into "res" array.
! Note: Applies boundary conditions for Phi: Dirichlet boundary conditions.
SUBROUTINE nablaphi(res, phiv)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(nptx, npty, nptz), INTENT(IN)     :: phiv ! Input array.

   REAL(rk) :: hx, hy, hz    ! Node spacings,
   REAL(rk) :: x05, y05, z05 ! Halved recipocals of node spacings.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   !INTEGER :: l       ! N/A.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of halved recipocals of node spacings.
   x05 = 0.5D0/hx
   y05 = 0.5D0/hy
   z05 = 0.5D0/hz

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1

            ! Discrete (dPhi/dx, dPhi/dy, dPhi/dz) all using CFD.
            res(1, i, j, k) = (phiv(i + 1, j, k) - phiv(i - 1, j, k))*x05
            res(2, i, j, k) = (phiv(i, j + 1, k) - phiv(i, j - 1, k))*y05
            res(3, i, j, k) = (phiv(i, j, k + 1) - phiv(i, j, k - 1))*z05

         END DO
      END DO
   END DO

   ! Boundary nodes for x = 0.
   i = 1
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete (dPhi/dx, 0, 0) using FFD for x coordinate.
         res(1, i, j, k) = (phiv(i + 1, j, k) - phiv(i, j, k))*x05
         res(2, i, j, k) = (0.0D0, 0.0D0)
         res(3, i, j, k) = (0.0D0, 0.0D0)

      END DO
   END DO

   ! Boundary nodes for x = Lx.
   i = nptx
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete (dPhi/dx, 0, 0) using BFD for x coordinate.
         res(1, i, j, k) = (phiv(i, j, k) - phiv(i - 1, j, k))*x05
         res(2, i, j, k) = (0.0D0, 0.0D0)
         res(3, i, j, k) = (0.0D0, 0.0D0)

      END DO
   END DO

   ! Boundary nodes for y = 0.
   j = 1
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete (0, dPhi/dy, 0) using FFD for y coordinate.
         res(1, i, j, k) = (0.0D0, 0.0D0)
         res(2, i, j, k) = (phiv(i, j + 1, k) - phiv(i, j, k))*y05
         res(3, i, j, k) = (0.0D0, 0.0D0)

      END DO
   END DO

   ! Boundary nodes for y = Ly.
   j = npty
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete (0, dPhi/dy, 0) using BFD for y coordinate.
         res(1, i, j, k) = (0.0D0, 0.0D0)
         res(2, i, j, k) = (phiv(i, j, k) - phiv(i, j - 1, k))*y05
         res(3, i, j, k) = (0.0D0, 0.0D0)

      END DO
   END DO

   ! Boundary nodes for z = 0.
   k = 1
   DO j = 2, npty - 1
      DO i = 2, nptx - 1

         ! Discrete (0, 0, dPhi/dz) using FFD for z coordinate.
         res(1, i, j, k) = (0.0D0, 0.0D0)
         res(2, i, j, k) = (0.0D0, 0.0D0)
         res(3, i, j, k) = (phiv(i, j, k + 1) - phiv(i, j, k))*z05

      END DO
   END DO

   ! Boundary nodes for z = Lz.
   k = nptz
   DO j = 2, npty - 1
      DO i = 2, nptx - 1

         ! Discrete (0, 0, dPhi/dz) using BFD for z coordinate.
         res(1, i, j, k) = (0.0D0, 0.0D0)
         res(2, i, j, k) = (0.0D0, 0.0D0)
         res(3, i, j, k) = (phiv(i, j, k) - phiv(i, j, k - 1))*z05

      END DO
   END DO

END SUBROUTINE nablaphi

! 3.8. A subroutine for calculating a curl from "mfld" array into "res" array.
! Note: Applies boundary conditions for H: Normal components have Dirichlet boundary conditions, tangent components have Neumann boundary conditions.
SUBROUTINE roth(res, mfld)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: mfld ! Input array.

   REAL(rk) :: hx, hy, hz    ! Node spacings,
   REAL(rk) :: x05, y05, z05 ! Halved recipocals of node spacings.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   !INTEGER :: l       ! N/A.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of halved recipocals of node spacings.
   x05 = 0.5D0/hx
   y05 = 0.5D0/hy
   z05 = 0.5D0/hz

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1

            ! Discrete (dHz/dy - dHy/dz, dHx/dz - dHz/dx, dHy/dx - dHx/dy) all using CFDs.
            res(1, i, j, k) = (mfld(3, i, j + 1, k) - mfld(3, i, j - 1, k))*y05 - (mfld(2, i, j, k + 1) - mfld(2, i, j, k - 1))*z05
            res(2, i, j, k) = (mfld(1, i, j, k + 1) - mfld(1, i, j, k - 1))*z05 - (mfld(3, i + 1, j, k) - mfld(3, i - 1, j, k))*x05
            res(3, i, j, k) = (mfld(2, i + 1, j, k) - mfld(2, i - 1, j, k))*x05 - (mfld(1, i, j + 1, k) - mfld(1, i, j - 1, k))*y05

         END DO
      END DO
   END DO

   ! Boundary nodes for x = 0.
   i = 1
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete dHz/dy - dHy/dz all using CFD.
     res(1, i, j, k) = ((mfld(3, i, j + 1, k) - mfld(3, i, j - 1, k))*y05 - (mfld(2, i, j, k + 1) - mfld(2, i, j, k - 1))*z05)/2.0D0

      END DO
   END DO

   ! Boundary nodes for x = Lx.
   i = nptx
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete dHz/dy - dHy/dz all using CFD.
     res(1, i, j, k) = ((mfld(3, i, j + 1, k) - mfld(3, i, j - 1, k))*y05 - (mfld(2, i, j, k + 1) - mfld(2, i, j, k - 1))*z05)/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = 0.
   j = 1
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete dHx/dz - dHz/dx all using CFDs.
     res(2, i, j, k) = ((mfld(1, i, j, k + 1) - mfld(1, i, j, k - 1))*z05 - (mfld(3, i + 1, j, k) - mfld(3, i - 1, j, k))*x05)/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = Ly.
   j = npty
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete dHx/dz - dHz/dx all using CFDs.
     res(2, i, j, k) = ((mfld(1, i, j, k + 1) - mfld(1, i, j, k - 1))*z05 - (mfld(3, i + 1, j, k) - mfld(3, i - 1, j, k))*x05)/2.0D0

      END DO
   END DO

   ! Boundary nodes for z = 0.
   k = 1
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete dHy/dx - dHx/dy all using CFDs.
     res(3, i, j, k) = ((mfld(2, i + 1, j, k) - mfld(2, i - 1, j, k))*x05 - (mfld(1, i, j + 1, k) - mfld(1, i, j - 1, k))*y05)/2.0D0

      END DO
   END DO

   ! Boundary nodes for z = Lz.
   k = nptz
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete dHy/dx - dHx/dy all using CFDs.
     res(3, i, j, k) = ((mfld(2, i + 1, j, k) - mfld(2, i - 1, j, k))*x05 - (mfld(1, i, j + 1, k) - mfld(1, i, j - 1, k))*y05)/2.0D0

      END DO
   END DO

END SUBROUTINE roth

! 3.9. A subroutine for calculating a curl from "dfld" array into "res" array.
! Note: Applies boundary conditions for D: Normal components have Neumann boundary conditions, tangent components have Dirichlet boundary conditions.
SUBROUTINE rotd(res, dfld)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT) :: res  ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: dfld ! Input array.

   REAL(rk) :: hx, hy, hz    ! Node spacings,
   REAL(rk) :: x05, y05, z05 ! Halved recipocals of node spacings.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   !INTEGER :: l       ! N/A.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of halved recipocals of node spacings.
   x05 = 0.5D0/hx
   y05 = 0.5D0/hy
   z05 = 0.5D0/hz

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1

            ! Discrete (dDz/dy - dDy/dz, dDx/dz - dDz/dx, dDy/dx - dDx/dy) all using CFDs.
            res(1, i, j, k) = (dfld(3, i, j + 1, k) - dfld(3, i, j - 1, k))*y05 - (dfld(2, i, j, k + 1) - dfld(2, i, j, k - 1))*z05
            res(2, i, j, k) = (dfld(1, i, j, k + 1) - dfld(1, i, j, k - 1))*z05 - (dfld(3, i + 1, j, k) - dfld(3, i - 1, j, k))*x05
            res(3, i, j, k) = (dfld(2, i + 1, j, k) - dfld(2, i - 1, j, k))*x05 - (dfld(1, i, j + 1, k) - dfld(1, i, j - 1, k))*y05

         END DO
      END DO
   END DO

   ! Boundary nodes for x = 0.
   i = 1
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete dDx/dz - dDz/dx using FFD for x coordinate and CFD for z coordinate.
         res(2, i, j, k) = (dfld(1, i, j, k + 1) - dfld(1, i, j, k - 1))*z05/2.0D0 - dfld(3, i + 1, j, k)*x05

         ! Discrete dDy/dx - dDx/dy using FFD for x coordinate and CFD for y coordinate.
         res(3, i, j, k) = dfld(2, i + 1, j, k)*x05 - (dfld(1, i, j + 1, k) - dfld(1, i, j - 1, k))*y05/2.0D0

      END DO
   END DO

   ! Boundary nodes for x = Lx.
   i = nptx
   DO j = 2, npty - 1
      DO k = 2, nptz - 1

         ! Discrete dDx/dz - dDz/dx using BFD for x coordinate and CFD for z coordinate.
         res(2, i, j, k) = (dfld(1, i, j, k + 1) - dfld(1, i, j, k - 1))*z05/2.0D0 + dfld(3, i - 1, j, k)*x05

         ! Discrete dDy/dx - dDx/dy using BFD for x coordinate and CFD for y coordinate.
         res(3, i, j, k) = -dfld(2, i - 1, j, k)*x05 - (dfld(1, i, j + 1, k) - dfld(1, i, j - 1, k))*y05/2.0D0

      END DO
   END DO

   ! Boundary nodes for y = 0.
   j = 1
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete dDz/dy - dDy/dz using FFD for y coordinate and CFD for z coordinate.
         res(1, i, j, k) = dfld(3, i, j + 1, k)*y05 - (dfld(2, i, j, k + 1) - dfld(2, i, j, k - 1))*z05/2.0D0

         ! Discrete dDy/dx - dDx/dy using FFD for y coordinate and CFD for x coordinate.
         res(3, i, j, k) = (dfld(2, i + 1, j, k) - dfld(2, i - 1, j, k))*x05/2.0D0 - dfld(1, i, j + 1, k)*y05

      END DO
   END DO

   ! Boundary nodes for y = Ly.
   j = npty
   DO i = 2, nptx - 1
      DO k = 2, nptz - 1

         ! Discrete dDz/dy - dDy/dz using BFD for y coordinate and CFD for z coordinate.
         res(1, i, j, k) = -dfld(3, i, j - 1, k)*y05 - (dfld(2, i, j, k + 1) - dfld(2, i, j, k - 1))*z05/2.0D0

         ! Discrete dDy/dx - dDx/dy using BFD for y coordinate and CFD for x coordinate.
         res(3, i, j, k) = (dfld(2, i + 1, j, k) - dfld(2, i - 1, j, k))*x05/2.0D0 + dfld(1, i, j - 1, k)*y05

      END DO
   END DO

   ! Boundary nodes for z = 0.
   k = 1
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete dDz/dy - dDy/dz using FFD for z coordinate and CFD for y coordinate.
         res(1, i, j, k) = (dfld(3, i, j + 1, k) - dfld(3, i, j - 1, k))*y05/2.0D0 - dfld(2, i, j, k + 1)*z05

         ! Discrete dDx/dz - dDz/dx using FFD for z coordinate and CFD for x coordinate.
         res(2, i, j, k) = dfld(1, i, j, k + 1)*z05 - (dfld(3, i + 1, j, k) - dfld(3, i - 1, j, k))*x05/2.0D0

      END DO
   END DO

   ! Boundary nodes for z = Lz.
   k = nptz
   DO i = 2, nptx - 1
      DO j = 2, npty - 1

         ! Discrete dDz/dy - dDy/dz using BFD for z coordinate and CFD for y coordinate.
         res(1, i, j, k) = (dfld(3, i, j + 1, k) - dfld(3, i, j - 1, k))*y05/2.0D0 + dfld(2, i, j, k - 1)*z05

         ! Discrete dDx/dz - dDz/dx using BFD for z coordinate and CFD for x coordinate.
         res(2, i, j, k) = -dfld(1, i, j, k - 1)*z05 - (dfld(3, i + 1, j, k) - dfld(3, i - 1, j, k))*x05/2.0D0

      END DO
   END DO

   ! Boundary nodes for x = 0 and y = 0.
   i = 1
   j = 1
   DO k = 2, nptz - 1

      ! Discrete dDy/dx - dDx/dy using FFDs for x and y coordinates.
      res(3, i, j, k) = dfld(2, i + 1, j, k)*x05 - dfld(1, i, j + 1, k)*y05

   END DO

   ! Boundary nodes for x = 0 and y = Ly.
   i = 1
   j = npty
   DO k = 2, nptz - 1

      ! Discrete dDy/dx - dDx/dy FFD for x coordinate and BFD for y coordinate.
      res(3, i, j, k) = dfld(2, i + 1, j, k)*x05 + dfld(1, i, j - 1, k)*y05

   END DO

   ! Boundary nodes for x = Lx and y = 0.
   i = nptx
   j = 1
   DO k = 2, nptz - 1

      ! Discrete dDy/dx - dDx/dy using BFD for x coordinate and FFD for y coordinate.
      res(3, i, j, k) = -dfld(2, i - 1, j, k)*x05 - dfld(1, i, j + 1, k)*y05

   END DO

   ! Boundary nodes for x = Lx and y = Ly.
   i = nptx
   j = npty
   DO k = 2, nptz - 1

      ! Discrete dDy/dx - dDx/dy using BFDs for x and y coordinates.
      res(3, i, j, k) = -dfld(2, i - 1, j, k)*x05 + dfld(1, i, j - 1, k)*y05

   END DO

   ! Boundary nodes for x = 0 and z = 0.
   i = 1
   k = 1
   DO j = 2, npty - 1

      ! Discrete dDx/dz - dDz/dx using FFDs for x and z coordinates.
      res(2, i, j, k) = dfld(1, i, j, k + 1)*z05 - dfld(3, i + 1, j, k)*x05

   END DO

   ! Boundary nodes for x = 0 and z = Lz.
   i = 1
   k = nptz
   DO j = 2, npty - 1

      ! Discrete dDx/dz - dDz/dx using FFD for x coordinate and BFD for z coordinate.
      res(2, i, j, k) = -dfld(1, i, j, k - 1)*z05 - dfld(3, i + 1, j, k)*x05

   END DO

   ! Boundary nodes for x = Lx and z = 0.
   i = nptx
   k = 1
   DO j = 2, npty - 1

      ! Discrete dDx/dz - dDz/dx using BFD for x coordinate and FFD for z coordinate.
      res(2, i, j, k) = dfld(1, i, j, k + 1)*z05 + dfld(3, i - 1, j, k)*x05

   END DO

   ! Boundary nodes for x = Lx and z = Lz.
   i = nptx
   k = nptz
   DO j = 2, npty - 1

      ! Discrete dDx/dz - dDz/dx using BFDs for x and z coordinates.
      res(2, i, j, k) = -dfld(1, i, j, k - 1)*z05 + dfld(3, i - 1, j, k)*x05

   END DO

   ! Boundary nodes for y = 0 and z = 0.
   j = 1
   k = 1
   DO i = 2, nptx - 1

      ! Discrete dDz/dy - dDy/dz using FFDs for y and z coordinates.
      res(1, i, j, k) = dfld(3, i, j + 1, k)*y05 - dfld(2, i, j, k + 1)*z05

   END DO

   ! Boundary nodes for y = 0 and z = Lz.
   j = 1
   k = nptz
   DO i = 2, nptx - 1

      ! Discrete dDz/dy - dDy/dz using FFD for y coordinate and BFD for z coordinate.
      res(1, i, j, k) = dfld(3, i, j + 1, k)*y05 + dfld(2, i, j, k - 1)*z05

   END DO

   ! Boundary nodes for y = 0 and z = Lz.
   j = npty
   k = 1
   DO i = 2, nptx - 1

      ! Discrete dDz/dy - dDy/dz using BFD for y coordinate and FFD for z coordinate.
      res(1, i, j, k) = -dfld(3, i, j - 1, k)*y05 - dfld(2, i, j, k + 1)*z05

   END DO

   ! Boundary nodes for y = Ly and z = Lz.
   j = npty
   k = nptz
   DO i = 2, nptx - 1

      ! Discrete dDz/dy - dDy/dz using BFDs for y and z coordinates.
      res(1, i, j, k) = -dfld(3, i, j - 1, k)*y05 + dfld(2, i, j, k - 1)*z05

   END DO

END SUBROUTINE rotd

! 3.10. A subroutine for calculating a divergence from "dvec" array into "res" array.
! Note: Applies boundary conditions for D: Normal components have Neumann boundary conditions, tangent components have Dirichlet boundary conditions.
SUBROUTINE divd(res, dvec)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(nptx, npty, nptz), INTENT(OUT)   :: res  ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN) :: dvec ! Input array.

   REAL(rk) :: hx, hy, hz    ! Node spacings,
   REAL(rk) :: x05, y05, z05 ! Halved recipocals of node spacings.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   !INTEGER :: l       ! N/A.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Calculation of node spacings.
   hx = sizex/(nptx - 1.0D0)
   hy = sizey/(npty - 1.0D0)
   hz = sizez/(nptz - 1.0D0)

   ! Calculation of halved recipocals of node spacings.
   x05 = 0.5D0/hx
   y05 = 0.5D0/hy
   z05 = 0.5D0/hz

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   !===================!
   ! Array assignment. !
   !===================!

   ! Interior nodes.
   DO i = 2, nptx - 1
      DO j = 2, npty - 1
         DO k = 2, nptz - 1

            ! Discrete dDx/dx + dDy/dy + dDz/dz all using CFDs.
            res(i, j, k) = (dvec(1, i + 1, j, k) - dvec(1, i - 1, j, k))*x05 &
                           + (dvec(2, i, j + 1, k) - dvec(2, i, j - 1, k))*y05 &
                           + (dvec(3, i, j, k + 1) - dvec(3, i, j, k - 1))*z05

         END DO
      END DO
   END DO

END SUBROUTINE divd

! 3.11. A subroutine for calculating the operator action U(x), where x = (Dx, Dy, Dz, Phi, Hx, Hy, Hz, Ax, Ay, Az) is the unknown vector.
SUBROUTINE mulmat(arg, res, tempr)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN)   :: arg   ! Input array,
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(OUT)  :: res   ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(INOUT) :: tempr ! Temporary array.

   COMPLEX(rk) :: ik0 ! Wavenumber of the free space multiplied by an imaginary unit.
   REAL(rk)    :: k02 ! Wavenumber of the free space squared.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Wavevector.
   ik0 = (0.0D0, 1.0D0)*(pi/size0)
   k02 = (pi/size0)**2

   !=================================================================!
   ! A set of subroutines to calculate the governing equation for D. !
   !=================================================================!

   CALL graddiv(res(1:3, :, :, :), arg(1:3, :, :, :)) ! U(x){1:3} = grad(div(D))
   CALL nablaphi(tempr, arg(4, :, :, :))              ! tempr     = grad(Phi)
   res(1:3, :, :, :) = res(1:3, :, :, :) + ik0*tempr  ! U(x){1:3} = grad(div(D)) + i.k0.grad(Phi)
   CALL roth(tempr, arg(5:7, :, :, :))                ! tempr     = curl(H)
   res(1:3, :, :, :) = res(1:3, :, :, :) + ik0*tempr  ! U(x){1:3} = grad(div(D)) + i.k0.curl(H) + i.k0.grad(Phi)
   tempr = -res(1:3, :, :, :) &                      ! tempr = - grad(div(D)) + 2.k0^2.D - k0^2.A - i.k0.curl(H) - i.k0.grad(Phi)
           + 2.0D0*k02*arg(1:3, :, :, :) - k02*arg(8:10, :, :, :)
   res(1:3, :, :, :) = tempr                          ! U(x){1:3} = - grad(div(D)) + 2.k0^2.D - k0^2.A - i.k0.curl(H) - i.k0.grad(Phi)

   ! A halving adjustment of non-zero mass term Dx at boundary nodes for x = 0 and x = Lx.
   res(1, 1, :, :) = res(1, 1, :, :) - k02*arg(1, 1, :, :) + 0.5D0*k02*arg(8, 1, :, :)
   res(1, nptx, :, :) = res(1, nptx, :, :) - k02*arg(1, nptx, :, :) + 0.5D0*k02*arg(8, nptx, :, :)

   ! A halving adjustment of non-zero mass term Dy at boundary nodes for y = 0 and y = Ly.
   res(2, :, 1, :) = res(2, :, 1, :) - k02*arg(2, :, 1, :) + 0.5D0*k02*arg(9, :, 1, :)
   res(2, :, npty, :) = res(2, :, npty, :) - k02*arg(2, :, npty, :) + 0.5D0*k02*arg(9, :, npty, :)

   ! A halving adjustment of non-zero mass term Dz at boundary nodes for z = 0 and z = Lz.
   res(3, :, :, 1) = res(3, :, :, 1) - k02*arg(3, :, :, 1) + 0.5D0*k02*arg(10, :, :, 1)
   res(3, :, :, nptz) = res(3, :, :, nptz) - k02*arg(3, :, :, nptz) + 0.5D0*k02*arg(10, :, :, nptz)

   !===================================================================!
   ! A set of subroutines to calculate the governing equation for Phi. !
   !===================================================================!

   CALL laplas(res(4, :, :, :), arg(4, :, :, :))             ! U(x){4}  = laplacian(Phi)
   CALL divd(tempr(1, :, :, :), arg(1:3, :, :, :))           ! tempr{1} = div(D)
   res(4, :, :, :) = res(4, :, :, :) + ik0*tempr(1, :, :, :) ! U(x){4}  = laplacian(Phi) + i.k0.div(D)
   CALL divd(tempr(1, :, :, :), arg(8:10, :, :, :))          ! tempr{1} = div(A)
   res(4, :, :, :) = res(4, :, :, :) - ik0*tempr(1, :, :, :) ! U(x){4}  = laplacian(Phi) - i.k0.div(A) + i.k0.div(D)
   tempr(1, :, :, :) = -res(4, :, :, :)                     ! tempr{1} = - laplacian(Phi) + i.k0.div(A) - i.k0.div(D)
   res(4, :, :, :) = tempr(1, :, :, :)                       ! U(x){4}  = - laplacian(Phi) + i.k0.div(A) - i.k0.div(D)

   !=================================================================!
   ! A set of subroutines to calculate the governing equation for H. !
   !=================================================================!

   CALL laplash(res(5:7, :, :, :), arg(5:7, :, :, :))  ! U(x){5:7} = laplacian(H)
   CALL rotd(tempr, arg(1:3, :, :, :))                 ! tempr     = curl(D)
   res(5:7, :, :, :) = res(5:7, :, :, :) - ik0*tempr   ! U(x){5:7} = laplacian(H) - i.k0.curl(D)
   CALL rotd(tempr, arg(8:10, :, :, :))                ! tempr     = curl(A)
   res(5:7, :, :, :) = res(5:7, :, :, :) - ik0*tempr   ! U(x){5:7} = laplacian(H) - i.k0.curl(A) - i.k0.curl(D)
   tempr = -res(5:7, :, :, :) + k02*arg(5:7, :, :, :) ! tempr     = - laplacian(H) + k0^2.H + i.k0.curl(A) + i.k0.curl(D)
   res(5:7, :, :, :) = tempr                           ! U(x){5:7} = - laplacian(H) + k0^2.H + i.k0.curl(A) + i.k0.curl(D)

   ! A halving adjustment of non-zero mass term Hx in H governing equation at boundary nodes for y = 0 and y = Ly.
   res(5, :, 1, :) = res(5, :, 1, :) - 0.5D0*k02*arg(5, :, 1, :)
   res(5, :, npty, :) = res(5, :, npty, :) - 0.5D0*k02*arg(5, :, npty, :)

   ! A halving adjustment of non-zero mass term Hx in H governing equation at boundary nodes for z = 0 and z = Ly.
   res(5, :, :, 1) = res(5, :, :, 1) - 0.5D0*k02*arg(5, :, :, 1)
   res(5, :, :, nptz) = res(5, :, :, nptz) - 0.5D0*k02*arg(5, :, :, nptz)

   ! A halving adjustment of non-zero mass term Hy in H governing equation at boundary nodes for x = 0 and x = Lx.
   res(6, 1, :, :) = res(6, 1, :, :) - 0.5D0*k02*arg(6, 1, :, :)
   res(6, nptx, :, :) = res(6, nptx, :, :) - 0.5D0*k02*arg(6, nptx, :, :)

   ! A halving adjustment of non-zero mass term Hy in H governing equation at boundary nodes for z = 0 and z = Lz.
   res(6, :, :, 1) = res(6, :, :, 1) - 0.5D0*k02*arg(6, :, :, 1)
   res(6, :, :, nptz) = res(6, :, :, nptz) - 0.5D0*k02*arg(6, :, :, nptz)

   ! A halving adjustment of non-zero mass term Hz in H governing equation at boundary nodes for x = 0 and x = Lx.
   res(7, 1, :, :) = res(7, 1, :, :) - 0.5D0*k02*arg(7, 1, :, :)
   res(7, nptx, :, :) = res(7, nptx, :, :) - 0.5D0*k02*arg(7, nptx, :, :)

   ! A halving adjustment of non-zero mass term Hz in H governing equation at boundary nodes for y = 0 and y = Ly.
   res(7, :, 1, :) = res(7, :, 1, :) - 0.5D0*k02*arg(7, :, 1, :)
   res(7, :, npty, :) = res(7, :, npty, :) - 0.5D0*k02*arg(7, :, npty, :)

   ! A halving adjustment of x-th component of H governing equation. ! <--- Might be a bug a BUG.
   !res(5, :, 1, 1) = res(5, :, 1, 1)*0.5D0             ! .. at boundary nodes for y = 0, z = 0.
   !res(5, :, npty, 1) = res(5, :, npty, 1)*0.5D0       ! .. at boundary nodes for y = Ly, z = 0.
   !res(5, :, 1, nptz) = res(5, :, 1, nptz)*0.5D0       ! .. at boundary nodes for y = 0, z = Lz.
   !res(5, :, npty, nptz) = res(5, :, npty, nptz)*0.5D0 ! .. at boundary nodes for y = Ly, z = Lz.

   ! A halving adjustment of y-th component of H governing equation. ! <--- Might be a bug a BUG.
   !res(6, 1, :, 1) = res(6, 1, :, 1)*0.5D0             ! .. at boundary nodes for x = 0, z = 0.
   !res(6, nptx, :, 1) = res(6, nptx, :, 1)*0.5D0       ! .. at boundary nodes for x = Lx, z = 0.
   !res(6, 1, :, nptz) = res(6, 1, :, nptz)*0.5D0       ! .. at boundary nodes for x = 0, z = Lz.
   !res(6, nptx, :, nptz) = res(6, nptx, :, nptz)*0.5D0 ! .. at boundary nodes for x = Lx, z = Lz.

   ! A halving adjustment of z-th component of H governing equation. ! <--- Might be a bug a BUG.
   !res(7, 1, 1, :) = res(7, 1, 1, :)*0.5D0             ! .. at boundary nodes for x = 0, y = 0.
   !res(7, nptx, 1, :) = res(7, nptx, 1, :)*0.5D0       ! .. at boundary nodes for x = Lx, y = 0.
   !res(7, 1, npty, :) = res(7, 1, npty, :)*0.5D0       ! .. at boundary nodes for x = 0, y = Ly.
   !res(7, nptx, npty, :) = res(7, nptx, npty, :)*0.5D0 ! .. at boundary nodes for x = Lx, y = Ly.

   ! A halving adjustment of x-th component of H governing equation.
   res(5, :, 1, 1) = res(5, :, 1, 1) + 0.25D0*k02*arg(5, :, 1, 1)                   ! .. at boundary nodes for y = 0, z = 0.
   res(5, :, npty, 1) = res(5, :, npty, 1) + 0.25D0*k02*arg(5, :, npty, 1)          ! .. at boundary nodes for y = Ly, z = 0.
   res(5, :, 1, nptz) = res(5, :, 1, nptz) + 0.25D0*k02*arg(5, :, 1, nptz)          ! .. at boundary nodes for y = 0, z = Lz.
   res(5, :, npty, nptz) = res(5, :, npty, nptz) + 0.25D0*k02*arg(5, :, npty, nptz) ! .. at boundary nodes for y = Ly, z = Lz.

   ! A halving adjustment of y-th component of H governing equation.
   res(6, 1, :, 1) = res(6, 1, :, 1) + 0.25D0*k02*arg(6, 1, :, 1)                   ! .. at boundary nodes for x = 0, z = 0.
   res(6, nptx, :, 1) = res(6, nptx, :, 1) + 0.25D0*k02*arg(6, nptx, :, 1)          ! .. at boundary nodes for x = Lx, z = 0.
   res(6, 1, :, nptz) = res(6, 1, :, nptz) + 0.25D0*k02*arg(6, 1, :, nptz)          ! .. at boundary nodes for x = 0, z = Lz.
   res(6, nptx, :, nptz) = res(6, nptx, :, nptz) + 0.25D0*k02*arg(6, nptx, :, nptz) ! .. at boundary nodes for x = Lx, z = Lz.

   ! A halving adjustment of z-th component of H governing equation.
   res(7, 1, 1, :) = res(7, 1, 1, :) + 0.25D0*k02*arg(7, 1, 1, :)                   ! .. at boundary nodes for x = 0, y = 0.
   res(7, nptx, 1, :) = res(7, nptx, 1, :) + 0.25D0*k02*arg(7, nptx, 1, :)          ! .. at boundary nodes for x = Lx, y = 0.
   res(7, 1, npty, :) = res(7, 1, npty, :) + 0.25D0*k02*arg(7, 1, npty, :)          ! .. at boundary nodes for x = 0, y = Ly.
   res(7, nptx, npty, :) = res(7, nptx, npty, :) + 0.25D0*k02*arg(7, nptx, npty, :) ! .. at boundary nodes for x = Lx, y = Ly.

   !=================================================================!
   ! A set of subroutines to calculate the governing equation for A. !
   !=================================================================!

   CALL laplasa(res(8:10, :, :, :), arg(8:10, :, :, :))    ! U(x){8:10} = laplacian(A)
   CALL roth(tempr, arg(5:7, :, :, :))                     ! tempr      = curl(H)
   res(8:10, :, :, :) = res(8:10, :, :, :) + ik0*tempr     ! U(x){8:10} = laplacian(A) + i.k0.curl(H)
   CALL nablaphi(tempr, arg(4, :, :, :))                   ! tempr      = grad(Phi)
   res(8:10, :, :, :) = res(8:10, :, :, :) - ik0*tempr     ! U(x){8:10} = laplacian(A) + i.k0.curl(H) - i.k0.grad(Phi)
   tempr = -res(8:10, :, :, :) - k02*arg(1:3, :, :, :) &  ! tempr      = - laplacian(A) + k0^2.A - i.k0.curl(H) - k0^2.D + i.k0.grad(Phi)
           + k02*arg(8:10, :, :, :)
   res(8:10, :, :, :) = tempr                              ! U(x){8:10} = - laplacian(A) + k0^2.A - i.k0.curl(H) - k0^2.D + i.k0.grad(Phi)

   ! A halving adjustment of non-zero mass term Dx - Ax at boundary nodes for x = 0 and x = Lx.
   res(8, 1, :, :) = res(8, 1, :, :) - 0.5D0*k02*arg(8, 1, :, :) + 0.5D0*k02*arg(1, 1, :, :)
   res(8, nptx, :, :) = res(8, nptx, :, :) - 0.5D0*k02*arg(8, nptx, :, :) + 0.5D0*k02*arg(1, nptx, :, :)

   ! A halving adjustment of non-zero mass term Dy - Ay at boundary nodes for y = 0 and y = Ly.
   res(9, :, 1, :) = res(9, :, 1, :) - 0.5D0*k02*arg(9, :, 1, :) + 0.5D0*k02*arg(2, :, 1, :)
   res(9, :, npty, :) = res(9, :, npty, :) - 0.5D0*k02*arg(9, :, npty, :) + 0.5D0*k02*arg(2, :, npty, :)

   ! A halving adjustment of non-zero mass term Dz - Az at boundary nodes for z = 0 and z = Lz.
   res(10, :, :, 1) = res(10, :, :, 1) - 0.5D0*k02*arg(10, :, :, 1) + 0.5D0*k02*arg(3, :, :, 1)
   res(10, :, :, nptz) = res(10, :, :, nptz) - 0.5D0*k02*arg(10, :, :, nptz) + 0.5D0*k02*arg(3, :, :, nptz)

   !res(2,:,:,:)=(0.0D0,0.0D0)
   !res(1:3,:,:,:)=(0.0D0,0.0D0)
   !res(4,:,:,:)=(0.0D0,0.0D0)

END SUBROUTINE mulmat

! A SUBROUTINE TO CALCULATE U(X), note: checks on mulmat revealed some issues and since i cannot directly access the stack, i'll have to rewrite to ensure it is correct.
SUBROUTINE mul_l_vacuum(res, tmp, arg)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   ! Interface variables.
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(OUT)   :: res
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz),  INTENT(INOUT) :: tmp
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN)    :: arg

   ! Local variables.
   COMPLEX(rk) :: ik0
   REAL(rk)    :: k02

   ! Assigning runtime variables.
   ik0 = (0.0D0, 1.0D0)*pi/size0
   k02 = ABS(ik0)**2

   ! Initial assignment.
   res = (0.0D0, 0.0D0)

   ! Governing equation for D field.
   CALL graddiv(tmp, arg(1:3, :, :, :))                                            ! tmp = grad(div(D))
   res(1:3, :, :, :) = -tmp + 2.0D0*k02*arg(1:3, :, :, :) - k02*arg(8:10, :, :, :) ! res(1:3) = -grad(div(D)) + 2k0^2*D - k0^2*A
   CALL roth(tmp, arg(5:7, :, :, :))                                               ! tmp = curl(H)
   res(1:3, :, :, :) = res(1:3, :, :, :) - ik0*tmp                                 ! res(1:3) = -grad(div(D)) + 2k0^2*D - k0^2*A - ik0*curl(H)
   CALL nablaphi(tmp, arg(4, :, :, :))                                             ! tmp = grad(Phi)
   res(1:3, :, :, :) = res(1:3, :, :, :) - ik0*tmp                                 ! res(1:3) = -grad(div(D)) + 2k0^2*D - k0^2*A - ik0*curl(H) - ik0*grad(Phi)

   ! Correction of mass terms near boundaries.
   res(1, 1, :, :)    = res(1, 1, :, :)    - k02*arg(1, 1, :, :)    + 0.5D0*k02*arg(8, 1, :, :)
   res(1, nptx, :, :) = res(1, nptx, :, :) - k02*arg(1, nptx, :, :) + 0.5D0*k02*arg(8, nptx, :, :)
   res(2, :, 1, :)    = res(2, :, 1, :)    - k02*arg(2, :, 1, :)    + 0.5D0*k02*arg(9, :, 1, :)
   res(2, :, npty, :) = res(2, :, npty, :) - k02*arg(2, :, npty, :) + 0.5D0*k02*arg(9, :, npty, :)
   res(3, :, :, 1)    = res(3, :, :, 1)    - k02*arg(3, :, :, 1)    + 0.5D0*k02*arg(10, :, :, 1)
   res(3, :, :, nptz) = res(3, :, :, nptz) - k02*arg(3, :, :, nptz) + 0.5D0*k02*arg(10, :, :, nptz)

   ! Governing equation for Phi field.
   CALL laplas(tmp(1, :, :, :), arg(4, :, :, :))                                  ! tmp(1) = laplacian(Phi)
   CALL divd(tmp(2, :, :, :), arg(8:10, :, :, :))                                 ! tmp(2) = div(A)
   CALL divd(tmp(3, :, :, :), arg(1:3, :, :, :))                                  ! tmp(3) = div(D)
   res(4, :, :, :) = -tmp(1, :, :, :) + ik0*tmp(2, :, :, :) - ik0*tmp(3, :, :, :) ! res(4) = -laplacian(Phi) + ik0*div(A) - ik0*div(D)

   ! Governing equation for H field.
   CALL laplash(tmp, arg(5:7, :, :, :))             ! tmp = laplacian(H)
   res(5:7, :, :, :) = -tmp + k02*arg(5:7, :, :, :) ! res(5:7) = -laplacian(H) + k0^2*H
   CALL rotd(tmp, arg(8:10, :, :, :))               ! tmp = curl(A)
   res(5:7, :, :, :) = res(5:7, :, :, :) + ik0*tmp  ! res(5:7) = -laplacian(H) + k0^2*H + ik0*curl(A)
   CALL rotd(tmp, arg(1:3, :, :, :))                ! tmp = curl(D)
   res(5:7, :, :, :) = res(5:7, :, :, :) + ik0*tmp  ! res(5:7) = -laplacian(H) + k0^2*H + ik0*curl(A) + ik0*curl(D)

   ! Correction of mass terms near boundaries.
   res(5, :, 1, :)    = res(5, :, 1, :)    - 0.5D0*k02*arg(5, :, 1, :)
   res(5, :, npty, :) = res(5, :, npty, :) - 0.5D0*k02*arg(5, :, npty, :)
   res(5, :, :, 1)    = res(5, :, :, 1)    - 0.5D0*k02*arg(5, :, :, 1)
   res(5, :, :, nptz) = res(5, :, :, nptz) - 0.5D0*k02*arg(5, :, :, nptz)
   res(6, 1, :, :)    = res(6, 1, :, :)    - 0.5D0*k02*arg(6, 1, :, :)
   res(6, nptx, :, :) = res(6, nptx, :, :) - 0.5D0*k02*arg(6, nptx, :, :)
   res(6, :, :, 1)    = res(6, :, :, 1)    - 0.5D0*k02*arg(6, :, :, 1)
   res(6, :, :, nptz) = res(6, :, :, nptz) - 0.5D0*k02*arg(6, :, :, nptz)
   res(7, 1, :, :)    = res(7, 1, :, :)    - 0.5D0*k02*arg(7, 1, :, :)
   res(7, nptx, :, :) = res(7, nptx, :, :) - 0.5D0*k02*arg(7, nptx, :, :)
   res(7, :, 1, :)    = res(7, :, 1, :)    - 0.5D0*k02*arg(7, :, 1, :)
   res(7, :, npty, :) = res(7, :, npty, :) - 0.5D0*k02*arg(7, :, npty, :)

   res(5, :, 1, 1)       = res(5, :, 1, 1)       + 0.25D0*k02*arg(5, :, 1, 1)
   res(5, :, npty, 1)    = res(5, :, npty, 1)    + 0.25D0*k02*arg(5, :, npty, 1)
   res(5, :, 1, nptz)    = res(5, :, 1, nptz)    + 0.25D0*k02*arg(5, :, 1, nptz)
   res(5, :, npty, nptz) = res(5, :, npty, nptz) + 0.25D0*k02*arg(5, :, npty, nptz)
   res(6, 1, :, 1)       = res(6, 1, :, 1)       + 0.25D0*k02*arg(6, 1, :, 1)
   res(6, nptx, :, 1)    = res(6, nptx, :, 1)    + 0.25D0*k02*arg(6, nptx, :, 1)
   res(6, 1, :, nptz)    = res(6, 1, :, nptz)    + 0.25D0*k02*arg(6, 1, :, nptz)
   res(6, nptx, :, nptz) = res(6, nptx, :, nptz) + 0.25D0*k02*arg(6, nptx, :, nptz)
   res(7, 1, 1, :)       = res(7, 1, 1, :)       + 0.25D0*k02*arg(7, 1, 1, :)
   res(7, nptx, 1, :)    = res(7, nptx, 1, :)    + 0.25D0*k02*arg(7, nptx, 1, :)
   res(7, 1, npty, :)    = res(7, 1, npty, :)    + 0.25D0*k02*arg(7, 1, npty, :)
   res(7, nptx, npty, :) = res(7, nptx, npty, :) + 0.25D0*k02*arg(7, nptx, npty, :)

   ! Governing equation for A field.
   CALL laplasa(tmp, arg(8:10, :, :, :))                           ! tmp = laplacian(A)
   res(8:10, :, :, :) = -tmp + k02*arg(8:10, :, :, :)              ! res(8:10) = -laplacian(A) + k0^2*A
   CALL roth(tmp, arg(5:7, :, :, :))                               ! tmp = curl(H)
   res(8:10, :, :, :) = res(8:10, :, :, :) - ik0*tmp               ! res(8:10) = -laplacian(A) + k0^2*A - ik0*curl(H)
   res(8:10, :, :, :) = res(8:10, :, :, :) - k02*arg(1:3, :, :, :) ! res(8:10) = -laplacian(A) + k0^2*A - ik0*curl(H) - k0^2*D
   CALL nablaphi(tmp, arg(4, :, :, :))                             ! tmp = grad(Phi)
   res(8:10, :, :, :) = res(8:10, :, :, :) + ik0*tmp               ! res(8:10) = -laplacian(A) + k0^2*A - ik0*curl(H) - k0^2*D +ik*grad(Phi) 

   ! Correction of mass terms near boundaries.
   res(8, 1, :, :)     = res(8, 1, :, :)     - 0.5D0*k02*arg(8, 1, :, :)     + 0.5D0*k02*arg(1, 1, :, :)
   res(8, nptx, :, :)  = res(8, nptx, :, :)  - 0.5D0*k02*arg(8, nptx, :, :)  + 0.5D0*k02*arg(1, nptx, :, :)
   res(9, :, 1, :)     = res(9, :, 1, :)     - 0.5D0*k02*arg(9, :, 1, :)     + 0.5D0*k02*arg(2, :, 1, :)
   res(9, :, npty, :)  = res(9, :, npty, :)  - 0.5D0*k02*arg(9, :, npty, :)  + 0.5D0*k02*arg(2, :, npty, :)
   res(10, :, :, 1)    = res(10, :, :, 1)    - 0.5D0*k02*arg(10, :, :, 1)    + 0.5D0*k02*arg(3, :, :, 1)
   res(10, :, :, nptz) = res(10, :, :, nptz) - 0.5D0*k02*arg(10, :, :, nptz) + 0.5D0*k02*arg(3, :, :, nptz)

END SUBROUTINE mul_l_vacuum

! 3.12. A subroutine for calculating the discrete inner product of two 10D vectors.
! Note: Normally a tripple looped sum is defined with a single DOT_PRODUCT instruction, however here it is defined as a double looped sum with 10 DOT_PRODUCT instructions.
!       Subroutine is equivalent to standard summation in this sense but has a downside: DOT_PRODUCT works with larger arrays, which for fine mesh settings makes it computationally more expensive as referencing 20.Nx elements of two Nx.Ny.Nz arrays at a time is slower than referencing just 20 elements of two Nx.Ny.Nz arrays at a time.
!       In theory if compiler vectorizes DOT_PRODUCT by default this would boost computation but in practice this does not scale with finer mesh configurations and only provides an advantage for coarse meshes. TODO: Consider implementing the usual way of computing the discrete inner product.
FUNCTION scprod(left, right)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: left, right ! Input arrays.

   !REAL(rk) :: scprod ! Return variable. ! <--- BUG.
   COMPLEX(rk) :: scprod ! Return variable.

   !INTEGER :: i, j, k ! Iterator variables for x, y, and z. ! <--- Added a skeleton of a more efficient implementation.
   INTEGER :: j, k ! Iterator variables for y and z.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Array initialization.
   scprod = (0.0D0, 0.0D0)

   !===========================!
   ! Accumulative calculation. !
   !===========================!

   ! Double looped summation. (Ineficient for fine meshes.)
   !DO i = 1, nptx ! <--- Added a skeleton of a more efficient implementation.
   DO j = 1, npty
      DO k = 1, nptz

         ! Discrete inner product of x.x = Dx*.Dx + Dy*.Dy + Dz*.Dz + Phi*.Phi + Hx*.Hx + Hy*.Hy + Hz*.Hz + Ax*.Ax + Ay*.Ay + Az*.Az.
         !scprod = scprod(left(:, i, j, k), right(:, i, j, k)) ! <--- Added a skeleton of a more efficient implementation.
         scprod = scprod &
                  + DOT_PRODUCT(left(1, :, j, k), right(1, :, j, k)) &
                  + DOT_PRODUCT(left(2, :, j, k), right(2, :, j, k)) &
                  + DOT_PRODUCT(left(3, :, j, k), right(3, :, j, k)) &
                  + DOT_PRODUCT(left(4, :, j, k), right(4, :, j, k)) &
                  + DOT_PRODUCT(left(5, :, j, k), right(5, :, j, k)) &
                  + DOT_PRODUCT(left(6, :, j, k), right(6, :, j, k)) &
                  + DOT_PRODUCT(left(7, :, j, k), right(7, :, j, k)) &
                  + DOT_PRODUCT(left(8, :, j, k), right(8, :, j, k)) &
                  + DOT_PRODUCT(left(9, :, j, k), right(9, :, j, k)) &
                  + DOT_PRODUCT(left(10, :, j, k), right(10, :, j, k))

      END DO
   END DO
   !END DO ! <--- Added a skeleton of a more efficient implementation.

END FUNCTION scprod

! 3.13. A subroutine for calculating the error between approximate ("xx" array) and analytical ("bmid" array) solution.
FUNCTION solerr(xx, bmid)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: xx, bmid ! Input arrays.

   REAL(rk) :: solerr ! Return variable.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z,
   INTEGER :: l       ! Component variable of 10D vector.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Array initialization.
   solerr = 0.0D0

   !===========================!
   ! Accumulative calculation. !
   !===========================!

   ! Quadruple looped summation.
   DO i = 1, nptx
      DO j = 1, npty
         DO k = 1, nptz
            DO l = 1, 10

               ! Discrete sum of absolute values of the differences between approximate and analytical solution.
               solerr = solerr + ABS((xx(l, i, j, k) - bmid(l, i, j, k)))**2
               !solerr = solerr + ABS((xx(l, i, j, k) - bmid(l, i, j, k))**2) ! <--- BUG.

            END DO
         END DO
      END DO
   END DO

END FUNCTION solerr

! New code:
! Note: The following modules and subroutines are responsible for computing positive operator action of C for a single iteration and diagnostics.
! **TODO: PRINT TENSORS AND VERIFY {Z.D}, {dag(Z).D} for tests

! 4.1. A subroutine for multiplying a 3x3D tensor "tensor" array to a 3D vector "vec" array.
SUBROUTINE multen(res, tensor, vec)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT)   :: res    ! Output array,
   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN) :: tensor ! Input array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)    :: vec    ! Input array.

   INTEGER :: i, j, k, l ! Iterator variables for x, y, and z.
   INTEGER :: ic, jc, kc ! Central point indices for debug output.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   ! Central point of the domain (integer division truncates, so this
   ! picks the nearest point at/after the true center for even sizes).
   ic = nptx/2 + 1
   jc = npty/2 + 1
   kc = nptz/2 + 1

   !===================!
   ! Array assignment. !
   !===================!

   ! Multiplying a vector by a tensor.
   DO i = 1, nptx
      DO j = 1, npty
         DO k = 1, nptz

            ! Assigning the vector.
            res(:, i, j, k) = MATMUL(tensor(:, :, i, j, k), vec(:, i, j, k))

            ! Debug: only print at the single central point of the domain.
            !IF (i == ic .AND. j == jc .AND. k == kc) THEN
            !   DO l = 1, 3
            !      PRINT *, "[DBG] tensor at (", i, ", ", j, ", ", k, ") = ", tensor(l, :, i, j, k)
            !   END DO
            !   PRINT *, "[DBG] vector at (", i, ", ", j, ", ", k, ") = ", vec(:, i, j, k)
            !   PRINT *, "[DBG] tensor*vector at (", i, ", ", j, ", ", k, ") = ", res(:, i, j, k)
            !END IF

         END DO
      END DO
   END DO

END SUBROUTINE multen

! 4.2. A subroutine for multiplying a 3x3D dagger tensor from "tensor" array to a 3D vector "vec" array.
SUBROUTINE muldagten(res, tensor, vec)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(OUT)   :: res    ! Output array,
   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN) :: tensor ! Input array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)    :: vec    ! Input array.

   COMPLEX(rk), DIMENSION(3, 3) :: dagger ! Temporary dagger-ed tensor.

   INTEGER :: i, j, k, l ! Iterator variables for x, y, and z.
   INTEGER :: ic, jc, kc ! Central point indices for debug output.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Array initialization.
   res = (0.0D0, 0.0D0)

   ! Central point of the domain.
   ic = nptx/2 + 1
   jc = npty/2 + 1
   kc = nptz/2 + 1

   !===================!
   ! Array assignment. !
   !===================!

   ! Multiplying a vector by a tensor.
   DO i = 1, nptx
      DO j = 1, npty
         DO k = 1, nptz

            ! Calculating dagger of a tensor.
            dagger = CONJG(TRANSPOSE(tensor(:, :, i, j, k)))

            ! Assigning the vector.
            res(:, i, j, k) = MATMUL(dagger, vec(:, i, j, k))

            ! Debug: only print at the single central point of the domain.
            !IF (i == ic .AND. j == jc .AND. k == kc) THEN
            !   DO l = 1, 3
            !      PRINT *, "[DBG] daggered tensor at (", i, ", ", j, ", ", k, ") = ", dagger(l, :)
            !   END DO
            !   PRINT *, "[DBG] vector at (", i, ", ", j, ", ", k, ") = ", vec(:, i, j, k)
            !   PRINT *, "[DBG] (daggered tensor)*vector at (", i, ", ", j, ", ", k, ") = ", res(:, i, j, k)
            !END IF

         END DO
      END DO
   END DO

END SUBROUTINE muldagten

! 4.3. A subroutine for calculating the operator action C(x), where x = (Dx, Dy, Dz, Phi, Hx, Hy, Hz, Ax, Ay, Az) is the unknown vector.
SUBROUTINE muladd(zeta, arg, res, tempr)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN)  :: zeta  ! Input array,
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN)    :: arg   ! Input array,
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(INOUT) :: res   ! Input/Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(INOUT)  :: tempr ! Temporary array.

   COMPLEX(rk) :: ik0 ! Wavenumber of the free space multiplied by an imaginary unit.
   REAL(rk)    :: k02 ! Wavenumber of the free space squared.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Wavevector.
   ik0 = (0.0D0, 1.0D0)*(pi/size0)
   k02 = (pi/size0)**2

   !================================================================================!
   ! A set of subroutines to calculate corrections to the governing equation for D. !
   !================================================================================!

   CALL multen(tempr, zeta, arg(1:3, :, :, :))                     ! tempr     = Z.D
   CALL muldagten(res(5:7, :, :, :), zeta, tempr)                  ! C(x){5:7} = dagger(Z).Z.D [temporary storage]
   res(1:3, :, :, :) = k02*(res(5:7, :, :, :) - arg(1:3, :, :, :)) ! C(x){1:3} = k0^2.(dagger(Z).Z.D - D)

   ! A halving adjustment of non-zero mass Dx-like terms at boundary nodes for x = 0 and x = Lx.
   res(1, 1, :, :) = res(1, 1, :, :) - 0.5D0*k02*(res(5, 1, :, :) - arg(1, 1, :, :))
   res(1, nptx, :, :) = res(1, nptx, :, :) - 0.5D0*k02*(res(5, nptx, :, :) - arg(1, nptx, :, :))

   ! A halving adjustment of non-zero mass Dy-like terms at boundary nodes for y = 0 and y = Ly.
   res(2, :, 1, :) = res(2, :, 1, :) - 0.5D0*k02*(res(6, :, 1, :) - arg(2, :, 1, :))
   res(2, :, npty, :) = res(2, :, npty, :) - 0.5D0*k02*(res(6, :, npty, :) - arg(2, :, npty, :))

   ! A halving adjustment of non-zero mass Dz-like terms at boundary nodes for z = 0 and z = Lz.
   res(3, :, :, 1) = res(3, :, :, 1) - 0.5D0*k02*(res(7, :, :, 1) - arg(3, :, :, 1))
   res(3, :, :, nptz) = res(3, :, :, nptz) - 0.5D0*k02*(res(7, :, :, nptz) - arg(3, :, :, nptz))

   CALL muldagten(tempr, zeta, arg(8:10, :, :, :))                          ! tempr     = dagger(Z).A
   res(1:3, :, :, :) = res(1:3, :, :, :) + k02*(arg(8:10, :, :, :) - tempr) ! C(x){1:3} = k0^2.(dagger(Z).Z.D - D) + k0^2.(A - dagger(Z).A)

   ! A halving adjustment of non-zero mass Ax-like terms at boundary nodes for x = 0 and x = Lx.
   res(1, 1, :, :) = res(1, 1, :, :) - 0.5D0*k02*(arg(8, 1, :, :) - tempr(1, 1, :, :))
   res(1, nptx, :, :) = res(1, nptx, :, :) - 0.5D0*k02*(arg(8, nptx, :, :) - tempr(1, nptx, :, :))

   ! A halving adjustment of non-zero mass Ay-like terms at boundary nodes for y = 0 and y = Ly.
   res(2, :, 1, :) = res(2, :, 1, :) - 0.5D0*k02*(arg(9, :, 1, :) - tempr(2, :, 1, :))
   res(2, :, npty, :) = res(2, :, npty, :) - 0.5D0*k02*(arg(9, :, npty, :) - tempr(2, :, npty, :))

   ! A halving adjustment of non-zero mass Az-like terms at boundary nodes for z = 0 and z = Lz.
   res(3, :, :, 1) = res(3, :, :, 1) - 0.5D0*k02*(arg(10, :, :, 1) - tempr(3, :, :, 1))
   res(3, :, :, nptz) = res(3, :, :, nptz) - 0.5D0*k02*(arg(10, :, :, nptz) - tempr(3, :, :, nptz))

   CALL nablaphi(res(5:7, :, :, :), arg(4, :, :, :))                       ! C(x){5:7} = grad(Phi) [temporary storage]
   CALL muldagten(tempr, zeta, res(5:7, :, :, :))                          ! tempr     = dagger(Z).grad(Phi)
   res(1:3, :, :, :) = res(1:3, :, :, :) + ik0*(res(5:7, :, :, :) - tempr) ! C(x){1:3} = k0^2.(A - dagger(Z).A) + k0^2.(A - dagger(Z).A) + i.k0.(grad(Phi) - dagger(Z).grad(Phi))

   !==================================================================================!
   ! A set of subroutines to calculate corrections to the governing equation for Phi. !
   !==================================================================================!

   CALL multen(tempr, zeta, arg(1:3, :, :, :))   ! tempr     = Z.D
   res(5:7, :, :, :) = arg(1:3, :, :, :) - tempr ! C(x){5:7} = D - Z.D [temporary storage]
   CALL divd(res(4, :, :, :), res(5:7, :, :, :)) ! C(x){4}   = div(D - Z.D)
   res(4, :, :, :) = ik0*res(4, :, :, :)         ! C(x){4}   = i.k0.div(D - Z.D)

   !================================================================================!
   ! A set of subroutines to calculate corrections to the governing equation for A. !
   !================================================================================!

   res(8:10, :, :, :) = k02*res(5:7, :, :, :) ! C(x){8:10} = k0^2.(D - Z.D)

   ! A halving adjustment of non-zero mass Dx-like terms at boundary nodes for x = 0 and x = Lx.
   res(8, 1, :, :) = res(8, 1, :, :) - 0.5D0*k02*res(5, 1, :, :)
   res(8, nptx, :, :) = res(8, nptx, :, :) - 0.5D0*k02*res(5, nptx, :, :)

   ! A halving adjustment of non-zero mass Dy-like terms at boundary nodes for y = 0 and y = Ly.
   res(9, :, 1, :) = res(9, :, 1, :) - 0.5D0*k02*res(6, :, 1, :)
   res(9, :, npty, :) = res(9, :, npty, :) - 0.5D0*k02*res(6, :, npty, :)

   ! A halving adjustment of non-zero mass Dz-like terms at boundary nodes for z = 0 and z = Lz.
   res(10, :, :, 1) = res(10, :, :, 1) - 0.5D0*k02*res(7, :, :, 1)
   res(10, :, :, nptz) = res(10, :, :, nptz) - 0.5D0*k02*res(7, :, :, nptz)

   !================================================================================!
   ! A set of subroutines to calculate corrections to the governing equation for H. !
   !================================================================================!

   res(5:7, :, :, :) = (0.0D0, 0.0D0) ! C(x){5:7} = 0

END SUBROUTINE muladd

! A SUBROUTINE TO CALCULATE C(X)
SUBROUTINE mul_l_correction(res, tmp, arg, zeta)

   USE numberformat
   USE indata
   USE constants

   IMPLICIT NONE

   ! Interface variables.
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz),   INTENT(INOUT) :: res
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz),    INTENT(INOUT) :: tmp
   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz),   INTENT(IN)    :: arg
   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN)    :: zeta

   ! Local variables.
   COMPLEX(rk) :: ik0
   REAL(rk)    :: k02

   ! Assigning runtime variables.
   ik0 = (0.0D0, 1.0D0)*pi/size0
   k02 = ABS(ik0)**2

   ! Initial assignment.
   res = (0.0D0, 0.0D0)

   ! Correction to governing equation for D field. (part 1)
   CALL multen(res(5:7, :, :, :), zeta, arg(1:3, :, :, :)) ! res(5:7) = Z*D
   CALL muldagten(tmp, zeta, res(5:7, :, :, :))            ! tmp = dagger(Z)*Z*D
   res(1:3, :, :, :) = k02*(tmp - arg(1:3, :, :, :))       ! res(1:3) = k0^2*(dagger(Z)*Z*D - D)

   ! Correction of mass terms near boundaries.
   res(1, 1, :, :)    = res(1, 1, :, :)    - 0.5D0*k02*(tmp(1, 1, :, :)    - arg(1, 1, :, :))
   res(1, nptx, :, :) = res(1, nptx, :, :) - 0.5D0*k02*(tmp(1, nptx, :, :) - arg(1, nptx, :, :))
   res(2, :, 1, :)    = res(2, :, 1, :)    - 0.5D0*k02*(tmp(2, :, 1, :)    - arg(2, :, 1, :))
   res(2, :, npty, :) = res(2, :, npty, :) - 0.5D0*k02*(tmp(2, :, npty, :) - arg(2, :, npty, :))
   res(3, :, :, 1)    = res(3, :, :, 1)    - 0.5D0*k02*(tmp(3, :, :, 1)    - arg(3, :, :, 1))
   res(3, :, :, nptz) = res(3, :, :, nptz) - 0.5D0*k02*(tmp(3, :, :, nptz) - arg(3, :, :, nptz))

   ! Correction to governing equation for D field. (part 2)
   CALL muldagten(tmp, zeta, arg(8:10, :, :, :))                          ! tmp = dagger(Z)*A
   res(1:3, :, :, :) = res(1:3, :, :, :) + k02*(arg(8:10, :, :, :) - tmp) ! res(1:3) = k0^2*(dagger(Z)*Z*D - D) + k0^2*(A - dagger(Z)*A)

   ! Correction of mass terms near boundaries.
   res(1, 1, :, :)    = res(1, 1, :, :)    - 0.5D0*k02*(arg(8, 1, :, :)     - tmp(1, 1, :, :))
   res(1, nptx, :, :) = res(1, nptx, :, :) - 0.5D0*k02*(arg(8, nptx, :, :)  - tmp(1, nptx, :, :))
   res(2, :, 1, :)    = res(2, :, 1, :)    - 0.5D0*k02*(arg(9, :, 1, :)     - tmp(2, :, 1, :))
   res(2, :, npty, :) = res(2, :, npty, :) - 0.5D0*k02*(arg(9, :, npty, :)  - tmp(2, :, npty, :))
   res(3, :, :, 1)    = res(3, :, :, 1)    - 0.5D0*k02*(arg(10, :, :, 1)    - tmp(3, :, :, 1))
   res(3, :, :, nptz) = res(3, :, :, nptz) - 0.5D0*k02*(arg(10, :, :, nptz) - tmp(3, :, :, nptz))

   ! Correction to governing equation for D field. (part 2)
   CALL nablaphi(res(5:7, :, :, :), arg(4, :, :, :))                     ! res(5:7) = grad(Phi)
   CALL muldagten(tmp, zeta, res(5:7, :, :, :))                          ! tmp = dagger(Z)*grad(Phi)
   res(1:3, :, :, :) = res(1:3, :, :, :) + ik0*(res(5:7, :, :, :) - tmp) ! res(1:3) = k0^2*(dagger(Z)*Z*D - D) + k0^2*(A - dagger(Z)*A) + ik0*(grad(Phi) - dagger(Z)*grad(Phi))

   ! Correction to governing equation for Phi field.
   CALL multen(res(5:7, :, :, :), zeta, arg(1:3, :, :, :))   ! res(5:7) = Z*D
   CALL divd(tmp(1, :, :, :), res(5:7, :, :, :))             ! tmp(1) = div(Z*D)
   CALL divd(tmp(2, :, :, :), arg(1:3, :, :, :))             ! tmp(2) = div(D)
   res(4, :, :, :) = ik0*(tmp(2, :, :, :) - tmp(1, :, :, :)) ! res(4) = ik0*div(D) - ik0*div(Z*D)

   ! Correction to governing equation for H field.
   res(5:7, :, :, :) = (0.0D0, 0.0D0)

   ! Correction to governing equation for A field.
   CALL multen(tmp, zeta, arg(1:3, :, :, :))          ! tmp = Z•D
   res(8:10, :, :, :) = k02*(arg(1:3, :, :, :) - tmp) ! res(8:10) = k₀²(D - Z•D)

   ! Correction of mass terms near boundaries.
   res(8, 1, :, :)     = res(8, 1, :, :)     - 0.5D0*k02*(arg(1, 1, :, :)    - tmp(1, 1, :, :))
   res(8, nptx, :, :)  = res(8, nptx, :, :)  - 0.5D0*k02*(arg(1, nptx, :, :) - tmp(1, nptx, :, :))
   res(9, :, 1, :)     = res(9, :, 1, :)     - 0.5D0*k02*(arg(2, :, 1, :)    - tmp(2, :, 1, :))
   res(9, :, npty, :)  = res(9, :, npty, :)  - 0.5D0*k02*(arg(2, :, npty, :) - tmp(2, :, npty, :))
   res(10, :, :, 1)    = res(10, :, :, 1)    - 0.5D0*k02*(arg(3, :, :, 1)    - tmp(3, :, :, 1))
   res(10, :, :, nptz) = res(10, :, :, nptz) - 0.5D0*k02*(arg(3, :, :, nptz) - tmp(3, :, :, nptz))

END SUBROUTINE mul_l_correction

! 4.4. A subroutine for calculating the right-hand-side of positive definite equations using just "Jext" array.
! Note: Applies boundary conditions for Jext: Normal components have Neumann boundary conditions, tangent components have Dirichlet boundary conditions.
SUBROUTINE setrhs(rhs, jext, tempr)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(OUT)  :: rhs   ! Output array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)    :: jext  ! Input array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(INOUT) :: tempr ! Temporary array.

   COMPLEX(rk) :: ik0 ! Wavenumber of the free space multiplied by an imaginary unit.
   REAL(rk)    :: k02 ! Wavenumber of the free space squared.
   REAL(rk)    :: cns ! A constant factor: cns = 4.pi/c

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Wavevector.
   ik0 = (0.0D0, 1.0D0)*(pi/size0)
   k02 = (pi/size0)**2

   ! Constant factor.
   cns = 4.0D0*pi/lightspeed

   !========================================================================!
   ! A set of subroutines to calculate rhs to the governing equation for D. !
   !========================================================================!

   CALL graddiv(tempr, jext)                        ! tempr        = grad(div[Jext]),
   rhs(1:3, :, :, :) = -cns*(ik0*jext + tempr/ik0) ! R(Jext){1:3} = - 4.pi.(i.k0.Jext + grad[div{Jext}]/[i.k0])/c.

   ! A halving adjustment of non-zero mass term Jextx at boundary nodes for x = 0 and x = Lx.
   rhs(1, 1, :, :) = rhs(1, 1, :, :) + 0.5D0*cns*ik0*jext(1, 1, :, :)
   rhs(1, nptx, :, :) = rhs(1, nptx, :, :) + 0.5D0*cns*ik0*jext(1, nptx, :, :)

   ! A halving adjustment of non-zero mass term Jexty at boundary nodes for y = 0 and y = Ly.
   rhs(2, :, 1, :) = rhs(2, :, 1, :) + 0.5D0*cns*ik0*jext(2, :, 1, :)
   rhs(2, :, npty, :) = rhs(2, :, npty, :) + 0.5D0*cns*ik0*jext(2, :, npty, :)

   ! A halving adjustment of non-zero mass term Jextz at boundary nodes for z = 0 and z = Lz.
   rhs(3, :, :, 1) = rhs(3, :, :, 1) + 0.5D0*cns*ik0*jext(3, :, :, 1)
   rhs(3, :, :, nptz) = rhs(3, :, :, nptz) + 0.5D0*cns*ik0*jext(3, :, :, nptz)

   !==========================================================================!
   ! A set of subroutines to calculate rhs to the governing equation for Phi. !
   !==========================================================================!

   rhs(4, :, :, :) = (0.0D0, 0.0D0) ! R(Jext){4} = 0.

   !========================================================================!
   ! A set of subroutines to calculate rhs to the governing equation for H. !
   !========================================================================!

   CALL rotd(tempr, jext)        ! tempr = curl(Jext),
   rhs(5:7, :, :, :) = cns*tempr ! R(Jext){5:7} = 4.pi.curl(Jext)/c.

   !========================================================================!
   ! A set of subroutines to calculate rhs to the governing equation for A. !
   !========================================================================!

   rhs(8:10, :, :, :) = (0.0D0, 0.0D0) ! R(Jext){8:10} = 0.

END SUBROUTINE setrhs

! 5.1. A function for computing a normalized residual of div(D) = rhoext from approximate solution in "xx" array and "Jext" array.
FUNCTION resdivd(xx, Jext)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: xx   ! Input array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: Jext ! Input array.

   REAL(rk) :: resdivd ! Return value.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: ihx, ihy, ihz ! Recipocals of spacings between nodes.
   REAL(rk) :: wx, wy, wz ! Quadrature.

   COMPLEX(rk), DIMENSION(3) :: lhst ! Individual terms of div(D) (i.e. dDx/dx, dDy/dy, and dDz/dz),
   COMPLEX(rk)               :: lhs  ! div(D) value,
   COMPLEX(rk)               :: rhs  ! rhoext value.

   REAL(rk) :: num, den ! Numerator and denominator.

   LOGICAL, PARAMETER :: docomputemax = .TRUE.              ! Debug flag.
   REAL(rk)           :: candidate                          ! Current residual.
   REAL(rk)           :: maxresdivdi, maxresdivdo           ! Maximum value of the residual of resdivd
   INTEGER            :: maxii = -1, maxij = -1, maxik = -1 ! Interior node.
   INTEGER            :: maxoi = -1, maxoj = -1, maxok = -1 ! Outer node.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Accumulation initialization.
   resdivd = 0.0D0
   num = 0.0D0
   den = 0.0D0

   ! Calculation of recipocals of spacings.
   ihx = (nptx - 1.0D0)/sizex
   ihy = (npty - 1.0D0)/sizey
   ihz = (nptz - 1.0D0)/sizez

   !===========================!
   ! Accumulative calculation. !
   !===========================!

   maxresdivdi = 0.0D0
   maxresdivdo = 0.0D0

   ! Tripple looped summation.
   DO i = 1, nptx

      wx = 1.0D0
      IF (i == 1 .OR. i == nptx) wx = 0.5D0

      DO j = 1, npty

         wy = 1.0D0
         IF (j == 1 .OR. j == npty) wy = 0.5D0

         DO k = 1, nptz

            wz = 1.0D0
            IF (k == 1 .OR. k == nptz) wz = 0.5D0

            ! Initialization for incremental accumulation.
            lhst = (0.0D0, 0.0D0)
            lhs = (0.0D0, 0.0D0)
            rhs = (0.0D0, 0.0D0)

            !=================================!
            ! Computation of individual term. !
            !=================================!

            ! dDx/dx and dJextx/dx.
            IF (i == 1) THEN
               lhst(1) = (xx(1, i + 1, j, k) - xx(1, i, j, k))*ihx
               rhs = rhs + 4.0D0*pi*(jext(1, i + 1, j, k) - jext(1, i, j, k))*ihx
            ELSEIF (i == nptx) THEN
               lhst(1) = (xx(1, i, j, k) - xx(1, i - 1, j, k))*ihx
               rhs = rhs + 4.0D0*pi*(jext(1, i, j, k) - jext(1, i - 1, j, k))*ihx
            ELSE
               lhst(1) = 0.5D0*(xx(1, i + 1, j, k) - xx(1, i - 1, j, k))*ihx
               rhs = rhs + 2.0D0*pi*(jext(1, i + 1, j, k) - jext(1, i - 1, j, k))*ihx
            END IF

            ! dDy/dy and dJexty/dy.
            IF (j == 1) THEN
               lhst(2) = (xx(2, i, j + 1, k) - xx(2, i, j, k))*ihy
               rhs = rhs + 4.0D0*pi*(jext(2, i, j + 1, k) - jext(2, i, j, k))*ihy
            ELSEIF (j == npty) THEN
               lhst(2) = (xx(2, i, j, k) - xx(2, i, j - 1, k))*ihy
               rhs = rhs + 4.0D0*pi*(jext(2, i, j, k) - jext(2, i, j - 1, k))*ihy
            ELSE
               lhst(2) = 0.5D0*(xx(2, i, j + 1, k) - xx(2, i, j - 1, k))*ihy
               rhs = rhs + 2.0D0*pi*(jext(2, i, j + 1, k) - jext(2, i, j - 1, k))*ihy
            END IF

            ! dDz/dz and dJextz/dz.
            IF (k == 1) THEN
               lhst(3) = (xx(3, i, j, k + 1) - xx(3, i, j, k))*ihz
               rhs = rhs + 4.0D0*pi*(jext(3, i, j, k + 1) - jext(3, i, j, k))*ihz
            ELSEIF (k == nptz) THEN
               lhst(3) = (xx(3, i, j, k) - xx(3, i, j, k - 1))*ihz
               rhs = rhs + 4.0D0*pi*(jext(3, i, j, k) - jext(3, i, j, k - 1))*ihz
            ELSE
               lhst(3) = 0.5D0*(xx(3, i, j, k + 1) - xx(3, i, j, k - 1))*ihz
               rhs = rhs + 2.0D0*pi*(jext(3, i, j, k + 1) - jext(3, i, j, k - 1))*ihz
            END IF

            ! div(D) = dDx/dx + dDy/dy + dDz/dz
            lhs = SUM(lhst)

            !============================================!
            ! Calculation of fraction under the radical. !
            !============================================!

            num = num + wx*wy*wz*ABS(lhs - rhs)**2
            den = den + wx*wy*wz*SUM(ABS(lhst)**2)

            !==============================================!
            ! Calculation of the residual at current node. !
            !==============================================!

            IF (docomputemax) THEN

               candidate = SQRT(wx*wy*wz*(ABS(lhs - rhs)**2))

               IF (i > 1 .AND. i < nptx .AND. j > 1 .AND. j < npty .AND. k > 1 .AND. k < nptz) THEN

                  IF (maxresdivdi < candidate) THEN

                     maxresdivdi = candidate
                     maxii = i
                     maxij = j
                     maxik = k

                  END IF

               ELSE

                  IF (maxresdivdo < candidate) THEN

                     maxresdivdo = candidate
                     maxoi = i
                     maxoj = j
                     maxok = k

                  END IF

               END IF

            END IF

         END DO
      END DO
   END DO

  IF (docomputemax) PRINT *, "[DBG] Inner maxresdivd = ", maxresdivdi*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxii, ", ", maxij, ", ", maxik, "."
  IF (docomputemax) PRINT *, "[DBG] Outer maxresdivd = ", maxresdivdo*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxoi, ", ", maxoj, ", ", maxok, "."

   ! Taking a square root.
   resdivd = SQRT(num/den)

END FUNCTION resdivd

! 5.2. A function for computing a normalized residual of div(H) = 0 from approximate solution in "xx" array.
FUNCTION resdivh(xx)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: xx ! Input array.

   REAL(rk) :: resdivh ! Return value5

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: ihx, ihy, ihz ! Recipocals of spacings between nodes.
   REAL(rk) :: wx, wy, wz ! Quadrature.

   COMPLEX(rk), DIMENSION(3) :: lhst ! Individual terms of div(H) (i.e. dHx/dx, dHy/dy, and dHz/dz),
   COMPLEX(rk)               :: lhs  ! div(H) value.

   REAL(rk) :: num, den ! Numerator and denominator.

   LOGICAL, PARAMETER :: docomputemax = .TRUE.              ! Debug flag.
   REAL(rk)           :: candidate                          ! Current residual.
   REAL(rk)           :: maxresdivhi, maxresdivho           ! Maximum value of the residual of resdivd
   INTEGER            :: maxii = -1, maxij = -1, maxik = -1 ! Inner node.
   INTEGER            :: maxoi = -1, maxoj = -1, maxok = -1 ! Outer node.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Accumulation initialization.
   resdivh = 0.0D0
   num = 0.0D0
   den = 0.0D0

   ! Calculation of recipocals of spacings.
   ihx = (nptx - 1.0D0)/sizex
   ihy = (npty - 1.0D0)/sizey
   ihz = (nptz - 1.0D0)/sizez

   !===========================!
   ! Accumulative calculation. !
   !===========================!

   maxresdivhi = 0.0D0
   maxresdivho = 0.0D0

   ! Tripple looped summation.
   DO i = 1, nptx

      wx = 1.0D0
      IF (i == 1 .OR. i == nptx) wx = 0.5D0

      DO j = 1, npty

         wy = 1.0D0
         IF (j == 1 .OR. j == npty) wy = 0.5D0

         DO k = 1, nptz

            wz = 1.0D0
            IF (k == 1 .OR. k == nptz) wz = 0.5D0

            ! Initialization for incremental accumulation.
            lhst = (0.0D0, 0.0D0)
            lhs = (0.0D0, 0.0D0)

            !=================================!
            ! Computation of individual term. !
            !=================================!

            ! dHx/dx.
            IF (i == 1) THEN
               lhst(1) = (xx(5, i + 1, j, k) - xx(5, i, j, k))*ihx
            ELSEIF (i == nptx) THEN
               lhst(1) = (xx(5, i, j, k) - xx(5, i - 1, j, k))*ihx
            ELSE
               lhst(1) = 0.5D0*(xx(5, i + 1, j, k) - xx(5, i - 1, j, k))*ihx
            END IF

            ! dHy/dy.
            IF (j == 1) THEN
               lhst(2) = (xx(6, i, j + 1, k) - xx(6, i, j, k))*ihy
            ELSEIF (j == npty) THEN
               lhst(2) = (xx(6, i, j, k) - xx(6, i, j - 1, k))*ihy
            ELSE
               lhst(2) = 0.5D0*(xx(6, i, j + 1, k) - xx(6, i, j - 1, k))*ihy
            END IF

            ! dHz/dz.
            IF (k == 1) THEN
               lhst(3) = (xx(7, i, j, k + 1) - xx(7, i, j, k))*ihz
            ELSEIF (k == nptz) THEN
               lhst(3) = (xx(7, i, j, k) - xx(7, i, j, k - 1))*ihz
            ELSE
               lhst(3) = 0.5D0*(xx(7, i, j, k + 1) - xx(7, i, j, k - 1))*ihz
            END IF

            ! div(H) = dHx/dx + dHy/dy + dHz/dz
            lhs = SUM(lhst)

            !=======================================!
            ! Calculation of sum under the radical. !
            !=======================================!

            num = num + wx*wy*wz*ABS(lhs)**2
            den = den + wx*wy*wz*SUM(ABS(lhst)**2)

            !==============================================!
            ! Calculation of the residual at current node. !
            !==============================================!

            IF (docomputemax) THEN

               candidate = SQRT(wx*wy*wz*ABS(lhs)**2)
               IF (i > 1 .AND. i < nptx .AND. j > 1 .AND. j < npty .AND. k > 1 .AND. k < nptz) THEN

                  IF (maxresdivhi < candidate) THEN

                     maxresdivhi = candidate
                     maxii = i
                     maxij = j
                     maxik = k

                  END IF

               ELSE

                  IF (maxresdivho < candidate) THEN

                     maxresdivho = candidate
                     maxoi = i
                     maxoj = j
                     maxok = k

                  END IF

               END IF

            END IF

         END DO
      END DO
   END DO

  IF (docomputemax) PRINT *, "[DBG] Inner maxresdivh = ", maxresdivhi*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxii, ", ", maxij, ", ", maxik, "."
  IF (docomputemax) PRINT *, "[DBG] Outer maxresdivh = ", maxresdivho*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxoi, ", ", maxoj, ", ", maxok, "."

   ! Taking a square root.
   resdivh = SQRT(num/den)

END FUNCTION resdivh

! 5.3. A function for computing a normalized residual of div(A) = 0 from approximate solution in "xx" array.
FUNCTION resdiva(xx)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: xx ! Input array.

   REAL(rk) :: resdiva ! Return value

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: ihx, ihy, ihz ! Recipocals of spacings between nodes.
   REAL(rk) :: wx, wy, wz ! Quadrature.

   COMPLEX(rk), DIMENSION(3) :: lhst ! Individual terms of div(A) (i.e. dAx/dx, dAy/dy, and dAz/dz),
   COMPLEX(rk)               :: lhs  ! div(A) value.

   REAL(rk) :: num, den ! Numerator and denominator.

   LOGICAL, PARAMETER :: docomputemax = .TRUE.              ! Debug flag.
   REAL(rk)           :: candidate                          ! Current residual.
   REAL(rk)           :: maxresdivai, maxresdivao           ! Maximum value of the residual of resdivd
   INTEGER            :: maxii = -1, maxij = -1, maxik = -1 ! Inner node.
   INTEGER            :: maxoi = -1, maxoj = -1, maxok = -1 ! Outer node.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   maxresdivai = 0.0D0
   maxresdivao = 0.0D0

   ! Accumulation initialization.
   resdiva = 0.0D0
   num = 0.0D0
   den = 0.0D0

   ! Calculation of recipocals of spacings.
   ihx = (nptx - 1.0D0)/sizex
   ihy = (npty - 1.0D0)/sizey
   ihz = (nptz - 1.0D0)/sizez

   !===========================!
   ! Accumulative calculation. !
   !===========================!

   ! Tripple looped summation.
   DO i = 1, nptx

      wx = 1.0D0
      IF (i == 1 .OR. i == nptx) wx = 0.5D0

      DO j = 1, npty

         wy = 1.0D0
         IF (j == 1 .OR. j == npty) wy = 0.5D0

         DO k = 1, nptz

            wz = 1.0D0
            IF (k == 1 .OR. k == nptz) wz = 0.5D0

            ! Initialization for incremental accumulation.
            lhst = (0.0D0, 0.0D0)
            lhs = (0.0D0, 0.0D0)

            !=================================!
            ! Computation of individual term. !
            !=================================!

            ! dAx/dx.
            IF (i == 1) THEN
               lhst(1) = (xx(8, i + 1, j, k) - xx(8, i, j, k))*ihx
            ELSEIF (i == nptx) THEN
               lhst(1) = (xx(8, i, j, k) - xx(8, i - 1, j, k))*ihx
            ELSE
               lhst(1) = 0.5D0*(xx(8, i + 1, j, k) - xx(8, i - 1, j, k))*ihx
            END IF

            ! dAy/dy.
            IF (j == 1) THEN
               lhst(2) = (xx(9, i, j + 1, k) - xx(9, i, j, k))*ihy
            ELSEIF (j == npty) THEN
               lhst(2) = (xx(9, i, j, k) - xx(9, i, j - 1, k))*ihy
            ELSE
               lhst(2) = 0.5D0*(xx(9, i, j + 1, k) - xx(9, i, j - 1, k))*ihy
            END IF

            ! dAz/dz.
            IF (k == 1) THEN
               lhst(3) = (xx(10, i, j, k + 1) - xx(10, i, j, k))*ihz
            ELSEIF (k == nptz) THEN
               lhst(3) = (xx(10, i, j, k) - xx(10, i, j, k - 1))*ihz
            ELSE
               lhst(3) = 0.5D0*(xx(10, i, j, k + 1) - xx(10, i, j, k - 1))*ihz
            END IF

            ! div(A) = dAx/dx + dAy/dy + dAz/dz
            lhs = SUM(lhst)

            !=======================================!
            ! Calculation of sum under the radical. !
            !=======================================!

            num = num + wx*wy*wz*ABS(lhs)**2
            den = den + wx*wy*wz*SUM(ABS(lhst)**2)

            !==============================================!
            ! Calculation of the residual at current node. !
            !==============================================!

            IF (docomputemax) THEN

               candidate = SQRT(wx*wy*wz*ABS(lhs)**2)
               IF (i > 1 .AND. i < nptx .AND. j > 1 .AND. j < npty .AND. k > 1 .AND. k < nptz) THEN

                  IF (maxresdivai < candidate) THEN

                     maxresdivai = candidate
                     maxii = i
                     maxij = j
                     maxik = k

                  END IF

               ELSE

                  IF (maxresdivao < candidate) THEN

                     maxresdivao = candidate
                     maxoi = i
                     maxoj = j
                     maxok = k

                  END IF

               END IF

            END IF

         END DO
      END DO
   END DO

  IF (docomputemax) PRINT *, "[DBG] Inner maxresdiva = ", maxresdivai*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxii, ", ", maxij, ", ", maxik, "."
  IF (docomputemax) PRINT *, "[DBG] Outer maxresdiva = ", maxresdivao*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxoi, ", ", maxoj, ", ", maxok, "."

   ! Taking a square root.
   resdiva = SQRT(num/den)

END FUNCTION resdiva

! 5.4. A function for computing a normalized residual of curl(A) = i.k0.H and curl(H) = 4.pi.Jext/c - i.k0.D from approximate solution in "xx" array and "Jext" array.
FUNCTION rescurlsrc(xx, Jext)

   USE numberformat
   USE constants
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: xx   ! Input array,
   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN)  :: Jext ! Input array.

   REAL(rk) :: rescurlsrc ! Return value.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: ihx, ihy, ihz ! Recipocals of spacings between nodes.
   REAL(rk) :: wx, wy, wz ! Quadrature.

   COMPLEX(rk), DIMENSION(12) :: pdtm ! 6 individual values for temporary calculation of partial derivatives.
   COMPLEX(rk), DIMENSION(6)  :: lhst ! 6 individual values, 3 for curl(A) and 3 for curl(H) (both of the form dFz/dy - dFy/dz, dFx/dz - dFz/dx, and dFy/dx - dFx/dy, where F can be H or A),
   COMPLEX(rk), DIMENSION(6)  :: rhst ! 6 individual values, 3 for i.k0.H and 3 individual for 4.pi.Jext/c - i.k0.D.

   COMPLEX(rk) :: ik0 ! Wavenumber of the free space multiplied by an imaginary unit,
   REAL(rk)    :: cns ! A constant factor: cns = 4.pi/c

   REAL(rk) :: num, den ! Numerator and denominator.

   LOGICAL, PARAMETER :: docomputemax = .TRUE.              ! Debug flag.
   REAL(rk)           :: candidate                          ! Current residual.
   REAL(rk)           :: maxrescurlsrci, maxrescurlsrco     ! Maximum value of the residual of resdivd
   INTEGER            :: maxii = -1, maxij = -1, maxik = -1 ! Inner node.
   INTEGER            :: maxoi = -1, maxoj = -1, maxok = -1 ! Outer node.

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Wavevector.
   ik0 = (0.0D0, 1.0D0)*(pi/size0)

   ! Constant factor.
   cns = 4.0D0*pi/lightspeed

   !===========================!
   ! Calculation of constants. !
   !===========================!

   ! Accumulation initialization.
   rescurlsrc = 0.0D0
   num = 0.0D0
   den = 0.0D0

   ! Calculation of recipocals of spacings.
   ihx = (nptx - 1.0D0)/sizex
   ihy = (npty - 1.0D0)/sizey
   ihz = (nptz - 1.0D0)/sizez

   !===========================!
   ! Accumulative calculation. !
   !===========================!

   maxrescurlsrci = 0.0D0
   maxrescurlsrco = 0.0D0

   ! Tripple looped summation.
   DO i = 1, nptx

      wx = 1.0D0
      IF (i == 1 .OR. i == nptx) wx = 0.5D0

      DO j = 1, npty

         wy = 1.0D0
         IF (j == 1 .OR. j == npty) wy = 0.5D0

         DO k = 1, nptz

            wz = 1.0D0
            IF (k == 1 .OR. k == nptz) wz = 0.5D0

            !==================================!
            ! Computation of individual terms. !
            !==================================!

            ! Partial derivatives with respect to x.
            IF (i == 1) THEN
               pdtm(1) = (xx(6, i + 1, j, k) - xx(6, i, j, k))*ihx ! dHy/dx
               pdtm(2) = (xx(7, i + 1, j, k) - xx(7, i, j, k))*ihx ! dHz/dx
               pdtm(3) = (xx(9, i + 1, j, k) - xx(9, i, j, k))*ihx ! dAy/dx
               pdtm(4) = (xx(10, i + 1, j, k) - xx(10, i, j, k))*ihx ! dAz/dx
            ELSEIF (i == nptx) THEN
               pdtm(1) = (xx(6, i, j, k) - xx(6, i - 1, j, k))*ihx ! dHy/dx
               pdtm(2) = (xx(7, i, j, k) - xx(7, i - 1, j, k))*ihx ! dHz/dx
               pdtm(3) = (xx(9, i, j, k) - xx(9, i - 1, j, k))*ihx ! dAy/dx
               pdtm(4) = (xx(10, i, j, k) - xx(10, i - 1, j, k))*ihx ! dAz/dx
            ELSE
               pdtm(1) = 0.5D0*(xx(6, i + 1, j, k) - xx(6, i - 1, j, k))*ihx ! dHy/dx
               pdtm(2) = 0.5D0*(xx(7, i + 1, j, k) - xx(7, i - 1, j, k))*ihx ! dHz/dx
               pdtm(3) = 0.5D0*(xx(9, i + 1, j, k) - xx(9, i - 1, j, k))*ihx ! dAy/dx
               pdtm(4) = 0.5D0*(xx(10, i + 1, j, k) - xx(10, i - 1, j, k))*ihx ! dAz/dx
            END IF

            ! Partial derivatives with respect to y.
            IF (j == 1) THEN
               pdtm(5) = (xx(5, i, j + 1, k) - xx(5, i, j, k))*ihy ! dHx/dy
               pdtm(6) = (xx(7, i, j + 1, k) - xx(7, i, j, k))*ihy ! dHz/dy
               pdtm(7) = (xx(8, i, j + 1, k) - xx(8, i, j, k))*ihy ! dAx/dy
               pdtm(8) = (xx(10, i, j + 1, k) - xx(10, i, j, k))*ihy ! dAz/dy
            ELSEIF (j == npty) THEN
               pdtm(5) = (xx(5, i, j, k) - xx(5, i, j - 1, k))*ihy ! dHx/dy
               pdtm(6) = (xx(7, i, j, k) - xx(7, i, j - 1, k))*ihy ! dHz/dy
               pdtm(7) = (xx(8, i, j, k) - xx(8, i, j - 1, k))*ihy ! dAx/dy
               pdtm(8) = (xx(10, i, j, k) - xx(10, i, j - 1, k))*ihy ! dAz/dy
            ELSE
               pdtm(5) = 0.5D0*(xx(5, i, j + 1, k) - xx(5, i, j - 1, k))*ihy ! dHx/dy
               pdtm(6) = 0.5D0*(xx(7, i, j + 1, k) - xx(7, i, j - 1, k))*ihy ! dHz/dy
               pdtm(7) = 0.5D0*(xx(8, i, j + 1, k) - xx(8, i, j - 1, k))*ihy ! dAx/dy
               pdtm(8) = 0.5D0*(xx(10, i, j + 1, k) - xx(10, i, j - 1, k))*ihy ! dAz/dy
            END IF

            ! Partial derivatives with respect to z.
            IF (k == 1) THEN
               pdtm(9) = (xx(5, i, j, k + 1) - xx(5, i, j, k))*ihz ! dHx/dz
               pdtm(10) = (xx(6, i, j, k + 1) - xx(6, i, j, k))*ihz ! dHy/dz
               pdtm(11) = (xx(8, i, j, k + 1) - xx(8, i, j, k))*ihz ! dAx/dz
               pdtm(12) = (xx(9, i, j, k + 1) - xx(9, i, j, k))*ihz ! dAy/dz
            ELSEIF (k == nptz) THEN
               pdtm(9) = (xx(5, i, j, k) - xx(5, i, j, k - 1))*ihz ! dHx/dz
               pdtm(10) = (xx(6, i, j, k) - xx(6, i, j, k - 1))*ihz ! dHy/dz
               pdtm(11) = (xx(8, i, j, k) - xx(8, i, j, k - 1))*ihz ! dAx/dz
               pdtm(12) = (xx(9, i, j, k) - xx(9, i, j, k - 1))*ihz ! dAy/dz
            ELSE
               pdtm(9) = 0.5D0*(xx(5, i, j, k + 1) - xx(5, i, j, k - 1))*ihz ! dHx/dz
               pdtm(10) = 0.5D0*(xx(6, i, j, k + 1) - xx(6, i, j, k - 1))*ihz ! dHy/dz
               pdtm(11) = 0.5D0*(xx(8, i, j, k + 1) - xx(8, i, j, k - 1))*ihz ! dAx/dz
               pdtm(12) = 0.5D0*(xx(9, i, j, k + 1) - xx(9, i, j, k - 1))*ihz ! dAy/dz
            END IF

            !=======================================!
            ! Assigning left-, and right-hand-side. !
            !=======================================!

            ! Left-hand-side.
            lhst(1) = pdtm(8) - pdtm(12) ! dAz/dy - dAy/dz
            lhst(2) = pdtm(11) - pdtm(4) ! dAx/dz - dAz/dx
            lhst(3) = pdtm(3) - pdtm(7)  ! dAy/dx - dAx/dy
            lhst(4) = pdtm(6) - pdtm(10) ! dHz/dy - dHy/dz
            lhst(5) = pdtm(9) - pdtm(2)  ! dHx/dz - dHz/dx
            lhst(6) = pdtm(1) - pdtm(5)  ! dHy/dx - dHx/dy

            ! Right-hand-side.
            rhst(1:3) = ik0*xx(5:7, i, j, k)
            rhst(4:6) = cns*jext(1:3, i, j, k) - ik0*xx(1:3, i, j, k)

            !=======================================!
            ! Calculation of sum under the radical. !
            !=======================================!

            num = num + wx*wy*wz*SUM(ABS(lhst(1:3) - rhst(1:3))**2 + ABS(lhst(4:6) - rhst(4:6))**2)
            den = den + wx*wy*wz*SUM(ABS(rhst)**2)

            !==============================================!
            ! Calculation of the residual at current node. !
            !==============================================!

            IF (docomputemax) THEN

               candidate = SQRT(wx*wy*wz*SUM(ABS(lhst(1:3) - rhst(1:3))**2 + ABS(lhst(4:6) - rhst(4:6))**2))
               IF (i > 1 .AND. i < nptx .AND. j > 1 .AND. j < npty .AND. k > 1 .AND. k < nptz) THEN

                  IF (maxrescurlsrci < candidate) THEN

                     maxrescurlsrci = candidate
                     maxii = i
                     maxij = j
                     maxik = k

                  END IF

               ELSE

                  IF (maxrescurlsrco < candidate) THEN

                     maxrescurlsrco = candidate
                     maxoi = i
                     maxoj = j
                     maxok = k

                  END IF

               END IF

            END IF

         END DO
      END DO
   END DO

  IF (docomputemax) PRINT *, "[DBG] Inner maxrescurlsrc = ", maxrescurlsrci*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxii, ", ", maxij, ", ", maxik, "."
  IF (docomputemax) PRINT *, "[DBG] Outer maxrescurlsrc = ", maxrescurlsrco*SQRT(nptx*npty*nptz/den), ", at node i, j, k = ", maxoi, ", ", maxoj, ", ", maxok, "."

   ! Taking a square root.
   rescurlsrc = SQRT(num/den)

END FUNCTION rescurlsrc

! 5.5. A subroutine for writing output from a 3D vector field array.
SUBROUTINE write3dvector(vec, filename)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, nptx, npty, nptz), INTENT(IN) :: vec      ! Input array,
   CHARACTER(LEN=*)                                      :: filename ! Input string.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: argx, argy, argz ! Actual values of x, y, and z.

   OPEN (16, FILE=filename, STATUS='replace')

   WRITE (16, *) "X, Y, Z, RE(Vx), IM(Vx), RE(Vy), IM(Vy), RE(Vz), IM(Vz)" ! Writing the header.

   DO i = 1, nptx

      ! Computing x.
      argx = sizex*(i - 1.0D0)/(nptx - 1.0D0)

      DO j = 1, npty

         ! Computing y.
         argy = sizey*(j - 1.0D0)/(npty - 1.0D0)

         DO k = 1, nptz

            ! Computing z.
            argz = sizez*(k - 1.0D0)/(nptz - 1.0D0)

            WRITE (16, '(8(ES24.16,", "), ES24.16)') &
               argx, argy, argz, &
               REAL(vec(1, i, j, k)), AIMAG(vec(1, i, j, k)), &
               REAL(vec(2, i, j, k)), AIMAG(vec(2, i, j, k)), &
               REAL(vec(3, i, j, k)), AIMAG(vec(3, i, j, k))

         END DO
      END DO
   END DO

   CLOSE (16)

END SUBROUTINE write3dvector

! 5.5. A subroutine for writing output from a 3x3D tensor field array.
SUBROUTINE write3x3dtensor(tensor, filename)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(3, 3, nptx, npty, nptz), INTENT(IN) :: tensor   ! Input array,
   CHARACTER(LEN=*)                                         :: filename ! Input string.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: argx, argy, argz ! Actual values of x, y, and z.

   OPEN (17, FILE=filename, STATUS='replace')

  WRITE(17, *) "X, Y, Z, RE(Txx), IM(Txx), RE(Txy), IM(Txy), RE(Txz), IM(Txz), RE(Tyx), IM(Tyx), RE(Tyy), IM(Tyy), RE(Tyz), IM(Tyz), RE(Tzx), IM(Tzx), RE(Tzy), IM(Tzy), RE(Tzz), IM(Tzz)" ! Writing the header.

   DO i = 1, nptx

      ! Computing x.
      argx = sizex*(i - 1.0D0)/(nptx - 1.0D0)

      DO j = 1, npty

         ! Computing y.
         argy = sizey*(j - 1.0D0)/(npty - 1.0D0)

         DO k = 1, nptz

            ! Computing z.
            argz = sizez*(k - 1.0D0)/(nptz - 1.0D0)

            WRITE (17, '(20(ES24.16,", "), ES24.16)') &
               argx, argy, argz, &
          REAL(tensor(1, 1, i, j, k)), AIMAG(tensor(1, 1, i, j, k)), REAL(tensor(1, 2, i, j, k)), AIMAG(tensor(1, 2, i, j, k)), REAL(tensor(1, 3, i, j, k)), AIMAG(tensor(1, 3, i, j, k)), &
          REAL(tensor(2, 1, i, j, k)), AIMAG(tensor(2, 1, i, j, k)), REAL(tensor(2, 2, i, j, k)), AIMAG(tensor(2, 2, i, j, k)), REAL(tensor(2, 3, i, j, k)), AIMAG(tensor(2, 3, i, j, k)), &
          REAL(tensor(3, 1, i, j, k)), AIMAG(tensor(3, 1, i, j, k)), REAL(tensor(3, 2, i, j, k)), AIMAG(tensor(3, 2, i, j, k)), REAL(tensor(3, 3, i, j, k)), AIMAG(tensor(3, 3, i, j, k))

         END DO
      END DO
   END DO

   CLOSE (17)

END SUBROUTINE write3x3dtensor

! 5.5. A subroutine for writing output from a 10D vector field array.
SUBROUTINE write10dvector(xx, filename)

   USE numberformat
   USE indata

   IMPLICIT NONE

   COMPLEX(rk), DIMENSION(10, nptx, npty, nptz), INTENT(IN) :: xx       ! Input array,
   CHARACTER(LEN=*)                                       :: filename ! Input string.

   INTEGER :: i, j, k ! Iterator variables for x, y, and z.

   REAL(rk) :: argx, argy, argz ! Actual values of x, y, and z.

   OPEN (18, FILE=filename, STATUS='replace')

  WRITE(18, *) "X, Y, Z, RE(Dx), IM(Dx), RE(Dy), IM(Dy), RE(Dz), IM(Dz), RE(Phi), IM(Phi), RE(Hx), IM(Hx), RE(Hy), IM(Hy), RE(Hz), IM(Hz), RE(Ax), IM(Ax), RE(Ay), IM(Ay), RE(Az), IM(Az)" ! Writing the header.

   DO i = 1, nptx

      ! Computing x.
      argx = sizex*(i - 1.0D0)/(nptx - 1.0D0)

      DO j = 1, npty

         ! Computing y.
         argy = sizey*(j - 1.0D0)/(npty - 1.0D0)

         DO k = 1, nptz

            ! Computing z.
            argz = sizez*(k - 1.0D0)/(nptz - 1.0D0)

            WRITE (18, '(22(ES24.16,", "), ES24.16)') &
               argx, argy, argz, &
               REAL(xx(1, i, j, k)), AIMAG(xx(1, i, j, k)), &
               REAL(xx(2, i, j, k)), AIMAG(xx(2, i, j, k)), &
               REAL(xx(3, i, j, k)), AIMAG(xx(3, i, j, k)), &
               REAL(xx(4, i, j, k)), AIMAG(xx(4, i, j, k)), &
               REAL(xx(5, i, j, k)), AIMAG(xx(5, i, j, k)), &
               REAL(xx(6, i, j, k)), AIMAG(xx(6, i, j, k)), &
               REAL(xx(7, i, j, k)), AIMAG(xx(7, i, j, k)), &
               REAL(xx(8, i, j, k)), AIMAG(xx(8, i, j, k)), &
               REAL(xx(9, i, j, k)), AIMAG(xx(9, i, j, k)), &
               REAL(xx(10, i, j, k)), AIMAG(xx(10, i, j, k))

         END DO
      END DO
   END DO

   CLOSE (18)

END SUBROUTINE write10dvector
