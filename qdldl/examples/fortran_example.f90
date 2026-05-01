! Standalone Fortran example for the QDLDL Fortran bindings.
!
! Solves the 4x4 tridiagonal system:
!
!   A = [ 4 -1  0  0 ]   b = [3, 1, 5, 6]
!       [-1  4 -1  0 ]
!       [ 0 -1  4 -1 ]
!       [ 0  0 -1  4 ]
!
! b = A * [1, 1, 2, 2], so the exact solution is x = [1, 1, 2, 2].
!
! The matrix is stored in Fortran 1-based CRS format (full symmetric).
! Three calls to C_FORTRAN_QDLDL demonstrate the full workflow:
!   iopt=1  LDL^T factorization
!   iopt=2  triangular solve (factors reused)
!   iopt=3  free all storage

      PROGRAM QDLDL_FORTRAN_EXAMPLE

      IMPLICIT NONE

      ! Interface to the C bridge (no Fortran module needed)
      EXTERNAL  C_FORTRAN_QDLDL
      INTEGER   C_QDLDL_GET_NEG_COUNT
      EXTERNAL  C_QDLDL_GET_NEG_COUNT

      INTEGER, PARAMETER :: N    = 4
      INTEGER, PARAMETER :: NNZ  = 10   ! full symmetric (includes both triangles)
      INTEGER, PARAMETER :: NRHS = 1

      ! CRS arrays (1-based)
      INTEGER :: COLPTR(N+1) = (/ 1, 3, 6, 9, 11 /)
      INTEGER :: ROWIND(NNZ) = (/ 1, 2,  1, 2, 3,  2, 3, 4,  3, 4 /)
      REAL(8) :: VALUES(NNZ) = (/ 4.0D0, -1.0D0, &
                                   -1.0D0,  4.0D0, -1.0D0, &
                                   -1.0D0,  4.0D0, -1.0D0, &
                                   -1.0D0,  4.0D0 /)

      REAL(8) :: RHS(N) = (/ 3.0D0, 1.0D0, 5.0D0, 6.0D0 /)
      REAL(8) :: EXACT(N) = (/ 1.0D0, 1.0D0, 2.0D0, 2.0D0 /)

      INTEGER(8) :: F_FACTORS
      INTEGER    :: IOPT, INFO, LDB
      INTEGER    :: NEG_COUNT
      REAL(8)    :: RESID, MAXRESID
      INTEGER    :: I

      F_FACTORS = 0
      LDB       = N

      ! --- factorize ---
      IOPT = 1
      CALL C_FORTRAN_QDLDL( IOPT, N, NNZ, NRHS, VALUES, ROWIND, COLPTR, &
                             RHS, LDB, F_FACTORS, INFO )
      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'QDLDL factorization failed, INFO =', INFO
        CALL EXIT(1)
      END IF

      NEG_COUNT = C_QDLDL_GET_NEG_COUNT()
      WRITE(*,'(A,I0)') 'Negative D entries (inertia): ', NEG_COUNT

      ! --- solve (RHS overwritten with solution) ---
      IOPT = 2
      CALL C_FORTRAN_QDLDL( IOPT, N, NNZ, NRHS, VALUES, ROWIND, COLPTR, &
                             RHS, LDB, F_FACTORS, INFO )
      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'QDLDL solve failed, INFO =', INFO
        CALL EXIT(1)
      END IF

      ! --- print and verify ---
      WRITE(*,'(A)') 'Solution x:'
      MAXRESID = 0.0D0
      DO I = 1, N
        RESID = ABS(RHS(I) - EXACT(I))
        IF (RESID .GT. MAXRESID) MAXRESID = RESID
        WRITE(*,'(A,I1,A,F12.8,A,F6.1,A,ES10.2,A)') &
          '  x(', I, ') = ', RHS(I), '   exact =', EXACT(I), &
          '   |err| =', RESID, ''
      END DO
      WRITE(*,'(A,ES10.2)') 'Max absolute error: ', MAXRESID

      ! --- free storage ---
      IOPT = 3
      CALL C_FORTRAN_QDLDL( IOPT, N, NNZ, NRHS, VALUES, ROWIND, COLPTR, &
                             RHS, LDB, F_FACTORS, INFO )

      IF (MAXRESID .GT. 1.0D-12) THEN
        WRITE(*,*) 'FAIL: residual exceeds tolerance'
        CALL EXIT(1)
      END IF
      WRITE(*,'(A)') 'PASS'

      END PROGRAM QDLDL_FORTRAN_EXAMPLE
