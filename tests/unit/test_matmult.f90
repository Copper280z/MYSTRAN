! Unit tests for SUBROUTINE MATMULT_FFF and SUBROUTINE MATMULT_FFF_T
! (Source/UTIL/MATMULT_FFF.f90 and MATMULT_FFF_T.f90).
!
! MATMULT_FFF  computes C    = A  × B
! MATMULT_FFF_T computes C   = A' × B   (A is transposed)
!
! The logging block is gated by WRT_LOG=0 (stub), so OURTIM is never called.
!
! Tests cover:
!   1. Square identity: A × I₃ = A
!   2. Non-square product: 2×3 × 3×2 = 2×2 with known result
!   3. Result is zero when one input matrix is zero
!   4. MATMULT_FFF_T vs explicit transpose: (A'×B) = (A_transposed × B)
!   5. MATMULT_FFF_T non-square: 3×2 matrix A (so A' is 2×3) times 3×2 B → 2×2

      PROGRAM test_matmult

      USE PENTIUM_II_KIND, ONLY  :  LONG, DOUBLE

      IMPLICIT NONE

      INTEGER(LONG) :: NROWA, NCOLA, NCOLB
      INTEGER(LONG) :: I, J, K
      INTEGER       :: n_failed

      REAL(DOUBLE), PARAMETER :: TOL = 1.0D-12

      REAL(DOUBLE) :: A33(3,3), B33(3,3), C33(3,3), I33(3,3)
      REAL(DOUBLE) :: A23(2,3), B32(3,2), C22(2,2), expected22(2,2)
      REAL(DOUBLE) :: Z33(3,3), R33(3,3)
      REAL(DOUBLE) :: A32(3,2), B32b(3,2), Ct22(2,2), At23(2,3)

      n_failed = 0

! --------------------------------------------------------------------
! Test 1: A × I = A  (3×3 square, identity multiplication)
! A = [[1,2,3],[4,5,6],[7,8,9]]
! --------------------------------------------------------------------
      A33(1,:) = [1.0D0, 2.0D0, 3.0D0]
      A33(2,:) = [4.0D0, 5.0D0, 6.0D0]
      A33(3,:) = [7.0D0, 8.0D0, 9.0D0]

      ! Build identity
      I33 = 0.0D0
      DO I = 1, 3
        I33(I,I) = 1.0D0
      END DO

      NROWA = 3;  NCOLA = 3;  NCOLB = 3
      CALL MATMULT_FFF(A33, I33, NROWA, NCOLA, NCOLB, C33)

      DO I = 1, 3
        DO J = 1, 3
          IF (ABS(C33(I,J) - A33(I,J)) > TOL) THEN
            WRITE(*,'(A)') 'FAIL test_matmult: Test 1 (A × I ≠ A)'
            n_failed = n_failed + 1
            GOTO 100
          END IF
        END DO
      END DO
 100  CONTINUE

! --------------------------------------------------------------------
! Test 2: Non-square product 2×3 × 3×2 = 2×2
! A = [[1,0,2],[0,3,1]],  B = [[1,1],[0,1],[2,0]]
! Expected C = A*B:
!   C(1,1) = 1*1+0*0+2*2 = 5,  C(1,2) = 1*1+0*1+2*0 = 1
!   C(2,1) = 0*1+3*0+1*2 = 2,  C(2,2) = 0*1+3*1+1*0 = 3
! --------------------------------------------------------------------
      A23(1,:) = [1.0D0, 0.0D0, 2.0D0]
      A23(2,:) = [0.0D0, 3.0D0, 1.0D0]
      B32(1,:) = [1.0D0, 1.0D0]
      B32(2,:) = [0.0D0, 1.0D0]
      B32(3,:) = [2.0D0, 0.0D0]
      expected22(1,:) = [5.0D0, 1.0D0]
      expected22(2,:) = [2.0D0, 3.0D0]

      NROWA = 2;  NCOLA = 3;  NCOLB = 2
      CALL MATMULT_FFF(A23, B32, NROWA, NCOLA, NCOLB, C22)

      DO I = 1, 2
        DO J = 1, 2
          IF (ABS(C22(I,J) - expected22(I,J)) > TOL) THEN
            WRITE(*,'(A,I0,A,I0,A,F10.4,A,F10.4,A)') &
              'FAIL test_matmult: Test 2 element (', I, ',', J, &
              ') = ', C22(I,J), ' expected ', expected22(I,J), ''
            n_failed = n_failed + 1
          END IF
        END DO
      END DO

! --------------------------------------------------------------------
! Test 3: Result is zero when one factor is the zero matrix
! --------------------------------------------------------------------
      Z33 = 0.0D0
      NROWA = 3;  NCOLA = 3;  NCOLB = 3
      CALL MATMULT_FFF(A33, Z33, NROWA, NCOLA, NCOLB, R33)

      IF (ANY(ABS(R33) > TOL)) THEN
        WRITE(*,'(A)') 'FAIL test_matmult: Test 3 (A × 0 ≠ 0)'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
! Test 4: MATMULT_FFF_T vs explicit transpose
! Use the same A23 (2×3) and B32 (3×2) from Test 2.
! MATMULT_FFF_T(A, B, NROWA=2, NCOLA=3, NCOLB=2, C) computes A'(3×2) × B(2×2)?
! Wait — MATMULT_FFF_T signature: A(NROWA,NCOLA), B(NROWA,NCOLB), C(NCOLA,NCOLB)
! so it computes A'(NCOLA×NROWA) × B(NROWA×NCOLB) → C(NCOLA×NCOLB)
!
! Use A32 (3×2) and B32b (3×2):
! MATMULT_FFF_T(A32, B32b, NROWA=3, NCOLA=2, NCOLB=2, C22) computes
!   A32'(2×3) × B32b(3×2) → C22(2×2)
!
! A32 = [[1,4],[2,5],[3,6]], B32b = [[1,1],[0,1],[2,0]]
! A32' = [[1,2,3],[4,5,6]]
! C = A32' × B32b:
!   C(1,1) = 1*1+2*0+3*2 = 7,  C(1,2) = 1*1+2*1+3*0 = 3
!   C(2,1) = 4*1+5*0+6*2 = 16, C(2,2) = 4*1+5*1+6*0 = 9
! --------------------------------------------------------------------
      A32(1,:) = [1.0D0, 4.0D0]
      A32(2,:) = [2.0D0, 5.0D0]
      A32(3,:) = [3.0D0, 6.0D0]
      B32b(1,:) = [1.0D0, 1.0D0]
      B32b(2,:) = [0.0D0, 1.0D0]
      B32b(3,:) = [2.0D0, 0.0D0]

      NROWA = 3;  NCOLA = 2;  NCOLB = 2
      CALL MATMULT_FFF_T(A32, B32b, NROWA, NCOLA, NCOLB, Ct22)

      expected22(1,:) = [7.0D0,  3.0D0]
      expected22(2,:) = [16.0D0, 9.0D0]

      DO I = 1, 2
        DO J = 1, 2
          IF (ABS(Ct22(I,J) - expected22(I,J)) > TOL) THEN
            WRITE(*,'(A,I0,A,I0,A,F10.4,A,F10.4,A)') &
              'FAIL test_matmult: Test 4 MATMULT_FFF_T element (', I, ',', J, &
              ') = ', Ct22(I,J), ' expected ', expected22(I,J), ''
            n_failed = n_failed + 1
          END IF
        END DO
      END DO

! --------------------------------------------------------------------
! Test 5: MATMULT_FFF_T vs explicit MATMULT_FFF with manually transposed A
! Confirms the two routines agree on a 3×3 case
! A = [[2,3,1],[0,4,5],[1,2,6]]
! B = [[1,0,2],[3,1,0],[0,2,1]]
! Compute C1 = MATMULT_FFF_T(A, B, 3, 3, 3)  (A'×B)
! Compute A_T manually, then C2 = MATMULT_FFF(A_T, B, 3, 3, 3)
! C1 should equal C2 to within tolerance
! --------------------------------------------------------------------
      BLOCK
        REAL(DOUBLE) :: A3x3(3,3), B3x3(3,3), AT3x3(3,3)
        REAL(DOUBLE) :: C1_33(3,3), C2_33(3,3)
        INTEGER(LONG) :: ii, jj

        A3x3(1,:) = [2.0D0, 3.0D0, 1.0D0]
        A3x3(2,:) = [0.0D0, 4.0D0, 5.0D0]
        A3x3(3,:) = [1.0D0, 2.0D0, 6.0D0]
        B3x3(1,:) = [1.0D0, 0.0D0, 2.0D0]
        B3x3(2,:) = [3.0D0, 1.0D0, 0.0D0]
        B3x3(3,:) = [0.0D0, 2.0D0, 1.0D0]

        ! Transpose A manually
        DO ii = 1, 3
          DO jj = 1, 3
            AT3x3(ii,jj) = A3x3(jj,ii)
          END DO
        END DO

        CALL MATMULT_FFF_T(A3x3, B3x3, 3_LONG, 3_LONG, 3_LONG, C1_33)
        CALL MATMULT_FFF  (AT3x3, B3x3, 3_LONG, 3_LONG, 3_LONG, C2_33)

        DO ii = 1, 3
          DO jj = 1, 3
            IF (ABS(C1_33(ii,jj) - C2_33(ii,jj)) > TOL) THEN
              WRITE(*,'(A)') &
                'FAIL test_matmult: Test 5 MATMULT_FFF_T disagrees with explicit transpose'
              n_failed = n_failed + 1
              GOTO 200
            END IF
          END DO
        END DO
 200    CONTINUE
      END BLOCK

! --------------------------------------------------------------------
      IF (n_failed > 0) THEN
        WRITE(*,'(A,I0,A)') 'test_matmult: ', n_failed, ' test(s) FAILED'
        STOP 1
      END IF

      WRITE(*,'(A)') 'PASS test_matmult'

      END PROGRAM test_matmult
