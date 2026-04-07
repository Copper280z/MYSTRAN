! Unit tests for SUBROUTINE OUTER_PRODUCT (Source/UTIL/OUTER_PRODUCT.f90).
!
! OUTER_PRODUCT has no logging or I/O dependencies — it only uses PENTIUM_II_KIND.
! Tests verify: correct values for a general 3×2 case, the rank-1 self-outer-product
! of a vector, and that every element satisfies C(i,j) = A(i)*B(j).

      PROGRAM test_outer_product

      USE PENTIUM_II_KIND, ONLY  :  LONG, DOUBLE

      IMPLICIT NONE

      INTEGER(LONG) :: NA, NB, I, J
      REAL(DOUBLE)  :: A3(3), B2(2), C32(3,2)
      REAL(DOUBLE)  :: A4(4), C44(4,4)
      INTEGER       :: n_failed

      REAL(DOUBLE), PARAMETER :: TOL = 1.0D-14

      n_failed = 0

! --------------------------------------------------------------------
! Test 1: 3×2 outer product — explicit expected values
! A = [1, 2, 3],  B = [4, 5]
! Expected C:
!   row 1: [4,  5]
!   row 2: [8, 10]
!   row 3: [12, 15]
! --------------------------------------------------------------------
      NA = 3;  NB = 2
      A3 = [1.0D0, 2.0D0, 3.0D0]
      B2 = [4.0D0, 5.0D0]

      CALL OUTER_PRODUCT(A3, B2, NA, NB, C32)

      IF (ABS(C32(1,1) -  4.0D0) > TOL .OR. &
          ABS(C32(1,2) -  5.0D0) > TOL .OR. &
          ABS(C32(2,1) -  8.0D0) > TOL .OR. &
          ABS(C32(2,2) - 10.0D0) > TOL .OR. &
          ABS(C32(3,1) - 12.0D0) > TOL .OR. &
          ABS(C32(3,2) - 15.0D0) > TOL) THEN
        WRITE(*,'(A)') 'FAIL test_outer_product: Test 1 (3x2 product) wrong values'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
! Test 2: Every element satisfies C(i,j) = A(i)*B(j)  (property check)
! Use vectors with non-trivial entries
! --------------------------------------------------------------------
      NA = 3;  NB = 2
      A3 = [-1.5D0, 0.0D0, 3.7D0]
      B2 = [ 2.0D0, -4.0D0]

      CALL OUTER_PRODUCT(A3, B2, NA, NB, C32)

      DO I = 1, NA
        DO J = 1, NB
          IF (ABS(C32(I,J) - A3(I)*B2(J)) > TOL) THEN
            WRITE(*,'(A,I0,A,I0,A)') &
              'FAIL test_outer_product: Test 2 element (', I, ',', J, ') wrong'
            n_failed = n_failed + 1
          END IF
        END DO
      END DO

! --------------------------------------------------------------------
! Test 3: Self outer product — result must be symmetric
! A = [1, -2, 3, 0],  C = A⊗A  →  C(i,j) = C(j,i)
! --------------------------------------------------------------------
      NA = 4;  NB = 4
      A4 = [1.0D0, -2.0D0, 3.0D0, 0.0D0]

      CALL OUTER_PRODUCT(A4, A4, NA, NB, C44)

      DO I = 1, NA
        DO J = 1, NB
          IF (ABS(C44(I,J) - C44(J,I)) > TOL) THEN
            WRITE(*,'(A)') &
              'FAIL test_outer_product: Test 3 self-outer-product not symmetric'
            n_failed = n_failed + 1
            EXIT
          END IF
        END DO
      END DO

! --------------------------------------------------------------------
! Test 4: Zero vector — all outputs must be zero
! --------------------------------------------------------------------
      NA = 3;  NB = 2
      A3 = [0.0D0, 0.0D0, 0.0D0]
      B2 = [7.0D0, -3.0D0]

      CALL OUTER_PRODUCT(A3, B2, NA, NB, C32)

      IF (ANY(ABS(C32) > TOL)) THEN
        WRITE(*,'(A)') 'FAIL test_outer_product: Test 4 zero-A product not zero'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
      IF (n_failed > 0) THEN
        WRITE(*,'(A,I0,A)') 'test_outer_product: ', n_failed, ' test(s) FAILED'
        STOP 1
      END IF

      WRITE(*,'(A)') 'PASS test_outer_product'

      END PROGRAM test_outer_product
