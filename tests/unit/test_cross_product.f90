! Unit tests for SUBROUTINE CROSS (Source/UTIL/CROSS.f90).
!
! CROSS computes the vector cross product C = A × B.  The logging block
! (IF WRT_LOG >= SUBR_BEGEND) is never entered because the stub sets WRT_LOG = 0
! and CROSS_BEGEND = 11 > 0.
!
! Tests cover:
!   1. Unit-axis pairs: x×y=z, y×z=x, z×x=y
!   2. Anti-commutativity: A×B = -(B×A)
!   3. Self cross-product: A×A = 0
!   4. Scalar linearity: (α·A)×B = α·(A×B)

      PROGRAM test_cross_product

      USE PENTIUM_II_KIND, ONLY  :  LONG, DOUBLE

      IMPLICIT NONE

      REAL(DOUBLE) :: A(3), B(3), C(3), C2(3), D(3)
      INTEGER      :: n_failed, I

      REAL(DOUBLE), PARAMETER :: TOL = 1.0D-14

      n_failed = 0

! --------------------------------------------------------------------
! Test 1a: x × y = z   →  [1,0,0] × [0,1,0] = [0,0,1]
! --------------------------------------------------------------------
      A = [1.0D0, 0.0D0, 0.0D0]
      B = [0.0D0, 1.0D0, 0.0D0]
      CALL CROSS(A, B, C)
      IF (ABS(C(1) - 0.0D0) > TOL .OR. &
          ABS(C(2) - 0.0D0) > TOL .OR. &
          ABS(C(3) - 1.0D0) > TOL) THEN
        WRITE(*,'(A)') 'FAIL test_cross_product: Test 1a (x cross y) wrong'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
! Test 1b: y × z = x   →  [0,1,0] × [0,0,1] = [1,0,0]
! --------------------------------------------------------------------
      A = [0.0D0, 1.0D0, 0.0D0]
      B = [0.0D0, 0.0D0, 1.0D0]
      CALL CROSS(A, B, C)
      IF (ABS(C(1) - 1.0D0) > TOL .OR. &
          ABS(C(2) - 0.0D0) > TOL .OR. &
          ABS(C(3) - 0.0D0) > TOL) THEN
        WRITE(*,'(A)') 'FAIL test_cross_product: Test 1b (y cross z) wrong'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
! Test 1c: z × x = y   →  [0,0,1] × [1,0,0] = [0,1,0]
! --------------------------------------------------------------------
      A = [0.0D0, 0.0D0, 1.0D0]
      B = [1.0D0, 0.0D0, 0.0D0]
      CALL CROSS(A, B, C)
      IF (ABS(C(1) - 0.0D0) > TOL .OR. &
          ABS(C(2) - 1.0D0) > TOL .OR. &
          ABS(C(3) - 0.0D0) > TOL) THEN
        WRITE(*,'(A)') 'FAIL test_cross_product: Test 1c (z cross x) wrong'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
! Test 2: Anti-commutativity — A×B = -(B×A) for a general pair
! A = [1, 2, 3], B = [4, -1, 0.5]
! --------------------------------------------------------------------
      A = [1.0D0,  2.0D0, 3.0D0]
      B = [4.0D0, -1.0D0, 0.5D0]
      CALL CROSS(A, B, C)   ! C  = A × B
      CALL CROSS(B, A, C2)  ! C2 = B × A  (should be -C)
      DO I = 1, 3
        IF (ABS(C(I) + C2(I)) > TOL) THEN
          WRITE(*,'(A,I0,A)') &
            'FAIL test_cross_product: Test 2 anti-commutativity, component ', I, ' wrong'
          n_failed = n_failed + 1
        END IF
      END DO

! --------------------------------------------------------------------
! Test 3: Self cross-product — A×A = 0 for any A
! A = [1.5, -3.0, 2.7]
! --------------------------------------------------------------------
      A = [1.5D0, -3.0D0, 2.7D0]
      CALL CROSS(A, A, C)
      IF (ABS(C(1)) > TOL .OR. ABS(C(2)) > TOL .OR. ABS(C(3)) > TOL) THEN
        WRITE(*,'(A)') 'FAIL test_cross_product: Test 3 (A cross A) not zero'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
! Test 4: Scalar linearity — (α·A) × B = α · (A × B)
! A = [1, 0, -1], B = [0, 1, 1], α = 3.5
! --------------------------------------------------------------------
      A  = [1.0D0, 0.0D0, -1.0D0]
      B  = [0.0D0, 1.0D0,  1.0D0]
      D  = 3.5D0 * A   ! scaled A

      CALL CROSS(A, B, C)   ! C  = A × B
      CALL CROSS(D, B, C2)  ! C2 = (3.5·A) × B

      DO I = 1, 3
        IF (ABS(C2(I) - 3.5D0*C(I)) > TOL) THEN
          WRITE(*,'(A,I0,A)') &
            'FAIL test_cross_product: Test 4 scalar linearity, component ', I, ' wrong'
          n_failed = n_failed + 1
        END IF
      END DO

! --------------------------------------------------------------------
! Test 5: Known numerical result
! A = [2, 3, 4], B = [5, 6, 7]
! A×B = [3*7-4*6, 4*5-2*7, 2*6-3*5] = [21-24, 20-14, 12-15] = [-3, 6, -3]
! --------------------------------------------------------------------
      A = [2.0D0, 3.0D0, 4.0D0]
      B = [5.0D0, 6.0D0, 7.0D0]
      CALL CROSS(A, B, C)
      IF (ABS(C(1) - (-3.0D0)) > TOL .OR. &
          ABS(C(2) -   6.0D0)  > TOL .OR. &
          ABS(C(3) - (-3.0D0)) > TOL) THEN
        WRITE(*,'(A)') 'FAIL test_cross_product: Test 5 known numerical result wrong'
        n_failed = n_failed + 1
      END IF

! --------------------------------------------------------------------
      IF (n_failed > 0) THEN
        WRITE(*,'(A,I0,A)') 'test_cross_product: ', n_failed, ' test(s) FAILED'
        STOP 1
      END IF

      WRITE(*,'(A)') 'PASS test_cross_product'

      END PROGRAM test_cross_product
