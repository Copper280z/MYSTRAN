! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

      SUBROUTINE SYM_MAT_DECOMP_QDLDL ( CALLING_SUBR, MATIN_NAME, MATIN_SET, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, INFO )

! Decomposes a symmetric matrix into LDL^T triangular factors using QDLDL with METIS ordering.
! The input matrix, MATIN, is stored in CRS sparse format (full, non-symmetric storage).
! After factorization:
!   - QDLDL_NEG_COUNT (in SuperLU_STUF) holds the number of negative D diagonal entries,
!     which equals the number of eigenvalues below the current shift (Sylvester's law of inertia).
!   - MAXRATIO is printed to stdout by the C bridge.
! Shares the SLU_FACTORS handle with the SuperLU path (same storage, different C routine).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  BAILOUT
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, QDLDL_NEG_COUNT

      USE SYM_MAT_DECOMP_QDLDL_USE_IFs

      IMPLICIT NONE

      INTEGER(LONG), EXTERNAL         :: C_QDLDL_GET_NEG_COUNT

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SYM_MAT_DECOMP_QDLDL'

      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_SET

      INTEGER(LONG)                   :: I
      INTEGER(LONG), INTENT(IN)       :: NROWS
      INTEGER(LONG), INTENT(IN)       :: NTERMS
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)
      INTEGER(LONG), INTENT(INOUT)    :: INFO

      INTEGER(LONG)                   :: COMPV
      INTEGER(LONG)                   :: GRIDV

      REAL(DOUBLE),  INTENT(IN)       :: MATIN(NTERMS)
      REAL(DOUBLE)                    :: DUM_COL(NROWS)

! **********************************************************************************************************************************

      DO I=1,NROWS
         DUM_COL(I) = ZERO
      ENDDO

      CALL C_FORTRAN_QDLDL( 1, NROWS, NTERMS, 1, MATIN, J_MATIN, I_MATIN, DUM_COL, NROWS, SLU_FACTORS, INFO )

      ! Retrieve inertia (number of eigenvalues below current shift)
      QDLDL_NEG_COUNT = C_QDLDL_GET_NEG_COUNT()

      IF (INFO == 0) THEN

         WRITE (SC1,9902) MATIN_NAME, SUBR_NAME
         WRITE (F06,9902) MATIN_NAME, SUBR_NAME

      ELSE IF (INFO < 0) THEN

         WRITE(SC1,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         WRITE(ERR,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         WRITE(F06,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE IF (INFO > 0) THEN

         CALL GET_GRID_AND_COMP ( MATIN_SET, INFO, GRIDV, COMPV  )

         WRITE(ERR,981) MATIN_NAME, INFO
         WRITE(F06,981) MATIN_NAME, INFO
         IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
            WRITE(ERR,9811) GRIDV, COMPV, CALLING_SUBR
            WRITE(F06,9811) GRIDV, COMPV, CALLING_SUBR
         ELSE
            WRITE(ERR,9812) INFO, CALLING_SUBR
            WRITE(F06,9812) INFO, CALLING_SUBR
         ENDIF

         IF (BAILOUT >= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,99999) BAILOUT
            WRITE(F06,99999) BAILOUT
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ENDIF

      RETURN

!***********************************************************************************************************************************
  981 FORMAT(' *ERROR   981: THE FACTORIZATION OF THE MATRIX ',A,' BY QDLDL HAD ERROR WITH INFO = ', I12, '.')

 9811 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMP ',I3,'. THE CALLING SUBR WAS: ',A,/)

 9812 FORMAT('               THIS IS FOR ROW AND COL ',I8,' IN THE MATRIX. THE CALLING SUBR WAS: ',A,/)

 9902 FORMAT(' QDLDL FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)

 9903 FORMAT(' *ERROR  9903: QDLDL SPARSE SOLVER HAS FAILED WITH INFO = ', I12,' IN SUBR ', A, ' CALLED BY SUBR ', A)

99999 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE MESSAGES AND BULK DATA PARAMETER BAILOUT = ',I7)

!***********************************************************************************************************************************

      END SUBROUTINE SYM_MAT_DECOMP_QDLDL
