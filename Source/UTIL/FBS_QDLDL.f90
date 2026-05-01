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

      SUBROUTINE FBS_QDLDL ( CALLING_SUBR, MATIN_NAME, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, ICOL, RHS_COL, INFO )

! Forward-backward substitution using pre-computed QDLDL LDL^T factors.
! Solves A*x = b using the factorization stored by SYM_MAT_DECOMP_QDLDL.
! RHS_COL is overwritten with the solution on return.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS
      USE FBS_QDLDL_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FBS_QDLDL'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME

      INTEGER(LONG), INTENT(IN)       :: ICOL
      INTEGER(LONG), INTENT(IN)       :: NROWS
      INTEGER(LONG), INTENT(IN)       :: NTERMS
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)
      INTEGER(LONG), INTENT(INOUT)    :: INFO

      REAL(DOUBLE),  INTENT(IN)       :: MATIN(NTERMS)
      REAL(DOUBLE),  INTENT(IN)       :: RHS_COL(NROWS)

! **********************************************************************************************************************************

      CALL C_FORTRAN_QDLDL( 2, NROWS, NTERMS, 1, MATIN, J_MATIN, I_MATIN, RHS_COL, NROWS, SLU_FACTORS, INFO )

      IF (INFO .NE. 0) THEN
         WRITE( * ,*) 'INFO FROM QDLDL TRIANGULAR SOLVE = ', INFO
         WRITE(ERR,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         WRITE(F06,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      RETURN

!***********************************************************************************************************************************
 9903 FORMAT(' *ERROR  9903: QDLDL SPARSE SOLVER HAS FAILED WITH INFO = ', I11,' IN SUBR ', A, ' CALLED BY SUBR ', A)

!***********************************************************************************************************************************

      END SUBROUTINE FBS_QDLDL
