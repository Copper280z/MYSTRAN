! ###############################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

   MODULE FBS_CHOLMOD_Interface

   INTERFACE

      SUBROUTINE FBS_CHOLMOD ( CALLING_SUBR, MATIN_NAME, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, ICOL, RHS_COL, INFO )

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_LOG, ERR, F04, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  FBS_CHOLMOD_BEGEND
      USE CHOLMOD_STUF, ONLY          :  CHM_FACTORS

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME

      INTEGER(LONG), INTENT(IN)       :: ICOL
      INTEGER(LONG), INTENT(IN)       :: NROWS
      INTEGER(LONG), INTENT(IN)       :: NTERMS
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)
      INTEGER(LONG), INTENT(INOUT)    :: INFO
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = FBS_CHOLMOD_BEGEND

      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERMS)
      REAL(DOUBLE) , INTENT(INOUT)    :: RHS_COL(NROWS)

      END SUBROUTINE FBS_CHOLMOD

   END INTERFACE

   END MODULE FBS_CHOLMOD_Interface
