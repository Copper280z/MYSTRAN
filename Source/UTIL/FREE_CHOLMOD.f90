! ##################################################################################################################################
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

      SUBROUTINE FREE_CHOLMOD ( CALLING_SUBR, MATIN_NAME, INFO )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  WRT_LOG, ERR, F04, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CHOLMOD_STUF, ONLY          :  CHM_FACTORS
      USE FREE_CHOLMOD_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FREE_CHOLMOD'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME
      INTEGER(LONG), INTENT(INOUT)    :: INFO

! **********************************************************************************************************************************
      IF (WRT_LOG > 0) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

      CALL C_FORTRAN_DCHOLMOD_FREE ( CHM_FACTORS, INFO )

      IF (INFO /= 0) THEN
         WRITE(ERR,9903) INFO, TRIM(MATIN_NAME), TRIM(CALLING_SUBR)
         WRITE(F06,9903) INFO, TRIM(MATIN_NAME), TRIM(CALLING_SUBR)
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      IF (WRT_LOG > 0) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

 9903 FORMAT(' *ERROR  9903: CHOLMOD FREE STORAGE FAILED WITH INFO = ', I12,' FOR MATRIX ', A, ' IN SUBR ', A)

      END SUBROUTINE FREE_CHOLMOD
