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

      SUBROUTINE SYM_MAT_DECOMP_CHOLMOD ( CALLING_SUBR, MATIN_NAME, MATIN_SET, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_LOG, ERR, F04, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  BAILOUT
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1
      USE SUBR_BEGEND_LEVELS, ONLY    :  SYM_MAT_DECOMP_CHOLMOD_BEGEND
      USE CHOLMOD_STUF, ONLY          :  CHM_FACTORS
      USE SYM_MAT_DECOMP_CHOLMOD_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SYM_MAT_DECOMP_CHOLMOD'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_SET

      INTEGER(LONG), INTENT(IN)       :: NROWS
      INTEGER(LONG), INTENT(IN)       :: NTERMS
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)
      INTEGER(LONG), INTENT(INOUT)    :: INFO

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: POS
      INTEGER(LONG)                   :: GRIDV
      INTEGER(LONG)                   :: COMPV
      INTEGER(LONG)                   :: NTERM_CHOLMOD
      INTEGER(LONG), ALLOCATABLE      :: COL_COUNTS(:)
      INTEGER(LONG), ALLOCATABLE      :: COL_NEXT(:)
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SYM_MAT_DECOMP_CHOLMOD_BEGEND
      INTEGER(LONG)                   :: HS_SLOT
      REAL(DOUBLE)                    :: HS_T0
      REAL(DOUBLE)                    :: HS_PHASE_T0

      REAL(DOUBLE), INTENT(IN)        :: MATIN(NTERMS)

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

      CALL HOTSPOT_TIMER_BEGIN ( 'SYM_MAT_DECOMP_CHOLMOD', HS_SLOT, HS_T0 )

      IF (ALLOCATED(CCS1)) CALL DEALLOCATE_SCR_MAT ( 'CCS1' )
      ALLOCATE ( COL_COUNTS(NROWS), COL_NEXT(NROWS) )

      DO I=1,NROWS
         COL_COUNTS(I) = 0
         COL_NEXT(I)   = 0
      ENDDO

      NTERM_CHOLMOD = 0
      DO I=1,NROWS
         DO K=I_MATIN(I),I_MATIN(I+1)-1
            J = J_MATIN(K)
            IF (J >= I) THEN
               COL_COUNTS(J) = COL_COUNTS(J) + 1
               NTERM_CHOLMOD = NTERM_CHOLMOD + 1
            ENDIF
         ENDDO
      ENDDO

      CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NROWS, NTERM_CHOLMOD, SUBR_NAME )

      J_CCS1(1) = 1
      DO J=1,NROWS
         J_CCS1(J+1) = J_CCS1(J) + COL_COUNTS(J)
         COL_NEXT(J) = J_CCS1(J)
      ENDDO

      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      DO I=1,NROWS
         DO K=I_MATIN(I),I_MATIN(I+1)-1
            J = J_MATIN(K)
            IF (J >= I) THEN
               POS = COL_NEXT(J)
               I_CCS1(POS) = I
                 CCS1(POS) = MATIN(K)
               COL_NEXT(J) = POS + 1
            ENDIF
         ENDDO
      ENDDO
      CALL HOTSPOT_TIMER_ADD ( 'SYM_MAT_DECOMP_CHOLMOD/CRS_TO_CCS', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      CALL C_FORTRAN_DCHOLMOD_FACTOR ( NROWS, NTERM_CHOLMOD, CCS1, I_CCS1, J_CCS1, CHM_FACTORS, INFO )
      CALL HOTSPOT_TIMER_ADD ( 'SYM_MAT_DECOMP_CHOLMOD/CHOLMOD_CALL', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

      DEALLOCATE ( COL_COUNTS, COL_NEXT )
      CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

      IF (INFO == 0) THEN

         WRITE (SC1,9902) MATIN_NAME, SUBR_NAME
         WRITE (F06,9902) MATIN_NAME, SUBR_NAME

      ELSE IF (INFO < 0) THEN

         WRITE(ERR,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         WRITE(F06,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE

         CALL GET_GRID_AND_COMP ( MATIN_SET, INFO, GRIDV, COMPV )

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

      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      CALL HOTSPOT_TIMER_END ( HS_SLOT, HS_T0 )

      RETURN

  981 FORMAT(' *ERROR   981: THE FACTORIZATION OF THE MATRIX ',A,' BY CHOLMOD FAILED AT ROW/COL = ', I12, '.')

 9811 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMP ',I3,'. THE CALLING SUBR WAS: ',A,/)

 9812 FORMAT('               THIS IS FOR ROW AND COL ',I8,' IN THE MATRIX. THE CALLING SUBR WAS: ',A,/)

 9902 FORMAT(' CHOLMOD FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)

 9903 FORMAT(' *ERROR  9903: CHOLMOD SPARSE SOLVER HAS FAILED WITH INFO = ', I12,' IN SUBR ', A, ' CALLED BY SUBR ', A)

99999 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE MESSAGES AND BULK DATA PARAMETER BAILOUT = ',I7)

      END SUBROUTINE SYM_MAT_DECOMP_CHOLMOD
