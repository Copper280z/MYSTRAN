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

      SUBROUTINE SPARSE_KGGD

! Converts the hash-assembled KGGD differential stiffness matrix to CRS format.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE, DBL_LONG
      USE IOUNT1, ONLY                :  ERR, F04, F06, SC1, WRT_ERR, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NGRID, NTERM_KGGD, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  SPARSE_KGGD_BEGEND
      USE PARAMS, ONLY                :  EPSIL, PRTSTIFF, SUPINFO
      USE STF_ARRAYS, ONLY            :  STF_ROW_HM, STF_COL_HM, STF_VAL_HM
      USE SPARSE_MATRICES, ONLY       :  I_KGGD, J_KGGD, KGGD

      USE SPARSE_KGGD_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_KGGD'

      INTEGER(LONG)                   :: I, J
      INTEGER(LONG)                   :: KTERM_KGGD
      INTEGER(LONG)                   :: NUM_MAX = 0
      INTEGER(LONG)                   :: NUM_NONZERO_IN_ROW
      INTEGER(LONG)                   :: NZERO
      INTEGER(LONG)                   :: POS
      INTEGER(LONG)                   :: RAW_NTERM
      INTEGER(LONG), ALLOCATABLE      :: ROW_NEXT(:)
      INTEGER(LONG)                   :: ROW_END
      INTEGER(LONG)                   :: ROW_START
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SPARSE_KGGD_BEGEND
      REAL(DOUBLE)                    :: EPS1

      INTRINSIC                       :: DABS

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
      RAW_NTERM = NTERM_KGGD

      NZERO = 0
      POS = 0
      DO I=1,RAW_NTERM
         IF (DABS(STF_VAL_HM(I)) < EPS1) THEN
            NZERO = NZERO + 1
         ELSE
            POS = POS + 1
            STF_ROW_HM(POS) = STF_ROW_HM(I)
            STF_COL_HM(POS) = STF_COL_HM(I)
            STF_VAL_HM(POS) = STF_VAL_HM(I)
         ENDIF
      ENDDO
      NTERM_KGGD = POS

      WRITE(ERR,146) NTERM_KGGD
      IF (SUPINFO == 'N') THEN
         WRITE(F06,146) NTERM_KGGD
      ENDIF

      IF (NTERM_KGGD <= 0) THEN
         WRITE(ERR,1611) NTERM_KGGD
         WRITE(F06,1611) NTERM_KGGD
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL ALLOCATE_SPARSE_MAT ( 'KGGD', NDOFG, NTERM_KGGD, SUBR_NAME )

      DO I=1,NDOFG+1
         I_KGGD(I) = 0
      ENDDO
      DO I=1,NTERM_KGGD
         I_KGGD(STF_ROW_HM(I)+1) = I_KGGD(STF_ROW_HM(I)+1) + 1
      ENDDO
      I_KGGD(1) = 1
      DO I=1,NDOFG
         I_KGGD(I+1) = I_KGGD(I+1) + I_KGGD(I)
      ENDDO
      ALLOCATE ( ROW_NEXT(NDOFG) )
      DO I=1,NDOFG
         ROW_NEXT(I) = I_KGGD(I)
      ENDDO
      DO I=1,NTERM_KGGD
         POS = ROW_NEXT(STF_ROW_HM(I))
         J_KGGD(POS) = STF_COL_HM(I)
           KGGD(POS) = STF_VAL_HM(I)
         ROW_NEXT(STF_ROW_HM(I)) = POS + 1
      ENDDO
      DO I=1,NDOFG
         ROW_START = I_KGGD(I)
         ROW_END   = I_KGGD(I+1) - 1
         NUM_NONZERO_IN_ROW = ROW_END - ROW_START + 1
         IF (NUM_NONZERO_IN_ROW > 1) THEN
            CALL SORT_INT1_REAL1 ( SUBR_NAME, 'KGGD row cols', NUM_NONZERO_IN_ROW, J_KGGD(ROW_START:ROW_END),                   &
                                   KGGD(ROW_START:ROW_END) )
         ENDIF
      ENDDO
      DEALLOCATE ( ROW_NEXT )

      KTERM_KGGD = 0
      CALL COUNTER_INIT('     Working on grid ', NGRID)
      DO I=1,NDOFG
         ROW_START = I_KGGD(I)
         ROW_END   = I_KGGD(I+1) - 1
         NUM_NONZERO_IN_ROW = MAX(0_LONG, ROW_END - ROW_START + 1)
         IF (NUM_NONZERO_IN_ROW > NUM_MAX) NUM_MAX = NUM_NONZERO_IN_ROW
         DO J=ROW_START,ROW_END
            KTERM_KGGD = KTERM_KGGD + 1
         ENDDO
         IF (I <= NGRID) CALL COUNTER_PROGRESS(I)
      ENDDO

      WRITE(SC1,*) CR13

      IF (PRTSTIFF(1) >= 1) THEN
         CALL WRITE_SPARSE_CRS ( 'STIFFNESS MATRIX KGGD', 'G ', 'G ', NTERM_KGGD, NDOFG, I_KGGD, J_KGGD, KGGD )
      ENDIF

      WRITE(ERR,101) NUM_MAX
      IF (SUPINFO == 'N') THEN
         WRITE(F06,101) NUM_MAX
      ENDIF

      IF (ALLOCATED(STF_ROW_HM)) DEALLOCATE ( STF_ROW_HM )
      IF (ALLOCATED(STF_COL_HM)) DEALLOCATE ( STF_COL_HM )
      IF (ALLOCATED(STF_VAL_HM)) DEALLOCATE ( STF_VAL_HM )

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: MAX NUMBER OF NONZERO TERMS IN A ROW OF THE G-SET STIFFNESS MATRIX     = ',I12,/)

  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE KGGD STIFFNESS MATRIX IS                 = ',I12,/)

 1611 FORMAT(' *ERROR  1611: THE G-SET DIFFERENTIAL STIFF MATRIX, KGGD, MUST HAVE SOME NONZERO TERMS. HOWEVER IT HAS ',I12,' TERMS')

      END SUBROUTINE SPARSE_KGGD
