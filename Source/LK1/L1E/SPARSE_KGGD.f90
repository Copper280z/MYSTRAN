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

      IF (NTERM_KGGD > 1) THEN
         CALL SORT_INT2_REAL1 ( SUBR_NAME, 'KGGD hash triplets', NTERM_KGGD, STF_ROW_HM(1:NTERM_KGGD), STF_COL_HM(1:NTERM_KGGD),  &
                                STF_VAL_HM(1:NTERM_KGGD) )
         POS = 1
         DO WHILE (POS <= NTERM_KGGD)
            ROW_START = POS
            DO WHILE ((POS <= NTERM_KGGD) .AND. (STF_ROW_HM(POS) == STF_ROW_HM(ROW_START)))
               POS = POS + 1
            ENDDO
            ROW_END = POS - 1
            NUM_NONZERO_IN_ROW = ROW_END - ROW_START + 1
            IF (NUM_NONZERO_IN_ROW > 1) THEN
               CALL SORT_INT1_REAL1 ( SUBR_NAME, 'KGGD row cols', NUM_NONZERO_IN_ROW, STF_COL_HM(ROW_START:ROW_END),             &
                                      STF_VAL_HM(ROW_START:ROW_END) )
            ENDIF
         ENDDO
      ENDIF

      KTERM_KGGD = 0
      POS = 1
      I_KGGD(1) = 1
      CALL COUNTER_INIT('     Working on grid ', NGRID)
i_do: DO I = 1,NGRID

         DO K=1,6                                          ! Make KGGD_II 6x6 even though for SPOINT's we only use 1-1 term
            DO L=1,6
               KGGD_II = ZERO
            ENDDO
         ENDDO 

!xx      CALL CALC_TDOF_ROW_NUM ( GRID_ID(INV_GRID_SEQ(I)), IROW_START, 'N' )
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(INV_GRID_SEQ(I)), IGRID )
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         KGGD_COL_NUM = TDOF(ROW_NUM_START,G_SET_COL)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
k_do:    DO K=1,NUM_COMPS

            KGGD_ROW_NUM = KGGD_ROW_NUM + 1
            IS = STFKEY(KGGD_ROW_NUM)

            IF (IS == 0) THEN                              ! Check for null row in stiffness matrix
               I_KGGD(KGGD_ROW_NUM+1) = I_KGGD(KGGD_ROW_NUM)
               CYCLE k_do
            ENDIF

            NUM_NONZERO_IN_ROW = 0                         ! Form row of non-zero's in arrays RJ, RSTF
j_do1:      DO J=1,NDOFG
               IF (DABS(STF3(IS)%Col_3) >= EPS1) THEN
                  NUM_NONZERO_IN_ROW = NUM_NONZERO_IN_ROW + 1
                  RSTF(NUM_NONZERO_IN_ROW) = STF3(IS)%Col_3
                  RJ(NUM_NONZERO_IN_ROW)   = STF3(IS)%Col_1
               ENDIF
               IS = STF3(IS)%Col_2
               IF (IS == 0) THEN
                  EXIT j_do1
               ENDIF
            ENDDO j_do1

            IF (NUM_NONZERO_IN_ROW > NUM_MAX) THEN
               NUM_MAX = NUM_NONZERO_IN_ROW
            ENDIF   
            IF (IS /= 0) THEN
               WRITE(ERR,1625) SUBR_NAME,I
               WRITE(F06,1625) SUBR_NAME,I
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                       ! Coding error, so quit
            ENDIF
 
            IF (NUM_NONZERO_IN_ROW /= 1) THEN               ! Sort row by the shell method so that RJ is in numerical order
               CALL SORT_INT1_REAL1 ( SUBR_NAME, 'RJ, RSTF', NUM_NONZERO_IN_ROW, RJ, RSTF )
            ENDIF   


n_do:       DO N=1,NUM_NONZERO_IN_ROW                      ! Formulate the K-th row of KGGD_II
               IF ((RJ(N) >= KGGD_COL_NUM) .AND. (RJ(N) <= KGGD_COL_NUM+NUM_COMPS-1)) THEN
                  KGGD_II_COL_NUM = RJ(N) - (KGGD_COL_NUM - 1)
                  KGGD_II(K,KGGD_II_COL_NUM) = RSTF(N)
               ENDIF
            ENDDO n_do
            
j_do3:      DO J=1,NUM_NONZERO_IN_ROW
               KTERM_KGGD = KTERM_KGGD + 1                    ! KTERM_KGGD is a count on the no. records written
               J_KGGD(KTERM_KGGD) = RJ(J)
                 KGGD(KTERM_KGGD) = RSTF(J)
            ENDDO j_do3

            I_KGGD(KGGD_ROW_NUM+1) = I_KGGD(KGGD_ROW_NUM) + NUM_NONZERO_IN_ROW

         ENDDO k_do

         DO K=1,6                                           ! Set lower portion of KGGD_II to be symmetric
            DO J=1,K-1
               KGGD_II(K,J) = KGGD_II(J,K)
            ENDDO
         ENDDO 

         CALL COUNTER_PROGRESS(I)
      ENDDO i_do

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
