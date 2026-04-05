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

      SUBROUTINE SPARSE_MGG

! Add sparse arrays for concentrated masses (MGGC), scalar masses (MGGS) and hash-assembled element mass (MGGE)
! to get the final sparse G-set mass matrix, MGG.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE, DBL_LONG
      USE IOUNT1, ONLY                :  ERR, F04, F06, L1R, L1R_MSG, LINK1R, SC1, WRT_ERR, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCMASS, NDOFG, NGRID, NTERM_MGG, NTERM_MGGC, NTERM_MGGE,         &
                                         NTERM_MGGS, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  SPARSE_MGG_BEGEND
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES,ONLY             :  TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID
      USE PARAMS, ONLY                :  EPSIL, PRTMASS, SUPINFO, WTMASS
      USE EMS_ARRAYS, ONLY            :  EMS_ROW_HM, EMS_COL_HM, EMS_VAL_HM
      USE SPARSE_MATRICES, ONLY       :  I2_MGG, I_MGG, J_MGG, MGG, I_MGGC, J_MGGC, MGGC, I_MGGE, J_MGGE, MGGE,                    &
                                         I_MGGS, J_MGGS, MGGS,  SYM_MGGC, SYM_MGGE, SYM_MGGS
      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1
      USE HOTSPOT_PROFILER, ONLY      :  HOTSPOT_COUNTER_ADD, HOTSPOT_TIMER_ADD, HOTSPOT_TIMER_BEGIN, HOTSPOT_TIMER_END,         &
                                         HOTSPOT_VALUE_ADD, HOTSPOT_WALL_TIME

      USE SPARSE_MGG_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_MGG'
      CHARACTER(  1*BYTE)             :: FOUND
      CHARACTER(LEN=LEN(SYM_MGGE))    :: SYM_CRS1

      INTEGER(LONG)                   :: GRID_NUM
      INTEGER(LONG)                   :: I, J, K
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: IGRID
      INTEGER(LONG)                   :: IK
      INTEGER(LONG)                   :: KTERM_MGGE
      INTEGER(LONG)                   :: MAX_NUM_IN_ROW
      INTEGER(LONG)                   :: NTERM_CRS1
      INTEGER(LONG)                   :: NUM_IN_ROW_I
      INTEGER(LONG)                   :: NUM_COMPS
      INTEGER(LONG)                   :: NZERO
      INTEGER(LONG)                   :: OUNT(2)
      INTEGER(LONG)                   :: POS
      INTEGER(LONG)                   :: RAW_NTERM
      INTEGER(LONG)                   :: ROW_END
      INTEGER(LONG)                   :: ROW_START
      INTEGER(LONG)                   :: ROW_NUM_START
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SPARSE_MGG_BEGEND
      INTEGER(LONG)                   :: HS_SLOT
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: HS_T0
      REAL(DOUBLE)                    :: HS_PHASE_T0
      REAL(DOUBLE)                    :: GRID_MGG(6,6)

      INTRINSIC                       :: DABS

! **********************************************************************************************************************************
      CALL HOTSPOT_TIMER_BEGIN ( 'SPARSE_MGG', HS_SLOT, HS_T0 )

      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
      RAW_NTERM = NTERM_MGGE

! Compact out exact-zero terms after duplicate accumulation in the hash-backed build.

      NZERO = 0
      POS = 0
      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      DO I=1,RAW_NTERM
         IF (DABS(EMS_VAL_HM(I)) < EPS1) THEN
            NZERO = NZERO + 1
            CALL HOTSPOT_COUNTER_ADD ( 'MGGE_ZERO_DROPS', INT(1,DBL_LONG) )
         ELSE
            POS = POS + 1
            EMS_ROW_HM(POS) = EMS_ROW_HM(I)
            EMS_COL_HM(POS) = EMS_COL_HM(I)
            EMS_VAL_HM(POS) = EMS_VAL_HM(I)
         ENDIF
      ENDDO
      NTERM_MGGE = POS
      CALL HOTSPOT_TIMER_ADD ( 'SPARSE_MGG/ZERO_STRIP', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

      WRITE(ERR,146) NTERM_MGGE
      IF (SUPINFO == 'N') THEN
         WRITE(F06,146) NTERM_MGGE
      ENDIF

! Open L1R to write the final mass matrix later in the routine.

      OUNT(1) = ERR
      OUNT(2) = F06
      CALL FILE_OPEN ( L1R, LINK1R, OUNT, 'REPLACE', L1R_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N', 'Y' )

! Build MGGE in CRS directly from the sorted unique triplets.

      KTERM_MGGE = 0
      I_MGGE(1) = 1
      WRITE(SC1, * )
      CALL COUNTER_INIT('     Working on grid ', NGRID)
      HS_LOOP_T0 = HOTSPOT_WALL_TIME()
i_do: DO I = 1,NGRID

         GRID_NUM = GRID_ID(I)
         IGRID = I
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
k_do:    DO K=1,NUM_COMPS

            IK = ROW_NUM_START + K - 1
            IS = EMSKEY(IK)

            IF (IS == 0) THEN                              ! Check for null row in mass matrix
               I_MGGE(IK+1) = I_MGGE(IK)
               CYCLE k_do
            ENDIF
         ENDDO
         CALL HOTSPOT_TIMER_ADD ( 'SPARSE_MGG/ROW_SORT', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
      ENDIF

      POS = 1
      DO I=1,NDOFG
         ROW_START = POS
         DO WHILE ((POS <= NTERM_MGGE) .AND. (EMS_ROW_HM(POS) == I))
            POS = POS + 1
         ENDDO
         ROW_END = POS - 1
         NUM_IN_ROW_I = MAX(0_LONG, ROW_END - ROW_START + 1)
         CALL HOTSPOT_VALUE_ADD ( 'MGGE_ROWLEN_FINAL', DBLE(NUM_IN_ROW_I) )
         DO J=ROW_START,ROW_END
            KTERM_MGGE = KTERM_MGGE + 1
            IF (KTERM_MGGE > NTERM_MGGE) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERM_MGGE, 'MGGE' )
            J_MGGE(KTERM_MGGE) = EMS_COL_HM(J)
              MGGE(KTERM_MGGE) = EMS_VAL_HM(J)
         ENDDO
         I_MGGE(I+1) = I_MGGE(I) + NUM_IN_ROW_I
      ENDDO

      IF (KTERM_MGGE /= NTERM_MGGE) THEN
         WRITE(ERR,1614) SUBR_NAME,LINK1R,KTERM_MGGE,NTERM_MGGE
         WRITE(F06,1614) SUBR_NAME,LINK1R,KTERM_MGGE,NTERM_MGGE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (ALLOCATED(EMS_ROW_HM)) DEALLOCATE ( EMS_ROW_HM )
      IF (ALLOCATED(EMS_COL_HM)) DEALLOCATE ( EMS_COL_HM )
      IF (ALLOCATED(EMS_VAL_HM)) DEALLOCATE ( EMS_VAL_HM )

! *********************************************************************************************************************************
! Call subr to calc MGGS matrix of scalar masses

      IF (NCMASS > 0) THEN
         CALL MGGS_MASS_MATRIX
      ENDIF

! Add MGGC, MGGE and MGGS to get MGG. This is done in 2 steps: add MGGC and MGGE to get temporary CRS1 then add CRS1 to MGGS

!  (1) add MGGC and MGGE to get CRS1 (do not mult by WTMASS here)
!  --------------------------------------------------------------

      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      CALL MATADD_SSS_NTERM ( NDOFG, 'MGGC', NTERM_MGGC, I_MGGC, J_MGGC, SYM_MGGC, 'MGGE', NTERM_MGGE, I_MGGE, J_MGGE, SYM_MGGE,&
                                     'CRS1' , NTERM_CRS1 )

      CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFG, NTERM_CRS1, SUBR_NAME )

      CALL MATADD_SSS ( NDOFG, 'MGGC', NTERM_MGGC, I_MGGC, J_MGGC, MGGC, ONE, 'MGGE', NTERM_MGGE, I_MGGE, J_MGGE, MGGE,            &
                        ONE, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )

!  (2) add CRS1 = MGGC + MGGE and MGGS to get CMGG (mult by WTMASS here)
!  ---------------------------------------------------------------------

      IF (NTERM_MGGS > 0) THEN

         SYM_CRS1 = SYM_MGGS
         CALL MATADD_SSS_NTERM ( NDOFG, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, SYM_CRS1, 'MGGS', NTERM_MGGS, I_MGGS, J_MGGS, SYM_MGGS,&
                                        'MGG' , NTERM_MGG )

         CALL ALLOCATE_L1_MGG ( 'I2_MGG', SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'MGG', NDOFG, NDOFG, SUBR_NAME )

         CALL MATADD_SSS ( NDOFG, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, WTMASS, 'MGGS', NTERM_MGGS, I_MGGS, J_MGGS, MGGS,      &
                           WTMASS, 'MGG', NTERM_MGG, I_MGG, J_MGG, MGG )

      ELSE

         NTERM_MGG = NTERM_CRS1
         CALL ALLOCATE_L1_MGG ( 'I2_MGG', SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'MGG', NDOFG, NTERM_MGG, SUBR_NAME )
         DO I=1,NDOFG+1
            I_MGG(I) = I_CRS1(I)
         ENDDO
         DO I=1,NTERM_MGG
            J_MGG(I) = J_CRS1(I)
              MGG(I) = WTMASS*CRS1(I)
         ENDDO

      ENDIF
      CALL HOTSPOT_TIMER_ADD ( 'SPARSE_MGG/MATADD_MERGES', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

! Deallocate CRS1

      CALL DEALLOCATE_SCR_MAT ( 'CRS1' )

      IF (PRTMASS(1) >= 2) THEN
         IF (NTERM_MGG > 0) THEN
            IF (ALLOCATED(MGGC)) THEN
               CALL WRITE_SPARSE_CRS (  'Conc mass matrix, MGGC', 'G ', 'G ', NTERM_MGGC, NDOFG, I_MGGC, J_MGGC, MGGC )
            ENDIF
            IF (ALLOCATED(MGGE)) THEN
               CALL WRITE_SPARSE_CRS (  'Elem mass matrix, MGGE', 'G ', 'G ', NTERM_MGGE, NDOFG, I_MGGE, J_MGGE, MGGE )
            ENDIF
            IF (ALLOCATED(MGGS)) THEN
               CALL WRITE_SPARSE_CRS ('Scalar mass matrix, MGGS', 'G ', 'G ', NTERM_MGGS, NDOFG, I_MGGS, J_MGGS, MGGS )
            ENDIF
         ENDIF
      ENDIF

! *********************************************************************************************************************************
! Write row, col, value to L1R for matrix MGG

      K = 0
      WRITE(L1R) NTERM_MGG
      IF (NTERM_MGG > 0) THEN
         DO I=1,NDOFG
            NUM_IN_ROW_I = I_MGG(I+1) - I_MGG(I)
            CALL HOTSPOT_VALUE_ADD ( 'MGG_ROWLEN_FINAL', DBLE(NUM_IN_ROW_I) )
            DO J=1,NUM_IN_ROW_I
               K = K + 1
               IF (K > NTERM_MGG)  CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, K, 'MGG' )
               I2_MGG(K) = I
               WRITE(L1R) I2_MGG(K), J_MGG(K), MGG(K)
            ENDDO
         ENDDO
      ENDIF

      CALL FILE_CLOSE ( L1R, LINK1R, 'KEEP', 'Y' )

! Get stats on MGG to write to F06

      IF (NTERM_MGG > 0) THEN

         MAX_NUM_IN_ROW = 0
         DO I=1,NDOFG
            IK = I_MGG(I+1) - I_MGG(I)
            IF (IK > MAX_NUM_IN_ROW) THEN
               MAX_NUM_IN_ROW = IK
            ENDIF
         ENDDO

         WRITE(ERR,147) NTERM_MGG
         WRITE(ERR,101) MAX_NUM_IN_ROW
         IF (SUPINFO == 'N') THEN
            WRITE(F06,147) NTERM_MGG
            WRITE(F06,101) MAX_NUM_IN_ROW
         ENDIF

      ENDIF

! Debug output (print grid 6x6 mass for every grid)

      IERR = 0
      IF (DEBUG(36) > 0) THEN
         WRITE(F06,1101)
         DO K=1,NGRID
            CALL GET_GRID_NUM_COMPS ( K, NUM_COMPS, SUBR_NAME )
            IF (NUM_COMPS == 6) THEN                       ! Only do output for actual grids, not SPOINT's
               IGRID = K
               CALL GET_GRID_6X6_MASS (  GRID_ID(K), IGRID, FOUND, GRID_MGG )
               WRITE(F06,1102) GRID_ID(K)
               DO I=1,3
                  WRITE(F06,1103) (GRID_MGG(I,J),J=1,6)
               ENDDO
               WRITE(F06,*)
               DO I=4,6
                  WRITE(F06,1103) (GRID_MGG(I,J),J=1,6)
               ENDDO
               WRITE(F06,*)
            ENDIF
         ENDDO
         WRITE(F06,1104)
      ENDIF

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      CALL HOTSPOT_TIMER_END ( HS_SLOT, HS_T0 )

      RETURN

! **********************************************************************************************************************************
  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE MGGE MASS MATRIX (ELEMS) IS             = ',I12,/)

  147 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE MGG MASS MATRIX (ELEMS + CONM) IS       = ',I12,/)

  101 FORMAT(' *INFORMATION: MAX NUMBER OF NONZERO TERMS IN A ROW OF THE G-SET MASS MATRIX          = ',I12,/)

 1101 FORMAT(' ___________________________________________________________________________________________________________________'&
            ,'________________'                                                                                                ,//,&
             ' ::::::::::::::::::::::::::::::::::::::::START DEBUG(36) OUTPUT FROM SUBROUTINE SPARSE_MGG:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1102 FORMAT('6 x 6 mass matrix for grid ',I8,/,'-----------------------------------')

 1103 FORMAT(3(1ES14.6),2X,3(1ES14.6))

 1104 FORMAT(' :::::::::::::::::::::::::::::::::::::::::END DEBUG(36) OUTPUT FROM SUBROUTINE SPARSE_MGG::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' ___________________________________________________________________________________________________________________'&
            ,'________________',/)

 1614 FORMAT(' *ERROR  1614: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF G-SET MASS MATRIX RECORDS WRITTEN TO FILE:'                                             &
                    ,/,15X,A                                                                                                       &
                    ,/,14X,' WAS KTERM_MGGE = ',I12,'. IT SHOULD HAVE BEEN NTERM_MGGE = ',I12)

      END SUBROUTINE SPARSE_MGG
