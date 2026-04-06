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

      SUBROUTINE SPARSE_KGG

! (1) Converts the hash-assembled system KGG matrix into CRS format and writes it to LINK1L.
! (2) Calls KGG_SINGULARITY_PROC on each grid's 6x6 diagonal partition.
! (3) Regenerates TDOF/TDOFI if AUTOSPC changed the TSET table.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE, DBL_LONG
      USE IOUNT1, ONLY                :  ERR, F04, F06, L1L, L1L_MSG, LINK1L, SC1, SPCFIL, SPC, WRT_ERR, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NGRID, NIND_GRDS_MPCS,                                    &
                                         NTERM_KGG, NUM_PCHD_SPC1, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  SPARSE_KGG_BEGEND
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_RAT, EPSIL, PRTTSET, PRTSTIFF, SPARSTOR, SPC1QUIT, SUPINFO, SUPWARN
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, GRID_SEQ, MPC_IND_GRIDS, INV_GRID_SEQ
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START, TDOFI, TSET
      USE STF_ARRAYS, ONLY            :  STF_ROW_HM, STF_COL_HM, STF_VAL_HM
      USE SPARSE_MATRICES, ONLY       :  I_KGG, J_KGG, KGG, SYM_KGG
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE HOTSPOT_PROFILER, ONLY      :  HOTSPOT_COUNTER_ADD, HOTSPOT_TIMER_ADD, HOTSPOT_TIMER_BEGIN, HOTSPOT_TIMER_END,         &
                                         HOTSPOT_VALUE_ADD, HOTSPOT_WALL_TIME

      USE SPARSE_KGG_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_KGG'
      CHARACTER(  7*BYTE)             :: ASPC_SUM_MSG1
      CHARACTER(100*BYTE)             :: ASPC_SUM_MSG2
      CHARACTER( 13*BYTE)             :: ASPC_SUM_MSG3
      CHARACTER(132*BYTE)             :: TDOF_MSG
      CHARACTER( 1*BYTE)              :: SKIPIT = 'N'

      INTEGER(LONG)                   :: AGRIDI
      INTEGER(LONG)                   :: G_SET_COL
      INTEGER(LONG)                   :: I, J, K, N
      INTEGER(LONG)                   :: IGRID
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: KGG_COL_NUM
      INTEGER(LONG)                   :: KGG_ROW_NUM
      INTEGER(LONG)                   :: KGG_II_COL_NUM
      INTEGER(LONG)                   :: KGG_NUM_ASPC
      INTEGER(LONG)                   :: KTERM_KGG
      INTEGER(LONG)                   :: NUM_ASPC_BY_COMP(6)
      INTEGER(LONG)                   :: NUM_MAX = 0
      INTEGER(LONG)                   :: NUM_COMPS
      INTEGER(LONG)                   :: NUM_NONZERO_IN_ROW
      INTEGER(LONG)                   :: NZERO
      INTEGER(LONG)                   :: OUNT(2)
      INTEGER(LONG)                   :: POS
      INTEGER(LONG)                   :: RAW_NTERM
      INTEGER(LONG), ALLOCATABLE      :: ROW_NEXT(:)
      INTEGER(LONG)                   :: ROW_END
      INTEGER(LONG)                   :: ROW_NUM_START
      INTEGER(LONG)                   :: ROW_START
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SPARSE_KGG_BEGEND
      INTEGER(LONG)                   :: HS_SLOT

      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: HS_T0
      REAL(DOUBLE)                    :: HS_PHASE_T0
      REAL(DOUBLE)                    :: HS_LOOP_T0
      REAL(DOUBLE)                    :: KGG_II(6,6)
      INTRINSIC                       :: DABS

! **********************************************************************************************************************************
      CALL HOTSPOT_TIMER_BEGIN ( 'SPARSE_KGG', HS_SLOT, HS_T0 )

      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
      RAW_NTERM = NTERM_KGG
      IF (SPARSTOR == 'SYM') THEN
         SYM_KGG = 'Y'
      ELSE
         SYM_KGG = 'N'
      ENDIF

! Compact out exact-zero terms after duplicate accumulation in the hash-backed build.

      NZERO = 0
      POS = 0
      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      DO I=1,RAW_NTERM
         IF (DABS(STF_VAL_HM(I)) < EPS1) THEN
            NZERO = NZERO + 1
            CALL HOTSPOT_COUNTER_ADD ( 'KGG_ZERO_DROPS', INT(1,DBL_LONG) )
         ELSE
            POS = POS + 1
            STF_ROW_HM(POS) = STF_ROW_HM(I)
            STF_COL_HM(POS) = STF_COL_HM(I)
            STF_VAL_HM(POS) = STF_VAL_HM(I)
         ENDIF
      ENDDO
      NTERM_KGG = POS
      CALL HOTSPOT_TIMER_ADD ( 'SPARSE_KGG/ZERO_STRIP', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

      WRITE(ERR,146) NTERM_KGG
      IF (SUPINFO == 'N') THEN
         WRITE(F06,146) NTERM_KGG
      ENDIF

      IF (NTERM_KGG <= 0) THEN
         WRITE(ERR,1611) NTERM_KGG
         WRITE(F06,1611) NTERM_KGG
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL ALLOCATE_SPARSE_MAT ( 'KGG', NDOFG, NTERM_KGG, SUBR_NAME )

! Build CRS directly from the hash-backed triplets, then sort within each row only.

      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      DO I=1,NDOFG+1
         I_KGG(I) = 0
      ENDDO
      DO I=1,NTERM_KGG
         I_KGG(STF_ROW_HM(I)+1) = I_KGG(STF_ROW_HM(I)+1) + 1
      ENDDO
      I_KGG(1) = 1
      DO I=1,NDOFG
         I_KGG(I+1) = I_KGG(I+1) + I_KGG(I)
      ENDDO
      ALLOCATE ( ROW_NEXT(NDOFG) )
      DO I=1,NDOFG
         ROW_NEXT(I) = I_KGG(I)
      ENDDO
      DO I=1,NTERM_KGG
         POS = ROW_NEXT(STF_ROW_HM(I))
         J_KGG(POS) = STF_COL_HM(I)
           KGG(POS) = STF_VAL_HM(I)
         ROW_NEXT(STF_ROW_HM(I)) = POS + 1
      ENDDO
      CALL HOTSPOT_TIMER_ADD ( 'SPARSE_KGG/CRS_BUILD', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

      HS_PHASE_T0 = HOTSPOT_WALL_TIME()
      DO I=1,NDOFG
         ROW_START = I_KGG(I)
         ROW_END   = I_KGG(I+1) - 1
         NUM_NONZERO_IN_ROW = ROW_END - ROW_START + 1
         IF (NUM_NONZERO_IN_ROW > 1) THEN
            CALL SORT_INT1_REAL1 ( SUBR_NAME, 'KGG row cols', NUM_NONZERO_IN_ROW, J_KGG(ROW_START:ROW_END), KGG(ROW_START:ROW_END) )
         ENDIF
      ENDDO
      CALL HOTSPOT_TIMER_ADD ( 'SPARSE_KGG/ROW_SORT', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
      DEALLOCATE ( ROW_NEXT )

! Open L1L to write stiffness.

      OUNT(1) = ERR
      OUNT(2) = F06
      CALL FILE_OPEN ( L1L, LINK1L, OUNT, 'REPLACE', L1L_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N', 'Y' )
      WRITE(L1L) NTERM_KGG

! Open SPC to write SPC1 records if KGG_SINGULARITY_PROC finds singularities

      OPEN (SPC,FILE=SPCFIL,STATUS='REPLACE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SPCFIL, OUNT, 'Y')
         CALL FILERR ( OUNT, 'Y' )
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      DO I=1,6
         NUM_ASPC_BY_COMP(I) = 0
      ENDDO

      IF (DEBUG(17) > 0) THEN
         WRITE(F06,9901) AUTOSPC_RAT
      ENDIF

! Write the finished CRS arrays and LINK1L records.

      KTERM_KGG = 0
      DO I=1,NDOFG
         ROW_START = I_KGG(I)
         ROW_END   = I_KGG(I+1) - 1
         NUM_NONZERO_IN_ROW = MAX(0_LONG, ROW_END - ROW_START + 1)
         CALL HOTSPOT_VALUE_ADD ( 'KGG_ROWLEN_FINAL', DBLE(NUM_NONZERO_IN_ROW) )
         IF (NUM_NONZERO_IN_ROW > NUM_MAX) NUM_MAX = NUM_NONZERO_IN_ROW
         DO J=ROW_START,ROW_END
            KTERM_KGG = KTERM_KGG + 1
            WRITE(L1L) I, J_KGG(J), KGG(J)
         ENDDO
      ENDDO

! Call singularity processor using the finished CRS matrix.

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
      KGG_ROW_NUM = 0
      CALL COUNTER_INIT('     Working on grid ', NGRID)
      HS_LOOP_T0 = HOTSPOT_WALL_TIME()
i_do: DO I = 1,NGRID
         SKIPIT = 'N'
         KGG_II = ZERO

         IGRID = INV_GRID_SEQ(I)
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         KGG_COL_NUM = TDOF(ROW_NUM_START,G_SET_COL)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         DO K=1,NUM_COMPS

            KGG_ROW_NUM = KGG_ROW_NUM + 1
            ROW_START = I_KGG(KGG_ROW_NUM)
            ROW_END   = I_KGG(KGG_ROW_NUM+1) - 1
            DO N=ROW_START,ROW_END
               IF ((J_KGG(N) >= KGG_COL_NUM) .AND. (J_KGG(N) <= KGG_COL_NUM + NUM_COMPS - 1)) THEN
                  KGG_II_COL_NUM = J_KGG(N) - (KGG_COL_NUM - 1)
                  KGG_II(K,KGG_II_COL_NUM) = KGG(N)
               ENDIF
            ENDDO
         ENDDO

         DO K=1,6
            DO J=1,K-1
               KGG_II(K,J) = KGG_II(J,K)
            ENDDO
         ENDDO

         AGRIDI = GRID_ID(INV_GRID_SEQ(I))
         DO J=1,NIND_GRDS_MPCS
            IF (AGRIDI == MPC_IND_GRIDS(J)) THEN
               SKIPIT = 'Y'
               EXIT
            ENDIF
         ENDDO

         IF (SKIPIT == 'N') THEN
            HS_PHASE_T0 = HOTSPOT_WALL_TIME()
            CALL KGG_SINGULARITY_PROC ( AGRIDI, KGG_II, NUM_ASPC_BY_COMP )
            CALL HOTSPOT_TIMER_ADD ( 'SPARSE_KGG/SINGULARITY_PROC', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         ENDIF
         CALL COUNTER_PROGRESS(I)
      ENDDO i_do
      CALL HOTSPOT_TIMER_ADD ( 'SPARSE_KGG/ROW_EXTRACT', HOTSPOT_WALL_TIME() - HS_LOOP_T0 )

      IF (DEBUG(17) > 0) THEN
         WRITE(F06,9902)
      ENDIF
      WRITE(F06,*)
      WRITE(SC1,*) CR13

      IF (PRTSTIFF(1) >= 1) THEN
         CALL WRITE_SPARSE_CRS ( 'STIFFNESS MATRIX KGG' , 'G ', 'G ', NTERM_KGG, NDOFG, I_KGG, J_KGG, KGG )
      ENDIF

      IF (AUTOSPC == 'Y') THEN
         KGG_NUM_ASPC = 0
         DO I=1,6
            KGG_NUM_ASPC = KGG_NUM_ASPC + NUM_ASPC_BY_COMP(I)
         ENDDO

         IF (KGG_NUM_ASPC > 0) THEN
            IF (PRTTSET > 0) THEN
               WRITE(F06,56)
               WRITE(F06,57)
               DO J = 1,NGRID
                  WRITE(F06,58) GRID(J,1), GRID_SEQ(J), (TSET(J,K),K = 1,6)
               ENDDO
               WRITE(F06,'(//)')
            ENDIF

            ASPC_SUM_MSG1(1:) = 'Stage 1:'
            ASPC_SUM_MSG2(1:) = 'after identification of AUTOSPC''s at the grid level'
            ASPC_SUM_MSG3(1:) = 'in this stage'
            CALL AUTOSPC_SUMMARY_MSGS ( ASPC_SUM_MSG1, ASPC_SUM_MSG2, ASPC_SUM_MSG3, 'Y', NUM_ASPC_BY_COMP )

            TDOF_MSG(1:)  = ' '
            TDOF_MSG(39:) = ASPC_SUM_MSG2(1:)
            CALL TDOF_PROC ( TDOF_MSG )
         ENDIF
      ENDIF

      IF (NUM_PCHD_SPC1 > 0) THEN
         CALL FILE_CLOSE ( SPC, SPCFIL, 'KEEP', 'Y' )
         IF (SPC1QUIT == 'Y') THEN
            WRITE(ERR,9991) SUBR_NAME, SPC1QUIT
            WRITE(F06,9991) SUBR_NAME, SPC1QUIT
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ELSE
         CALL FILE_CLOSE ( SPC, SPCFIL, 'DELETE', 'Y' )
      ENDIF

      IF (KTERM_KGG /= NTERM_KGG) THEN
         WRITE(ERR,1623) SUBR_NAME, LINK1L, KTERM_KGG, NTERM_KGG
         WRITE(F06,1623) SUBR_NAME, LINK1L, KTERM_KGG, NTERM_KGG
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL FILE_CLOSE ( L1L, LINK1L, 'KEEP', 'Y' )
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

      CALL HOTSPOT_TIMER_END ( HS_SLOT, HS_T0 )

      RETURN

! **********************************************************************************************************************************
   56 FORMAT(64X,'DEGREE OF FREEDOM SET TABLE (TSET)')

   57 FORMAT(33x,'     GRID SEQUENCE       T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(33x,2(1X,I8),6(7X,A2))

  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE KGG STIFFNESS MATRIX IS                 = ',I12,/)

  101 FORMAT(' *INFORMATION: MAX NUMBER OF NONZERO TERMS IN A ROW OF THE G-SET STIFFNESS MATRIX     = ',I12,/)

 1611 FORMAT(' *ERROR  1611: THE G-SET STIFFNESS MATRIX, KGG, MUST HAVE SOME NONZERO TERMS. HOWEVER IT HAS ',I12,' TERMS')

 1623 FORMAT(' *ERROR  1623: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF G-SET STIFFNESS MATRIX RECORDS WRITTEN TO FILE:'                                        &
                    ,/,15X,A                                                                                                       &
                    ,/,14X,' WAS KTERM_KGG = ',I12,'. IT SHOULD HAVE BEEN NTERM_KGG = ',I12)

 9901 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(17) OUTPUT FROM SUBROUTINE KGG_SINGULARITY_PROC:::::::::::::::::::',&
             ':::::::::::::::::',//,54X,'AUTOSPC_RAT = ',1ES13.6,/)

 9902 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(17) OUTPUT FROM SUBROUTINE KGG_SINGULARITY_PROC::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

 9991 FORMAT(' PROCESSING ABORTED IN SUBR ',A,' BASED ON PARAMETER SPC1QUIT = ',A)

      END SUBROUTINE SPARSE_KGG
