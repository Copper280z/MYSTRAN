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
 
      SUBROUTINE EMP
 
! Element mass processor
 
! EMP generates the portion of the G-set mass matrix due to element mass and stores unique nonzero terms in
! hash-backed row/col/value arrays for later conversion to CRS in SPARSE_MGG.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE, DBL_LONG
      USE ISO_FORTRAN_ENV, ONLY       :  INT64
      USE IOUNT1, ONLY                :  ERR, F04, F06, F22, F22FIL, F22_MSG, SC1, WRT_BUG, WRT_ERR, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_ME_BIT, ELDT_F22_ME_BIT, FATAL_ERR, IBIT, LINKNO, LTERM_MGGE,   &
                                         MBUG, MELDOF, NDOFG, NELE, NGRID, NTERM_MGGE, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR
      USE SUBR_BEGEND_LEVELS, ONLY    :  EMP_BEGEND
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  BGRID, ELDT, ELDOF, ELGP, GRID, NUM_EMG_FATAL_ERRS, ME, OELDT, PLY_NUM, TYPE
      USE EMS_ARRAYS, ONLY            :  EMS_ROW_HM, EMS_COL_HM, EMS_VAL_HM
      USE HOTSPOT_PROFILER, ONLY      :  HOTSPOT_COUNTER_ADD, HOTSPOT_TIMER_ADD, HOTSPOT_TIMER_BEGIN,                           &
                                         HOTSPOT_TIMER_END, HOTSPOT_WALL_TIME
      USE MATRIX_ASSEMBLY_FFHASH, ONLY:  FFH_T
 
      USE EMP_USE_IFs

      IMPLICIT NONE
 
      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EMP'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)
 
      INTEGER(LONG)                   :: EDOF(MELDOF)      ! A list of the G-set DOF's for an elem
      INTEGER(LONG)                   :: EDOF_ROW_NUM      ! Row number in array EDOF
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept 
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Intermediate variable resulting from an IAND operation
      INTEGER(LONG)                   :: IDUM              ! Dummy variable used when flipping DOF's
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: KSTART            ! Used in deciding whether to process all elem mass terms or only
!                                                            the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: MGG_ROW           ! A row no. in MGG
      INTEGER(LONG)                   :: MGG_ROWJ          ! Another row no. in EMS
      INTEGER(LONG)                   :: MGG_COL           ! A col no. in MGG
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR  
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW_NUM      ! Row number in array TDOF
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: ALLOC_TERMS
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = EMP_BEGEND
      INTEGER(LONG)                   :: HS_SLOT
      TYPE(FFH_T)                     :: ENTRY_MAP
 
      REAL(DOUBLE)                    :: DQE(MELDOF,NSUB)  ! Dummy array in call to ELEM_TRANSFORM_LBG
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: HS_T0
      REAL(DOUBLE)                    :: HS_PHASE_T0
      CHARACTER(LEN=128)              :: HS_PHASE_NAME
 
      INTRINSIC                       :: DABS, IAND

! **********************************************************************************************************************************
      CALL HOTSPOT_TIMER_BEGIN ( 'EMP', HS_SLOT, HS_T0 )

      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
      NTERM_MGGE = 0
      ALLOC_TERMS = MAX(1_LONG, LTERM_MGGE)

      IF (ALLOCATED(EMS_ROW_HM)) DEALLOCATE ( EMS_ROW_HM )
      IF (ALLOCATED(EMS_COL_HM)) DEALLOCATE ( EMS_COL_HM )
      IF (ALLOCATED(EMS_VAL_HM)) DEALLOCATE ( EMS_VAL_HM )
      ALLOCATE ( EMS_ROW_HM(ALLOC_TERMS), EMS_COL_HM(ALLOC_TERMS), EMS_VAL_HM(ALLOC_TERMS) )
      EMS_ROW_HM = 0
      EMS_COL_HM = 0
      EMS_VAL_HM = ZERO
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

!! Null dummy array DQE used in call to ELEM_TRANSFORM_LBG

      DO I=1,MELDOF
         DO J=1,NSUB
            DQE(I,J) = ZERO
         ENDDO
      ENDDO

! Set up the option flags for EMG:
 
      OPT(1) = 'Y'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
 
! Process the elements:
 
      IERROR = 0
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages         
      CALL COUNTER_INIT('     Calculating mass matrix. Process elem   ', NELE)
      elems:DO I=1,NELE

         DO J=0,MBUG-1
            WRT_BUG(J) = 0
         ENDDO

         IF (LINKNO == 1) THEN                             ! Only want element mass matrix, ME, written to BUG file in LINK 1
            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_ME_BIT))       ! WRT_BUG(3): printed output of ME
            IF (I1 > 0) THEN
               WRT_BUG(3) = 1
            ENDIF
         ENDIF

         IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
            WRITE(F06,14001)
         ENDIF

         NUM_EMG_FATAL_ERRS = 0
         PLY_NUM = 0
         HS_PHASE_T0 = HOTSPOT_WALL_TIME()
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'Y' )       ! 'Y' means write to BUG file
         HS_PHASE_NAME = 'EMP/EMG/' // TRIM(TYPE)
         CALL HOTSPOT_TIMER_ADD ( HS_PHASE_NAME, HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE elems
         ENDIF 

         I1 = IAND(OELDT,IBIT(ELDT_F22_ME_BIT))            ! Do we need to write elem mass matrices to F22 files
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 2, 0 )
         ENDIF

         EDOF_ROW_NUM = 0                                  ! Generate element DOF'S
         DO J = 1,ELGP
            CALL HOTSPOT_COUNTER_ADD ( 'GRID_LOOKUP/REPLACEABLE/TOTAL', INT(1,DBL_LONG) )
            CALL HOTSPOT_COUNTER_ADD ( 'GRID_LOOKUP/REPLACEABLE/EMP'  , INT(1,DBL_LONG) )
            IGRID = BGRID(J)
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, SUBR_NAME )
            DO K = 1,NUM_COMPS
               TDOF_ROW_NUM       = ROW_NUM_START + K - 1
               EDOF_ROW_NUM       = EDOF_ROW_NUM + 1
               EDOF(EDOF_ROW_NUM) = TDOF(TDOF_ROW_NUM, G_SET_COL_NUM)
            ENDDO 
         ENDDO
 
! Transform ME from local at the elem ends to basic at elem ends to global at elem ends to global at grids.

                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN

            HS_PHASE_T0 = HOTSPOT_WALL_TIME()
            CALL ELEM_TRANSFORM_LBG ( 'ME', ME, DQE )
            CALL HOTSPOT_TIMER_ADD ( 'EMP/ELEM_TRANSFORM_LBG', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
24357 format(6(1es14.6))

         ENDIF 

! Put the element mass matrix, ME, into EMS array. J ranges over rows, K over cols of elem mass matrix, ME 
 
         HS_PHASE_T0 = HOTSPOT_WALL_TIME()
mgg_rows:DO J = 1,ELDOF
            MGG_ROWJ  = EDOF(J)
            IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
               WRITE(F06,*)
            ENDIF

            IF (SPARSTOR == 'SYM') THEN                    ! Set KSTART depending on SPARSTOR
               KSTART = J                                  ! Process only upper right portion of ME
            ELSE
               KSTART = 1                                  ! Process all of ME
            ENDIF

mgg_cols:   DO K = KSTART,ELDOF
               MGG_ROW  = MGG_ROWJ                         ! Make sure we have correct row num. It may have been flipped w/ col
               MGG_COL  = EDOF(K)
               IF (DABS(ME(J,K)) < EPS1) THEN
                  CYCLE mgg_cols
               ENDIF
 
               IF (SPARSTOR == 'SYM') THEN                 ! If 'SYM', Flip MGG_COL,MGG_ROW if MGG_COL < MGG_ROW
                  IF (MGG_COL < MGG_ROW) THEN
                     IDUM    = MGG_ROW
                     MGG_ROW = MGG_COL
                     MGG_COL = IDUM
                  ENDIF
               ENDIF

               CALL ADD_HASH_MASS_TERM ( MGG_ROW, MGG_COL, ME(J,K) )

            ENDDO mgg_cols 

         ENDDO mgg_rows 
         CALL HOTSPOT_TIMER_ADD ( 'EMP/MASS_INSERT', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         CALL COUNTER_PROGRESS(I)

      ENDDO elems 

      WRITE(SC1,*) CR13

! Debug output:

      CALL ENTRY_MAP%RESET()
  
      IF((DEBUG(10) == 21) .OR. (DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR.                                                     &
         (DEBUG(10) == 31) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
         WRITE(F06,1260)
         DO I=1,NTERM_MGGE
            WRITE(F06,1261) I, EMS_ROW_HM(I), EMS_COL_HM(I), EMS_VAL_HM(I)
         ENDDO
         WRITE(F06,*)
      ENDIF
 
! Reset subr EMG option flags:
 
      OPT(3) = 'N'
      OPT(4) = 'N'
 
! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                                    ! IERROR is count of all subr EMG errors, so quit
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
 1260 FORMAT(/,'            I   EMS_ROW_HM   EMS_COL_HM        EMS_VAL_HM')

 1261 FORMAT(1X,I12,I12,I12,3X,1ES21.14)

 1624 FORMAT(' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY NON-ZERO TERMS IN THE ',A,' MATRIX. LIMIT IS ',A,' = ',I12)

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

14001 FORMAT(' ******************************************************************************************************************&
&******')

! **********************************************************************************************************************************
 
! ##################################################################################################################################
 
      CONTAINS
 
      SUBROUTINE ADD_HASH_MASS_TERM ( ROW_NUM, COL_NUM, VALUE )

      INTEGER(LONG), INTENT(IN)       :: ROW_NUM
      INTEGER(LONG), INTENT(IN)       :: COL_NUM
      REAL(DOUBLE) , INTENT(IN)       :: VALUE

      INTEGER(LONG)                   :: IDX
      INTEGER(INT64)                  :: KEY

      KEY = PACK_MATRIX_KEY(ROW_NUM, COL_NUM)
      IDX = ENTRY_MAP%GET_INDEX(KEY)

      IF (IDX >= 0) THEN
         CALL HOTSPOT_COUNTER_ADD ( 'MGGE_DUPLICATE_INSERTS', INT(1,DBL_LONG) )
         EMS_VAL_HM(ENTRY_MAP%VALS(IDX)) = EMS_VAL_HM(ENTRY_MAP%VALS(IDX)) + VALUE
      ELSE
         NTERM_MGGE = NTERM_MGGE + 1
         IF (TYPE == 'QUAD4   ') THEN
            CALL HOTSPOT_COUNTER_ADD ( 'MGGE_NEW_TERMS/CQUAD4', INT(1,DBL_LONG) )
         ELSE
            CALL HOTSPOT_COUNTER_ADD ( 'MGGE_NEW_TERMS/OTHER' , INT(1,DBL_LONG) )
         ENDIF
         IF (NTERM_MGGE > LTERM_MGGE) THEN
            WRITE(ERR,'(A,A,/,A,A,A,I12)') &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME), &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'MASS', ' MATRIX. LIMIT IS LTERM_MGGE = ', LTERM_MGGE
            WRITE(F06,'(A,A,/,A,A,A,I12)') &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME), &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'MASS', ' MATRIX. LIMIT IS LTERM_MGGE = ', LTERM_MGGE
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         EMS_ROW_HM(NTERM_MGGE) = ROW_NUM
         EMS_COL_HM(NTERM_MGGE) = COL_NUM
         EMS_VAL_HM(NTERM_MGGE) = VALUE
         CALL ENTRY_MAP%STORE_VALUE ( KEY, NTERM_MGGE, IDX )
         IF (IDX < 0) THEN
            WRITE(ERR,'(A,A,/,A,A,A,I12)') &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME), &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'MASS', ' MATRIX. LIMIT IS HASH INDEX = ', IDX
            WRITE(F06,'(A,A,/,A,A,A,I12)') &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME), &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'MASS', ' MATRIX. LIMIT IS HASH INDEX = ', IDX
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

      END SUBROUTINE ADD_HASH_MASS_TERM

! ##################################################################################################################################

      PURE INTEGER(INT64) FUNCTION PACK_MATRIX_KEY ( ROW_NUM, COL_NUM )

      INTEGER(LONG), INTENT(IN)       :: ROW_NUM
      INTEGER(LONG), INTENT(IN)       :: COL_NUM

      PACK_MATRIX_KEY = IOR( SHIFTL(INT(ROW_NUM,INT64), 32), INT(COL_NUM,INT64) )

      END FUNCTION PACK_MATRIX_KEY

      END SUBROUTINE EMP
