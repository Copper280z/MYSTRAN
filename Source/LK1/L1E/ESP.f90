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
 
      SUBROUTINE ESP
 
! Element stiffness processor
 
! ESP generates the G-set stiffness matrix and puts it into the 1D array STF of nonzero stiffness terms above the
! diagonal.
 
! ESP processes the elements sequentially to generate element KE matrix using the EMG set of routines. The element
! stiffness are transformed from local to basic to global coords for each grid and then merged into
! hash-backed row/col/value arrays for later conversion to CRS.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE, DBL_LONG
      USE ISO_FORTRAN_ENV, ONLY       :  INT64
      USE IOUNT1, ONLY                :  ERR, F04, F06, F23, F23FIL, F23_MSG, F24, F24FIL, F24_MSG, SC1,                           &
                                         WRT_BUG, WRT_ERR, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_KE_BIT, ELDT_BUG_SE_BIT,                                           &
                                         ELDT_F23_KE_BIT, ELDT_F24_SE_BIT, ELDT_BUG_BCHK_BIT, ELDT_BUG_BMAT_BIT, ELDT_BUG_SHPJ_BIT,&
                                         FATAL_ERR, IBIT, LINKNO, LTERM_KGG, LTERM_KGGD, MBUG, MELDOF, NDOFG, NELE, NGRID,         &
                                         NTERM_KGG, NTERM_KGGD, NSUB, SOL_NAME
      USE PARAMS, ONLY                :  EPSIL
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SUBR_BEGEND_LEVELS, ONLY    :  ESP_BEGEND
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  BGRID, ELDT, ELDOF, ELGP, GRID, NUM_EMG_FATAL_ERRS, PLY_NUM, OELDT, KE, KED, TYPE
      USE STF_ARRAYS, ONLY            :  STF_ROW_HM, STF_COL_HM, STF_VAL_HM
      USE STF_TEMPLATE_ARRAYS, ONLY   :  CROW, TEMPLATE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE HOTSPOT_PROFILER, ONLY      :  HOTSPOT_COUNTER_ADD, HOTSPOT_TIMER_ADD, HOTSPOT_TIMER_BEGIN,                           &
                                         HOTSPOT_TIMER_END, HOTSPOT_WALL_TIME
      USE MATRIX_ASSEMBLY_FFHASH, ONLY:  FFH_T
 
      USE ESP_USE_IFs

      IMPLICIT NONE
 
      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)
      CHARACTER(24*BYTE)              :: NAME              ! Name for output error purposes
      INTEGER(LONG), PARAMETER        :: DEB_NUM   = 46    ! Debug number for output error message
      INTEGER(LONG)                   :: EDOF(MELDOF)      ! A list of the G-set DOF's for an elem
      INTEGER(LONG)                   :: EDOF_ROW_NUM      ! Row number in array EDOF
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept 
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Intermediate variable resulting from an IAND operation
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: KGG_ROW           ! A row no. in KGG or KGGD
      INTEGER(LONG)                   :: KGG_ROWJ          ! Another row no. in KGG or KGGD
      INTEGER(LONG)                   :: KGG_COL           ! A col no. in KGG or KGGD
      INTEGER(LONG)                   :: NTERM             ! Either NTERM_KGGD (BUCKLING) or NTERM_KGG otherwise
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.   
      INTEGER(LONG)                   :: PKTERM            ! Count of the terms in TEMPLATE for nonzero stiffness terms
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW_NUM      ! Row number in array TDOF
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: LTERM             ! Either LTERM_KGGD (BUCKLING) or LTERM_KGG otherwise
      INTEGER(LONG)                   :: ALLOC_TERMS
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = ESP_BEGEND
      INTEGER(LONG)                   :: HS_SLOT
      TYPE(FFH_T)                     :: ENTRY_MAP

      REAL(DOUBLE)                    :: DQE(MELDOF,NSUB)  ! Dummy array in call to ELEM_TRANSFORM_LBG
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: HS_T0
      REAL(DOUBLE)                    :: HS_PHASE_T0
      REAL(DOUBLE)                    :: TERM_VALUE
 
      INTRINSIC                       :: DABS
      INTRINSIC                       :: IAND
      INTRINSIC                       :: MAX

! **********************************************************************************************************************************
      CALL HOTSPOT_TIMER_BEGIN ( 'ESP', HS_SLOT, HS_T0 )

      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      EPS1  = EPSIL(1)
      NTERM = 0

! Make units for writing errors the error file and output file
 
      OUNT(1) = ERR
      OUNT(2) = F06
 
! Null dummy array DQE used in call to ELEM_TRANSFORM_LBG

      DO I=1,MELDOF
         DO J=1,NSUB
            DQE(I,J) = ZERO
         ENDDO
      ENDDO

! LTERM is used in this subr to make sure we do not try to use more array dimension than was allocated. So here, set LTERM to be
! either LTERM_KGG or LTERM_KGGD depending on BUCKLING

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         LTERM = LTERM_KGGD
      ELSE
         LTERM = LTERM_KGG
      ENDIF
      ALLOC_TERMS = MAX(1_LONG, LTERM)

      IF (ALLOCATED(STF_ROW_HM)) DEALLOCATE ( STF_ROW_HM )
      IF (ALLOCATED(STF_COL_HM)) DEALLOCATE ( STF_COL_HM )
      IF (ALLOCATED(STF_VAL_HM)) DEALLOCATE ( STF_VAL_HM )
      ALLOCATE ( STF_ROW_HM(ALLOC_TERMS), STF_COL_HM(ALLOC_TERMS), STF_VAL_HM(ALLOC_TERMS) )
      STF_ROW_HM = 0
      STF_COL_HM = 0
      STF_VAL_HM = ZERO

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )


! DEBUG(10) = 13 or 33 requests that array TEMPLATE be printed

      IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN
         CALL ALLOCATE_TEMPLATE ( SUBR_NAME )
      ENDIF

! Set up the option flags for EMG:

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
 
      IF      ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         OPT(4) = 'N'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'Y'                                      ! OPT(6) is for calc of KE-nonlinear
      ELSE IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         OPT(4) = 'Y'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'Y'                                      ! OPT(6) is for calc of KE-nonlinear
      ELSE
         OPT(4) = 'Y'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'N'                                      ! OPT(6) is for calc of KE-nonlinear
      ENDIF

! Process the elements:
 
      IERROR = 0
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages         
      CALL COUNTER_INIT('     Calculating stiff matrix. Process elem  ', NELE)
      elems:DO I=1,NELE

         IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
            WRITE(F06,14001)
         ENDIF

         DO J=0,MBUG-1
            WRT_BUG(J) = 0
         ENDDO

         IF (LINKNO == 1) THEN                             ! Only want element stiff matrix, KE, written to BUG file in LINK 1

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_KE_BIT))       ! WRT_BUG(4): printed output of KE
            IF (I1 > 0) THEN
               WRT_BUG(4) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_SE_BIT))       ! WRT_BUG(5): printed output of SEi, STEi
            IF (I1 > 0) THEN
               WRT_BUG(5) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_SHPJ_BIT))     ! WRT_BUG(7): printed output of shape fcns and Jacobians for some elems
            IF (I1 > 0) THEN
               WRT_BUG(7) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_BMAT_BIT))     ! WRT_BUG(8): printed output of strain-displ matrices for some elements
            IF (I1 > 0) THEN
               WRT_BUG(8) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_BCHK_BIT))     ! WRT_BUG(9): printed output of R.B., const strain checks for some elems
            IF (I1 > 0) THEN
               WRT_BUG(9) = 1
            ENDIF

         ENDIF

         OPT(3) = 'N'                                      ! OPT(3) is for calc of SEi, STEi
         IF ((WRT_BUG(4) == 1) .OR. (WRT_BUG(5) == 1)) THEN
            OPT(3) = 'Y'
         ENDIF
 
         PLY_NUM = 0
         HS_PHASE_T0 = HOTSPOT_WALL_TIME()
         CALL EMG ( I   , OPT, 'Y', SUBR_NAME, 'Y' )       ! 'N' means do not write to BUG file
         CALL HOTSPOT_TIMER_ADD ( 'ESP/EMG', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )

         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE elems
         ENDIF 

         I1 = IAND(OELDT,IBIT(ELDT_F23_KE_BIT))           ! Do we need to write elem stiff matrices to F23 files?
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 3, 0 )
         ENDIF

         I1 = IAND(OELDT,IBIT(ELDT_F24_SE_BIT))           ! Do we need to write elem stress recovery matrices to F24 files?
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 4, 0 )
         ENDIF

         EDOF_ROW_NUM = 0                                  ! Generate element DOF'S
         DO J = 1,ELGP
            CALL HOTSPOT_COUNTER_ADD ( 'GRID_LOOKUP/REPLACEABLE/TOTAL', INT(1,DBL_LONG) )
            CALL HOTSPOT_COUNTER_ADD ( 'GRID_LOOKUP/REPLACEABLE/ESP'  , INT(1,DBL_LONG) )
            IGRID = BGRID(J)
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, SUBR_NAME )
            DO K = 1,NUM_COMPS
               TDOF_ROW_NUM       = ROW_NUM_START + K - 1
               EDOF_ROW_NUM       = EDOF_ROW_NUM + 1
               EDOF(EDOF_ROW_NUM) = TDOF(TDOF_ROW_NUM, G_SET_COL_NUM)
            ENDDO 
         ENDDO 
 
! Write diagonostics on negative diag stiffness before transformation to global

         IF ((DEBUG(189) == 1) .OR. (DEBUG(189) == 3)) THEN
            CALL WRITE_NEG_DIAG_STIFFNESS ( 1 )
         ENDIF

! Transform KE from local at the elem ends to basic at elem ends to global at elem ends to global at grids.
                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  ')) THEN
            HS_PHASE_T0 = HOTSPOT_WALL_TIME()
            IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
               CALL ELEM_TRANSFORM_LBG ( 'KED', KED, DQE )
            ELSE
               CALL ELEM_TRANSFORM_LBG ( 'KE' , KE , DQE )
            ENDIF
            CALL HOTSPOT_TIMER_ADD ( 'ESP/ELEM_TRANSFORM_LBG', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         ENDIF 

! Write diagonostics on negative diag stiffness after  transformation to global

         IF ((DEBUG(189) == 2) .OR. (DEBUG(189) == 3)) THEN
            CALL WRITE_NEG_DIAG_STIFFNESS ( 2 )
         ENDIF

! Put the elem stiff matrix, KE, or KED (now in global coords), into STF array. J ranges over rows, K over cols of elem stiff mat 
 
         HS_PHASE_T0 = HOTSPOT_WALL_TIME()
kgg_rows:DO J = 1,ELDOF
            KGG_ROWJ  = EDOF(J)
            IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
               WRITE(F06,*)
            ENDIF

kgg_cols:   DO K = 1,ELDOF
               KGG_ROW  = KGG_ROWJ                         ! Make sure we have correct row num. It may have been flipped w/ col
               KGG_COL  = EDOF(K)
               IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
                  TERM_VALUE = KED(J,K)
                  IF (DABS(TERM_VALUE) < EPS1) THEN
                     CYCLE kgg_cols
                  ENDIF
               ELSE
                  TERM_VALUE = KE(J,K)
                  IF (DABS(TERM_VALUE) < EPS1) THEN
                     CYCLE kgg_cols
                  ENDIF
               ENDIF
 
               IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN
                  IF (ALLOCATED(TEMPLATE)) THEN
                     TEMPLATE(KGG_ROW,KGG_COL) = .TRUE.
                  ELSE
                     NAME = 'TEMPLATE                '
                     WRITE(ERR,1628) SUBR_NAME,DEB_NUM,NAME
                     WRITE(F06,1628) SUBR_NAME,DEB_NUM,NAME
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ENDIF
               CALL ADD_HASH_STIFF_TERM ( KGG_ROW, KGG_COL, TERM_VALUE )
 
            ENDDO kgg_cols 

         ENDDO kgg_rows 
         CALL HOTSPOT_TIMER_ADD ( 'ESP/STIFF_INSERT', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         CALL COUNTER_PROGRESS(I)
      ENDDO elems 
      WRITE(SC1,*) CR13
 
! Reset subr EMG option flags:
 
      OPT(3) = 'N'
      OPT(4) = 'N'
      OPT(6) = 'N'
 
! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! IERROR is count of all subr EMG errors, so quit
      ENDIF

! Print out TEMPLATE which shows where the nonzero values are in the upper triangle of the stiffness matrix

      IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN 

         IF (ALLOCATED(TEMPLATE)) THEN

            PKTERM = 0                                     ! Count nonzero terms in K based on TEMPLATE array.
            DO I=1,NDOFG                                   ! Call this PKTERM and it should be same as LTERM
               DO J=I,NDOFG
                  IF (TEMPLATE(I,J)) THEN
                     PKTERM = PKTERM + 1
                  ENDIF
               ENDDO 
            ENDDO

            WRITE(F06,14002) PKTERM 
            WRITE(F06,*)
            DO I=1,NDOFG
               DO J=1,NDOFG
                  CROW(J) = ' '
               ENDDO
               DO J=I,NDOFG
                  IF (TEMPLATE(I,J)) THEN
                     CROW(J) = 'K'
                  ELSE
                     CROW(J) = '_'
                  ENDIF
               ENDDO 
               WRITE(F06,*) (CROW(J),J=1,NDOFG)
            ENDDO
            WRITE(F06,*)

         ELSE

            WRITE(ERR,1628) SUBR_NAME,DEB_NUM,NAME
            WRITE(F06,1628) SUBR_NAME,DEB_NUM,NAME
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Coding error (TEMPLATE should be allocated), so quit
            ENDIF

      ENDIF
   
! Deallocate TEMPLATE and CROW arrays  

      IF ((ALLOCATED(TEMPLATE)) .OR. (ALLOCATED(CROW))) THEN
         CALL DEALLOCATE_TEMPLATE
      ENDIF 

! Reset LTERM and NTERM to appropriate values

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         NTERM_KGGD = NTERM
         LTERM_KGGD = NTERM_KGGD                           ! reset LTERM now that we have det. actual number of terms in KGGD
      ELSE
         NTERM_KGG  = NTERM
         LTERM_KGG  = NTERM_KGG                            ! reset LTERM now that we have det. actual number of terms in KGG
      ENDIF

      CALL ENTRY_MAP%RESET()


! **********************************************************************************************************************************
! Debug output:
  
      IF((DEBUG(10) == 11) .OR. (DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR.                                                     &
         (DEBUG(10) == 31) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
         WRITE(F06,1260)
         DO I=1,NTERM
            WRITE(F06,1261) I,STF_ROW_HM(I),STF_COL_HM(I),STF_VAL_HM(I)
         ENDDO 
         WRITE(F06,*)
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
 1260 FORMAT(/,'            I   STF_ROW_HM   STF_COL_HM        STF_VAL_HM')

 1261 FORMAT(1X,I12,I12,I12,3X,1ES21.14)

 1624 FORMAT(' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY NON-ZERO TERMS IN THE ',A,' MATRIX. LIMIT IS ',A,'    = ',I12)

 1628 FORMAT(' *ERROR  1628: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' BASED ON DEBUG ',I3,' VALUE, ARRAY ',A,' SHOULD BE ALLOCATED BUT IT IS NOT')

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

14001 FORMAT(' ********************************************************************************************************************&
&******')

14002 FORMAT('From subr ESP: TEMPLATE array showing ',I12,' actual nonzero terms in upper triangle of K')

35791 format(' In ESP: I, EID, J, K, BGRID(J), TDOF_ROW_NUM, 6*(BGRID(J)-1)+K, Diff = ',8i8)

88770 format(' In ESP:                         KGG_MAX_DIAG_TERM  = ',47X,1ES15.6)

88771 format(' In ESP #1: J, K, IS, NTERM_KGG, KGG(diagonal term) = ',4i8,1es15.6)

88772 format(' In ESP #2: J, K, IS, NTERM_KGG, KGG(diagonal term) = ',4i8,1es15.6)

88773 format(' In ESP #3: J, K, IS, NTERM_KGG, KGG(diagonal term) = ',4i8,1es15.6)

! **********************************************************************************************************************************
 
! ##################################################################################################################################
 
      CONTAINS
 
! ##################################################################################################################################

      SUBROUTINE ADD_HASH_STIFF_TERM ( ROW_NUM, COL_NUM, VALUE )

      INTEGER(LONG), INTENT(IN)       :: ROW_NUM
      INTEGER(LONG), INTENT(IN)       :: COL_NUM
      REAL(DOUBLE) , INTENT(IN)       :: VALUE

      INTEGER(LONG)                   :: IDX
      INTEGER(INT64)                  :: KEY

      KEY = PACK_MATRIX_KEY(ROW_NUM, COL_NUM)
      IDX = ENTRY_MAP%GET_INDEX(KEY)

      IF (IDX >= 0) THEN
         CALL HOTSPOT_COUNTER_ADD ( 'KGG_DUPLICATE_INSERTS', INT(1,DBL_LONG) )
         STF_VAL_HM(ENTRY_MAP%VALS(IDX)) = STF_VAL_HM(ENTRY_MAP%VALS(IDX)) + VALUE
      ELSE
         NTERM = NTERM + 1
         IF (NTERM > LTERM) THEN
            WRITE(ERR,'(A,A,/,A,A,A,I12)')                                                                                           &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME),                                               &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'STIFFNESS', ' MATRIX. LIMIT IS LTERM = ', LTERM
            WRITE(F06,'(A,A,/,A,A,A,I12)')                                                                                           &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME),                                               &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'STIFFNESS', ' MATRIX. LIMIT IS LTERM = ', LTERM
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         STF_ROW_HM(NTERM) = ROW_NUM
         STF_COL_HM(NTERM) = COL_NUM
         STF_VAL_HM(NTERM) = VALUE
         CALL ENTRY_MAP%STORE_VALUE ( KEY, NTERM, IDX )
         IF (IDX < 0) THEN
            WRITE(ERR,'(A,A,/,A,A,A,I12)')                                                                                           &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME),                                               &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'STIFFNESS', ' MATRIX. LIMIT IS HASH IDX = ', IDX
            WRITE(F06,'(A,A,/,A,A,A,I12)')                                                                                           &
                  ' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ', TRIM(SUBR_NAME),                                               &
                  '              TOO MANY NON-ZERO TERMS IN THE ', 'STIFFNESS', ' MATRIX. LIMIT IS HASH IDX = ', IDX
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

      END SUBROUTINE ADD_HASH_STIFF_TERM

! ##################################################################################################################################

      PURE INTEGER(INT64) FUNCTION PACK_MATRIX_KEY ( ROW_NUM, COL_NUM )

      INTEGER(LONG), INTENT(IN)       :: ROW_NUM
      INTEGER(LONG), INTENT(IN)       :: COL_NUM

      PACK_MATRIX_KEY = IOR( SHIFTL(INT(ROW_NUM,INT64), 32), INT(COL_NUM,INT64) )

      END FUNCTION PACK_MATRIX_KEY

      SUBROUTINE WRITE_NEG_DIAG_STIFFNESS ( WHAT )

      USE CONSTANTS_1,ONLY            :  ZERO
      USE MODEL_STUF, ONLY            :  AGRID, BGRID, EID, ELGP, ELDOF, TYPE

      IMPLICIT NONE

      INTEGER(LONG)                   :: II,JJ,KK,LL,MM    ! DO loop indices or counters
      INTEGER(LONG)                   :: NUM_DIAG_NEGS     ! Number of neg values on the diag of the quad element stiffness matrix
      INTEGER(LONG)                   :: WHAT              ! What header to write out

      REAL(DOUBLE)                    :: MAX_ABS_DIAG      ! Max absolute diagonal term
      REAL(DOUBLE)                    :: RATIO             ! Ratio of diagonal term to MAX_ABS_DIAG
      REAL(DOUBLE)                    :: ZE(ELDOF,ELDOF)   ! Either KE or KED (KED for BUCKLING)

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         DO II=1,ELDOF
            DO JJ=1,ELDOF
               ZE(II,JJ) = KED(II,JJ)
            ENDDO
         ENDDO
      ELSE
         DO II=1,ELDOF
            DO JJ=1,ELDOF
               ZE(II,JJ) = KE(II,JJ)
            ENDDO
         ENDDO
      ENDIF


      MAX_ABS_DIAG = ZERO
      NUM_DIAG_NEGS = 0
      DO II=1,ELDOF
         IF (DABS(ZE(II,II)) > MAX_ABS_DIAG) THEN
            MAX_ABS_DIAG = ZE(II,II)
         ENDIF
         IF (ZE(II,II) < 0.D0) THEN
            NUM_DIAG_NEGS = NUM_DIAG_NEGS + 1               
         ENDIF
      ENDDO

      IF (NUM_DIAG_NEGS > 0) THEN

         IF (WHAT == 1) THEN
            WRITE(F06,97531) TYPE, EID, NUM_DIAG_NEGS
         ELSE IF (WHAT == 2) THEN
            WRITE(F06,97532) TYPE, EID, NUM_DIAG_NEGS
         ENDIF

         WRITE(F06,97533)
         KK=0
         DO LL=1,ELGP
            CALL GET_GRID_NUM_COMPS ( BGRID(LL), NUM_COMPS, SUBR_NAME )
            DO MM=1,NUM_COMPS
               KK = KK + 1
               RATIO = ZERO
               IF (MAX_ABS_DIAG > ZERO) THEN
                  RATIO = ZE(KK,KK)/MAX_ABS_DIAG
               ENDIF
               IF (ZE(KK,KK) >= ZERO) THEN
                  IF (MM == 1) THEN
                     WRITE(F06,90869) AGRID(LL), MM, ZE(KK,KK), RATIO
                  ELSE
                     WRITE(F06,90870) MM, ZE(KK,KK), RATIO
                  ENDIF
               ELSE
                  IF (MM == 1) THEN
                     WRITE(F06,90871) AGRID(LL), MM, ZE(KK,KK), RATIO
                  ELSE
                     WRITE(F06,90872) MM, ZE(KK,KK), RATIO
                  ENDIF
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         WRITE(F06,*)
      ENDIF

97531 FORMAT(' Diagonal stiffnesses for ',A,' element ',I8, ' in local elem coords. Element has ', I3,' negative diagonal terms')

97532 FORMAT(' Diagonal stiffnesses for ',A,' element ',I8, ' in global     coords. Element has ', I3,' negative diagonal terms')

97533 FORMAT( '    Grid  Comp  Diagonal stiffness      Ratio')

90869 FORMAT( I8,I6,1ES17.6,1ES17.6)

90870 FORMAT( 8X,I6,1ES17.6,1ES17.6)

90871 FORMAT( I8,I6,1ES17.6,1ES17.6,' *****NEGATIVE STIFFNESS*****')

90872 FORMAT( 8X,I6,1ES17.6,1ES17.6,' *****NEGATIVE STIFFNESS*****')

      END SUBROUTINE WRITE_NEG_DIAG_STIFFNESS

      END SUBROUTINE ESP
