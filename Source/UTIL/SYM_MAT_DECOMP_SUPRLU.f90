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
 
      SUBROUTINE SYM_MAT_DECOMP_SUPRLU ( CALLING_SUBR, MATIN_NAME, MATIN_SET, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, INFO )

! Decomposes a symmetric band matrix into triangular factors. The input matrix, MATIN, is stored in CRS sparse format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_LOG, ERR, F04, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC       
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  CRS_CCS, SPARSTOR, BAILOUT
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1, I_CRS1, J_CRS1, CRS1
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS
      USE SUBR_BEGEND_LEVELS, ONLY    :  SYM_MAT_DECOMP_SUPRLU_BEGEND
      USE HOTSPOT_PROFILER, ONLY      :  HOTSPOT_TIMER_ADD, HOTSPOT_TIMER_BEGIN, HOTSPOT_TIMER_END, HOTSPOT_WALL_TIME

      USE SYM_MAT_DECOMP_SUPRLU_USE_IFs
                      
      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SYM_MAT_DECOMP_SUPRLU'

      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! The subr that called this subr (used for output error purposes)
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME        ! Name of matrix to be decomposed
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_SET         ! Set designator for the input matrix. If it corresponds to a MYSTRAN
!                                                            displ set (e.g. 'L ' set) then error messages about singulatities
!                                                            can reference the grid/comp that is singular (otherwise the row/col
!                                                            where the singularity occurs is referenced). If it is not a MYSTRAN
!                                                            set designator it should be blank
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in sparse matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of nonzeros in sparse matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)  ! Indicators of number of nonzero terms in rows of matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)   ! Col numberts of nonzero terms in matrix MATIN

      INTEGER(LONG), INTENT(INOUT)    :: INFO              ! Output from SuperLU routine

      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SYM_MAT_DECOMP_SUPRLU_BEGEND
      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: HS_SLOT
      INTEGER(LONG)                   :: NTERM_SUPRLU      ! Number of terms handed to SuperLU
      INTEGER(LONG)                   :: NUM_MATIN_DIAG_0  ! Number of zero diagonal terms in the input matrix
      LOGICAL                         :: MATIN_IS_FULL

      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERMS)
      REAL(DOUBLE)                    :: DUM_COL(NROWS)    ! Temp variable for solving equations
      REAL(DOUBLE)                    :: HS_PHASE_T0
      REAL(DOUBLE)                    :: HS_T0

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

      CALL HOTSPOT_TIMER_BEGIN ( 'SYM_MAT_DECOMP_SUPRLU', HS_SLOT, HS_T0 )

! **********************************************************************************************************************************

      DO I=1,NROWS                                         ! Need a null col of loads when SuperLU is called to factor KLL
         DUM_COL(I) = ZERO                                 ! (only because it appears in the calling list)
      ENDDO

      IF (SPARSTOR == 'SYM   ') THEN

         MATIN_IS_FULL = .FALSE.
         DO I=1,NROWS
            DO J=I_MATIN(I),I_MATIN(I+1)-1
               IF (J_MATIN(J) < I) THEN
                  MATIN_IS_FULL = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (MATIN_IS_FULL) EXIT
         ENDDO

         IF (.NOT.MATIN_IS_FULL) THEN
            HS_PHASE_T0 = HOTSPOT_WALL_TIME()
            CALL SPARSE_MAT_DIAG_ZEROS ( MATIN_NAME, NROWS, NTERMS, I_MATIN, J_MATIN, NUM_MATIN_DIAG_0 )
            NTERM_SUPRLU = 2*NTERMS - (NROWS - NUM_MATIN_DIAG_0)
            IF (ALLOCATED(CRS1)) CALL DEALLOCATE_SCR_MAT ( 'CRS1' )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS, NTERM_SUPRLU, SUBR_NAME )
            CALL CRS_SYM_TO_CRS_NONSYM ( MATIN_NAME, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, 'CRS1', NTERM_SUPRLU,                 &
                                         I_CRS1, J_CRS1, CRS1, 'N' )
            CALL HOTSPOT_TIMER_ADD ( 'SYM_MAT_DECOMP_SUPRLU/EXPAND_SYM_TO_NONSYM', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         ELSE
            NTERM_SUPRLU = NTERMS
         ENDIF

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         MATIN_IS_FULL = .TRUE.
         NTERM_SUPRLU = NTERMS

      ELSE                                              ! Error - incorrect CRS_CCS 

         WRITE(ERR,932) SUBR_NAME, 'SPARSTOR'
         WRITE(F06,932) SUBR_NAME, 'SPARSTOR'
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF      (CRS_CCS == 'CRS') THEN                   ! Use MATIN stored in Compressed Row Storage (CRS) format

         IF (MATIN_IS_FULL) THEN
            HS_PHASE_T0 = HOTSPOT_WALL_TIME()
            CALL C_FORTRAN_DGSSV( 1, NROWS, NTERM_SUPRLU, 1, MATIN, J_MATIN, I_MATIN, DUM_COL, NROWS, SLU_FACTORS, INFO )
            CALL HOTSPOT_TIMER_ADD ( 'SYM_MAT_DECOMP_SUPRLU/SUPERLU_CALL', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         ELSE
            HS_PHASE_T0 = HOTSPOT_WALL_TIME()
            CALL C_FORTRAN_DGSSV( 1, NROWS, NTERM_SUPRLU, 1, CRS1, J_CRS1, I_CRS1, DUM_COL, NROWS, SLU_FACTORS, INFO )
            CALL HOTSPOT_TIMER_ADD ( 'SYM_MAT_DECOMP_SUPRLU/SUPERLU_CALL', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         ENDIF

      ELSE IF (CRS_CCS == 'CCS') THEN                   ! Convert CRS input to Compressed Col Storage (CCS) format

         IF (ALLOCATED(CCS1)) CALL DEALLOCATE_SCR_MAT ( 'CCS1' )
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NROWS, NTERM_SUPRLU, SUBR_NAME )
         IF (MATIN_IS_FULL) THEN
            CALL SPARSE_CRS_SPARSE_CCS ( NROWS, NROWS, NTERM_SUPRLU, MATIN_NAME, I_MATIN, J_MATIN, MATIN, 'CCS1', J_CCS1, I_CCS1, &
                                        CCS1, 'Y' )
         ELSE
            CALL SPARSE_CRS_SPARSE_CCS ( NROWS, NROWS, NTERM_SUPRLU, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CCS1', J_CCS1, I_CCS1, CCS1, &
                                        'Y' )
         ENDIF
         HS_PHASE_T0 = HOTSPOT_WALL_TIME()
         CALL C_FORTRAN_DGSSV( 1, NROWS, NTERM_SUPRLU, 1, CCS1, I_CCS1, J_CCS1, DUM_COL, NROWS, SLU_FACTORS, INFO )
         CALL HOTSPOT_TIMER_ADD ( 'SYM_MAT_DECOMP_SUPRLU/SUPERLU_CALL', HOTSPOT_WALL_TIME() - HS_PHASE_T0 )
         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

      ELSE

         WRITE(ERR,933) SUBR_NAME, 'CRS_CCS'
         WRITE(F06,933) SUBR_NAME, 'CRS_CCS'
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF (.NOT.MATIN_IS_FULL) THEN
         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )
      ENDIF


      IF (INFO == 0) THEN

         WRITE (SC1,9902) MATIN_NAME, SUBR_NAME 
         WRITE (F06,9902) MATIN_NAME, SUBR_NAME

      ELSE IF (INFO < 0) THEN                              ! Illegal value of an argument to SuperLU

         WRITE(SC1,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         WRITE(ERR,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         WRITE(F06,9903) INFO, TRIM(SUBR_NAME), TRIM(CALLING_SUBR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE IF (INFO > 0) THEN                              ! Singular matrix, memory error, or other.

        CALL GET_GRID_AND_COMP ( MATIN_SET, INFO, GRIDV, COMPV  )

        WRITE(ERR,981) MATIN_NAME, INFO
        WRITE(F06,981) MATIN_NAME, INFO
        IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
          WRITE(ERR,9811) GRIDV, COMPV, CALLING_SUBR 
          WRITE(F06,9811) GRIDV, COMPV, CALLING_SUBR
        ELSE 
          WRITE(ERR,9812) INFO, CALLING_SUBR 
          WRITE(F06,9812) INFO, CALLING_SUBR
        ENDIF
         
      ENDIF

      ! This should also use BAILOUT_CHECK like SYM_MAT_DECOMP_LAPACK does, however we need to 
      ! build an array of the diagonal values of the U factor. SLU_FACTORS is a pointer to the 
      ! structure containing L, U, perm_c, perm_r. See superlu/FORTRAN/c_fortran_dgssv.c. We
      ! may need to look up the permutations to find the corresponding rows of U.
      !
      ! It still respects the BAILOUT parameter when SuperLU reports an error while processing 
      ! the matrix which seems to work in a similar way.

      IF ((INFO > 0)) THEN
                                                           ! If BAILOUT >= 0 then quit. Otherwise, continue processing.
        IF (BAILOUT >= 0) THEN
          FATAL_ERR = FATAL_ERR + 1
          WRITE(ERR,99999) BAILOUT
          WRITE(F06,99999) BAILOUT
          CALL OUTA_HERE ( 'Y' )
        ENDIF

      ENDIF


! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
  9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      CALL HOTSPOT_TIMER_END ( HS_SLOT, HS_T0 )

      RETURN

!***********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER ', A, ' MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  933 FORMAT(' *ERROR   933: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER ', A, ' MUST BE EITHER "CRS" OR "CCS" BUT VALUE IS ',A)

  981 FORMAT(' *ERROR   981: THE FACTORIZATION OF THE MATRIX ',A,' BY SUPERLU HAD ERROR WITH INFO = ', I12, '.')

 9811 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMP ',I3,'. THE CALLING SUBR WAS: ',A,/)

 9812 FORMAT('               THIS IS FOR ROW AND COL ',I8,' IN THE MATRIX. THE CALLING SUBR WAS: ',A,/)

 9902 FORMAT(' SUPERLU FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)

 9903 FORMAT(' *ERROR  9903: SUPERLU SPARSE SOLVER HAS FAILED WITH INFO = ', I12,' IN SUBR ', A, ' CALLED BY SUBR ', A)

99999 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE MESSAGES AND BULK DATA PARAMETER BAILOUT = ',I7)

!***********************************************************************************************************************************

      END SUBROUTINE SYM_MAT_DECOMP_SUPRLU

