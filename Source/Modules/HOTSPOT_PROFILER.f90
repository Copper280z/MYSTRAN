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

      MODULE HOTSPOT_PROFILER

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DBL_LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), PARAMETER        :: MAX_TIMERS   = 512
      INTEGER(LONG), PARAMETER        :: MAX_COUNTERS = 512
      INTEGER(LONG), PARAMETER        :: MAX_VALUES   = 256
      INTEGER(LONG), PARAMETER        :: NAME_LEN     = 128

      CHARACTER(NAME_LEN*BYTE)        :: TIMER_NAMES  (MAX_TIMERS)
      CHARACTER(NAME_LEN*BYTE)        :: COUNTER_NAMES(MAX_COUNTERS)
      CHARACTER(NAME_LEN*BYTE)        :: VALUE_NAMES  (MAX_VALUES)

      REAL(DOUBLE)                    :: TIMER_TOTALS (MAX_TIMERS)
      REAL(DOUBLE)                    :: VALUE_SUMS   (MAX_VALUES)
      REAL(DOUBLE)                    :: VALUE_MAXIMA (MAX_VALUES)

      INTEGER(DBL_LONG)               :: TIMER_CALLS  (MAX_TIMERS)
      INTEGER(DBL_LONG)               :: COUNTER_VALS (MAX_COUNTERS)
      INTEGER(DBL_LONG)               :: VALUE_COUNTS (MAX_VALUES)

      INTEGER(LONG)                   :: NUM_TIMERS   = 0
      INTEGER(LONG)                   :: NUM_COUNTERS = 0
      INTEGER(LONG)                   :: NUM_VALUES   = 0

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_PROFILER_INIT

      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      NUM_TIMERS   = 0
      NUM_COUNTERS = 0
      NUM_VALUES   = 0

      DO I=1,MAX_TIMERS
         TIMER_NAMES (I)(1:) = ' '
         TIMER_TOTALS(I)     = 0.0D0
         TIMER_CALLS (I)     = 0
      ENDDO

      DO I=1,MAX_COUNTERS
         COUNTER_NAMES(I)(1:) = ' '
         COUNTER_VALS (I)     = 0
      ENDDO

      DO I=1,MAX_VALUES
         VALUE_NAMES (I)(1:) = ' '
         VALUE_SUMS  (I)     = 0.0D0
         VALUE_MAXIMA(I)     = 0.0D0
         VALUE_COUNTS(I)     = 0
      ENDDO

      END SUBROUTINE HOTSPOT_PROFILER_INIT

! ##################################################################################################################################

      REAL(DOUBLE) FUNCTION HOTSPOT_WALL_TIME ()

      INTEGER(DBL_LONG)               :: COUNT
      INTEGER(DBL_LONG)               :: COUNT_RATE
      INTEGER(DBL_LONG)               :: COUNT_MAX

! **********************************************************************************************************************************
      CALL SYSTEM_CLOCK ( COUNT, COUNT_RATE, COUNT_MAX )

      IF (COUNT_RATE > 0) THEN
         HOTSPOT_WALL_TIME = DBLE(COUNT)/DBLE(COUNT_RATE)
      ELSE
         HOTSPOT_WALL_TIME = 0.0D0
      ENDIF

      END FUNCTION HOTSPOT_WALL_TIME

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_TIMER_BEGIN ( NAME, SLOT, T0 )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)  , INTENT(OUT)    :: SLOT
      REAL(DOUBLE)   , INTENT(OUT)    :: T0

! **********************************************************************************************************************************
      SLOT = GET_TIMER_SLOT ( NAME )
      T0   = HOTSPOT_WALL_TIME()

      END SUBROUTINE HOTSPOT_TIMER_BEGIN

! ##################################################################################################################################

      CHARACTER(6*BYTE) FUNCTION HOTSPOT_OPT_MASK ( OPT )

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      DO I=1,6
         HOTSPOT_OPT_MASK(I:I) = OPT(I)
      ENDDO

      END FUNCTION HOTSPOT_OPT_MASK

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_TIMER_END ( SLOT, T0 )

      INTEGER(LONG), INTENT(IN)       :: SLOT
      REAL(DOUBLE) , INTENT(IN)       :: T0

      REAL(DOUBLE)                    :: DT

! **********************************************************************************************************************************
      IF ((SLOT < 1) .OR. (SLOT > NUM_TIMERS)) RETURN

      DT = HOTSPOT_WALL_TIME() - T0
      IF (DT < 0.0D0) DT = 0.0D0

      TIMER_TOTALS(SLOT) = TIMER_TOTALS(SLOT) + DT
      TIMER_CALLS (SLOT) = TIMER_CALLS (SLOT) + 1

      END SUBROUTINE HOTSPOT_TIMER_END

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_TIMER_ADD ( NAME, DT )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      REAL(DOUBLE)   , INTENT(IN)     :: DT

      INTEGER(LONG)                   :: SLOT
      REAL(DOUBLE)                    :: DT_USE

! **********************************************************************************************************************************
      SLOT = GET_TIMER_SLOT ( NAME )

      DT_USE = DT
      IF (DT_USE < 0.0D0) DT_USE = 0.0D0

      TIMER_TOTALS(SLOT) = TIMER_TOTALS(SLOT) + DT_USE
      TIMER_CALLS (SLOT) = TIMER_CALLS (SLOT) + 1

      END SUBROUTINE HOTSPOT_TIMER_ADD

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_COUNTER_ADD ( NAME, DELTA )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(DBL_LONG), INTENT(IN)   :: DELTA

      INTEGER(LONG)                   :: SLOT

! **********************************************************************************************************************************
      SLOT = GET_COUNTER_SLOT ( NAME )
      COUNTER_VALS(SLOT) = COUNTER_VALS(SLOT) + DELTA

      END SUBROUTINE HOTSPOT_COUNTER_ADD

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_VALUE_ADD ( NAME, VALUE )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      REAL(DOUBLE)   , INTENT(IN)     :: VALUE

      INTEGER(LONG)                   :: SLOT

! **********************************************************************************************************************************
      SLOT = GET_VALUE_SLOT ( NAME )

      VALUE_SUMS  (SLOT) = VALUE_SUMS  (SLOT) + VALUE
      VALUE_COUNTS(SLOT) = VALUE_COUNTS(SLOT) + 1

      IF (VALUE_COUNTS(SLOT) == 1) THEN
         VALUE_MAXIMA(SLOT) = VALUE
      ELSE IF (VALUE > VALUE_MAXIMA(SLOT)) THEN
         VALUE_MAXIMA(SLOT) = VALUE
      ENDIF

      END SUBROUTINE HOTSPOT_VALUE_ADD

! ##################################################################################################################################

      REAL(DOUBLE) FUNCTION HOTSPOT_TIMER_TOTAL ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      HOTSPOT_TIMER_TOTAL = 0.0D0

      DO I=1,NUM_TIMERS
         IF (TRIM(TIMER_NAMES(I)) == TRIM(NAME)) THEN
            HOTSPOT_TIMER_TOTAL = TIMER_TOTALS(I)
            RETURN
         ENDIF
      ENDDO

      END FUNCTION HOTSPOT_TIMER_TOTAL

! ##################################################################################################################################

      INTEGER(DBL_LONG) FUNCTION HOTSPOT_TIMER_CALLS ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      HOTSPOT_TIMER_CALLS = 0

      DO I=1,NUM_TIMERS
         IF (TRIM(TIMER_NAMES(I)) == TRIM(NAME)) THEN
            HOTSPOT_TIMER_CALLS = TIMER_CALLS(I)
            RETURN
         ENDIF
      ENDDO

      END FUNCTION HOTSPOT_TIMER_CALLS

! ##################################################################################################################################

      INTEGER(DBL_LONG) FUNCTION HOTSPOT_COUNTER_VALUE ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      HOTSPOT_COUNTER_VALUE = 0

      DO I=1,NUM_COUNTERS
         IF (TRIM(COUNTER_NAMES(I)) == TRIM(NAME)) THEN
            HOTSPOT_COUNTER_VALUE = COUNTER_VALS(I)
            RETURN
         ENDIF
      ENDDO

      END FUNCTION HOTSPOT_COUNTER_VALUE

! ##################################################################################################################################

      REAL(DOUBLE) FUNCTION HOTSPOT_VALUE_SUM ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      HOTSPOT_VALUE_SUM = 0.0D0

      DO I=1,NUM_VALUES
         IF (TRIM(VALUE_NAMES(I)) == TRIM(NAME)) THEN
            HOTSPOT_VALUE_SUM = VALUE_SUMS(I)
            RETURN
         ENDIF
      ENDDO

      END FUNCTION HOTSPOT_VALUE_SUM

! ##################################################################################################################################

      REAL(DOUBLE) FUNCTION HOTSPOT_VALUE_MAX ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      HOTSPOT_VALUE_MAX = 0.0D0

      DO I=1,NUM_VALUES
         IF (TRIM(VALUE_NAMES(I)) == TRIM(NAME)) THEN
            HOTSPOT_VALUE_MAX = VALUE_MAXIMA(I)
            RETURN
         ENDIF
      ENDDO

      END FUNCTION HOTSPOT_VALUE_MAX

! ##################################################################################################################################

      INTEGER(DBL_LONG) FUNCTION HOTSPOT_VALUE_COUNT ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      HOTSPOT_VALUE_COUNT = 0

      DO I=1,NUM_VALUES
         IF (TRIM(VALUE_NAMES(I)) == TRIM(NAME)) THEN
            HOTSPOT_VALUE_COUNT = VALUE_COUNTS(I)
            RETURN
         ENDIF
      ENDDO

      END FUNCTION HOTSPOT_VALUE_COUNT

! ##################################################################################################################################

      SUBROUTINE HOTSPOT_WRITE_REPORT ( UNT )

      INTEGER(LONG), INTENT(IN)       :: UNT

      REAL(DOUBLE)                    :: AVG_ROW_LEN
      INTEGER(DBL_LONG)               :: NROW

! **********************************************************************************************************************************
      WRITE(UNT,1000)

      CALL WRITE_TIMER_LINE   ( UNT, 'OFP3_ELFE_2D',          'OFP3_ELFE_2D' )
      CALL WRITE_TIMER_LINE   ( UNT, 'OFP3_STRE_NO_PCOMP',    'OFP3_STRE_NO_PCOMP' )
      CALL WRITE_TIMER_LINE   ( UNT, 'OFP3_STRN_NO_PCOMP',    'OFP3_STRN_NO_PCOMP' )
      CALL WRITE_TIMER_LINE   ( UNT, 'ESP',                   'ESP' )
      CALL WRITE_TIMER_LINE   ( UNT, 'REDUCE_G_NM',           'REDUCE_G_NM' )
      CALL WRITE_TIMER_LINE   ( UNT, 'ELEM_TRANSFORM_LBG',    'ELEM_TRANSFORM_LBG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_KGG',            'SPARSE_KGG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'N_SET_AUTOSPC_PROC_1',  'N_SET_AUTOSPC_PROC_1' )
      CALL WRITE_TIMER_LINE   ( UNT, 'N_SET_AUTOSPC_PROC_2',  'N_SET_AUTOSPC_PROC_2' )
      CALL WRITE_TIMER_LINE   ( UNT, 'TDOF_PROC',             'TDOF_PROC' )
      CALL WRITE_TIMER_LINE   ( UNT, 'WRITE_DOF_TABLES',      'WRITE_DOF_TABLES' )
      CALL WRITE_TIMER_LINE   ( UNT, 'EMP',                   'EMP' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_MGG',            'SPARSE_MGG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'QMEM1',                 'QMEM1' )
      CALL WRITE_TIMER_LINE   ( UNT, 'QPLT2',                 'QPLT2' )
      CALL WRITE_TIMER_LINE   ( UNT, 'POLYNOM_FIT_STRE_STRN', 'POLYNOM_FIT_STRE_STRN' )
      CALL WRITE_TIMER_LINE   ( UNT, 'EMG',                   'EMG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'GET_ELEM_AGRID_BGRID',  'GET_ELEM_AGRID_BGRID' )
      CALL WRITE_TIMER_LINE   ( UNT, 'ASSERT_ARRAY_SORTED',   'ASSERT_ARRAY_SORTED' )
      CALL WRITE_TIMER_LINE   ( UNT, 'GET_ARRAY_ROW_NUM',     'GET_ARRAY_ROW_NUM' )

      WRITE(UNT,1010)
      CALL WRITE_TIMER_LINE   ( UNT, 'ESP/EMG',                         'ESP internal: EMG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'ESP/ELEM_TRANSFORM_LBG',          'ESP internal: ELEM_TRANSFORM_LBG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'ESP/STIFF_INSERT',                'ESP internal: stiffness insertion' )
      CALL WRITE_TIMER_LINE   ( UNT, 'EMP/EMG',                         'EMP internal: EMG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'EMP/ELEM_TRANSFORM_LBG',          'EMP internal: ELEM_TRANSFORM_LBG' )
      CALL WRITE_TIMER_LINE   ( UNT, 'EMP/MASS_INSERT',                 'EMP internal: mass insertion' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_KGG/ZERO_STRIP',           'SPARSE_KGG: recount / zero strip' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_KGG/ROW_EXTRACT',          'SPARSE_KGG: row extraction' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_KGG/ROW_SORT',             'SPARSE_KGG: row sort' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_KGG/SINGULARITY_PROC',     'SPARSE_KGG: singularity processing' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_MGG/ZERO_STRIP',           'SPARSE_MGG: recount / zero strip' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_MGG/ROW_EXTRACT',          'SPARSE_MGG: row extraction' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_MGG/ROW_SORT',             'SPARSE_MGG: row sort' )
      CALL WRITE_TIMER_LINE   ( UNT, 'SPARSE_MGG/MATADD_MERGES',        'SPARSE_MGG: MATADD merges' )
      CALL WRITE_TIMER_LINE   ( UNT, 'GET_ELEM_AGRID_BGRID/CONNECT',    'GET_ELEM_AGRID_BGRID: connectivity lookups' )
      CALL WRITE_TIMER_LINE   ( UNT, 'GET_ELEM_AGRID_BGRID/CHECK_AGRID','GET_ELEM_AGRID_BGRID: CHECK_AGRID branch' )
      CALL WRITE_TIMER_LINE   ( UNT, 'AUTOSPC/TDOFI_SEARCH/PROC1',      'AUTOSPC proc1: TDOFI search' )
      CALL WRITE_TIMER_LINE   ( UNT, 'AUTOSPC/TDOFI_SEARCH/PROC2',      'AUTOSPC proc2: TDOFI search' )

      WRITE(UNT,1020)
      CALL WRITE_PREFIX_TIMERS ( UNT, 'GET_ARRAY_ROW_NUM/CALLER/', 'GET_ARRAY_ROW_NUM by caller' )
      CALL WRITE_PREFIX_TIMERS ( UNT, 'EMG/CALL/',                  'EMG calls by caller / element / OPT' )
      CALL WRITE_PREFIX_TIMERS ( UNT, 'EMP/EMG/',                  'EMP EMG calls by element type' )
      CALL WRITE_PREFIX_TIMERS ( UNT, 'QMEM1/OPT/',                 'QMEM1 calls by OPT mask' )
      CALL WRITE_PREFIX_TIMERS ( UNT, 'QPLT2/OPT/',                 'QPLT2 calls by OPT mask' )
      CALL WRITE_PREFIX_TIMERS ( UNT, 'HELPER/',                    'Helper kernel timers' )

      WRITE(UNT,1030)
      CALL WRITE_COUNTER_LINE ( UNT, 'LINK9_REQUEST/ELEMENTS/ELFE',                'Elements requesting ELFE' )
      CALL WRITE_COUNTER_LINE ( UNT, 'LINK9_REQUEST/ELEMENTS/STRE',                'Elements requesting STRE' )
      CALL WRITE_COUNTER_LINE ( UNT, 'LINK9_REQUEST/ELEMENTS/STRN',                'Elements requesting STRN' )
      CALL WRITE_COUNTER_LINE ( UNT, 'GRID_LOOKUP/REPLACEABLE/TOTAL',              'GRID_ID lookups replaceable with BGRID' )
      CALL WRITE_COUNTER_LINE ( UNT, 'GRID_LOOKUP/REQUIRED/TOTAL',                 'GRID_ID lookups requiring a real search' )
      CALL WRITE_COUNTER_LINE ( UNT, 'AUTOSPC/NULL_ROWS_FOUND',                    'AUTOSPC null rows found' )
      CALL WRITE_COUNTER_LINE ( UNT, 'AUTOSPC/SMALL_DIAG_ROWS_FOUND',              'AUTOSPC small-diagonal rows found' )
      CALL WRITE_COUNTER_LINE ( UNT, 'KGG_DUPLICATE_INSERTS',                      'KGG duplicate inserts combined' )
      CALL WRITE_COUNTER_LINE ( UNT, 'KGG_ZERO_DROPS',                             'KGG entries dropped as near-zero' )
      CALL WRITE_COUNTER_LINE ( UNT, 'MGGE_DUPLICATE_INSERTS',                     'MGGE duplicate inserts combined' )
      CALL WRITE_COUNTER_LINE ( UNT, 'MGGE_ZERO_DROPS',                            'MGGE entries dropped as near-zero' )
      CALL WRITE_COUNTER_LINE ( UNT, 'MGGE_NEW_TERMS/CQUAD4',                      'MGGE new linked-list terms from CQUAD4' )
      CALL WRITE_COUNTER_LINE ( UNT, 'MGGE_NEW_TERMS/OTHER',                       'MGGE new linked-list terms from other elems' )

      NROW = HOTSPOT_VALUE_COUNT ( 'KGG_ROWLEN_RAW' )
      IF (NROW > 0) THEN
         AVG_ROW_LEN = HOTSPOT_VALUE_SUM ( 'KGG_ROWLEN_RAW' )/DBLE(NROW)
         WRITE(UNT,1040) 'KGG raw row length avg/max', AVG_ROW_LEN, HOTSPOT_VALUE_MAX('KGG_ROWLEN_RAW')
      ENDIF

      NROW = HOTSPOT_VALUE_COUNT ( 'KGG_ROWLEN_FINAL' )
      IF (NROW > 0) THEN
         AVG_ROW_LEN = HOTSPOT_VALUE_SUM ( 'KGG_ROWLEN_FINAL' )/DBLE(NROW)
         WRITE(UNT,1040) 'KGG final row length avg/max', AVG_ROW_LEN, HOTSPOT_VALUE_MAX('KGG_ROWLEN_FINAL')
      ENDIF

      NROW = HOTSPOT_VALUE_COUNT ( 'MGGE_ROWLEN_RAW' )
      IF (NROW > 0) THEN
         AVG_ROW_LEN = HOTSPOT_VALUE_SUM ( 'MGGE_ROWLEN_RAW' )/DBLE(NROW)
         WRITE(UNT,1040) 'MGGE raw row length avg/max', AVG_ROW_LEN, HOTSPOT_VALUE_MAX('MGGE_ROWLEN_RAW')
      ENDIF

      NROW = HOTSPOT_VALUE_COUNT ( 'MGGE_ROWLEN_FINAL' )
      IF (NROW > 0) THEN
         AVG_ROW_LEN = HOTSPOT_VALUE_SUM ( 'MGGE_ROWLEN_FINAL' )/DBLE(NROW)
         WRITE(UNT,1040) 'MGGE final row length avg/max', AVG_ROW_LEN, HOTSPOT_VALUE_MAX('MGGE_ROWLEN_FINAL')
      ENDIF

      NROW = HOTSPOT_VALUE_COUNT ( 'MGG_ROWLEN_FINAL' )
      IF (NROW > 0) THEN
         AVG_ROW_LEN = HOTSPOT_VALUE_SUM ( 'MGG_ROWLEN_FINAL' )/DBLE(NROW)
         WRITE(UNT,1040) 'MGG final row length avg/max', AVG_ROW_LEN, HOTSPOT_VALUE_MAX('MGG_ROWLEN_FINAL')
      ENDIF

 1000 FORMAT(/,1X,'HOTSPOT PROFILER SUMMARY',/,1X,'========================')
 1010 FORMAT(/,1X,'Internal phase timers:')
 1020 FORMAT(/,1X,'Call-shape summaries:')
 1030 FORMAT(/,1X,'Key counters:')
 1040 FORMAT(1X,'  ',A40,' : ',F12.3,' / ',F12.3)

      END SUBROUTINE HOTSPOT_WRITE_REPORT

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION GET_TIMER_SLOT ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: LNAME

! **********************************************************************************************************************************
      DO I=1,NUM_TIMERS
         IF (TRIM(TIMER_NAMES(I)) == TRIM(NAME)) THEN
            GET_TIMER_SLOT = I
            RETURN
         ENDIF
      ENDDO

      NUM_TIMERS = NUM_TIMERS + 1
      IF (NUM_TIMERS > MAX_TIMERS) THEN
         GET_TIMER_SLOT = MAX_TIMERS
         RETURN
      ENDIF

      TIMER_NAMES (NUM_TIMERS)(1:) = ' '
      LNAME = MIN(LEN_TRIM(NAME),NAME_LEN)
      IF (LNAME > 0) THEN
         TIMER_NAMES (NUM_TIMERS)(1:LNAME) = NAME(1:LNAME)
      ENDIF
      TIMER_TOTALS(NUM_TIMERS)      = 0.0D0
      TIMER_CALLS (NUM_TIMERS)      = 0

      GET_TIMER_SLOT = NUM_TIMERS

      END FUNCTION GET_TIMER_SLOT

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION GET_COUNTER_SLOT ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: LNAME

! **********************************************************************************************************************************
      DO I=1,NUM_COUNTERS
         IF (TRIM(COUNTER_NAMES(I)) == TRIM(NAME)) THEN
            GET_COUNTER_SLOT = I
            RETURN
         ENDIF
      ENDDO

      NUM_COUNTERS = NUM_COUNTERS + 1
      IF (NUM_COUNTERS > MAX_COUNTERS) THEN
         GET_COUNTER_SLOT = MAX_COUNTERS
         RETURN
      ENDIF

      COUNTER_NAMES(NUM_COUNTERS)(1:) = ' '
      LNAME = MIN(LEN_TRIM(NAME),NAME_LEN)
      IF (LNAME > 0) THEN
         COUNTER_NAMES(NUM_COUNTERS)(1:LNAME) = NAME(1:LNAME)
      ENDIF
      COUNTER_VALS (NUM_COUNTERS) = 0

      GET_COUNTER_SLOT = NUM_COUNTERS

      END FUNCTION GET_COUNTER_SLOT

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION GET_VALUE_SLOT ( NAME )

      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: LNAME

! **********************************************************************************************************************************
      DO I=1,NUM_VALUES
         IF (TRIM(VALUE_NAMES(I)) == TRIM(NAME)) THEN
            GET_VALUE_SLOT = I
            RETURN
         ENDIF
      ENDDO

      NUM_VALUES = NUM_VALUES + 1
      IF (NUM_VALUES > MAX_VALUES) THEN
         GET_VALUE_SLOT = MAX_VALUES
         RETURN
      ENDIF

      VALUE_NAMES (NUM_VALUES)(1:) = ' '
      LNAME = MIN(LEN_TRIM(NAME),NAME_LEN)
      IF (LNAME > 0) THEN
         VALUE_NAMES (NUM_VALUES)(1:LNAME) = NAME(1:LNAME)
      ENDIF
      VALUE_SUMS  (NUM_VALUES) = 0.0D0
      VALUE_MAXIMA(NUM_VALUES) = 0.0D0
      VALUE_COUNTS(NUM_VALUES) = 0

      GET_VALUE_SLOT = NUM_VALUES

      END FUNCTION GET_VALUE_SLOT

! ##################################################################################################################################

      SUBROUTINE WRITE_TIMER_LINE ( UNT, NAME, LABEL )

      INTEGER(LONG), INTENT(IN)       :: UNT
      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      CHARACTER(LEN=*), INTENT(IN)    :: LABEL

      REAL(DOUBLE)                    :: TOTAL
      INTEGER(DBL_LONG)               :: NCALLS

! **********************************************************************************************************************************
      TOTAL  = HOTSPOT_TIMER_TOTAL ( NAME )
      NCALLS = HOTSPOT_TIMER_CALLS ( NAME )

      IF (NCALLS > 0) THEN
         WRITE(UNT,2000) LABEL, TOTAL, NCALLS
      ENDIF

 2000 FORMAT(1X,'  ',A40,' : ',F12.3,' s in ',I12,' call(s)')

      END SUBROUTINE WRITE_TIMER_LINE

! ##################################################################################################################################

      SUBROUTINE WRITE_COUNTER_LINE ( UNT, NAME, LABEL )

      INTEGER(LONG), INTENT(IN)       :: UNT
      CHARACTER(LEN=*), INTENT(IN)    :: NAME
      CHARACTER(LEN=*), INTENT(IN)    :: LABEL

      INTEGER(DBL_LONG)               :: VALUE

! **********************************************************************************************************************************
      VALUE = HOTSPOT_COUNTER_VALUE ( NAME )
      IF (VALUE > 0) THEN
         WRITE(UNT,2100) LABEL, VALUE
      ENDIF

 2100 FORMAT(1X,'  ',A40,' : ',I18)

      END SUBROUTINE WRITE_COUNTER_LINE

! ##################################################################################################################################

      SUBROUTINE WRITE_PREFIX_TIMERS ( UNT, PREFIX, TITLE )

      INTEGER(LONG), INTENT(IN)       :: UNT
      CHARACTER(LEN=*), INTENT(IN)    :: PREFIX
      CHARACTER(LEN=*), INTENT(IN)    :: TITLE

      LOGICAL                         :: ANY_WRITTEN
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: LP
      CHARACTER(NAME_LEN*BYTE)        :: LABEL

! **********************************************************************************************************************************
      ANY_WRITTEN = .FALSE.
      LP = LEN_TRIM(PREFIX)

      DO I=1,NUM_TIMERS
         IF (INDEX(TRIM(TIMER_NAMES(I)), TRIM(PREFIX)) == 1) THEN
            IF (.NOT. ANY_WRITTEN) THEN
               WRITE(UNT,2200) TITLE
               ANY_WRITTEN = .TRUE.
            ENDIF
            LABEL(1:) = ' '
            IF (LEN_TRIM(TIMER_NAMES(I)) > LP) THEN
               LABEL(1:LEN_TRIM(TIMER_NAMES(I))-LP) = TIMER_NAMES(I)(LP+1:LEN_TRIM(TIMER_NAMES(I)))
            ELSE
               LABEL(1:1) = '?'
            ENDIF
            WRITE(UNT,2210) TRIM(LABEL), TIMER_TOTALS(I), TIMER_CALLS(I)
         ENDIF
      ENDDO

 2200 FORMAT(1X,'  ',A,':')
 2210 FORMAT(1X,'    ',A48,' : ',F12.3,' s in ',I12,' call(s)')

      END SUBROUTINE WRITE_PREFIX_TIMERS

! ##################################################################################################################################

      END MODULE HOTSPOT_PROFILER
