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

      SUBROUTINE OFP3_SHELL_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Combined non-PCOMP shell output driver for force, stress, and strain output.
! The fast path collapses the repeated EMG traversal for the standard F06/OP2 path.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE, DBL_LONG
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_LOG, ERR, F04, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_ELFE_BIT, ELOUT_STRE_BIT, ELOUT_STRN_BIT, FATAL_ERR, IBIT,           &
                                         INT_SC_NUM, MBUG, MOGEL, NELE, NCQUAD4, NCQUAD4K, NCSHEAR, NCTRIA3, NCTRIA3K, SOL_NAME, &
                                         MAX_STRESS_POINTS, MELGP
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  FORC_LOC, STRE_LOC, STRN_LOC
      USE MODEL_STUF, ONLY            :  AGRID, EDAT, EPNT, ETYPE, EID, ELGP, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS,           &
                                         NUM_SEi, PCOMP_PROPS, PLY_NUM, STRESS, STRAIN, TYPE, SHELL_STR_ANGLE
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL, POLY_FIT_ERR, POLY_FIT_ERR_INDEX
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRE, TXT_STRE, OTM_STRN, TXT_STRN, OTM_ELFE, TXT_ELFE
      USE SUBR_BEGEND_LEVELS, ONLY    :  OFP3_BEGEND
      USE OFP3_ELFE_2D_USE_IFs
      USE OFP3_STRE_NO_PCOMP_USE_IFs
      USE OFP3_STRN_NO_PCOMP_USE_IFs
      USE CALC_ELEM_STRESSES_USE_IFs
      USE CALC_ELEM_STRAINS_USE_IFs
      USE SHELL_ENGR_FORCE_OGEL_USE_IFs
      USE PLANE_COORD_TRANS_21_Interface
      USE TRANSFORM_SHELL_STR_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_SHELL_NO_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'
      CHARACTER( 1*BYTE)              :: OPT(6)
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR
      CHARACTER(20*BYTE)              :: FORCE_ITEM(8)
      CHARACTER(20*BYTE)              :: STRESS_ITEM(11)
      CHARACTER(20*BYTE)              :: STRAIN_ITEM(11)

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID
      INTEGER(LONG), INTENT(IN)       :: ITE
      INTEGER(LONG), INTENT(IN)       :: JVEC
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = OFP3_BEGEND

      INTEGER(LONG)                   :: I,J,K,M,R
      INTEGER(LONG)                   :: IERROR = 0
      INTEGER(LONG)                   :: NUM_PTS_ELFE(METYPE), NUM_PTS_STRE(METYPE), NUM_PTS_STRN(METYPE)
      INTEGER(LONG)                   :: NELREQ_ELFE(METYPE), NELREQ_STRE(METYPE), NELREQ_STRN(METYPE)
      INTEGER(LONG)                   :: NUM_ELFE_ROWS, NUM_STRE_ROWS, NUM_STRN_ROWS
      INTEGER(LONG)                   :: NUM_ELFE_PREV, NUM_STRE_PREV, NUM_STRN_PREV
      INTEGER(LONG)                   :: NUM_OGEL_ELFE, NUM_OGEL_STRE, NUM_OGEL_STRN
      INTEGER(LONG)                   :: NUM_OGEL_TMP
      INTEGER(LONG)                   :: NUM_FROWS
      INTEGER(LONG)                   :: NUM_LOCAL
      INTEGER(LONG)                   :: NROWS_FORCE, NROWS_STRE, NROWS_STRN
      INTEGER(LONG)                   :: STRESS_OUT_ERR_INDEX(MAX_STRESS_POINTS)
      INTEGER(LONG)                   :: STRAIN_OUT_ERR_INDEX(MAX_STRESS_POINTS)
      INTEGER(LONG)                   :: ITABLE
      CHARACTER(8*BYTE)               :: TABLE_NAME
      LOGICAL                         :: WRITE_NEU
      LOGICAL                         :: FAST_PATH
      REAL(DOUBLE)                    :: PCT_ERR_MAX
      REAL(DOUBLE)                    :: STRESS_OUT_PCT_ERR(MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRAIN_OUT_PCT_ERR(MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRESS_RAW(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRAIN_RAW(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRESS_OUT(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRAIN_OUT(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: TEL(3,3)

      REAL(DOUBLE), ALLOCATABLE       :: OGEL_ELFE_LOC(:,:), OGEL_STRE_LOC(:,:), OGEL_STRN_LOC(:,:)
      INTEGER(LONG), ALLOCATABLE      :: EID_ELFE_LOC(:,:), EID_STRE_LOC(:,:), EID_STRN_LOC(:,:)
      INTEGER(LONG), ALLOCATABLE      :: GID_ELFE_LOC(:,:), GID_STRE_LOC(:,:), GID_STRN_LOC(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: POLY_FIT_ERR_STRE_LOC(:), POLY_FIT_ERR_STRN_LOC(:)
      INTEGER(LONG), ALLOCATABLE      :: POLY_FIT_ERR_INDEX_STRE_LOC(:), POLY_FIT_ERR_INDEX_STRN_LOC(:)

      INTEGER(LONG), PARAMETER        :: ZERO_I = 0
      INTRINSIC IAND

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

      WRITE_NEU = (PRTNEU == 'Y')
      FAST_PATH = (PRTNEU /= 'Y') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')

      IF (.NOT. FAST_PATH) THEN
         CALL OFP3_ELFE_2D       ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )
         CALL OFP3_STRE_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )
         CALL OFP3_STRN_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )
         GOTO 9999
      ENDIF

      OPT(1) = 'N'
      OPT(2) = 'N'
      OPT(3) = 'Y'
      OPT(4) = 'N'
      OPT(5) = 'N'
      OPT(6) = 'N'

      FORCE_ITEM(1) = 'Nxx: Normal x Force '
      FORCE_ITEM(2) = 'Nyy: Normal y Force '
      FORCE_ITEM(3) = 'Nxy: Shear xy Force '
      FORCE_ITEM(4) = 'Mxx: Moment x Plane '
      FORCE_ITEM(5) = 'Myy: Moment y Plane '
      FORCE_ITEM(6) = 'Mxy: Twist Mom xy   '
      FORCE_ITEM(7) = 'Qx : Transv Shear x '
      FORCE_ITEM(8) = 'Qy : Transv Shear y '

      ALLOCATE ( OGEL_ELFE_LOC(MAXREQ,MOGEL) )
      ALLOCATE ( OGEL_STRE_LOC(MAXREQ,MOGEL) )
      ALLOCATE ( OGEL_STRN_LOC(MAXREQ,MOGEL) )
      ALLOCATE ( EID_ELFE_LOC(MAXREQ,2), EID_STRE_LOC(MAXREQ,2), EID_STRN_LOC(MAXREQ,2) )
      ALLOCATE ( GID_ELFE_LOC(MAXREQ,MELGP+1), GID_STRE_LOC(MAXREQ,MELGP+1), GID_STRN_LOC(MAXREQ,MELGP+1) )
      ALLOCATE ( POLY_FIT_ERR_STRE_LOC(MAXREQ), POLY_FIT_ERR_STRN_LOC(MAXREQ) )
      ALLOCATE ( POLY_FIT_ERR_INDEX_STRE_LOC(MAXREQ), POLY_FIT_ERR_INDEX_STRN_LOC(MAXREQ) )

      OGEL_ELFE_LOC           = ZERO
      OGEL_STRE_LOC           = ZERO
      OGEL_STRN_LOC           = ZERO
      EID_ELFE_LOC            = ZERO_I
      EID_STRE_LOC            = ZERO_I
      EID_STRN_LOC            = ZERO_I
      GID_ELFE_LOC            = ZERO_I
      GID_STRE_LOC            = ZERO_I
      GID_STRN_LOC            = ZERO_I
      POLY_FIT_ERR_STRE_LOC   = ZERO
      POLY_FIT_ERR_STRN_LOC   = ZERO
      POLY_FIT_ERR_INDEX_STRE_LOC = ZERO_I
      POLY_FIT_ERR_INDEX_STRN_LOC = ZERO_I

      DO I=1,METYPE
         NELREQ_ELFE(I) = 0
         NELREQ_STRE(I) = 0
         NELREQ_STRN(I) = 0
         NUM_PTS_ELFE(I) = 0
         NUM_PTS_STRE(I) = 0
         NUM_PTS_STRN(I) = 0
      ENDDO

      DO I=1,METYPE
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  IF ((ELMTYP(I)(1:5) == 'TRIA3') .OR. (ELMTYP(I)(1:5) == 'QUAD4') .OR. (ELMTYP(I)(1:5) == 'QUAD8') .OR.        &
                      (ELMTYP(I)(1:5) == 'SHEAR') .OR. (ELMTYP(I)(1:6) == 'USERIN')) THEN
                     IF ((FORC_LOC == 'CORNER  ') .OR. (FORC_LOC == 'GAUSS   ') .OR. (ELMTYP(I)(1:5) == 'QUAD8')) THEN
                        NUM_PTS_ELFE(I) = NUM_SEi(I)
                     ELSE
                        NUM_PTS_ELFE(I) = 1
                     ENDIF
                     IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT)) > 0) NELREQ_ELFE(I) = NELREQ_ELFE(I) + NUM_PTS_ELFE(I)
                  ENDIF

                  IF ((ELMTYP(I)(1:5) == 'TRIA3') .OR. (ELMTYP(I)(1:5) == 'QUAD4') .OR. (ELMTYP(I)(1:5) == 'QUAD8') .OR.        &
                      (ELMTYP(I)(1:5) == 'SHEAR')) THEN
                     IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ') .OR. (ELMTYP(I)(1:4) == 'HEXA') .OR.            &
                         (ELMTYP(I)(1:5) == 'PENTA') .OR. (ELMTYP(I)(1:5) == 'TETRA') .OR. (ELMTYP(I)(1:5) == 'QUAD8')) THEN
                        NUM_PTS_STRE(I) = NUM_SEi(I)
                     ELSE
                        NUM_PTS_STRE(I) = 1
                     ENDIF
                     IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT)) > 0) NELREQ_STRE(I) = NELREQ_STRE(I) + NUM_PTS_STRE(I)
                  ENDIF

                  IF ((ELMTYP(I)(1:5) == 'TRIA3') .OR. (ELMTYP(I)(1:5) == 'QUAD4') .OR. (ELMTYP(I)(1:5) == 'QUAD8') .OR.        &
                      (ELMTYP(I)(1:5) == 'SHEAR')) THEN
                     IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ') .OR. (ELMTYP(I)(1:4) == 'HEXA') .OR.            &
                         (ELMTYP(I)(1:5) == 'PENTA') .OR. (ELMTYP(I)(1:5) == 'TETRA') .OR. (ELMTYP(I)(1:5) == 'QUAD8')) THEN
                        NUM_PTS_STRN(I) = NUM_SEi(I)
                     ELSE
                        NUM_PTS_STRN(I) = 1
                     ENDIF
                     IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT)) > 0) NELREQ_STRN(I) = NELREQ_STRN(I) + NUM_PTS_STRN(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO I=1,METYPE
         IF ((NELREQ_ELFE(I) == 0) .AND. (NELREQ_STRE(I) == 0) .AND. (NELREQ_STRN(I) == 0)) CYCLE

         NUM_OGEL_ELFE = 0
         NUM_OGEL_STRE = 0
         NUM_OGEL_STRN = 0

         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS /= 'N') CYCLE
            IF (ETYPE(J) /= ELMTYP(I)) CYCLE

            TYPE = ETYPE(J)
            EID  = EDAT(EPNT(J))

            IF ((IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT)) > 0) .OR.                                       &
                (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT)) > 0) .OR.                                       &
                (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT)) > 0)) THEN

               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0

               IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT)) > 0) THEN
                  OPT(4) = 'Y'
               ELSE
                  OPT(4) = 'N'
               ENDIF

               CALL EMG ( J, OPT, 'N', SUBR_NAME, 'N' )
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS

               IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT)) > 0) THEN
                  CALL HANDLE_FORCE_ELEMENT ( J, I, NUM_OGEL_ELFE, NUM_PTS_ELFE, OGEL_ELFE_LOC, EID_ELFE_LOC, GID_ELFE_LOC )
               ENDIF

               IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT)) > 0) THEN
                  CALL HANDLE_STRESS_ELEMENT ( J, I, NUM_OGEL_STRE, NUM_PTS_STRE, OGEL_STRE_LOC, EID_STRE_LOC, GID_STRE_LOC,     &
                                               POLY_FIT_ERR_STRE_LOC, POLY_FIT_ERR_INDEX_STRE_LOC )
               ENDIF

               IF (IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT)) > 0) THEN
                  CALL HANDLE_STRAIN_ELEMENT ( J, I, NUM_OGEL_STRN, NUM_PTS_STRN, OGEL_STRN_LOC, EID_STRN_LOC, GID_STRN_LOC,     &
                                               POLY_FIT_ERR_STRN_LOC, POLY_FIT_ERR_INDEX_STRN_LOC )
               ENDIF
            ENDIF
         ENDDO

         TYPE = ELMTYP(I)

         IF (NUM_OGEL_ELFE > 0) THEN
            OGEL(1:NUM_OGEL_ELFE,1:MOGEL) = OGEL_ELFE_LOC(1:NUM_OGEL_ELFE,1:MOGEL)
            EID_OUT_ARRAY(1:NUM_OGEL_ELFE,1:2) = EID_ELFE_LOC(1:NUM_OGEL_ELFE,1:2)
            GID_OUT_ARRAY(1:NUM_OGEL_ELFE,1:MELGP+1) = GID_ELFE_LOC(1:NUM_OGEL_ELFE,1:MELGP+1)
            CALL WRITE_ELEM_ENGR_FORCE ( JVEC, NUM_OGEL_ELFE, IHDR, NUM_PTS_ELFE(I), ITABLE )
         ENDIF

         IF (NUM_OGEL_STRE > 0) THEN
            OGEL(1:NUM_OGEL_STRE,1:MOGEL) = OGEL_STRE_LOC(1:NUM_OGEL_STRE,1:MOGEL)
            EID_OUT_ARRAY(1:NUM_OGEL_STRE,1:2) = EID_STRE_LOC(1:NUM_OGEL_STRE,1:2)
            GID_OUT_ARRAY(1:NUM_OGEL_STRE,1:MELGP+1) = GID_STRE_LOC(1:NUM_OGEL_STRE,1:MELGP+1)
            POLY_FIT_ERR(1:NUM_OGEL_STRE) = POLY_FIT_ERR_STRE_LOC(1:NUM_OGEL_STRE)
            POLY_FIT_ERR_INDEX(1:NUM_OGEL_STRE) = POLY_FIT_ERR_INDEX_STRE_LOC(1:NUM_OGEL_STRE)
            CALL WRITE_ELEM_STRESSES ( JVEC, NUM_OGEL_STRE, IHDR, NUM_PTS_STRE(I), ITABLE )
         ENDIF

         IF (NUM_OGEL_STRN > 0) THEN
            OGEL(1:NUM_OGEL_STRN,1:MOGEL) = OGEL_STRN_LOC(1:NUM_OGEL_STRN,1:MOGEL)
            EID_OUT_ARRAY(1:NUM_OGEL_STRN,1:2) = EID_STRN_LOC(1:NUM_OGEL_STRN,1:2)
            GID_OUT_ARRAY(1:NUM_OGEL_STRN,1:MELGP+1) = GID_STRN_LOC(1:NUM_OGEL_STRN,1:MELGP+1)
            POLY_FIT_ERR(1:NUM_OGEL_STRN) = POLY_FIT_ERR_STRN_LOC(1:NUM_OGEL_STRN)
            POLY_FIT_ERR_INDEX(1:NUM_OGEL_STRN) = POLY_FIT_ERR_INDEX_STRN_LOC(1:NUM_OGEL_STRN)
            CALL WRITE_ELEM_STRAINS ( JVEC, NUM_OGEL_STRN, IHDR, NUM_PTS_STRN(I), ITABLE )
         ENDIF
      ENDDO

      DEALLOCATE ( OGEL_ELFE_LOC, OGEL_STRE_LOC, OGEL_STRN_LOC )
      DEALLOCATE ( EID_ELFE_LOC, EID_STRE_LOC, EID_STRN_LOC )
      DEALLOCATE ( GID_ELFE_LOC, GID_STRE_LOC, GID_STRN_LOC )
      DEALLOCATE ( POLY_FIT_ERR_STRE_LOC, POLY_FIT_ERR_STRN_LOC )
      DEALLOCATE ( POLY_FIT_ERR_INDEX_STRE_LOC, POLY_FIT_ERR_INDEX_STRN_LOC )

 9999 CONTINUE

      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
      CONTAINS

      SUBROUTINE HANDLE_FORCE_ELEMENT ( J, I, NUM_OGEL, NUM_PTS, OGEL_LOC, EID_LOC, GID_LOC )

      INTEGER(LONG), INTENT(IN)       :: J, I
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS(METYPE)
      INTEGER(LONG), INTENT(INOUT)    :: NUM_OGEL
      REAL(DOUBLE), INTENT(INOUT)     :: OGEL_LOC(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: EID_LOC(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: GID_LOC(:,:)

      INTEGER(LONG)                   :: M, K, NUM_PREV, NUM_TMP
      REAL(DOUBLE)                    :: STRESS_RAW_LOC(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRESS_OUT_LOC(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: PCT_ERR_MAX_LOC, TEL_LOC(3,3)

      NUM_PREV = NUM_OGEL
      STRESS_RAW_LOC = ZERO
      STRESS_OUT_LOC = ZERO
      PCT_ERR_MAX_LOC = ZERO

      IF ((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J)(1:5) == 'QUAD8') .OR.                  &
          (ETYPE(J)(1:5) == 'SHEAR')) THEN

         DO M=1,NUM_PTS(I)
            CALL ELEM_STRE_STRN_ARRAYS ( M )
            STRESS_RAW_LOC(:,M) = STRESS(:)
         ENDDO

         STRESS_OUT_LOC(:,1) = STRESS_RAW_LOC(:,1)

         IF ((FORC_LOC == 'CORNER  ') .OR. (ETYPE(J)(1:5) == 'QUAD8')) THEN
            IF (TYPE(1:5) == 'QUAD4') THEN
               CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW_LOC, 9, NUM_PTS(I), STRESS_OUT_LOC, STRESS_OUT_PCT_ERR,                &
                                            STRESS_OUT_ERR_INDEX, PCT_ERR_MAX_LOC )
            ELSE IF (ETYPE(J)(1:5) == 'QUAD8') THEN
               CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW_LOC, 9, NUM_PTS(I), STRESS_OUT_LOC, STRESS_OUT_PCT_ERR,                &
                                            STRESS_OUT_ERR_INDEX, PCT_ERR_MAX_LOC )
               DO M=2,NUM_PTS(I)
                  CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE(M), TEL_LOC, '' )
                  CALL TRANSFORM_SHELL_STR( TEL_LOC, STRESS_OUT_LOC(:,M), ONE )
               ENDDO
               STRESS_OUT_LOC(:,1) = (STRESS_OUT_LOC(:,2) + STRESS_OUT_LOC(:,3) + STRESS_OUT_LOC(:,4) + STRESS_OUT_LOC(:,5)) / FOUR
            ENDIF
         ENDIF

         DO M=1,NUM_PTS(I)
            STRESS(:) = STRESS_OUT_LOC(:,M)
            NUM_TMP = NUM_OGEL
            CALL SHELL_ENGR_FORCE_OGEL ( NUM_OGEL )
            OGEL_LOC(NUM_OGEL,:) = OGEL(NUM_OGEL,:)
            EID_LOC(NUM_OGEL,1) = EID
            EID_LOC(NUM_OGEL,2) = ZERO_I
            GID_LOC(NUM_OGEL,1) = ZERO_I
            DO K=1,ELGP
               GID_LOC(NUM_OGEL,K+1) = AGRID(K)
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE HANDLE_FORCE_ELEMENT

      SUBROUTINE HANDLE_STRESS_ELEMENT ( J, I, NUM_OGEL, NUM_PTS, OGEL_LOC, EID_LOC, GID_LOC, POLY_LOC, POLY_IDX_LOC )

      INTEGER(LONG), INTENT(IN)       :: J, I
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS(METYPE)
      INTEGER(LONG), INTENT(INOUT)    :: NUM_OGEL
      REAL(DOUBLE), INTENT(INOUT)     :: OGEL_LOC(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: EID_LOC(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: GID_LOC(:,:)
      REAL(DOUBLE), INTENT(INOUT)     :: POLY_LOC(:)
      INTEGER(LONG), INTENT(INOUT)    :: POLY_IDX_LOC(:)

      INTEGER(LONG)                   :: M, K, NUM_PREV, NUM_TMP
      REAL(DOUBLE)                    :: STRESS_RAW_LOC(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRESS_OUT_LOC(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: PCT_ERR_MAX_LOC, TEL_LOC(3,3)

      NUM_PREV = NUM_OGEL
      STRESS_RAW_LOC = ZERO
      STRESS_OUT_LOC = ZERO
      PCT_ERR_MAX_LOC = ZERO

      DO M=1,NUM_PTS(I)
         CALL ELEM_STRE_STRN_ARRAYS ( M )
         STRESS_RAW_LOC(:,M) = STRESS(:)
      ENDDO

      STRESS_OUT_LOC(:,1) = STRESS_RAW_LOC(:,1)

      IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ') .OR. (ETYPE(J)(1:4) == 'HEXA') .OR.                    &
          (ETYPE(J)(1:5) == 'PENTA') .OR. (ETYPE(J)(1:5) == 'TETRA') .OR. (ETYPE(J)(1:5) == 'QUAD8')) THEN
         IF (TYPE(1:5) == 'QUAD4') THEN
            CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW_LOC, 9, NUM_PTS(I), STRESS_OUT_LOC, STRESS_OUT_PCT_ERR,                &
                                         STRESS_OUT_ERR_INDEX, PCT_ERR_MAX_LOC )
         ELSE IF (ETYPE(J)(1:5) == 'QUAD8') THEN
            CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW_LOC, 9, NUM_PTS(I), STRESS_OUT_LOC, STRESS_OUT_PCT_ERR,                &
                                         STRESS_OUT_ERR_INDEX, PCT_ERR_MAX_LOC )
            DO M=2,NUM_PTS(I)
               CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE(M), TEL_LOC, '' )
               CALL TRANSFORM_SHELL_STR( TEL_LOC, STRESS_OUT_LOC(:,M), ONE )
            ENDDO
            STRESS_OUT_LOC(:,1) = (STRESS_OUT_LOC(:,2) + STRESS_OUT_LOC(:,3) + STRESS_OUT_LOC(:,4) + STRESS_OUT_LOC(:,5)) / FOUR
         ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
            STRESS_OUT_LOC(:,:) = STRESS_RAW_LOC(:,:)
         ENDIF
      ENDIF

      DO M=1,NUM_PTS(I)
         STRESS(:) = STRESS_OUT_LOC(:,M)
         NUM_TMP = NUM_OGEL
         CALL CALC_ELEM_STRESSES ( MAXREQ, NUM_OGEL, ZERO_I, 'Y', 'N' )
         DO K=NUM_PREV+1,NUM_OGEL
            OGEL_LOC(K,:) = OGEL(K,:)
            EID_LOC(K,1) = EID
            EID_LOC(K,2) = ZERO_I
            GID_LOC(K,1) = ZERO_I
            DO R=1,ELGP
               GID_LOC(K,R+1) = AGRID(R)
            ENDDO
            POLY_LOC(K) = ZERO
            POLY_IDX_LOC(K) = ZERO_I
            IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) THEN
               IF (TYPE(1:5) == 'QUAD4') THEN
                  POLY_LOC(K) = STRESS_OUT_PCT_ERR(M)
                  POLY_IDX_LOC(K) = STRESS_OUT_ERR_INDEX(M)
               ENDIF
            ENDIF
         ENDDO
         NUM_PREV = NUM_OGEL
      ENDDO

      END SUBROUTINE HANDLE_STRESS_ELEMENT

      SUBROUTINE HANDLE_STRAIN_ELEMENT ( J, I, NUM_OGEL, NUM_PTS, OGEL_LOC, EID_LOC, GID_LOC, POLY_LOC, POLY_IDX_LOC )

      INTEGER(LONG), INTENT(IN)       :: J, I
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS(METYPE)
      INTEGER(LONG), INTENT(INOUT)    :: NUM_OGEL
      REAL(DOUBLE), INTENT(INOUT)     :: OGEL_LOC(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: EID_LOC(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: GID_LOC(:,:)
      REAL(DOUBLE), INTENT(INOUT)     :: POLY_LOC(:)
      INTEGER(LONG), INTENT(INOUT)    :: POLY_IDX_LOC(:)

      INTEGER(LONG)                   :: M, K, NUM_PREV, NUM_TMP
      REAL(DOUBLE)                    :: STRAIN_RAW_LOC(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: STRAIN_OUT_LOC(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: PCT_ERR_MAX_LOC, TEL_LOC(3,3)

      NUM_PREV = NUM_OGEL
      STRAIN_RAW_LOC = ZERO
      STRAIN_OUT_LOC = ZERO
      PCT_ERR_MAX_LOC = ZERO

      DO M=1,NUM_PTS(I)
         CALL ELEM_STRE_STRN_ARRAYS ( M )
         STRAIN_RAW_LOC(:,M) = STRAIN(:)
      ENDDO

      STRAIN_OUT_LOC(:,1) = STRAIN_RAW_LOC(:,1)

      IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ') .OR. (ETYPE(J)(1:4) == 'HEXA') .OR.                    &
          (ETYPE(J)(1:5) == 'PENTA') .OR. (ETYPE(J)(1:5) == 'TETRA') .OR. (ETYPE(J)(1:5) == 'QUAD8')) THEN
         IF (TYPE(1:5) == 'QUAD4') THEN
            CALL POLYNOM_FIT_STRE_STRN ( STRAIN_RAW_LOC, 9, NUM_PTS(I), STRAIN_OUT_LOC, STRAIN_OUT_PCT_ERR,                &
                                         STRAIN_OUT_ERR_INDEX, PCT_ERR_MAX_LOC )
         ELSE IF (ETYPE(J)(1:5) == 'QUAD8') THEN
            CALL POLYNOM_FIT_STRE_STRN ( STRAIN_RAW_LOC, 9, NUM_PTS(I), STRAIN_OUT_LOC, STRAIN_OUT_PCT_ERR,                &
                                         STRAIN_OUT_ERR_INDEX, PCT_ERR_MAX_LOC )
            DO M=2,NUM_PTS(I)
               CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE(M), TEL_LOC, '' )
               CALL TRANSFORM_SHELL_STR( TEL_LOC, STRAIN_OUT_LOC(:,M), TWO )
            ENDDO
            STRAIN_OUT_LOC(:,1) = (STRAIN_OUT_LOC(:,2) + STRAIN_OUT_LOC(:,3) + STRAIN_OUT_LOC(:,4) + STRAIN_OUT_LOC(:,5)) / FOUR
         ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
            STRAIN_OUT_LOC(:,:) = STRAIN_RAW_LOC(:,:)
         ENDIF
      ENDIF

      DO M=1,NUM_PTS(I)
         STRAIN(:) = STRAIN_OUT_LOC(:,M)
         NUM_TMP = NUM_OGEL
         CALL CALC_ELEM_STRAINS ( MAXREQ, NUM_OGEL, ZERO_I, 'Y', 'N' )
         DO K=NUM_PREV+1,NUM_OGEL
            OGEL_LOC(K,:) = OGEL(K,:)
            EID_LOC(K,1) = EID
            EID_LOC(K,2) = ZERO_I
            GID_LOC(K,1) = ZERO_I
            DO R=1,ELGP
               GID_LOC(K,R+1) = AGRID(R)
            ENDDO
            POLY_LOC(K) = ZERO
            POLY_IDX_LOC(K) = ZERO_I
            IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
               IF (TYPE(1:5) == 'QUAD4') THEN
                  POLY_LOC(K) = STRAIN_OUT_PCT_ERR(M)
                  POLY_IDX_LOC(K) = STRAIN_OUT_ERR_INDEX(M)
               ENDIF
            ENDIF
         ENDDO
         NUM_PREV = NUM_OGEL
      ENDDO

      END SUBROUTINE HANDLE_STRAIN_ELEMENT

      END SUBROUTINE OFP3_SHELL_NO_PCOMP
