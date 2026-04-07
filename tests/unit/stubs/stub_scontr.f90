! Minimal stub for MODULE SCONTR providing only what unit-tested subroutines need.
! The real SCONTR has hundreds of variables; this stub provides just the symbols
! required so that CROSS, MATMULT_FFF, and MATMULT_FFF_T compile cleanly.

      MODULE SCONTR

      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG

      IMPLICIT NONE

      SAVE

      ! Used as the length of SUBR_NAME character variables
      CHARACTER(31*BYTE), PARAMETER :: BLNK_SUB_NAM = '                               '

      ! Used by IOUNT1 as array dimensions
      INTEGER(LONG), PARAMETER :: MBUG =  10
      INTEGER(LONG), PARAMETER :: MFIJ =   5

      END MODULE SCONTR
