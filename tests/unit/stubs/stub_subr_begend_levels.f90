! Minimal stub for MODULE SUBR_BEGEND_LEVELS providing only the constants needed
! by unit-tested subroutines.  Values match the real module (LINK_BEGEND + N).

      MODULE SUBR_BEGEND_LEVELS

      USE PENTIUM_II_KIND, ONLY  :  LONG

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER :: LINK_BEGEND           =  1
      INTEGER(LONG), PARAMETER :: CROSS_BEGEND          = LINK_BEGEND + 10  ! = 11
      INTEGER(LONG), PARAMETER :: MATMULT_FFF_BEGEND    = LINK_BEGEND +  8  ! =  9
      INTEGER(LONG), PARAMETER :: MATMULT_FFF_T_BEGEND  = LINK_BEGEND +  8  ! =  9

      END MODULE SUBR_BEGEND_LEVELS
