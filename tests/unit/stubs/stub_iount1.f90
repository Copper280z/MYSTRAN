! Minimal stub for MODULE IOUNT1 providing only what unit-tested subroutines need.
! WRT_LOG = 0 means all logging guards (IF WRT_LOG >= SUBR_BEGEND) evaluate false,
! so F04 is never written to and OURTIM is never called.

      MODULE IOUNT1

      USE PENTIUM_II_KIND, ONLY  :  LONG

      IMPLICIT NONE

      SAVE

      ! File unit for the log (F04) — stdout used here so any accidental writes are visible
      INTEGER(LONG) :: F04     = 6

      ! Logging verbosity level; 0 disables all subroutine entry/exit logging
      INTEGER(LONG) :: WRT_LOG = 0

      END MODULE IOUNT1
