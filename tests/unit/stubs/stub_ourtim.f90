! No-op stub for SUBROUTINE OURTIM.
! The real OURTIM records wall-clock time into TIMDAT module variables.
! In unit tests WRT_LOG = 0, so OURTIM is never called; this stub only
! satisfies the linker.

      SUBROUTINE OURTIM

      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG
      USE TIMDAT, ONLY           :  HOUR, MINUTE, SEC, SFRAC, TSEC, DSEC

      IMPLICIT NONE

      ! Nothing to do — logging is disabled in unit tests
      RETURN

      END SUBROUTINE OURTIM
