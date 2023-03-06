C ============================================================================
C This dummy subroutine is here if the AVSDX library is not used
C ============================================================================

      subroutine avsdxconvert(x,w,MDAT,WORKLOVER)

      REAL*8 x,w,MDAT,WORKLOVER

      write(*,*)'Transformation to AVS and DX formats is switched off.'
      write(*,*)'Edit Makefile and recompile convertdata.'      

      stop
      end
