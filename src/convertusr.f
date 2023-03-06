C ============================================================================

      subroutine savespecial(x,w)

C     Cut out a rectangular block from a data file.
C     The block may have lower dimensionality than the original data.

      include 'convertdef.f'

      REAL*8 x(nx1,nx2,nx3,ndim), w(nx1,nx2,nx3,nw)
      LOGICAL fileopen
      INTEGER idim,ix1,ix2,ix3,iw,ieqpar,i,l
      INTEGER ncutdim,icutdim(3),ixcutmin(3),ixcutmax(3)
      CHARACTER*79 fileheadout,varnamesout
 
      data ncutdim,icutdim /4*0/, ixcutmin,ixcutmax /6*1/
      save fileheadout,varnamesout
C ----------------------------------------------------------------------------

      if(ncutdim.eq.0)then

C        Read index limits for cutting block

         write(*,*)'Give ixcutmin,ixcutmax for each dimensions:'
         read(*,*)(ixcutmin(idim),ixcutmax(idim),idim=1,ndim)

C        Check index limits and dimensionality of block

         ncutdim=0
         do idim=1,ndim
            if(ixcutmin(idim).lt.1) stop 'ixcutmin must be >= 1 !'
            if(ixcutmax(idim).gt.nx(idim)) stop 'ixcutmax must be <= nx'
            if(ixcutmin(idim).gt.ixcutmax(idim)) 
     &           stop 'ixcutmin must be <= ixcutmax'

            if(ixcutmin(idim).lt.ixcutmax(idim))then

C              Non-degenerate dimension

               ncutdim=ncutdim+1
               icutdim(ncutdim)=idim
            endif
         enddo
         
         if(ncutdim.eq.0)
     &   stop 'ixcutmax must be > ixcutmin at least in one dimension!'
         
         write(*,*)'ncutdim, icutdim:',
     &        ncutdim,(icutdim(idim),idim=1,ncutdim)

C        Determine fileheadout and varnamesout

         if(ncutdim.eq.ndim)then
            fileheadout=fileheadini
            varnamesout=varnames
         else
C           Correct number of dimensions for fileheadout

            i=index(fileheadini,'    ')
            write(fileheadout,'(a,i1,a)')
     &           fileheadini(1:i-3),ncutdim,fileheadini(i-1:i-1)


C           Remove coordinate names from varnamesini
C           and append to varnamesout for non-degenerate dimensions only

            l=0
            do idim=1,ndim

C              Get rid of leading spaces and 
C              find position after the first coordinate name
 10            i=index(varnames,' ')
               if(i.eq.1)then
                  varnames=varnames(2:79)
                  goto 10
               endif

               if(ixcutmin(idim).lt.ixcutmax(idim))then

C                 Append coordinate name to varnamesout
                  if(l.gt.0)then
                      varnamesout=varnamesout(1:l)//varnames(1:i)
                  else
                      varnamesout=varnames(1:i)
                  endif
                  l=l+i
               endif

C              Remove coordinate name
               varnames=varnames(i+1:79)
            enddo

C           Append rest of variable names
            varnamesout=varnamesout(1:l)//varnames(1:79-l)

            write(*,*)'fileheadout = ',fileheadout(1:54)
            write(*,*)'varnamesout = ',varnamesout(1:54)
         endif
      endif

C Check if the file needs to be opened

      inquire(unitout,opened=fileopen)

C Open an unformatted file

      if(.not.fileopen)
     &   open(unitout,FILE=filenameout,STATUS='unknown',
     &        FORM='unformatted')

C Save header info

      write(unitout)fileheadout
      write(unitout)it,t,ncutdim,neqpar,nw
      write(unitout)
     &     (ixcutmax(icutdim(idim))-ixcutmin(icutdim(idim))+1,
     &     idim=1,ncutdim)
      write(unitout)(eqpar(ieqpar),ieqpar=1,neqpar)
      write(unitout)varnamesout

C Save non-degenerate coordinates
      write(unitout)((((x(ix1,ix2,ix3,icutdim(idim)),
     &     ix1=ixcutmin(1),ixcutmax(1)),
     &     ix2=ixcutmin(2),ixcutmax(2)),
     &     ix3=ixcutmin(3),ixcutmax(3)),
     &     idim=1,ncutdim)


C Save conservative variables
      do iw=1,nw
         write(unitout)((((w(ix1,ix2,ix3,iw)),
     &     ix1=ixcutmin(1),ixcutmax(1)),
     &     ix2=ixcutmin(2),ixcutmax(2)),
     &     ix3=ixcutmin(3),ixcutmax(3))
      enddo

      return
      end

C ============================================================================
