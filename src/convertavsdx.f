C ============================================================================

      subroutine avsdxconvert(x,w,MDAT,WORKLOVER)

      INCLUDE '../AVSDX/include/lover.h'
      include 'convertdef.f'

      REAL*8 x(nx1,nx2,nx3,3), w(nx1,nx2,nx3,nw)
      REAL*8 MDAT(nw,nx1,nx2,nx3),WORKLOVER(4*nw*nxs)

      REAL*8     TIME(MAXNTIME)

      INTEGER    IWORK(MAXNTIME*MAXNFIELDS)                                
      INTEGER    RI(4), CO(4), SL(4), ST(MAXNFIELDS+1)                     
      INTEGER    FORM                                             

      CHARACTER  FILENAME*23, FORWHICH
      CHARACTER  FN(MAXNFIELDS+1)*23, TY(MAXNFIELDS+1)*8

      integer ix1,ix2,ix3,idim,iw,isave,nsave

      character*10 physics,typephys

      integer i,irho,ib,im,ie,nfields,ndims

C ---------------------------------------------------------------------------

C     Determine physics from headline containing sg. like '_mhd12'

      i=index(fileheadini,'_')
      if(i.le.0)then
         write(*,*)'Physics is not defined in file header'
         write(*,*)'Physics (rho/hdadiab/hd/mhdiso/mhd):'
         read(*,*)typephys
      else
         physics=fileheadini(i+1:i+10)
         i=index(physics,' ')
         typephys=physics(1:i-3)
         write(*,'(a,a)')' physics = ',typephys
      endif

      irho=1
      im=1
      ie=0
      ib=0
      if(typephys.eq.'rho')then
         im=0
      else if(typephys.eq.'hdadiab')then
         im=1
      else if(typephys.eq.'hd')then
         ie=1
      else if(typephys.eq.'mhdiso')then
         ib=1
      else if(typephys.eq.'mhd')then
         ie=1
         ib=1
      else
         stop 'Error: Unknown physics'
      endif

      if(typephys.ne.'rho')ndir=(nw-irho-ie)/(im+ib)
      write(*,'(a,i4,a,i4,a,i4,a,i1,a,i1)')
     &     ' nx1=',nx1,' nx2=',nx2,' nx3=',nx3,
     &     ' ndim=',ndim,' ndir=',ndir

      nfields = irho+im+ie+ib

      if(typefileini.eq.'ascii')then
         open(unitini,FILE=filenameini,STATUS='old')
      else
         open(unitini,FILE=filenameini,STATUS='old',
     &        FORM='unformatted')
      endif

      FILENAME  = filenameout

      if(typefileout.eq.'AVSbin')then
         FORWHICH='A'
         FORM    = VISFORMBIN
         NDIMS=NDIM+1
      else if (typefileout.eq.'AVSxdr')then
         FORWHICH='A'
         FORM    = VISFORMXDR
         NDIMS=NDIM+1
      else if (typefileout.eq.'DXbin')then
         FORWHICH='D'
         FORM    = VISFORMBIN
         NDIMS=NDIM
         typefileout='binary'
      else if (typefileout.eq.'DXasc')then
         FORWHICH='D'
         FORM    = VISFORMASC
         NDIMS=NDIM
      else
         stop 'Unknown file format'
      endif

C     * On the first pass through the raw data file,                   
C       the number of time steps and the time labels are read.          

      nsave = 0                                                         
 100  continue
         call readfileini(x,w)
         if(ios.lt.0) goto 200
         nsave = nsave + 1                                              
         TIME(nsave)=t
         goto 100
 200  continue

C     To be safe we put 0-s in unused components of x

      do idim=ndim+1,3
         do ix3=1,nx3
            do ix2=1,nx2
               do ix1=1,nx1
                  x(ix1,ix2,ix3,idim)=0.D0
               enddo
            enddo
         enddo
      enddo

C     * end-of-file                                                     

      RI(1) = nx1
      RI(2) = nx2
      RI(3) = nx3
      RI(4) = nsave                                                     
      CO(1) = 1
      CO(2) = 2
      CO(3) = 3
      CO(4) = 4
      IF (FORWHICH.eq.'A' .and. nsave.GT.1) CO(NDIMS) = 4

C     Here iw is an index for scalar/vector fields

      iw = 0                                                          
      if (irho.ne.0) then
         iw = iw + 1                                                
         ST(iw) = VISRANKSCA                                          
         TY(iw) = VISTYPEDOU                                          
         FN(iw) = 'rho'                                               
      endif
      if (im.ne.0) then
         iw = iw + 1                                                
         ST(iw) = ndir
         TY(iw) = VISTYPEDOU                                          
         FN(iw) = 'M'                                                 
      endif                                      
      if (ie.ne.0) then
         iw = iw + 1                                                
         ST(iw) = VISRANKSCA                                          
         TY(iw) = VISTYPEDOU                                          
         FN(iw) = 'e'                                                 
      endif
      if (ib.ne.0) then                                                 
         iw = iw + 1                                                
         ST(iw) = ndir
         TY(iw) = VISTYPEDOU                                          
         FN(iw) = 'B'                                                 
      endif

C     * Setting arguments                                               
      CALL VISSETARGS(nw,nx1,nx2,nx3,nsave,RI,CO,SL,NDIMS,
     &                .TRUE.,NFIELDS,FORM,VISGEOMAUX)

C     * We will output one timestep a time                              
      CALL VISONESTEP                                                   

C     * Producing DX/AVS internal file (header part)                    
      CALL VISGOFORIT(FORWHICH,VISINTERNA,
     &     MDAT,x(1,1,1,1),x(1,1,1,2),x(1,1,1,3),
     &     TIME,FILENAME,FN,ST,TY,WORKLOVER,IWORK)

C     * On the second pass through the raw data file,                  
C       the actual data are written in DX/AVS format.                   

      rewind(unitini)

      DO isave = 1, nsave

         call readfileini(x,w)
         write(*,'(a,i4)')' Converting snapshot ',ipict

C        Reordering w of VAC to MDAT of LOVER
         do iw=1,nw
            do ix3=1,nx3
               do ix2=1,nx2
                  do ix1=1,nx1
                     MDAT(iw,ix1,ix2,ix3)=w(ix1,ix2,ix3,iw)
                  enddo
               enddo
            enddo
         enddo

C        * Producing DX/AVS internal file                   

         CALL VISGOFORIT(FORWHICH,VISINTERNA,
     &        MDAT,x(1,1,1,1),x(1,1,1,2),x(1,1,1,3),
     &        TIME,FILENAME,FN,ST,TY,WORKLOVER,IWORK)

      enddo                                                             


      return
      end
