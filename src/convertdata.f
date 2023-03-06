      PROGRAM CONVERT

C     Reads an ASCII or BINARY VAC data file and converts it to 
C
C     1. VAC ASCII format
C     2. VAC BINARY format
C     3. AVS internal file format
C        ascii/xdr
C     4. DX  internal file format
C        ascii/binary
C
C     Original version of a program, which converts to AVS and DX was provided 
C     by A.J.C. Belien.
C
C     Further adaptation to VAC by G. Toth with A.J.C. Belien's permission.
C
C     The LOVER library itself is A.J.C. Belien's work.

C ----------------------------------------------------------------------------

      include 'convertdef.f'

      integer i,neededsize
      character*79 pictstr

C ----------------------------------------------------------------------------

      write(*,'(a)')'Input file name:'
      read(*,'(a)')filenameini

      call getfiletype

      if(verbose.ge.2)then
         write(*,*)'headline      =',fileheadini(1:40)
         write(*,*)'it, t         =',it,t
         write(*,*)'ndim,neqpar,nw=',signedndim,neqpar,nw
         write(*,*)'nx            =',nx
         write(*,*)'eqpar         =',(eqpar(i),i=1,neqpar)
         write(*,*)'varnames      =',varnames(1:40)
         write(*,*)
      endif
      if(verbose.ge.1)then
         write(*,*)'Input file type = ',typefileini
         write(*,*)'Snapshot size   =',pictsize,' bytes'
         write(*,*)
      endif

      if(pictsize .gt. worksize*8)then
         write(*,*)'Error: input snapshot size is too big!'
         write(*,*)'Adjust worksize in src/convertdef.f to at least',
     &        (pictsize+7)/8
         stop
      endif

      write(*,'(a)')'Which snapshots should be read?'
      if(verbose.ge.1)then
         write(*,*)'example string: 1,3,4,6:10:2,13:100'
         write(*,*)'hit RETURN for the whole file'
      endif

      read(*,'(a)')pictstr
      call parsepictstr(pictstr,npictmax,readpict,npict)

      if(verbose.gt.1.and.npict.lt.20)then
         write(*,*)'npict    = ',npict
         write(*,*)'readpict = ',(readpict(ipict),ipict=1,20)
      endif

      write(*,'(a)')'Output file format ',
     &   '(ascii/binary/real4/AVSbin/AVSxdr/DXbin/DXasc/special):'
      read(*,'(a)')typefileout

C     Check if we have enough work space

      if(typefileout.eq.'ascii'.or.typefileout.eq.'binary')then
         neededsize=(ndim+nw)*nxs + 1
      else if(typefileout.eq.'real4')then
         neededsize=(ndim+nw)*nxs + (ndim*nxs)/2
      else if(typefileout.eq.'AVSbin'.or.typefileout.eq.'AVSxdr'.or.
     &        typefileout.eq.'DXbin' .or.typefileout.eq.'DXasc')then
         neededsize=(3+nw+nw+4*nw)*nxs
      else if(typefileout.eq.'special') then
         neededsize=1
      else
         stop 'Error: Unknown file format'
      endif

      if(neededsize .gt. worksize)then
         write(*,*)'Error: output snapshot size is too big!'
         write(*,*)'Adjust worksize in src/convertdef.f to at least',
     &        neededsize
         stop
      endif

      ipict=0

      if(typefileout.eq.'ascii'.or.typefileout.eq.'binary'.or.
     &   typefileout.eq.'real4'.or.typefileout.eq.'special')then

         write(*,'(a)')
     &        'Output file name (different from input file name):'
         read(*,'(a)')filenameout

         if(filenameini .eq. filenameout)
     &        stop 'Error: Input and output file names have to differ!'

C        Use work array to contain x, w, xreal4
         call vacconvert(work(1),work(ndim*nxs+1),work((ndim+nw)*nxs+1))
      else
         write(*,'(a)')
     &        'Output file name (without the .dx or .fld extension):'
         read(*,'(a)')filenameout

C        Use work array to contain x, w, MDAT, WORKLOVER
         call avsdxconvert(work(1),work(3*nxs+1),
     &        work((3+nw)*nxs+1),work((3+nw+nw)*nxs+1))
      endif

      write(*,'(a)')'Data conversion finished.'

      stop
      end

C ============================================================================

      subroutine getfiletype

C     Check if the file is ASCII, BINARY, or REAL4
C     Read in header information and calculate size of a snapshot

      include 'convertdef.f'

      logical fileexist
      integer i,iosauto
      character*91 test

C ----------------------------------------------------------------------------

      inquire(FILE=filenameini,EXIST=fileexist)
      if(.not.fileexist) stop 'Stop: Input file does not exist'

C     Figure out the input file type based on the first characters. 
C     The header consists of 79 alphabetic characters. The files contain:
C
C     ascii:  header \n ...
C     binary: 0 0 0 79 header 0 0 0 79 0 0 0 24 it t ndim neqpar nw 0 0 0 24
C             79 0 0 0 header 79 0 0 0 24 0 0 0 it t ndim neqpar nw 24 0 0 0
C     real4:  0 0 0 79 header 0 0 0 79 0 0 0 24 it t ndim neqpar nw 0 0 0 24 
C             79 0 0 0 header 79 0 0 0 20 0 0 0 it t ndim neqpar nw 20 0 0 0
C
C     The order of the zeros and the 79/24/20 depends on the endianness

      open(unitini,FILE=filenameini,STATUS='old')
      read(unitini,'(a91)',iostat=iosauto)test
      close(unitini)

      if(iosauto.ne.0)then
         typefileini='unknown'
      else
         typefileini='ascii'
         do i=1,79
            if(ichar(test(i:i)).eq.0)typefileini='unknown'
         enddo
         if(typefileini.eq.'unknown')then
            do i=88,91
               if(ichar(test(i:i)).eq.24)typefileini='binary'
               if(ichar(test(i:i)).eq.20)typefileini='real4'
            enddo
         endif
      endif

      if(typefileini.eq.'unknown')then
          write(*,*)'Could not determine file type automatically...'
          write(*,*)'typefileini (ascii/binary/real4)?'
          read(*,'(a)')typefileini
      endif

      if(typefileini.ne.'ascii'.and.typefileini.ne.'binary'.and.
     &   typefileini.ne.'real4')stop 'Unknown input file type!'

      do i=1,3
         nx(i)=1
      enddo

      if(typefileini.eq.'ascii')then
         open(unitini,FILE=filenameini,STATUS='old')
      else
         open(unitini,FILE=filenameini,STATUS='old',FORM='unformatted')
      endif

      call readheader

      close(unitini)

      nx1=nx(1)
      nx2=nx(2)
      nx3=nx(3)
      nxs=nx1*nx2*nx3

C     Calculate size of one snapshot

      if(typefileini.eq.'ascii')then
         pictsize=1+79 + 1+7+13+9 + 1+ndim*4 + 1+neqpar*13 + 1+79 + 
     &            (18*(ndim+nw)+1)*nxs
      else if(typefileini.eq.'binary')then
         pictsize=8+79 + 8+4*4+8  + 8+ndim*4 + 8+neqpar*8  + 8+79 +
     &            8*(1+nw)+8*(ndim+nw)*nxs
      else if(typefileini.eq.'real4')then
         pictsize=8+79 + 5*4+8  + 8+ndim*4 + 8+neqpar*4  + 8+79 +
     &            8*(1+nw)+4*(ndim+nw)*nxs
      endif

      return
      end

C ============================================================================

      subroutine readheader

C     Read in header information

      include 'convertdef.f'

      integer idim,ieqpar
      real tmpreal4, eqparreal4(100)
C ----------------------------------------------------------------------------

      if(typefileini.eq.'ascii')then

         read(unitini,'(a)',iostat=ios)fileheadini
         if(ios .lt. 0)return
         read(unitini,*)it,t,signedndim,neqpar,nw
         ndim=abs(signedndim)
         read(unitini,*)(nx(idim),idim=1,ndim)
         read(unitini,*)(eqpar(ieqpar),ieqpar=1,neqpar)
         read(unitini,'(a)')varnames

      else if(typefileini.eq.'binary')then

         read(unitini,iostat=ios)fileheadini
         if(ios .lt. 0)return
         read(unitini)it,t,signedndim,neqpar,nw
         ndim=abs(signedndim)
         read(unitini)(nx(idim),idim=1,ndim)
         read(unitini)(eqpar(ieqpar),ieqpar=1,neqpar)
         read(unitini)varnames

      else if(typefileini.eq.'real4')then

         read(unitini,iostat=ios)fileheadini
         if(ios .lt. 0)return
         read(unitini)it,tmpreal4,signedndim,neqpar,nw
         t=dble(tmpreal4)
         ndim=abs(signedndim)
         read(unitini)(nx(idim),idim=1,ndim)
         read(unitini)(eqparreal4(ieqpar),ieqpar=1,neqpar)
         do ieqpar=1,neqpar
            eqpar(ieqpar)=dble(eqparreal4(ieqpar))
         enddo
         read(unitini)varnames

      else
         stop 'Unknown input file type!'
      endif

      return
      end

C ============================================================================

      subroutine vacconvert(x,w,xreal4)

C     Convert between ASCII, BINARY, and REAL4 formats for VAC files

      include 'convertdef.f'

      REAL*8 x(nx1,nx2,nx3,ndim), w(nx1,nx2,nx3,nw)
      REAL   xreal4(nx1,nx2,nx3,ndim)
      integer ix1,ix2,ix3,iw

C ----------------------------------------------------------------------------

      if(typefileini.eq.'ascii')then
         open(unitini,FILE=filenameini,STATUS='old')
      else
         open(unitini,FILE=filenameini,STATUS='old',
     &        FORM='unformatted')
      endif

      if(typefileout.eq.'ascii')then
         open(unitout,FILE=filenameout,STATUS='unknown')
      else if(typefileout.eq.'binary'.or.typefileout.eq.'real4')then
         open(unitout,FILE=filenameout,STATUS='unknown',
     &        FORM='unformatted')
      endif

C     Loop for reading and writing

 100  continue
         call readfileini(x,w,xreal4)
         if(ios.lt.0)goto 200
         write(*,'(a,i5)')' Converting snapshot',ipict


C        iw=1 
C	do ix3=1,nx3
C            do ix2=1,nx2
C               do ix1=1,nx1
C                  if(w(ix1,ix2,ix3,iw).lt.-1.0e16) w(ix1,ix2,ix3,iw)=-1.0e16
C                  if(w(ix1,ix2,ix3,iw).gt.1.0e16) w(ix1,ix2,ix3,iw)=1.0e16
C               enddo
C            enddo
C         enddo



         if(typefileout.ne.'special')then
            call savefileout(x,w)
         else
            call savespecial(x,w)
         endif
         goto 100
 200  continue

      return
      end

C ============================================================================

      subroutine readfileini(x,w,xreal4)

C     Read a snapshot from a VAC ASCII or BINARY data file
C     Only return if readpict is true for the picture
C     Set ios to negative if eof is encountered or npict is reached

      include 'convertdef.f'

      REAL*8 x(nx1,nx2,nx3,ndim), w(nx1,nx2,nx3,nw)
      REAL   xreal4(nx1,nx2,nx3,ndim)
      integer ix1,ix2,ix3,idim,iw,ieqpar

C ----------------------------------------------------------------------------

C     Pretend end of file if the last picture has been read
      if(ipict.eq.npict)then
         ios=-1
         return
      endif

 100  continue

      call readheader
      if(ios .lt. 0)return

      if(typefileini.eq.'ascii')then
         do ix3=1,nx3
            do ix2=1,nx2
               do ix1=1,nx1
                  read(unitini,*)
     &                 (x(ix1,ix2,ix3,idim),idim=1,ndim),
     &                 (w(ix1,ix2,ix3,iw),iw=1,nw)
               enddo
            enddo
         enddo
      else if(typefileini.eq.'binary')then
         read(unitini)x
         do iw=1,nw
            read(unitini)
     &          (((w(ix1,ix2,ix3,iw),ix1=1,nx1),ix2=1,nx2),ix3=1,nx3)
         enddo
      else if(typefileini.eq.'real4')then
         read(unitini)xreal4
         do idim=1,ndim
            do ix3=1,nx3
               do ix2=1,nx2
                  do ix1=1,nx1
                     x(ix1,ix2,ix3,idim)=dble(xreal4(ix1,ix2,ix3,idim))
                  enddo
               enddo
            enddo
         end do
         do iw=1,nw
            read(unitini)
     &      (((xreal4(ix1,ix2,ix3,1),ix1=1,nx1),ix2=1,nx2),ix3=1,nx3)
            do ix3=1,nx3
               do ix2=1,nx2
                  do ix1=1,nx1
                     w(ix1,ix2,ix3,iw)=dble(xreal4(ix1,ix2,ix3,1))
                  enddo
               enddo
            enddo
         end do
      else
         stop 'Unknown input file type!'
      endif

      ipict=ipict+1

      if(readpict(ipict))then
C        Return with read data if picture is to be read
         return
      else
C        Continue reading the file if picture is not to be read
         goto 100
      endif

      end

C ============================================================================

      subroutine savefileout(x,w)

C     Save a snapshot to a VAC ASCII or BINARY data file

      include 'convertdef.f'

      REAL*8 x(nx1,nx2,nx3,ndim), w(nx1,nx2,nx3,nw)
      integer ix1,ix2,ix3,idim,iw,ieqpar

C ----------------------------------------------------------------------------




      if(typefileout.eq.'ascii')then
         write(unitout,'(a)')
     &        fileheadini
         write(unitout,'(i7,1pe13.5,3i3)')
     &        it,t,signedndim,neqpar,nw
         write(unitout,'(3i4)')
     &        (nx(idim),idim=1,ndim)
         write(unitout,'(100(1pe13.5))')
     &        (eqpar(ieqpar),ieqpar=1,neqpar)
         write(unitout,'(a)')
     &        varnames
         do ix3=1,nx3
            do ix2=1,nx2
               do ix1=1,nx1
                  write(unitout,'(100(1pe18.10))')
     &                 (x(ix1,ix2,ix3,idim),idim=1,ndim),
     &                 (w(ix1,ix2,ix3,iw),iw=1,nw)
               enddo
            enddo
         enddo
      else if(typefileout.eq.'binary')then
         write(unitout)fileheadini
         write(unitout)it,t,signedndim,neqpar,nw
         write(unitout)(nx(idim),idim=1,ndim)
         write(unitout)(eqpar(ieqpar),ieqpar=1,neqpar)
         write(unitout)varnames
         write(unitout)x
         do iw=1,nw
            write(unitout)
     &           (((w(ix1,ix2,ix3,iw),ix1=1,nx1),ix2=1,nx2),ix3=1,nx3)
         enddo
      else if(typefileout.eq.'real4')then
         write(unitout)fileheadini
         write(unitout)it,real(t),signedndim,neqpar,nw
         write(unitout)(nx(idim),idim=1,ndim)
         write(unitout)(real(eqpar(ieqpar)),ieqpar=1,neqpar)
         write(unitout)varnames
         write(unitout)((((real(x(ix1,ix2,ix3,idim)),
     &      ix1=1,nx1),ix2=1,nx2),ix3=1,nx3),idim=1,ndim)
         do iw=1,nw
            write(unitout)
     &      (((real(w(ix1,ix2,ix3,iw)),ix1=1,nx1),ix2=1,nx2),ix3=1,nx3)
         enddo
      else
         
      endif

      return
      end

C ============================================================================

      subroutine parsepictstr(pictstr,npictmax,readpict,npict)

C     Extract information from the pictstr string. Determine which
C     snapshots should be read, and the number of the final snapshot.
C
C     Syntax for pictstr: [INT[:INT[:INT]]][, INT[:INT[:INT]]]...
C
C     The string consists of SEGMENTS separated by commas. 
C     
C     Each segment consist of 1, 2, or 3 integers separated by colons.
C
C     1 number is the snapshot number in the file               (e.g. 5).
C     2 numbers separated by a colon give a range of snapshots  (e.g. 5:9).
C     3 numbers separated by a colon give a range with a stride (e.g. 5:9:2).
C
C     The order of the segments is arbitrary. References to snapshot numbers
C     that exceed the number of frames in the file are ignored.
C
C     An empty string means that all pictures are read, just like 1:1000.

      character*79 pictstr
      integer npictmax,npict
      logical readpict(npictmax) 

      character*79 segment,str1,str2,str3
      integer      comma,colon,ipict0,ipict1,dpict,ipict

C --------------------------------------------------------------

      if(pictstr(1:5).eq.'     ')then
         npict=npictmax
         do ipict=1,npictmax
            readpict(ipict)=.true.
         enddo
         return
      endif

      do ipict=1,npictmax
         readpict(ipict)=.false.
      enddo
      npict=0

 100  continue
C        Extract first segment left in pictstr

         comma=index(pictstr,',')
         if(comma.le.0)then
            segment=pictstr
         else
            segment=pictstr(1:comma-1)
            pictstr=pictstr(comma+1:79)
         endif

C        Extract numbers in segment
         colon=index(segment,':')
         if(colon.le.0)then
C           Single number
            read(segment,*)ipict0
            ipict1=ipict0
            dpict=1
         else
            str1=segment(1:colon-1)
            read(str1,*)ipict0
            segment=segment(colon+1:79)
            colon=index(segment,':')
            if(colon.le.0)then
C              Range "ipict0:ipict1"
               read(segment,*)ipict1
               dpict=1
            else
C              Range witha stride "ipict0:ipict1:dpict"
               str2=segment(1:colon-1)
               read(str2,*)ipict1
               str3=segment(colon+1:79)
               read(str3,*)dpict
            endif
         endif

C        Check consistency
         if(ipict1.lt.ipict0.or.ipict0.lt.1.or.ipict1.lt.1
     &      .or.dpict.lt.1)then
            write(*,*)'Incorrect segment:',segment(1:40)
            ipict0=max(1,abs(ipict0))
            ipict1=ipict0
            dpict=1
         endif

         if(ipict0.gt.npictmax.or.ipict1.gt.npictmax)
     &        write(*,*)'Reducing numbers to NPICTMAX=',npictmax,
     &        ' in segment ',segment(1:20)

         ipict0=min(ipict0,npictmax)
         ipict1=min(ipict1,npictmax)

CTEST       write(*,*)'Segment: ',segment
CTEST       write(*,*)'Loop: ',ipict0,ipict1,dpict

C        Assign frames to be read
         do ipict=ipict0,ipict1,dpict
            readpict(ipict)=.true.
            npict=max(ipict,npict)
         enddo

C        Goto next segment if needed
      if(comma.ge.1)goto 100

      return
      end
C ============================================================================
