!=============================================================================
program convertbinary

! Converts between different binary formats for VAC data files.
! The file may contain more snapshots, they are all read and converted. 

integer:: it,ndimgen,neqpar,nw,xsize
integer, allocatable:: nx(:)
double precision:: t
double precision, allocatable:: eqpar(:),xbuf(:),wbuf(:,:)
character*79 :: filehead,varnames

logical:: oktest
integer:: ios                           ! 0 if not EOF, -1 if EOF, >0 if error
integer:: ipict,idim,ndim,iw,ieqpar
!-----------------------------------------------------------------------------

oktest=.true.

ipict=0
do
    read(*,iostat=ios)filehead
        if(ios<0)exit                ! Cycle until the last recorded state
        if(ios>0) stop 'Error reading filehead'
        if(oktest) write(99,"('filehead=',a30)")filehead

    ipict=ipict+1

    read(*,iostat=ios)it,t,ndimgen,neqpar,nw
        if(ios>0) stop 'Error reading it,t,ndimgen,neqpar,nw'
        if(oktest) write(99, &
           "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")&
           it,t,ndimgen,neqpar,nw
        ndim=abs(ndimgen)

    if(ipict==1) allocate(nx(ndim),eqpar(neqpar))

    read(*,iostat=ios)nx
        if(ios>0) stop 'Error reading nx'
        if(oktest) write(99,"('nx =',3i4)")nx
    read(*,iostat=ios)eqpar
        if(ios>0) stop 'Error reading eqpar'
      if(oktest) write(99,"('eqpar =',10g14.7)")eqpar
    read(*,iostat=ios)varnames
        if(ios>0) stop 'Error reading varnames'
        if(oktest) write(99,"('varnames=',a30)")varnames

    xsize=product(nx)
    if(ipict==1)allocate(xbuf(xsize*ndim),wbuf(xsize,nw))

    read(*,iostat=ios)xbuf
        if(ios>0) stop 'Error reading xbuf'
        if(oktest) write(99,*)'xbuf read'

    ! To conform savefileout_bin we use loop for iw
    do iw=1,nw
       read(*,iostat=ios)wbuf(:,iw)
          if(ios>0) stop 'Error reading wbuf'
          if(oktest) write(99,*)'wbuf',iw
    end do

    call savefile(filehead,it,t,ndim,ndimgen,neqpar,nw,nx,eqpar,varnames,&
                  xsize,xbuf,wbuf)
end do

stop
end
!=============================================================================
