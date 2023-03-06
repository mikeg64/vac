!==============================================================================
!
!    THE FOLLOWING SUBROUTINES IMPLEMENT REAL*4 INPUT AND OUTPUT
!
!------------------------------------------------------------------------------
!    To include these subroutines into a VACUSR module, at the top write
!
!INCLUDE:vacusr.real4io.t
!
!    If no VACUSR module is needed, use the src/vacusr.t.real4 module with
!
!setvac -u=real4
!
!=============================================================================
subroutine readfileini_special(w)

! Reads from unitini,filenameini in real format.
!
! The file may contain more than one snapshots, in which case the last set is 
! read unless snapshotini is set. 
! The compatibility of initial data with internal parameters is checked.

include 'vacdef.f'

double precision:: w(ixG^T,nw)

logical:: fileexist
integer:: ios                           ! 0 if not EOF, -1 if EOF, >0 if error
integer:: ndimini,neqparini,neqparin,nwini,nwin ! values describing input data
integer:: ix^L,iw,ieqpar,snapshot
double precision:: eqparextra
character*^LENNAME :: varnamesini

! These arrays are needed for temporary arrays
real:: tmpreal4
real:: eqparreal4(neqpar+nspecialpar)
real:: xreal4(ixG^T,ndim)
!-----------------------------------------------------------------------------

oktest=index(teststr,'readfileini')>=1

write(unitterm,*)'ReadFileIni_real4'

inquire(file=filenameini,exist=fileexist)
if(.not.fileexist) call die('Stop: file does not exist, filenameini='// &
   filenameini)
open(unitini,file=filenameini,status='old',form='unformatted')

snapshot=0
do
    read(unitini,iostat=ios,end=100)fileheadini
        if(ios<0)exit                ! Cycle until the last recorded state
        if(oktest) write(unitterm,*)'fileheadini=',fileheadini(1:30)
    read(unitini,iostat=ios)it,tmpreal4,ndimini,neqparini,nwini
    t=dble(tmpreal4)
        if(oktest) write(unitterm, &
"('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")&
           it,t,ndimini,neqparini,nwini
        gencoord= ndimini<0
        call checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)
    read(unitini,iostat=ios)nx
        if(oktest) write(unitterm,"('nx =',3i4)")nx
        call setixGixMix(ix^L)
    read(unitini,iostat=ios)(eqparreal4(ieqpar),ieqpar=1,neqparin),&
                            (tmpreal4,ieqpar=neqparin+1,neqparini)
    do ieqpar=1,neqparin
       eqpar(ieqpar)=dble(eqparreal4(ieqpar))
    enddo
    read(unitini,iostat=ios)varnamesini
    if(varnames=='default')varnames=varnamesini

    read(unitini,iostat=ios)xreal4(ix^S,1:ndim)
    x(ix^S,1:ndim)=dble(xreal4(ix^S,1:ndim))
    ! To conform savefileout_bin we use loop for iw
    do iw=1,nwin
       read(unitini,iostat=ios)xreal4(ix^S,1)
       w(ix^S,iw)=dble(xreal4(ix^S,1))
    end do
    if(ios/=0)then
        write(uniterr,*)'Error: iostat=',ios
        call die('Error in ReadFileIni')
    end if
    snapshot=snapshot+1
    if(snapshot==snapshotini)exit
end do

100 continue

close(unitini)

if(oktest) write(*,*)'x,w:',&
  x(ixtest^D,idimtest),w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine savefileout_special(qunit,w,ix^L)

! This version saves into filename(fileout_) real binary data at every 
! save time in full accordance with the ReadFileini_special subroutine, 
! except that the first line is fileheadout and not fileheadini.

include 'vacdef.f'

integer:: qunit,ix^L,iw,ndimout,ieqpar
double precision:: w(ixG^T,nw)
logical:: fileopen
real:: eqparreal4(neqpar+nspecialpar)
real:: xreal4(ixG^T,ndim)
!-----------------------------------------------------------------------------

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   write(*,*)'SaveFileout_real4'
   open(qunit,file=filenameout,status='unknown',form='unformatted')
endif

if(gencoord)then
   ndimout= -ndim
else
   ndimout= ndim
endif

write(qunit)fileheadout
!write(*,*)fileheadout
write(qunit)it,real(t),ndimout,neqpar+nspecialpar,nw
!write(*,*)it,real(t),ndimout,neqpar+nspecialpar,nw
write(qunit) ixmax^D-ixmin^D+1
!write(*,*) ixmax^D-ixmin^D+1
do ieqpar=1,neqpar+nspecialpar
   eqparreal4(ieqpar)=real(eqpar(ieqpar))
enddo
write(qunit)eqparreal4
!write(*,*)eqparreal4
write(qunit)varnames
!write(*,*)varnames
xreal4(ix^S,1:ndim)=real(x(ix^S,1:ndim))
write(qunit)xreal4(ix^S,1:ndim)
!write(*,*)xreal4(ixtest^D,1:ndim)

do iw=1,nw
   xreal4(ix^S,1)=real(w(ix^S,iw))
   write(qunit)xreal4(ix^S,1)
   !write(*,*)iw,xreal4(ixtest^D,1)
end do

call flushunit(qunit)

return 
end

!=============================================================================
subroutine savefilelog_special(qunit,w,ix^L)

! Save user-defined log data into filename(filelog_) in user-defined format.
! Check savefilelog_default on opening the file etc.

include 'vacdef.f'

integer:: qunit,ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

call die('Special savefilelog is not defined')
end
!=============================================================================
