!=============================================================================
subroutine savefile(filehead,it,t,ndim,ndimgen,neqpar,nw,nx,eqpar,varnames,&
                    xsize,xbuf,wbuf)

! Saves a binary VAC file to STDOUT

integer:: it,ndim,ndimgen,neqpar,nw,nx(ndim),xsize
double precision:: t,eqpar(neqpar),xbuf(ndim*xsize),wbuf(xsize,nw)
character*79 :: filehead,varnames

integer:: iw
logical:: oktest
!-----------------------------------------------------------------------------

oktest=.false.

ndim=abs(ndimgen)

write(*)filehead
   if(oktest)write(98,*)filehead
write(*)it,t,ndimgen,neqpar,nw
   if(oktest)write(98,*)it,t,ndimgen,neqpar,nw
write(*)nx
   if(oktest)write(98,*)nx
write(*)eqpar
   if(oktest)write(98,*),eqpar
write(*)varnames
   if(oktest)write(98,*)varnames

write(*)xbuf
    if(oktest)write(98,*)'xbuf written'

do iw=1,nw
   write(*)wbuf(:,iw)
   if(oktest)write(98,*)'wbuf',iw
end do

return 
end
!=============================================================================
