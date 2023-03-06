!#############################################################################
! module vaccd
! Centered difference scheme
!=============================================================================
subroutine centdiff(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,w)

! Advance the iws flow variables from t to t+qdt within ixO^L by centered 
! differencing in space the dw/dt+dF_i(w)/dx_i=S type equation. 
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'vacdef.f'

double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
integer:: ixI^L,ixO^L,iws(niw_),idim^LIM
logical :: transport,tvdpreproc

double precision:: v(ixG^T,ndim),f(ixG^T),fC(ixG^T)
integer:: iiw,iw,ix^L,hxO^L,ixC^L,jxC^L,idim,idir

logical:: firstsweep,lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

oktest=index(teststr,'centdiff')>=1
if(oktest)write(*,*)'CentDiff wCT,w:',&
   wCT(ixtest^D,iwtest),w(ixtest^D,iwtest)

! An extra layer is needed in each direction for which fluxes are added.
ix^L=ixO^L;
do idim= idim^LIM
   ix^L=ix^L^LADDkr(idim,^D);
enddo
if(ixI^L^LTix^L|.or.|.or.) call die( &
   'Error in CentDiff: Non-conforming input limits')

tvdpreproc=(typefull1=='tvd').and.(typeadvance=='onestep').and.&
           (typeaxial=='cylinder'.or.(typeaxial=='slab'.and.sourceunsplit))

if(index(teststr,'noryuprep')>=1)tvdpreproc=.false.

if(oktest.and.tvdpreproc)write(*,*)'TVD preproc is done:',firstsweep

if(tvdpreproc.and.firstsweep)then
   ! Add half of curveture and source terms following Ryu, ApJ 452,364, eq.(A9)
   if(typeaxial=='cylinder')then
      ! Apply curvature correction
      ! wCT=w-(dt/2)*F_r(w)/r
      ! In case of angmomfix
      ! wCT=w-dt*F_r(w)/r
      ! We use fC for velocity in the r_ direction
      call getv(w,ixM^L,r_,fC)
      do iiw=1,iws(niw_); iw=iws(iiw)
         call getflux(w,ixM^L,iw,r_,f,transport)
         if(transport)f(ixM^S)=f(ixM^S)+fC(ixM^S)*w(ixM^S,iw)
         if(angmomfix.and.iw==mphi_)then
            wCT(ixM^S,iw)=wCT(ixM^S,iw)-qdt*f(ixM^S)/x(ixM^S,r_)
         else
            wCT(ixM^S,iw)=wCT(ixM^S,iw)-(qdt/2)*f(ixM^S)/x(ixM^S,r_)
         endif
      enddo
      ! Add half of geometrical sources as well
      call addgeometry(qdt/2,ixM^L,iws,w,wCT)
   endif
   ! Add half of unsplit sources
   if(sourceunsplit)call addsource2(qdt/2,ixG^L,ixM^L,iws,qt,w,qt,wCT)

   call getboundary(qt+qdt/2,1,nw,1,ndim,wCT)
!!! This is just a partial step, so processing is probably not useful here
!   if(nproc(2)/=0.and.implpar<=zero)then
!      if(nproc(2)<0.or.(it-itmin==((it-itmin)/nproc(2))*nproc(2)))then
!         call process(2,idim^LIM,wCT)
!      endif
!!! endif
endif

! Add fluxes to w
do idim= idim^LIM
   ix^L=ixO^L^LADDkr(idim,^D); ixCmin^D=ixmin^D; ixCmax^D=ixOmax^D;
   !SHIFT
   jxC^L=ixC^L+kr(idim,^D); 
   hxO^L=ixO^L-kr(idim,^D);

   if(gencoord)then
      do idir=1,ndim
         call getv(wCT,ix^L,idir,f)
         v(ix^S,idir)=f(ix^S)
      enddo
   else
      call getv(wCT,ix^L,idim,f)
      v(ix^S,idim)=f(ix^S)
   endif
   do iiw=1,iws(niw_); iw=iws(iiw)
      if(gencoord)then
         fC(ixC^S)=zero
         do idir=1,ndim
            call getflux(wCT,ix^L,iw,idir,f,transport)
            ! Add transport flux
            if(transport)f(ix^S)=f(ix^S)+v(ix^S,idir)*wCT(ix^S,iw)
            ! Center flux to interface 
            !SHIFT BEGIN
            tmp(ixC^S)=half*(f(ixC^S)+f(jxC^S))
            !SHIFT END

            ! Store flux for fluxCT or fluxCD schemes
            ! Check idir==idim since we do not project on the normal vector
            if((typeconstrain.eq.'fluxCT'.or.typeconstrain.eq.'fluxCD')&
               .and.idim==idir.and.iw/=b0_+idim &
               .and.iw>b0_.and.iw<=b0_+ndim.and.istep==nstep) &
                call storeflux(qdt,tmp,ixC^L,idim,iw)

            ! Project to normal vector
            fC(ixC^S)=fC(ixC^S)+normalC(ixC^S,idim,idir)*tmp(ixC^S)
         enddo
      else
         ! Get non-transported flux
         call getflux(wCT,ix^L,iw,idim,f,transport)
         if(oktest.and.iw==iwtest)write(*,*)'  fj,fi,fh:',&
            f(ixtest^D+kr(idim,^D)),f(ixtest^D),f(ixtest^D-kr(idim,^D))
         ! Add transport flux
         if(transport)f(ix^S)=f(ix^S)+v(ix^S,idim)*wCT(ix^S,iw)
         if(oktest.and.iw==iwtest)write(*,*)'v+fj,fi,fh:',&
            f(ixtest^D+kr(idim,^D)),f(ixtest^D),f(ixtest^D-kr(idim,^D))
         ! Center flux to interface
         !SHIFT BEGIN
         fC(ixC^S)=half*(f(ixC^S)+f(jxC^S))
         !SHIFT END
      endif

      call addflux(qdt,ixO^L,iw,idim,fC,ixO^L,fC,hxO^L,w)

      if(oktest.and.iw==iwtest)write(*,*)'fCi,fCh:',&
            fC(ixtest^D),fC(ixtest^D-kr(idim,^D))
      if(oktest.and.iw==iwtest)write(*,*)'idim,wnew:',idim,w(ixtest^D,iw)
   end do    !next iw
end do       !next idim

if(.not.tvdpreproc)then
   if(typeaxial/='slab'.and.idimmin==r_) &
      call addgeometry(qdt,ixO^L,iws,wCT,w)
   if(sourceunsplit) &
      call addsource2(qdt*(idimmax-idimmin+one)/ndim, &
         ixI^L,ixO^L,iws,qtC,wCT,qt,w)
endif

return
end

!=============================================================================
subroutine centdiff4(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,w)

! Advance the iws flow variables from t to t+qdt within ixO^L by 
! fourth order centered  differencing in space the dw/dt+dF_i(w)/dx_i=S 
! type equation. 
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'vacdef.f'

double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
integer:: ixI^L,ixO^L,iws(niw_),idim^LIM
logical :: transport

double precision:: v(ixG^T),f(ixG^T)
integer:: iiw,iw,ix^L,idim,idir
!-----------------------------------------------------------------------------

oktest=index(teststr,'centdiff')>=1
if(oktest)write(*,*)'CentDiff4 wCT,w:',&
   wCT(ixtest^D,iwtest),w(ixtest^D,iwtest)

if(gencoord)call die('CentDiff4 is not implemented for gen.coords yet!')

! Two extra layers are needed in each direction for which fluxes are added.
ix^L=ixO^L;
do idim= idim^LIM
   ix^L=ix^L^LADD2*kr(idim,^D);
enddo
if(ixI^L^LTix^L|.or.|.or.) call die( &
   'Error in CentDiff4: Non-conforming input limits')

! Add fluxes to w
do idim= idim^LIM
   ix^L=ixO^L^LADD2*kr(idim,^D);

   call getv(wCT,ix^L,idim,v)

   do iiw=1,iws(niw_); iw=iws(iiw)

      ! Get non-transported flux
      call getflux(wCT,ix^L,iw,idim,f,transport)

      if(oktest.and.iw==iwtest)write(*,*)'  fj,fi,fh:',&
         f(ixtest^D+kr(idim,^D)),f(ixtest^D),f(ixtest^D-kr(idim,^D))

      ! Add transport flux
      if(transport)f(ix^S)=f(ix^S)+v(ix^S)*wCT(ix^S,iw)
      if(oktest.and.iw==iwtest)write(*,*)'v+fj,fi,fh:',&
         f(ixtest^D+kr(idim,^D)),f(ixtest^D),f(ixtest^D-kr(idim,^D))

      ! Add divergence of flux
      call gradient4(.false.,f,ixO^L,idim,tmp)
      w(ixO^S,iw)=w(ixO^S,iw)-qdt*tmp(ixO^S)

      if(oktest.and.iw==iwtest)write(*,*)'fk,fj,fh,fg:',&
         f(ixtest^D+2*kr(idim,^D)),f(ixtest^D+kr(idim,^D)),&
         f(ixtest^D-kr(idim,^D)),f(ixtest^D-2*kr(idim,^D))
      if(oktest.and.iw==iwtest)write(*,*)'idim,wnew:',idim,w(ixtest^D,iw)

   end do    !next iw
end do       !next idim

! Add sources
if(typeaxial/='slab'.and.idimmin==r_) &
   call addgeometry(qdt,ixO^L,iws,wCT,w)
if(sourceunsplit) &
   call addsource2(qdt*(idimmax-idimmin+one)/ndim, &
      ixI^L,ixO^L,iws,qtC,wCT,qt,w)

return
end

!=============================================================================
! end module vaccd
!#############################################################################
