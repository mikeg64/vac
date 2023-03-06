!##############################################################################
! module vactvdlf 
! Subroutines for TVDLF with Hancock predictor. Also needed for TVD-MUSCL.
!=============================================================================
subroutine hancock(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,wnew)

! Predictor for TVDLF and TVD-MUSCL schemes due to Hancock 
! with flux limiting over conservative variables.

include 'vacdef.f'

double precision:: qdt,qtC,qt
integer:: ixI^L,ixO^L,iws(niw_),idim^LIM
double precision, dimension(ixG^T,nw):: wCT,wnew,wLC,wRC
double precision, dimension(ixG^T):: fLC,fRC,vLC,vRC
integer:: ix^L,hxO^L,iiw,iw,idim
logical:: transport
!-----------------------------------------------------------------------------

oktest=index(teststr,'hancock')>=1
if(oktest)write(*,*)'Hancock: wCT,wnew', &
   wCT(ixtest^D,iwtest),wnew(ixtest^D,iwtest)

! Expand limits in each idim direction in which fluxes are added
do idim= idim^LIM
   ix^L=ixO^L^LADDkr(idim,^D);
enddo
if(ixI^L^LTix^L|.or.|.or.) call die( &
   'Error in Hancock: Nonconforming input limits')

do idim= idim^LIM
   ! Calculate w_j+g_j/2 and w_j-g_j/2
   ! First copy all variables, then upwind wLC and wRC for the iws.
   ! wLC is to the left of ixO, wRC is to the right of wCT.
   !SHIFT
   hxO^L=ixO^L-kr(idim,^D);
   !SHIFT BEGIN
   wRC(hxO^S,1:nw)=wCT(ixO^S,1:nw)
   !SHIFT END
   wLC(ixO^S,1:nw)=wCT(ixO^S,1:nw)
   call upwindLR(ixO^L,hxO^L,iws,idim,wCT,wLC,wRC)

   {^IFGEN
   if(gencoord)then
      call rotatew(hxO^L,idim,wRC)
      call rotatew(ixO^L,idim,wLC)
   endif
   }

   if(oktest)write(*,*)'idim,wLC,wRC:', &
      idim,wLC(ixtest^D,iwtest),wRC(ixtest^D,iwtest)

   ! Calculate fLC=f(w_j+g_j/2) and fRC=f(w_j-g_j/2) 
   ! Add fluxes (eq 4.38c with the choice W~=U, P=I, qdt is not halved here)

   ! Calculate vLC and vRC velocities
   call getv(wRC,hxO^L,idim,vRC)
   call getv(wLC,ixO^L,idim,vLC)
   if(oktest)write(*,*)'vLC,vRC',&
      vLC(ixtest^D),vRC(ixtest^D)

   ! Advect w(iws)
   do iiw=1,iws(niw_); iw=iws(iiw)
      ! Calculate the fLC and fRC fluxes
      call getflux(wRC,hxO^L,iw,idim,fRC,transport)
      call getflux(wLC,ixO^L,iw,idim,fLC,transport)
      if(transport)then
         fRC(hxO^S)=fRC(hxO^S)+vRC(hxO^S)*wRC(hxO^S,iw)
         fLC(ixO^S)=fLC(ixO^S)+vLC(ixO^S)*wLC(ixO^S,iw)
      endif
      ! Advect w(iw)
      
      if(gencoord.and.vectoriw(iw)>=0)then
      {^IFGEN call addflux_rotate(qdt,ixO^L,iw,idim,fLC,ixO^L,fRC,hxO^L,wnew)}
      {^NOGEN call die('Error: Gencoord is off')}
      else
          call addflux(qdt,ixO^L,iw,idim,fLC,ixO^L,fRC,hxO^L,wnew)
      endif
   end do
   if(oktest)write(*,*)'wnew:',wnew(ixtest^D,iwtest)
end do ! next idim

if(typeaxial/='slab'.and.idimmin==1) &
   call addgeometry(qdt,ixO^L,iws,wCT,wnew)
if(sourceunsplit)call addsource2(qdt*(idimmax-idimmin+one)/ndim,&
   ixI^L,ixO^L,iws,qtC,wCT,qt,wnew)

return
end

!=============================================================================
subroutine tvdmusclf(addfluxok,method,qdt,ixII^L,ixO^L,iws,idim^LIM,&
                     qtC,wCT,qt,wnew)

! method=='tvdmu'  --> 2nd oreder (may be 3rd order in 1D) TVD-MUSCL scheme. 
! method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.).
! method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
! method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.
!
! addfluxok=.false. is used for dissipative filter
!
! Flux limiting is over primitive or conservative variables. 
! Based on section 4.4.2. 
!
! Equation 4.39b should read j+1/2 !
! Equation BARMIN(10) should read uL_i+1/2=u_n+1/2_i+...

include 'vacdef.f'

logical:: addfluxok
character*^LENTYPE :: method
double precision:: qdt,qtC,qt
integer:: ixII^L,ixO^L,iws(niw_),idim^LIM,ix
double precision, dimension(ixG^T,nw):: wCT,wnew,wLC,wRC
double precision, dimension(ixG^T):: fLC,fRC,vLC,vRC,cmaxC
double precision:: courantmax
integer:: ixC^L,ixI^L,hxO^L,jxC^L,iiw,iw,idim
logical:: transport,new_cmax
!-----------------------------------------------------------------------------

oktest=index(teststr,'tvd')>=1
if(oktest)write(*,*)'TVDMUSCLF: wCT,wnew', &
   wCT(ixtest^D,iwtest),wnew(ixtest^D,iwtest)

if(idimmax>idimmin .and. typelimited=='original' .and. &
   method/='tvdlf1' .and. method/='tvdmu1') call die( &
   'Error in TVDMUSCLF: Unsplit dim. and original is limited')

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixI^L=ixO^L;
do idim= idim^LIM
   ixI^L=ixI^L^LADD2*kr(idim,^D);
enddo
if(ixII^L^LTixI^L|.or.|.or.)  &
   call die('Error in TVDMUSCLF: Nonconforming input limits')

do idim= idim^LIM
   hxO^L=ixO^L-kr(idim,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 (eq.4.38a,b)
   !SHIFT
   jxC^L=ixC^L+kr(idim,^D);
   !SHIFT BEGIN
   wRC(ixC^S,1:nw)=wCT(jxC^S,1:nw)
   !SHIFT END
   wLC(ixC^S,1:nw)=wCT(ixC^S,1:nw)
   if(method/='tvdlf1'.and.method/='tvdmu1')then
      select case(typelimited)
      case('previous')
         call upwindLR(ixC^L,ixC^L,iws,idim,wold,wLC,wRC)
      case('predictor')
         call upwindLR(ixC^L,ixC^L,iws,idim,wCT,wLC,wRC)
      case('original')
         call upwindLR(ixC^L,ixC^L,iws,idim,wnew,wLC,wRC)
      case default
         call die(&
             'Error in TVDMUSCLF: No such base for limiter:'//typelimited)
      end select
   endif

   {^IFGEN
   if(gencoord)then
      call rotatew(ixC^L,idim,wLC)
      call rotatew(ixC^L,idim,wRC)
   endif
   }

   if(oktest)write(*,*)'idim,wLC,wRC:', &
      idim,wLC(ixtest^D,iwtest),wRC(ixtest^D,iwtest)

   if(method=='tvdlf'.or.method=='tvdlf1')then
      ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
      ! the maximum eigenvalue, it is calculated in advance.
      ! The spectral radius of the Jacobian BARMIN(eq.12)=cmaxC.
      ! Bernardo&Lin's local Lax-Fri (J.Comp.Phys 1989, 84, 90), though
      ! they propose max(cmaxL,cmaxR) instead of cmax((uR+uL)/2) that I use.
      ! To save memory we use wLC temporarily to store the mean 0.5*(uR+uL)
      wLC(ixC^S,1:nw)=half*(wLC(ixC^S,1:nw)+wRC(ixC^S,1:nw))
      new_cmax=.true.
      call getcmax(new_cmax,wLC,ixC^L,idim,cmaxC)
      if(courantpar>zero.and.istep==nstep.and.implpar<zero)then
         ! Calculate dtcourant(idim) for next time step
         if(gencoord)then
            courantmax=maxval(cmaxC(ixC^S)*&
               surfaceC(ixC^S,idim)*two/(dvolume(ixC^S)+dvolume(jxC^S)))
         else
            courantmax=maxval(cmaxC(ixC^S)*two/(dx(ixC^S,idim)+dx(jxC^S,idim)))
         endif
         {^IFMPI call mpiallreduce(courantmax,MPI_MAX)}
         if(qdt*courantmax>one)then
            nerror(couranterr_)=nerror(couranterr_)+1
            if(nerror(couranterr_)==1)write(*,*)&
               'Courant condition error (code=',couranterr_,&
               ') at it=',it,' for direction ',idim,' CFL=',qdt*courantmax
         endif
         if(courantmax>smalldouble)&
            dtcourant(idim)=min(dtcourant(idim),courantpar/courantmax)
      endif
      if(oktest)write(*,*)',cmaxC,wLR:',&
         cmaxC(ixtest^D),wLC(ixtest^D,iwtest)
      ! We regain wLC for further use
      wLC(ixC^S,1:nw)=2*wLC(ixC^S,1:nw)-wRC(ixC^S,1:nw)
   endif ! if TVDLF

   if(addfluxok)then
      ! Calculate velocities for transport fluxes
      call getv(wLC,ixC^L,idim,vLC)
      call getv(wRC,ixC^L,idim,vRC)
   endif

   if(addfluxok.or.method=='tvdlf'.or.method=='tvdlf1')then
      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iws
      do iiw=1,iws(niw_); iw=iws(iiw)
         if(addfluxok)then
            call getflux(wLC,ixC^L,iw,idim,fLC,transport)
            call getflux(wRC,ixC^L,iw,idim,fRC,transport)
            if(transport)then
               fLC(ixC^S)=fLC(ixC^S)+vLC(ixC^S)*wLC(ixC^S,iw)
               fRC(ixC^S)=fRC(ixC^S)+vRC(ixC^S)*wRC(ixC^S,iw)
            endif
            if(oktest.and.iw==iwtest)write(*,*)'fLC,fRC:', &
               fLC(ixtest^D),fRC(ixtest^D)
            ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
            fLC(ixC^S)=half*(fLC(ixC^S)+fRC(ixC^S))
         endif

         ! Add TVDLF dissipation to the flux
         if(method=='tvdlf'.or.method=='tvdlf1')then
            ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
            fRC(ixC^S)=-cmaxC(ixC^S)*half*(wRC(ixC^S,iw)-wLC(ixC^S,iw))

            ! Reduce dissipation by acmcoef if set
            if(acmcoef(iw)>=zero)fRC(ixC^S)=acmcoef(iw)*fRC(ixC^S)

            ! Reduce dissipation by ACM switch if exponent is set
            if(acmexpo>zero)then
                ! We use tmp for the jump.
                if(acmnolim)then
                   !SHIFT BEGIN
                   tmp(ixC^S)=abs(wCT(jxC^S,iw)-wCT(ixC^S,iw))
                   !SHIFT END
                else
                   tmp(ixC^S)=abs(wRC(ixC^S,iw)-wLC(ixC^S,iw))
                endif
                ! acmswitch will overwrite tmp but that is OK.
                call acmswitch(tmp,ixC^L,idim,fRC)
            endif

            if(addfluxok)then
               ! fLC contains physical+dissipative fluxes
               fLC(ixC^S)=fLC(ixC^S)+fRC(ixC^S)
            else
               ! fLC contains dissipative flux only
               fLC(ixC^S)=fRC(ixC^S)
            endif
         endif

         ! Add high order flux_idim to wnew=U_n+1 (eq.4.9)
         if(gencoord.and.vectoriw(iw)>=0)then
            {^IFGEN 
            call addflux_rotate(qdt,ixO^L,iw,idim,fLC,ixO^L,fLC,hxO^L,wnew)}
            {^NOGEN call die('Error: gencoord is off')}
         else
            call addflux(qdt,ixO^L,iw,idim,fLC,ixO^L,fLC,hxO^L,wnew)
         endif

         if(oktest.and.iw==iwtest)write(*,*)'fLR(ix),fLR(hx):', &
            fLC(ixtest^D),fLC(ixtest^D-kr(idim,^D))
      end do ! Next iw

      if(oktest)write(*,*)'wnew:',wnew(ixtest^D,iwtest)

   endif

   ! For the MUSCL scheme apply the characteristic based limiter
   if(method=='tvdmu'.or.method=='tvdmu1')then
      {^IFTVD
      call tvdlimit2(method,qdt,ixC^L,ixO^L,iws,idim,wLC,wRC,wnew)
      if(.false.)}call die('TVD-MUSCL requires the TVD module to be on!')
   endif
enddo ! Next idim

if(oktest)write(*,*)'wnew with fluxes:',wnew(ixtest^D,iwtest)

! Add geometrical and physical sources if physical fluxes were added
if(addfluxok)then
   if(typeaxial/='slab'.and.idimmin==1) &
      call addgeometry(qdt,ixO^L,iws,wCT,wnew)
   if(sourceunsplit)call addsource2(qdt*(idimmax-idimmin+one)/ndim,&
      ixI^L,ixO^L,iws,qtC,wCT,qt,wnew)

   if(oktest)write(*,*)'wnew with source:',wnew(ixtest^D,iwtest)
endif

return
end

!=============================================================================
subroutine upwindLR(ixL^L,ixR^L,iws,idim,w,wLC,wRC)

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. For musclomega=1 
! any typelimiter can be used, for musclomega>1 the 'muscl1','muscl2' limiters,
! minmod(x,omega*y) and minmod(omega*x,y) respectively, are used.

include 'vacdef.f'

double precision:: w(ixG^T,nw),wLC(ixG^T,nw),wRC(ixG^T,nw)
double precision:: ldw(ixG^T),dwC(ixG^T)
integer:: ixL^L,ixR^L,jxR^L,ixC^L,jxC^L,iws(niw_),iiw,iw,idim
!-----------------------------------------------------------------------------

oktest=index(teststr,'upwindlr')>=1
if(oktest)write(*,*)'UpwindLR omega:',musclomega

! Eqs.3.38a and 3.38b

if(oktest.and.idim==idimtest)write(*,*)'wLC,wRC:',&
   wLC(ixtest^D,iwtest),wRC(ixtest^D,iwtest)
if(oktest.and.idim==idimtest)write(*,*)'wg,wh,wi,wj,wk:',&
  w(ixtest^D-2*kr(idim,^D),iwtest),w(ixtest^D-kr(idim,^D),iwtest),&
  w(ixtest^D,iwtest),&
  w(ixtest^D+kr(idim,^D),iwtest),w(ixtest^D+2*kr(idim,^D),iwtest)

! Transform w,wL,wR to primitive variables 
if(useprimitive)then
   call primitive(ixG^L,w)
   call primitive(ixL^L,wLC)
   call primitive(ixR^L,wRC)
endif

!SHIFT
jxR^L=ixR^L+kr(idim,^D);
ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idim,^D); 
!SHIFT MORE
jxC^L=ixC^L+kr(idim,^D);
!SHIFT BEGIN
do iiw=1,iws(niw_); iw=iws(iiw)
   dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
   if(musclomega>one)then
      typelimiter(iw)='muscl1'
      call dwlimiter2(dwC,ixC^L,iw,idim,ldw)
      wRC(ixR^S,iw)=wRC(ixR^S,iw)-muscleta1*ldw(jxR^S)
      wLC(ixL^S,iw)=wLC(ixL^S,iw)+muscleta2*ldw(ixL^S)
      typelimiter(iw)='muscl2'
      call dwlimiter2(dwC,ixC^L,iw,idim,ldw)
      wRC(ixR^S,iw)=wRC(ixR^S,iw)-muscleta2*ldw(jxR^S)
      wLC(ixL^S,iw)=wLC(ixL^S,iw)+muscleta1*ldw(ixL^S)
   else
      call dwlimiter2(dwC,ixC^L,iw,idim,ldw)
      wRC(ixR^S,iw)=wRC(ixR^S,iw)-half*ldw(jxR^S)
      wLC(ixL^S,iw)=wLC(ixL^S,iw)+half*ldw(ixL^S)
   endif
end do
!SHIFT END

! Transform w,wL,wR back to conservative variables 
if(useprimitive)then
   call conserve(ixG^L,w)
   call conserve(ixL^L,wLC)
   call conserve(ixR^L,wRC)
endif

return
end

!=============================================================================
! end module vactvdlf
!##############################################################################
