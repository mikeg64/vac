!##############################################################################
! module vactvd
! Subroutines for TVD, TVD-MacCormack and TVD-MUSCL schemes
!=============================================================================
subroutine tvdlimit(method,qdt,ixI^L,ixO^L,iws,idim^LIM,w,qt,wnew)

! Limit the iws flow variables in wnew. 
! Limiting is based on w or wnew depending of "typelimited".
! Create the left and right values by shifting then call tvdlimit2.
! "method" can be temporally second order 'tvd' or first order 'tvd1'
! For the 2nd order 'tvd' source terms are also added to second order accuracy

include 'vacdef.f'

character*^LENTYPE :: method
double precision:: qdt,qt
integer:: ixI^L,ixO^L,iws(niw_),idim^LIM
double precision, dimension(ixG^T,nw):: w,wnew,wR
!!! wL
integer:: ixIC^L,jxIC^L,iiw,iw,idim

logical:: firstsweep,lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

oktest=index(teststr,'tvd')>=1
if(oktest)write(*,*)'TVDLimit w,wnew,method:',&
   w(ixtest^D,iwtest),wnew(ixtest^D,iwtest),method

if(typelimited=='predictor')then
   !!!Side effect, w is overwritten, but after TVDlimit w is not needed anymore
   call getboundary(qt,1,nw,idim^LIM,wnew)
   do iiw=1,iws(niw_); iw=iws(iiw)
      w(ixI^S,iw)=wnew(ixI^S,iw)
   end do
end if

!!!Side effect, w is overwritten
!!!if(typelimited=='previous')w(ixI^S,1:nw)=wold(ixI^S,1:nw)

do idim= idim^LIM
   ixICmax^D=ixOmax^D+kr(idim,^D); ixICmin^D=ixOmin^D-2*kr(idim,^D);

   !!! wL(ixIC^S,1:nw)=w(ixIC^S,1:nw)
   !SHIFT
   jxIC^L=ixIC^L+kr(idim,^D);
   !SHIFT BEGIN
   wR(ixIC^S,1:nw)=w(jxIC^S,1:nw)
   !SHIFT END
   {^IFGEN
   if(gencoord)then
       if(idim<idimmax)call die('TVD does not work for unsplit gen.coord.!')
       call rotatew(ixIC^L,idim,w)  !!! w-->wL
       call rotatew(ixIC^L,idim,wR)
   endif
   }
   call tvdlimit2(method,qdt,ixIC^L,ixO^L,iws,idim,w,wR,wnew) !!!w-->wL
enddo

if(method=='tvd'.and.lastsweep.and.index(teststr,'noryuprep')<1.and.&
   (typeaxial=='cylinder'.or.(typeaxial=='slab'.and.sourceunsplit)))then

   if(oktest)write(*,*)'Adding Ryu style (geometrical) sources'

   ! Following Ryu et al (APJ 452,364) (geometrical) sources are added 
   ! with second order accuracy:
   ! whalf=(wold+whydro+S)/2     eq.A22 corrected by including S
   ! wnew=whydro+S(whalf)        eq.A21

   !!!Side effect: We overwrite w, but that is OK, it is not needed anymore
   w(ixO^S,1:nw)=wnew(ixO^S,1:nw)

   if(typeaxial/='slab')call addgeometry(qdt,ixO^L,iws,wold,w)
   if(sourceunsplit)call addsource2(qdt,ixG^L,ixO^L,iws,t+qdt,wold,t+qdt,w)

   w(ixO^S,1:nw)=half*(wold(ixO^S,1:nw)+w(ixO^S,1:nw))

   if(typeaxial/='slab')call addgeometry(qdt,ixO^L,iws,w,wnew)
   if(sourceunsplit)call addsource2(qdt,ixO^L,ixO^L,iws,t+qdt,w,t+qdt,wnew)
endif

return
end

!=============================================================================
subroutine tvdlimit2(method,qdt,ixIC^L,ixO^L,iws,idim,wL,wR,wnew)

! Limit the iws flow variables in wnew according to typetvd. 
! wroeC is based on wL and wR.
! If method=='tvd' an extra adtdx**2*jumpC is added to phiC for 2nd order
! accuracy in time.

include 'vacdef.f'

character*^LENTYPE :: method
double precision:: qdt
integer:: ixIC^L,ixO^L,iws(niw_),idim
double precision, dimension(ixG^T,nw):: wL,wR,wnew,wroeC
double precision, dimension(ixG^T)   :: phiC,rphiC,jumpC,adtdxC,smallaC
double precision:: courantmax,wtest(nw)
integer:: hxO^L,ixC^L,jxC^L,jxIC^L,iiw,iw,il
!-----------------------------------------------------------------------------

oktest=index(teststr,'tvdlimit')>=1
if(oktest)write(*,*)'TVDLimit2 wL,wR,wnew:',&
   wL(ixtest^D,iwtest),wR(ixtest^D,iwtest),wnew(ixtest^D,iwtest)
if(typetvd=='testing')wtest(1:nw)=0.

hxO^L=ixO^L-kr(idim,^D);
ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D; 

jxC^L=ixC^L+kr(idim,^D);
jxIC^L=ixIC^L+kr(idim,^D);

call average(wL,wR,ixIC^L,iws,idim,wroeC)

! A loop on characteristic variables to calculate the dissipative flux phiC.
do il=1,nw
   !Calculate the jump in the il-th characteristic variable: L(wroe)*dw
   call geteigenjump(wL,wR,wroeC,ixIC^L,il,idim,smallaC,adtdxC,jumpC)

   if(oktest)write(*,*)'wL(ix,il),wR(ix,il):',&
      wL(ixtest^D,il),wR(ixtest^D,il)
   if(oktest)write(*,*)'L*dwC,a:',&
      jumpC(ixtest^D),adtdxC(ixtest^D)

   ! Normalize the eigenvalue "a" (and its limit "smalla" if needed):
   if(gencoord)then
      adtdxC(ixIC^S)=adtdxC(ixIC^S)*qdt*surfaceC(ixIC^S,idim)*&
         two/(dvolume(ixIC^S)+dvolume(jxIC^S))
      if(typeentropy(il)=='harten' .or. typeentropy(il)=='powell')&
      smallaC(ixIC^S)=smallaC(ixIC^S)*qdt*surfaceC(ixIC^S,idim)*&
         two/(dvolume(ixIC^S)+dvolume(jxIC^S))
   else
      adtdxC(ixIC^S)=adtdxC(ixIC^S)*qdt*&
         two/(dx(ixIC^S,idim)+dx(jxIC^S,idim))
      if(typeentropy(il)=='harten' .or. typeentropy(il)=='powell')&
      smallaC(ixIC^S)=smallaC(ixIC^S)*qdt*&
         two/(dx(ixIC^S,idim)+dx(jxIC^S,idim))
   endif

   if(oktest)write(*,*)'qdt,adtdxC:',qdt,adtdxC(ixtest^D)

   ! For potentially extreme eigenvalues calculate dtcourant for next step
   if((il==extremeLW_.or.il==extremeRW_).and.&
      courantpar>zero.and.istep==nstep.and.implpar<zero)then
      courantmax=maxval(abs(adtdxC(ixIC^S)))
      {^IFMPI call mpiallreduce(courantmax,MPI_MAX)}
      if(courantmax>one) then
         nerror(couranterr_)=nerror(couranterr_)+1
         if(nerror(couranterr_)==1)write(*,*)&
            'Courant condition error (code=',couranterr_,&
            ') at it=',it,' for direction',idim,' CFL=',courantmax
      endif
      if(courantmax>smalldouble)&
         dtcourant(idim)=min(dtcourant(idim),qdt*courantpar/courantmax)
   endif

   ! Calculate the flux limiter function phi
   call getphi(method,jumpC,adtdxC,smallaC,ixIC^L,ixC^L,il,idim,phiC)

   if(oktest)write(*,*)'il,phiC:',il,phiC(ixtest^D)

   !Multiply phiC by 2>acmcoef>0 (kappa in eq.2.21 Yee, Sandham, Djomehri)
   if(acmcoef(il)>=zero)phiC(ixC^S)=acmcoef(il)*phiC(ixC^S)

   !Multiply by a term of Artificial Compression Method (ACM)
   !(theta_j+1/2 in eq.2.21 Yee, Sandham, Djomehri) if axmexpo>0 is set

   if(acmexpo>zero)call acmswitch(jumpC,ixC^L,idim,phiC)

   if(oktest)write(*,*)'phiCstar:',phiC(ixtest^D)

   !Add R(iw,il)*phiC(il) to each variable iw in wnew
   do iiw=1,iws(niw_); iw=iws(iiw)
      call rtimes(phiC,wroeC,ixC^L,iw,il,idim,rphiC)

      !Based on Eqs 5.6a,b in Yee
      rphiC(ixC^S)=rphiC(ixC^S)*quarter*(dvolume(jxC^S)+dvolume(ixC^S))
      if(gencoord.and.vectoriw(iw)>=0)then
         ! The dt=-one tells addflux_rotate that this is an artificial flux
         {^IFGEN 
         call addflux_rotate(-one,ixO^L,iw,idim,rphiC,ixO^L,rphiC,hxO^L,wnew)}
         {^NOGEN call die('Error: gencoord is off')}
      else
         ! The dt=-one tells addflux that this is an artificial flux
         call addflux(-one,ixO^L,iw,idim,rphiC,ixO^L,rphiC,hxO^L,wnew)
      endif

      if(oktest)write(*,*)'rphiCi,rphiCh:',&
         rphiC(ixtest^D)/dvolume(ixtest^D),&
         rphiC(ixtest^D-kr(idim,^D))/dvolume(ixtest^D)
      if(typetvd=='testing')wtest(iw)=wtest(iw)+rphiC(ixtest^D)
   end do  !iw
   if(oktest)write(*,*)'wnew:',wnew(ixtest^D,iwtest)
end do     !il

if(typetvd=='testing')then
   write(*,*)'wtest                   wR-wL'
   do iiw=1,iws(niw_); iw=iws(iiw)
      write(*,*)wtest(iw),wR(ixtest^D,iw)-wL(ixtest^D,iw)
   end do
endif

return
end

!=============================================================================
subroutine getphi(method,jumpC,adtdxC,smallaC,ixIC^L,ixC^L,il,idim,phiC)

! Calculate the dissipative flux from jumpC=L*dw and adtdx=eigenvalue*dt/dx.
! Add Lax-Wendroff type correction if method=='tvd'.
! Limit according to method and typetvd.

include 'vacdef.f'

character*^LENTYPE :: method
integer:: ixIC^L,ixC^L
double precision, dimension(ixG^T):: jumpC,ljumpC,adtdxC,smallaC,phiC
integer:: jxC^L,ix^L,hx^L,il,idim
!-----------------------------------------------------------------------------

oktest=index(teststr,'getphi')>=1
if(oktest)write(*,*)'Getphi jumpC,adtdxC:',&
   jumpC(ixtest^D),adtdxC(ixtest^D)

if(method=='tvdmu'.or.method=='tvdmu1')then
   ! In the MUSCL scheme phi=|a|*jump, apply entropy fix to it
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixC^S)=abs(adtdxC(ixC^S))*jumpC(ixC^S)
   else
      where(abs(adtdxC(ixC^S))>=smallaC(ixC^S))
         phiC(ixC^S)=abs(adtdxC(ixC^S))*jumpC(ixC^S)
      elsewhere
         phiC(ixC^S)=half*(smallaC(ixC^S)+adtdxC(ixC^S)**2/smallaC(ixC^S))&
            *jumpC(ixC^S)
      endwhere
   endif
   ! That's all for the MUSCL scheme
   if(oktest)write(*,*)'GetPhi Final phiC:',phiC(ixtest^D)
   return
endif

if(method=='tvd')then
   !Entropy fix to |a|-a**2
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
     phiC(ixIC^S)=abs(adtdxC(ixIC^S))-adtdxC(ixIC^S)**2
   else
      where(abs(adtdxC(ixIC^S))>=smallaC(ixIC^S))
         phiC(ixIC^S)=abs(adtdxC(ixIC^S))-adtdxC(ixIC^S)**2
      elsewhere
         phiC(ixIC^S)=half*smallaC(ixIC^S)+&
            (half/smallaC(ixIC^S)-one)*adtdxC(ixIC^S)**2
      endwhere
   endif
   if(oktest)write(*,*)'abs(nu)-nu**2:',phiC(ixtest^D)
else
   !Entropy fix to |a|
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixIC^S)=abs(adtdxC(ixIC^S))
   else
      where(abs(adtdxC(ixIC^S))>=smallaC(ixIC^S))
         phiC(ixIC^S)=abs(adtdxC(ixIC^S))
      elsewhere
         phiC(ixIC^S)=half*smallaC(ixIC^S)+&
            half/smallaC(ixIC^S)*adtdxC(ixIC^S)**2
      endwhere
   endif
   if(oktest)write(*,*)'abs(nu)-nu**2:',phiC(ixtest^D)
endif

!SHIFT
jxC^L=ixC^L+kr(idim,^D);

hxmin^D=ixICmin^D; hxmax^D=ixICmax^D-kr(idim,^D);
!SHIFT MORE
ix^L=hx^L+kr(idim,^D);

!SHIFT BEGIN
select case(typetvd)
   case('symmetric')
      !eq.3.53 and eq.3.69d
      call dwlimiter3(jumpC,ixIC^L,il,idim,ljumpC)
      phiC(ixC^S)=phiC(ixC^S)*(jumpC(ixC^S)-ljumpC(ixC^S))
      !extra (a*lambda)**2*delta
      if(method=='tvd')phiC(ixC^S)=phiC(ixC^S)+adtdxC(ixC^S)**2*jumpC(ixC^S)
   case('roe')
      !eq.3.52 with correction of sign of a
      select case(typelimiter(il))
      case('roe','superroe')
         call dwlimiterroe(adtdxC,jumpC,ixC^L,il,idim,ljumpC)
         phiC(ixC^S)=phiC(ixC^S)*(jumpC(ixC^S)-ljumpC(ixC^S))
      case default
         call dwlimiter2(jumpC,ixIC^L,il,idim,ljumpC)
         where(adtdxC(ixC^S)<=0)
            phiC(ixC^S)=phiC(ixC^S)*(jumpC(ixC^S)-ljumpC(jxC^S))
         elsewhere
            phiC(ixC^S)=phiC(ixC^S)*(jumpC(ixC^S)-ljumpC(ixC^S))
         end where
      end select
      !extra (a*lambda)**2*delta
      if(method=='tvd')phiC(ixC^S)=phiC(ixC^S)+adtdxC(ixC^S)**2*jumpC(ixC^S)
   case('sweby')
      !Sweby eqs.4.11-4.15, but no 0.5 ?!
      phiC(ixIC^S)=phiC(ixIC^S)*jumpC(ixIC^S)
      select case(typelimiter(il))
      case('roe','superroe')
         call dwlimiterroe(adtdxC,phiC,ixC^L,il,idim,ljumpC)
         phiC(ixC^S)=phiC(ixC^S)-ljumpC(ixC^S)
      case default
         call dwlimiter2(phiC,ixIC^L,il,idim,ljumpC)
         where(adtdxC(ixC^S)<=0)
            phiC(ixC^S)=phiC(ixC^S)-ljumpC(jxC^S)
         elsewhere
            phiC(ixC^S)=phiC(ixC^S)-ljumpC(ixC^S)
         end where
      end select
      !extra (a*lambda)**2*delta
      if(method=='tvd')phiC(ixC^S)=phiC(ixC^S)+adtdxC(ixC^S)**2*jumpC(ixC^S)
   case('yee')
      !eq.3.51 with correction
      call dwlimiter2(jumpC,ixIC^L,il,idim,ljumpC)

      if(oktest)write(*,*)'Limiter jumpC,ljumpC(ix,jx):',&
          jumpC(ixtest^D),ljumpC(ixtest^D),ljumpC(ixtest^D+kr(idim,^D))

      !Use phiC as 0.5*(|nu|-nu**2) eq.3.45e for tvd otherwise 0.5*|nu|
      phiC(ixC^S)=half*phiC(ixC^S)
      !gamma*lambda eq.3.51d, use tmp to store agdtdxC
      where(abs(jumpC(ixC^S))>smalldouble)
         tmp(ixC^S)=adtdxC(ixC^S)+phiC(ixC^S)*&
         (ljumpC(jxC^S)-ljumpC(ixC^S))/jumpC(ixC^S)
      elsewhere
         tmp(ixC^S)=adtdxC(ixC^S)
      end where

      if(oktest)write(*,*)'agdtdxC:',tmp(ixtest^D)

      !eq.3.51a
      if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
         phiC(ixC^S)=-phiC(ixC^S)*(ljumpC(jxC^S)+ljumpC(ixC^S))+&
            abs(tmp(ixC^S))*jumpC(ixC^S)
      else
         where(abs(tmp(ixC^S))>=smallaC(ixC^S))
            phiC(ixC^S)=-phiC(ixC^S)*(ljumpC(jxC^S)+ljumpC(ixC^S))+&
               abs(tmp(ixC^S))*jumpC(ixC^S)
         elsewhere
            phiC(ixC^S)=-phiC(ixC^S)*(ljumpC(jxC^S)+ljumpC(ixC^S))+&
               (half*smallaC(ixC^S)+half/smallaC(ixC^S)*&
               tmp(ixC^S)**2)*jumpC(ixC^S)
         endwhere
      endif
   case('harten')
      !See Ryu, section 2.3
      !Use phiC as 0.5*(|nu|-nu**2)*jumpC eq.3.45b,e
      phiC(ixIC^S)=half*phiC(ixIC^S)*jumpC(ixIC^S)
      call dwlimiter2(phiC,ixIC^L,il,idim,ljumpC)

      if(artcomp(il))then
         ! ARTIFICIAL COMPRESSION, ORIGINALLY IMPLEMENTED BY M. NAUTA
         !
         ! Equations 5.8b, 5.7c, 5.8a, 5.7a of Harten change ljumpC,
         ! or Ryu's  2.98, 2.97, 2.96:
         !
         ! sigma=(1-a*|dt/dx|)/2
         ! gbar_i=minmod(sigma*alfa_i-1/2,sigma*alfa_i+1/2)
         ! theta_i=(|alfa_i+1/2 - alfa_i-1/2|)/(|alfa_i+1/2|+|alfa_i-1/2|)
         ! ljumpC=ljumpC+theta*gbar
         !
         ! To save memory tmp is used for sigma & theta, tmp2 is used for gbar

         tmp(ixIC^S)=half*(one-abs(adtdxC(ixIC^S)))*jumpC(ixIC^S)
         tmp2(ix^S)=sign(one,jumpC(ix^S))*&
             max(zero,min(sign(one,jumpC(ix^S))*tmp(hx^S),abs(tmp(ix^S))))
         where(abs(jumpC(ix^S))+abs(jumpC(hx^S))>smalldouble)
            tmp(ix^S)=abs(jumpC(ix^S)-jumpC(hx^S))/&
                  (abs(jumpC(ix^S))+abs(jumpC(hx^S)))
         elsewhere
            tmp(ix^S)=zero
         end where
         ljumpC(ix^S)=ljumpC(ix^S)+tmp(ix^S)*tmp2(ix^S)
      endif

      !gamma*lambda eq.3.45d, use tmp as agdtdxC
      where(abs(jumpC(ixC^S))>smalldouble)
         tmp(ixC^S)=adtdxC(ixC^S)+&
         (ljumpC(jxC^S)-ljumpC(ixC^S))/jumpC(ixC^S)
      elsewhere
         tmp(ixC^S)=adtdxC(ixC^S)
      end where
      if(oktest)write(*,*)'agdtdxC',tmp(ixtest^D)
      !eq.3.45a with correction
      if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
         phiC(ixC^S)=-ljumpC(jxC^S)-ljumpC(ixC^S)+jumpC(ixC^S)*&
            abs(tmp(ixC^S))
      else
         where(abs(tmp(ixC^S))>=smallaC(ixC^S))
            phiC(ixC^S)=-ljumpC(jxC^S)-ljumpC(ixC^S)+jumpC(ixC^S)*&
               abs(tmp(ixC^S))
         elsewhere
            phiC(ixC^S)=-ljumpC(jxC^S)-ljumpC(ixC^S)+jumpC(ixC^S)*&
               (half*smallaC(ixC^S)+half/smallaC(ixC^S)*tmp(ixC^S)**2)
         endwhere
      endif
      !extra -(a*lambda)**2*delta
   case('testing')
      !phiC:=jumpC, thus R*L*(w(jx)-w(ix))==w(jx)-w(ix) can be tested
      !Use tvdmc, otherwise the second order correction is applied
      phiC(ixC^S)=jumpC(ixC^S)
   case default
      call die('Error in TVDLimit: Unknown TVD type='//typetvd)
end select
!SHIFT END

if(oktest)write(*,*)'GetPhi Final phiC:',&
    phiC(ixtest^D)

return
end

!=============================================================================
subroutine entropyfix(ix^L,il,aL,aR,a,smalla)

! Apply entropyfix based on typeentropy(il),aL,aR, and a
! Calculate "smalla" (Harten,Powell) or modify "a" (ratio)
!!! tmp and tmp2 are not to be used in this subroutine

include 'vacdef.f'

integer:: ix^L,il
double precision, dimension(ixG^T):: aL,aR,a,smalla
!-----------------------------------------------------------------------------

select case(typeentropy(il))
case('harten')
   smalla(ix^S)=max(zero,a(ix^S)-aL(ix^S),aR(ix^S)-a(ix^S))
case('powell')
   smalla(ix^S)=max(zero,two*(aR(ix^S)-aL(ix^S)))
case('ratio')
   where(aL(ix^S)<zero .and. aR(ix^S)>zero)&
      a(ix^S)=a(ix^S)-2*aR(ix^S)*aL(ix^S)/(aR(ix^S)-aL(ix^S))
case('yee')
   ! This has been done in geteigenjump already
case('nul')
   ! No entropyfix is applied
case default
   call die('No such type of entropy fix:'//typeentropy(il))
end select

return
end

!=============================================================================
! end module vactvd
!##############################################################################
