!##############################################################################
! module vactvd
! Subroutines for TVD, TVD-MacCormack and TVD-MUSCL schemes
!=============================================================================
subroutine tvdlimit(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iws,idimmin,idimmax,w,qt,wnew)

! Limit the iws flow variables in wnew. 
! Limiting is based on w or wnew depending of "typelimited".
! Create the left and right values by shifting then call tvdlimit2.
! "method" can be temporally second order 'tvd' or first order 'tvd1'
! For the 2nd order 'tvd' source terms are also added to second order accuracy

include 'vacdef.f'

character*10 :: method
double precision:: qdt,qt
integer:: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   iws(niw_),idimmin,idimmax
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: w,wnew,wR
!!! wL
integer:: ixICmin1,ixICmin2,ixICmax1,ixICmax2,jxICmin1,jxICmin2,jxICmax1,&
   jxICmax2,iiw,iw,idim

logical:: firstsweep,lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

oktest=index(teststr,'tvd')>=1
if(oktest)write(*,*)'TVDLimit w,wnew,method:',w(ixtest1,ixtest2,iwtest),&
   wnew(ixtest1,ixtest2,iwtest),method

if(typelimited=='predictor')then
   !!!Side effect, w is overwritten, but after TVDlimit w is not needed anymore
   call getboundary(qt,1,nw,idimmin,idimmax,wnew)
   do iiw=1,iws(niw_); iw=iws(iiw)
      w(ixImin1:ixImax1,ixImin2:ixImax2,iw)=wnew(ixImin1:ixImax1,&
         ixImin2:ixImax2,iw)
   end do
end if

!!!Side effect, w is overwritten
!!!if(typelimited=='previous')w(ixI^S,1:nw)=wold(ixI^S,1:nw)

do idim= idimmin,idimmax
   ixICmax1=ixOmax1+kr(idim,1);ixICmax2=ixOmax2+kr(idim,2)
   ixICmin1=ixOmin1-2*kr(idim,1);ixICmin2=ixOmin2-2*kr(idim,2);

   !!! wL(ixIC^S,1:nw)=w(ixIC^S,1:nw)
   !SHIFT
   jxICmin1=ixICmin1+kr(idim,1);jxICmin2=ixICmin2+kr(idim,2)
   jxICmax1=ixICmax1+kr(idim,1);jxICmax2=ixICmax2+kr(idim,2);
   !SHIFT BEGIN
   wR(ixICmin1:ixICmax1,ixICmin2:ixICmax2,1:nw)=w(jxICmin1:jxICmax1,&
      jxICmin2:jxICmax2,1:nw)
   !SHIFT END
   
   call tvdlimit2(method,qdt,ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,iws,idim,w,wR,wnew) !!!w-->wL
enddo

if(method=='tvd'.and.lastsweep.and.index(teststr,'noryuprep')<1&
   .and.(typeaxial=='cylinder'.or.(typeaxial=='slab'.and.sourceunsplit)))then

   if(oktest)write(*,*)'Adding Ryu style (geometrical) sources'

   ! Following Ryu et al (APJ 452,364) (geometrical) sources are added 
   ! with second order accuracy:
   ! whalf=(wold+whydro+S)/2     eq.A22 corrected by including S
   ! wnew=whydro+S(whalf)        eq.A21

   !!!Side effect: We overwrite w, but that is OK, it is not needed anymore
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw)=wnew(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nw)

   if(typeaxial/='slab')call addgeometry(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
      iws,wold,w)
   if(sourceunsplit)call addsource2(qdt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,t+qdt,wold,t+qdt,w)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw)=half*(wold(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nw)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw))

   if(typeaxial/='slab')call addgeometry(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
      iws,w,wnew)
   if(sourceunsplit)call addsource2(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,t+qdt,w,t+qdt,wnew)
endif

return
end

!=============================================================================
subroutine tvdlimit2(method,qdt,ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iws,idim,wL,wR,wnew)

! Limit the iws flow variables in wnew according to typetvd. 
! wroeC is based on wL and wR.
! If method=='tvd' an extra adtdx**2*jumpC is added to phiC for 2nd order
! accuracy in time.

include 'vacdef.f'

character*10 :: method
double precision:: qdt
integer:: ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   iws(niw_),idim
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wnew,wroeC
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   :: phiC,rphiC,&
   jumpC,adtdxC,smallaC
double precision:: courantmax,wtest(nw)
integer:: hxOmin1,hxOmin2,hxOmax1,hxOmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
   jxCmin1,jxCmin2,jxCmax1,jxCmax2,jxICmin1,jxICmin2,jxICmax1,jxICmax2,iiw,iw,&
   il
!-----------------------------------------------------------------------------

oktest=index(teststr,'tvdlimit')>=1
if(oktest)write(*,*)'TVDLimit2 wL,wR,wnew:',wL(ixtest1,ixtest2,iwtest),&
   wR(ixtest1,ixtest2,iwtest),wnew(ixtest1,ixtest2,iwtest)
if(typetvd=='testing')wtest(1:nw)=0.

hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2; 

jxCmin1=ixCmin1+kr(idim,1);jxCmin2=ixCmin2+kr(idim,2)
jxCmax1=ixCmax1+kr(idim,1);jxCmax2=ixCmax2+kr(idim,2);
jxICmin1=ixICmin1+kr(idim,1);jxICmin2=ixICmin2+kr(idim,2)
jxICmax1=ixICmax1+kr(idim,1);jxICmax2=ixICmax2+kr(idim,2);

call average(wL,wR,ixICmin1,ixICmin2,ixICmax1,ixICmax2,iws,idim,wroeC)

! A loop on characteristic variables to calculate the dissipative flux phiC.
do il=1,nw
   !Calculate the jump in the il-th characteristic variable: L(wroe)*dw
   call geteigenjump(wL,wR,wroeC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,&
      smallaC,adtdxC,jumpC)

   if(oktest)write(*,*)'wL(ix,il),wR(ix,il):',wL(ixtest1,ixtest2,il),&
      wR(ixtest1,ixtest2,il)
   if(oktest)write(*,*)'L*dwC,a:',jumpC(ixtest1,ixtest2),adtdxC(ixtest1,&
      ixtest2)

   ! Normalize the eigenvalue "a" (and its limit "smalla" if needed):
   if(gencoord)then
      adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=adtdxC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2)*qdt*surfaceC(ixICmin1:ixICmax1,ixICmin2:ixICmax2,&
         idim)*two/(dvolume(ixICmin1:ixICmax1,ixICmin2:ixICmax2)&
         +dvolume(jxICmin1:jxICmax1,jxICmin2:jxICmax2))
      if(typeentropy(il)=='harten' .or. typeentropy(il)=='powell')smallaC&
         (ixICmin1:ixICmax1,ixICmin2:ixICmax2)=smallaC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2)*qdt*surfaceC(ixICmin1:ixICmax1,ixICmin2:ixICmax2,&
         idim)*two/(dvolume(ixICmin1:ixICmax1,ixICmin2:ixICmax2)&
         +dvolume(jxICmin1:jxICmax1,jxICmin2:jxICmax2))
   else
      adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=adtdxC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2)*qdt*two/(dx(ixICmin1:ixICmax1,ixICmin2:ixICmax2,&
         idim)+dx(jxICmin1:jxICmax1,jxICmin2:jxICmax2,idim))
      if(typeentropy(il)=='harten' .or. typeentropy(il)=='powell')smallaC&
         (ixICmin1:ixICmax1,ixICmin2:ixICmax2)=smallaC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2)*qdt*two/(dx(ixICmin1:ixICmax1,ixICmin2:ixICmax2,&
         idim)+dx(jxICmin1:jxICmax1,jxICmin2:jxICmax2,idim))
   endif

   if(oktest)write(*,*)'qdt,adtdxC:',qdt,adtdxC(ixtest1,ixtest2)

   ! For potentially extreme eigenvalues calculate dtcourant for next step
   if((il==extremeLW_.or.il==extremeRW_).and.courantpar>zero&
      .and.istep==nstep.and.implpar<zero)then
      courantmax=maxval(abs(adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)))
      
      if(courantmax>one) then
         nerror(couranterr_)=nerror(couranterr_)+1
         if(nerror(couranterr_)==1)write(*,*)'Courant condition error (code=',&
            couranterr_,') at it=',it,' for direction',idim,' CFL=',courantmax
      endif
      if(courantmax>smalldouble)dtcourant(idim)=min(dtcourant(idim),qdt&
         *courantpar/courantmax)
   endif

   ! Calculate the flux limiter function phi
   call getphi(method,jumpC,adtdxC,smallaC,ixICmin1,ixICmin2,ixICmax1,&
      ixICmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,il,idim,phiC)

   if(oktest)write(*,*)'il,phiC:',il,phiC(ixtest1,ixtest2)

   !Multiply phiC by 2>acmcoef>0 (kappa in eq.2.21 Yee, Sandham, Djomehri)
   if(acmcoef(il)>=zero)phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=acmcoef(il)&
      *phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

   !Multiply by a term of Artificial Compression Method (ACM)
   !(theta_j+1/2 in eq.2.21 Yee, Sandham, Djomehri) if axmexpo>0 is set

   if(acmexpo>zero)call acmswitch(jumpC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,&
      phiC)

   if(oktest)write(*,*)'phiCstar:',phiC(ixtest1,ixtest2)

   !Add R(iw,il)*phiC(il) to each variable iw in wnew
   do iiw=1,iws(niw_); iw=iws(iiw)
      call rtimes(phiC,wroeC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,il,idim,rphiC)

      !Based on Eqs 5.6a,b in Yee
      rphiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=rphiC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*quarter*(dvolume(jxCmin1:jxCmax1,jxCmin2:jxCmax2)&
         +dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      if(gencoord.and.vectoriw(iw)>=0)then
         ! The dt=-one tells addflux_rotate that this is an artificial flux
         
                call die('Error: gencoord is off')
      else
         ! The dt=-one tells addflux that this is an artificial flux
         call addflux(-one,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,idim,rphiC,&
            ixOmin1,ixOmin2,ixOmax1,ixOmax2,rphiC,hxOmin1,hxOmin2,hxOmax1,&
            hxOmax2,wnew)
      endif

      if(oktest)write(*,*)'rphiCi,rphiCh:',rphiC(ixtest1,ixtest2)&
         /dvolume(ixtest1,ixtest2),rphiC(ixtest1-kr(idim,1),ixtest2&
         -kr(idim,2))/dvolume(ixtest1,ixtest2)
      if(typetvd=='testing')wtest(iw)=wtest(iw)+rphiC(ixtest1,ixtest2)
   end do  !iw
   if(oktest)write(*,*)'wnew:',wnew(ixtest1,ixtest2,iwtest)
end do     !il

if(typetvd=='testing')then
   write(*,*)'wtest                   wR-wL'
   do iiw=1,iws(niw_); iw=iws(iiw)
      write(*,*)wtest(iw),wR(ixtest1,ixtest2,iw)-wL(ixtest1,ixtest2,iw)
   end do
endif

return
end

!=============================================================================
subroutine getphi(method,jumpC,adtdxC,smallaC,ixICmin1,ixICmin2,ixICmax1,&
   ixICmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,il,idim,phiC)

! Calculate the dissipative flux from jumpC=L*dw and adtdx=eigenvalue*dt/dx.
! Add Lax-Wendroff type correction if method=='tvd'.
! Limit according to method and typetvd.

include 'vacdef.f'

character*10 :: method
integer:: ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: jumpC,ljumpC,&
   adtdxC,smallaC,phiC
integer:: jxCmin1,jxCmin2,jxCmax1,jxCmax2,ixmin1,ixmin2,ixmax1,ixmax2,hxmin1,&
   hxmin2,hxmax1,hxmax2,il,idim
!-----------------------------------------------------------------------------

oktest=index(teststr,'getphi')>=1
if(oktest)write(*,*)'Getphi jumpC,adtdxC:',jumpC(ixtest1,ixtest2),&
   adtdxC(ixtest1,ixtest2)

if(method=='tvdmu'.or.method=='tvdmu1')then
   ! In the MUSCL scheme phi=|a|*jump, apply entropy fix to it
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=abs(adtdxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   else
      where(abs(adtdxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>=smallaC&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2))
         phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=abs(adtdxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      elsewhere
         phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=half*(smallaC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2&
            /smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
      endwhere
   endif
   ! That's all for the MUSCL scheme
   if(oktest)write(*,*)'GetPhi Final phiC:',phiC(ixtest1,ixtest2)
   return
endif

if(method=='tvd')then
   !Entropy fix to |a|-a**2
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
     phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC(ixICmin1:ixICmax1,&
        ixICmin2:ixICmax2))-adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)**2
   else
      where(abs(adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2))&
         >=smallaC(ixICmin1:ixICmax1,ixICmin2:ixICmax2))
         phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2))-adtdxC(ixICmin1:ixICmax1,&
            ixICmin2:ixICmax2)**2
      elsewhere
         phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*smallaC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2)+(half/smallaC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2)-one)*adtdxC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2)**2
      endwhere
   endif
   if(oktest)write(*,*)'abs(nu)-nu**2:',phiC(ixtest1,ixtest2)
else
   !Entropy fix to |a|
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2))
   else
      where(abs(adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2))&
         >=smallaC(ixICmin1:ixICmax1,ixICmin2:ixICmax2))
         phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2))
      elsewhere
         phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*smallaC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2)+half/smallaC&
            (ixICmin1:ixICmax1,ixICmin2:ixICmax2)*adtdxC(ixICmin1:ixICmax1,&
            ixICmin2:ixICmax2)**2
      endwhere
   endif
   if(oktest)write(*,*)'abs(nu)-nu**2:',phiC(ixtest1,ixtest2)
endif

!SHIFT
jxCmin1=ixCmin1+kr(idim,1);jxCmin2=ixCmin2+kr(idim,2)
jxCmax1=ixCmax1+kr(idim,1);jxCmax2=ixCmax2+kr(idim,2);

hxmin1=ixICmin1;hxmin2=ixICmin2; hxmax1=ixICmax1-kr(idim,1)
hxmax2=ixICmax2-kr(idim,2);
!SHIFT MORE
ixmin1=hxmin1+kr(idim,1);ixmin2=hxmin2+kr(idim,2);ixmax1=hxmax1+kr(idim,1)
ixmax2=hxmax2+kr(idim,2);

!SHIFT BEGIN
select case(typetvd)
   case('symmetric')
      !eq.3.53 and eq.3.69d
      call dwlimiter3(jumpC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,&
         ljumpC)
      phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
         -ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      !extra (a*lambda)**2*delta
      if(method=='tvd')phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
         =phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   case('roe')
      !eq.3.52 with correction of sign of a
      select case(typelimiter(il))
      case('roe','superroe')
         call dwlimiterroe(adtdxC,jumpC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,il,&
            idim,ljumpC)
         phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)*(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            -ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      case default
         call dwlimiter2(jumpC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,&
            ljumpC)
         where(adtdxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)<=0)
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)*(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
               -ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
         elsewhere
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)*(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
               -ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
         end where
      end select
      !extra (a*lambda)**2*delta
      if(method=='tvd')phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
         =phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   case('sweby')
      !Sweby eqs.4.11-4.15, but no 0.5 ?!
      phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=phiC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2)*jumpC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)
      select case(typelimiter(il))
      case('roe','superroe')
         call dwlimiterroe(adtdxC,phiC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,il,&
            idim,ljumpC)
         phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      case default
         call dwlimiter2(phiC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,&
            ljumpC)
         where(adtdxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)<=0)
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)-ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
         elsewhere
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
         end where
      end select
      !extra (a*lambda)**2*delta
      if(method=='tvd')phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
         =phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   case('yee')
      !eq.3.51 with correction
      call dwlimiter2(jumpC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,&
         ljumpC)

      if(oktest)write(*,*)'Limiter jumpC,ljumpC(ix,jx):',jumpC(ixtest1,&
         ixtest2),ljumpC(ixtest1,ixtest2),ljumpC(ixtest1+kr(idim,1),ixtest2&
         +kr(idim,2))

      !Use phiC as 0.5*(|nu|-nu**2) eq.3.45e for tvd otherwise 0.5*|nu|
      phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=half*phiC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      !gamma*lambda eq.3.51d, use tmp to store agdtdxC
      where(abs(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>smalldouble)
         tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            *(ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2))/jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      elsewhere
         tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
      end where

      if(oktest)write(*,*)'agdtdxC:',tmp(ixtest1,ixtest2)

      !eq.3.51a
      if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
         phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-phiC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)&
            +ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))+abs(tmp(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      else
         where(abs(tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>=smallaC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2))
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-phiC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)&
               +ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))+abs(tmp&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
         elsewhere
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-phiC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)&
               +ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))+(half&
               *smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+half&
               /smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*tmp(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)**2)*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
         endwhere
      endif
   case('harten')
      !See Ryu, section 2.3
      !Use phiC as 0.5*(|nu|-nu**2)*jumpC eq.3.45b,e
      phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*phiC(ixICmin1:ixICmax1,&
         ixICmin2:ixICmax2)*jumpC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)
      call dwlimiter2(phiC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,ljumpC)

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

         tmp(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*(one&
            -abs(adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)))&
            *jumpC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)
         tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=sign(one,jumpC(ixmin1:ixmax1,&
            ixmin2:ixmax2))*max(zero,min(sign(one,jumpC(ixmin1:ixmax1,&
            ixmin2:ixmax2))*tmp(hxmin1:hxmax1,hxmin2:hxmax2),&
            abs(tmp(ixmin1:ixmax1,ixmin2:ixmax2))))
         where(abs(jumpC(ixmin1:ixmax1,ixmin2:ixmax2))+abs(jumpC&
            (hxmin1:hxmax1,hxmin2:hxmax2))>smalldouble)
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=abs(jumpC(ixmin1:ixmax1,&
               ixmin2:ixmax2)-jumpC(hxmin1:hxmax1,hxmin2:hxmax2))&
               /(abs(jumpC(ixmin1:ixmax1,ixmin2:ixmax2))+abs(jumpC&
               (hxmin1:hxmax1,hxmin2:hxmax2)))
         elsewhere
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=zero
         end where
         ljumpC(ixmin1:ixmax1,ixmin2:ixmax2)=ljumpC(ixmin1:ixmax1,&
            ixmin2:ixmax2)+tmp(ixmin1:ixmax1,ixmin2:ixmax2)&
            *tmp2(ixmin1:ixmax1,ixmin2:ixmax2)
      endif

      !gamma*lambda eq.3.45d, use tmp as agdtdxC
      where(abs(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>smalldouble)
         tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+(ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)&
            -ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))/jumpC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
      elsewhere
         tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
      end where
      if(oktest)write(*,*)'agdtdxC',tmp(ixtest1,ixtest2)
      !eq.3.45a with correction
      if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
         phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-ljumpC(jxCmin1:jxCmax1,&
            jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            +jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*abs(tmp(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2))
      else
         where(abs(tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>=smallaC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2))
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-ljumpC(jxCmin1:jxCmax1,&
               jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
               +jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*abs(tmp&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2))
         elsewhere
            phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-ljumpC(jxCmin1:jxCmax1,&
               jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
               +jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*(half&
               *smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+half&
               /smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*tmp(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)**2)
         endwhere
      endif
      !extra -(a*lambda)**2*delta
   case('testing')
      !phiC:=jumpC, thus R*L*(w(jx)-w(ix))==w(jx)-w(ix) can be tested
      !Use tvdmc, otherwise the second order correction is applied
      phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=jumpC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
   case default
      call die('Error in TVDLimit: Unknown TVD type='//typetvd)
end select
!SHIFT END

if(oktest)write(*,*)'GetPhi Final phiC:',phiC(ixtest1,ixtest2)

return
end

!=============================================================================
subroutine entropyfix(ixmin1,ixmin2,ixmax1,ixmax2,il,aL,aR,a,smalla)

! Apply entropyfix based on typeentropy(il),aL,aR, and a
! Calculate "smalla" (Harten,Powell) or modify "a" (ratio)
!!! tmp and tmp2 are not to be used in this subroutine

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,il
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: aL,aR,a,smalla
!-----------------------------------------------------------------------------

select case(typeentropy(il))
case('harten')
   smalla(ixmin1:ixmax1,ixmin2:ixmax2)=max(zero,a(ixmin1:ixmax1,&
      ixmin2:ixmax2)-aL(ixmin1:ixmax1,ixmin2:ixmax2),aR(ixmin1:ixmax1,&
      ixmin2:ixmax2)-a(ixmin1:ixmax1,ixmin2:ixmax2))
case('powell')
   smalla(ixmin1:ixmax1,ixmin2:ixmax2)=max(zero,two*(aR(ixmin1:ixmax1,&
      ixmin2:ixmax2)-aL(ixmin1:ixmax1,ixmin2:ixmax2)))
case('ratio')
   where(aL(ixmin1:ixmax1,ixmin2:ixmax2)<zero .and. aR(ixmin1:ixmax1,&
      ixmin2:ixmax2)>zero)a(ixmin1:ixmax1,ixmin2:ixmax2)=a(ixmin1:ixmax1,&
      ixmin2:ixmax2)-2*aR(ixmin1:ixmax1,ixmin2:ixmax2)*aL(ixmin1:ixmax1,&
      ixmin2:ixmax2)/(aR(ixmin1:ixmax1,ixmin2:ixmax2)-aL(ixmin1:ixmax1,&
      ixmin2:ixmax2))
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
