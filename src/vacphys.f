!##############################################################################
! module vacphys - hdadiab

!=============================================================================
subroutine addsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------

return
end
!=============================================================================
!=============================================================================
subroutine process(count,idimmin,idimmax,w)

! Process w before it is advected in directions idim^LIM, or before save
! count=1 and 2 for first and second (half step) processing during advection
! count=ifile+2 for saving results into the file indexed by ifile

include 'vacdef.f'

integer:: count,idimmin,idimmax
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

!-----------------------------------------------------------------------------

return
end
!=============================================================================
!##############################################################################
! module vacphys.hdadroe - subroutines for Roe-type Rieman solver for adiab HD
!=============================================================================
subroutine average(wL,wR,ixmin1,ixmin2,ixmax1,ixmax2,iws,idim,wroe)

! This subroutine has been added by M. Nauta
! Calculate the mean as indicated in: A Riemann solver for barotropic flow,
! by P. Glaister in JCP 93, 477-480 (1991)
! rho -> rho, m -> v

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim,iws(niw_)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wroe
integer:: idir
double precision:: csound(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
common /roe/ csound
!!!There is a similar qsimple parameter in the next geteigenjump subroutine
LOGICAL,PARAMETER:: qsimple=.true.
DOUBLE PRECISION,PARAMETER:: qsmall=1.D-6
!-----------------------------------------------------------------------------

if(qsimple)then
   ! This is the simple arithmetic average

   wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=half*(wL(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
    wroe(ixmin1:ixmax1,ixmin2:ixmax2,m1_)=half*(wL(ixmin1:ixmax1,&
       ixmin2:ixmax2,m1_)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
       +wR(ixmin1:ixmax1,ixmin2:ixmax2,m1_)/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
       rho_)) 
     wroe(ixmin1:ixmax1,ixmin2:ixmax2,m2_)=half*(wL(ixmin1:ixmax1,&
        ixmin2:ixmax2,m2_)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
        +wR(ixmin1:ixmax1,ixmin2:ixmax2,m2_)/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
        rho_)) 

   csound(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(eqpar(adiab_)*eqpar(gamma_)&
      *wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)-1))

   return
endif

! Calculate the Roe-average, use tmp,tmp2 to hold sqrt(rhoL),sqrt(rhoR)
tmp(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
!The averaged density is sqrt(rhoL*rhoR)
wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=tmp(ixmin1:ixmax1,ixmin2:ixmax2)&
   *tmp2(ixmin1:ixmax1,ixmin2:ixmax2)

!Now the more useful ratio is put into tmp
tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(ixmin1:ixmax1,ixmin2:ixmax2)&
   /tmp2(ixmin1:ixmax1,ixmin2:ixmax2)
! Roe-average velocities
do idir=1,ndir
   wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)=(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
      m0_+idir)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)*tmp(ixmin1:ixmax1,&
      ixmin2:ixmax2)+wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)&
      /wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))/(one+tmp(ixmin1:ixmax1,&
      ixmin2:ixmax2))
end do

! To avoid division by 0 check if densities are different enough
where(abs(wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-wR(ixmin1:ixmax1,ixmin2:ixmax2,&
   rho_))<=qsmall*(wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+wR(ixmin1:ixmax1,&
   ixmin2:ixmax2,rho_)))
   csound(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(eqpar(adiab_)*eqpar(gamma_)&
      *wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)-1))
elsewhere
   csound(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(eqpar(adiab_)*(wR(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)**eqpar(gamma_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
      **eqpar(gamma_))/(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
      -wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)))
end where

return
end

!=============================================================================
subroutine geteigenjump(wL,wR,wroe,ixmin1,ixmin2,ixmax1,ixmax2,il,idim,smalla,&
   a,jump)

! This subroutine has been added by M. Nauta
! based on P. Glaister JCP 93, 477-480 (1991)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated at wroe. 
! jump(il)=Sum_il L(il,iw)*(w(jx,iw)-w(ix,iw))

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,il,idim,idir
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wroe
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   :: smalla,a,jump,&
   csound
common /roe/ csound
!!!There is a similar qsimple parameter in the previous average subroutine
LOGICAL,PARAMETER:: qsimple=.true.
!-----------------------------------------------------------------------------

oktest=index(teststr,'geteigenjump')>=1
if(oktest)write(*,*)'GetEigenJump wL,wR:',wL(ixtest1,ixtest2,iwtest),&
   wR(ixtest1,ixtest2,iwtest)

if(qsimple)then
   ! This is the original simple Roe-solver

   select case(il)
   case(soundRW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)+csound(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half*((one-wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,m0_+idim)/csound(ixmin1:ixmax1,ixmin2:ixmax2))&
         *(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_))+(wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim))/csound&
         (ixmin1:ixmax1,ixmin2:ixmax2))
   case(soundLW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-csound(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half*((one+wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,m0_+idim)/csound(ixmin1:ixmax1,ixmin2:ixmax2))&
         *(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_))-(wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim))/csound&
         (ixmin1:ixmax1,ixmin2:ixmax2))
   case default
      ! Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=-wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idir)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_))+(wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idir)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir))
   end select
else
   ! This is the Roe solver by Glaister 
   select case(il)
   case(soundRW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)+csound(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half*((wR(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))+wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_)/csound(ixmin1:ixmax1,ixmin2:ixmax2)&
         *(wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)))
   case(soundLW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-csound(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half*((wR(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))-wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_)/csound(ixmin1:ixmax1,ixmin2:ixmax2)&
         *(wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)))
   case default
      ! Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)/wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idir)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
   end select
endif

! Calculate "smalla" or modify "a" based on the "typeentropy" switch
! Use tmp and tmp2 for the left and right eigenvalues if needed
select case(typeentropy(il))
case('yee')
   ! Based on Yee JCP 68,151 eq 3.23
   smalla(ixmin1:ixmax1,ixmin2:ixmax2)=entropycoef(il)
case('harten','powell')
   ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
   select case(il)
   case(soundRW_)
      tmp(ixmin1:ixmax1,ixmin2:ixmax2) =wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+ sqrt(eqpar(adiab_)&
         *eqpar(gamma_)*wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)&
         -1))
      tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+ sqrt(eqpar(adiab_)&
         *eqpar(gamma_)*wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)&
         -1))
   case(soundLW_)
      tmp(ixmin1:ixmax1,ixmin2:ixmax2) =wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)- sqrt(eqpar(adiab_)&
         *eqpar(gamma_)*wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)&
         -1))
      tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)- sqrt(eqpar(adiab_)&
         *eqpar(gamma_)*wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)&
         -1))
   case default
      tmp(ixmin1:ixmax1,ixmin2:ixmax2) =wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
      tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
   end select
end select

call entropyfix(ixmin1,ixmin2,ixmax1,ixmax2,il,tmp,tmp2,a,smalla)

return
end

!=============================================================================
subroutine rtimes(q,wroe,ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,rq)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,idir
double precision:: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: q,rq,csound
common /roe/ csound
!-----------------------------------------------------------------------------

oktest=index(teststr,'rtimes')>=1

if(iw==rho_)then
   select case(il)
   case(soundRW_,soundLW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)
   case default
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
   end select
else if(iw==m0_+idim)then
   select case(il)
   case(soundRW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
         *(wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)+csound(ixmin1:ixmax1,&
         ixmin2:ixmax2))
   case(soundLW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
         *(wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)-csound(ixmin1:ixmax1,&
         ixmin2:ixmax2))
   case default
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
   end select
else
   select case(il)
   case(soundRW_,soundLW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
         *wroe(ixmin1:ixmax1,ixmin2:ixmax2,iw)
   case default
      !Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      if(iw==m0_+idir) then
         rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)
      else
         rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
      endif
   end select
endif

return
end

!=============================================================================
! end module vacphys.hdadroe
!##############################################################################
!=============================================================================
subroutine physini

! Tell VAC which variables are vectors
! Set default values for entropycoef for typeentropy='yee'

include 'vacdef.f'

integer:: il
!-----------------------------------------------------------------------------

iw_vector(1)=m0_

! The values of the constants are taken from Ryu & Jones ApJ 442, 228
! assuming that it works for HD just as well as for MHD.
do il=1,nw
   select case(il)
   case(soundRW_,soundLW_)
      entropycoef(il)= 0.2
   case default
      entropycoef(il)= -one
   end select
end do

return
end

!=============================================================================
subroutine getdt(w,ixmin1,ixmin2,ixmax1,ixmax2)

! No extra restriction on dt

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
!-----------------------------------------------------------------------------

return
end

!=============================================================================
subroutine getflux(w,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,f,transport)

! Calculate non-transport flux f_idim[iw] within ix^L.
! Set transport=.true. if a transport flux should be added

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2,iw,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),f(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
logical::          transport
!-----------------------------------------------------------------------------

if(iw==m0_+idim)then
   ! f_i[m_i]=v_i*m_i+p and p=adiab*rho**gamma
   f(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(adiab_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
      rho_)**eqpar(gamma_)
else
   f(ixmin1:ixmax1,ixmin2:ixmax2)=zero
endif

transport=.true.

return
end

!=============================================================================
subroutine addgeometry(qdt,ixmin1,ixmin2,ixmax1,ixmax2,iws,w,wnew)

! Add geometrical source terms to wnew

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2,iws(niw_),ix,iiw,iw,idir
double precision:: qdt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'addgeometry')>=1

select case(typeaxial)
   case('slab','test')
       ! No source terms in slab symmetry
   case('nozzle')
      if(ndir>1)call die('Nozzle geometry for ndir>1 is not implemented')
      do iiw=1,iws(niw_); iw=iws(iiw)
         if(iw==mr_)then
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(adiab_)*w(ixmin1:ixmax1,&
               ixmin2:ixmax2,rho_)**eqpar(gamma_)
            forall(ix= ixmin1:ixmax1) wnew(ix,ixmin2:ixmax2,iw)&
               =wnew(ix,ixmin2:ixmax2,iw) +qdt*tmp(ix,ixmin2:ixmax2)&
               *areaside(ix)
         endif
      enddo
   case('sphere')
      do iiw=1,iws(niw_); iw=iws(iiw)
         select case(iw)
         ! s[mr]=pthermal*2/r+(mphi**2/rho+mtheta**2/rho)/r
         case(mr_)
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(adiab_)*w(ixmin1:ixmax1,&
               ixmin2:ixmax2,rho_)**eqpar(gamma_)
            ! This discretization maintains an exact equilibrium for p=const.
            ! areaside=(areaCi-areaCh)/areadx
            forall(ix= ixmin1:ixmax1) wnew(ix,ixmin2:ixmax2,iw)&
               =wnew(ix,ixmin2:ixmax2,iw) +qdt*tmp(ix,ixmin2:ixmax2)&
               *areaside(ix)
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=zero
            do idir=2,ndir
               tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(ixmin1:ixmax1,&
                  ixmin2:ixmax2)+w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
                  +idir)**2/w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
            end do
         ! s[mphi]=(-mphi*mr/rho)/radius and phi-->theta
         case(m2_)
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)= -w(ixmin1:ixmax1,ixmin2:ixmax2,&
               m2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,mr_)/w(ixmin1:ixmax1,&
               ixmin2:ixmax2,rho_) 
         end select
         ! Divide by radius and add to wnew for variables other than rho
         if(iw/=rho_)wnew(ixmin1:ixmax1,ixmin2:ixmax2,iw)=wnew(ixmin1:ixmax1,&
            ixmin2:ixmax2,iw)+qdt*tmp(ixmin1:ixmax1,ixmin2:ixmax2)&
            /x(ixmin1:ixmax1,ixmin2:ixmax2,r_)
         if(oktest.and.iw==iwtest)write(*,*)'Geometrical source:',tmp(ixtest1,&
            ixtest2)
      end do
   case('cylinder')
      do iiw=1,iws(niw_); iw=iws(iiw)
         select case(iw)
         ! s[mr]=(pthermal+mphi**2/rho)/radius
         case(mr_)
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(adiab_)*w(ixmin1:ixmax1,&
               ixmin2:ixmax2,rho_)**eqpar(gamma_)
            if(.not.gencoord)then
               ! For nonuniform Cartesian grid this provides hydrostatic equil.
               forall(ix= ixmin1:ixmax1) wnew(ix,ixmin2:ixmax2,iw)&
                  =wnew(ix,ixmin2:ixmax2,iw) +qdt*tmp(ix,ixmin2:ixmax2)&
                  *areaside(ix)
               tmp(ixmin1:ixmax1,ixmin2:ixmax2)=zero
            endif

         end select
         ! Divide by radius and add to wnew
         if(iw==mr_.or.(iw==mphi_.and..not.angmomfix)) wnew(ixmin1:ixmax1,&
            ixmin2:ixmax2,iw)=wnew(ixmin1:ixmax1,ixmin2:ixmax2,iw)&
            +qdt*tmp(ixmin1:ixmax1,ixmin2:ixmax2)/x(ixmin1:ixmax1,&
            ixmin2:ixmax2,r_)
         if(oktest.and.iw==iwtest)write(*,*)'Geometrical source:',tmp(ixtest1,&
            ixtest2)
      end do
end select

return
end

!=============================================================================

!=============================================================================
subroutine keeppositive_rho(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Keep density positive

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
logical:: toosmallr
!-----------------------------------------------------------------------------

! Overwrite smallrho if it is not yet defined
if(smallrho== -one)then
    smallrho=minval(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_))*smallrhocoeff
    
endif

! Check for very small density
toosmallr=any(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)<max(zero,smallrho))

if(toosmallr)then
   nerror(toosmallr_)=nerror(toosmallr_)+1
   if(nerror(toosmallr_)==1)then
      write(*,'(a,i2,a,i7)')'Too small density (code=',toosmallr_,') at it=',&
         it
      write(*,*)'Value < smallrho: ',minval(w(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)),smallrho
!     write(*,*)'Location: ',minloc(w(ix^S,rho_)) !F77_
      
   endif
   if(smallrho>zero)w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=max(w(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_),smallrho)
endif

return
end
!=============================================================================

subroutine keeppositive(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Keep density (and consequently pressure) positive

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------

if(vacuumrho<zero)then
   ! Keep density positive
   call keeppositive_rho(ixmin1,ixmin2,ixmax1,ixmax2,w)
else
   ! Where rho is small use vacuum state: rho=vacuumrho, v=0
   where(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)<smallrho)
      w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=vacuumrho
      w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)=zero
      w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)=zero;
   endwhere
endif

return
end

!=============================================================================
! end module vacphys - hdadiab
!##############################################################################


