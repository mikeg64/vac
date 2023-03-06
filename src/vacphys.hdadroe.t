!##############################################################################
! module vacphys.hdadroe - subroutines for Roe-type Rieman solver for adiab HD
!=============================================================================
subroutine average(wL,wR,ix^L,iws,idim,wroe)

! This subroutine has been added by M. Nauta
! Calculate the mean as indicated in: A Riemann solver for barotropic flow,
! by P. Glaister in JCP 93, 477-480 (1991)
! rho -> rho, m -> v

include 'vacdef.f'

integer:: ix^L,idim,iws(niw_)
double precision, dimension(ixG^T,nw):: wL,wR,wroe
integer:: idir
double precision:: csound(ixG^T)
common /roe/ csound
!!!There is a similar qsimple parameter in the next geteigenjump subroutine
LOGICAL,PARAMETER:: qsimple=.true.
DOUBLE PRECISION,PARAMETER:: qsmall=1.D-6
!-----------------------------------------------------------------------------

if(qsimple)then
   ! This is the simple arithmetic average

   wroe(ix^S,rho_)=half*(wL(ix^S,rho_)+wR(ix^S,rho_))
   {^C& wroe(ix^S,m^C_)=&
      half*(wL(ix^S,m^C_)/wL(ix^S,rho_)+wR(ix^S,m^C_)/wR(ix^S,rho_)) \}

   csound(ix^S)=sqrt(eqpar(adiab_)*eqpar(gamma_)*&
      wroe(ix^S,rho_)**(eqpar(gamma_)-1))

   return
endif

! Calculate the Roe-average, use tmp,tmp2 to hold sqrt(rhoL),sqrt(rhoR)
tmp(ix^S)=sqrt(wL(ix^S,rho_))
tmp2(ix^S)=sqrt(wR(ix^S,rho_))
!The averaged density is sqrt(rhoL*rhoR)
wroe(ix^S,rho_)=tmp(ix^S)*tmp2(ix^S)

!Now the more useful ratio is put into tmp
tmp(ix^S)=tmp(ix^S)/tmp2(ix^S)
! Roe-average velocities
do idir=1,ndir
   wroe(ix^S,m0_+idir)=(wL(ix^S,m0_+idir)/wL(ix^S,rho_)*tmp(ix^S)+&
                        wR(ix^S,m0_+idir)/wR(ix^S,rho_))/(one+tmp(ix^S))
end do

! To avoid division by 0 check if densities are different enough
where(abs(wL(ix^S,rho_)-wR(ix^S,rho_))<=qsmall*(wL(ix^S,rho_)+wR(ix^S,rho_)))
   csound(ix^S)=sqrt(eqpar(adiab_)*eqpar(gamma_)*&
         wroe(ix^S,rho_)**(eqpar(gamma_)-1))
elsewhere
   csound(ix^S)=sqrt(eqpar(adiab_)*(wR(ix^S,rho_)**eqpar(gamma_)-&
              wL(ix^S,rho_)**eqpar(gamma_))/(wR(ix^S,rho_)-wL(ix^S,rho_)))
end where

return
end

!=============================================================================
subroutine geteigenjump(wL,wR,wroe,ix^L,il,idim,smalla,a,jump)

! This subroutine has been added by M. Nauta
! based on P. Glaister JCP 93, 477-480 (1991)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated at wroe. 
! jump(il)=Sum_il L(il,iw)*(w(jx,iw)-w(ix,iw))

include 'vacdef.f'

integer:: ix^L,il,idim,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T)   :: smalla,a,jump,csound
common /roe/ csound
!!!There is a similar qsimple parameter in the previous average subroutine
LOGICAL,PARAMETER:: qsimple=.true.
!-----------------------------------------------------------------------------

oktest=index(teststr,'geteigenjump')>=1
if(oktest)write(*,*)'GetEigenJump wL,wR:',&
   wL(ixtest^D,iwtest),wR(ixtest^D,iwtest)

if(qsimple)then
   ! This is the original simple Roe-solver

   select case(il)
   case(soundRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+csound(ix^S)
      jump(ix^S)=half*((one-wroe(ix^S,m0_+idim)/csound(ix^S))*&
               (wR(ix^S,rho_)-wL(ix^S,rho_))&
              +(wR(ix^S,m0_+idim)-wL(ix^S,m0_+idim))/csound(ix^S))
   case(soundLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-csound(ix^S)
      jump(ix^S)=half*((one+wroe(ix^S,m0_+idim)/csound(ix^S))*&
               (wR(ix^S,rho_)-wL(ix^S,rho_))&
               -(wR(ix^S,m0_+idim)-wL(ix^S,m0_+idim))/csound(ix^S))
   case default
      ! Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=-wroe(ix^S,m0_+idir)*(wR(ix^S,rho_)-wL(ix^S,rho_))&
               +(wR(ix^S,m0_+idir)-wL(ix^S,m0_+idir))
   end select
else
   ! This is the Roe solver by Glaister 
   select case(il)
   case(soundRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+csound(ix^S)
      jump(ix^S)=half*((wR(ix^S,rho_)-wL(ix^S,rho_))+&
              wroe(ix^S,rho_)/csound(ix^S)*(wR(ix^S,m0_+idim)/wR(ix^S,rho_)-&
              wL(ix^S,m0_+idim)/wL(ix^S,rho_)))
   case(soundLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-csound(ix^S)
      jump(ix^S)=half*((wR(ix^S,rho_)-wL(ix^S,rho_))-&
              wroe(ix^S,rho_)/csound(ix^S)*(wR(ix^S,m0_+idim)/wR(ix^S,rho_)-&
              wL(ix^S,m0_+idim)/wL(ix^S,rho_)))
   case default
      ! Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=wroe(ix^S,rho_)*(wR(ix^S,m0_+idir)/wR(ix^S,rho_)-&
               wL(ix^S,m0_+idir)/wL(ix^S,rho_))
   end select
endif

! Calculate "smalla" or modify "a" based on the "typeentropy" switch
! Use tmp and tmp2 for the left and right eigenvalues if needed
select case(typeentropy(il))
case('yee')
   ! Based on Yee JCP 68,151 eq 3.23
   smalla(ix^S)=entropycoef(il)
case('harten','powell')
   ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
   select case(il)
   case(soundRW_)
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)&
         + sqrt(eqpar(adiab_)*eqpar(gamma_)*wL(ix^S,rho_)**(eqpar(gamma_)-1))
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)&
         + sqrt(eqpar(adiab_)*eqpar(gamma_)*wR(ix^S,rho_)**(eqpar(gamma_)-1))
   case(soundLW_)
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)&
         - sqrt(eqpar(adiab_)*eqpar(gamma_)*wL(ix^S,rho_)**(eqpar(gamma_)-1))
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)&
         - sqrt(eqpar(adiab_)*eqpar(gamma_)*wR(ix^S,rho_)**(eqpar(gamma_)-1))
   case default
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)
   end select
end select

call entropyfix(ix^L,il,tmp,tmp2,a,smalla)

return
end

!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'vacdef.f'

integer::          ix^L,iw,il,idim,idir
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq,csound
common /roe/ csound
!-----------------------------------------------------------------------------

oktest=index(teststr,'rtimes')>=1

if(iw==rho_)then
   select case(il)
   case(soundRW_,soundLW_)
      rq(ix^S)=q(ix^S)
   case default
      rq(ix^S)=zero
   end select
else if(iw==m0_+idim)then
   select case(il)
   case(soundRW_)
      rq(ix^S)=q(ix^S)*(wroe(ix^S,m0_+idim)+csound(ix^S))
   case(soundLW_)
      rq(ix^S)=q(ix^S)*(wroe(ix^S,m0_+idim)-csound(ix^S))
   case default
      rq(ix^S)=zero
   end select
else
   select case(il)
   case(soundRW_,soundLW_)
      rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
   case default
      !Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      if(iw==m0_+idir) then
         rq(ix^S)=q(ix^S)
      else
         rq(ix^S)=zero
      endif
   end select
endif

return
end

!=============================================================================
! end module vacphys.hdadroe
!##############################################################################
