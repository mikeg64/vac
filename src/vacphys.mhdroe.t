!##############################################################################
! module vacphys.mhdroe - subroutines for Roe-type Riemann solver for MHD
!=============================================================================
subroutine average(wL,wR,ix^L,iws,idim,wroe)

! Eight-wave MHD Riemann solver. See Powell, Notes on the eigeinsystem, Gombosi
! Calculate the wroe average of primitive variables in wL and wR, assignment:
! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
!
! wL,wR,wroe are all interface centered quantities

include 'vacdef.f'

integer:: ix^L,idim,iws(niw_),idir,jdir,iw
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T):: cfast,cslow,afast,aslow,csound2,dp,rhodv
common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv

!-----------------------------------------------------------------------------

if(ndir==1)call die('MHD with d=11 is the same as HD')

oktest=index(teststr,'average')>=1
if(oktest)write(*,*)'Average wL,wR:',&
   wL(ixtest^D,iwtest),wR(ixtest^D,iwtest)

!Averaging primitive variables
wroe(ix^S,rho_)=half*(wL(ix^S,rho_)+wR(ix^S,rho_))
{^C&
wroe(ix^S,v^C_)=half*(wL(ix^S,m^C_)/wL(ix^S,rho_)+wR(ix^S,m^C_)/wR(ix^S,rho_))
wroe(ix^S,b^C_)=half*(wL(ix^S,b^C_)+wR(ix^S,b^C_))
\}

! Use afast and aslow for pressures pL and pR
call getpthermal(.true.,wL,ix^L,afast)
call getpthermal(.true.,wR,ix^L,aslow)
if(p_>0)wroe(ix^S,pp_)=half*(afast(ix^S)+aslow(ix^S))

if(oktest)write(*,*)'Calculate saved variables'

if(useprimitive.or.p_<1.or.eqpar(gamma_)<=zero)then
   ! dp=pR-pL
   dp(ix^S)=aslow(ix^S)-afast(ix^S)
else
   !CONSERVATIVE dp=(g-1)*(de-v*dm+0.5*v**2*drho-0.5*d(B**2))
   dp(ix^S)=(eqpar(gamma_)-1)*(wR(ix^S,ee_)-wL(ix^S,ee_)&
      -(^C&wroe(ix^S,m^C_)*(wR(ix^S,m^C_)-wL(ix^S,m^C_))+)&
      +half*(^C&wroe(ix^S,m^C_)**2+)*(wR(ix^S,rho_)-wL(ix^S,rho_))&
      -half*(^C&wR(ix^S,b^C_)**2-wL(ix^S,b^C_)**2+))
endif

!CONSERVATIVE rho*dv_idim=dm_idim-v_idim*drho
rhodv(ix^S)=wR(ix^S,m0_+idim)-wL(ix^S,m0_+idim)-&
          wroe(ix^S,m0_+idim)*(wR(ix^S,rho_)-wL(ix^S,rho_))

!Calculate csound2,cfast,cslow,alphafast and alphaslow

! get csound**2
call getcsound2prim(wroe,ix^L,csound2)

! aa=B**2/rho+a**2
cfast(ix^S)=( ^C&wroe(ix^S,b^C_)**2+ )/wroe(ix^S,rho_)+csound2(ix^S)

! cs**2=0.5*(aa+sqrt(aa**2-4*a**2*(b_i**2/rho)))
cslow(ix^S)=max(zero, half*(cfast(ix^S)-sqrt(cfast(ix^S)**2-&
   4d0*csound2(ix^S)*wroe(ix^S,b0_+idim)**2/wroe(ix^S,rho_))))

! cf**2=aa-cs**2
cfast(ix^S)=cfast(ix^S)-cslow(ix^S)

! alpha_f**2=(a**2-cs**2)/(cf**2-cs**2)
afast(ix^S)=(csound2(ix^S)-cslow(ix^S))/(cfast(ix^S)-cslow(ix^S))
afast(ix^S)=min(one,max(afast(ix^S),zero))

! alpha_s=sqrt(1-alpha_f**2)
aslow(ix^S)=sqrt(one-afast(ix^S))

! alpha_f=sqrt(alpha_f**2)
afast(ix^S)=sqrt(afast(ix^S))

! cf=sqrt(cf**2)
cfast(ix^S)=sqrt(cfast(ix^S))

! cs=sqrt(cs**2)
cslow(ix^S)=sqrt(cslow(ix^S))

if(oktest)write(*,*)'Average:rho,csound2,dp,rhodv',&
   wroe(ixtest^D,rho_),csound2(ixtest^D),dp(ixtest^D),rhodv(ixtest^D)
if(oktest)write(*,*)'Average:cf,cs,af,as',&
   cfast(ixtest^D),cslow(ixtest^D),afast(ixtest^D),aslow(ixtest^D)

!Replace the primitive variables with more useful quantities:
! rho -> sqrt(rho)
wroe(ix^S,rho_)=sqrt(wroe(ix^S,rho_))

! Avoid sgn(b_idim)==0
where(abs(wroe(ix^S,b0_+idim))<smalldouble)&
   wroe(ix^S,b0_+idim)=smalldouble
! B_idir,jdir -> beta_idir,jdir
idir=idim+1-ndir*(idim/ndir)
if(ndir==2)then
    where(wroe(ix^S,b0_+idir)>=zero)
       wroe(ix^S,b0_+idir)=one
    elsewhere
       wroe(ix^S,b0_+idir)=-one
    end where
else
    !beta_j=B_j/sqrt(B_i**2+B_j**2); beta_i=B_i/sqrt(B_i**2+B_j**2)
    jdir=idir+1-ndir*(idir/ndir)
    tmp(ix^S)=sqrt(wroe(ix^S,b0_+idir)**2+wroe(ix^S,b0_+jdir)**2)
    where(tmp(ix^S)>smalldouble)
       wroe(ix^S,b0_+idir)=wroe(ix^S,b0_+idir)/tmp(ix^S)
       wroe(ix^S,b0_+jdir)=wroe(ix^S,b0_+jdir)/tmp(ix^S)
    elsewhere
       wroe(ix^S,b0_+idir)=sqrt(half)
       wroe(ix^S,b0_+jdir)=sqrt(half)
    end where
endif

return
end

!=============================================================================
subroutine geteigenjump(wL,wR,wroe,ix^L,il,idim,smalla,a,jump)

! Calculate the il-th characteristic speed and the jump in the il-th
! characteristic variable in the idim direction within ixL.
! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
! variables. However part of the summation is done in advance and saved into
! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
! used in the entropy fix.
!
! All the variables are centered on the cell interface, thus the 
! "*C" notation is omitted for sake of brevity.

include 'vacdef.f'

integer:: ix^L,il,idim,idir,jdir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T)   :: smalla,a,jump
double precision, dimension(ixG^T)   ::cfast,cslow,afast,aslow,csound2,dp,rhodv
common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
double precision, dimension(ixG^T)   :: bdv,bdb
save bdv,bdb
double precision, dimension(ixG^T)   :: aL,aR,cs2L,cs2R,cs2ca2L,cs2ca2R
save cs2L,cs2R,cs2ca2L,cs2ca2R
!-----------------------------------------------------------------------------

oktest=index(teststr,'geteigenjump')>=1

idir=idim+1-ndir*(idim/ndir)
jdir=idir+1-ndir*(idir/ndir)

if(il==fastRW_)then
   !Fast and slow waves use bdv=sqrho**2*sign(bx)*(betay*dvy+betaz*dvz)
   !                        bdb=sqrho*a*          (betay*dBy+betaz*dBz)
   bdv(ix^S)=wroe(ix^S,b0_+idir)* &
      (wR(ix^S,m0_+idir)/wR(ix^S,rho_)-wL(ix^S,m0_+idir)/wL(ix^S,rho_))
   if(ndir==3)bdv(ix^S)=bdv(ix^S)+wroe(ix^S,b0_+jdir)* &
      (wR(ix^S,m0_+jdir)/wR(ix^S,rho_)-wL(ix^S,m0_+jdir)/wL(ix^S,rho_))
   bdv(ix^S)=bdv(ix^S)*sign(wroe(ix^S,rho_)**2,wroe(ix^S,b0_+idim))

   bdb(ix^S)=wroe(ix^S,b0_+idir)*(wR(ix^S,b0_+idir)-wL(ix^S,b0_+idir))
   if(ndir==3)bdb(ix^S)=bdb(ix^S)+&
              wroe(ix^S,b0_+jdir)*(wR(ix^S,b0_+jdir)-wL(ix^S,b0_+jdir))
   bdb(ix^S)=bdb(ix^S)*sqrt(csound2(ix^S))*wroe(ix^S,rho_)
   if(oktest)write(*,*)'rhobetadv,sqrhoabetadb:',&
      bdv(ixtest^D),bdb(ixtest^D)
endif

if(il==alfvRW_)then
   !Alfven waves use      bdv=0.5*sqrho**2*      (betaz*dvy-betay*dvz)
   !                      bdb=0.5*sqrho*sign(bx)*(betaz*dBy-betay*dBz)
   bdv(ix^S)=wroe(ix^S,b0_+jdir)* &
       (wR(ix^S,m0_+idir)/wR(ix^S,rho_)-wL(ix^S,m0_+idir)/wL(ix^S,rho_)) &
               -wroe(ix^S,b0_+idir)* &
       (wR(ix^S,m0_+jdir)/wR(ix^S,rho_)-wL(ix^S,m0_+jdir)/wL(ix^S,rho_))
   bdb(ix^S)=wroe(ix^S,b0_+jdir)*(wR(ix^S,b0_+idir)-wL(ix^S,b0_+idir)) &
               -wroe(ix^S,b0_+idir)*(wR(ix^S,b0_+jdir)-wL(ix^S,b0_+jdir))
   bdv(ix^S)=bdv(ix^S)*half*wroe(ix^S,rho_)**2
   bdb(ix^S)=bdb(ix^S)*half*sign(wroe(ix^S,rho_),wroe(ix^S,b0_+idim))
   if(oktest)write(*,*)'rhobetaXdv/2,sqrhobetaXdb/2:',&
      bdv(ixtest^D),bdb(ixtest^D)
endif

select case(il)
   case(fastRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+cfast(ix^S)
      jump(ix^S)=half/csound2(ix^S)*(&
         afast(ix^S)*(+cfast(ix^S)*rhodv(ix^S)+dp(ix^S))&
        +aslow(ix^S)*(-cslow(ix^S)*bdv(ix^S)+bdb(ix^S)))
   case(fastLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-cfast(ix^S)
      jump(ix^S)=half/csound2(ix^S)*(&
         afast(ix^S)*(-cfast(ix^S)*rhodv(ix^S)+dp(ix^S))&
        +aslow(ix^S)*(+cslow(ix^S)*bdv(ix^S)+bdb(ix^S)))
   case(slowRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+cslow(ix^S)
      jump(ix^S)=half/csound2(ix^S)*(&
         aslow(ix^S)*(+cslow(ix^S)*rhodv(ix^S)+dp(ix^S))&
        +afast(ix^S)*(+cfast(ix^S)*bdv(ix^S)-bdb(ix^S)))
   case(slowLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-cslow(ix^S)
      jump(ix^S)=half/csound2(ix^S)*(&
         aslow(ix^S)*(-cslow(ix^S)*rhodv(ix^S)+dp(ix^S))&
        +afast(ix^S)*(-cfast(ix^S)*bdv(ix^S)-bdb(ix^S)))
   case(entroW_)
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=wR(ix^S,rho_)-wL(ix^S,rho_)-dp(ix^S)/csound2(ix^S)
   case(diverW_)
      if(divbwave)then
         a(ix^S)=wroe(ix^S,m0_+idim)
         jump(ix^S)=wR(ix^S,b0_+idim)-wL(ix^S,b0_+idim)
      else
         a(ix^S)=zero
         jump(ix^S)=zero
      endif
   case(alfvRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+abs(wroe(ix^S,b0_+idim))/wroe(ix^S,rho_)
      jump(ix^S)=+bdv(ix^S)-bdb(ix^S)
   case(alfvLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-abs(wroe(ix^S,b0_+idim))/wroe(ix^S,rho_)
      jump(ix^S)=-bdv(ix^S)-bdb(ix^S)
end select

! Calculate "smalla" or modify "a" based on the "typeentropy" switch

select case(typeentropy(il))
case('yee')
   ! Based on Yee JCP 68,151 eq 3.23
   smalla(ix^S)=entropycoef(il)
case('harten','powell', 'ratio')
   ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
   ! Initialize left and right eigenvalues by velocities
   aL(ix^S)= wL(ix^S,m0_+idim)/wL(ix^S,rho_)
   aR(ix^S)= wR(ix^S,m0_+idim)/wR(ix^S,rho_)
   ! Calculate the final "aL" and "aR"
   select case(il)
   case(fastRW_)
      ! These quantities will be used for all the fast and slow waves
      ! Calculate soundspeed**2 and cs**2+ca**2.
      call getcsound2(wL,ix^L,cs2L)
      call getcsound2(wR,ix^L,cs2R)
      cs2ca2L(ix^S)=cs2L(ix^S)+(^C&wL(ix^S,b^C_)**2+)/wL(ix^S,rho_)
      cs2ca2R(ix^S)=cs2R(ix^S)+(^C&wR(ix^S,b^C_)**2+)/wR(ix^S,rho_)
      ! Save the discriminants into cs2L and cs2R
      cs2L(ix^S)=&
         sqrt(cs2ca2L(ix^S)**2-4*cs2L(ix^S)*wL(ix^S,b0_+idim)**2/wL(ix^S,rho_))
      cs2R(ix^S)=&
         sqrt(cs2ca2R(ix^S)**2-4*cs2R(ix^S)*wR(ix^S,b0_+idim)**2/wR(ix^S,rho_))

      ! The left and right eigenvalues for the fast wave going to right
      aL(ix^S)=aL(ix^S) + sqrt(half*(cs2ca2L(ix^S) + cs2L(ix^S)))
      aR(ix^S)=aR(ix^S) + sqrt(half*(cs2ca2R(ix^S) + cs2R(ix^S)))
   case(fastLW_)
      aL(ix^S)=aL(ix^S) - sqrt(half*(cs2ca2L(ix^S) + cs2L(ix^S)))
      aR(ix^S)=aR(ix^S) - sqrt(half*(cs2ca2R(ix^S) + cs2R(ix^S)))
   case(slowRW_)
      aL(ix^S)=aL(ix^S) + sqrt(half*(cs2ca2L(ix^S) - cs2L(ix^S)))
      aR(ix^S)=aR(ix^S) + sqrt(half*(cs2ca2R(ix^S) - cs2R(ix^S)))
   case(slowLW_)
      aL(ix^S)=aL(ix^S) - sqrt(half*(cs2ca2L(ix^S) - cs2L(ix^S)))
      aR(ix^S)=aR(ix^S) - sqrt(half*(cs2ca2R(ix^S) - cs2R(ix^S)))
   case(entroW_,diverW_)
      ! These propagate by the velocity
   case(alfvRW_)
      ! Store the Alfven speeds into cs2ca2L and cs2ca2R
      cs2ca2L(ix^S)=abs(wL(ix^S,b0_+idim))/sqrt(wL(ix^S,rho_))
      cs2ca2R(ix^S)=abs(wR(ix^S,b0_+idim))/sqrt(wR(ix^S,rho_))

      aL(ix^S)=aL(ix^S) + cs2ca2L(ix^S)
      aR(ix^S)=aR(ix^S) + cs2ca2R(ix^S)
   case(alfvLW_)
      aL(ix^S)=aL(ix^S) - cs2ca2L(ix^S)
      aR(ix^S)=aR(ix^S) - cs2ca2R(ix^S)
   end select
end select

call entropyfix(ix^L,il,aL,aR,a,smalla)

return
end

!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'vacdef.f'

integer::          ix^L,iw,il,idim,idir,jdir
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq
double precision, dimension(ixG^T):: cfast,cslow,afast,aslow,csound2,dp,rhodv
common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
double precision, dimension(ixG^T):: bv,v2a2
save:: bv,v2a2
!-----------------------------------------------------------------------------

oktest=index(teststr,'rtimes')>=1

idir=idim+1-ndir*(idim/ndir)
jdir=idir+1-ndir*(idir/ndir)

select case(iw)
  case(rho_)
    select case(il)
      case(fastRW_,fastLW_)
        rq(ix^S)=q(ix^S)*afast(ix^S)
      case(slowRW_,slowLW_)
        rq(ix^S)=q(ix^S)*aslow(ix^S)
      case(entroW_)
        rq(ix^S)=q(ix^S)
      case(diverW_,alfvRW_,alfvLW_)
        rq(ix^S)=zero
    end select
  case(e_)
    if(il==fastRW_)then
      if(eqpar(gamma_)>zero)then
         !!!IDEAL GAS
         ! Store 0.5*v**2+(2-gamma)/(gamma-1)*a**2
         v2a2(ix^S)=half*(^C&wroe(ix^S,m^C_)**2+)+ &
            (2-eqpar(gamma_))/(eqpar(gamma_)-1)*csound2(ix^S)
      else
         call die('Correct rTimes for NONIDEAL gas in src/vacphys.mhdroe.t')
         ! Express the partial derivative de/dp using wroe
         !v2a2(ix^S)=half*(^C&wroe(ix^S,m^C_)**2+)+ &
         !  (??dedp(ix^S)??-1)*csound2(ix^S)
      endif
      ! Store sgn(bx)*(betay*vy+betaz*vz) in bv
      bv(ix^S)=wroe(ix^S,b0_+idir)*wroe(ix^S,m0_+idir)
      if(ndir==3)bv(ix^S)=bv(ix^S)+&
                wroe(ix^S,b0_+jdir)*wroe(ix^S,m0_+jdir)
      bv(ix^S)=bv(ix^S)*sign(one,wroe(ix^S,b0_+idim))
      if(oktest)write(*,*)'v2/2+(2-g)/(g-1)a2,betav:',&
         v2a2(ixtest^D),bv(ixtest^D)
    else if(il==alfvRW_)then
      !Store betaz*vy-betay*vz in bv
      bv(ix^S)=(wroe(ix^S,b0_+jdir)*wroe(ix^S,m0_+idir)-&
                 wroe(ix^S,b0_+idir)*wroe(ix^S,m0_+jdir))
      if(oktest)write(*,*)'betaXv:',&
         bv(ixtest^D)
    endif
    select case(il)
      case(fastRW_)
        rq(ix^S)=q(ix^S)*(-aslow(ix^S)*cslow(ix^S)*bv(ix^S)+afast(ix^S)*&
              (v2a2(ix^S)+cfast(ix^S)*(cfast(ix^S)+wroe(ix^S,m0_+idim))))
      case(fastLW_)
        rq(ix^S)=q(ix^S)*(+aslow(ix^S)*cslow(ix^S)*bv(ix^S)+afast(ix^S)*&
              (v2a2(ix^S)+cfast(ix^S)*(cfast(ix^S)-wroe(ix^S,m0_+idim))))
      case(slowRW_)
        rq(ix^S)=q(ix^S)*(+afast(ix^S)*cfast(ix^S)*bv(ix^S)+aslow(ix^S)*&
              (v2a2(ix^S)+cslow(ix^S)*(cslow(ix^S)+wroe(ix^S,m0_+idim))))
      case(slowLW_)
        rq(ix^S)=q(ix^S)*(-afast(ix^S)*cfast(ix^S)*bv(ix^S)+aslow(ix^S)*&
              (v2a2(ix^S)+cslow(ix^S)*(cslow(ix^S)-wroe(ix^S,m0_+idim))))
      case(entroW_)
        rq(ix^S)= q(ix^S)*half*(^C&wroe(ix^S,m^C_)**2+)
      case(diverW_)
        if(divbwave)then
           rq(ix^S)= q(ix^S)*wroe(ix^S,b0_+idim)
        else
           rq(ix^S)= zero
        endif
      case(alfvRW_)
        rq(ix^S)=+q(ix^S)*bv(ix^S)
      case(alfvLW_)
        rq(ix^S)=-q(ix^S)*bv(ix^S)
    end select
  case(m^C_)
    if(iw==m0_+idim)then
      select case(il)
        case(fastRW_)
          rq(ix^S)=q(ix^S)*afast(ix^S)*(wroe(ix^S,iw)+cfast(ix^S))
        case(fastLW_)
          rq(ix^S)=q(ix^S)*afast(ix^S)*(wroe(ix^S,iw)-cfast(ix^S))
        case(slowRW_)
          rq(ix^S)=q(ix^S)*aslow(ix^S)*(wroe(ix^S,iw)+cslow(ix^S))
        case(slowLW_)
          rq(ix^S)=q(ix^S)*aslow(ix^S)*(wroe(ix^S,iw)-cslow(ix^S))
        case(entroW_)
          rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
        case(diverW_,alfvLW_,alfvRW_)
          rq(ix^S)=zero
      end select
    else
      select case(il)
        case(fastRW_)
          rq(ix^S)=q(ix^S)*(afast(ix^S)*wroe(ix^S,iw)-aslow(ix^S)*&
            cslow(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
        case(fastLW_)
          rq(ix^S)=q(ix^S)*(afast(ix^S)*wroe(ix^S,iw)+aslow(ix^S)*&
            cslow(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
        case(slowRW_)
          rq(ix^S)=q(ix^S)*(aslow(ix^S)*wroe(ix^S,iw)+afast(ix^S)*&
            cfast(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
        case(slowLW_)
          rq(ix^S)=q(ix^S)*(aslow(ix^S)*wroe(ix^S,iw)-afast(ix^S)*&
            cfast(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
        case(entroW_)
          rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
        case(diverW_)
          rq(ix^S)=zero
        case(alfvRW_)
          if(iw==m0_+idir)then
            rq(ix^S)=+q(ix^S)*wroe(ix^S,b0_+jdir)
          else
            rq(ix^S)=-q(ix^S)*wroe(ix^S,b0_+idir)
          endif
        case(alfvLW_)
          if(iw==m0_+idir)then
            rq(ix^S)=-q(ix^S)*wroe(ix^S,b0_+jdir)
          else
            rq(ix^S)=+q(ix^S)*wroe(ix^S,b0_+idir)
          endif        
      end select
    end if ! iw=m_idir,m_jdir
  case(b^C_)
    if(iw==b0_+idim)then
      if(il==diverW_ .and. divbwave)then
        rq(ix^S)=q(ix^S)
      else
        rq(ix^S)=zero
      endif
    else
      select case(il)
        case(fastRW_,fastLW_)
          rq(ix^S)=+q(ix^S)*aslow(ix^S)*sqrt(csound2(ix^S))*wroe(ix^S,iw)&
                    /wroe(ix^S,rho_)
        case(slowRW_,slowLW_)
          rq(ix^S)=-q(ix^S)*afast(ix^S)*sqrt(csound2(ix^S))*wroe(ix^S,iw)&
                    /wroe(ix^S,rho_)
        case(entroW_,diverW_)
          rq(ix^S)=zero
        case(alfvRW_,alfvLW_)
          if(iw==b0_+idir)then
             rq(ix^S)=-q(ix^S)*wroe(ix^S,b0_+jdir)&
                  /sign(wroe(ix^S,rho_),wroe(ix^S,b0_+idim))
          else
             rq(ix^S)=+q(ix^S)*wroe(ix^S,b0_+idir)&
                  /sign(wroe(ix^S,rho_),wroe(ix^S,b0_+idim))
          end if
      end select
    end if ! iw=b_idir,b_jdir
end select

return
end
!=============================================================================
! end module vacphys.mhdroe
!##############################################################################
