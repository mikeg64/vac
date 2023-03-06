!##############################################################################
! module vacphys0 - hdadiab

!=============================================================================
subroutine conserve(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Transform primitive variables into conservative ones

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'conserv')>=1
if(oktest)write(*,*)'Conserve w:',w(ixtest1,ixtest2,iwtest)

! Convert velocity to momentum
w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)=w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   *w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)
w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)=w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   *w(ixmin1:ixmax1,ixmin2:ixmax2,v2_);

if(oktest)write(*,*)'New w:',w(ixtest1,ixtest2,iwtest)

return
end

!=============================================================================
subroutine primitive(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Transform conservative variables into primitive ones

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'primitive')>=1
if(oktest)write(*,*)'Primitive w:',w(ixtest1,ixtest2,iwtest)

! Convert momentum to velocity
w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)=w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)&
   /w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)=w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)&
   /w(ixmin1:ixmax1,ixmin2:ixmax2,rho_);

if(oktest)write(*,*)'New w:',w(ixtest1,ixtest2,iwtest)

return
end

!=============================================================================
subroutine getv(w,ixmin1,ixmin2,ixmax1,ixmax2,idim,v)

! Calculate v_idim=m_idim/rho within ix^L

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),v(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

v(ixmin1:ixmax1,ixmin2:ixmax2)=w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
   +idim)/w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)

return
end

!=============================================================================
subroutine getcmax(new_cmax,w,ixmin1,ixmin2,ixmax1,ixmax2,idim,cmax)

! Calculate cmax_idim=csound+abs(v_idim) within ix^L

include 'vacdef.f'

logical:: new_cmax
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),cmax(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),csound(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
save csound
!-----------------------------------------------------------------------------

!Direction independent part of getcmax:
if(new_cmax)then
   new_cmax=.false.
   csound(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(eqpar(adiab_)*eqpar(gamma_)&
      *w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**(eqpar(gamma_)-1))
endif

cmax(ixmin1:ixmax1,ixmin2:ixmax2)=csound(ixmin1:ixmax1,ixmin2:ixmax2)&
   +abs(w(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/w(ixmin1:ixmax1,ixmin2:ixmax2,&
   rho_))

return 
end

!=============================================================================
! end module vacphys0 - hdadiab
!##############################################################################
