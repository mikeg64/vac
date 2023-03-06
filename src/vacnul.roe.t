! Subroutines related to Roe-type Riemann solvers

!=============================================================================
subroutine average(wL,wR,ix^L,iws,idim,wroe)

! Calculate some symmetric average of wL and wR and some useful quantities 
! for geteigenjump and rtimes

include 'vacdef.f'

integer:: ix^L,iws(niw_),idim
double precision, dimension(ixG^T,nw):: wL,wR,wroe
!-----------------------------------------------------------------------------

call die('vacnul.roe:Average is Undefined')
end

!=============================================================================
subroutine geteigenjump(wL,wR,wroe,ix^L,il,idim,smalla,a,jump)

! Calculate the il-th characteristic speed "a" and the "jump" in the il-th
! characteristic variable in the idim direction within ixL.
! Return the smalla array for the entropyfix based on "typeentropy(il)" switch.
! The "a" eigenvalues and the l=r**(-1) matrix is calculated from wroe.
! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw))

include 'vacdef.f'

integer:: ix^L,il,idim
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T)   :: smalla,a,jump
!-----------------------------------------------------------------------------

call die('vacnul.roe:GetEigenJump is Undefined')
end

!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'vacdef.f'

integer::          ix^L,iw,il,idim
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq
!-----------------------------------------------------------------------------

call die('vacnul.roe:Rtimes is Undefined')
end
!=============================================================================
