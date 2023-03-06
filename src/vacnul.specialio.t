!=============================================================================
subroutine readfileini_special(w)

! Reads from unitini,filenameini in user-defined format.
! Check readfileini_asc and readfileini_bin in vacio.t on what should be done.

include 'vacdef.f'

double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

call die('Special readfileini is not defined')
end
!=============================================================================
subroutine savefileout_special(qunit,w,ix^L)

! Save current results into filenameout in user-defined format.
! Check savefileout_asc and savefileout_bin in vacio.t on what should be done.

include 'vacdef.f'

integer:: qunit,ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

call die('Special savefileout is not defined')
end
!=============================================================================
subroutine savefilelog_special(qunit,w,ix^L)

! Save user-defined log data into filename(filelog_) in user-defined format.
! Check savefilelog_default on opening the file etc.

include 'vacdef.f'

integer:: qunit,ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

call die('Special savefilelog is not defined')
end
!=============================================================================
