!##############################################################################
! include vacpar - hdadiab

! For adiabatic HD: density,momentum + adiabatic_index, adiabatic_constant
!                                      p=adiab*rho**gamma

CHARACTER*7,PARAMETER::  typephys='hdadiab'       ! VACPHYS module name
CHARACTER*11,PARAMETER:: eqparname='gamma adiab'  ! Equation parameter names

INTEGER,PARAMETER:: rho_=1,m0_=rho_,m1_=m0_+1,m2_=m0_+2,nw=m2_ !flow variables
INTEGER,PARAMETER:: e_=0, ee_=1                           ! No energy
INTEGER,PARAMETER:: b0_=0                                 ! No magnetic field

INTEGER,PARAMETER:: v0_=m0_, v1_=m1_,v2_=m2_, p_=e_, pp_=ee_ !Primitive variables

INTEGER,PARAMETER:: soundRW_=1,soundLW_=2,shearW0_=2        ! Character. waves
INTEGER,PARAMETER:: extremeRW_=soundRW_,extremeLW_=soundLW_ !Potential extrema

INTEGER,PARAMETER:: mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_

INTEGER,PARAMETER:: nvector=1                             ! No. vector vars.

INTEGER,PARAMETER:: gamma_=1,adiab_=2,neqpar=2            ! equation params.
INTEGER,PARAMETER:: nprocpar=1                            ! processing params.

! end include vacpar - hdadiab
!##############################################################################
