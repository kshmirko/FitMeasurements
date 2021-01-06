MODULE mo_DLS
! 	use iso_c_binding
	use mo_par_DLS, only  : KNpar, KRpar, KMpar, KR1par, KMD 
	implicit none
	integer :: 	key, key_RD, keyEL, keySUB, keyLS, key_org, key_fx, key_grid1, &
						 	&key_RD1, KN, KM, KR, NRATN, NDP          
	real		::	WL, RN, RK, pomin, pomax, xext, xabs, xsca, albedo
	real*4, dimension (KRpar)	::	R

	real, dimension	(KNpar) :: grid,SD
	real, dimension (KRpar) :: RD
	real, dimension (KMpar) :: f11, f12, f22, f33, f34, f44
	real, dimension (KMpar) :: ANGLE
	real	::	XBLR,XLDR

	character (len=255)                    ::	distname_O, distname_F, distname_N   
	character (len=255), dimension(KR1par) ::	comm_name   
	integer																 ::	key_SD, ID, NMD, NSD
	! number of aerosol component, need development for KSD>1
  integer,parameter :: KSD=1 
	real, dimension(KMD)   :: CM, SM, RMM     
  real, dimension(KNpar) :: RRR, AR, xgrid
  real									 :: AC
END MODULE mo_DLS

      