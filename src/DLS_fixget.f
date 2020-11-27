!      SUBROUTINE MATRIX_FIX(key,key_RD,keyEL,keySUB
!     &                     ,key_org,key_fx,key_grid1
!     &                     ,KR,R,RD,KN1,grid1,KM,ANGLE
!     &                     ,WAVEL,KRE,KIM,ARE,AIM,pomin,pomax
!     &                     ,NRATN,RATIO
!     &                     ,distname_O,distname_F,distname_N, NDP)
      SUBROUTINE MATRIX_FIX(KN1,grid1
     &                     ,WAVEL,KRE,KIM,ARE,AIM
     &                     ,RATIO
     &                     ,NDP)
!c ** matrix_fixget_S.f previous file name 30/06/2011
!c ** ANGLE SPLINE interpolation version
!c ** rootdir is defined in "mo_par_DLS.f90"
!c ** 12/04/03 f22 interpolation is logarithmic
!c ** 05/05/03 this version can be used to retrieve an aspect ratio
!c ** 13/10/03 IF(ANGLE<=40.)ln(f33)-interpol.
!c **          IF(ANGLE<=50.)ln(f44)-interpol.
!c ** 16/07/18 IF(ANGLE<=20.)ln(f33)-interpol.
!c **          IF(ANGLE<=20.)ln(f44)-interpol.

!c **************************************************************** c
!c **   Subroutine gets original and calculates fixed kernel     ** c
!c **   matrices for given                                       ** c 
!c **   aspect ratio distribution and scattering angles          ** c
!c **                                                            ** c
!c **************************************************************** c
!c **                                                            ** c
!c ** INPUT:                                                     ** c
!c **                                                            ** c
!c **   key  = 1 - create fixed kernels (for fixed axis          ** c 
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **          2 - read fixed kernels from 'Rke...fix...'  files ** c
!c **          3 - create fixed kernels but don't save them      ** c
!c **          4 - don't create fixed kernels, calculate         ** c
!c **              opt.characteristics from original kernels     ** c
!c **   key_RD =1 - volume mixture of spheroids                  ** c
!c **           2 - surface area  mixture of spheroids           ** c
!c **   keyEL=  1 - calculate F11                                ** c
!c **           2 - F11,F12                                      ** c
!c **           3 - F11,F12,F22                                  ** c
!c **           4 - F11,F12,F22,F33                              ** c
!c **           5 - F11,F12,F22,F33,F34                          ** c
!c **           6 - F11,F12,F22,F33,F34,F44                      ** c
!c **   key_org=0 - read original kernels simultaneously make    ** c 
!c **               angle interpolation and calculate opt.char.  ** c
!c **           1 -  -"-, save new kernels in                    ** c
!c **                /NEWORG/ directory, STOP                    ** c
!c **   key_fx works when key=1                                  ** c
!c **        =0 -  create fixed kernels (for fixed axis          ** c 
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **         1 - save fixed kernels with original kernel format ** c
!c **             in order to be used as input kernels;          ** c
!c **             'Rke...fix' kernels have to be renamed and moved* c
!c **             into directory 'dir_name'(see 'matrix_fixget.f')* c
!c **             The files can not be used if key=2.            ** c
!c **                                                            ** c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.txt'                                  ** c
!c **           1 - 'grid1_fix.txt'                              ** c
!c **   KR  - number of aspect ratios                            ** c
!c **   R(KR)  - grid ratios                                     ** c
!c **   RD(KR) - aspect ratio distribution for grid aspect ratios** c
!c **   KM   - number of scattering angles for fixed kernels     ** c
!c **   ANGLE(KM)- scattering angles for fixed kernels           ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory               ** c
!c **                                      name (key_org=1)      ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   UF...(KMpar,KN1par,KIMpar,KREpar) - kernels for          ** c
!c **           given aspect ratio distribution                  ** c
!c **                                                            ** c
!c **   UFEA(2,KN1par,KIMpar,KREpar) - extinction and absorption ** c
!c **                                                    kernels ** c
!c **              1 - extinction                                ** c
!c **              2 - absorption                                ** c
!c **   KN1 - number of grid radii from original or fixed kernels** c
!c **   grid1(KN1) - grid radii                                  ** c
!c **   WAVEL - wavelength from original or fixed kernels        ** c
!c **   KRE   - number of real parts of refr.ind.                ** c
!c **   KIM   - number of imaginary parts of refr.ind.           ** c
!c **   ARE(KRE) - real parts of refr.ind.                       ** c
!c **   AIM(KIM) - imaginary parts of refr.ind.                  ** c
!c **************************************************************** c
!c **************************************************************** c
!c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.txt'                                  ** c
!c **           1 - 'grid1_fix.txt'                              ** c
!c

      use, intrinsic :: IEEE_ARITHMETIC ! fd
      use alloc      
      use alloc1
      USE mo_par_DLS
      use mo_DLS, only : key,key_RD,keyEL,keySUB
     &                  ,key_org,key_fx,key_grid1
     &                  ,KR,R,RD,KM,ANGLE
     &                  ,pomin,pomax
     &                  ,distname_O,distname_F,distname_N
     &                  ,NRATN,comm_name
!      use mod_os
      use mo_intrpl_linear
	  use mo_intrpl_spline

      integer          ::  NEL (0:6)
      character(len=2) ::  NELC(0:6)

      real rmin1,rmax1,rgrid1,rgrid2
      real RATIO(KR1par),RRATN
      dimension grid1(KN1par) 
     &         ,ANGLE1(KM1par)
     &         ,ANGLE2(KM1par)
     &         ,RAB1(KM1par)
     &         ,RAB2(KM1par,KN1par,KIMpar,KREpar)
      dimension RDc(KR), X(2), Y(2)
     &         ,ARE(KREpar), AIM(KIMpar)
      real WAVEL
!      real LINEAR		!AH original
	  integer NDP
!cl     &, LINEAR_LN
      CHARACTER(255) name,dir_name_O,dir_name_F,dir_name_N
      CHARACTER(255) full_name
      CHARACTER(255) comm_name1	  

!c      save NDP
      INTEGER                      :: key_spln, NNEL, KKEL
      double precision, dimension(KM1par)    :: XXS2,YYS2
      double precision                       :: XARG,YFIT
      double precision, dimension(KM1par+4)  :: KS1,CS1
!c ----- Timer ------
      !real*4, external :: etime,dtime

      real ::  tarray(2),T_RFM,T_CFM,T_CFM0
      real   ::  UTEST, tiny = 1e-29
      real :: xx1, xx2
      real :: ang_f33, ang_f44

      ang_f33 = 20.0
      ang_f44 = 20.0

!cl    write(*,*)'distname_O=',TRIM(distname_O)
!cl	  write(*,*)'distname_F=',TRIM(distname_F)
!cl	  write(*,*)'distname_N=',TRIM(distname_N)

	    !dir_name_O=TRIM(rootdir)//TRIM(distname_O)//'/'
	    !dir_name_F=TRIM(rootdir)//TRIM(distname_F)//'/'
	    !dir_name_N=TRIM(rootdir)//TRIM(distname_N)//'/'

	    dir_name_O = TRIM(distname_O)//'/'
	    dir_name_F = TRIM(distname_F)//'/'
	    dir_name_N = TRIM(distname_N)//'/'
       
      PI2=2.*ACOS(-1.)
      NEL=(/0, 11, 12, 22, 33, 34, 44/)
      NELC = (/'00','11','12','22','33','34','44'/)
        T_RFM=0.
        T_CFM=0.
       write(*,*) 'key_grid1=',key_grid1
	       
      if(key_grid1.eq.0) then
       full_name=TRIM(dir_name_O)//'grid1.txt'
      else
       full_name=TRIM(dir_name_F)//'grid1_fix.txt'
      endif ! key_grid1
      write(*,*)'file=',TRIM(full_name)
      OPEN (10,file=full_name,status='unknown')
      READ (10,*) KN1, XWL

      if(KN1 .gt. KN1par) then
      write(*,*) 'in GET_MATRIX: KN1=',KN1,' KN1par=',KN1par
      stop ' !!! KN1.ne.KN1par'
      endif
      DO I=1,KN1
      READ (10,*) grid1(I)
      ENDDO ! I 
      rgrid1=grid1(1)
      rgrid2=grid1(KN1)
      READ (10,*) KM1
      if(KM1.gt.KM1par) then
      write(*,*) 'in GET_MATRIX: KM1=',KM1,' KM1par=',KM1par
      stop ' !!! KM1.ne.KM1par'
      endif
      DO J=1,KM1
      READ (10,*) ANGLE1(J)
      ENDDO ! J
      CLOSE(10)
      if((key .eq. 4) .and. (KM .ne. KM1)) then
      write(*,*) 'if key=4, in input.dat KM=',KM,
     &  ' should be equal to KM1=',KM1,' in grid1.txt'
      stop 'STOP in MATRIX_FIX (matrix_fixget.f)'      
	endif ! key
!c
!c ** Redefine grids for fixed or NEWORG kernels
!c
      NDPP=0
      if(key.eq.1.or.key_org.eq.1) then
  	rgmin=pomin*XWL/pi2
	rgmax=pomax*XWL/pi2 
!tl	write(*,*) 'pomin=',pomin,' XWL=',XWL,' pi2=',pi2,rgmin
!tl	stop
!c
      xx1 = rgmin
      xx2 = grid1(1)
      if(xx1 .lt. xx2) then
        write(*,*) 'XWL=',XWL,' IS wave length in kernels equal to XWL?'
        write(*,*) 'check input.dat: rgmin=',rgmin,' < ',
     &     'grid1(1)=',grid1(1)
        STOP 'STOP in MARTRIX_FIX (matrix_fixget.f)'
      endif
      xx1 = rgmax
      xx2 = grid1(KN1)
      if(xx1 .gt. xx2) then
        write(*,*) 'XWL=',XWL,' IS wave length in kernels equal to XWL?'
        write(*,*) 'check input.dat: rgmax=',rgmax,' > ',
     &      'grid1(KN1)=',grid1(KN1)
        STOP 'STOP in MARTRIX_FIX (matrix_fixget.f)'
      endif

      ind=0
	do i=1,KN1-1
!c	write(*,*) 'i=',i,' rgmin=',rgmin,' grid1(i+1)=',grid1(i+1)
      if(rgmin.le.grid1(i+1)) then
!c	xr1=grid1(i)
	nn1=i
	ind=ind+1
	endif
      if(ind.eq.1) EXIT
	enddo ! i
!c	 
      ind=0
	do i=KN1,2,-1
!c	write(*,*) 'i=',i,' rgmax=',rgmax,' grid1(i-1)=',grid1(i-1)
      if(rgmax.ge.grid1(i-1)) then
!c	xr2=grid1(i)
	nn2=i
	ind=ind+1
	endif
      if(ind.eq.1) EXIT
	enddo ! i
!c
      if(rgmin.eq.grid1(1)) then
!c	xr1=grid1(1)
	nn1=1
	endif ! rgmin 
      if(rgmax.eq.grid1(KN1)) then
!c	xr2=grid1(KN1)
	nn2=KN1
	endif ! rgmax
      nn3=nn2-nn1+1
!cl	write(*,*) 'nn1=',nn1,' nn2=',nn2,' nn3=',nn3
!cl	write(*,*) 'xr1=',xr1,' xr2=',xr2,' KN1=',KN1

!cl	STOP 'TEST STOP'
!c	 	  
!      write(*,*) 'before grid1_fix.txt' 
      if(key_org.eq.0) then 
       full_name=TRIM(dir_name_F)//'grid1_fix.txt'
	else
       full_name=TRIM(dir_name_N)//'grid1.txt'
	endif
	write(*,*)'file=',TRIM(full_name)
      open(77,file=full_name, status='replace')
      WRITE(77,'(i4,f8.3)') nn3,XWL
      DO I=nn1,nn2
      WRITE(77,'(e15.7)') grid1(I)
      ENDDO ! I 
      WRITE(77,*) KM
      DO J=1,KM
      WRITE(77,'(f6.2)') ANGLE(J)
      ENDDO ! J      
	close(77)
      WRITE(*,*) 
      WRITE(*,*) '  ATTENTION:'
      WRITE(*,*) '    New input file grid1.txt (key=1,key_org=1) or' 
	WRITE(*,*) '    grid1_fix.txt (key=1,key_org=0) has been created'
	WRITE(*,*) '    for your further calculations !!!'      
      WRITE(*,*) 
	endif ! key&key_org
!c
!c ** Define directories
!c
!c ** for AERONET
!cl      NNEL=2
!cl      if(key_f11.eq.0) NNEL=7
	  NNEL=keyEL+1
	  KKEL=6-keyEL
	  
      IF(key .EQ. 2) THEN
								if(keySUB .eq. 0) then
									!T_RFM=dtime(tarray) !+++
									!T_RFM=tarray(1)+tarray(2)
						        endif

        open(10,file=TRIM(TRIM(dir_name_F)//"Rkernel1_fix_00.txt"), 
     &                                                 status="old")
	  IF(keyEL.gt.0) THEN
        open(11,file=TRIM(TRIM(dir_name_F)//"Rkernel1_fix_11.txt"),
     &                                                 status="old")
      if(keyEL.gt.1)
     &  open(12,file=TRIM(dir_name_F)//"Rkernel1_fix_12.txt",
     &                                                 status="old")
      if(keyEL.gt.2)	 
     &  open(13,file=TRIM(dir_name_F)//"Rkernel1_fix_22.txt",
     &                                                 status="old")
      if(keyEL.gt.3)	 
     &  open(14,file=TRIM(dir_name_F)//"Rkernel1_fix_33.txt",
     &                                                 status="old")
	  if(keyEL.gt.4) 
     &  open(15,file=TRIM(dir_name_F)//"Rkernel1_fix_34.txt",
     &                                                 status="old")
	  if(keyEL.gt.5) 
     &  open(16,file=TRIM(dir_name_F)//"Rkernel1_fix_44.txt",
     &                                                 status="old")
      ENDIF ! keyEL
        open(20,file='CHECK.dat',       status='unknown')

        DO II=1,NNEL
        WRITE(20,*) 'INPUT CHECK: Nelement=',NEL(II-1)
        READ(10+II-1,*) key_RD1
        IF(key_RD.ne.key_RD1) then
        WRITE(20,14) II,key_RD,key_RD1
        WRITE(20,*) 'STOP: key_RD.ne.keyRD1 in Rke...fix'
        WRITE(*,14) II,key_RD,key_RD1
        WRITE(*,*) 'STOP: key_RD.ne.keyRD1 in Rke...fix'
        STOP
        ENDIF ! key_RD
        READ(10+II-1,*) rmin,rmax
          rmin1=rmin
	    rmax1=rmax
          IF(rmin1.ne.rgrid1.or.rmax1.ne.rgrid2) THEN
          WRITE(20,*) 'grid1(1)=',rgrid1,' rmin=',rmin1,
     &    ' rmax=',rmax1,' grid1(KN1)=',rgrid2
          WRITE(*,*) 'grid1(1)=',rgrid1,' rmin=',rmin1,
     &    ' rmax=',rmax1,' grid1(KN1)=',rgrid2
	    WRITE(*,*) 'Compare1: grid1(1) in grid1.dat',
     &                    ' with RMIN in Rke...fix'
	    WRITE(*,*) 'Compare2: grid1(KN1) in grid1.dat',
     &                    ' with RMAX in Rke...fix'
          STOP 'STOP: grid1(1)/grid(KN1).ne.rmin/rmax'
          ENDIF
        READ(10+II-1,*) kk
        WRITE(20,*) 'KR=',KR,' KR=',kk
        READ(10+II-1,*) 
        DO I=1,KR
        READ(10+II-1,*) xx,yy
        WRITE(20,*) R(I),RD(I),' R,RD',xx,yy,' R,RD' 
        ENDDO ! I
        READ(10+II-1,*) kk
        WRITE(20,*) 'KN1=',KN1,' KN1=',-kk
	   IF(KN1 .ne. -kk) THEN
         WRITE(*,*) 'in grid1.dat KN1=',KN1,
     &                 ' .ne. KN1=',-kk,' in Rke...fix'
         STOP 'STOP in matrix_fixget'
	   ENDIF
        ENDDO ! II

        DO II=2,NNEL
        WRITE(20,*) 'INPUT CHECK: Nelement=',NEL(II-1)
        READ(10+II-1,*) kk
        WRITE(20,*) 'KM=',KM,' KM=',kk 
        READ(10+II-1,*) ANGLE2(:KM)
        WRITE(20,*)  'SCATTERING ANGLES2:'
        WRITE(20,15) ANGLE2(:KM)
        WRITE(20,*)  'SCATTERING ANGLES:'
        WRITE(20,15) ANGLE(:KM)
        ENDDO ! II

        DO II=1,NNEL
        WRITE(20,*) 'INPUT CHECK: Nelement=',NEL(II-1)
        READ(10+II-1,*) xx,yy
        WRITE(20,*) 'ARE(1),ARE(KRE)',xx,yy
        READ(10+II-1,*) xx,yy
        WRITE(20,*) 'AIM(1),AIM(KIM)',xx,yy 
        READ(10+II-1,*) KRE,KIM
        WRITE(20,*) 'KRE,KIM',KRE,KIM
        IF(KIM.lt.0) KIM=-KIM
        ENDDO ! II
        WRITE(*,*) 'Fixed kernel matrices have been read :'
        DO IRE=1,KRE
        DO IIM=1,KIM
        DO II=1,NNEL
        if(key_fx.eq.1) then
        READ(10+II-1,*) kk,aa
	  else ! key_fx=0
        READ(10+II-1,*) kk  
	  endif ! key_fx    
        READ(10+II-1,*) WAVEL,ARE(IRE),AIM(IIM)
        IF(AIM(IIM).LT.0) AIM(IIM)=-AIM(IIM) 
!cl        WRITE(20,*) 'WAVEL,ARE,AIM,NEL',
!cl     &                WAVEL,ARE(IRE),-AIM(IIM),NEL(II-1) 
        ENDDO ! II

        READ(10,*)
        READ(10,11) UFEA(1,:KN1,IIM,IRE)
        READ(10,*)
        READ(10,11) UFEA(2,:KN1,IIM,IRE)

        if(keyEL .gt. 0) then
          DO I=1,KN1
          READ(11,11) UF11(1:KM,I,IIM,IRE)
          ENDDO ! I
        endif
        if(keyEL .gt. 1) then
          DO I=1,KN1
          READ(12,11) UF12(1:KM,I,IIM,IRE)
          ENDDO ! I
        endif
        if(keyEL .gt. 2) then
          DO I=1,KN1
          READ(13,11) UF22(1:KM,I,IIM,IRE)
          ENDDO ! I
        endif
        if(keyEL .gt. 3) then
          DO I=1,KN1
          READ(14,11) UF33(1:KM,I,IIM,IRE)
          ENDDO ! I
        endif
        if(keyEL.gt.4) then
          DO I=1,KN1
          READ(15,11) UF34(1:KM,I,IIM,IRE)
          do j=1,KM
          UTEST = UF34(j,I,IIM,IRE)
          if((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0)) 
     &    UF34(j,I,IIM,IRE) = 0.0 
          enddo ! j
          ENDDO ! I
        endif
        if(keyEL.gt.5) then
          DO I=1,KN1
          READ(16,11) UF44(1:KM,I,IIM,IRE)
          ENDDO ! I
        endif

        ENDDO ! IIM
        WRITE(*,16) WAVEL,ARE(IRE)
        ENDDO ! IRE
      close(10)
      if(keyEL .gt. 0) close(11)
      if(keyEL .gt. 1) close(12)
      if(keyEL .gt. 2) close(13)
      if(keyEL .gt. 3) close(14)
      if(keyEL .gt. 4) close(15)
      if(keyEL .gt. 5) close(16)
      close (20)

	  IF(keyEL.gt.0) THEN
      UF11(1:KM,1:KN1,1:KIM,1:KRE)=
     &     LOG( UF11(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL )
      if(keyEL.gt.1)
     &  UF12(1:KM,1:KN1,1:KIM,1:KRE)=
     &         UF12(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL 
      if(keyEL.gt.2)
     &  UF22(1:KM,1:KN1,1:KIM,1:KRE)=
     &     LOG( UF22(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL ) 
	  if(keyEL.gt.3) then
	  DO J=1,KM
		UF33(J,1:KN1,1:KIM,1:KRE)=
     &      UF33(J,1:KN1,1:KIM,1:KRE)*WAVEL
       IF(ANGLE(J) .LE. ang_f33)
     & UF33(J,1:KN1,1:KIM,1:KRE)=
     &      LOG( UF33(J,1:KN1,1:KIM,1:KRE) )
      ENDDO ! J
	  endif
      if(keyEL.gt.4)
     &  UF34(1:KM,1:KN1,1:KIM,1:KRE)=
     &         UF34(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL
	  if(keyEL.gt.5) then
      DO J=1,KM	  	  
        UF44(J,1:KN1,1:KIM,1:KRE)=
     &         UF44(J,1:KN1,1:KIM,1:KRE)*WAVEL
      IF(ANGLE(J) .LE. ang_f44)
     &     UF44(J,1:KN1,1:KIM,1:KRE)=
     &         LOG( UF44(J,1:KN1,1:KIM,1:KRE) )
      ENDDO ! J
	  endif
	  ENDIF ! keyEL

13    FORMAT(3E12.4,I4)                    
14    FORMAT('II=',i2,' key_RD=',i2,' key_RD1=',i2)
16    FORMAT(12x,'wl=',f5.2,2x,'n=',f8.5)
!cl      WRITE(*,*) 'Fixed kernel matrices have been read'
      IF(key_RD.EQ.1) WRITE(*,*) 'Volume mixture of spheroids'
      IF(key_RD.EQ.2) WRITE(*,*) 'Surface area mixture of spheroids'

                                        if(keySUB.eq.0) then
                                            !T_RFM=dtime(tarray) !+++
                                            !T_RFM=tarray(1)+tarray(2)
      !write(6,*) 
      !write(6,*) 
!     &        '------------------ T I M I N G ------------------' 
      !WRITE(*,61) T_RFM/60.     
										endif
										
   61 format('  Read fixed matr. ........ ',f8.3,' min.')

      ENDIF ! key=2

      IF(key .NE. 2) THEN
      IF(key .EQ. 1) THEN
!c
!c *** ALLOCATE and INITIALIZE ARRAYS
!c
        ALLOCATE(UEA(KR1par,2,KN1par,KIMpar,KREpar),    stat=ierr)
        if(ierr/=0) stop 'Can not allocate UEA array'
        UEA=0.
        IF(keyEL.gt.0) THEN
          ALLOCATE(U11(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U11 array'
          U11=0. 
        if(keyEL.gt.1) then
          ALLOCATE(U12(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U12 array'
          U12=0.
        endif
        if(keyEL.gt.2) then		
          ALLOCATE(U22(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U22 array'
          U22=0.
        endif
        if(keyEL.gt.3) then		
          ALLOCATE(U33(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U33 array'
          U33=0.
        endif
        if(keyEL.gt.4) then				  
          ALLOCATE(U34(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U34 array'
          U34=0.
        endif
        if(keyEL.gt.5) then		
          ALLOCATE(U44(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U44 array'
          U44=0.
        endif
      ENDIF ! keyEL=0
      ENDIF ! key=1
      IF(NDP.EQ.0) THEN
      T_CFM=0.
!c
!c ** READ ORIGINAL kernels 
!c
!c      OPEN (10,file=TRIM(rootdir)//'name.dat',status='old')
      !OPEN (10,file='name.dat',status='old')
      !READ(10,*) NRATN	  
      IF(NRATN.gt.KR1par) 
     &    STOP ' in GET_MATRIX 1: NRATN.gt.KR1par !!!'

      DO IRATN=1,NRATN
!c **
!c ** READ U11, U12, U22, U33, U34, U44 matrices
!c **
       DO KEL=1,NNEL-1
        !READ(10,*) name
!cl        OPEN(11,FILE=name,status='old')    
		   
        comm_name1=trim(comm_name(IRATN))//'_'//NELC(KEL)//'.txt'
        full_name=TRIM(dir_name_O)//trim(comm_name1)
        write(*,*) trim(full_name)
        OPEN(11,FILE=trim(full_name),status='old')
      if(key_org.eq.1) then
        full_name=TRIM(dir_name_N)//trim(comm_name1)
!cl        full_name=TRIM(name)
        OPEN(21,FILE=trim(full_name),status='unknown')
        write(*,*) trim(full_name)
      endif ! key_org 
        READ(11,*) rmin,rmax,RATIO(iratn)
        READ(11,*) KNN
        IF(KNN.LT.0) KNN=-KNN
        READ(11,*) KMM
        IF(KNN.ne.KN1) STOP ' in GET_MATRIX 1: KNN.ne.KN1 !!!'
        IF(KMM.ne.KM1) STOP ' in GET_MATRIX 1: KMM.ne.KM1 !!!'
        READ(11,*) ANGLE2(1:KMM)
!cl        WRITE(*,*)  'ANGLE2 in GET_MATRIX NEL=',NEL(KEL)
!cl        WRITE(*,11) ANGLE2(:KM1)
!cl        WRITE(*,*)  'ANGLE1  in GET_MATRIX NEL=',NEL(KEL)
!cl        WRITE(*,11) ANGLE1(:KM1)

        WRITE(*,*) 'READ matrix U',NEL(KEL)
        READ(11,*) RREMIN, RREMAX
!cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'   
        READ(11,*) RIMMIN, RIMMAX
!cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'        
        READ(11,*) KRE, KIM
!cl        WRITE(*,*) KRE, KIM,' KRE, KIM'
      if(key_org.eq.1) then
        WRITE(21,*) grid1(nn1),grid1(nn2),RATIO(iratn),
     &                      '  rmin, rmax, RATIO'
        WRITE(21,*) -nn3,' number of intervals'
        WRITE(21,*) KM,'  number of angles'
        WRITE(21,27) ANGLE(1:KM)

        WRITE(*,*) 'WRITE matrix U',NEL(KEL)
        WRITE(21,*) RREMIN, RREMAX,' real refr. indices'
!cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'   
        WRITE(21,*) RIMMIN, RIMMAX,' imag refr. indices'
!cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'        
        WRITE(21,*) KRE, KIM,' number of intervals for opt. const'
!cl        WRITE(*,*) KRE, KIM,' KRE, KIM'
       endif ! key_org=1

        IF(KRE.LT.0) KRE=-KRE
        IF(KIM.LT.0) KIM=-KIM 
       DO  IRE=1,KRE
       DO  IIM=1,KIM
        READ(11,*) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
        READ(11,*) WAVEL,ARE(IRE),AIM(IIM) 
        AIM(IIM)=-AIM(IIM)

      if(key_org.eq.1) then
      if(NDPP.eq.0) then
        if(XWL.ne.WAVEL) then
        write(*,*) 'in grid1.dat: XWL=',XWL,'.NE.',
     &                         ' in Rke...: WAVEL=',WAVEL
      STOP 'STOP in MATRIX_FIX (matrix_fixget.f)'
      endif ! XWL
        NDPP=1
      endif ! NDPP
        WRITE(21,21) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
        WRITE(21,20) WAVEL,ARE(IRE),-AIM(IIM)
      endif ! key_org=1 

! key .lt. 4
      IF(key .lt. 4) THEN
        IF(KEL .EQ. 1) then
          DO  I=1,KN1
          READ(11,*) RAB1(1:KM1)
          RAB2(1:KM1,I,IIM,IRE)=RAB1(1:KM1)
          XXS2(1:KM1)=ANGLE1(1:KM1)
          YYS2(1:KM1)=LOG(RAB1(1:KM1))
          key_spln=0
          DO J=1,KM
          XARG=ANGLE(J)
          CALL intrpl_spline(KM1,XXS2(1:KM1),YYS2(1:KM1)
     &          ,XARG,YFIT,key_spln,KS1(1:KM1+4),CS1(1:KM1+4))
          U11(IRATN,J,I,IIM,IRE)=EXP(YFIT)
          key_spln=1
          ENDDO ! J
          ENDDO ! I
          if(key_org.eq.1) then
            DO  I=nn1,nn2
            WRITE(21,27) U11(IRATN,1:KM,I,IIM,IRE)
            ENDDO ! I
          endif ! key_org=1
        ELSE IF(KEL .EQ. 2) THEN
          DO  I=1,KN1
          READ(11,*) RAB1(1:KM1)
           if(ANGLE1(1).eq.0.)     RAB1(1)  = 0.
           if(ANGLE1(KM1).eq.180.) RAB1(KM1)= 0.
          XXS2(1:KM1)=ANGLE1(1:KM1)
          YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
          key_spln=0
          DO J=1,KM
          XARG=ANGLE(J)
          CALL intrpl_spline(KM1,XXS2(1:KM1),YYS2(1:KM1)
     &          ,XARG,YFIT,key_spln,KS1(1:KM1+4),CS1(1:KM1+4))
           if(J.eq.1.and.ANGLE1(J).eq.0.)       YFIT = 0.
           if(J.eq.KM1.and.ANGLE1(KM1).eq.180.) YFIT = 0.
          U12(IRATN,J,I,IIM,IRE)=YFIT*U11(IRATN,J,I,IIM,IRE)
          key_spln=1
          ENDDO ! J
          ENDDO ! I
          if(key_org.eq.1) then
            DO  I=nn1,nn2
            WRITE(21,27) U12(IRATN,1:KM,I,IIM,IRE)
            ENDDO ! I
          endif ! key_org=1
        ELSE IF(KEL .EQ. 3) THEN
          DO  I=1,KN1
          READ(11,*) RAB1(1:KM1)
          XXS2(1:KM1)=ANGLE1(1:KM1)
          YYS2(1:KM1)=LOG(RAB1(1:KM1))
          key_spln=0
          DO J=1,KM
          XARG=ANGLE(J)
          CALL intrpl_spline(KM1,XXS2(1:KM1),YYS2(1:KM1)
     &          ,XARG,YFIT,key_spln,KS1(1:KM1+4),CS1(1:KM1+4))
          U22(IRATN,J,I,IIM,IRE)=EXP(YFIT)
          key_spln=1
          ENDDO ! J
          ENDDO ! I
          if(key_org.eq.1) then
            DO  I=nn1,nn2
            WRITE(21,27) U22(IRATN,1:KM,I,IIM,IRE)
            ENDDO ! I
          endif ! key_org=1
        ELSE IF(KEL .EQ. 4) THEN
          DO  I=1,KN1
          READ(11,*) RAB1(1:KM1)
          XXS2(1:KM1)=ANGLE1(1:KM1)
          YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
          key_spln=0
          DO J=1,KM
          XARG=ANGLE(J)
          CALL intrpl_spline(KM1,XXS2(1:KM1),YYS2(1:KM1)
     &          ,XARG,YFIT,key_spln,KS1(1:KM1+4),CS1(1:KM1+4))
          U33(IRATN,J,I,IIM,IRE)=YFIT*U11(IRATN,J,I,IIM,IRE)
          key_spln=1
          ENDDO ! J
          ENDDO ! I
          if(key_org.eq.1) then
            DO  I=nn1,nn2
            WRITE(21,27) U33(IRATN,1:KM,I,IIM,IRE)
            ENDDO ! I
          endif ! key_org=1
        ELSE IF(KEL .EQ. 5) THEN
          DO  I=1,KN1
          READ(11,*) RAB1(1:KM1) 
           do j=1,KM1
             UTEST = RAB1(j)
             if((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0))
     &       RAB1(j) = 0.0
           enddo ! j
           if(ANGLE1(1).eq.0.)     RAB1(1)  = 0.
           if(ANGLE1(KM1).eq.180.) RAB1(KM1)= 0.
           XXS2(1:KM1)=ANGLE1(1:KM1)
           YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
           key_spln=0
	         DO J=1,KM
           XARG=ANGLE(J)
           CALL intrpl_spline(KM1,XXS2(1:KM1),YYS2(1:KM1)
     &          ,XARG,YFIT,key_spln,KS1(1:KM1+4),CS1(1:KM1+4))
           if(J.eq.1.and.ANGLE1(J).eq.0.)     YFIT = 0.
           if(J.eq.KM1.and.ANGLE1(J).eq.180.) YFIT = 0.
           if(J .eq. -180) then
           write(*,*) 'U34,YFIT,U11 - ',  
     &     U34(IRATN,J,I,IIM,IRE),YFIT,U11(IRATN,J,I,IIM,IRE)
           endif                     
           if((abs(YFIT) .lt. tiny) .and. (abs(YFIT) .gt. 0.0)) 
     &        YFIT = 0.0
           U34(IRATN,J,I,IIM,IRE)=YFIT*U11(IRATN,J,I,IIM,IRE)
           if((abs(U34(IRATN,J,I,IIM,IRE)) .lt. tiny) .and.
     &        (abs(U34(IRATN,J,I,IIM,IRE)) .gt. 0.0))
     &        U34(IRATN,J,I,IIM,IRE) = 0.0
           key_spln=1
          ENDDO ! J
          ENDDO ! I
          if(key_org.eq.1) then
            DO  I=nn1,nn2
            WRITE(21,27) U34(IRATN,1:KM,I,IIM,IRE)
            ENDDO ! I
          endif ! key_org=1
        ELSE IF(KEL .EQ. 6) THEN
          DO  I=1,KN1
          READ(11,*) RAB1(1:KM1)
          XXS2(1:KM1)=ANGLE1(1:KM1)
          YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
          key_spln=0
          DO J=1,KM
          XARG=ANGLE(J)
          CALL intrpl_spline(KM1,XXS2(1:KM1),YYS2(1:KM1)
     &          ,XARG,YFIT,key_spln,KS1(1:KM1+4),CS1(1:KM1+4))
          U44(IRATN,J,I,IIM,IRE)=YFIT*U11(IRATN,J,I,IIM,IRE)
          key_spln=1
          ENDDO ! J
          ENDDO ! I
          if(key_org.eq.1) then
            DO  I=nn1,nn2
            WRITE(21,27) U44(IRATN,1:KM,I,IIM,IRE)
            ENDDO ! I
          endif ! key_org=1
        ENDIF ! KEL
! key .eq. 4
      ELSE ! key=4
        IF(KEL .EQ. 1) then
          DO  I=1,KN1
          READ(11,*) U11(IRATN,1:KM1,I,IIM,IRE)
          ENDDO ! I
        ELSE IF(KEL .EQ. 2) THEN
          DO  I=1,KN1
          READ(11,*) U12(IRATN,1:KM1,I,IIM,IRE)
          ENDDO ! I
          if(ANGLE1(1)  .eq.  0.) U12(IRATN,1,  1:KN1,IIM,IRE) = 0.
          if(ANGLE1(KM1).eq.180.) U12(IRATN,KM1,1:KN1,IIM,IRE) = 0.
        ELSE IF(KEL .EQ. 3) THEN
          DO  I=1,KN1
          READ(11,*) U22(IRATN,1:KM1,I,IIM,IRE)
          ENDDO ! I
        ELSE IF(KEL .EQ. 4) THEN
          DO  I=1,KN1
          READ(11,*) U33(IRATN,1:KM1,I,IIM,IRE)
          ENDDO ! I
        ELSE IF(KEL .EQ. 5) THEN
          DO  I=1,KN1
          READ(11,*) U34(IRATN,1:KM1,I,IIM,IRE)
            do j=1,KM1
              UTEST = U34(IRATN,j,I,IIM,IRE)
              if((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0)) 
     &        U34(IRATN,j,I,IIM,IRE) = 0.0 
            enddo ! j    
          ENDDO ! I
          if(ANGLE1(1)  .eq.  0.) U34(IRATN,1,  1:KN1,IIM,IRE) = 0.
          if(ANGLE1(KM1).eq.180.) U34(IRATN,KM1,1:KN1,IIM,IRE) = 0.
        ELSE IF(KEL .EQ. 6) THEN
          DO  I=1,KN1
          READ(11,*) U44(IRATN,1:KM1,I,IIM,IRE)
          ENDDO ! I
        ENDIF ! KEL
      ENDIF ! key
        ENDDO ! IIM
        ENDDO ! IRE
      CLOSE (11)
      if(key_org.eq.1) CLOSE (21)
        ENDDO ! KEL	  	
!c **
!c ** READ UEA MATRIX
!c **
      !if(keyEL.lt.6) then
      !do KEL=1,KKEL
      !READ(10,*) name
      !enddo ! KEL
      !endif ! key_EL
      !READ(10,*) name

        comm_name1=trim(comm_name(IRATN))//'_'//NELC(0)//'.txt'
	    full_name=TRIM(dir_name_O)//trim(comm_name1)
      write(*,*) trim(full_name)
      OPEN(11,FILE=trim(full_name),status='old')

      if(key_org.eq.1) then
        full_name=TRIM(dir_name_N)//trim(comm_name1)
        OPEN(21,FILE=trim(full_name),status='unknown')
        write(*,*) trim(full_name)
      endif ! key_org 

      READ(11,*) rmin,rmax,RATIO(iratn)
!cl      WRITE(*,*) rmin,rmax,RATIO(iratn),' rmin,rmax,ratio'
      READ(11,*) KNN
      IF(KNN.LT.0) KNN=-KNN
      IF(KNN.ne.KN1) THEN
      WRITE(*,*) KNN,KN1,' KNN,KN1'
      STOP ' in GET_MATRIX 2: KNN.ne.KN1 !!!'
      ENDIF
      WRITE(*,*) 'READ matrix U',NEL(0)
      READ(11,*) RREMIN, RREMAX
!cl      WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'  
      READ(11,*) RIMMIN, RIMMAX
!cl      WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
      READ(11,*) KRE, KIM
      if(key_org.eq.1) then
        WRITE(21,*) grid1(nn1),grid1(nn2),RATIO(iratn),
     &                           '  rmin, rmax, RATIO'
        WRITE(21,*) -nn3,' number of intervals'
        WRITE(*,*) 'WRITE matrix U',NEL(KEL)
        WRITE(21,*) RREMIN, RREMAX,' real refr. indices'
!cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'   
        WRITE(21,*) RIMMIN, RIMMAX,' imag refr. indices'
!cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'        
        WRITE(21,*) KRE, KIM,' number of intervals for opt. const'
!cl        WRITE(*,*) KRE, KIM,' KRE, KIM'
       endif ! key_org=1

      IF(KRE.LT.0) KRE=-KRE
      IF(KIM.LT.0) KIM=-KIM 
      DO IRE=1,KRE
      DO IIM=1,KIM
      READ(11,*) nn,xx
!cl      WRITE(*,*) nn,xx,' KEL,RATIO'
      READ(11,*) WAVEL1,ARE(IRE),AIM(IIM)
      IF(keyEL.EQ.0) THEN
      WAVEL=WAVEL1
      ELSE
      IF(WAVEL.ne.WAVEL1) THEN
      WRITE(*,*) WAVEL, WAVEL1,' WAVEL,WAVEL1'
      STOP 'in GET_MATRIX: WAVEL.ne.WAVEL1'
      ENDIF 
      ENDIF ! keyEL=0
      AIM(IIM)=-AIM(IIM)             
      if(key_org.eq.1) then
        WRITE(21,21) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
        WRITE(21,20) WAVEL,ARE(IRE),-AIM(IIM)
      endif ! key_org=1 

      READ(11,*)           
      READ(11,*) (UEA(IRATN,1,I,IIM,IRE),I=1,KN1)
	    if(key_org.eq.1) then     
          WRITE(21,*) ' EXTINCTION (1/km, for d()/dlnr m3/m3*km):'
          WRITE(21,27) (UEA(IRATN,1,I,IIM,IRE),I=nn1,nn2)
	    endif ! key_org=1
      READ(11,*)           
      READ(11,*) (UEA(IRATN,2,I,IIM,IRE),I=1,KN1)
	    if(key_org.eq.1) then     
          WRITE(21,*) ' ABSORPTION (1/km, for dv/dlnr m3/m3*km):'
          WRITE(21,27) (UEA(IRATN,2,I,IIM,IRE),I=nn1,nn2)
	    endif ! key_org=1
      ENDDO ! IIM
      ENDDO ! IRE
      CLOSE (11)
!cl	endif ! key_org=0
      ENDDO ! IRATN
      IF(NRATN.lt.KR.and.(R(1).lt.RATIO(1).or.R(KR).lt.RATIO(NRATN))) 
     &    STOP ' in GET_MATRIX: NRATN.lt.KR' 
	 !check input.dat and name.dat)!!!
      !CLOSE(10)

   27 FORMAT (7E16.7)
   20 FORMAT (F14.5,2E16.7,'  wavel,rreal,rimag')
   21 FORMAT (I3,F14.5,'  element,ratn')

      if(key_org.eq.1) 
     & STOP 'STOP: key_org=1, new kernels have been saved'

      NDP=1

      ENDIF ! NDP

!cl ** TEST for R=1 (sphere) write kernel elements

!cl      write(*,*) 'EXT, ABS : n=',ARE(13),' k=',AIM(8)
!cl      do i=17,38
!cl      write(*,'(3e17.5)') grid1(i),UEA(1,1,I,8,13),UEA(1,2,I,8,13)
!cl	enddo 

      IF(key .EQ. 4) then
      goto 555
      write(iu_output,'(a,i0)') 'test read kernels, key = ',key
      write(iu_output,'(es12.4,a,es12.4)')
     &              minval(UEA(1,1,1:KN1,1:KIM,1:KRE)),
     &' <= UEXT <=',maxval(UEA(1,1,1:KN1,1:KIM,1:KRE))
      write(iu_output,'(es12.4,a,es12.4)')
     &              minval(UEA(1,2,1:KN1,1:KIM,1:KRE)),
     &' <= UABS <=',maxval(UEA(1,2,1:KN1,1:KIM,1:KRE))

      if(keyEL .gt. 0)
     &write(iu_output,'(es12.4,a,es12.4)')
     &             minval(U11(1,1:KM,1:KN1,1:KIM,1:KRE)),
     &' <= U11 <=',maxval(U11(1,1:KM,1:KN1,1:KIM,1:KRE))
      if(keyEL .gt. 1)
     &write(iu_output,'(es12.4,a,es12.4)')
     &             minval(U12(1,1:KM,1:KN1,1:KIM,1:KRE)),
     &' <= U12 <=',maxval(U12(1,1:KM,1:KN1,1:KIM,1:KRE))
      if(keyEL .gt. 2)
     &write(iu_output,'(es12.4,a,es12.4)')
     &             minval(U22(1,1:KM,1:KN1,1:KIM,1:KRE)),
     &' <= U22 <=',maxval(U22(1,1:KM,1:KN1,1:KIM,1:KRE))
      if(keyEL .gt. 3)
     &write(iu_output,'(es12.4,a,es12.4)')
     &             minval(U33(1,1:KM,1:KN1,1:KIM,1:KRE)),
     &' <= U33 <=',maxval(U33(1,1:KM,1:KN1,1:KIM,1:KRE))
      if(keyEL .gt. 4)
     &write(iu_output,'(es12.4,a,es12.4)')
     &             minval(U34(1,1:KM,1:KN1,1:KIM,1:KRE)),
     &' <= U34 <=',maxval(U34(1,1:KM,1:KN1,1:KIM,1:KRE))
      if(keyEL .gt. 5)
     &write(iu_output,'(es12.4,a,es12.4)')
     &             minval(U44(1,1:KM,1:KN1,1:KIM,1:KRE)),
     &' <= U44 <=',maxval(U44(1,1:KM,1:KN1,1:KIM,1:KRE))
555   continue

!c ** LOG(U...)
!c
        if(keyEL .gt. 0)
     &     U11(1:KR,1:KM,1:KN1,1:KIM,1:KRE)=
     &     LOG( U11(1:KR,1:KM,1:KN1,1:KIM,1:KRE)*WAVEL )
        !write(*,*) 'after u11'
        if(keyEL .gt. 1)
     &     U12(1:KR,1:KM,1:KN1,1:KIM,1:KRE)=
     &     U12(1:KR,1:KM,1:KN1,1:KIM,1:KRE)*WAVEL 
        !write(*,*) 'i am here after u12'
        if(keyEL .gt. 2)
     &         U22(1:KR,1:KM,1:KN1,1:KIM,1:KRE)=
     &         LOG( U22(1:KR,1:KM,1:KN1,1:KIM,1:KRE)*WAVEL )
        !write(*,*) 'after u22'
        if(keyEL .gt. 3) then
        DO J=1,KM
        U33(1:KR,J,1:KN1,1:KIM,1:KRE)=
     &        U33(1:KR,J,1:KN1,1:KIM,1:KRE)*WAVEL
        IF(ANGLE(J) .LE. ang_f33) THEN
!         do l=1,kr
!         do n=1,kre
!         do k=1,kim
!         do i=1,kn1
!         if(U33(l,J,i,k,n) .le. 0.0)
!     &     write(*,*) l,n,k,i,J,ANGLE(J),U33(l,J,i,k,n)/WAVEL,
!     &     '  l,n,k,i,J,ANGLE(J),U33(l,J,i,k,n)'
!         enddo
!         enddo
!         enddo
!         enddo

        U33(1:KR,J,1:KN1,1:KIM,1:KRE)=
     &  LOG( U33(1:KR,J,1:KN1,1:KIM,1:KRE) )
        ENDIF
        ENDDO ! J
      do l=1,KR
      do j=1,KM
      do i=1,KN1
      do k=1,KIM
      do n=1,KRE
!fd
      if (.not. ieee_is_normal(U33(l,j,i,k,n))) then
        write(*,'(a,5i5,a,e16.6)') 'l,j,i,k,n: ',l,j,i,k,n,
     &  '  U33=',U33(l,j,i,k,n)
      end if
      enddo
      enddo
      enddo
      enddo
      enddo
        endif
        !write(*,*) 'after u33'

        if(keyEL .gt. 4) then
        U34(1:KR,1:KM,1:KN1,1:KIM,1:KRE) =
     &        U34(1:KR,1:KM,1:KN1,1:KIM,1:KRE)*WAVEL 
      do l=1,KR
      do j=1,KM
      do i=1,KN1
      do k=1,KIM
      do n=1,KRE
!fd
      if (.not. ieee_is_normal(U34(l,j,i,k,n))) then
        write(*,'(a,5i5,a,e16.6)') 'l,j,i,k,n: ',l,j,i,k,n,
     &  '  U34=',U34(l,j,i,k,n)
      end if
      enddo
      enddo
      enddo
      enddo
      enddo
        endif
      !write(*,*) 'after u34'

        if(keyEL .gt. 5) then
        DO J=1,KM
        U44(1:KR,J,1:KN1,1:KIM,1:KRE)=
     &        U44(1:KR,J,1:KN1,1:KIM,1:KRE)*WAVEL
        IF(ANGLE(J) .LE. ang_f44) THEN
      do l=1,kr
      do n=1,kre
      do k=1,kim
      do i=1,kn1
      if(U44(l,J,i,k,n) .le. 0.0)
     & write(*,*) l,n,k,i,J,ANGLE(J),U44(l,J,i,k,n)/WAVEL,
     & '  l,n,k,i,J,ANGLE(J),U44(l,J,i,k,n)'
      enddo
      enddo
      enddo
      enddo
        U44(1:KR,J,1:KN1,1:KIM,1:KRE)=
     &  LOG( U44(1:KR,J,1:KN1,1:KIM,1:KRE) )
!cl       write(*,*) 'U44 before intrpl orgnl:',U44(1,j,1,1,1)
        ENDIF ! ANGLE(J) .LE. ang_f44
        ENDDO ! J
        endif
!      write(*,*) 'after u44'
      RETURN
      ENDIF ! key .EQ. 4
										if(keySUB.eq.0) then
										 !T_CFM0=dtime(tarray)!+++
										 !T_CFM0=tarray(1)+tarray(2)
										endif

      IF(key_RD .eq. 2) then
! ** RECALCULATE ASPECT RATIO DISTRIBUTION (RDc()=SAREA/VOLUME)
! ** RDc()=RD()/RDc(); sumRD=sum(RDc())
!
! ** OBLATE
        do IR=1,KR
        if(R(IR).lt.1.) then              
        E=SQRT( 1.-R(IR)*R(IR) )
        xa1=LOG((1.+E)/(1.-E))/E
        RDc(IR)=1.5*(R(IR)**(-2./3.)+
     &            0.5*xa1*R(IR)**(4./3.))
!c ** PROLATE            
        elseif(R(IR).gt.1.) then          
        E=SQRT( 1.-1./R(IR)/R(IR) )
        xa2=ASIN(E)/E
        RDc(IR)=1.5*(R(IR)**(-2./3.)+
     &               xa2*R(IR)**(1./3.))
!c ** SPHERE
        elseif(R(IR).eq.1.) then             
        RDc(IR)=3.
        endif ! R()
!c ** WRITE ASPECT RATIO DISTRIBUTION 
!c          write(*,*) 'R=',R(IR),' B=',RDc(IR),
!c     &  ' 1/B=',1./RDc(IR),' RD=',RD(IR)

        enddo ! IR
        RDc(:KR)=RD(:KR)/RDc(:KR)
      ENDIF ! key_RD

      IF(key_RD .eq. 1) RDc(:KR)=RD(:KR)
      sumRD=sum(RDc(:KR))
      if(keySUB.eq.0) then
       write(*,*)
       do IR=1,KR 
       write(*,'(''R='',f8.4,4x,'' RDc='',e13.5)') R(IR),RDc(IR)
       enddo ! IR
       write(*,'(''sumRD='',e13.5)') sumRD
      endif

!c ** ALL SCATTERING MATRIX ELEMENTS
      DO IRE=1,KRE
      DO IIM=1,KIM
!c
!c **  SCATTERING MATRIX ELEMENTS
!c
      DO I=1,KN1
      IF(keyEL.gt.0) THEN
      DO J=1,KM

      DO IR=1,KR
      IF(NRATN.EQ.1) then
        RRATN=RATIO(1)
        IF(RRATN.NE.R(IR)) THEN
        WRITE(*,*) 'R=',R(IR),' .NE. RATIO=',RATIO(1)
        STOP 'in subroutine MATRIX_FIX 1'
!c        WRITE(*,*) 'R has been changed',
!c     &                  R(IR),' => ',RATIO(1)
!c        R(IR)=RRATN
        ENDIF
      ELSE
        RRATN=R(IR)
      ENDIF
      IF(NRATN.NE.1) THEN
      DO IRATN=1,NRATN-1
      IF(RRATN.GE.RATIO(IRATN).AND.RRATN.LE.RATIO(IRATN+1)) THEN
        L0=IRATN
        L1=IRATN+1
      ENDIF
      ENDDO
      IF(RRATN.LE.RATIO(1)) THEN
        IF(RRATN.LT.RATIO(1)) THEN
        WRITE(*,*) 'R=',R(IR),' is out of the range:',
     &                  RATIO(1),'< R <',RATIO(NRATN)
        STOP 'in subroutine MATRIX_FIX 2'
!c        WRITE(*,*) 'R has been changed',R(IR),' => ',RATIO(1)
        ENDIF
        L0=1
        L1=2
        R(IR)=RATIO(1)
        RRATN=RATIO(1)
      ENDIF
      IF(RRATN.GE.RATIO(NRATN)) THEN
        IF(RRATN.GT.RATIO(NRATN)) THEN
        WRITE(*,*) 'R=',R(IR),' is out of the range:',
     &                  RATIO(1),'< R <',RATIO(NRATN)
        STOP 'in subroutine MATRIX_FIX 3'
!c        WRITE(*,*) 'R has been changed',R(IR),
!c     &                            ' => ',RATIO(NRATN)
        ENDIF
        L0=NRATN-1
        L1=NRATN
        R(IR)=RATIO(NRATN)
        RRATN=RATIO(NRATN)
      ENDIF
      ELSE
        L0=1
        L1=1
      ENDIF
     
      RRATN1=RRATN

      X(1)=RATIO(L0)
      X(2)=RATIO(L1)
!c
!c ** U11 ->    UF11
!c
      if(keyEL .gt. 0) then
      Y(1)=U11(L0,J,I,IIM,IRE)
      Y(2)=U11(L1,J,I,IIM,IRE)
      UF11(J,I,IIM,IRE)=UF11(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)
	    endif
!      write(*,*) 'after UF11',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
!c
!c ** U12 ->    UF12
!c
      if(keyEL .gt. 1) then
      Y(1)=U12(L0,J,I,IIM,IRE)
      Y(2)=U12(L1,J,I,IIM,IRE)
      UF12(J,I,IIM,IRE)=UF12(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)
	    endif
!      write(*,*) 'after UF12',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
!c
!c ** U22 ->    UF22
!c
      if(keyEL .gt. 2) then
      Y(1)=U22(L0,J,I,IIM,IRE)
      Y(2)=U22(L1,J,I,IIM,IRE)
      UF22(J,I,IIM,IRE)=UF22(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)
      endif
!      write(*,*) 'after UF22',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
!c
!c ** U33 ->    UF33
!c
      if(keyEL .gt. 3) then
      Y(1)=U33(L0,J,I,IIM,IRE)
      Y(2)=U33(L1,J,I,IIM,IRE)
      UF33(J,I,IIM,IRE)=UF33(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)
      endif
!      write(*,*) 'after UF33',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
!c
!c ** U34 ->    UF34
!c
      if(keyEL .gt. 4) then
      Y(1)=U34(L0,J,I,IIM,IRE)
      Y(2)=U34(L1,J,I,IIM,IRE)
!      UTEST = LINEAR(X,Y,2,RRATN1)
!      if(ir .eq. 2 .and. j .eq. 180 .and. i .eq. 5 .and.
!     &   iim .eq. 1 .and. ire .eq. 1) then
!        write(*,*) X(1),X(2),Y(1),Y(2),RRATN1,RDc(IR),UTEST,
!     &  '  X(1),X(2),Y(1),Y(2),RRATN1,RDc(IR),UTEST'
!      endif
      UF34(J,I,IIM,IRE)=UF34(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)
      endif
!      write(*,*) 'after UF34',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
!c
!c ** U44 ->    UF44
!c
      if(keyEL .gt. 5) then
      Y(1)=U44(L0,J,I,IIM,IRE)
      Y(2)=U44(L1,J,I,IIM,IRE)
      UF44(J,I,IIM,IRE)=UF44(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)

      endif
!      write(*,*) 'after UF44',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
      ENDDO ! IR RD
      if(keyEL.gt.0) UF11(J,I,IIM,IRE)=UF11(J,I,IIM,IRE)/sumRD
      if(keyEL.gt.1) UF12(J,I,IIM,IRE)=UF12(J,I,IIM,IRE)/sumRD
      if(keyEL.gt.2) UF22(J,I,IIM,IRE)=UF22(J,I,IIM,IRE)/sumRD
      if(keyEL.gt.3) UF33(J,I,IIM,IRE)=UF33(J,I,IIM,IRE)/sumRD
      if(keyEL.gt.4) UF34(J,I,IIM,IRE)=UF34(J,I,IIM,IRE)/sumRD
      if(keyEL.gt.5) UF44(J,I,IIM,IRE)=UF44(J,I,IIM,IRE)/sumRD

      ENDDO   ! J KM
      ENDIF ! keyEL

!c
!c ** EXTINCTION & ABSORPTION
!c
      DO J=1,2
      DO IR=1,KR
      IF(NRATN.EQ.1) then
      RRATN=RATIO(NRATN)
      ELSE
      RRATN=R(IR)
      ENDIF
      IF(NRATN.NE.1) THEN
      DO IRATN=1,NRATN-1
      IF(RRATN.GE.RATIO(IRATN).AND.RRATN.LE.RATIO(IRATN+1)) THEN
          L0=IRATN
          L1=IRATN+1
      ENDIF
      ENDDO

      IF(RRATN.GE.RATIO(NRATN)) THEN
        L0=NRATN-1
        L1=NRATN
        R(IR)=RATIO(NRATN)
        RRATN=RATIO(NRATN)
      ENDIF
      IF(RRATN.LE.RATIO(1)) THEN
        L0=1
        L1=2
        R(IR)=RATIO(1)
        RRATN=RATIO(1)
      ENDIF
      ELSE
        L0=1
        L1=1
      ENDIF

      RRATN1=RRATN

      X(1)=RATIO(L0)
      X(2)=RATIO(L1)

      Y(1)=UEA(L0,J,I,IIM,IRE)
      Y(2)=UEA(L1,J,I,IIM,IRE)
      UFEA(J,I,IIM,IRE)=UFEA(J,I,IIM,IRE)+
     &                      LINEAR(X,Y,2,RRATN1)*RDc(IR)

!      write(*,*) 'after UFEA',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
      ENDDO ! IR KR

      UFEA(J,I,IIM,IRE)=UFEA(J,I,IIM,IRE)/sumRD     
      
      ENDDO ! J 2

      ENDDO ! I KN1
     
      ENDDO ! IIM
      ENDDO ! IRE

      !goto 554
      write(iu_output,'(a)') 'in kernels_tb_shd_equal'
      write(iu_output,'(es12.4,a,es12.4)')
     &             minval(UFEA(1,1:KN1,1:KIM,1:KRE)),
     &' <= UE  <=',maxval(UFEA(1,1:KN1,1:KIM,1:KRE))
      write(iu_output,'(es12.4,a,es12.4)')
     &             minval(UFEA(2,1:KN1,1:KIM,1:KRE)),
     &' <= UA  <=',maxval(UFEA(2,1:KN1,1:KIM,1:KRE))

554   continue

                                         !T_CFM0=dtime(tarray)
                                         !T_CFM0=tarray(1)+tarray(2)
                                         !T_CFM=T_CFM+T_CFM0 !+++
      IF(key.EQ.1) THEN
!c
!c *** SAVE FIXED matrices
!c


      open(10,file=TRIM(dir_name_F)//'Rkernel1_fix_00.txt',
     &                                      status='unknown')
      IF(keyEL.gt.0) THEN
       open(11,file=TRIM(dir_name_F)//'Rkernel1_fix_11.txt',
     &                                      status='unknown')
      if(keyEL.gt.1) 
     & open(12,file=TRIM(dir_name_F)//'Rkernel1_fix_12.txt', 
     &                                      status='unknown')
      if(keyEL.gt.2)	  
     & open(13,file=TRIM(dir_name_F)//'Rkernel1_fix_22.txt', 
     &                                      status='unknown')
      if(keyEL.gt.3)
     & open(14,file=TRIM(dir_name_F)//'Rkernel1_fix_33.txt', 
     &                                      status='unknown')
      if(keyEL.gt.4)     
     & open(15,file=TRIM(dir_name_F)//'Rkernel1_fix_34.txt', 
     &                                      status='unknown')
      if(keyEL.gt.5)
     & open(16,file=TRIM(dir_name_F)//'Rkernel1_fix_44.txt', 
     &                                      status='unknown')
      ENDIF
        DO II=1,NNEL
        if(key_fx.eq.1) then
          WRITE(10+II-1,10) grid1(nn1),grid1(nn2),R(KR)
        else ! key_fx=0 
          WRITE(10+II-1,*) key_RD,                               
     &  '  key_RD=1-volume mixture, 2-surface area mixture'      
          WRITE(10+II-1,10) grid1(nn1),grid1(nn2)
          WRITE(10+II-1,*) KR,' a number of grid aspect ratios'  
          WRITE(10+II-1,*) 'aspect ratio distribution'           
          DO I=1,KR                                              
          WRITE(10+II-1,11) R(I),RDc(I)                          
          ENDDO ! I 
        endif ! key_fx                                               
        WRITE(10+II-1,*) -nn3,' a number of grid radii'
        ENDDO ! II
        DO II=2,NNEL
        WRITE(10+II-1,*) KM,' a number of scattering angles'
        WRITE(10+II-1,15) ANGLE(1:KM)
        ENDDO ! II
        DO II=1,NNEL
        WRITE(10+II-1,'(2e13.5,'' real  refr. index'')') 
     &                ARE(1),ARE(KRE)
        WRITE(10+II-1,'(2e13.5,'' imag refr. index'')') 
     &                                -AIM(1),-AIM(KIM)
        WRITE(10+II-1,*) KRE,-KIM,' number of intervals for opt. const'       
        ENDDO ! II

        DO IRE=1,KRE
        DO IIM=1,KIM
        DO II=1,NNEL
        if(key_fx.eq.1) then
        WRITE(10+II-1,'(i3,E13.5,'' element number,ratn'')') 
     &                  NEL(II-1),R(KR)  
        else ! key_fx=0
        WRITE(10+II-1,'(i3,E13.5,'' element number'')') 
     &                  NEL(II-1)  
        endif ! key_fx       
        WRITE(10+II-1,'(E11.4,2E15.7,'' wavel,rreal,rimag'')')
     &         WAVEL,ARE(IRE),-AIM(IIM)
        ENDDO ! II

      IF(keyEL.gt.0) THEN
        DO I=nn1,nn2
        WRITE(11,11) UF11(1:KM,I,IIM,IRE)
        ENDDO ! I
        if(keyEL.gt.1) then
        DO I=nn1,nn2
        WRITE(12,11) UF12(1:KM,I,IIM,IRE)
      ENDDO ! I
      endif
        if(keyEL.gt.2) then	
        DO I=nn1,nn2
        WRITE(13,11) UF22(1:KM,I,IIM,IRE)
      ENDDO
      endif
        if(keyEL.gt.3) then	
        DO I=nn1,nn2		
        WRITE(14,11) UF33(1:KM,I,IIM,IRE)
	    ENDDO ! I
      endif
	    if(keyEL.gt.4) then
        DO I=nn1,nn2
        WRITE(15,11) UF34(1:KM,I,IIM,IRE)
        ENDDO ! I
      endif
	    if(keyEL.gt.5) then
        DO I=nn1,nn2		
        WRITE(16,11) UF44(1:KM,I,IIM,IRE)
        ENDDO ! I
        endif 
        ENDIF

        WRITE (10,*)'EXTINCTION (1/km, for d()/dlnr m3/m3*km)'
        WRITE (10,11) UFEA(1,nn1:nn2,IIM,IRE)
        WRITE (10,*)'ABSORPTION (1/km, for d()/dlnr m3/m3*km)'
        WRITE (10,11) UFEA(2,nn1:nn2,IIM,IRE)
        ENDDO ! IRE
        ENDDO ! IIM
      close(10)
      IF(keyEL.gt.0) THEN
      close(11)
      if(keyEL.gt.1) close(12)
      if(keyEL.gt.2) close(13)
      if(keyEL.gt.3) close(14)
      if(keyEL.gt.4) close(15)
      if(keyEL.gt.5) close(16)
      ENDIF ! keyEL

      WRITE(*,*) 'Fixed kernels have been calculated and saved'
      IF(key_RD.EQ.1) WRITE(*,*) 'Volume mixture of spheroids'
      IF(key_RD.EQ.2) WRITE(*,*) 'Surface area mixture of spheroids'
!c
!c *** DEALLOCATE ARRAYS
!c
      IF(keyEL.gt.0) THEN
      DEALLOCATE(U11,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate U11 array'
      DEALLOCATE(UEA,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate UEA array'
      if(keyEL.gt.1) then
        DEALLOCATE(U12,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U12 array'
      endif
      if(keyEL.gt.2) then	  
        DEALLOCATE(U22,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U22 array'
      endif
      if(keyEL.gt.3) then
        DEALLOCATE(U33,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U33 array'
      endif
      if(keyEL.gt.4) then
        DEALLOCATE(U34,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U34 array'
      endif
      if(keyEL.gt.5) then
        DEALLOCATE(U44,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U44 array'
      endif 
      ENDIF ! keyEL
      ENDIF ! key=1

      IF(keyEL.gt.0) THEN
      UF11(1:KM,1:KN1,1:KIM,1:KRE)=
     &     LOG( UF11(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL )
      if(keyEL.gt.1)
     &  UF12(1:KM,1:KN1,1:KIM,1:KRE)=
     &         UF12(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL
      if(keyEL.gt.2)
     &  UF22(1:KM,1:KN1,1:KIM,1:KRE)=
     &         LOG( UF22(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL )
      if(keyEL.gt.3) then
      DO J=1,KM
        UF33(J,1:KN1,1:KIM,1:KRE)=
     &         UF33(J,1:KN1,1:KIM,1:KRE)*WAVEL
      IF(ANGLE(J) .LE. ang_f33)
     &     UF33(J,1:KN1,1:KIM,1:KRE)=
     &         LOG( UF33(J,1:KN1,1:KIM,1:KRE) )
      ENDDO ! J
      endif
      if(keyEL.gt.4)
     &  UF34(1:KM,1:KN1,1:KIM,1:KRE)=
     &         UF34(1:KM,1:KN1,1:KIM,1:KRE)*WAVEL	  
      if(keyEL.gt.5) then
      DO J=1,KM	  	  
        UF44(J,1:KN1,1:KIM,1:KRE)=
     &         UF44(J,1:KN1,1:KIM,1:KRE)*WAVEL
      IF(ANGLE(J) .LE. ang_f44)
     &     UF44(J,1:KN1,1:KIM,1:KRE)=
     &         LOG( UF44(J,1:KN1,1:KIM,1:KRE) )
      ENDDO ! J
      endif
      ENDIF ! keyEL
                                       
      !if(keySUB .eq. 0) then
      !write(6,*) 
      !write(6,*) 
      !&        '------------------ T I M I N G ------------------' 
      !WRITE(*,62) T_CFM/60.     
      !endif
62    format('  Calcul. fixed kernels ... ',f8.3,' min.')

      ENDIF ! key.NE.2
    
10    FORMAT(3E15.7,' rmin, rmax, RATIO')                    
11    FORMAT(7E15.7)                    
12    FORMAT(3E15.7,I4,'    wavelength, n, k, NEL')                    
15    FORMAT(7F12.2)                    
       
      RETURN
      END SUBROUTINE MATRIX_FIX

!c **************************************************************** 
