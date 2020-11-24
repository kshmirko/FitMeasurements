program FitMeasurements
	use mo_DLS
	use ObjFuncMod
	USE DE_mod
! 	use f90getopt
	
	implicit None
	integer  I, nop, NP, T, II
	real v, Cr, Fl, Fh
	character(len=255) fname, data_fname, iniFname
!   type(option_s)	::	opts(2)
!
! 	opts(1) = option_s( "func_type", .true., 'f' )
! 	opts(2) = option_s( "help", .false., 'h')
!
! 	if (command_argument_count() .eq. 0 ) then
! 		write (*,*) "Available Options: --func_type=0|1 -f0|1 --help -h"
! 		STOP
! 	end if
!
!
! 	do
! 		select case( getopt( "f:h", opts ) ) ! opts is optional (for longopts only)
! 	  	case( char(0) )
! 	    	exit
! 	    case( 'f' )
! 				print *, 'func_type/f=', optarg
! 				STOP
! 	    case( 'h' )
! 				print *, 'help-screen'
! 				STOP
! 	  end select
! 	end do
			
	fname = "./input.dat"
	data_fname = "./input.lst"
	iniFname = "./inifile.ini"
	
	call InitObjFuncVariables()
	
	
	! Call init variables from file
	call ReadIniFile(iniFname)
	call DLS_read_input(fname)
	
	! NDP = 0, then berfore calculations we load tables from file
	NDP=0
	
	! Allocate internal arrays last parameter = 1
	call alloc_DLS_array(key,keyEL,1)
	
	!Wvl = [0.355, 0.532, 1.064]
	!Ym = [0.037000,0.042000,0.029000,1.414000,1.218000]
	!LoParamVal = (/0.01, 0.15, 0.05, 1.3, 0.00000001/)
	!UpParamVal = (/10.0, 0.48, 1.2 , 1.56, 0.05/)
	!X =  (/5.0, 0.45, 0.8, 1.45, 0.01/)

	!KN=22
	!RMin = 0.05
	!RMax = 15.0
	!discrKind = DiscrKindRMS
	
	! Number of population vectors:
	NP = 10 * SearchParamsCount + 10
	! Number of generations
	T = 1000
	! Crossover variable: Cr \in [0,1]
	Cr = 0.85D0 ! manually set

	! Scale factor determined by dither:
	Fl = SQRT(1.0D0 - 0.50D0* Cr )
	Fh = 0.950D0 ! L,H order should be fine as long as Cr > 0.2
	
	!call powerlaw(KN,1,1.0,-4.3,0.1,1.0,RRR,AR,AC,KNpar)
	
	DO II=1, InputVectorsCount
		Ym(1:InpParamsCount) = datum(1:InpParamsCount, II)
		
		
		
		CALL Diff_Evol( SearchParamsCount, NP, LoParamVal, UpParamVal, T, Fl, Fh, Cr, 1, X, v, 100, ObjectiveFunction )
  	CALL ObjectiveFunction(SearchParamsCount, X, v)
		
		! REPORT RESULTS
		PRINT '("REPORT:")'
		PRINT *
		PRINT *, "The program found a solution to the inverse scattering problem in the spheroid"
		PRINT *, "approximation, using a log-normal function as the distribution. The search for"
		PRINT *, "a solution consisted in the selection of the parameters of this distribution "
		PRINT *, "that minimize residual function."
		PRINT *
		PRINT *, "The original solution vector contains five components: 3 backscatter coefficients"
		PRINT *, "at wavelengths (0.355, 0.532, and 1.064 μm) and 2 attenuation coefficients "
		PRINT *, "at wavelengths (0.355 and 0.532 μm)."
		PRINT *
	
	
	
	
		PRINT '("GOODNESS OF FIT:", F13.2, "%",/)', v
	
		PRINT '("                 ", 5A13)', "CM", "SM", "RMM", "RN", "RK"
		PRINT '("                 ", 5A13)', "==========", "==========", "==========", "==========", "=========="
		PRINT '("Final solution = ", 5F13.6)', X
		PRINT *, "                   ----------   ----------   ----------   ----------   ----------"
		PRINT *
	
		PRINT *, "MEASURED AND CALCULATED DATA"
		PRINT *, "============================"
		PRINT *
		PRINT '(A3,3A13," %")', "#", "Y meas", "Y calc", "Error,"
	
		DO I=1, InpParamsCount
			PRINT '(I3, 2E13.3, F13.2, " %")', I, Ym(I), Yc(I), ((Ym(I)-Yc(I))/Ym(I))*100
		END DO
	
		PRINT *
		PRINT *, "OPTICAL PROPERTIES (CALC)"
		PRINT '(5A13)', "WVL", "Ext", "Bsc", "BLR", "LDR"
		DO I=1, WavelengthCount 
			PRINT '(F13.3,2E13.3,2F13.3)', Wvl(I), Ext(I), Bsc(I), Blr(I), Ldr(I)
		END DO
		PRINT *
		PRINT *, "SIZE DISTRIBUTION"
		PRINT '(A13,A13)', "R, mkm", "SD"
		
		DO I=1, KN
			PRINT '(3E13.3)', RRR(I), AR(I), SD(I)
		END DO
	END DO
		
	
	
	
	
	! Deallocate internal arrays last parameter = 2
	call alloc_DLS_array(key,keyEL,2)
	
contains
	subroutine ReadIniFile(iniFname)
		use ObjFuncMod
		implicit none
		character(len=255), intent(in) :: iniFname
		INTEGER I, J
		
		open(101, FILE=trim(iniFname), status='old')
		READ(101, *) WavelengthCount
		READ(101, *) AlphaInpParamsCount
		READ(101, *) DepolInpParamsCount
		READ(101, *) (Wvl(I),I=1,WavelengthCount)
		READ(101, *) discrKind
		READ(101, *) func_type
		READ(101, *) KN
		READ(101, *) Rmin, Rmax
		
		if(func_type .eq. FunctLogNormal) then
			SearchParamsCount = LnParamsCount
		else if(func_type .eq. FunctPowerLaw) then
			SearchParamsCount = PowParamsCount
		end if
		
		InpParamsCount = WavelengthCount+AlphaInpParamsCount+DepolInpParamsCount
		READ(101, *) (LoParamVal(I),I=1,SearchParamsCount)
		READ(101, *) (UpParamVal(I),I=1,SearchParamsCount)
		READ(101, *) InputVectorsCount
		DO I=1, InputVectorsCount
			READ(101, *) (datum(J, I),J=1,InpParamsCount)
		END DO
		close(101)
		
		
	end subroutine ReadIniFile
end program FitMeasurements