!
!  main.f95
!  FortranFindSolution
!
!  Created by Константин Шмирко on 2020-12-30.
!  Copyright 2020 Константин Шмирко. All rights reserved.
!
program FitMeasurements
	use mo_DLS
	use ObjFuncMod
	
	implicit None
	integer  I, nop, NP, T, II, arg_count, ftype, JJ, ires
	integer*8 opt
	real v, Cr, Fl, Fh
	character(len=255) fname, data_fname, iniFname
	external ObjectiveFunction, fmin
	include 'nlopt.inc'
	
	
	! Проверяем наличие аргументов командной строки
	arg_count = command_argument_count()
	if ( arg_count .ne. 1 ) then
		write (*, *) "./FitMeasurements iniFilename.ini"
		STOP
	endif
	
			
	fname    = "./input.dat"
	iniFname = "./inifile.ini"
	
	! Get first command line argument - init file name
	call get_command_argument(arg_count, iniFname)
	
	! Initialize all module variables for solving inversion poblem
	call InitObjFuncVariables()
	
	! Initialize all shared variables from the ini file
	call ReadIniFile(iniFname, fname)
	
	call DLS_read_input(fname)
	
	! NDP = 0, then berfore calculations we load tables from file
	NDP=0
	
	! Allocate internal arrays last parameter = 1
	call alloc_DLS_array(key,keyEL,1)
	
	
	! Number of population vectors:
	NP = 10 * SearchParamsCount + 10
	! Number of generations
	T = NGEN
	! Crossover variable: Cr \in [0,1]
	Cr = 0.85D0 ! manually set

	! Scale factor determined by dither:
	Fl = SQRT(1.0D0 - 0.50D0* Cr )
	Fh = 0.950D0 ! L,H order should be fine as long as Cr > 0.2
	
	!call powerlaw(KN,1,1.0,-4.3,0.1,1.0,RRR,AR,AC,KNpar)
	
	
	X = 0.5*(LoParamVal+UpParamVal)
	
	
	
	
	DO II=1, InputVectorsCount	
		! Инициализация солвера с методом глобальной оптимизации
		if (func_type .eq. FunctLogNormal) then
			call nlo_create(opt, NLOPT_GN_ISRES, LnParamsCount)
		else if (func_type .eq. FunctPowerLaw) then
			call nlo_create(opt, NLOPT_GN_ISRES, PowParamsCount)
		else
			stop 'Неподдерживаемый тип аппроксимирующей функции'
		endif

		call nlo_set_lower_bounds(ires, opt, LoParamVal)
		call nlo_set_upper_bounds(ires, opt, UpParamVal)
		call nlo_set_min_objective(ires, opt, fmin, 0)
		call nlo_set_xtol_rel(ires, opt, 1.0e-5)
		call nlo_set_ftol_rel(ires, opt, 1.0e-5)
		call nlo_set_maxeval(ires, opt, 100000)
		
		
		! Готовим вектор значений и абсолютной ошибки для вычислений
		DO JJ = 1, InpParamsCount
			Ym(JJ) = datum(JJ, II)*1e-9
			Yerr(JJ) = daterr(JJ, II) * Ym(JJ) / 100.0
		END DO
		
		call nlo_optimize(ires, opt, X, v) 
		
		if (func_type .eq. FunctLogNormal) then
			print '(I4, "   RETCODE: ", I4, " SOLUTION: ", E13.3, 2F7.3, F7.2, F7.4, " FIT: ", F5.1, $)', II, ires, X(1:5), v
		elseif (func_type .eq. FunctPowerLaw) then
			print '(I4, "   RETCODE: ", I4, " SOLUTION: ", E13.3, 2F7.2, F7.4, " FIT: ", F5.1, $)', II, ires, X(1:4), v		
		endif
		
		call nlo_destroy(opt)
		Call PrintReport3(v, X)
		
		
	END DO
	
	
	
	
	! Deallocate internal arrays last parameter = 2
	call alloc_DLS_array(key,keyEL,2)
	
contains
	subroutine ReadIniFile(iniFname, fname)
		use ObjFuncMod
		implicit none
		character(len=255), intent(in) 		:: iniFname
		character(len=255), intent(inout) :: fname
		INTEGER I, J
		
		open(101, FILE=trim(iniFname), status='old')
		READ(101, *) fname
		READ(101, *) WavelengthCount
		READ(101, *) AlphaInpParamsCount
		READ(101, *) DepolInpParamsCount
		READ(101, *) (Wvl(I),I=1,WavelengthCount)
		READ(101, *) discrKind
		READ(101, *) func_type
		READ(101, *) KN
		READ(101, *) Rmin, Rmax
		READ(101, *) threshold
		READ(101, *) NGEN
		
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
			READ(101, *) (daterr(J, I),J=1,InpParamsCount)
		END DO
		close(101)
		
		
	end subroutine ReadIniFile
	
	subroutine PrintReport(gof, fit)
		use mo_DLS
		use ObjFuncMod
		use mo_par_DLS, only : KN1par
		implicit none
		
		
		real, intent(in) :: gof, fit(:)
		real						 :: tmpr(KN1par), a, b, reff, dlogr, Ntot
		
		! REPORT RESULTS
		PRINT *
		PRINT '("REPORT:")'
		PRINT '("=======")'
		PRINT *
		PRINT *, "The program found a solution to the inverse scattering problem in the spheroid"
		PRINT *, "approximation, using a log-normal function as the distribution'sshape. The search"
		PRINT *, "for a solution consisted in the selection of the parameters of this distribution "
		PRINT *, "that minimize residual function."
		PRINT *
		PRINT *, "The original solution vector contains five components: 3 backscatter coefficients"
		PRINT *, "at wavelengths (0.355, 0.532, and 1.064 mkm) and 2 attenuation coefficients "
		PRINT *, "at wavelengths (0.355 and 0.532 mkm)."
		PRINT *
	
	
	
	
		PRINT '("GOODNESS OF FIT:", F13.2, "%",/)', gof
	
		PRINT '("                 ", 5A13)', "C1", "C2", "C3", "C4", "C5"
		PRINT '("                 ", 5A13)', "==========", "==========", "==========", "==========", "=========="
		PRINT '("Final solution = ", 5E13.6)', fit(1:SearchParamsCount)
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
		PRINT *, "VOLUME SIZE DISTRIBUTION"
		PRINT '(A13,A13)', "R, mkm", "SD"
		
		A = 0.0
		Ntot = 0.0
		DO I=1, KN
			PRINT '(2E13.3)', RRR(I), AR(I)
			tmpr(I) = AR(I) / (4.1888 * RRR(I)**3)
			Ntot = Ntot + tmpr(I) * RRR(I)
		END DO
		
		dlogr = LOG(RRR(2)/RRR(1))
		Ntot = Ntot * dlogr
		print '("Ntot = ", F10.4, " cm-3")', Ntot*1e12
		! Вычислим эффективный радиус частиц, для этого возьмем и переведем
		! нашу концентрацию в 
		
		reff = 0.0
		a= 0.0
		b= 0.0
		
		do I=1, KN
			a = a + (RRR(I)**4 * tmpr(I))
			b = b + (RRR(I)**3 * tmpr(I))
		end do
		reff = a / b
		!A = X(1)*0.238732414637843/(EXP(3*X(3)+4.5*X(2)**2))
		A = Yc(2)/(Ntot*reff**2)
		
		
		print '("REFF = ", F7.3)', reff
		print '("A    = ", F7.3, " LG(A) = ", F7.3)', A, LOG10(A)
	end subroutine PrintReport
	
	subroutine PrintReport1(gof)
		use mo_DLS
		use mo_par_DLS, only : KN1par
		use ObjFuncMod
		implicit none
		
		
		real, intent(in) :: gof
		real						 :: dlogr, reff, a, b
		real						 :: tmpr(KN1par)
		
		! REPORT RESULTS
		PRINT *
		PRINT '("REPORT:")'
		PRINT '("=======")'
		PRINT *
		PRINT *, "The program found a solution to the inverse scattering problem in the spheroid"
		PRINT *, "approximation, using a log-normal function as the distribution'sshape. The search"
		PRINT *, "for a solution consisted in the selection of the parameters of this distribution "
		PRINT *, "that minimize residual function."
		PRINT *
		PRINT *, "The original solution vector contains five components: 3 backscatter coefficients"
		PRINT *, "at wavelengths (0.355, 0.532, and 1.064 mkm) and 2 attenuation coefficients "
		PRINT *, "at wavelengths (0.355 and 0.532 mkm)."
		PRINT *
	
	
	
	
		PRINT '("GOODNESS OF FIT:", F13.2, "%",/)', gof
	
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
		PRINT *, "VOLUME SIZE DISTRIBUTION"
		PRINT '(A13,A13)', "R, mkm", "SD"
		
		DO I=1, KN
			PRINT '(2E13.3)', RRR(I), AR(I)
			tmpr(I) = AR(I) / (4.1888 * RRR(I)**3)
		END DO
		
		! Вычислим эффективный радиус частиц, для этого возьмем и переведем
		! нашу концентрацию в 
		dlogr = LOG(RRR(2)/RRR(1))
		reff = 0.0
		a= 0.0
		b= 0.0
		
		do I=1, KN
			a = a + (RRR(I)**4 * tmpr(I))
			b = b + (RRR(I)**3 * tmpr(I))
		end do
		reff = a / b
		A = Yc(2)/(X(1)*reff**2)
		
		
		print '("REFF = ", F7.3)', reff
		print '("A    = ", F7.3, " LG(A) = ", F7.3)', A, LOG10(A)
	end subroutine PrintReport1
	
	subroutine PlotSizeDistribution(fname)
		implicit none
		character(len=*), intent(in)	::	fname
		
	end subroutine PlotSizeDistribution
	
	subroutine PrintReport3(gof, fit)
		use mo_DLS
		use ObjFuncMod
		use mo_par_DLS, only : KN1par
		implicit none
		
		
		real, intent(in) :: gof, fit(:)
		real						 :: tmpr(KN1par), a, b, reff, dlogr, Ntot
		
		! REPORT RESULTS
		A = 0.0
		Ntot = 0.0
		DO I=1, KN
			tmpr(I) = AR(I) / (4.1888 * RRR(I)**3)
			Ntot = Ntot + tmpr(I) * RRR(I)
		END DO
		
		dlogr = LOG(RRR(2)/RRR(1))
		Ntot = Ntot * dlogr
		print '(" Ntot = ", F8.2, " cm-3 ", $)', Ntot*1e12
		! Вычислим эффективный радиус частиц, для этого возьмем и переведем
		! нашу концентрацию в 
		
		reff = 0.0
		a= 0.0
		b= 0.0
		
		do I=1, KN
			a = a + (RRR(I)**4 * tmpr(I))
			b = b + (RRR(I)**3 * tmpr(I))
		end do
		reff = a / b
		!A = X(1)*0.238732414637843/(EXP(3*X(3)+4.5*X(2)**2))
		A = Yc(2)/(Ntot*reff**2)
		
		
		print '("REFF = ", F5.2, $)', reff
		print '(" A    = ", F5.3, " LG(A) = ", F6.3)', A, LOG10(A)
	end subroutine PrintReport3
end program FitMeasurements
