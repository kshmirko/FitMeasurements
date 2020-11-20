program FitMeasurements
	use mo_DLS
	use ObjFuncMod
	USE DE_mod
	
	
	implicit None
	integer  I, nop, NP, T
	real v, Cr, Fl, Fh
	character(len=255) fname
	fname = "./input.dat"
	
	call InitObjFuncVariables()
	
	call DLS_read_input(fname)
	
	print '("Wl=", F7.3)', WL
	
	! NDP = 0, then berfore calculations we load tables from file
	NDP=0
	
	! Allocate internal arrays last parameter = 1
	call alloc_DLS_array(key,keyEL,1)
	
	Wvl = [0.355, 0.532, 1.064]
	Ym = [0.037000,0.042000,0.029000,1.414000,1.218000]
	LoParamVal = (/0.01, 0.10, 0.01, 1.3, 0.00000001/)
	UpParamVal = (/10.0, 0.48, 1.2 , 1.56, 0.05/)
	X =  (/0.05, 0.45, 0.8, 1.45, 0.01/)

	KN=22
	RMin = 0.05
	RMax = 15.0
	nop = 5
	discrKind = DiscrKindRMS
	
	! Number of population vectors:
	NP = 10 * nop + 10
	! Number of generations
	T = 1000
	! Crossover variable: Cr \in [0,1]
	Cr = 0.85D0 ! manually set

	! Scale factor determined by dither:
	Fl = SQRT(1.0D0 - 0.50D0* Cr )
	Fh = 0.950D0 ! L,H order should be fine as long as Cr > 0.2
	
	
	CALL Diff_Evol( nop, NP, LoParamVal, UpParamVal, T, Fl, Fh, Cr, 1, X, v, 100, ObjectiveFunction )
  CALL ObjectiveFunction(nop, X, v)
	
	
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
	PRINT '("                 ", 5A13)', "----------", "----------", "----------", "----------", "----------"
	PRINT '("Final solution = ", 5F13.6)', X
	PRINT *, "                   ----------   ----------   ----------   ----------   ----------"
	PRINT *
	
	PRINT *, "MEASUREMENT AND CALCULATED DATA"
	PRINT *, "==============================="
	PRINT *
	PRINT '(A3,3A13," %")', "#", "Y meas", "Y calc", "Error,"
	
	DO I=1, SearchParamsCount
		PRINT '(I3, 2E13.3, F13.2, " %")', I, Ym(I), Yc(I), ((Ym(I)-Yc(I))/Ym(I))*100
	END DO
	
	PRINT *
	PRINT *, "OPTICAL PROPERTIES (CALC)"
	PRINT '(5A13)', "WVL", "Ext", "Bsc", "BLR", "LDR"
	DO I=1, WavelengthCount 
		PRINT '(F13.3,2E13.3,2F13.3)', Wvl(I), Ext(I), Bsc(I), Blr(I), Ldr(I)
	END DO
	
	! Deallocate internal arrays last parameter = 2
	call alloc_DLS_array(key,keyEL,2)

end program FitMeasurements