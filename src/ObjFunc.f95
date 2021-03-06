module ObjFuncMod
	! при работе с модулями в фортране важно специально инициализировать
	! все его переменные, особенно в случае масивов
	! Основные используемые переменные:
	! Ym(MaxInpParamsCount)		- вектор измерений
	! Yc(MaxInpParamsCount)		- вектор расчетных величин
	! Yerr(MaxInpParamsCount)	-	вектор ошибок измерений в процентах
	! 													в дальнейшем, он преобразуется в вектор
	!														абсолютных ошибок Yerr = Yerr * Ym
	!	Добавил вариант DiscrKindChiSqr - невязка считается как хи-квадрат
	
	implicit none
	public
	integer, parameter	:: MaxInputVectorsCount = 100
	integer, parameter	:: MaxWavelengthsCount = 3
	
	
	integer, parameter	:: MaxSearchParamsCount = 10
	integer, parameter	:: MaxInpParamsCount    = 10
	
	integer, parameter	:: DiscrKindRMS			 = 0
	integer, parameter	:: DiscrKindMAXABS   = 1
	integer, parameter	:: DiscrKindChiSqr	 = 2
	
	integer, parameter	:: FunctLogNormal		 = 0
	integer, parameter	:: FunctPowerLaw		 = 1
	integer, parameter	:: FunctManualInput	 = 2
	
	integer, parameter	:: LnParamsCount		 = 5
	integer, parameter	:: PowParamsCount		 = 4
	  
	
	real, dimension(MaxInpParamsCount, MaxInputVectorsCount)	:: datum, &
																															 &daterr
	real, dimension(MaxInpParamsCount)		:: Ym,	Yc, Yerr
	
	real, dimension(MaxSearchParamsCount)	:: LoParamVal, UpParamVal, X
	real, dimension(MaxWavelengthsCount)		:: Bsc, Ext, Wvl, Absb, Ldr, Blr
	
	real	::	Rmin, Rmax, threshold
	Integer:: discrKind, func_type
	integer	:: AlphaInpParamsCount, DepolInpParamsCount, InputVectorsCount
	integer :: WavelengthCount 
	integer	:: InpParamsCount, SearchParamsCount, NGEN
	
	
	! Для случая использования внешней процепуры поиска решения создадим
	! отдельный ип данных, который бы мог хранить используемые параметры
	
	type SupplementalParams
		real, dimension(MaxInpParamsCount, MaxInputVectorsCount)	:: datum, &
																															 &daterr
		real, dimension(MaxInpParamsCount)		:: Ym,	Yc, Yerr
   	real, dimension(MaxSearchParamsCount)	:: LoParamVal, UpParamVal, X
		real, dimension(MaxWavelengthsCount)		:: Bsc, Ext, Wvl, Absb, Ldr, Blr
		
		real	::	Rmin, Rmax, threshold
		Integer:: discrKind, func_type
		integer	:: AlphaInpParamsCount, DepolInpParamsCount, InputVectorsCount
		integer :: WavelengthCount 
		integer	:: InpParamsCount, SearchParamsCount, NGEN																												
	end type
	
contains
	
	! Специальная функция инициализации переменных модуля
	subroutine InitObjFuncVariables()
		NGEN = 1000
		datum = 0.0
		daterr = 0.0
		Ym = 0.0
		Yc = 0.0
		Yerr = 0.0
		LoParamVal = 0.0
		UpParamVal = 0.0
		Bsc = 0.0
		Ext = 0.0
		Wvl = 0.0
		Absb = 0.0
		Ldr = 0.0
		Blr = 0.0
		X= 0.0
		Rmin = 0.05
		Rmax = 15.0
		discrKind = DiscrKindRMS
		func_type = FunctLogNormal
		WavelengthCount = 3
		InpParamsCount = 5
		SearchParamsCount =5
	end subroutine InitObjFuncVariables
	
	! subroutine ObjectiveFunction(NP, X, func_val)
! 		use mo_DLS
! 		implicit none
! 		integer, intent(IN)	::	NP
! 		real, intent(IN) 		::	X(NP)
! 		real, intent(inout) ::	func_val
! 		real tmp
! 		logical b_print
!
! 		integer I
!
! 		! Chreck bounds
! 		DO I=1, NP
! 			if ( ( X(I).lt.LoParamVal(I) ) .or.  ( X(I).gt.UpParamVal(I) ) ) then
! 				func_val = 99999.0
! 				return
! 			end if
! 		END DO
!
!
!
! 		if ( (func_type .eq. FunctLogNormal) .and. (NP .eq. LnParamsCount) ) then
! 			call SIZEDIS2(-KN,1,X(1),X(2),X(3),RMIN,RMAX,RRR,AR,AC)
! 			! Last two items in X are always RN and RK
! 			RN = X(NP-1)
! 			RK = X(NP)
! 			SD(1:KN) = AR(1:KN)
! 		else if ( (func_type .eq. FunctPowerLaw) .and. (NP .eq. PowParamsCount) ) then
! 			call powerlaw(KN,X(1),X(2),RMIN,RMAX,RRR,AR,AC, KNpar)
! 			! Last two items in X are always RN and RK
! 			RN = X(NP-1)
! 			RK = X(NP)
! 			SD(1:KN) = AR(1:KN)
! 		else if ( func_type .eq. FunctManualInput ) then
! 			print *, "USER FUNC"
! 		else
! 			print *, "WRONG PARAMETERS COUNT OF FITTING FUNCTION TYPE"
! 			STOP
! 		endif
!
!
!
!
! 		DO I=1, WavelengthCount
! 			IF (NDP==0) then
! 				print *, "LOADING DATABASE OF SPHEROIDS..."
! 				b_print = .true.
! 			endif
!
! 			WL = WVL(I)
!
! 			CALL OPTCHAR(NDP)
!
! 			if(b_print .eqv. .true.) then
! 				print *, "DONE."
! 				b_print = .false.
! 			endif
!
!
! 			Ext(I) = xext
! 			Absb(I) = xabs
! 			Bsc(I) = xext / xblr
! 			Ldr(I) = xldr
! 			Yc(I) = Bsc(I)
! 			Blr(I) = xblr
!
! 		END DO
!
! 		DO I=1, AlphaInpParamsCount
! 			Yc(I+WavelengthCount) = Ext(I)
! 		END DO
!
! 		! В зависимости от того, какой свособ расчета невязки выбран
! 		! производим вычисления
! 		if(discrKind==DiscrKindRMS) then
! 			!случай невязки в виде root mean square
! 			func_val = 0.0
!
! 			! Переписали код короче
! 			!func_val = sum( ((Ym - Yc)/Ym)**2 )
! 			DO I=1, SearchParamsCount
! 				func_val = func_val + (( Ym(I)-Yc(I) )/Ym(I))**2.0
! 			END DO
!
! 			func_val = SQRT(func_val/SearchParamsCount)
! 			! переводим в проценты полученную величину
! 			func_val = func_val * 100
! 		else if (discrKind==DiscrKindMAXABS) then
! 			! случай максимального абсолютного отклонения
! 			func_val=0.0
!
! 			! Перевисали код короче
! 			!func_val = MAXVAL( (Ym-Yc)/Ym )
!  			DO I=1, SearchParamsCount
!  				tmp = ABS(( Ym(I)-Yc(I) )/Ym(I))
!  				if (func_val .gt. tmp) then
!  					func_val = tmp
!  				end if
!  			END DO
! 			! переводим в проценты полученную величину
! 			func_val = func_val * 100
! 		else if (discrKind==DiscrKindChiSqr) then
! 			! случай максимального абсолютного отклонения
! 			func_val=0.0
!
!  			DO I=1, SearchParamsCount
!  				tmp = (( Ym(I)-Yc(I) )/Yerr(I))**2.0
!  				func_val = func_val + tmp
!  			END DO
! 		endif
!
!
! 	end subroutine ObjectiveFunction
	
end module ObjFuncMod