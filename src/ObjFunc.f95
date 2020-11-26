module ObjFuncMod
	! при работе с модулями в фортране важно спечиально инициализировать все 
	! его переменные, особенно в случае масивов
	implicit none
	public
	integer, parameter	:: MaxInputVectorsCount = 100
	integer, parameter	:: MaxWavelengthsCount = 3
	
	
	integer, parameter	:: MaxSearchParamsCount = 10
	integer, parameter	:: MaxInpParamsCount    = 10
	
	integer, parameter	:: DiscrKindRMS			 = 0
	integer, parameter	:: DiscrKindMAXABS   = 1
	integer, parameter	:: FunctLogNormal		 = 0
	integer, parameter	:: FunctPowerLaw		 = 1
	integer, parameter	:: LnParamsCount		 = 5
	integer, parameter	:: PowParamsCount		 = 4
	  
	
	real, dimension(MaxInpParamsCount, MaxInputVectorsCount)	:: datum
	real, dimension(MaxInpParamsCount)		:: Ym,	 Yc
	real, dimension(MaxSearchParamsCount)	:: LoParamVal, UpParamVal, X
	real, dimension(MaxWavelengthsCount)		:: Bsc, Ext, Wvl, Absb, Ldr, Blr
	real	::	Rmin, Rmax
	Integer:: discrKind, func_type
	integer	:: AlphaInpParamsCount, DepolInpParamsCount, InputVectorsCount
	integer :: WavelengthCount 
	integer	:: InpParamsCount, SearchParamsCount, NGEN
contains
	
	! Специальная функция инициализации переменных модуля
	subroutine InitObjFuncVariables()
		NGEN = 1000
		datum = 0.0
		Ym = 0.0
		Yc = 0.0
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
	
	subroutine ObjectiveFunction(NP, X, func_val)
		use mo_DLS
		implicit none
		integer, intent(IN)	::	NP
		real, intent(IN) 		::	X(NP)
		real, intent(inout) ::	func_val
		real tmp
		logical b_print
		
		integer I
		
		! Chreck bounds
		DO I=1, NP
			if ( ( X(I).lt.LoParamVal(I) ) .or.  ( X(I).gt.UpParamVal(I) ) ) then
				func_val = 99999.0
				return
			end if
		END DO
		
		! Last two items in X are always RN and RK
		RN = X(NP-1)
		RK = X(NP)
		
		if ( (func_type .eq. FunctLogNormal) .and. (NP .eq. LnParamsCount) ) then
			call SIZEDIS2(-KN,1,X(1),X(2),X(3),RMIN,RMAX,RRR,AR,AC) 
		else if ( (func_type .eq. FunctPowerLaw) .and. (NP .eq. PowParamsCount) ) then
			call powerlaw(KN,X(1),X(2),RMIN,RMAX,RRR,AR,AC, KNpar)
		else
			print *, "WRONG PARAMETERS COUNT OF FITTING FUNCTION TYPE"
			STOP
		endif		
		
		SD(1:KN) = AR(1:KN)
		
		
		DO I=1, WavelengthCount
			IF (NDP==0) then
				print *, "LOADING DATABASE OF SPHEROIDS..."
				b_print = .true.
			endif
			
			WL = WVL(I)
			
			CALL OPTCHAR(NDP)
			
			if(b_print .eqv. .true.) then
				print *, "DONE."
				b_print = .false.
			endif
			
			
			Ext(I) = xext
			Absb(I) = xabs
			Bsc(I) = xext / xblr
			Ldr(I) = xldr
			Yc(I) = Bsc(I)
			Blr(I) = xblr
			
		END DO
		
		DO I=1, AlphaInpParamsCount
			Yc(I+WavelengthCount) = Ext(I)
		END DO
		
		! В зависимости от того, какой свособ расчета невязки выбран
		! роизводим ычисления
		if(discrKind==DiscrKindRMS) then
			!случай невязки в виде root mean square 
			func_val = 0.0
			
			! Переписали код короче
			!func_val = sum( ((Ym - Yc)/Ym)**2 )
			DO I=1, SearchParamsCount
				func_val = func_val + (( Ym(I)-Yc(I) )/Ym(I))**2.0
			END DO
			
			func_val = SQRT(func_val/SearchParamsCount)
			
		else if (discrKind==DiscrKindMAXABS) then
			! случай максимального абсолютного отклонения
			func_val=0.0
			
			! Перевисали код короче
			!func_val = MAXVAL( (Ym-Yc)/Ym )
 			DO I=1, SearchParamsCount
 				tmp = ABS(( Ym(I)-Yc(I) )/Ym(I))
 				if (func_val .gt. tmp) then
 					func_val = tmp
 				end if
 			END DO
		
		endif
		
		! переводим в проценты полученную величину
		func_val = func_val * 100
	end subroutine ObjectiveFunction
	
	subroutine powerlaw(KN,A,GAMMA,RMIN,RMAX,RRR,AR,AC,KNPar)
		implicit none
		integer, intent(IN) :: KN, KNpar
		real, intent(IN)		:: A, GAMMA, RMIN, RMAX
		real, intent(INOUT)	:: RRR(KNpar), AR(KNpar), AC
		integer							:: KNN, I
		REAL								:: DR, NORMAL
		
		KNN = KN
		
		DR = (LOG(RMAX)-LOG(RMIN))/REAL(KN-1)
		NORMAL = 4.188*(RMAX**(GAMMA+1.0+3.0) - RMIN**(GAMMA+1.0+3.0))/(GAMMA+1.0+3.0)
		DO I=1, KN
			RRR(I)=EXP(LOG(RMIN)+(I-1)*DR)
			AR(I) = A/NORMAL*RRR(I)**gamma
		END DO
		AR = 4.188*(RRR**3)*AR
		AC = 0.0
		DO I=2, KN
			AC = AC + 0.5*(AR(I-1)+AR(I))*(RRR(I)-RRR(I-1))
		END DO
	end subroutine powerlaw
	
end module ObjFuncMod