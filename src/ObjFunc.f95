module ObjFuncMod
	! при работе с модулями в фортране важно спечиально инициализировать все 
	! его переменные, особенно в случае масивов
	implicit none
	public
	integer, parameter	:: WavelengthCount   = 3
	integer, parameter	:: SearchParamsCount = 5
	integer, parameter	:: InpParamsCount    = 5
	integer, parameter	:: DiscrKindRMS			 = 0
	integer, parameter	:: DiscrKindMAXABS   = 1
	

	real	:: Ym(SearchParamsCount)=(/0.0,0.0,0.0,0.0,0.0/),	 Yc(SearchParamsCount)
	real	:: LoParamVal(1:SearchParamsCount), UpParamVal(1:SearchParamsCount)
	real  :: Bsc(1:WavelengthCount), Ext(1:WavelengthCount), Wvl(1:WavelengthCount), &
					&Absb(1:WavelengthCount), Ldr(1:WavelengthCount), X(1:SearchParamsCount),&
					&Blr(1:WavelengthCount)
	real	::	Rmin, Rmax
	Integer:: discrKind
contains
	
	! Специальная функция инициализации переменных модуля
	subroutine InitObjFuncVariables()
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
		discrKind = 0
	end subroutine InitObjFuncVariables
	
	subroutine ObjectiveFunction(NP, X, func_val)
		use mo_DLS
		integer, intent(IN)	::	NP
		real, intent(IN) 		:: X(NP)
		real, intent(inout) :: func_val
		real tmp
		logical b_print
		
		integer I
		
		DO I=1, SearchParamsCount
			if ( ( X(I).lt.LoParamVal(I) ) .or.  ( X(I).gt.UpParamVal(I) ) ) then
				func_val = 99999.0
				return
			end if
		END DO
		RN = X(4)
		RK = X(5)
		
		call SIZEDIS2(-KN,1,X(1),X(2),X(3),RMIN,RMAX,RRR,AR,AC) 
		
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
		
		DO I=1, SearchParamsCount-WavelengthCount
			Yc(I+WavelengthCount) = Ext(I)
		END DO
		
		! В зависимости от того, какой свособ расчета невязки выбран
		! роизводим ычисления
		if(discrKind==DiscrKindRMS) then
			!случай невязки в виде root mean square 
			func_val = 0.0
			
			! Переписали код короче
			func_val = sum( ((Ym - Yc)/Ym)**2 )
			!DO I=1, SearchParamsCount
			!	func_val = func_val + (( Ym(I)-Yc(I) )/Ym(I))**2.0
			!END DO
			
			func_val = SQRT(func_val/SearchParamsCount)
			
		else if (discrKind==DiscrKindMAXABS) then
			! случай максимального абсолютного отклонения
			func_val=0.0
			
			! Перевисали код короче
			func_val = MAXVAL( (Ym-Yc)/Ym )
! 			DO I=1, SearchParamsCount
! 				tmp = ABS(( Ym(I)-Yc(I) )/Ym(I))
! 				if (func_val .gt. tmp) then
! 					func_val = tmp
! 				end if
! 			END DO
		
		endif
		
		! переводим в проценты полученную величину
		func_val = func_val * 100
	end subroutine ObjectiveFunction
	
	subroutine powerlaw(KN,ID,A,GAMMA,RMIN,RMAX,RRR,AR,AC,KNpar)
		implicit none
		integer, intent(IN) :: KN, ID, KNpar
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
		print *, AC, NORMAL
	end subroutine powerlaw
	
end module ObjFuncMod