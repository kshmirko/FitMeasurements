subroutine ObjectiveFunction(NP, AX, func_val)
	use mo_DLS
	use ObjFuncMod, only : LoParamVal, UpParamVal, FunctLogNormal, LnParamsCount,&
	&FunctPowerLaw, PowParamsCount, Wvl, Ext, Bsc, Absb, Ldr, Blr, Ym, &
	&Yc, Yerr, FunctManualInput, AlphaInpParamsCount, SearchParamsCount, &
	&DiscrKindRMS, DiscrKindMAXABS, discrKind, func_type, Rmin, Rmax, &
	&DiscrKindChiSqr, WavelengthCount
	
	implicit none
	integer, intent(INOUT)	::	NP
	real, intent(INOUT) 		::	AX(NP)
	real, intent(INOUT) 		::	func_val
	real tmp
	logical b_print
	
	integer I
	
	! Chreck bounds
	DO I=1, NP
		if ( ( AX(I).lt.LoParamVal(I) ) .or.  ( AX(I).gt.UpParamVal(I) ) ) then
			func_val = 99999.0
			return
		end if
	END DO
	
	
	
	if ( (func_type .eq. FunctLogNormal) .and. (NP .eq. LnParamsCount) ) then
		call SIZEDIS2(-KN,1,AX(1),AX(2),AX(3),RMIN,RMAX,RRR,AR,AC) 
		! Last two items in X are always RN and RK
		RN = AX(NP-1)
		RK = AX(NP)
		SD(1:KN) = AR(1:KN)
	else if ( (func_type .eq. FunctPowerLaw) .and. (NP .eq. PowParamsCount) ) then
		call powerlaw(KN,AX(1),AX(2),RMIN,RMAX,RRR,AR,AC, KNpar)
		! Last two items in X are always RN and RK
		RN = AX(NP-1)
		RK = AX(NP)
		SD(1:KN) = AR(1:KN)
	else if ( func_type .eq. FunctManualInput ) then
		print *, "USER FUNC"
	else
		print *, "WRONG PARAMETERS COUNT OF FITTING FUNCTION TYPE"
		STOP
	endif		
	
	
	
	
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
	! производим вычисления
	if(discrKind==DiscrKindRMS) then
		!случай невязки в виде root mean square 
		func_val = 0.0
		
		! Переписали код короче
		!func_val = sum( ((Ym - Yc)/Ym)**2 )
		DO I=1, SearchParamsCount
			func_val = func_val + (( Ym(I)-Yc(I) )/Ym(I))**2.0
		END DO
		
		func_val = SQRT(func_val/SearchParamsCount)
		! переводим в проценты полученную величину
		func_val = func_val * 100
	else if (discrKind==DiscrKindMAXABS) then
		! случай максимального абсолютного отклонения
		func_val=0.0
		
		! Переписали код короче
		!func_val = MAXVAL( (Ym-Yc)/Ym )
		DO I=1, SearchParamsCount
			tmp = ABS(( Ym(I)-Yc(I) )/Ym(I))
			if (func_val .gt. tmp) then
				func_val = tmp
			end if
		END DO
		! переводим в проценты полученную величину
		func_val = func_val * 100
	else if (discrKind==DiscrKindChiSqr) then
		! случай максимального абсолютного отклонения
		func_val=0.0
	
		DO I=1, SearchParamsCount
			tmp = (( Ym(I)-Yc(I) )/Yerr(I))**2.0
			func_val = func_val + tmp
		END DO
	endif
	
	!print '( F10.4)', func_val
end subroutine ObjectiveFunction


subroutine fmin(result, np, ax, grad, need_gradient, f_data)
	implicit none
	integer need_gradient, np, f_data
	real result, ax(np), grad(np)
	
	Call ObjectiveFunction(np, ax, result)
end subroutine fmin