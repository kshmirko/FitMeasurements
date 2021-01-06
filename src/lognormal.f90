!
!  lognormal.f90
!  FortranFindSolution
!
!  Created by Константин Шмирко on 2020-12-30.
!  Copyright 2020 Константин Шмирко. All rights reserved.
!
! вычисление значений логнормального распределения 
! KN	- число отсчетов на распределение
! V		- суммарный объем частиц
! S		- полуширина распределения
! RM	- медианный радиус
! Важно: мы ищем решение в виде объемного распрелеления, поэтому 
! параметры RM - будут соответствовать медиане объемного распределения


subroutine lognormal(KN,V,S,RM,RMIN,RMAX,RRR,AR,AC,KNPar)
	implicit none
	integer, intent(in)		::	KN, KNPar
	real, 	 intent(in)		::	A, V, S, RM, RMIN, RMAX
	real,		 intent(inout)::	RRR(KNPar), AR(KNPar), AC
	
end subroutine lognormal