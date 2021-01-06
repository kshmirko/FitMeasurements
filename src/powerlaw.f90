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
	