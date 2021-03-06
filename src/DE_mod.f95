!--------------------------------
!..   Fortran .f95 MODULE      ..
!--------------------------------
! MODULE for DE_main.f95
! License: https://github.com/ian-mmm/differential-evolution_f95/blob/master/LICENSE

MODULE DE_mod !% % % % % % % % %

IMPLICIT NONE
SAVE

!==GLOBAL-VARS==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Precision Parameters: - - - - - - - - - - - - - - - - - - - - - - - - - -
    INTEGER, PARAMETER                          :: dp = SELECTED_REAL_KIND(12, 100)
    INTEGER, PARAMETER                          :: sp = SELECTED_REAL_KIND(6, 35)
                                                    !(P,R) where P=precision and R=decimal exponent range
                                                    ! -- gfortran: dp=(15,307), sp=(6,37)
CONTAINS

!!======================================================

SUBROUTINE Diff_Evol( nop, NP, LB, UB, T, F_lo, F_hi, Cr, PD, BH_best, F_best, I_prn, Obj_Ftn, thres )
! Optimizes `functn' using a differential evolution optimization technique
! -- requires global vars: dp
! -- calls ftns/subroutines: Obj_Ftn
	use mo_DLS, only : RRR, AR, SD, RN, RK
	use mo_par_DLS, only : KN1par
	INTEGER,                     INTENT(IN)     :: nop, I_prn
	INTEGER,                     INTENT(IN)     :: NP
	REAL(dp),                    INTENT(IN)     :: LB(:)
	REAL(dp),                    INTENT(IN)     :: UB(:)
	INTEGER,                     INTENT(IN)     :: T, PD
	REAL(dp),                    INTENT(IN)     :: F_lo, F_hi, Cr
	REAL(dp),                    INTENT(OUT)    :: BH_best(:)
	REAL(dp),                    INTENT(OUT)    :: F_best
	REAL(dp),										 INTENT(IN)			:: thres

	! Local variables:
	INTEGER                                     :: nn, kk, kc, tt, ik
	INTEGER                                     :: IND
	INTEGER,  DIMENSION(1)                      :: f1_pos
	INTEGER,  DIMENSION(3)                      :: spouse
	REAL(dp), DIMENSION( nop, NP )              :: THETA
	REAL(dp), DIMENSION( NP )                   :: F_theta
	REAL(dp), DIMENSION( nop )                  :: range, Z0, theta_new, theta_prime
	REAL(dp), DIMENSION(3)                      :: Z2
	REAL(dp)                                    :: fval, f1, ftri, Z1, RNP, Fdither, &
	                                                RNoP, f1_turn
	real(dp), allocatable :: Sols(:, :)
	REAL(dp)																		:: Rtmp(KN1par), SDtmp(KN1par), RNtmp, RKtmp
	LOGICAL 																		:: bfirst
!	external Obj_Ftn
	 interface
 		subroutine Obj_Ftn(NP, X, func_val)
 			integer, intent(IN)	::	NP
 			real, intent(IN) 		:: X(NP)
 			real, intent(inout) :: func_val
 		end subroutine
	end interface
!	EXTERNAL Obj_Ftn
!--------------------
!   ARGUMENTS:-
!       nop     = (Intent IN) number of total parameters
!       NP      = (Intent IN) number of points in parameter grid
!       LB,UB   = (Intent IN) vector of [lower,upper] bounds for prarameters (beta's)
!       T       = (Intent IN) number of generations
!       PD      = (Intent IN) Indicator for pertubations: 0== off, 1== on
!       SMTH    = (Intent IN) Indicator for using smooth MSE: 0== off, 1== on
!   F_lo, F_hi  = (Intent IN) lower and upper bounds for F dither, "scale factor",
!       Cr      = (Intent IN) crossover value, for binomial where Cr= prob. of 
!									spouse/mutant gene
!       eta     = (Intent IN) minimum theshold for ftn value to keep theta vector
!       theta_prime     = "mutant vector" or "spouse vector"
!       theta_new       = "trial vector" or "offspring vector"
!   NOTES:- this diff evol method can be classified as DE/rand/1/bin


	allocate(Sols(nop+1, T))
	Sols = 0.0
	RNP     = REAL( NP, dp ) ! convert NP to real
	RNoP    = REAL( nop, dp )
	f1_turn = 0.0D0
	THETA   = 0.0D0

! Create initial candidate solution space:
! -- create unchanging variables beforehand outside NP loop
  range = UB - LB
  tt = 0 ! intial grid is generation zero

	DO nn = 1, NP ! -  - @ - - Grid Loop for initial generation - - @ - -
		! -- parameters uniformly drawn between given bounds
		CALL RANDOM_NUMBER( Z0 )
		THETA( :, nn ) = Z0 * range +  LB

		! Calculate the function value:
		CALL Obj_Ftn( nop, THETA(:, nn ), fval)
		! -- SUBROUTINE Obj_Ftn( N, BH, val)
		F_theta( nn )   = fval
	END DO ! -  - @ - - end grid Loop for initial gen - - @ - -

	F_best  = MINVAL( F_theta(:) )

	WRITE(*,'("SOLUTION SET:")')
	WRITE(*,'("=============",/)')
	gen_do: DO tt = 1, T !- - generation loop - - - - - - - - - - - - - - - - - - - -
		CALL RANDOM_NUMBER( Z1 )
		Fdither = Z1 * (F_hi - F_lo) + F_lo

		np_do: DO nn =1, NP !~ ~ ~ ~ ~ ~ NP loop ~ ~ ~ ~ ~ ~
			IND = 0
			! Need to create a SPOUSE for nn:
			! -- sample from {1,...,NP} without replacement (already removing nn)
			CALL RANDOM_NUMBER( Z2 )
			spouse = FLOOR( Z2 * RNP + 1.0D0)
			DO WHILE (( spouse(1) == nn ))
				CALL RANDOM_NUMBER( Z1 )
				spouse(1) = FLOOR( Z1 * RNP + 1.0D0)
			END DO
			DO WHILE (( spouse(2) == nn ) .AND. ( spouse(2) == spouse(1)))
				CALL RANDOM_NUMBER( Z1 )
				spouse(2) = FLOOR( Z1 * RNP + 1.0D0)
			END DO
			DO WHILE (( spouse(3) == nn ) .AND. ( spouse(3) == spouse(1)) .AND. ( spouse(3) == spouse(2)))
				CALL RANDOM_NUMBER( Z1 )
				spouse(3) = FLOOR( Z1 * RNP + 1.0D0)
			END DO
			! Standard spouse creation,
			theta_prime(1: nop ) = THETA(1: nop, spouse(1)) + Fdither *&
														&( THETA(1: nop, spouse(2)) - THETA(1: nop, spouse(3)) )

			! Gene transferring:
			CALL RANDOM_NUMBER( Z0 )
			! randomly choose first gene to have from spouse/mutatant
			kc = FLOOR( Z0( nop ) * ( nop + 1.0D0) + 1.0D0 ) 
			DO kk = 1, nop
				IF ( Z0( kk ) <= Cr .OR. kc == kk ) THEN
					theta_new( kk ) = theta_prime( kk )
				ELSE
					theta_new( kk ) = THETA( kk, nn )
				END IF
			END DO

			! Calculate the function value:
			CALL Obj_Ftn( nop, theta_new(:), ftri)
			! -- SUBROUTINE Obj_Ftn( N, BH, val)

			! Picking the Winner:
			IF ( ftri < F_theta( nn ) ) THEN
				THETA( :, nn ) = theta_new
				F_theta( nn ) = ftri
				IND = 1 ! record replacements
			END IF

			! Recording stats on all vectors,
			f1_turn = f1_turn + REAL(IND,dp)

		END DO np_do !~ ~ ~ ~ ~ ~ end NP loop ~ ~ ~ ~ ~ ~

		! Calculate averages and grid stats,
		f1              = MINVAL( F_theta ) ! ftn value of best theta for each generation
		f1_pos          = MINLOC( F_theta ) ! location of best vector
		f1_turn         = f1_turn / RNP

		! Live reporting
		if (mod(tt, I_prn) == 0) then
			!PRINT*, tt, f1, f1_pos, f1_turn
		endif

		IF ( f1 < F_best ) THEN
			BH_best(:)  = THETA(:, f1_pos(1))
			F_best  = f1
			! Print best function parameters at each generation
			
		END IF
		
		Sols(1,  tt)  = F_best
		Sols(2:, tt) = BH_best
		!WRITE(*,'(I5, E13.3, $)') tt, F_best
		!WRITE(*,*) (BH_best(ik),ik=1,nop) 
		
		
		! = = = = GRID LOCK EXIT = = = =
		IF( tt> 250 ) THEN
			IF (f1_turn <= 0.0000000000000001D0 .AND. PD==0) THEN
				EXIT
			END IF
		END IF

		IF (PD ==1) THEN ! = = = PERTURBATIONS = = = =
			kk = NINT( REAL(tt,dp) / 10.0D0) * 10
			IF ( tt == kk)  THEN
				IF ( f1_turn < 0.010D0 ) THEN
					PRINT*, "Pertubations! -- -- -- -- -- --"
					DO nn =1, NP !---
						IF ( nn /= f1_pos(1) ) THEN ! don't throw away the best
							CALL RANDOM_NUMBER( Z1 )
							IF ( Z1 > 0.1D0 ) THEN ! chance to keep some bad ones
								CALL RANDOM_NUMBER( Z0 )
								THETA( :, nn ) = Z0 * range + LB
								
								CALL Obj_Ftn( nop, THETA(:, nn ), fval)
	    					! -- SUBROUTINE Obj_Ftn( N, BH, val)
								F_theta( nn ) = fval
							END IF
						END IF
					END DO !---------
				END IF
			END IF
		END IF

	END DO gen_do 
!- - - - - - - - end generation loop  - - - - - - - - - - - - - - - - - 
	tt = 0
	bfirst = .true.
	do nn=1, T
		if (Sols(1, nn) .lt. thres) then
			CALL Obj_Ftn(nop, Sols(2:,nn), fval)
			if(bfirst .eqv. .true.) then
				SDtmp(1:KN1par) = SD(1:KN1par)
				RNtmp = Sols(nop, nn)
				RKtmp = Sols(nop+1, nn)
				bfirst = .false.
			else
				SDtmp(1:KN1par) = SDtmp(1:KN1par) + SD(1:KN1par)
				RNtmp = RNtmp + Sols(nop, nn)
				RKtmp = RKtmp + Sols(nop+1, nn)
			end if
			tt = tt+1
		endif
	end do
	
	SD = SDtmp / REAL(tt) 
	AR = SD
	RN = RNtmp / REAL(tt)
	RK = RKtmp / REAL(tt)
	print '("MEAN REF.IDX =", F7.2, "+", F7.5, "i")', RN, RK
	deallocate(Sols)
END SUBROUTINE Diff_Evol
!======================================================

END MODULE DE_mod !% % % % % % % % %