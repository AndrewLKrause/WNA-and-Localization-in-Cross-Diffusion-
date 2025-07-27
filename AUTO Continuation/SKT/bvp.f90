!Free parameter=NBC+NINT-NDIM+1

!unames = {1:'u', 2:'v', 3:'r', 4:'s'}
!parnames = {1:'alpha', 2:'beta', 3:'gamma', 4:'m', 5:'D1', 6:'D2'}


    SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)

		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		DIMENSION U(NDIM),PAR(*),F(NDIM)

		CALL FFFF(NDIM,U,ICP,PAR,IJAC,F)
		PERIOD=PAR(11)
        DO  I=1,NDIM
			F(I)=PERIOD*F(I)
		END DO

    END SUBROUTINE FUNC

	
	
    SUBROUTINE FFFF(NDM,U,ICP,PAR,IJAC,F)


		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		DIMENSION U(NDM),PAR(*),F(NDM),DFDU(NDM,NDM),DFDP(NDM,NDM)
		
		DOUBLE PRECISION r1, r2, a1, p, b1, b2, d1, d2, d11, d12, q, d22

		  r1 = PAR(1)
		  r2 = PAR(2)
		  a1 = PAR(3)
		  p = PAR(4)
		  b1 = PAR(5)
		  b2 = PAR(6)
		  d1 = PAR(7)
		  d2 = PAR(8)
		  d11 = PAR(9)
		  d12 = PAR(10)
		  q = PAR(12)
		  d22 = PAR(13)
		  
		  F(1) = U(3)
		  F(2) = U(4)
	!	  F(3) = - (((q*U(1) + d2 + U(2)*d22)*(U(1)*(- 1 + U(1)*a1 + U(2)*b1)*r1 - U(3)*(d11*U(3) + 2*d12*U(4))) &
	!	  			+ U(1)*d12*(- U(2)*(- 1 + p*U(2) + U(1)*b2)*r2 + U(4)*(2*q*U(3) + d22*U(4))))/(q*U(1)*U(2)*d12 &
	!				- (d1 + U(1)*d11 + U(2)*d12)*(q*U(1) + d2 + U(2)*d22)))
	!	  F(4) = (- q*U(1)*U(2)*(- 1 + U(1)*a1 + U(2)*b1)*r1 + U(2)*(- 1 + p*U(2) + U(1)*b2)*(d1 + U(1)*d11 + U(2)*d12)*r2 &
	!	  			+ q*U(2)*d11*U(3)**2.0 - 2*q*(d1 + U(1)*d11)*U(3)*U(4) - (d1 + U(1)*d11 + U(2)*d12)*d22*U(4)**2.0)/ &
	!				(U(2)*d12*(d2 + U(2)*d22) + d1*(q*U(1) + d2 + U(2)*d22) + U(1)*d11*(q*U(1) + d2 + U(2)*d22))
		  F(3) = U(1)*(- 1 + U(1)*a1 + U(2)*b1)*r1/d1
		  F(4) = (-q*U(1)*U(2)*(- 1 + U(1)*a1 + U(2)*b1)*r1 + d1*(U(2)*(- 1 + p*U(2) + U(1)*b2)*r2 - 2*q*U(3)*U(4)))/(d1*(q*U(1) + d2))

	END SUBROUTINE FFFF


 
	SUBROUTINE STPNT(NDIM,U,PAR,T)


		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		COMMON /FRST/ IFRST
		DIMENSION U(NDIM),PAR(*)

		IF(IFRST.NE.1)THEN
			IFRST = 1
			
			PERIOD = 4*8.224672274330526
			PAR(11)= PERIOD
			

			!parnames = {1:'r1', 2:'r2', 3:'a1', 4:'p', 5:'b1', 6:'b2', 7:'d1', 8:'d2', 9:'d11', 10: 'd12', 12: 'q', 13: 'd22'}
			
			PAR(1) = 1.0
      		PAR(2) = 2.0
      		PAR(3) = 0.9
			PAR(4) = 0.5
			PAR(5) = 0.6
			PAR(6) = 0.2
			PAR(7) = 1.0
			PAR(8) = 1.0
			PAR(9) = 0.0
			PAR(10) = 0.0
			PAR(12) = 52.86845318178516
			PAR(13) = 0.0


		ENDIF

		U(1) = (PAR(4) - PAR(5))/(PAR(4)*PAR(3) - PAR(5)*PAR(6))
		U(2) = (PAR(3) - PAR(6))/(PAR(4)*PAR(3) - PAR(5)*PAR(6))
		U(3) = 0.0
		U(4) = 0.0

		
    END SUBROUTINE STPNT

	
    SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)


		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC)

        PI = 4.D0*DATAN(1.D0)
		
		FB(1) = U0(3)
		FB(2) = U0(4)
		FB(3) = U1(3)
		FB(4) = U1(4)


	END SUBROUTINE BCND
	
    
	SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
		!     ---------- ----
	
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
		DOUBLE PRECISION, INTENT(IN) :: PAR(*)
		DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
		DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
		DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
	
		FI(1) = UPOLD(1)*(U(1)-UOLD(1)) ! phase condition
	
	END SUBROUTINE ICND


	
    SUBROUTINE FOPT
	END SUBROUTINE FOPT

	
	SUBROUTINE PVLS
	END SUBROUTINE PVLS