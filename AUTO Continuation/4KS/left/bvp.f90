!Free parameter=NBC+NINT-NDIM+1

!unames = {1:'u', 2:'v', 3:'r', 4:'s'}
!parnames = {1:'a', 2:'b', 3:'c', 4:'d', 5:'h', 6:'delta'}


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
		
		DOUBLE PRECISION a, b, c, D, r, r1, r2

		  a = PAR(1)
		  b = PAR(2)
		  c = PAR(3)
		  D = PAR(4)
		  r = PAR(5)
		  r1 = PAR(6)
		  r2 = PAR(7)
		  
		  F(1) = U(5)
		  F(2) = U(6)
		  F(3) = U(7)
		  F(4) = U(8)
		  F(5) = - r*U(1)*(1 - U(1)) + c*((1 - U(1)**2.0)*U(5)*U(6)/(1 + U(1)**2.0)**2.0 - (1/D)*U(1)/(1 + U(1)**2.0)* &
		  (U(1) - a*U(2) - r2*U(2)*U(3) + b*(U(1)**2.0)*U(4)))
		  F(6) = - (U(1) - a*U(2) - r2*U(2)*U(3) + b*(U(1)**2.0)*U(4))/D
		  F(7) = - (- r2*U(2)*U(3) + (r1 + b*(U(1)**2.0))*U(4) + U(1)*(1 - U(3)))/D
		  F(8) = - (r2*U(2)*U(3) - (r1 + b*(U(1)**2.0))*U(4))/D
		

	END SUBROUTINE FFFF


 
	SUBROUTINE STPNT(NDIM,U,PAR,T)


		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		COMMON /FRST/ IFRST
		DIMENSION U(NDIM),PAR(*)

		IF(IFRST.NE.1)THEN
			IFRST=1
			
			PERIOD = 55*5.183850242607483
			PAR(11) = PERIOD
			

			!parnames = {1:'a', 2:'b', 3:'c', 4:'D', 5:'r', 6:'r1', 7:'r2'}	
			
			PAR(1) = 0.9
			PAR(2) = 0.001
			PAR(3) = 4.32
			PAR(4) = 1.0
			PAR(5) = 0.1
			PAR(6) = 0.2
			PAR(7) = 0.5


		ENDIF
       
       	U(1) = 1.0
		U(2) = (PAR(2) + PAR(6))/(PAR(1)*(PAR(2) + PAR(6)) + PAR(6)*PAR(7))
		U(3) = 1.0
		U(4) = PAR(7)/(PAR(1)*(PAR(2) + PAR(6)) + PAR(6)*PAR(7))
		U(5) = 0.0
		U(6) = 0.0
		U(7) = 0.0
		U(8) = 0.0

		
    END SUBROUTINE STPNT

	
    SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)


		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC)

        PI = 4.D0*DATAN(1.D0)
		
		FB(1) = U0(5)
		FB(2) = U0(6)
		FB(3) = U0(7)
		FB(4) = U0(8)
		FB(5) = U1(5)
		FB(6) = U1(6)
		FB(7) = U1(7)
		FB(8) = U1(8)


	END SUBROUTINE BCND
	
    
	SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
	END SUBROUTINE ICND


	
    SUBROUTINE FOPT
	END SUBROUTINE FOPT

	
	SUBROUTINE PVLS
	END SUBROUTINE PVLS