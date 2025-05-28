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
		
		DOUBLE PRECISION sigma, lambda, d, p

		  sigma = PAR(1)
		  lambda = PAR(2)
		  d = PAR(3)
		  p = PAR(4)
		  
		  F(1) = U(3)
		  F(2) = U(4)
		  F(3) = U(1) - U(1)**2.0*U(2) - sigma*(U(1) - 1/U(2))**2.0 - p*(U(1)**2*U(2) + sigma*(U(1) - 1/U(2))**2.0 - lambda)/d
		  F(4) = (U(1)**2*U(2) + sigma*(U(1) - 1/U(2))**2.0 - lambda)/d
		

	END SUBROUTINE FFFF


 
	SUBROUTINE STPNT(NDIM,U,PAR,T)


		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		COMMON /FRST/ IFRST
		DIMENSION U(NDIM),PAR(*)

		IF(IFRST.NE.1)THEN
			IFRST=1
			
			PERIOD = 51*4.881324902151784
			PAR(11)= PERIOD
			

			!parnames = {1:'sigma', 2:'lambda', 3:'d', 4:'p'}	
			
			PAR(1) = - 1.5
      		PAR(2) = 1.5
      		PAR(3) = 5.82842712474619
			PAR(4) = 0.0


		ENDIF
       
       	U(1) = PAR(2)
		U(2) = 1/par(2)
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
	END SUBROUTINE ICND


	
    SUBROUTINE FOPT
	END SUBROUTINE FOPT

	
	SUBROUTINE PVLS
	END SUBROUTINE PVLS