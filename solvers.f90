! ##############################################################################

MODULE CAVITY_FLOW
IMPLICIT NONE

CONTAINS

! ==============================================================================

SUBROUTINE SET_PRESSURE_RHS(DX, DY, DT, RHO, U, V, P, RHS, M, N)
IMPLICIT NONE

INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, RHO
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: U, V, P
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: RHS

INTEGER :: I, J
DOUBLE PRECISION :: UDX, UDY, VDX, VDY, FDX, FDY, FACTOR

FDX = 1.0D0/(2.0D0*DX)
FDY = 1.0D0/(2.0D0*DY)
FACTOR = (RHO*DX*DX*DY*DY)/(2.0D0*(DX*DX + DY*DY))

DO J = 2, N-1
    DO I = 2, M-1
        UDX = (U(I+1,J) - U(I-1,J))*FDX
        UDY = (U(I,J+1) - U(I,J-1))*FDY
        VDX = (V(I+1,J) - V(I-1,J))*FDX
        VDY = (V(I,J+1) - V(I,J-1))*FDY
	    RHS(I,J) = FACTOR*((UDX + VDY)/DT - UDX*UDX - 2.0D0*UDY*VDX - VDY*VDY)
    END DO
END DO

END SUBROUTINE SET_PRESSURE_RHS

! ==============================================================================

SUBROUTINE SOLVE_PRESSURE(DX, DY, P, RHS, M, N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: RHS
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: P

INTEGER :: I, J
DOUBLE PRECISION :: DX2, DY2, DIV, PKI, PKJ
DOUBLE PRECISION, DIMENSION(M, N) :: PK
	
DX2 = DX*DX
DY2 = DY*DY
DIV = 1.0D0/(2.0D0*(DX2 + DY2))
PK = P

DO J = 2, N-1
    DO I = 2, M-1
        PKI = (PK(I+1,J) + PK(I-1,J))*DY2
        PKJ = (PK(I,J+1) + PK(I,J-1))*DX2
        P(I,J) = (PKI + PKJ)*DIV - RHS(I,J)
    END DO
END DO
END SUBROUTINE SOLVE_PRESSURE

! ==============================================================================

SUBROUTINE SOLVE_VELOCITY(DX, DY, DT, RHO, NU, U, V, P, M, N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, RHO, NU
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: P
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: U, V

INTEGER :: I, J
DOUBLE PRECISION, DIMENSION(M, N) :: UK, VK

UK = U
VK = V

DO J = 2, N-1
    DO I = 2, M-1
        U(I,J) = UK(I,J) &
            - (UK(I,J) - UK(I-1,J))*UK(I,J)*DT/DX &
            - (UK(I,J) - UK(I,J-1))*VK(I,J)*DT/DY &
            - (P(I+1,J) - P(I-1,J))*DT/(2.0D0*RHO*DX) &
            + (UK(I+1,J) - 2.0D0*UK(I,J) + UK(I-1,J))*NU*DT/DX/DX &
            + (UK(I,J+1) - 2.0D0*UK(I,J) + UK(I,J-1))*NU*DT/DY/DY
        
        V(I,J) = VK(I,J) &
            - (VK(I,J) - VK(I-1,J))*UK(I,J)*DT/DX &
            - (VK(I,J) - VK(I,J-1))*VK(I,J)*DT/DY &
            - (P(I,J+1) - P(I,J-1))*DT/(2.0D0*RHO*DY) &
            + (VK(I+1,J) - 2.0D0*VK(I,J) + VK(I-1,J))*NU*DT/DX/DX &
            + (VK(I,J+1) - 2.0D0*VK(I,J) + VK(I,J-1))*NU*DT/DY/DY
    END DO
END DO
END SUBROUTINE SOLVE_VELOCITY

END MODULE CAVITY_FLOW

! ##############################################################################

MODULE CHANNEL_FLOW
IMPLICIT NONE

CONTAINS

! ==============================================================================

SUBROUTINE SET_PRESSURE_RHS(DX, DY, DT, RHO, U, V, P, RHS, M, N)
IMPLICIT NONE

INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, RHO
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: U, V, P
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: RHS

INTEGER :: I, J
DOUBLE PRECISION :: UDX, UDY, VDX, VDY, FDX, FDY, FACTOR

FDX = 1.0D0/(2.0D0*DX)
FDY = 1.0D0/(2.0D0*DY)
FACTOR = (RHO*DX*DX*DY*DY)/(2.0D0*(DX*DX + DY*DY))

DO J = 2, N-1
    DO I = 2, M-1
        UDX = (U(I+1,J) - U(I-1,J))*FDX
        UDY = (U(I,J+1) - U(I,J-1))*FDY
        VDX = (V(I+1,J) - V(I-1,J))*FDX
        VDY = (V(I,J+1) - V(I,J-1))*FDY
	    RHS(I,J) = FACTOR*((UDX + VDY)/DT - UDX*UDX - 2.0D0*UDY*VDX - VDY*VDY)
    END DO
END DO

DO J = 2, N-1
    
    ! PERIODIC BOUNDARY CONDITION OUTLET (I=M).
    
    UDX = (U(1,J) - U(M-1,J))*FDX
    UDY = (U(M,J+1) - U(M,J-1))*FDY
    VDX = (V(1,J) - V(M-1,J))*FDX
    VDY = (V(M,J+1) - V(M,J-1))*FDY
    RHS(M,J) = FACTOR*((UDX + VDY)/DT - UDX*UDX - 2.0D0*UDY*VDX - VDY*VDY)

    ! PERIODIC BOUNDARY CONDITION INLET (I=1).

    UDX = (U(2,J) - U(M,J))*FDX
    UDY = (U(1,J+1) - U(1,J-1))*FDY
    VDX = (V(2,J) - V(M,J))*FDX
    VDY = (V(1,J+1) - V(1,J-1))*FDY
    RHS(I,J) = FACTOR*((UDX + VDY)/DT - UDX*UDX - 2.0D0*UDY*VDX - VDY*VDY)
    
END DO

END SUBROUTINE SET_PRESSURE_RHS

! ==============================================================================

SUBROUTINE SOLVE_PRESSURE(DX, DY, P, RHS, M, N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: RHS
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: P

INTEGER :: I, J
DOUBLE PRECISION :: DX2, DY2, DIV, PKI, PKJ
DOUBLE PRECISION, DIMENSION(M, N) :: PK

DX2 = DX*DX
DY2 = DY*DY
DIV = 1.0D0/(2.0D0*(DX2 + DY2))
PK = P

DO J = 2, N-1
    DO I = 2, M-1
        PKI = (PK(I+1,J) + PK(I-1,J))*DY2
        PKJ = (PK(I,J+1) + PK(I,J-1))*DX2
        P(I,J) = (PKI + PKJ)*DIV - RHS(I,J)
    END DO
END DO

DO J = 2, N-1

    ! PERIODIC BOUNDARY CONDITION OUTLET.

    PKI = (PK(1,J) + PK(M-1,J))*DY2
    PKJ = (PK(M,J+1) + PK(M,J-1))*DX2
    P(M,J) = (PKI + PKJ)*DIV - RHS(M,J)

    ! PERIODIC BOUNDARY CONDITION INLET.

    PKI = (PK(2,J) + PK(M,J))*DY2
    PKJ = (PK(1,J+1) + PK(1,J-1))*DX2
    P(1,J) = (PKI + PKJ)*DIV - RHS(1,J)
    
END DO

DO I = 1, M
    
    ! WALL BOUNDARY CONDITION DP/DY = 0.0
    
    P(I, 1) = P(I, 2)
    P(I, N) = P(I, N-1)
    
END DO

END SUBROUTINE SOLVE_PRESSURE

! ==============================================================================

SUBROUTINE SOLVE_VELOCITY(DX, DY, DT, RHO, NU, U, V, P, F, M, N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, RHO, NU
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: P, F
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: U, V

INTEGER :: I, J
DOUBLE PRECISION, DIMENSION(M, N) :: UK, VK

UK = U
VK = V

DO J = 2, N-1
    DO I = 2, M-1

        U(I,J) = UK(I,J) &
            - (UK(I,J) - UK(I-1,J))*UK(I,J)*DT/DX &
            - (UK(I,J) - UK(I,J-1))*VK(I,J)*DT/DY &
            - (P(I+1,J) - P(I-1,J))*DT/(2.0D0*RHO*DX) &
            + (UK(I+1,J) - 2.0D0*UK(I,J) + UK(I-1,J))*NU*DT/DX/DX &
            + (UK(I,J+1) - 2.0D0*UK(I,J) + UK(I,J-1))*NU*DT/DY/DY + F(I,J)*DT
        
        V(I,J) = VK(I,J) &
            - (VK(I,J) - VK(I-1,J))*UK(I,J)*DT/DX &
            - (VK(I,J) - VK(I,J-1))*VK(I,J)*DT/DY &
            - (P(I,J+1) - P(I,J-1))*DT/(2.0D0*RHO*DY) &
            + (VK(I+1,J) - 2.0D0*VK(I,J) + VK(I-1,J))*NU*DT/DX/DX &
            + (VK(I,J+1) - 2.0D0*VK(I,J) + VK(I,J-1))*NU*DT/DY/DY

    END DO
END DO

DO J = 2, N-1
        
    ! PERIODIC BOUNDARY CONDITION OUTLET.

    U(M,J) = UK(M,J) &
        - (UK(M,J) - UK(M-1,J))*UK(M,J)*DT/DX &
        - (UK(M,J) - UK(M,J-1))*VK(M,J)*DT/DY &
        - (P(1,J) - P(M-1,J))*DT/(2.0D0*RHO*DX) &
        + (UK(M-1,J) - 2.0D0*UK(M,J) + UK(1,J))*NU*DT/DX/DX &
        + (UK(M,J+1) - 2.0D0*UK(M,J) + UK(M,J-1))*NU*DT/DY/DY + F(M,J)*DT
    
    V(M,J) = VK(M,J) &
        - (VK(M,J) - VK(M-1,J))*UK(M,J)*DT/DX &
        - (VK(M,J) - VK(M,J-1))*VK(M,J)*DT/DY &
        - (P(M,J+1) - P(M,J-1))*DT/(2.0D0*RHO*DY) &
        + (VK(M-1,J) - 2.0D0*VK(M,J) + VK(M,J))*NU*DT/DX/DX &
        + (VK(M,J+1) - 2.0D0*VK(M,J) + VK(M,J-1))*NU*DT/DY/DY

    ! PERIODIC BOUNDARY CONDITION INLET.

    U(1,J) = UK(1,J) &
        - (UK(1,J) - UK(M,J))*UK(1,J)*DT/DX &
        - (UK(1,J) - UK(1,J-1))*VK(1,J)*DT/DY &
        - (P(2,J) - P(M,J))*DT/(2.0D0*RHO*DX) &
        + (UK(2,J) - 2.0D0*UK(1,J) + UK(M,J))*NU*DT/DX/DX &
        + (UK(1,J+1) - 2.0D0*UK(1,J) + UK(1,J-1))*NU*DT/DY/DY + F(1,J)*DT
    
    V(1,J) = VK(1,J) &
        - (VK(1,J) - VK(M,J))*UK(1,J)*DT/DX &
        - (VK(1,J) - VK(1,J-1))*VK(1,J)*DT/DY &
        - (P(1,J+1) - P(1,J-1))*DT/(2.0D0*RHO*DY) &
        + (VK(2,J) - 2.0D0*VK(1,J) + VK(M,J))*NU*DT/DX/DX &
        + (VK(1,J+1) - 2.0D0*VK(1,J) + VK(1,J-1))*NU*DT/DY/DY

END DO

DO I = 1, M
    
    ! NO SLIP WALL BOUNDARY CONDITION U = V = 0.0
    
    U(I,1) = 0.0D0
    V(I,1) = 0.0D0
    
    U(I,N) = 0.0D0
    V(I,N) = 0.0D0
    
END DO

END SUBROUTINE SOLVE_VELOCITY

END MODULE CHANNEL_FLOW
