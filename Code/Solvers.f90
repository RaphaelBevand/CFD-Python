SUBROUTINE STEP_01(DX, DT, U, M)
!
! SOLVE ONE DIMENSIONAL LINEAR CONVECTION EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M
DOUBLE PRECISION, INTENT(IN) :: DX, DT
DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: U

INTEGER :: I
DOUBLE PRECISION, DIMENSION(M) :: UN

UN = U

DO I = 2, M
    U(I) = UN(I) &
        - DT / DX * (UN(I) - UN(I-1))
END DO

END SUBROUTINE STEP_01

! ==============================================================================

SUBROUTINE STEP_02(DX, DT, U, M)
!
! SOLVE ONE DIMENSIONAL NONLINEAR CONVECTION EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M
DOUBLE PRECISION, INTENT(IN) :: DX, DT
DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: U

INTEGER :: I
DOUBLE PRECISION, DIMENSION(M) :: UN

UN = U

DO I = 2, M
    U(I) = UN(I) &
        - UN(I) * DT / DX * (UN(I) - UN(I-1))
END DO

END SUBROUTINE STEP_02

! ==============================================================================

SUBROUTINE STEP_03(NU, DX, DT, U, M)
!
! SOLVE ONE DIMENSIONAL DIFFUSION EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M
DOUBLE PRECISION, INTENT(IN) :: NU, DX, DT
DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: U

INTEGER :: I
DOUBLE PRECISION, DIMENSION(M) :: UN

UN = U

DO I = 2, M-1
    U(I) = UN(I) &
        + NU * DT / DX / DX * (UN(I+1) - 2.0D0*UN(I) + UN(I-1))
END DO

END SUBROUTINE STEP_03

! ==============================================================================

SUBROUTINE STEP_04(NU, DX, DT, U, M)
!
! SOLVE ONE DIMENSIONAL NONLINEAR BURGERS EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M
DOUBLE PRECISION, INTENT(IN) :: NU, DX, DT
DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: U

INTEGER :: I
DOUBLE PRECISION, DIMENSION(M) :: UN

UN = U

DO I = 2, M-1
    U(I) = UN(I) &
        - UN(I) * DT / DX * (UN(I) - UN(I-1)) &
        + NU * DT / DX / DX * (UN(I+1) - 2.0D0*UN(I) + UN(I-1))
END DO

! PERIODIC BOUNDARY CONDITION.

U(1) = UN(1) - UN(1) * DT/DX * (UN(1) - UN(M-1)) &
    + NU*DT/DX/DX * (UN(2) - 2.0D0*UN(1) + UN(M-1))
U(M) = U(1)

END SUBROUTINE STEP_04

! ==============================================================================

SUBROUTINE STEP_05(DX, DY, DT, U, M, N)
!
! SOLVE TWO DIMENSIONAL LINEAR CONVECTION EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: U

INTEGER :: I, J
DOUBLE PRECISION, DIMENSION(M, N) :: UN

DO I = 1, M
    U(I,1) = 1.0D0
    U(I,N) = 1.0D0
END DO

DO J = 1, N
    U(1,J) = 1.0D0
    U(M,J) = 1.0D0
END DO

UN = U

DO J = 2, N
    DO I = 2, M
        U(I,J) = UN(I,J) &
            - DT / DX * (UN(I,J) - UN(I-1,J)) &
            - DT / DY * (UN(I,J) - UN(I,J-1))
    END DO
END DO

END SUBROUTINE STEP_05

! ==============================================================================

SUBROUTINE STEP_06(DX, DY, DT, U, V, M, N)
!
! SOLVE TWO DIMENSIONAL NONLINEAR CONVECTION EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: U, V

INTEGER :: I, J
DOUBLE PRECISION, DIMENSION(M, N) :: UN, VN

DO I = 1, M
    U(I,1) = 1.0D0
    U(I,N) = 1.0D0
    V(I,1) = 1.0D0
    V(I,N) = 1.0D0    
END DO

DO J = 1, N
    U(1,J) = 1.0D0
    U(M,J) = 1.0D0
    V(1,J) = 1.0D0
    V(M,J) = 1.0D0
END DO

UN = U
VN = V

DO J = 2, N
    DO I = 2, M
        U(I,J) = UN(I,J) &
            - UN(I,J) * DT / DX * (UN(I,J) - UN(I-1,J)) &
            - VN(I,J) * DT / DY * (UN(I,J) - UN(I,J-1))

        V(I,J) = VN(I,J) &
            - UN(I,J) * DT / DX * (VN(I,J) - VN(I-1,J)) &
            - VN(I,J) * DT / DY * (VN(I,J) - VN(I,J-1))
    END DO
END DO

END SUBROUTINE STEP_06

! ==============================================================================

SUBROUTINE STEP_07(DX, DY, DT, NU, U, M, N)
!
! SOLVE TWO DIMENSIONAL DIFFUSION CONVECTION EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, NU
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: U

INTEGER :: I, J
DOUBLE PRECISION, DIMENSION(M, N) :: UN

DO I = 1, M
    U(I,1) = 1.0D0
    U(I,N) = 1.0D0   
END DO

DO J = 1, N
    U(1,J) = 1.0D0
    U(M,J) = 1.0D0
END DO

UN = U

DO J = 2, N-1
    DO I = 2, M-1
        U(I,J) = UN(I,J) &
            + NU * DT / DX / DX * (UN(I+1,J) - 2.0D0 * UN(I,J) + UN(I-1,J)) &
            + NU * DT / DY / DY * (UN(I,J+1) - 2.0D0 * UN(I,J) + UN(I,J-1))
    END DO
END DO

END SUBROUTINE STEP_07

! ==============================================================================

SUBROUTINE STEP_08(DX, DY, DT, NU, U, V, M, N)
!
! SOLVE TWO DIMENSIONAL NONLINEAR BURGERS EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, NU
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: U, V

INTEGER :: I, J
DOUBLE PRECISION, DIMENSION(M, N) :: UN, VN

DO I = 1, M
    U(I,1) = 1.0D0
    U(I,N) = 1.0D0
    V(I,1) = 1.0D0
    V(I,N) = 1.0D0    
END DO

DO J = 1, N
    U(1,J) = 1.0D0
    U(M,J) = 1.0D0
    V(1,J) = 1.0D0
    V(M,J) = 1.0D0
END DO

UN = U
VN = V

DO J = 2, N-1
    DO I = 2, M-1
        U(I,J) = UN(I,J) &
            - UN(I,J) * DT / DX * (UN(I,J) - UN(I-1,J)) &
            - VN(I,J) * DT / DY * (UN(I,J) - UN(I,J-1)) &
            + NU * DT / DX / DX * (UN(I+1,J) - 2.0D0 * UN(I,J) + UN(I-1,J)) &
            + NU * DT / DY / DY * (UN(I,J+1) - 2.0D0 * UN(I,J) + UN(I,J-1))

        V(I,J) = VN(I,J) &
            - UN(I,J) * DT / DX * (VN(I,J) - VN(I-1,J)) &
            - VN(I,J) * DT / DY * (VN(I,J) - VN(I,J-1)) &
            + NU * DT / DX / DX * (VN(I+1,J) - 2.0D0 * VN(I,J) + VN(I-1,J)) &
            + NU * DT / DY / DY * (VN(I,J+1) - 2.0D0 * VN(I,J) + VN(I,J-1))
    END DO
END DO

END SUBROUTINE STEP_08

! ==============================================================================

SUBROUTINE STEP_09(DX, DY, P, RESIDUAL, M, N)
!
! SOLVE TWO DIMENSIONAL LAPLACE EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY
DOUBLE PRECISION, INTENT(OUT) :: RESIDUAL
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: P

INTEGER :: I, J
DOUBLE PRECISION :: FACTOR
DOUBLE PRECISION, DIMENSION(M, N) :: PN

DO I = 1, M
    P(I,1) = P(I,2)
    P(I,N) = P(I,N-1)
END DO

DO J = 1, N
    P(1,J) = 0.0D0
    P(M,J) = (J - 1.0D0) * DY
END DO

PN = P
FACTOR = 1.0D0 / (2.0D0 * (DX * DX + DY * DY))
RESIDUAL = 0.0D0

DO J = 2, N-1
    DO I = 2, M-1
        P(I,J) =  FACTOR * ( &
            DY * DY * (PN(I+1,J) + PN(I-1,J)) + & 
            DX * DX * (PN(I,J+1) + PN(I,J-1)) )
        RESIDUAL = MAX(RESIDUAL, ABS(P(I,J) - PN(I,J)))
    END DO
END DO

END SUBROUTINE STEP_09

! ==============================================================================

SUBROUTINE STEP_10(DX, DY, RHS, P, RESIDUAL, M, N)
!
! SOLVE TWO DIMENSIONAL POISSON EQUATION.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: RHS
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: P
DOUBLE PRECISION, INTENT(OUT) :: RESIDUAL

INTEGER :: I, J
DOUBLE PRECISION :: FACTOR
DOUBLE PRECISION, DIMENSION(M, N) :: PN

DO I = 1, M
    P(I,1) = 0.0D0
    P(I,N) = 0.0D0
END DO

DO J = 1, N
    P(1,J) = 0.0D0
    P(M,J) = 0.0D0
END DO

PN = P
FACTOR = 1.0D0 / (2.0D0 * (DX * DX + DY * DY))
RESIDUAL = 0.0D0

DO J = 2, N-1
    DO I = 2, M-1
        P(I,J) =  FACTOR * ( &
            DY * DY * (PN(I+1,J) + PN(I-1,J)) + & 
            DX * DX * (PN(I,J+1) + PN(I,J-1)) - &
            DX * DX * DY * DY * RHS(I,J))
        RESIDUAL = MAX(RESIDUAL, ABS(P(I,J) - PN(I,J)))
    END DO
END DO

END SUBROUTINE STEP_10

! ==============================================================================

SUBROUTINE STEP_11(INNER, DX, DY, DT, RHO, NU, P, U, V, M, N)
!
! SOLVE THE TWO DIMENSIONAL NAVIER-STOKES CAVITY FLOW.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: INNER, M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, RHO, NU
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: P, U, V

INTEGER :: I, J, K
DOUBLE PRECISION :: UDX, UDY, VDX, VDY, FDX, FDY, FACTOR
DOUBLE PRECISION, DIMENSION(M, N) :: PN, UN, VN, RHS

DO I = 1, M
    U(I,1) = 0.0D0
    U(I,N) = 1.0D0
    V(I,1) = 0.0D0
    V(I,N) = 0.0D0    
END DO

DO J = 1, N
    U(1,J) = 0.0D0
    U(M,J) = 0.0D0
    V(1,J) = 0.0D0
    V(M,J) = 0.0D0
END DO

! SET PRESSURE POISSON EQUATION SOURCE TERM.

FDX = 1.0D0 / (2.0D0 * DX)
FDY = 1.0D0 / (2.0D0 * DY)
FACTOR = (RHO * DX * DX * DY * DY) / (2.0D0 * (DX * DX + DY * DY))

DO J = 2, N-1
    DO I = 2, M-1
        UDX = (U(I+1,J) - U(I-1,J)) * FDX
        UDY = (U(I,J+1) - U(I,J-1)) * FDY
        VDX = (V(I+1,J) - V(I-1,J)) * FDX
        VDY = (V(I,J+1) - V(I,J-1)) * FDY
        RHS(I,J) = FACTOR * ( &
            (UDX + VDY) / DT - UDX * UDX - 2.0D0 * UDY * VDX - VDY * VDY)
    END DO
END DO

FACTOR = 1.0D0 / (2.0D0 * (DX * DX + DY * DY))

DO K = 1, INNER

    ! SET BOUNDARY CONDITIONS OF PRESSURE.
    
    DO I = 1, M
        P(I,1) = P(I,2) ! y = ymin | dp/dy = 0.0
        P(I,N) = 0.0D0  ! y = ymax | p = 0.0
    END DO

    DO J = 2, N-1
        P(1,J) = P(2,J)   ! x = xmin | dp/dx = 0.0
        P(M,J) = P(M-1,J) ! x = xmax | dp/dx = 0.0
    END DO
    
    ! SOLVE PRESSURE POISSON EQUATION.
    
    PN = P

    DO J = 2, N-1
        DO I = 2, M-1
            P(I,J) = FACTOR * ( &
                (PN(I+1,J) + PN(I-1,J)) * DY * DY + &
                (PN(I,J+1) + PN(I,J-1)) * DX * DX ) - RHS(I,J)
        END DO
    END DO
END DO

! SOLVE VELOCITIES USING THE MOMENTUM EQUATIONS.

UN = U
VN = V

DO J = 2, N-1
    DO I = 2, M-1
        U(I,J) = UN(I,J) &
            - UN(I,J) * DT / DX * (UN(I,J) - UN(I-1,J)) &
            - VN(I,J) * DT / DY * (UN(I,J) - UN(I,J-1)) &
            - DT / (2.0D0 * RHO * DX) * (P(I+1,J) - P(I-1,J)) &
            + NU * DT / DX / DX * (UN(I+1,J) - 2.0D0*UN(I,J) + UN(I-1,J)) &
            + NU * DT / DY / DY * (UN(I,J+1) - 2.0D0*UN(I,J) + UN(I,J-1))
        
        V(I,J) = VN(I,J) &
            - UN(I,J) * DT / DX * (VN(I,J) - VN(I-1,J)) &
            - VN(I,J) * DT / DY * (VN(I,J) - VN(I,J-1)) &
            - DT / (2.0D0 * RHO * DY) * (P(I,J+1) - P(I,J-1)) &
            + NU * DT / DX / DX * (VN(I+1,J) - 2.0D0*VN(I,J) + VN(I-1,J)) &
            + NU * DT / DY / DY * (VN(I,J+1) - 2.0D0*VN(I,J) + VN(I,J-1))
    END DO
END DO

END SUBROUTINE STEP_11

! ==============================================================================

SUBROUTINE STEP_12(INNER, DX, DY, DT, RHO, NU, F, P, U, V, M, N)
!
! SOLVE THE TWO DIMENSIONAL NAVIER-STOKES CHANNEL FLOW.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: INNER, M, N
DOUBLE PRECISION, INTENT(IN) :: DX, DY, DT, RHO, NU
DOUBLE PRECISION, DIMENSION(M, N), INTENT(IN) :: F
DOUBLE PRECISION, DIMENSION(M, N), INTENT(INOUT) :: P, U, V

INTEGER :: I, J, K
DOUBLE PRECISION :: UDX, UDY, VDX, VDY, FDX, FDY, FACTOR
DOUBLE PRECISION, DIMENSION(M, N) :: PN, UN, VN, RHS

! SET PRESSURE POISSON EQUATION SOURCE TERM.

FDX = 1.0D0 / (2.0D0 * DX)
FDY = 1.0D0 / (2.0D0 * DY)
FACTOR = (RHO * DX * DX * DY * DY) / (2.0D0 * (DX * DX + DY * DY))

DO J = 2, N-1
    DO I = 2, M-1
        UDX = (U(I+1,J) - U(I-1,J)) * FDX
        UDY = (U(I,J+1) - U(I,J-1)) * FDY
        VDX = (V(I+1,J) - V(I-1,J)) * FDX
        VDY = (V(I,J+1) - V(I,J-1)) * FDY
        RHS(I,J) = FACTOR * ( &
            (UDX + VDY) / DT - UDX * UDX - 2.0D0 * UDY * VDX - VDY * VDY)
    END DO
END DO

DO J = 2, N-1
    ! PERIODIC BOUNDARY CONDITION OUTLET (I=M).
    
    UDX = (U(1,J) - U(M-1,J)) * FDX
    UDY = (U(M,J+1) - U(M,J-1)) * FDY
    VDX = (V(1,J) - V(M-1,J)) * FDX
    VDY = (V(M,J+1) - V(M,J-1)) * FDY
    RHS(M,J) = FACTOR * ( &
        (UDX + VDY) / DT - UDX * UDX - 2.0D0 * UDY * VDX - VDY * VDY)

    ! PERIODIC BOUNDARY CONDITION INLET (I=1).

    UDX = (U(2,J) - U(M,J)) * FDX
    UDY = (U(1,J+1) - U(1,J-1)) * FDY
    VDX = (V(2,J) - V(M,J)) * FDX
    VDY = (V(1,J+1) - V(1,J-1)) * FDY
    RHS(I,J) = FACTOR*( &
        (UDX + VDY) / DT - UDX * UDX - 2.0D0 * UDY * VDX - VDY * VDY)
END DO

FACTOR = 1.0D0 / (2.0D0 * (DX * DX + DY * DY))

DO K = 1, INNER
    
    ! SOLVE PRESSURE POISSON EQUATION.
    
    PN = P

    DO J = 2, N-1
        DO I = 2, M-1
            P(I,J) = FACTOR * ( &
                (PN(I+1,J) + PN(I-1,J)) * DY * DY + &
                (PN(I,J+1) + PN(I,J-1)) * DX * DX ) - RHS(I,J)
        END DO
    END DO

    ! WALL BOUNDARY CONDITION DP/DY = 0.0    

    DO I = 1, M
        P(I, 1) = P(I, 2)
        P(I, N) = P(I, N-1)
    END DO
    
    DO J = 2, N-1
        ! PERIODIC BOUNDARY CONDITION OUTLET.
        P(M,J) = FACTOR * ( &
            (PN(1,J) + PN(M-1,J)) * DY * DY + &
            (PN(M,J+1) + PN(M,J-1)) * DX * DX ) - RHS(M,J)

        ! PERIODIC BOUNDARY CONDITION INLET.
        P(1,J) = FACTOR * ( &
            (PN(2,J) + PN(M,J)) * DY * DY + &
            (PN(1,J+1) + PN(1,J-1)) * DX * DX ) - RHS(1,J)
    END DO
END DO

! SOLVE VELOCITIES USING THE MOMENTUM EQUATIONS.

UN = U
VN = V

DO J = 2, N-1
    DO I = 2, M-1
        U(I,J) = UN(I,J) &
            - (UN(I,J) - UN(I-1,J)) * UN(I,J) * DT / DX &
            - (UN(I,J) - UN(I,J-1)) * VN(I,J) * DT / DY &
            - (P(I+1,J) - P(I-1,J)) * DT / (2.0D0 * RHO * DX) &
            + (UN(I+1,J) - 2.0D0 * UN(I,J) + UN(I-1,J)) * NU * DT / DX / DX &
            + (UN(I,J+1) - 2.0D0 * UN(I,J) + UN(I,J-1)) * NU * DT / DY / DY &
            + F(I,J) * DT
        
        V(I,J) = VN(I,J) &
            - (VN(I,J) - VN(I-1,J)) * UN(I,J) * DT / DX &
            - (VN(I,J) - VN(I,J-1)) * VN(I,J) * DT / DY &
            - (P(I,J+1) - P(I,J-1)) * DT / (2.0D0 * RHO * DY) &
            + (VN(I+1,J) - 2.0D0 * VN(I,J) + VN(I-1,J)) * NU * DT / DX / DX &
            + (VN(I,J+1) - 2.0D0 * VN(I,J) + VN(I,J-1)) * NU * DT / DY / DY
    END DO
END DO

DO J = 2, N-1
    ! PERIODIC BOUNDARY CONDITION OUTLET.
    U(M,J) = UN(M,J) &
        - (UN(M,J) - UN(M-1,J)) * UN(M,J) * DT / DX &
        - (UN(M,J) - UN(M,J-1)) * VN(M,J) * DT / DY &
        - (P(1,J) - P(M-1,J)) * DT / (2.0D0 * RHO * DX) &
        + (UN(M-1,J) - 2.0D0 * UN(M,J) + UN(1,J)) * NU * DT / DX / DX &
        + (UN(M,J+1) - 2.0D0 * UN(M,J) + UN(M,J-1)) * NU * DT / DY / DY &
        + F(M,J) * DT
    
    V(M,J) = VN(M,J) &
        - (VN(M,J) - VN(M-1,J)) * UN(M,J) * DT / DX &
        - (VN(M,J) - VN(M,J-1)) * VN(M,J) * DT / DY &
        - (P(M,J+1) - P(M,J-1)) * DT / (2.0D0 * RHO * DY) &
        + (VN(M-1,J) - 2.0D0 * VN(M,J) + VN(M,J)) * NU * DT / DX / DX &
        + (VN(M,J+1) - 2.0D0 * VN(M,J) + VN(M,J-1)) * NU * DT / DY / DY

    ! PERIODIC BOUNDARY CONDITION INLET.
    U(1,J) = UN(1,J) &
        - (UN(1,J) - UN(M,J)) * UN(1,J) * DT / DX &
        - (UN(1,J) - UN(1,J-1)) * VN(1,J) * DT / DY &
        - (P(2,J) - P(M,J)) * DT / (2.0D0 * RHO * DX) &
        + (UN(2,J) - 2.0D0 * UN(1,J) + UN(M,J)) * NU * DT / DX / DX &
        + (UN(1,J+1) - 2.0D0 * UN(1,J) + UN(1,J-1)) * NU * DT / DY / DY &
        + F(1,J) * DT
    
    V(1,J) = VN(1,J) &
        - (VN(1,J) - VN(M,J)) * UN(1,J) * DT / DX &
        - (VN(1,J) - VN(1,J-1)) * VN(1,J) * DT / DY &
        - (P(1,J+1) - P(1,J-1)) * DT / (2.0D0 * RHO * DY) &
        + (VN(2,J) - 2.0D0 * VN(1,J) + VN(M,J)) * NU * DT / DX / DX &
        + (VN(1,J+1) - 2.0D0 * VN(1,J) + VN(1,J-1)) * NU * DT / DY / DY
END DO

! NO SLIP WALL BOUNDARY CONDITION U = V = 0.0

DO I = 1, M
    U(I,1) = 0.0D0
    V(I,1) = 0.0D0
    U(I,N) = 0.0D0
    V(I,N) = 0.0D0
END DO

END SUBROUTINE STEP_12
