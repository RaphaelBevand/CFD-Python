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

