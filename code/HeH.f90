!****************************************************************************
!
!  PROGRAM: HeHSCF
!
!****************************************************************************

program HeHSCF
    implicit none
    real*8 :: IOP, R, ZETA1, ZETA2, ZA, ZB
    integer :: N
    IOP=3
    N=3
    R=1.4632D0
    ZETA1=2.0925D0
    ZETA2=1.24D0
    ZA=2.0D0
    ZB=1.0D0
    CALL HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
end program HeHSCF

subroutine HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
    implicit none
    real*8 :: IOP, R, ZETA1, ZETA2, ZA, ZB
    integer :: N
    if (IOP /= 0) then
        write(*, 10) N, ZA, ZB
    endif
    call INTGRL(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
    call COLECT(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
    call SCF(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
10  format(' ',2X,'STO-',I1,'G FOR ATOMIC NUMBERS ',F5.2,' AND ',F5.2//) 
    return
end subroutine HFCALC
    
subroutine INTGRL(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
    implicit none
    real*8 :: IOP, R, ZETA1, ZETA2, ZA, ZB
    integer :: N
    real*8 :: R2,RAP,RBP,RAQ,RBQ,RPQ,RAP2,RBP2,RAQ2,RBQ2,RPQ2
    real*8 :: S12,T11,T12,T22,V11A,V12A,V22A,V11B,V12B,V22B,V1111,V2111,V2121,V2211,V2221,V2222
    common /INT/ S12,T11,T12,T22,V11A,V12A,V22A,V11B,V12B,V22B,V1111,V2111,V2121,V2211,V2221,V2222
    real*8 :: COEF(3,3),EXPON(3,3),D1(3),A1(3),D2(3),A2(3)
    real*8 :: PI=3.1415926535898D0
    integer I,J,K,L
    real*8, external :: S,T,V,TWOE
    DATA COEF,EXPON /1.0D0,2*0.0D0,0.678914D0,0.430129D0,0.0D0,0.444635D0,0.535328D0,0.154329D0,&
        0.270950D0,2*0.0D0,0.151623D0,0.851819D0,0.0D0,0.109818D0,0.405771D0,2.22766D0/
    R2=R*R
    do I=1,N
        A1(I)=EXPON(I,N)*(ZETA1**2)
        D1(I)=COEF(I,N)*((2.0D0*A1(I)/PI)**0.75D0)
        A2(I)=EXPON(I,N)*(ZETA2**2)
        D2(I)=COEF(I,N)*((2.0D0*A2(I)/PI)**0.75D0)
    end do
! Calculate electron integrals
    S12=0.0D0
    T11=0.0D0
    T12=0.0D0
    T22=0.0D0
    V11A=0.0D0
    V12A=0.0D0
    V22A=0.0D0
    V11B=0.0D0
    V12B=0.0D0
    V22B=0.0D0
    V1111=0.0D0
    V2111=0.0D0
    V2121=0.0D0
    V2211=0.0D0
    V2221=0.0D0
    V2222=0.0D0
! one-electron integrals I for atom1, J for atom2
    do I=1,N
        do J=1,N
            ! R1=0 R2=2
            RAP=A2(J)*R/(A1(I)+A2(J))
            RAP2=RAP**2
            RBP2=(R-RAP)**2
            S12=S12+S(A1(I),A2(J),R2)*D1(I)*D2(J)
            T11=T11+T(A1(I),A1(J),0.0D0)*D1(I)*D1(J)
            T12=T12+T(A1(I),A2(J),R2)*D1(I)*D2(J)
            T22=T22+T(A2(I),A2(J),0.0D0)*D2(I)*D2(J)
            V11A=V11A+V(A1(I),A1(J),0.0D0,0.0D0,ZA)*D1(I)*D1(J)
            V12A=V12A+V(A1(I),A2(J),R2,RAP2,ZA)*D1(I)*D2(J)
            V22A=V22A+V(A2(I),A2(J),0.0D0,R2,ZA)*D2(I)*D2(J)
            V11B=V11B+V(A1(I),A1(J),0.0D0,R2,ZB)*D1(I)*D1(J)
            V12B=V12B+V(A1(I),A2(J),R2,RBP2,ZB)*D1(I)*D2(J)
            V22B=V22B+V(A2(I),A2(J),0.0D0,0.0D0,ZB)*D2(I)*D2(J)
        end do
    end do
    do I=1,N
        do J=1,N
            do K=1,N
                do L=1,N
                    RAP=A2(I)*R/(A2(I)+A1(J))
                    RBP=R-RAP
                    RAQ=A2(K)*R/(A2(K)+A1(L))
                    RBQ=R-RAQ
                    RPQ=RAP-RAQ
                    RAP2=RAP*RAP
                    RBP2=RBP*RBP
                    RAQ2=RAQ*RAQ
                    RBQ2=RBQ*RBQ
                    RPQ2=RPQ*RPQ
                    V1111=V1111+TWOE(A1(I),A1(J),A1(K),A1(L),0.0D0,0.0D0,0.0D0)*D1(I)*D1(J)*D1(K)*D1(L)
                    V2111=V2111+TWOE(A2(I),A1(J),A1(K),A1(L),R2,0.0D0,RAP2)*D2(I)*D1(J)*D1(K)*D1(L)
                    V2121=V2121+TWOE(A2(I),A1(J),A2(K),A1(L),R2,R2,RPQ2)*D2(I)*D1(J)*D2(K)*D1(L)
                    V2211=V2211+TWOE(A2(I),A2(J),A1(K),A1(L),0.0D0,0.0D0,R2)*D2(I)*D2(J)*D1(K)*D1(L)
                    V2221=V2221+TWOE(A2(I),A2(J),A2(K),A1(L),0.0D0,R2,RBQ2)*D2(I)*D2(J)*D2(K)*D1(L)
                    V2222=V2222+TWOE(A2(I),A2(J),A2(K),A2(L),0.0D0,0.0D0,0.0D0)*D2(I)*D2(J)*D2(K)*D2(L)
                end do
            end do
        end do
    end do
    if (IOP .EQ. 0) return
    write(*,40)
40  FORMAT(3X,'R',10X,'ZETA1',6X,'ZETA2',6X,'S12',8X,'T11'/)
    write(*,50) R,ZETA1,ZETA2,S12,T11
50  FORMAT(5F11.6//)
    write(*,60)
60  FORMAT(3X,'T12',8X,'T22',8X,'V11A',7X,'V12A',7X,'V22A'/)
    write(*,50) T12,T22,V11A,V12A,V22A
    write(*,70)
70  FORMAT(3X,4HV11B,7X,4HV12B,7X,4HV22B,7X,'V1111',6X,'V2111'/)
    write(*,50) V11B,V12B,V22B,V1111,V2111
    write(*,80)
80  FORMAT(3X,5HV2121,6X,5HV2211,6X,5HV2221,6X,5HV2222/)
    write(*,50) V2121,V2211,V2221,V2222
    return
end subroutine

! Calculate the F FUNCTION F0 ONLY (s-type ORBITALS)
function F0(ARG)
    implicit none
    real*8 ::ARG
    real*8 :: PI=3.1415926535898D0
    real*8 :: F0
    if (ARG < 1.0D-6) then
        F0=1.0D0-ARG/3.0D0
    else
        F0=DSQRT(PI/ARG)*ERF(DSQRT(ARG))/2.0D0
    end if
    return
end function

! Calculate overlaps S
function S(A,B,RAB2)
    implicit none
    real*8 :: A,B,RAB2
    real*8 :: S
    real*8 :: PI=3.1415926535898D0
    S=(PI/(A+B))**1.5D0*DEXP(-A*B*RAB2/(A+B))
    return
end function

! Calculate kinetic energy integrals
function T(A,B,RAB2)
    implicit none
    real*8 :: A,B,RAB2
    real*8 :: T
    real*8 :: PI=3.1415926535898D0
    T=A*B/(A+B)*(3.0D0-2.0D0*A*B*RAB2/(A+B))*(PI/(A+B))**1.5D0*DEXP(-A*B*RAB2/(A+B))
    return
end function

! Calculate nuclear attraction integrals
function V(A,B,RAB2,RCP2,ZC)
    implicit none
    real*8 :: A,B,RAB2,RCP2,ZC
    real*8, external :: F0
    real*8 :: V
    real*8 :: PI=3.1415926535898D0
    V=2.0D0*PI/(A+B)*F0((A+B)*RCP2)*DEXP(-A*B*RAB2/(A+B))
    V=-V*ZC
    return
end function

! Calculate two-electron integrals
function TWOE(A,B,C,D,RAB2,RCD2,RPQ2)
    implicit none
    real*8 :: A,B,C,D,RAB2,RCD2,RPQ2
    real*8, external :: F0
    real*8 :: TWOE
    real*8 :: PI=3.1415926535898D0
    TWOE=2.0D0*(PI**2.5D0)/((A+B)*(C+D)*DSQRT(A+B+C+D))*F0((A+B)*(C+D)*RPQ2/(A+B+C+D))*DEXP(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D))
    return
end function

! THIS TAKES THE BASIC INTEGRALS FROM COMMON AND ASSEMBLES THE RELEVENT MATRICES
subroutine COLECT(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
    implicit none
    real*8 :: IOP,R,ZETA1,ZETA2,ZA,ZB
    integer :: N
    real*8 :: S(2,2),X(2,2),XT(2,2),H(2,2),F(2,2),G(2,2),C(2,2),FPRIME(2,2),CPRIME(2,2),P(2,2),OLDP(2,2),TT(2,2,2,2),E(2,2)
    common /MATRIX/ S,X,XT,H,F,G,C,FPRIME,CPRIME,P,OLDP,TT,E
    real*8 :: S12,T11,T12,T22,V11A,V12A,V22A,V11B,V12B,V22B,V1111,V2111,V2121,V2211,V2221,V2222
    common /INT/ S12,T11,T12,T22,V11A,V12A,V22A,V11B,V12B,V22B,V1111,V2111,V2121,V2211,V2221,V2222
    integer :: I,J,K,L
! core Hamiltonian
    H(1,1)=T11+V11A+V11B
    H(1,2)=T12+V12A+V12B
    H(2,1)=H(1,2)
    H(2,2)=T22+V22A+V22B
! S Matrix
    S(1,1)=1.0D0
    S(1,2)=S12
    S(2,1)=S(1,2)
    S(2,2)=1.0D0
! Schmidt orthogonalization for S Matrix
    X(1,1)=1.0D0/DSQRT(2.0D0*(1.0D0+S12))
    X(2,1)=X(1,1)
    X(1,2)=1.0D0/DSQRT(2.0D0*(1.0D0-S12))
    X(2,2)=-X(1,2)
! Transpose of X
    XT=transpose(X)
! Matrix of Two-electron integrals
    TT(1,1,1,1)=V1111
    TT(2,1,1,1)=V2111
    TT(1,2,1,1)=V2111
    TT(1,1,2,1)=V2111
    TT(1,1,1,2)=V2111
    TT(2,1,2,1)=V2121
    TT(1,2,2,1)=V2121
    TT(2,1,1,2)=V2121
    TT(1,2,1,2)=V2121
    TT(2,2,1,1)=V2211
    TT(1,1,2,2)=V2211
    TT(2,2,2,1)=V2221
    TT(2,2,1,2)=V2221
    TT(2,1,2,2)=V2221
    TT(1,2,2,2)=V2221
    TT(2,2,2,2)=V2222
    if (IOP /= 0) then
        CALL MATOUT(S,2,2,2,2,'S')
        CALL MATOUT(X,2,2,2,2,'X')
        CALL MATOUT(H,2,2,2,2,'H')
        write(*,10)
        do I=1,2
            do J=1,2
                do K=1,2
                    do L=1,2
                        write(*,20) I,J,K,L,TT(I,J,K,L)
                    end do
                end do
            end do
        end do
    end if
10  FORMAT(//)
20  FORMAT(3X,'(',4I2,')',F10.6)
end subroutine

subroutine SCF(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
    implicit none
    real*8 :: IOP,R,ZETA1,ZETA2,ZA,ZB
    integer :: N
    real*8 :: S(2,2),X(2,2),XT(2,2),H(2,2),F(2,2),G(2,2),C(2,2),FPRIME(2,2),CPRIME(2,2),P(2,2),OLDP(2,2),TT(2,2,2,2),E(2,2)
    common /MATRIX/ S,X,XT,H,F,G,C,FPRIME,CPRIME,P,OLDP,TT,E
    real*8 :: PI=3.1415926535898D0
    real*8 :: CRIT=1.0D-6
    integer :: max_iter=120,ITER,i,j,k,l,lwork,info
    real*8 :: energy,evals(2),dummy(1),ENT,delta
    real*8, allocatable :: work(:)
    ! Initial guess
    P = 0.0D0
    ! SCF cyclo
    do iter=1, max_iter
        write(*,'(///,4X,A,I2)') 'START OF ITERATION NUMBER =',iter
        ! Construct F Matrix
        ! F_uv = H_uv + Sum_ls P_ls * [ (uv|sl) - 0.5 * (ul|sv) ]
        F = H
        do i=1,2
            do j=1,2
                do k=1,2
                    do l=1,2
                        F(i,j)=F(i,j)+P(k,l)*(TT(i,j,l,k)-0.5D0*TT(i,k,l,j))
                    end do
                end do
            end do
        end do
        
        ! Calculate energy
        energy = 0.0D0
        do i = 1,2
            do j = 1,2
                energy = energy+0.5D0*P(i,j)*(H(i,j)+F(i,j))
            end do
        end do
        write(*,'(/,4X,A,D20.12)') 'ELECTRONIC ENERGY = ', energy
        
        ! Transition Fock matrix
        FPRIME=matmul(XT,matmul(F,X))
        
        ! DSYEV
        lwork=-1
        call dsyev('V','U',2,FPRIME,2,evals,dummy,lwork,info)
        lwork=int(dummy(1))
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',2,FPRIME,2,evals,work,lwork,info)
        E(1,1)=evals(1)
        E(2,2)=evals(2)
        ! r transition C
        C=matmul(X, FPRIME)
        
        ! Construct new P
        OLDP=P
        do i=1,2
            do j=1,2
                P(i,j)=0.0D0
                do k=1,1
                    P(i,j)=P(i,j)+2.0D0*C(i,k)*C(j,k)
                end do
            end do
        end do
        if (IOP > 2) then
            CALL MATOUT(F,2,2,2,2,'F')
            CALL MATOUT(E,2,2,2,2,'E')
            CALL MATOUT(C,2,2,2,2,'C')
            CALL MATOUT(P,2,2,2,2,'P')
        end if
        ! Calculate delta
        delta=0.0D0
        delta = sum((p(:,:) - oldp(:,:))**2)
        delta=DSQRT(delta/4.0D0)
        if (IOP > 0) then
            write(*,'(/,4X,A,F10.6,/)') 'DELTA(CONVERGENCE OF DENSITY MATRIX) =',delta
        end if
        if (delta<CRIT) exit
        if (iter==max_iter .and. delta>CRIT) then
            write(*,'(4X,A)') 'NO CONVERGENCE IN SCF'
        end if
    end do
    ENT=energy+ZA*ZB/R
    if (IOP /= 0) then
        write(*, '(//,4X,A,//,4X,A,D20.12,//,4X,A,D20.12)') &
        'CALCULATION CONVERGED', &
        'ELECTRONIC ENERGY = ', energy, &
        'TOTAL ENERGY = ', ENT
    end if
    if (IOP == 1) then
        call MATOUT(F, 2, 2, 2, 2, 'F')
        call MATOUT(C, 2, 2, 2, 2, 'C')
        call MATOUT(P, 2, 2, 2, 2, 'P')
    end if
        OLDP=matmul(P,S)
    if (IOP /= 0) then
        call MATOUT(OLDP, 2, 2, 2, 2, 'PS')
    end if
    return
end subroutine

! Print Matrix of size M by N
subroutine MATOUT(A,IM,INP,M,N,LABEL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IM, INP, M, N
    DOUBLE PRECISION, INTENT(IN) :: A(IM, INP)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER :: I, J, LOW, IHIGH
    CHARACTER(LEN=100) :: FMT_HEADER, FMT_ROW

    IHIGH = 0

    do
        LOW = IHIGH + 1
        IHIGH = IHIGH + 5
        IHIGH = MIN(IHIGH, N)
        write(*, '(/3X, "THE ", A, " ARRAY")') TRIM(LABEL)
        write(*, '(15X)', ADVANCE='NO')
        do I = LOW, IHIGH
            write(*, '(10X, I3, 6X)', ADVANCE='NO') I
        end do
        write(*, *)
  
        do I = 1, M
            write(*, '(I10, 5X)', ADVANCE='NO') I
            do J = LOW, IHIGH
                write(*, '(1X, D18.10)', ADVANCE='NO') A(I, J)
            end do
            write(*, *)
        end do
        IF (IHIGH >= N) EXIT
    end do
end subroutine