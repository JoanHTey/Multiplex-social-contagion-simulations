PROGRAM SIMULATION

    INTEGER :: N, NTOT, DT, STIME, I, K, M, K1, NEI, NNEI, KSUM, IOSTAT, F, SUM50, S, J
    INTEGER, ALLOCATABLE :: X(:), SAVEX(:,:), ADJA(:,:)

    REAL(8) :: NU, PR, PREP, EXPR, INPR, CS, ETA, GAMMA1, GAMMA2
    REAL(8), ALLOCATABLE :: R(:)

    ! Example content of 'INITIAL.txt':
    ! 10 3 0.01 0.5 1 1.0 0.1 0.2
    OPEN (1 , FILE = 'INITIAL.txt')
    READ (1,*) N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2
    CLOSE (1)

    NTOT = DIM * N
    ALLOCATE(X(NTOT), SAVEX(NTOT, M), ADJA(NTOT, N), R(NTOT))

    CALL setr1279(N*DIM+M)

    ! Open the file and read the matrix
    OPEN(UNIT=1, FILE='adjacency_matrix.txt', STATUS='old', ACTION='read', IOSTAT=IOSTAT)
    IF (IOSTAT /= 0) THEN
        PRINT *, 'Error opening file'
        STOP
    END IF

    DO I = 1, NTOT
        READ(1, *) (ADJA(I, J), J = 1, N)
    END DO

    CLOSE(1)

    ! Print the matrix to verify
    !DO I = 1, NTOT
    !    WRITE(*, *) (ADJA(I, J), J = 1, N)
    !END DO


    PR = 1.d0
    PREP = 0.01D0
    NSUM = 1000
    DT = 1

    DO I = 1, N
        R(I) = 1.D0-(1.D0-1.D0/DBLE(SUM(ADJA(I,:))))**GAMMA1
        R(I+N) = 1.D0-(1.D0-1.D0/DBLE(SUM(ADJA(I,:))))**GAMMA2
    END DO

    IF (ETA*INPR.LE.1) THEN
        NU = 0.5
        INPR = INPR * NU
        EXPR = ETA * INPR
    ELSE
        EXPR = 0.5
        NU = 1.D0/(ETA * INPR)
        INPR = 1.D0/ETA

    END IF

    DO I = 1, NTOT
        X(I)=0
    END DO

    X(1) = 1
    !X(N+1) = 1


    OPEN(UNIT=1, FILE='output.txt', STATUS='unknown', ACTION='write', IOSTAT=IOSTAT)

    DO WHILE (DT .LE. STIME)
        print *, DT
        F = MOD(DT,50)

        !Test to see if we are in the equilibrium
        !IF (DT.GE.M.AND.F.EQ.0.AND.STIME-NSUM.GE.DT) THEN
        !    CS = SUM(X)-SUM50
        !END IF

        !IF (0.025D0 .GT. CS) THEN
        !        DT=STIME-NSUM-1
        !END IF
    
        !IF (F.EQ.0) THEN
        !    SUM50 = SUM(X)
        !END IF

        IF (DT .LT. M) THEN
            DO I = 1, NTOT
                SAVEX(I,DT) = X(I)
            END DO
        END IF

        IF (DT .GE. M.AND. r1279().LE.PREP) THEN
            K = MOD(INT(M*r1279()),M)+1
            DO I = 1, NTOT
                SAVEX(I,K) = X(I)
            END DO
        END IF

        IF (SUM(X).EQ.0) THEN
            DO WHILE (S.EQ.0)
                K = MOD(INT(M*r1279()),M)+1
                IF (SUM(SAVEX(:,K)).NE.0) THEN
                    DO I = 1, NTOT
                        X(I) = SAVEX(I,K)
                    END DO
                    S = 1
                END IF
            END DO
            S = 0
        END IF

        DO I = 1, NTOT

            K = MOD(INT(NTOT*r1279()),NTOT)+1

            IF (X(K).EQ.1) THEN
                IF (r1279() .GT. NU) THEN
                    X(K) = 0
                END IF  
            END IF

            IF (X(K).EQ.0) THEN
                K1 = MOD(K-1+N,NTOT)+1
                PR = 1
                IF (X(K1).EQ.1) THEN
                    PR = PR * (1-EXPR)
                END IF
                NNEI = ADJA(K,1)
                KSUM = 0

                DO NEI = 2,NNEI+1
                    KSUM = KSUM + X(ADJA(K,NEI))
                END DO
                WRITE(1,*), KSUM
                PR = PR * (1-INPR*R(K))**KSUM
                WRITE(1,*) PR   
                IF (r1279() .GT. PR) THEN
                    X(K) = 1
                    WRITE(1,*), 'flip'
                END IF
            END IF
        END DO
        DT = DT + 1
        !WRITE(1,*) (X(I),I=1,NTOT)
        WRITE(1,*) SUM(X)
    END DO
    CLOSE(1)
END PROGRAM SIMULATION
