PROGRAM SIMULATION

    INTEGER :: N, NTOT, DT, STIME, I, K, M, K1, NEI, NNEI, KSUM, IOSTAT, F, SUM50, S, J,DISTIME,H,LAGS
    INTEGER, ALLOCATABLE :: X(:), SAVEX(:,:), ADJA(:,:), XP(:)

    REAL(8) :: NU, PR, PREP, EXPR, INPR, CS, ETA, GAMMA1, GAMMA2, BETAS(5)
    REAL(8), ALLOCATABLE :: R(:)

    ! Example content of 'INITIAL.txt':
    ! 10 3 0.01 0.5 1 1.0 0.1 0.2
    OPEN (1 , FILE = 'INITIAL.txt')
    READ (1,*) N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2, DISTIME, LAGS, NU
    CLOSE (1)
    
    OPEN (1, FILE = 'BETAS.txt', STATUS='old', ACTION='read', IOSTAT=IOSTAT)
    READ(1, *, IOSTAT=IOSTAT) (BETAS(I), I = 1, 5)
    CLOSE(1)
    WRITE(*,*) 'Betas:', (BETAS(I), I = 1, 5)

    NTOT = DIM * N
    ALLOCATE(X(NTOT), XP(NTOT), SAVEX(NTOT, M), ADJA(NTOT, N), R(NTOT))

    CALL setr1279(INT((N*DIM+M)*ETA*INPR + SYSTEM_CLOCK))

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


    PREP = 0.01D0
    DT = 1

    DO I = 1, N
        R(I) = 1.D0-(1.D0-1.D0/DBLE(ADJA(I,1)))**GAMMA1
        R(I+N) = 1.D0-(1.D0-1.D0/DBLE(ADJA(I+N,1)))**GAMMA2
    END DO
    
    DO I = 1, NTOT
        X(I)=0
    END DO


    DO I = 1, NTOT
        IF (r1279() < 0.05D0) THEN
            X(I) = 1
        END IF
    END DO

    DO I = 1, NTOT
        XP(I)=X(I)
    END DO

    OPEN(UNIT=1, FILE='output.bin', FORM='UNFORMATTED', ACCESS='SEQUENTIAL')

    DO H = 1, 5

        INPR = BETAS(H)

        !NU = 0.5D0
        INPR = INPR * NU
        EXPR = ETA * INPR
        DT = 1
        WRITE(*,*) 'YEY'
        DO WHILE (DT .LE. STIME)

            DO I = 1, NTOT
                X(I)=XP(I)
            END DO


            F = MOD(DT,50)

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

                K = I

                IF (X(K).EQ.1) THEN
                    IF (r1279() .LT. NU) THEN
                        XP(K) = 0
                    END IF  
                END IF
                    
                K1 = MOD(K-1+N,NTOT)+1
                PR = 1.d0
                    
                IF (X(K1).EQ.1) THEN
                    PR = PR * (1-EXPR)
                END IF

                NNEI = ADJA(K,1)
                KSUM = 0

                DO NEI = 2,NNEI+1
                    KSUM = KSUM + X(ADJA(K,NEI))
                END DO
                    
                PR = PR * (1-INPR*R(K))**KSUM 
                    
                IF (r1279() .GT. PR) THEN
                    XP(K) = 1              
                END IF

                

            END DO

            IF (MOD(DT-1, DISTIME) >= 0 .AND. MOD(DT-1, DISTIME) < LAGS) THEN
                WRITE(1) DT,(X(I), I=1, NTOT)
            END IF

            DT = DT + 1
        END DO
    END DO
    CLOSE(1)
END PROGRAM SIMULATION