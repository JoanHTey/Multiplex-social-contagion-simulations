PROGRAM MSTREQ

    INTEGER :: N, NTOT, DT, STIME, I, K, M, K1, NEI, NNEI, KSUM, IOSTAT, F, SUM50, S, J
    INTEGER, ALLOCATABLE ::  ADJA(:,:)
    CHARACTER(LEN=20) :: INPR_STR
    REAL(8) :: NU, PR, PREP, EXPR, INPR, CS, ETA, GAMMA1, GAMMA2, Q
    REAL(8), ALLOCATABLE :: R(:), X(:), XSAVE(:), QLIST(:)

    ! Example content of 'INITIAL.txt':
    ! 10 3 0.01 0.5 1 1.0 0.1 0.2
    OPEN (1 , FILE = 'INITIAL.txt')
    READ (1,*) N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2
    CLOSE (1)
    NTOT= DIM * N
    ALLOCATE(X(NTOT), XSAVE(NTOT), ADJA(NTOT, N), R(NTOT), QLIST(NTOT))
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
    
    DO I = 1, N
        R(I) = 1.D0-(1.D0-1.D0/DBLE(SUM(ADJA(I,:))))**GAMMA1
        R(I+N) = 1.D0-(1.D0-1.D0/DBLE(SUM(ADJA(I,:))))**GAMMA2
    END DO

    NU = 0.5D0
    INPR = INPR * NU
    EXPR = ETA * INPR
    !IF (ETA*INPR.LE.1) THEN
    !    NU = 0.5D0
    !    INPR = INPR * NU
    !    EXPR = ETA * INPR
    !ELSE
    !    EXPR = 1.D0
    !    NU = 1.D0/(ETA * INPR)
    !    INPR = 1.D0/ETA
    !END IF

    DO I = 1, NTOT
        X(I) = 0
    END DO
    X(1)=1


    !DO DT = 1,STIME
    DO WHILE (MAXVAL(ABS(X - XSAVE)) > 1E-10)
        DO K=1,NTOT
            XSAVE(K)=X(K)
        END DO
        DO K=1,NTOT 
            Q=1
            DO J=1,ADJA(K,1)
                Q=Q*(1-INPR*R(K)*XSAVE(ADJA(K,J+1)))
            END DO
            Q=Q*(1-EXPR*XSAVE(MOD(K-1+N,NTOT)+1))
            X(K) = (1-XSAVE(K))*(1-Q) + (1-NU)*XSAVE(K) + NU*(1-Q)*XSAVE(K)
        END DO
    END DO

   
    ! Open the file
    OPEN(UNIT=1, FILE='output2INPR.txt', STATUS='unknown', ACTION='write', IOSTAT=IOSTAT)
    
    WRITE(1,*)(X(j), j=1,NTOT)
        DO K=1,NTOT 
            Q=1
            DO J=1,ADJA(K,1)
                Q=Q*(1-INPR*R(K)*X(ADJA(K,J+1)))
            END DO
            Q=Q*(1-EXPR*X(MOD(K-1+N,NTOT)+1))
            QLIST(k)=Q
        END DO
    WRITE(1,*)(QLIST(j), j=1,NTOT)
    CLOSE(1)
 
    DEALLOCATE(X, XSAVE, ADJA, R)
END PROGRAM MSTREQ