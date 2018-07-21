      SUBROUTINE  ROTATE  (X,Y,CS,SN,N)
      INTEGER N
      REAL    X(N),Y(N),CS,SN
C
C
      REAL    XX
      INTEGER J
C
C
      DO 10 J = 1, N
         XX = X(J)
         X(J) = XX*CS + Y(J)*SN
10       Y(J) = Y(J)*CS - XX*SN
      RETURN
      END
