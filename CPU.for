      program TEST
      DIMENSION VAX(5)
      character*8 name
      DATA VAX / .830,.526,8.080,7.220,2.395 /
      open(10,FILE='test.out')
      NIT=500
      WRITE(*,*) 'LINPACK ............'
      CALL LINS(NIT,FLINS,DLINS,ELINS)
      NIT=500/2
      WRITE(*,*) 'DWHET ..............'
      CALL  DWHET(NIT,WHETS,DWHETS)
      NIT=500/8
      WRITE(*,*) 'SFLOPS .............'
      CALL  SFLOPS(NIT,FLOPS,DFLOPS,A)
      NIT=500
      WRITE(*,*) 'MATDOT .............'
      CALL  MATDOT(NIT,TDOT,DDOT,EDOT)
      NIT=500/5
      WRITE(*,*) 'MATINV .............'
      CALL MATINV(NIT,TINV,DINV,EINV)
      NIT=500/3
      WRITE(*,*) 'FFT ................'
      CALL FFT(NIT,TFFT,DFFT,EFFT)
      NIT=500/14
      WRITE(*,*) 'SIEVE ..............'
      CALL  SIEVE(NIT,TSIEVE,DSIEVE,ESIEVE)

      AVE1=WHETS/VAX(1)
      AVE2=FLOPS/VAX(2)
      AVE3=VAX(3)/TDOT
      AVE4=VAX(4)/TINV
      AVE5=VAX(5)/TFFT
      WRITE(*,*) AVE1,AVE2,AVE3,AVE4,AVE5

      AVE=WHETS/VAX(1)+FLOPS/VAX(2)+VAX(3)/TDOT+VAX(4)/TINV+VAX(5)/TFFT
      AVE=AVE/5
      WRITE(10,99)
      WRITE(10,100) NAME,AVE,FLINS,WHETS,FLOPS,TDOT,TINV,TFFT,TSIEVE
      WRITE(10,101) DLINS,DWHETS,DFLOPS,DDOT,DINV,DFFT,DSIEVE
      WRITE(10,102) ELINS,EDOT,EINV,EFFT,ESIEVE
      WRITE(*,99)
      WRITE(*,100) NAME,AVE,FLINS,WHETS,FLOPS,TDOT,TINV,TFFT,TSIEVE
      WRITE(*,101) DLINS,DWHETS,DFLOPS,DDOT,DINV,DFFT,DSIEVE
      WRITE(*,102) ELINS,EDOT,EINV,EFFT,ESIEVE
99    FORMAT(17X,'LINPACK  DWHETS  KFLOPS  MATDOT  MATINV',
     $ '     FFT   SIEVE')
100   FORMAT(1X,A7,4F15.0,4F15.6)
101   FORMAT(4X,'Times ',6X,7F11.2)
102   FORMAT(4X,'Errors',6X,1P1E11.1,16X,1P4E11.1)
      stop
      end

C *******************************************************************
      SUBROUTINE DWHET(NIT,WHETS,DIFTIM)
      DOUBLE PRECISION X1,X2,X3,X4,X,Y,Z,T,T1,T2,E1, begtim,endtim
      INTEGER   J,K,L,I, N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,ISAVE
      COMMON    T,T1,T2,E1(4),J,K,L
C
      I      =    40*NIT
      call dsec(BEGTIM)
C
C       ... the Whetstone code here ...
C
      ISAVE=I
      T=0.499975D00
      T1=0.50025D00
      T2=2.0D00
      N1=0
      N2=12*I
      N3=14*I
      N4=345*I
      N5=0
      N6=210*I
      N7=32*I
      N8=899*I
      N9=616*I
      N10=0
      N11=93*I
      N12=0
      X1=1.0D0
      X2=-1.0D0
      X3=-1.0D0
      X4=-1.0D0
      IF(N1)19,19,11
 11   DO 18 I=1,N1,1
      X1=(X1+X2+X3-X4)*T
      X2=(X1+X2-X3+X4)*T
      X4=(-X1+X2+X3+X4)*T
      X3=(X1-X2+X3+X4)*T
 18   CONTINUE
 19   CONTINUE
      E1(1)=1.0D0
      E1(2)=-1.0D0
      E1(3)=-1.0D0
      E1(4)=-1.0D0
      IF(N2)29,29,21
 21   DO 28 I=1,N2,1
      E1(1)=(E1(1)+E1(2)+E1(3)-E1(4))*T
      E1(2)=(E1(1)+E1(2)-E1(3)+E1(4))*T
      E1(3)=(E1(1)-E1(2)+E1(3)+E1(4))*T
      E1(4)=(-E1(1)+E1(2)+E1(3)+E1(4))*T
 28   CONTINUE
 29   CONTINUE
      IF(N3)39,39,31
 31   DO 38 I=1,N3,1
 38   CALL PA(E1)
 39   CONTINUE
      J=1
      IF(N4)49,49,41
 41   DO 48 I=1,N4,1
      IF(J-1)43,42,43
 42   J=2
      GOTO44
 43   J=3
 44   IF(J-2)46,46,45
 45   J=0
      GOTO47
 46   J=1
 47   IF(J-1)411,412,412
 411  J=1
      GOTO48
 412  J=0
 48   CONTINUE
 49   CONTINUE
      J=1
      K=2
      L=3
      IF(N6)69,69,61
 61   DO 68 I=1,N6,1
      J=J*(K-J)*(L-K)
      K=L*K-(L-J)*K
      L=(L-K)*(K+J)
      E1(L-1)=J+K+L
      E1(K-1)=J*K*L
 68   CONTINUE
 69   CONTINUE
      X=0.5D0
      Y=0.5D0
      IF(N7)79,79,71
 71   DO 78 I=1,N7,1
      X=T*DATAN(T2*DSIN(X)*DCOS(X)/(DCOS(X+Y)+DCOS(X-Y)-1.0D0))
      Y=T*DATAN(T2*DSIN(Y)*DCOS(Y)/(DCOS(X+Y)+DCOS(X-Y)-1.0D0))
 78   CONTINUE
 79   CONTINUE
      X=1.0D0
      Y=1.0D0
      Z=1.0D0
      IF(N8)89,89,81
 81   DO 88 I=1,N8,1
 88   CALL P3(X,Y,Z)
 89   CONTINUE
      J=1
      K=2
      L=3
      E1(1)=1.0D0
      E1(2)=2.0D0
      E1(3)=3.0D0
      IF(N9)99,99,91
 91   DO 98 I=1,N9,1
 98   CALL P0
 99   CONTINUE
      J=2
      K=3
      IF(N10)109,109,101
 101  DO 108 I=1,N10,1
      J=J+K
      K=J+K
      J=J-K
      K=K-J-J
 108  CONTINUE
 109  CONTINUE
      X=0.75D0
      IF(N11)119,119,111
 111  DO 118 I=1,N11,1
 118  X=DSQRT(DEXP(DLOG(X)/T1))
 119  CONTINUE
      I = ISAVE

C       ... the whetstone ends here

      call dsec(ENDTIM)
      DIFTIM = ENDTIM - BEGTIM
      WHETS = FLOAT(I)/DIFTIM/10
C        .....  RISULTATO IN MEGA WHETSTONE
      RETURN
      END
C  ************************************************************
      SUBROUTINE PA(E)
      DOUBLE PRECISION T,T1,T2,E,E1
      COMMON    T,T1,T2,E1(4),J,K,L
      DIMENSION E(4)
      J=0
 1    E(1)=(E(1)+E(2)+E(3)-E(4))*T
      E(2)=(E(1)+E(2)-E(3)+E(4))*T
      E(3)=(E(1)-E(2)+E(3)+E(4))*T
      E(4)=(-E(1)+E(2)+E(3)+E(4))/T2
      J=J+1
      IF(J-6)1,2,2
 2    CONTINUE
      RETURN
      END
C  ************************************************************
      SUBROUTINE P0
      DOUBLE PRECISION T,T1,T2,E1
      COMMON T,T1,T2,E1(4),J,K,L
      E1(J)=E1(K)
      E1(K)=E1(L)
      E1(L)=E1(J)
      RETURN
      END
C  ************************************************************
      SUBROUTINE P3(X,Y,Z)
      DOUBLE PRECISION T,T1,T2,X1,Y1,X,Y,Z
      COMMON    T,T1,T2,E1(4),J,K,L
      X1=X
      Y1=Y
      X1=T*(X1+Y1)
      Y1=T*(X1+Y1)
      Z=(X1+Y1)/T2
      RETURN
      END
C *******************************************************************
      SUBROUTINE SFLOPS(NIT,FLOPS,DIFTIM,A)
      doubleprecision t1,t2
      REAL X,Y,Z
      DATA X,Y,Z/1.,1.,1./
      CALL dsec(T1)
      DO 33 KK=1,NIT
      DO 1 J=1,500
      DO 1 I=1,500
      Z=X-Y*Z
      X=Z*X-Y
      Y=Z-X*Y
1     CONTINUE
33    CONTINUE
      CALL dsec(T2)
      DIFTIM=T2-T1
      FLOPS=6.*FLOAT(NIT)/(4.*DIFTIM)
C ...  RISULTATO IN MEGAFLOPS
      A=X+Y+Z
      RETURN
      END
C *******************************************************************
      SUBROUTINE FFT(NIT,Time,DIFTIM,ERROR)
      PARAMETER(NN=1024,NN2=2*NN)
      DIMENSION DATA(NN2),DCMP(NN2)
      doubleprecision t1,t2
      DO 1 I=1,2*NN-1,2
      DATA(I)=1.0/((0.5*(I-NN-1.0)/NN)**2 +1.0)
      DATA(I+1)=(0.25*(I-NN-1.0)/NN)*EXP(-(0.5*(I-NN-1.0)/NN)**2)
      DCMP(I)=DATA(I)
      DCMP(I+1)=DATA(I+1)
1     CONTINUE
      CALL dsec(T1)
      DO 33 KK=1,10*NIT
      ISIGN =1
      CALL FOUR1(DATA,NN,ISIGN)
      ISIGN =-1
      CALL FOUR1(DATA,NN,ISIGN)
      DO 2 I=1,NN2
2     DATA(I)=DATA(I)/NN
33    CONTINUE
      CALL dsec(T2)
      DIFTIM=T2-T1
      TIME=(T2-T1)/(2.*NIT)
      I=NN/2
      ERROR=ABS(DATA(I)-DCMP(I))
      RETURN
      END
C
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      DOUBLEPRECISION  WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
	IF(J.GT.I)THEN
	  TEMPR=DATA(J)
	  TEMPI=DATA(J+1)
	  DATA(J)=DATA(I)
	  DATA(J+1)=DATA(I+1)
	  DATA(I)=TEMPR
	  DATA(I+1)=TEMPI
	ENDIF
	M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
	  J=J-M
	  M=M/2
	GO TO 1
	ENDIF
	J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
	ISTEP=2*MMAX
	THETA=6.28318530717959D0/(ISIGN*MMAX)
	WPR=-2.D0*DSIN(0.5D0*THETA)**2
	WPI=DSIN(THETA)
	WR=1.D0
	WI=0.D0
	DO 13 M=1,MMAX,2
	  DO 12 I=M,N,ISTEP
	    J=I+MMAX
	    TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
	    TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
	    DATA(J)=DATA(I)-TEMPR
	    DATA(J+1)=DATA(I+1)-TEMPI
	    DATA(I)=DATA(I)+TEMPR
	    DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
	  WTEMP=WR
	  WR=WR*WPR-WI*WPI+WR
	  WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
	MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
C *******************************************************************
      SUBROUTINE MATDOT(NIT,Time,DIFTIM,ERROR)
      PARAMETER(N=100)
      doubleprecision t1,t2
      COMMON /UNO/ A(N,N),B(N,N),C(N,N)
      DO 1 I=1,N
      DO 1 J=1,N
	A(I,J)=2.75*i/100.
1       B(I,J)=3.25*J/75.
      call dsec(t1)
	do 33 k=1,nit
      DO 2 I=1,N
      DO 2 J=1,N
      CC=0.
      DO 3 L=1,N
3       CC=CC+A(I,L)*B(L,J)
2      C(I,J)=CC
33      continue
      call dsec(t2)
	DIFTIM=T2-T1
      Time=(t2-t1)/nit
      ERROR= ABS( C(20,20) - 47.666648)
      RETURN
      END
C *******************************************************************
      SUBROUTINE MATINV(NIT,TIME,DIFTIM,ERROR)
c ******************************************************
c ***   inversione di una matrice **********************
c ******************************************************
      PARAMETER(N=100)
      doubleprecision t1,t2
      DIMENSION INDX(N),WI(N,N)
      COMMON /UNO/ A(N,N),AI(N,N),W(N,N)
      DO 1 I=1,N
      DO 2 J=1,N
      WI(I,J)=0.
2     W(I,J)=1.
      WI(I,I)=1.
1     W(I,I)=2.
      call dsec(t1)
      DO 33 KK=1,NIT
      DO 3 I=1,N
      DO 3 J=1,N
      AI(I,J)=WI(I,J)
3     A(I,J)=W(I,J)
      CALL LUDCMP(A,N,N,INDX,D)
      DO 4 J=1,N
4     CALL LUBKSB(A,N,N,INDX,AI(1,J))
33    CONTINUE
      call dsec(t2)
      DIFTIM=T2-T1
      Time=(t2-t1)/NIT
      Error= ABS(AI(N,N)-.990098)
      RETURN
      END

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
	LL=INDX(I)
	SUM=B(LL)
	B(LL)=B(I)
	IF (II.NE.0)THEN
	  DO 11 J=II,I-1
	    SUM=SUM-A(I,J)*B(J)
11        CONTINUE
	ELSE IF (SUM.NE.0.) THEN
	  II=I
	ENDIF
	B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
	SUM=B(I)
	IF(I.LT.N)THEN
	  DO 13 J=I+1,N
	    SUM=SUM-A(I,J)*B(J)
13        CONTINUE
	ENDIF
	B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
	AAMAX=0.
	DO 11 J=1,N
	  IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
	IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
	VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
	IF (J.GT.1) THEN
	  DO 14 I=1,J-1
	    SUM=A(I,J)
	    IF (I.GT.1)THEN
	      DO 13 K=1,I-1
		SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
	      A(I,J)=SUM
	    ENDIF
14        CONTINUE
	ENDIF
	AAMAX=0.
	DO 16 I=J,N
	  SUM=A(I,J)
	  IF (J.GT.1)THEN
	    DO 15 K=1,J-1
	      SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
	    A(I,J)=SUM
	  ENDIF
	  DUM=VV(I)*ABS(SUM)
	  IF (DUM.GE.AAMAX) THEN
	    IMAX=I
	    AAMAX=DUM
	  ENDIF
16      CONTINUE
	IF (J.NE.IMAX)THEN
	  DO 17 K=1,N
	    DUM=A(IMAX,K)
	    A(IMAX,K)=A(J,K)
	    A(J,K)=DUM
17        CONTINUE
	  D=-D
	  VV(IMAX)=VV(J)
	ENDIF
	INDX(J)=IMAX
	IF(J.NE.N)THEN
	  IF(A(J,J).EQ.0.)A(J,J)=TINY
	  DUM=1./A(J,J)
	  DO 18 I=J+1,N
	    A(I,J)=A(I,J)*DUM
18        CONTINUE
	ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
C *******************************************************************
      SUBROUTINE SIEVE(NIT,TIME,DIFTIM,ERROR)
      doubleprecision t1,t2
      LOGICAL FLAGS(8191)
      INTEGER I,J,K,COUNT,ITER,PRIME
      call dsec(t1)
      DO 92 ITER=1,100*NIT
      COUNT=0
      I=0
      DO 10 I=1,8191
10      FLAGS(I)= .TRUE.
      DO 91 I=1,8191
      IF(.NOT. FLAGS(I)) GO TO 91
      PRIME=I+I+1
      COUNT=COUNT +1
      K=I+PRIME
      IF(K.GT.8191) GO TO 91
      DO 60 J=K,8191,PRIME
60      FLAGS(J)=.FALSE.
91      CONTINUE
92      CONTINUE
      call dsec(t2)
	DIFTIM=T2-T1
      Time = (t2-t1)/(NIT*10)
      ERROR=(PRIME-16381)**2+(COUNT-1899)**2
      RETURN
      END
C *******************************************************************
      subroutine lins(NTIMES,FLOPS,TIME,ERROR)
c     this program was updated on 10/12/92 to correct a
c     problem with the random number generator.
c     modified by sandro on 7/3/92 non era stabile
c     nel calcolo dei tempi sul pc
      parameter(n=100,lda=200)
      real a(lda,lda),b(lda),x(lda)
      real cray,ops,norma,normx
      real resid,residn,eps,epslon
      doubleprecision t1,t2,t3
      integer ipvt(lda)
      write(*,*) ' Times are reported for matrices of order ',n
      ops = (2.0e0*float(n)**3)/3.0e0 + 2.0e0*float(n)**2
      cray = .056
      eps = epslon(1.0)
      call matgen(a,lda,n,b,norma)
      call sgefa(a,lda,n,ipvt,info)
      call sgesl(a,lda,n,ipvt,b,0)
c ***  compute a residual to verify results.
      do 10 i = 1,n
   10    x(i) = b(i)
      call matgen(a,lda,n,b,norma)
      do 20 i = 1,n
   20    b(i) = -b(i)
      call smxpy(n,b,n,lda,x,a)
      resid = 0.0
      normx = 0.0
      error = 0.0
      do 30 i = 1,n
	 resid = amax1( resid, abs(b(i)) )
	 normx = amax1( normx, abs(X(I)) )
   30    error = amax1( error, abs(1-x(i)) )
      RESIDn = RESID/( N*NORMA*NORMX*EPS )
      call dsec(t1)
      do 40 k=1,ntimes
 40      call matgen(a,lda,n,b,norma)
      call dsec(t2)
      do 50 k=1,ntimes
	 call matgen(a,lda,n,b,norma)
	 call sgefa(a,lda,n,ipvt,info)
 50      call sgesl(a,lda,n,ipvt,b,0)
      call dsec(t3)
      write(*,*)' norm.resid       resid      machep       error',
     $ '          T1          T2'
      write(*,'(1p6e12.4,/)') residn,resid,eps,error,t3-t2,t2-t1
      time = ((t3-t2)-(t2-t1))
      flops = ntimes*ops/(1.0e6*time)
      return
      end
      subroutine matgen(a,lda,n,b,norma)
      integer lda,n,init(4),i,j
      real a(lda,1),b(1),norma,ran
      init(1) = 1
      init(2) = 2
      init(3) = 3
      init(4) = 1325
      norma = 0.0
      do 30 j = 1,n
	 do 20 i = 1,n
	    a(i,j) = rana(init) - .5
	    norma = amax1(abs(a(i,j)), norma)
   20    continue
   30 continue
      do 35 i = 1,n
	  b(i) = 0.0
   35 continue
      do 50 j = 1,n
	 do 40 i = 1,n
	    b(i) = b(i) + a(i,j)
   40    continue
   50 continue
      return
      end
      subroutine sgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      real a(lda,1)
      real t
      integer isamax,j,k,kp1,l,nm1
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
	 kp1 = k + 1
	 l = isamax(n-k+1,a(k,k),1) + k - 1
	 ipvt(k) = l
	 if (a(l,k) .eq. 0.0e0) go to 40
	    if (l .eq. k) go to 10
	       t = a(l,k)
	       a(l,k) = a(k,k)
	       a(k,k) = t
   10       continue
	    t = -1.0e0/a(k,k)
	    call sscal(n-k,t,a(k+1,k),1)
	    do 30 j = kp1, n
	       t = a(l,j)
	       if (l .eq. k) go to 20
		  a(l,j) = a(k,j)
		  a(k,j) = t
   20          continue
	       call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
	 go to 50
   40    continue
	    info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end
      subroutine sgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      real a(lda,1),b(1)
      real sdot,t
      integer k,kb,l,nm1
      nm1 = n - 1
      if (job .ne. 0) go to 50
	 if (nm1 .lt. 1) go to 30
	 do 20 k = 1, nm1
	    l = ipvt(k)
	    t = b(l)
	    if (l .eq. k) go to 10
	       b(l) = b(k)
	       b(k) = t
   10       continue
	    call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
	 do 40 kb = 1, n
	    k = n + 1 - kb
	    b(k) = b(k)/a(k,k)
	    t = -b(k)
	    call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
	 do 60 k = 1, n
	    t = sdot(k-1,a(1,k),1,b(1),1)
	    b(k) = (b(k) - t)/a(k,k)
   60    continue
	 if (nm1 .lt. 1) go to 90
	 do 80 kb = 1, nm1
	    k = n - kb
	    b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
	    l = ipvt(k)
	    if (l .eq. k) go to 70
	       t = b(l)
	       b(l) = b(k)
	       b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine saxpy(n,da,dx,incx,dy,incy)
      real dx(1),dy(1),da
      integer i,incx,incy,ix,iy,n
      if(n.le.0)return
      if (da .eq. 0.0e0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dy(iy) = dy(iy) + da*dx(ix)
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
   20 continue
      do 30 i = 1,n
	dy(i) = dy(i) + da*dx(i)
   30 continue
      return
      end
      real function sdot(n,dx,incx,dy,incy)
      real dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,n
      sdot = 0.0e0
      dtemp = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dtemp = dtemp + dx(ix)*dy(iy)
	ix = ix + incx
	iy = iy + incy
   10 continue
      sdot = dtemp
      return
   20 continue
      do 30 i = 1,n
	dtemp = dtemp + dx(i)*dy(i)
   30 continue
      sdot = dtemp
      return
      end
      subroutine  sscal(n,da,dx,incx)
      real da,dx(1)
      integer i,incx,n,nincx
      if(n.le.0)return
      if(incx.eq.1)go to 20
      nincx = n*incx
      do 10 i = 1,nincx,incx
	dx(i) = da*dx(i)
   10 continue
      return
   20 continue
      do 30 i = 1,n
	dx(i) = da*dx(i)
   30 continue
      return
      end
      integer function isamax(n,dx,incx)
      real dx(1),dmax
      integer i,incx,ix,n
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
      ix = 1
      dmax = abs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
	 if(abs(dx(ix)).le.dmax) go to 5
	 isamax = i
	 dmax = abs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 dmax = abs(dx(1))
      do 30 i = 2,n
	 if(abs(dx(i)).le.dmax) go to 30
	 isamax = i
	 dmax = abs(dx(i))
   30 continue
      return
      end
      REAL FUNCTION EPSLON (X)
      REAL X
      REAL A,B,C,EPS
      A = 4.0E0/3.0E0
   10 B = A - 1.0E0
      C = B + B + B
      EPS = ABS(C-1.0E0)
      IF (EPS .EQ. 0.0E0) GO TO 10
      EPSLON = EPS*ABS(X)
      RETURN
      END
      SUBROUTINE SMXPY (N1, Y, N2, LDM, X, M)
      REAL Y(*), X(*), M(LDM,*)
      J = MOD(N2,2)
      IF (J .GE. 1) THEN
	 DO 10 I = 1, N1
	    Y(I) = (Y(I)) + X(J)*M(I,J)
   10    CONTINUE
      ENDIF
      J = MOD(N2,4)
      IF (J .GE. 2) THEN
	 DO 20 I = 1, N1
	    Y(I) = ( (Y(I))
     $             + X(J-1)*M(I,J-1)) + X(J)*M(I,J)
   20    CONTINUE
      ENDIF
      J = MOD(N2,8)
      IF (J .GE. 4) THEN
	 DO 30 I = 1, N1
	    Y(I) = ((( (Y(I))
     $             + X(J-3)*M(I,J-3)) + X(J-2)*M(I,J-2))
     $             + X(J-1)*M(I,J-1)) + X(J)  *M(I,J)
   30    CONTINUE
      ENDIF
      J = MOD(N2,16)
      IF (J .GE. 8) THEN
	 DO 40 I = 1, N1
	    Y(I) = ((((((( (Y(I))
     $             + X(J-7)*M(I,J-7)) + X(J-6)*M(I,J-6))
     $             + X(J-5)*M(I,J-5)) + X(J-4)*M(I,J-4))
     $             + X(J-3)*M(I,J-3)) + X(J-2)*M(I,J-2))
     $             + X(J-1)*M(I,J-1)) + X(J)  *M(I,J)
   40    CONTINUE
      ENDIF
      JMIN = J+16
      DO 60 J = JMIN, N2, 16
	 DO 50 I = 1, N1
	    Y(I) = ((((((((((((((( (Y(I))
     $             + X(J-15)*M(I,J-15)) + X(J-14)*M(I,J-14))
     $             + X(J-13)*M(I,J-13)) + X(J-12)*M(I,J-12))
     $             + X(J-11)*M(I,J-11)) + X(J-10)*M(I,J-10))
     $             + X(J- 9)*M(I,J- 9)) + X(J- 8)*M(I,J- 8))
     $             + X(J- 7)*M(I,J- 7)) + X(J- 6)*M(I,J- 6))
     $             + X(J- 5)*M(I,J- 5)) + X(J- 4)*M(I,J- 4))
     $             + X(J- 3)*M(I,J- 3)) + X(J- 2)*M(I,J- 2))
     $             + X(J- 1)*M(I,J- 1)) + X(J)   *M(I,J)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
      REAL FUNCTION rana( ISEED )
      INTEGER            ISEED( 4 )
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            IPW2
      REAL               R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
      INTEGER            IT1, IT2, IT3, IT4
      INTRINSIC          MOD, REAL
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RANA = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R*
     $         ( REAL( IT4 ) ) ) ) )
      RETURN
      END

********************************************************************
      SUBROUTINE DSEC(DT)
!      INTEGER*2 IH,IM,IS,IHU
      REAL*8 DT
!      CALL GETTIM (IH,IM,IS,IHU)
!      DT = (IH*3600.D0+IM*60.D0+IS)+IHU/100.D0
      call cpu_time(DT)

      RETURN
      END
********************************************************************