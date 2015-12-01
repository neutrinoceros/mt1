C 
C  $$$$  RA15.FOR   !!! AJOUT d'UN parametre :  TD, pour Temps au Depart 
C 
      SUBROUTINE RA15M(X,V,TD,TF0,XL,LL,NV,NCLASS,NOR,nsor,FORCE,SORTIE) 
C  Integrator RADAU by E. Everhart, Physics Department, University of Denver 
C  This 15th-order version, called RA15, is written out for faster execution. 
C  y'=F(y,t) is  NCLASS=1,  y"=F(y,t) is NCLASS= -2,  y"=F(y',y,t) is NCLASS=2 
C  TF is t(final) - t(initial). It may be negative for backward integration. 
C  NV = the number of simultaneous differential equations. 
C  The dimensioning below assumes NV will not be larger than 50. 
C  LL controls sequence size. Thus SS=10**(-LL) controls the size of a term. 
C  A typical LL-value is in the range 6 to 12 for this order 11 program. 
C  However, if LL.LT.0 then XL is the constant sequence size used. 
C  X and V enter as the starting position-velocity vector, and are output as 
C  the final position-velocity vector. 
C  Integration is in double precision. A 64-bit double-word is assumed. 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *4 TVAL,PW 
      DIMENSION X(90000),V(90000),F1(90000),FJ(90000),C(21),D(21),
     & R(21),Y(90000),Z(90000), 
     A     B(7,90000),G(7,90000),E(7,90000),BD(7,90000),H(8),W(7),
     & U(7),NW(8) 
      LOGICAL NPQ,NSF,NPER,NCL,NES 
      EXTERNAL FORCE,SORTIE   
      DATA NW/0,0,1,3,6,10,15,21/ 
      DATA ZERO, HALF, ONE,SR/0.0D0, 0.5D0, 1.0D0,1.4D0/ 
C  These H values are the Gauss-Radau spacings, scaled to the range 0 to 1, 
C  for integrating to order 15. 
      DATA H/         0.D0, .05626256053692215D0, .18024069173689236D0, 
     A.35262471711316964D0, .54715362633055538D0, .73421017721541053D0, 
     B.88532094683909577D0, .97752061356128750D0/ 
C  The sum of the H-values should be 3.73333333333333333 
      NPER=.FALSE. 
      NSF=.FALSE. 
      NCL=NCLASS.EQ.1 
      NPQ=NCLASS.LT.2 
C  y'=F(y,t)  NCL=.TRUE.    y"=F(y,t)  NCL=.FALSE.   y"=F(y',y,t) NCL=.FALSE. 
C  NCLASS=1   NPQ=.TRUE.    NCLASS= -2 NPQ=.TRUE.    NCLASS= 2    NPQ=.FALSE. 
C  NSF is .FALSE. on starting sequence, otherwise .TRUE. 
C  NPER is .TRUE. only on last sequence of the integration. 
C  NES is .TRUE. only if LL is negative. Then the sequence size is XL. 
      DIR=ONE 
      TF=TF0-TD                !!! NOUVELLE INSTRUCTION  !!!!!!!!!!!!!!!!!!!!!! 
      IF(TF.LT.ZERO) DIR=-ONE   
      NES=LL.LT.0 
      XL=DIR*DABS(XL) 
      PW=1./9. 
C  Evaluate the constants in the W-, U-, C-, D-, and R-vectors 
      DO 14 N=2,8 
      WW=N+N*N 
      IF(NCL) WW=N 
      W(N-1)=ONE/WW 
      WW=N 
  14  U(N-1)=ONE/WW 
      DO 22 K=1,NV 
      IF(NCL) V(K)=ZERO 
      DO 22 L=1,7 
      BD(L,K)=ZERO 
  22  B(L,K)=ZERO 
      W1=HALF 
      IF(NCL) W1=ONE 
      C(1)=-H(2) 
      D(1)=H(2) 
      R(1)=ONE/(H(3)-H(2)) 
      LA=1 
      LC=1 
      DO 73 K=3,7 
      LB=LA 
      LA=LC+1 
      LC=NW(K+1) 
      C(LA)=-H(K)*C(LB) 
      C(LC)=C(LA-1)-H(K) 
      D(LA)=H(2)*D(LB) 
      D(LC)=-C(LC) 
      R(LA)=ONE/(H(K+1)-H(2)) 
      R(LC)=ONE/(H(K+1)-H(K)) 
      IF(K.EQ.3) GO TO 73 
      DO 72 L=4,K 
      LD=LA+L-3 
      LE=LB+L-4 
      C(LD)=C(LE)-H(K)*C(LE+1) 
      D(LD)=D(LE)+H(L-1)*D(LE+1) 
  72  R(LD)=ONE/(H(K+1)-H(L-1)) 
  73  CONTINUE 
      SS=10.**(-LL) 
C  The statements above are used only once in an integration to set up the 
C  constants. They use less than a second of execution time.  Next set in 
C  a reasonable estimate to TP based on experience. Same sign as DIR. 
C  An initial first sequence size can be set with XL even with LL positive. 
      TP=0.1D0*DIR 
      IF(XL.NE.ZERO) TP=XL 
      IF(NES) TP=XL 
      IF(TP/TF.GT.HALF) TP=HALF*TF 
      NCOUNT=0 
 
C  An * is the symbol for writing on the monitor. The printer is unit 4. 
C  Line 4000 is the starting place of the first sequence. 
4000  NS=0 
      NF=0 
      NI=6 
      TM=TD     !ZERO               !!! CORRECTION  ***************************** 
      CALL FORCE (X, V, TM, F1) 
      NF=NF+1 
C Line 722 is begins every sequence after the first. First find new beta- 
C  values from the predicted B-values, following Eq. (2.7) in text. 
 722  DO 58 K=1,NV 
      G(1,K)=B(1,K)+D(1)*B(2,K)+ 
     X  D(2)*B(3,K)+D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K) 
      G(2,K)=            B(2,K)+ 
     X  D(3)*B(3,K)+D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K) 
      G(3,K)=B(3,K)+D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K) 
      G(4,K)=            B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K) 
      G(5,K)=                         B(5,K)+D(15)*B(6,K)+D(20)*B(7,K) 
      G(6,K)=                                      B(6,K)+D(21)*B(7,K) 
  58  G(7,K)=                                                   B(7,K) 

      T=TP 
      T2=T*T 
      IF(NCL) T2=T 
      TVAL=DABS(T) 
      IF(NS/nsor*nsor.EQ.NS) then  
        call SORTIE(X,V,TM,VIP) 
 
        WRITE(116,*) TM,VIP 
      end if 
 
C  Loop 175 is 6 iterations on first sequence and two iterations therafter. 
      DO 175 M=1,NI 
C  Loop 174 is for each substep within a sequence. 
      DO 174 J=2,8 
      JD=J-1 
      JDM=J-2 
      S=H(J) 
      Q=S 
      IF(NCL) Q=ONE 
C  Use Eqs. (2.9) and (2.10) of text to predict positions at each aubstep. 
C  These collapsed series are broken into two parts because an otherwise 
C  excellent  compiler could not handle the complicated expression. 
      DO 130 K=1,NV 
      A=W(3)*B(3,K)+S*(W(4)*B(4,K)+S*(W(5)*B(5,K)+S*(W(6)*B(6,K)+ 
     V   S*W(7)*B(7,K)))) 
      Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*(W(1)*B(1,K)+S*(W(2)*B(2,K) 
     X  +S*A)))) 
      IF(NPQ) GO TO 130 
C  Next are calculated the velocity predictors need for general class II. 
      A=U(3)*B(3,K)+S*(U(4)*B(4,K)+S*(U(5)*B(5,K)+S*(U(6)*B(6,K)+ 
     T    S*U(7)*B(7,K)))) 
      Z(K)=V(K)+S*T*(F1(K)+S*(U(1)*B(1,K)+S*(U(2)*B(2,K)+S*A))) 
 130  CONTINUE 

C  Find forces at each substep. 
      CALL FORCE(Y,Z,TM+S*T,FJ) 

      NF=NF+1 
      DO 171 K=1,NV 
C  Find G-value for the force FJ found at the current substep. This 
C  section, including the many-branched GOTO, uses Eq. (2.4) of text. 
      TEMP=G(JD,K) 
      GK=(FJ(K)-F1(K))/S 
      GO TO (102,102,103,104,105,106,107,108),J 
 102  G(1,K)=GK 
      GO TO 160 
 103  G(2,K)=(GK-G(1,K))*R(1) 
      GO TO 160 
 104  G(3,K)=((GK-G(1,K))*R(2)-G(2,K))*R(3) 
      GO TO 160 
 105  G(4,K)=(((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*R(6) 
      GO TO 160 
 106  G(5,K)=((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*R(9)- 
     X     G(4,K))*R(10) 
      GO TO 160 
 107  G(6,K)=(((((GK-G(1,K))*R(11)-G(2,K))*R(12)-G(3,K))*R(13)- 
     X     G(4,K))*R(14)-G(5,K))*R(15) 
      GO TO 160 
 108  G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-G(3,K))*R(18)- 
     X     G(4,K))*R(19)-G(5,K))*R(20)-G(6,K))*R(21) 
C  Upgrade all B-values. 
 160  TEMP=G(JD,K)-TEMP 
      B(JD,K)=B(JD,K)+TEMP 
C  TEMP is now the improvement on G(JD,K) over its former value. 
C  Now we upgrade the B-value using this dfference in the one term. 
C  This section is based on Eq. (2.5). 
      GO TO (171,171,203,204,205,206,207,208),J 
 203  B(1,K)=B(1,K)+C(1)*TEMP 
      GO TO 171 
 204  B(1,K)=B(1,K)+C(2)*TEMP 
      B(2,K)=B(2,K)+C(3)*TEMP 
      GO TO 171 
 205  B(1,K)=B(1,K)+C(4)*TEMP 
      B(2,K)=B(2,K)+C(5)*TEMP 
      B(3,K)=B(3,K)+C(6)*TEMP 
      GO TO 171 
 206  B(1,K)=B(1,K)+C(7)*TEMP 
      B(2,K)=B(2,K)+C(8)*TEMP 
      B(3,K)=B(3,K)+C(9)*TEMP 
      B(4,K)=B(4,K)+C(10)*TEMP 
      GO TO 171 
 207  B(1,K)=B(1,K)+C(11)*TEMP 
      B(2,K)=B(2,K)+C(12)*TEMP 
      B(3,K)=B(3,K)+C(13)*TEMP 
      B(4,K)=B(4,K)+C(14)*TEMP 
      B(5,K)=B(5,K)+C(15)*TEMP 
      GO TO 171 
 208  B(1,K)=B(1,K)+C(16)*TEMP 
      B(2,K)=B(2,K)+C(17)*TEMP 
      B(3,K)=B(3,K)+C(18)*TEMP 
      B(4,K)=B(4,K)+C(19)*TEMP 
      B(5,K)=B(5,K)+C(20)*TEMP 
      B(6,K)=B(6,K)+C(21)*TEMP 
 171  CONTINUE 
 174  CONTINUE 

      IF(NES.OR.M.LT.NI) GO TO 175 
C  Integration of sequence is over. Next is sequence size control. 
      HV=ZERO 
      DO 635 K=1,NV 
 635  HV=DMAX1(HV,DABS(B(7,K))) 
      HV=HV*W(7)/TVAL**7 
 175  CONTINUE 
      IF (NSF) GO TO 180 
      IF(.NOT.NES) TP=(SS/HV)**PW*DIR 
      IF(NES) TP=XL 
      IF(NES) GO TO 170 
      IF(TP/T.GT.ONE) GO TO 170 
   8  FORMAT (2X,2I2,2D18.10) 
      TP=.8D0*TP 
      NCOUNT=NCOUNT+1 
      IF(NCOUNT.GT.10) RETURN 
c      IF(NCOUNT.GT.1) WRITE(14,8) NOR,NCOUNT,T,TP 
C  Restart with 0.8x sequence size if new size called for is smaller than 
C  originally chosen starting sequence size on first sequence. 
      GO TO 4000 
 170  NSF=.TRUE. 
C Loop 35 finds new X and V values at end of sequence using Eqs. (2.11),(2.12) 
 180  DO 35 K=1,NV 
      X(K)=X(K)+V(K)*T+T2*(F1(K)*W1+B(1,K)*W(1)+B(2,K)*W(2)+B(3,K)*W(3) 
     X    +B(4,K)*W(4)+B(5,K)*W(5)+B(6,K)*W(6)+B(7,K)*W(7)) 
      IF(NCL) GO TO 35 
      V(K)=V(K)+T*(F1(K)+B(1,K)*U(1)+B(2,K)*U(2)+B(3,K)*U(3) 
     V    +B(4,K)*U(4)+B(5,K)*U(5)+B(6,K)*U(6)+B(7,K)*U(7)) 
  35  CONTINUE 
      NS=NS+1 
      if(NES) then 
        TM=TD+NS*T     !!! NOUVELLE INSTRUCTION  !!!!!!!!!!!!!!!!!!!!!! 
**        TM=NS*T 
      else 
        TM=TM+T 
      end if 
C  Return if done. 
      IF(.NOT.NPER) GO TO 78 
        call SORTIE(X,V,TD+NS*T,VIP) 
 
      WRITE(116,*) TM,VIP 
cc        call SORTIE(X,V,TF0,VIP) 
 
**      WRITE(14,7) NF,NS 
      RETURN 
C  Control on size of next sequence and adjust last sequence to exactly 
C  cover the integration span. NPER=.TRUE. set on last sequence. 
78    CALL FORCE (X,V,TM,F1) 
      NF=NF+1 
      IF(NES) GO TO 341 
      TP=DIR*(SS/HV)**PW 
      IF(TP/T.GT.SR) TP=T*SR 
 341  IF(NES) TP=XL 
**      IF(DIR*(TM+TP).LT.DIR*TF-1.D-8) GO TO 77 
      IF(DIR*(TM+TP).LT.DIR*(TF0)-1.D-8) GO TO 77   !!!essai 
cc      IF(abs((TM+TP)-TF0) .LT. 1.D-8) GO TO 77      !!!!   CORRECTION  ************ 
      TP=TF0-TM 
      NPER=.TRUE. 
C  Now predict B-values for next step. The predicted values from the preceding 
C  sequence were saved in the E-matrix. Te correction BD between the actual 
C  B-values found and these predicted values is applied in advance to the 
C  next sequence. The gain in accuracy is significant. Using Eqs. (2.13): 
  77  Q=TP/T 
      DO 39 K=1,NV 
      IF(NS.EQ.1) GO TO 31
 
      DO 20 J=1,7 
  20  BD(J,K)=B(J,K)-E(J,K) 
  31  E(1,K)=      Q*(B(1,K)+ 2.D0*B(2,K)+ 3.D0*B(3,K)+ 
     X           4.D0*B(4,K)+ 5.D0*B(5,K)+ 6.D0*B(6,K)+ 7.D0*B(7,K)) 
      E(2,K)=                Q**2*(B(2,K)+ 3.D0*B(3,K)+ 
     Y           6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+21.D0*B(7,K)) 
      E(3,K)=                             Q**3*(B(3,K)+ 
     Z           4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+35.D0*B(7,K)) 
      E(4,K)=   Q**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+35.D0*B(7,K)) 
      E(5,K)=                Q**5*(B(5,K)+ 6.D0*B(6,K)+21.D0*B(7,K)) 
      E(6,K)=                             Q**6*(B(6,K)+ 7.D0*B(7,K)) 
      E(7,K)=                                           Q**7*B(7,K) 

      DO 39 L=1,7 
  39  B(L,K)=E(L,K)+BD(L,K) 

C  Two iterations for every sequence after the first. 
      NI=2 
      if(.true.) GO TO 722 
      END 
 
