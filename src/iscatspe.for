      SUBROUTINE ISCATR
C* ROUTINE FOR THEORETICALLY COMPUTING INCOHERENT POWER SPECTRA
C  AND THEIR DERIVATIVES WITH RESPECT TO ION AND ELECTRON TEMPERATURES,
C  ELECTRON DENSITY, COLLISION FREQUENCY AND ION COMPOSITION.
C  BY WES SWARTZ
C =============================================================
C
C  CALL "ISCATR" ONCE FOR SYSTEM INITIALIZATION.
C  CALL "ISFITC" ONCE FOR INITIALIZATION OF EACH FITTING TYPE.
C  CALL "SPECTR" FOR FULL SPECTRA OR ACF IN ONE CALL AFTER "ISCATR"
C  AND "ISFITC".
C
C  DEFINITION OF BASIC TERMS:
C  IONOSPHERIC PARAMETERS
C    "TI"=ION TEMPERATURE (K).
C    "TE"=ELECTRON TEMPERATURE (K).
C    "EN"=ELECTRON DENSITY (CM-3)
C    "CF"=COLLISION FREQENCY (KHZ)
C    "FRACTN(J-1)"=FRACTIONAL COMPOSITION OF J'TH ION SPECIES
C    "SYSTEM"=SYSTEM SCALING FACTOR
C    "TR"="TE/TI"
      COMMON /PARM/TI,TE,EN,CF,FRACTN(2),SYSTEM /CTR/TR
C
C  GEOGRAPHIC PARAMETERS:
C    "GFE"=ELECTRON GYRO FREQUENCY (KHZ)
C    "ALPHA"=ANGLE BETWEEN MAGNETIC FIELD AND SCATTERING VECTOR (DEG).
C    "NJ"=NUMBER OF ION SPECIES.
C    "AM(J)"=MOL MASS OF J'TH ION SPECIES
C           = 1. FOR H+
C           = 4. FOR HE+
C           =16. FOR O+
C           =30. FOR NO+
C           =32. FOR O2+
C           =31. FOR COMBINING NO+ AND O2+
C           =55.8 FOR FE+
C           =24.3 FOR MG+
C
C    "DTAU"=BASIC SAMPLE INTERVAL FOR ACF (MS).
C  SYSTEM PARAMETERS:
C    "F(K)"=FREQUENCY VALUES FOR SPECTRUM (KHZ).
C    "COSINE(L,K)"=ARRAY FOR FOURIER TRANSFORM
C    "TF"=TRANSMITTER FREQUENCY (MHZ).
C        =49.2 FOR JICAMARCA.
C        =430. FOR ARECIBO.
C        =1290. FOR CHATANIKA.
C        =933. FOR EISCAT.
C    "NDF"=NUMBER OF FREQUENCY INTERVALS.
C    "TAU(L)"=ACF LAGS (US).
      COMMON /PM/DTAU,TF0,TF,GFE,ALPHA,AM(3,9),NDF,NJ,IREGN
C
C  CONTROL PARAMETERS:
C    "IFIT(I)"=1 FOR FITTING I'TH PARAMETER.
C    "IFIT(I)"=4 FOR HOLDING I'TH PARAMETER CONSTANT.
C
C   SPECIAL CASES:
C    "IFIT(2)"=1 FOR FITTING "TE".
C             =2 FOR "TE/TI" FIXED.
C             =4 FOR "TE" FIXED.
C    "IFIT(5)"=1 FOR FITTING COMPSITION.
C             =2 FOR GIVEN NONZERO COMPOSITION WITHOUT FITTING.
      COMMON /FITC/IFITR(8,9),IFFIT(7),NP
C
C COMPLEX OPERATIONS HAVE BEEN ELIMINATED.
C CODE FOR SCATTERING PERPINDICULAR TO MAGNETIC FIELD IS NOT IMPLIMENTED
C IN THIS VERSION.
C
      REAL SQMI(3),YIR1(100),YII1(100),GP(7),S(100)
      DIMENSION COSINE(36,100),FUNC(100,6),DYIRDP(100,2),DYIIDP(100,2),
     >          IP(7),IFIT(8),DSDTI(100)
C
      COMMON /FUNCT/ACF(36,6),REFUNC(36,6),AIFUNC(36,6),TAU(36),NLAGS
     >       /CSPECT/YII(100),YIR(100),DYIRDT(100),DYIIDT(100),
     >               DYIRDC(100),DYIIDC(100)
      COMMON /PLAS/ZR(100),ZI(100),THETA(100),PSIION,F(120),SCALEF,NDX
     >       /FITP/ITI,ITE,INE,ICF,IP1,IP2,ISY
      EQUIVALENCE (GP(1),TI), (IP(1),ITI), (FUNC(1,1),YII(1),S(1)),
     >  (DYIRDP(1,2),YIR1(1)), (DYIIDP(1,2),YII1(1)), (DSDTI(1),YIR(1))
C======================================================================
      NDX    = NDF
      TR     = 1.
C  SET UP COSINE ARRAY FOR FOURIER TRANSFORM:
      NDFM   = NDF-1
      DO 2 K = 2,NDFM
    2 COSINE(1,K) = (F(K+1)-F(K-1))*0.5
      COSINE(1,1) = 0.5*F(2)
      COSINE(1,NDF) = 0.5*(F(NDF)-F(NDFM))
      DO 3 L = 2,36
      TP2    = TAU(L)*6.2831854
      TP2SQ  = TP2**2
      COSINE(L,1) = (1.-COS(TP2*F(2)))/(F(2)*TP2SQ)
      DO 3 K = 2,NDF
      COSINE(L,K) = ((COS(TP2*F(K))-COS(TP2*F(K-1)))/
     >                (F(K)-F(K-1))
     >               +(COS(TP2*F(K))-COS(TP2*F(K+1)))
     >               /(F(K+1)-F(K)))/TP2SQ
   3  CONTINUE
      RETURN
C =====================================================================
      ENTRY ISFITC
C
C NUMERICAL CONSTANT FOR FREQUENCY NORMALIZATION IS:
C          27226. = C*SQRT(ME/2./BOLTZK)/2.
      FNORMC = 27.226/TF
      PHIEL  = FNORMC*GFE
C NUMERICAL CONSTANT IS PIE/180.
      AL     = ALPHA*1.745329E-02
      SINFAC = .5*(SIN(AL)/PHIEL)**2
      COSAL  = COS(AL)
      COSAL2 = COSAL**2
C DEBYE LENGTH CONSTANT FACTOR IS:
C          8.36753E-12 = 16.*PIE**2*PERMITIVITY*BOLTZK/COUL**2/C**2
      DEBYEC = 8.36753E-06*TF**2
C  SET UP FITTING CONTROL PARAMETERS:
      NPP    = 1
      DO 4 I = 1,7
      IP(I)  = NPP
      IFIT(I) = IFITR(I,IREGN)
      IF( IFIT(I) .NE. 1 ) GO TO 4
      NP     = NPP
      NPP    = NPP+1
      IP(I)  = NPP
      IFFIT(NP) = I
    4 CONTINUE
C
      NJMF   = 0
      IF( IFIT(5) .LE. 2 ) NJMF = 1
      IF( IFIT(6) .LE. 2 ) NJMF = 2
      NJ     = IFITR(8,IREGN)
      IF( NJ      .EQ. 1 ) NJMF = 0
C
C NUMERICAL CONSTANT IN THE FOLLOWING STATEMENT IS:
C         42.86445=SQRT( (M(PROTON)+M(NEUTRON))/M(ELECTRON)/2.0)
      SQMIN  = 428.6445
      DO 5 J = 1,NJ
      SQMI(J) = SQRT(AM(J,IREGN))*42.86445
      IF( SQMI(J) .LT. SQMIN ) SQMIN = SQMI(J)
   5  CONTINUE
      RETURN
C ======================================================================
      ENTRY SPECTR
C ----------------------------------------------------------------------
C DEFINE VARIOUS FACTORS WITH REPETITIVE USAGE:
      IF( IFIT(2) .EQ. 2 ) TE = TI*TR
      SQTE   = SQRT(TE)
      SQTI   = SQRT(TI)
      R5TI   = 0.5/TI
C "DEBYE"=NORMALIZED DEBYE LENGTH CORRECTION FACTOR.
      DEBYE  = DEBYEC/EN
      TWODEB = 2.*DEBYE
      TWOBTE = 2./TE
      TWOBTI = 2./TI
      NJM    = NJ-1
C SET UP FRACTIONAL COMPOSITION FACTORS:
      FRCTN1 = 1.0
      IF( NJ .NE. 1 ) THEN
      DO 10 JM = 1,NJM
   10 FRCTN1 = FRCTN1-FRACTN(JM)
      ENDIF
      THETC  = FNORMC/SQTI
      SCALEF = THETC*SQMI(1)
C "CONST"=SPECTRUM SCALING FACTOR (PROPORTIONAL TO "EN")
      CONST  = SYSTEM*EN
C "TEX"=APPROXIMATION OF "EXP(-TE*SINFAC)" WHICH BREAKS DOWN FOR
C LARGE "TF" AND "TE".
      TEX    = 1.-TE*SINFAC
      TEFAC  = 0.5/TE+SINFAC/TEX
      TEXBTE = TEX/TE
      DPSIDC = FNORMC*SQMI(1)/SQTI
      IF( CF .EQ. 0.0 ) GO TO 50
C ======================================================================
C  I O N  A D M I T A N C E  FUNCTIONS INCLUDING  C O L L I S I O N S.
C ASSUME "PSI" IS INDEPENDENT OF SPECIES MASS.
      PSIION = CF*DPSIDC
      PSII2  = PSIION*2.
      DPSIDT =-PSIION*R5TI
C COMPUTE COMPLEX PLASMA DISPERSION FUNCTION FOR FIRST ION:
      CALL PLASMA
      THEMAX = THETA(NDX)
      THEMIN = THEMAX
C COMPUTE COMPLEX ADMITANCE FUNCTION AND ITS DERIVATIVE FOR FIRST ION:
      PSIZ   = 1.-PSIION*ZI(1)
      PSIZM  = SQMI(1)/PSIZ
C MAKE "YIR1(1)" THE LIMIT OF "YIR/THETA" AS THETA GOES TO ZERO:
      YIR1(1) = PSIZM*ZI(1)
      YIR(1)  = FRCTN1*YIR1(1)
      DYRDPS = FRCTN1*PSIZM*(ZI(1)*(ZI(1)+PSII2)-2.)/PSIZ
      DYIRDT(1) = DYRDPS*DPSIDT
      DYIRDC(1) = DYRDPS*DPSIDC
      FRCFQ  = FRCTN1*DPSIDC
      DO 39 K = 2,NDX
      PSIZM  = PSIION*(ZR(K)**2+ZI(K)**2)
      YD     = 1.-PSIION*(2.*ZI(K)-PSIZM)
      YIR1(K) = THETA(K)*(ZI(K)-PSIZM)/YD
      YII1(K) = 1.+THETA(K)*ZR(K)/YD
      YIR(K) = FRCTN1*YIR1(K)
      YII(K) = FRCTN1*YII1(K)
C
      TH     = 1.0/THETA(K)
      C1     = YIR1(K)**2+YII1(K)*(1.-YII1(K))
      C2     = TH-2.*THETA(K)
      C3     = YIR1(K)*(1.-2.*YII1(K))
      C4     = YIR1(K)*C2
      C5     = PSII2*C1
      DYRDTH = C4-C5
      DYIDTH = PSII2*C3-TH+YII1(K)*C2
      DYRDPS = TH*C1+DYIDTH
      DYIDPS = C5-TH*C3-C4
      C1     =-THETA(K)*R5TI
C
      DYIRDT(K) = FRCTN1*(C1*DYRDTH+DPSIDT*DYRDPS)
      DYIIDT(K) = FRCTN1*(C1*DYIDTH+DPSIDT*DYIDPS)
      DYIRDC(K) = FRCFQ*DYRDPS
      DYIIDC(K) = FRCFQ*DYIDPS
C
   39 CONTINUE
      IF( NJ .EQ. 1 ) GO TO 70
C ----------------------------------------------------------------------
      DO 49 J = 2,NJ
      JM     = J-1
C COMPUTE COMPLEX PLASMA DISPERSION FUNCTION FOR OTHER IONS:
      SCALEF = THETC*SQMI(J)
      CALL PLASMA
      IF( THETA(NDX) .GT. THEMAX ) THEMAX = THETA(NDX)
      IF( THETA(NDX) .LT. THEMIN ) THEMIN = THETA(NDX)
C COMPUTE COMPLEX ADMITANCE AND DERIVATIVES FUNCTIONS FOR OTHER IONS:
      PSIZ   = 1.-PSIION*ZI(1)
      PSIZM  = SQMI(J)/PSIZ
C MAKE "YIRJ(1)" THE LIMIT OF "YIR/THETA" AS THETA GOES TO ZERO:
      YR     = PSIZM*ZI(1)
      YIR(1) = YIR(1)+FRACTN(JM)*YR
      DYRDPS = FRACTN(JM)*PSIZM*(ZI(1)*(ZI(1)+PSII2)-2.)/PSIZ
      DYIRDT(1) = DYIRDT(1)+DYRDPS*DPSIDT
      DYIRDC(1) = DYIRDC(1)+DYRDPS*DPSIDC
      DYIRDP(1,JM) = YR-YIR1(1)
      FRCFQ  = FRACTN(JM)*DPSIDC
      DO 49 K = 2,NDX
      PSIZM  = PSIION*(ZR(K)**2+ZI(K)**2)
      YD     = 1.-PSIION*(2.*ZI(K)-PSIZM)
      YR     = THETA(K)*(ZI(K)-PSIZM)/YD
      YI     = 1.+THETA(K)*ZR(K)/YD
      YIR(K) = YIR(K)+FRACTN(JM)*YR
      YII(K) = YII(K)+FRACTN(JM)*YI
C
      TH     = 1.0/THETA(K)
      C1     = YR**2+YI*(1.-YI)
      C2     = TH-2.*THETA(K)
      C3     = YR*(1.-2.*YI)
      C4     = YR*C2
      C5     = PSII2*C1
      DYRDTH = C4-C5
      DYIDTH = PSII2*C3-TH+YI*C2
      DYRDPS = TH*C1+DYIDTH
      DYIDPS = C5-TH*C3-C4
      C1     =-THETA(K)*R5TI
C
      DYIRDT(K) = DYIRDT(K)+FRACTN(JM)*(C1*DYRDTH+DPSIDT*DYRDPS)
      DYIIDT(K) = DYIIDT(K)+FRACTN(JM)*(C1*DYIDTH+DPSIDT*DYIDPS)
      DYIRDC(K) = DYIRDC(K)+FRCFQ*DYRDPS
      DYIIDC(K) = DYIIDC(K)+FRCFQ*DYIDPS
      DYIRDP(K,JM) = YR-YIR1(K)
      DYIIDP(K,JM) = YI-YII1(K)
   49 CONTINUE
      GO TO 70
C ======================================================================
C  I O N  A D M I T A N C E  FUNCTIONS WITHOUT  C O L L I S I O N S.
   50 CF     = 0.0
C COMPUTE COMPLEX PLASMA DISPERSION FUNCTION FOR FIRST ION:
      PSIION = 0.0
      CALL PLASMA
      THEMAX = THETA(NDX)
      THEMIN = THEMAX
C COMPUTE COMPLEX ADMITANCE FUNCTION AND ITS DERIVATIVE FOR FIRST ION:
      YIR1(1) = ZI(1)*SQMI(1)
      YIR(1)  = FRCTN1*YIR1(1)
      FRCFQ  = FRCTN1*DPSIDC
      FRCFQ2 = FRCFQ*2.
      DYIRDT(1) = 0.0
      DYIRDC(1) = FRCFQ*(ZI(1)*YIR1(1)-2.*SQMI(1))
      DO 59 K = 2,NDX
      YIR1(K) = THETA(K)*ZI(K)
      YII1(K) = THETA(K)*ZR(K)+1.0
      YIR(K) = FRCTN1*YIR1(K)
      YII(K) = FRCTN1*YII1(K)
C
      TH     = 1.0/THETA(K)
      IF( IFIT(4) .EQ. 1 ) THEN
      DYIRDC(K) = FRCFQ*((YIR1(K)**2+YII1(K)*(1.-YII1(K)))*TH
     >                           +YII1(K)*(TH-2.*THETA(K))-TH)
      DYIIDC(K)=FRCFQ2*YIR1(K)*(THETA(K)+(YII1(K)-1.)*TH)
      ENDIF
      TEMP   = THETA(K)**2/TI-R5TI
      DYIRDT(K) = FRCTN1*TEMP*YIR1(K)
      DYIIDT(K) = FRCTN1*(TEMP*YII1(K)+R5TI)
   59 CONTINUE
      IF( NJ .EQ. 1 ) GO TO 70
C ----------------------------------------------------------------------
C COMPUTE COMPLEX PLASMA DISPERSION FUNCTION FOR OTHER IONS:
      DO 69 J = 2,NJ
      IF( IFIT(J+3) .GT. 2 .AND. FRACTN(J-1) .EQ. 0.0 ) GO TO 69
      JM     = J-1
      SCALEF = THETC*SQMI(J)
      CALL PLASMA
      IF( THETA(NDX) .GT. THEMAX ) THEMAX = THETA(NDX)
      IF( THETA(NDX) .LT. THEMIN ) THEMIN = THETA(NDX)
C COMPUTE COMPLEX ADMITANCE AND DERIVATIVES FUNCTIONS FOR OTHER IONS:
      YR     = ZI(1)*SQMI(J)
      YIR(1) = YIR(1)+FRACTN(JM)*YR
      DYIRDP(1,JM) = YR-YIR1(1)
      FRCFQ  = FRACTN(JM)*DPSIDC
      FRCFQ2 = FRCFQ*2.
      DYIRDC(1) = DYIRDC(1)+FRCFQ*(ZI(1)*YR-2.*SQMI(J))
      DO 69 K = 2,NDX
      YR     = THETA(K)*ZI(K)
      YI     = THETA(K)*ZR(K)+1.0
      YIR(K) = YIR(K)+FRACTN(JM)*YR
      YII(K) = YII(K)+FRACTN(JM)*YI
C
      IF( IFIT(4) .EQ. 1 ) THEN
      TH     = 1.0/THETA(K)
      DYIRDC(K) = DYIRDC(K)+FRCFQ*((YR**2+YI*(1.-YI))*TH
     >                           +YI*(TH-2.*THETA(K))-TH)
      DYIIDC(K) = DYIIDC(K)+FRCFQ2*YR*(THETA(K)+(YI-1.)*TH)
      ENDIF
      TEMP   = THETA(K)**2/TI-R5TI
      DYIRDT(K) = DYIRDT(K)+FRACTN(JM)*TEMP*YR
      DYIIDT(K) = DYIIDT(K)+FRACTN(JM)*(TEMP*YI+R5TI)
C  COMPUTE PARTIALS W.R.T. FRACTIONAL ION COMPOSITIONS:
      DYIRDP(K,JM) = YR-YIR1(K)
      DYIIDP(K,JM) = YI-YII1(K)
   69 CONTINUE
C
C ======================================================================
C COMPUTE COMPLEX PLASMA DISPERSION FUNCTION FOR ELECTRONS:
   70 SCALEF = FNORMC/SQTE/COSAL
      CALL PLASMA
C ======================================================================
C CENTER FREQUENCY SPECTRUM AND PARTIAL DERIVATIVES.
      SQTIN  = SQTI*EN**2
      QTINE  = SQTIN*YIR(1)
      TINEQ  = QTINE*TI
      TICD   = TI*DEBYEC
      CDTINE = TICD+EN
      SQTEC  = CDTINE*SQTE*ZI(1)/COSAL2
      STECX  = SQTEC*TEX*TE
      YN     = TINEQ+STECX*CDTINE
      YM     = TI*EN+TE*CDTINE
      YN2YM  = 2.*YN/YM
      TEMP   = FNORMC*CONST/YM**2
C ----------------------------------------------------------------------
C THEORETICAL SPECTRUM AT "F=0.0":
      S(1)   = TEMP*YN
C  PARTIAL DERIVATIVES W.R.T.  T E M P E R A T U R E S  AT "F=0.0":
      DSDTE  = TEMP*CDTINE*(SQTEC*(1.5-2.5*SINFAC*TE)-YN2YM)
      DSDTI(1) = TEMP*(1.5*QTINE+2.*DEBYEC*STECX-YN2YM*(EN+TE*DEBYEC)+
     >           SQTIN*TI*DYIRDT(1))
      IF( IFIT(2) .EQ. 1 ) FUNC(1,ITE) = DSDTE
      IF( IFIT(2) .EQ. 2 ) FUNC(1,2)   = FUNC(1,2)+DSDTE*TR
C  PARTIAL DERIVATIVES W.R.T. ELECTRON  D E N S I T Y  AT "F=0.0":
      IF( IFIT(3) .EQ. 1 ) FUNC(1,INE) =
     >       TEMP*(3.*TINEQ+STECX*(TICD+3.*EN)-YN2YM*EN*(TI+TE))/EN
      TEMP   = TEMP*SQTIN*TI
C PARTIAL DERIVATIVES W.R.T.  C O L L I S I O N  FREQ. AT "F=0.0":
      IF( IFIT(4) .EQ. 1 ) FUNC(1,ICF) = TEMP*DYIRDC(1)
      NI     = ICF
      IF( NJMF .NE. 0 ) THEN
C  PARTIAL DERIVATIVES W.R.T. ION  C O M P O S T I O N  AT "F=0.0":
      DO 75 J = 1,NJMF
      NI = NI+1
   75 FUNC(1,NI) = TEMP*DYIRDP(1,J)
      ENDIF
      DO 85 I = 1,NI
      DO 85 L = 1,NLAGS
   85 ACF(L,I) = FUNC(1,I)*COSINE(L,1)
C ======================================================================
      KSMAX  = 1
C LOOP OVER ALL FREQUENCIES.
      DO 200 K = 2,NDX
C COMPUTE COMPLEX ADMITANCE FUNCTION AND ITS DERIVATIVE FOR ELECTRONS:
      THE2   = THETA(K)**2
      TEMP   = TEX*THETA(K)
      YER    = TEMP*ZI(K)
      YEI    = 1.+TEMP*ZR(K)
C
      TEMP   = THE2/TE-TEFAC
      DYERDT = TEMP*YER
      DYEIDT = TEMP*(YEI-1.)+THE2*TEXBTE
C ----------------------------------------------------------------------
      YERT   = YER/TE
      YEIT   = YEI/TE
      YEM    = YERT**2+YEIT**2
      YIRT   = YIR(K)/TI
      YMRT   = YERT+YIRT
      YIIT   = YII(K)/TI
      YIIH   = YIIT+DEBYE
      YIM    = YIRT**2+YIIH**2
      YMIT   = YEIT+YIIH
      YM     = YMRT**2+YMIT**2
      YN     = YEM*YIR(K)+YIM*YER
      YNBYM  = YN/YM
      TEMP   = CONST/F(K)
C ----------------------------------------------------------------------
C THEORETICAL POWER SPECTRUM:
      S(K)   = TEMP*YNBYM
      IF( S(K) .GT. S(KSMAX) ) KSMAX = K
C  INTERMEDIATE DERIVATIVES W.R.T. TEMPERATURES:
      DMYIR  = DYIRDT(K)-YIRT
      DMYII  = DYIIDT(K)-YIIT
      DYIMDT = TWOBTI*(YIRT*DMYIR+YIIH*DMYII)
      DYMDTI = TWOBTI*(YMRT*DMYIR+YMIT*DMYII)
      DMYER  = DYERDT-YERT
      DMYEI  = DYEIDT-YEIT
      DYEMDT = TWOBTE*(YERT*DMYER+YEIT*DMYEI)
      DYMDTE = TWOBTE*(YMRT*DMYER+YMIT*DMYEI)
C  PARTIAL DERIVATIVES OF SPECTRUM W.R.T.  T E M P E R A T U R E S :
      DSDTE  = (TEMP*(YIM*DYERDT+YIR(K)*DYEMDT)-S(K)*DYMDTE)/YM
      DSDTI(K) = (TEMP*(YEM*DYIRDT(K)+YER*DYIMDT)-S(K)*DYMDTI)/YM
C  REDIFINE TEMPERATURE DERIVATIVES DEPENDING ON FIT CODE:
      IF( IFIT(2) .EQ. 1 ) FUNC(K,ITE) = DSDTE
      IF( IFIT(2) .EQ. 2 ) FUNC(K,2) = FUNC(K,2)+DSDTE*TR
C  PARTIAL DERIVATIVES OF SPECTRUM W.R.T. ELECTRON  D E N S I T Y :
      IF( IFIT(3) .EQ. 1 ) FUNC(K,INE) = SYSTEM*(YNBYM+(YMIT*YNBYM-
     >                     YER*YIIH)*TWODEB/YM)/F(K)
C  PARTIAL DERIVATIVES OF SPECTRUM W.R.T. C O L L I S I O N  FREQUENCY:
      IF( IFIT(4) .EQ. 1 ) FUNC(K,ICF) = (TEMP*(YEM*DYIRDC(K)+
     >      YER*TWOBTI*(YIRT*DYIRDC(K)+YIIH*DYIIDC(K)))-
     >      S(K)*TWOBTI*(YMRT*DYIRDC(K)+YMIT*DYIIDC(K)))/YM
      NI     = ICF
      IF( NJMF .NE. 0 ) THEN
C  PARTIAL DERIVATIVES OF SPECTRUM W.R.T. ION  C O M P O S I T I O N :
      DO 60 J = 1,NJMF
      NI     = NI+1
   60 FUNC(K,NI) = TEMP*(DYIRDP(K,J)*(YEM+TWOBTI*(YER*YIRT-YMRT*YNBYM)
     >               )+DYIIDP(K,J)*TWOBTI*(YER*YIIH-YMIT*YNBYM))/YM
      ENDIF
      DO 180 I = 1,NI
      ACF(1,I) = ACF(1,I)+FUNC(K,I)*COSINE(1,K)
      DO 180 L = 2,NLAGS
  180 ACF(L,I) = ACF(L,I)+FUNC(K,I)*COSINE(L,K)
  200 CONTINUE
      IF( IFIT(7) .NE. 1 )  RETURN
      DO 210 K = 1,NLAGS
  210 ACF(K,NPP) = ACF(K,1)/SYSTEM
      RETURN
      END
C
C
C
      SUBROUTINE PLASMA
C*THIS ROUTINE COMPUTES THE COMPLEX PLASMA DISPERSION FUNCTION
C GIVEN BY:
C             Z(S)=I*SQRT(PIE)*EXPC(-S**2)*(1.+ERFC(I*S)
C WHERE:
C         I=SQRT(-1.) ;  S=X+I*Y=COMPLEX ARGUMENT
C FOR ABS(Y).GT.1.0, THE CONTINUED FRACTION EXPANSION GIVEN BY FRIED
C AND CONTE (1961) IS USED; WHILE FOR ABS(Y).LE.1.0, THE FOLLOWING
C DIFFERENTIAL EQUATION IS SOLVED:
C           D Z(S)
C           ------ = -2.*(1.+S*Z(S))
C            D S
C SUBJECT TO Z(0)=I*SQRT(PIE)
C
C     "F(K)"=TRUE FREQUENCY.
C     "X(K)"=NORMALIZED FREQUENCY.
C     "SCALEF"=FREQUENCY SCALING FACTOR FOR NORMALIZATION.
C ----------------------------------------------------------------------
C BY WES SWARTZ
C WHEN "Y" IS ZERO OR VERY SMALL, AND "X"IS GREATER THAN 7., THEN
C "ZI(K)" IS SET ZERO AND "ZR(K)" IS COMPUTED FROM THE ASYMPTOTIC
C EXPANSION FOR LARGE REAL ARGUMENTS.
C WHEN "Y.GT.0.0 .AND. Y.LE.1.0" THEN CODE IS GOOD FOR "ABS(X).LE.15."
C WHILE FOR "ABS(X).GT.15." RESULTS STILL LOOK GOOD WHEN ZI IS ZEROED,
C BUT NO DEFINITIVE CHECKS HAVE BEEN MADE.
C ======================================================================
      COMMON /PLAS/ZR(100),ZI(100),X(100),Y,F(120),SCALEF,NX
      SR=0.0
      SI=ABS(Y)
      IF (SI.GE.1.0) GO TO 8
C CODE BELOW USES SOLUTION TO DIFFERENTIAL EQUATION WHERE
C                  "(CR,CI)=FUNF(H,S,Z)=H*(1.+S*Z)"
C IS THE MAIN FUNCTION FOR INTEGRATING THE DIFFERENTIAL EQUATION FOR Z.
C CASE FOR "0.0.LT.ABS(Y).AND.ABS(Y).LE.1.0":
C "Z=(0.,1.772454)". SQRT(PI)=1.7724539
      ZZR=0.0
      ZZI=1.772454
C IS ARGUMENT PURELY REAL? IF SO, SKIP "Y" INTEGRATION.
      IF (Y.EQ.0.0) GO TO 4
      SI=0.0
      NY=INT(Y*50.+1.1)
      DY=Y/FLOAT(NY)
      H=2.*DY
      H2=.5*DY
      DO 2 K=1,NY
        A0I=DY*(1.-SI*ZZI)
        SI=SI+H2
        A1I=DY*(1.-SI*(ZZI-A0I))
        B0I=H*(1.-SI*(ZZI-A1I))
        SI=SI+H2
        B1I=H*(1.-SI*(ZZI-B0I))
    2   ZZI=ZZI-(2.*(A0I+A1I+A1I+B0I)+B1I)/6.
C GENERAL INTEGRATION OVER "X":
    4 ZR(1)=ZZR
      ZI(1)=ZZI
      X(1)=SR
      XMAX=SCALEF*F(NX)
      XLIMIT=XMAX
      IF (XLIMIT.GT.7.2 .AND. SI.LT.1.E-05) XLIMIT=7.2
      IF (XLIMIT.GT.15.0 .AND. SI.GE.1.E-05) XLIMIT=15.
      K=2
      X(K)=SCALEF*F(K)
    5 DX=X(K)-SR
      NSKIP=INT(DX*9.1+1.1)
      DX=DX/FLOAT(NSKIP)
      H=2.*DX
      H2=.5*DX
      DO 6 KSKIP=1,NSKIP
        A0R=DX*(1.+SR*ZZR-SI*ZZI)
        A0I=DX*(SR*ZZI+SI*ZZR)
        AR=ZZR-A0R
        AI=ZZI-A0I
        SR=SR+H2
        A1R=DX*(1.+SR*AR-SI*AI)
        A1I=DX*(SR*AI+SI*AR)
        AR=ZZR-A1R
        AI=ZZI-A1I
        B0R=H*(1.+SR*AR-SI*AI)
        B0I=H*(SR*AI+SI*AR)
        AR=ZZR-B0R
        AI=ZZI-B0I
        SR=SR+H2
        B1R=H*(1.+SR*AR-SI*AI)
        B1I=H*(SR*AI+SI*AR)
        ZZR=ZZR-(2.*(A0R+A1R+A1R+B0R)+B1R)/6.
    6   ZZI=ZZI-(2.*(A0I+A1I+A1I+B0I)+B1I)/6.
      ZR(K)=ZZR
      ZI(K)=ZZI
      SR=X(K)
      K=K+1
      IF (K.GT.NX) RETURN
      X(K)=SCALEF*F(K)
      IF (X(K).LE.XLIMIT) GO TO 5
C USE ASYMPTOTIC EXPANSION FOR LARGE REAL ARGUMENTS:
      IF (SI.GT.1.E-05) WRITE(6,22)
   22 FORMAT(' ** *** WARNING, X OUT OF RANGE IN Z(X,Y) *** **')
      KP=K
      DO 7 K=KP,NX
        X(K)=SCALEF*F(K)
        X2=X(K)**(-2)
        ZR(K)=-(1.+X2*(.5+X2*(.75+X2*(.625+X2*.4375))))/X(K)
        ZI(K)=0.0
    7 CONTINUE
      RETURN
C CODE BELOW IS FOR "ABS(Y).GT.1" USING CONTINUED FRACTION EXPANSION:
    8 DO 30 K=1,NX
        SR=SCALEF*F(K)
        X(K)=SR
        A0R=0.0
        A0I=0.0
        B0R=1.0
        B0I=0.0
        RR=SI**2-SR**2+0.5
        RI=-2.*SR*SI
        A1R=SR
        A1I=SI
        B1R=RR
        B1I=RI
        BM=B1R**2+B1I**2
        FR=(A1R*B1R+A1I*B1I)/BM
        FI=(A1I*B1R-A1R*B1I)/BM
        DO 10 N=1,23
          AN=FLOAT(N)
          C=AN*(.5-AN)
          DR=RR+2.*AN
          AR=DR*A1R+C*A0R-RI*A1I
          AI=DR*A1I+C*A0I+RI*A1R
          BR=DR*B1R+C*B0R-RI*B1I
          BI=DR*B1I+C*B0I+RI*B1R
          BM=BR**2+BI**2
          GR=(AR*BR+AI*BI)/BM
          GI=(AI*BR-AR*BI)/BM
          IF (ABS(GI/FI-1.).LE.1.E-05) GO TO 20
          FR=GR
          FI=GI
          A0R=A1R
          A0I=A1I
          A1R=AR
          A1I=AI
          B0R=B1R
          B0I=B1I
          B1R=BR
   10     B1I=BI
   20   ZR(K)=GR
        ZI(K)=GI
   30 CONTINUE
      IF (Y.GE.0.0) RETURN
C IF "Y.LT.1.0", THEN CHANGE QUADRANT:
      DO 40 K=1,NX
        SR=X(K)
        SI=Y
C FOR "Y.LT.-1.0" THE FOLLOWING CARD AND COMPLEX EXP IS REQUIRED.
C       S=(0.,3.544908)*EXPC(-S**2)
        ZR(K)=SR
   40   ZI(K)=SI
      RETURN
      END
