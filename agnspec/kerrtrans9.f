
       PROGRAM KERRTRANS9
c       IMPLICIT NONE
C This program calculates the polarization and spectrum of an accretion
C disk in a Kerr spacetime, for arbitrary spin parameter.  
C
C----------------------- Variable definitions ---------------------------
C
      INTEGER MRE,NRO,NPHI,NK0,MK,NW,NWO,MKO,NDISK
      PARAMETER (MRE = 50, NRO = 128, NPHI = 128, NK0=0, MK = 23286, 
     %           NW  = 10, NWO =   1, MKO  = 531, NDISK =1)
      DOUBLE PRECISION A,AEQ,ALPHA,BETA,DEL,
     %       EPSINU,GAM,L,L2,ONE,OMEGAG,OMEGA,
     %       PHI,PI,Q2,RMS,ROUT,UE,VG,Z1,Z2
      DOUBLE PRECISION AMIN,AMAX,BBB,BMIN,BMAX,C,CCC,DNUE,DNUO,E,ENERGY,
     &       ENU,ETA,FLUX,GE,GN,INTG,INU(MRE,MK,NW+2),INUR(MK,NW+2),
     &       ISPEC,
     &       LEDD,LUM,MASS,MDOT,MP,ME,MSUN,MUE,MUO(NWO),NU(MK),NUR(MK),
     &       P(MRE,MK,NW+2),PR(MK,NW+2),PSI,
     &       PSIM,PSIP,QQQ,R1,R2,R(MRE),RCUT,RE(MRE),RE0,RO(NRO),
     &       RG,MU(NW+2),SIGMAB,SIGMAT,SM,SSQRT,SU,
     &       T,T1,T2,TEFF(MRE),TMIN,U,V,WRO(NRO),WR(MRE),WT(NWO),
     &       X,X0,X1,X2,X3,
     &       MUHUB(NW),MUHUBWT(NW),NMISS(MKO,NWO),WA,WB,JUNK
      DOUBLE PRECISION A1,A2,B1,B2,DELTA,FCTISP,LSPEC(3,MKO),NUE,
     &       PANGLE,SFAC(3),TWO,WPHI,THIRD,WK(MKO),TEFF0,TRES,LTOT(3),
     &       LAM,HNU,NUO(MKO),NUO1,NUO2,ZERO,FOUR,FOURTH,KMAX,
     &       A1T,A2T,B1T,B2T,HC,HK,SMALL,TEX,WAT,WBT
      INTEGER GR
      LOGICAL TPR,TPM
      CHARACTER*40 RADFILE
      EXTERNAL PANGLE,FCTISP
C
C*************************************************************** 
C Define HK=h/k=Planck's constant/Boltzmann's constant
C Define HC=2*H/C^2
C*************************************************************** 
      PARAMETER(HK = 4.79916D-11, HC=1.47455D-47)
      PARAMETER ( ZERO=0.D0, ONE=1.D0, TWO=2.D0, TRES=3.D0, FOUR= 4.D0,
     &             THIRD = 0.333333333333333D0, FOURTH= 0.25D0,
     &             C=2.9979D10, MP=1.6726D-24, ME=9.1095D-28, 
     &             E=4.8032D-10, 
     &             MSUN=1.989D33, SIGMAB=5.6703D-5, GN=6.672D-8)
C
C-  Set up parameters:
C
      PI=ACOS(-ONE) 
      SIGMAT=(8.D0*PI/TRES)*(E*E/ME/C/C)**2
C----------------------- Initializations -------------------------------
C
C First, read in the frequencies to be used at the disk:
C
      READ(5,*) NK
      IF(NK.GT.0) THEN
        OPEN(UNIT=10,FILE='freqrad.in',STATUS='old')
        READ(10,*) (NU(KK),KK=NK,1,-1)
C     if NK < 0, then use logarithmically spaced points:
      ELSE
        NK=ABS(NK)
        READ(5,*) NUO1,NUO2
        NUO1=LOG(NUO1)
        NUO2=LOG(NUO2)
        DNUO=(NUO2-NUO1)/DBLE(NK)
        DO I=1,NK
          NU(I)=EXP(NUO1+DNUO*(DBLE(I)-0.5D0))
        ENDDO
      ENDIF
C
C-  Read in the frequencies observed and weights:
C
      READ(5,*) NKO
      IF(NKO.GT.0) THEN
C-      read in the frequencies observed and weights:
        OPEN(UNIT=12,FILE='freq.in',STATUS='old')
        DO I=NKO,1,-1
          READ(12,*) NUO(I),WK(I)
          IF(I.LE.NKO-2) WK(I+1)=0.5D0*(NUO(I+2)-NUO(I))
        ENDDO
        WK(1)=0.5D0*(NUO(2)-NUO(1))
        WK(NKO)=0.5D0*(NUO(NKO)-NUO(NKO-1))
        CLOSE(UNIT=12)
C     if NKO is negative, use logarithmically spaced frequencies:
      ELSE
        NKO=ABS(NKO)
        READ(5,*) NUO1,NUO2
        NUO1=LOG(NUO1)
        NUO2=LOG(NUO2)
        DNUO=(NUO2-NUO1)/DBLE(NKO)
        DO I=1,NKO
          NUO(I)=EXP(NUO1+DNUO*(DBLE(I)-0.5D0))
          WK (I)=NUO(I)*DNUO
        ENDDO
      ENDIF
C
C     the observation angle:
C
      READ(5,*) MUO(1)
      WT(1)=1.D0
C
C     basic parameters of a disk:
C
C     Input mass of hole in solar masses
C     mass acretion rate (msun/year)
C     a is the spin of the black hole in geometrized units:
C
      READ(5,*) MASS,MDOT,A
C
      MASS=MASS*MSUN
      MDOT=MDOT*MSUN/365.25D0/24.D0/3600.D0
C-    Eddington Luminosity:
      LEDD=FOUR*PI*GN*MASS*MP*C/SIGMAT
C-    Gravitational radius:
      RG=GN*MASS/C/C
C
C     Marginally stable (direct) orbit radius, in units of mass
C
      Z1=ONE+(ONE-A*A)**THIRD*((ONE+A)**THIRD+(ONE-A)**THIRD)
      Z2=SSQRT(TRES*A*A+Z1*Z1)
      RMS=TRES+Z2-SSQRT((TRES-Z1)*(TRES+Z1+TWO*Z2))
      X0=SSQRT(RMS)
      X1=TWO*COS((ACOS(A)-PI)*THIRD)
      X2=TWO*COS((ACOS(A)+PI)*THIRD)
      X3=-TWO*COS(ACOS(A)*THIRD)
C
C-    Efficiency and accretion rate:
      ENERGY=(RMS*RMS-TWO*RMS+A*SSQRT(RMS))/RMS/
     %       SQRT(RMS*RMS-TRES*RMS+TWO*A*SSQRT(RMS))
      ETA=ONE-ENERGY
      LUM=ETA*MDOT*C*C
C
C-------------- Emission radius- grid of integration -------------
C     ROUT is the outer boundary of the disk:
C
      READ(5,*) ROUT,TMIN,NRE,NRETYPE, GR
      RCUT=ROUT
      IF(NRETYPE.EQ.1) THEN
C-  Logarithmic spacing:
C-  Range of integration with respect to log(re/rms):
        R1 = ZERO
C-  Radiation for re > rms (stable disk only):
        R2 = LOG( RCUT / RMS )
C-  Abscissas and weights for re-integration:
        CALL GAULEG(R1,R2,R,WR,NRE,MRE)
        DO I = 1, NRE
           RE(I) = RMS * EXP( R(I) )
        END DO
      ELSE IF(NRETYPE.EQ.2) THEN
C   Inverse spacing:
        R1 = ONE / RCUT
        R2 = ONE / RMS
C-  Radiation for re > rplus (from horizon):
C       R2 = ONE / RPLUS
C-  Abscissas and weights for re-integration:
        CALL GAULEG(R1,R2,R,WR,NRE,MRE)
        DO I=1,NRE
           RE(I) = ONE / R(I)
        END DO
       ELSE IF(NRETYPE.EQ.3) THEN
C   Third option is to read in the emission radii from a file:
C   The radii should be in ascending order
        OPEN(UNIT=13,FILE='emrad.in',STATUS='old')
        READ(13,*) (RE(I),I=1,NRE)
        iunin=11
       ELSE
        OPEN(UNIT=13,FILE='emrad.in',STATUS='old')
        iunin=13
      ENDIF
C
C-------------- Observation radius - grid of integration -----------
C
      NROTYPE=4
      A1=-ROUT
      A2=ROUT
      B1=-ROUT
      B2=ROUT
      WA=(A2-A1)/DBLE(NRO)
      WB=(B2-B1)/DBLE(NPHI)
      SMALL = RMS*0.05
      NS=MAX(INT(LOG((A2-A1)/SMALL/DBLE(NRO))/LOG(2.))+1,
     &   INT(LOG((B2-B1)/SMALL/DBLE(NPHI))/LOG(2.))+1)
c
c ******************************************************************
C--------------- Input atmosphere flux & polarization data ---------------
C
      DO 50 I = 1, NRE
         IF(NRETYPE.GT.3) READ(13,*) RE(I),NKR
C-  Calculate relativistic factors from Novikov & Thorne:
         BBB=ONE+A/RE(I)**1.5D0
         CCC=ONE-TRES/RE(I)+TWO*A/RE(I)**1.5D0
         X=SSQRT(RE(I))
         QQQ=(ONE+A/X**TRES)/SSQRT(ONE-TRES/X**TWO+TWO*A/X**TRES)/X*
     %       (X-X0-1.5D0*A*LOG(X/X0)-TRES*(X1-A)**TWO/X1/
     %       (X1-X2)/(X1-X3)*LOG((X-X1)/(X0-X1))-TRES*(X2-A)**TWO/
     %       X2/(X2-X1)/(X2-X3)*LOG((X-X2)/(X0-X2))-TRES*(X3-A)
     %       **TWO/X3/(X3-X1)/(X3-X2)*LOG((X-X3)/(X0-X3)))
C-  Calculate flux and effective temperature:
         FLUX=TRES*GN*MASS*MDOT*QQQ/8.D0/PI/RE(I)**TRES/
     %        BBB/SSQRT(CCC)/RG**TRES
         TEFF(I)=(FLUX/SIGMAB)**FOURTH
         IF(TEFF(I).GT.TMIN) THEN
C   Input the name of the file for the current radius:
           IF(NRETYPE.LE.3) THEN
              READ(13,*) RADFILE,NKR
              OPEN(UNIT=11,FILE=RADFILE,STATUS='old')
           END IF
C-----------------------------------------------------------------------------
C-  First, read the flux and polarization into an array
C-  NOTE: the index for mu ranges from 0 to nw+1, while
C-  the index for nu ranges from 1 to nk.
C            
           CALL GAULEG(ZERO,ONE,MUHUB,MUHUBWT,NW,NW)
           MU(1)=0.D0
           DO MM=2,NW+1
              MU(MM)=MUHUB(MM-1)
           END DO
           MU(NW+2)=1.D0
C  Hubeny has file ordered by increasing lambda, so reverse the reading:
           DO 10 KK = NKR, 1, -1
               READ(iunin,*) LAM, HNU
               READ(iunin,*) (INUR(KK,MM),PR(KK,MM),MM=2,NW+1)
C  Convert lambda (Angstroms) to frequency (Hertz)
               NUR(KK)=C/LAM*1.D8
C-  Linearly extrapolate to mu=0,1:
               INUR(KK,1)=INUR(KK,3)+MU(3)*(INUR(KK,2)-INUR(KK,3))/
     &                    (MU(3)-MU(2))
               INUR(KK,NW+2)=INUR(KK,NW)-(ONE-MU(NW))*
     &                       (INUR(KK,NW)-INUR(KK,NW+1))/
     &                       (MU(NW+1)-MU(NW))
               PR(KK,1)=PR(KK,3)+MU(3)*(PR(KK,2)-PR(KK,3))/(MU(3)-MU(2))
               PR(KK,NW+2)=0.D0
               DO MM=1,NW+2 
                  IF(INUR(KK,MM).LT.0.D0) INUR(KK,MM) = 0.D0
               ENDDO
 10        CONTINUE
C Now, we need to interpolate from the NUR grid to the NU grid:
           DO 20 K=1,NK
              CALL LOCATE(NUR,NKR,NU(K),KK,MK)
              DO MM=1,NW+2
                 IF(KK.EQ.0) THEN
                   TEX=LOG(HC/INUR(1,MM)*NUR(1)**3+1.D0)
                   IF(TEX.GT.0.D0) THEN
                     TEX=HK*NUR(1)/TEX
                    ELSE
                     TEX=HK/NUR(1)**2*INUR(1,MM)/HC
                   ENDIF
                   INU(I,K,MM)=FCTISP(NU(K),MU(MM),DELTA,TEX)
                   P(I,K,MM)=0.D0
                  ELSE IF (KK.EQ.NKR) THEN
                   TEX=LOG(HC/INUR(NKR,MM)*NUR(NKR)**3+1.D0)
                   IF(TEX.GT.0.D0) THEN
                     TEX=HK*NUR(NKR)/TEX
                    ELSE
                     TEX=HK/NUR(NKR)**2*INUR(NKR,MM)/HC
                   ENDIF
                   INU(I,K,MM)=FCTISP(NU(K),MU(MM),DELTA,TEX)
                   P(I,K,MM)=0.D0
                  ELSE
                   T=(NU(K)-NUR(KK))/(NUR(KK+1)-NUR(KK))
                   INU(I,K,MM)=INUR(KK,MM)*(1.D0-T)+INUR(KK+1,MM)*T
                   P(I,K,MM)=PR(KK,MM)*(1.D0-T)+PR(KK+1,MM)*T
                 ENDIF
              END DO
 20        CONTINUE
c ******************************************************************
          IF(NU(1).GT.NU(NK)) WRITE(6,*) 'Frequencies in wrong order'
c ******************************************************************
          IF(NRETYPE.LE.3) CLOSE(UNIT=iunin)
        ENDIF
   50 CONTINUE
      IF(NRETYPE.GT.3) CLOSE(UNIT=iunin)
C--------------------------------------------------------
C-  Set the total flux to zero:
      DO M=1,3
         LTOT(M)=ZERO
      END DO
c ******************************************************************
      II=1
      A1T=A1
      B1T=B1
      A2T=A2
      B2T=B2
      WAT=WA
      WBT=WB
c ******************************************************************
C     
C-  Initialize integral over frequencies:
      DO K=1,NKO
         DO J=1,3
            LSPEC(J,K) = ZERO
         END DO
      END DO
      AMIN = 0.D0
      BMIN = 0.D0
      AMAX = 0.D0
      BMAX = 0.D0
C
C-------------- Calculation of the spec. luminosity -----
c ******************************************************************
C Loop over nested grids:
      DO 100 NG = 1,NS
c ******************************************************************
C-  Integration with respect to ro:
         DO 90 I = 1, NRO
c ******************************************************************
            IF(NROTYPE.GT.1) ALPHA=A1T+(A2T-A1T)*(DBLE(I-1)+0.5D0)/
     &                             DBLE(NRO)
c ******************************************************************
C-  Set phi-integral to 0:
C-  Integration with respect to phi:
            DO 70 J = 1, NPHI
               IF(NROTYPE.EQ.1) THEN
                 IF(NPHI.NE.1.D0) THEN 
                   PHI=TWO*PI*DBLE(J-1)/DBLE(NPHI-1)
                  ELSE
                   PHI=0.D0
                 ENDIF
C Calculate impact parameters at infinity:
                 ALPHA=RO(I)*SIN(PHI)
                 BETA=-RO(I)*COS(PHI)
                ELSE
                 BETA=B1T+(B2T-B1T)*(DBLE(J-1)+0.5D0)/DBLE(NPHI)
               ENDIF
c ******************************************************************
C Make a hole in the center, which will be filled in with finer grids:
               IF(NG.LT.NS.AND.ABS(ALPHA).LT.A2T*0.5.AND.
     &           ABS(BETA).LT.B2T*0.5)  GOTO 70
c ******************************************************************
C Full general relativistic calculation?
               IF(GR.EQ.1) THEN
C Calculate the angular momentum and Carter's constant of motion:
               L=-ALPHA*SSQRT(ONE-MUO(II)*MUO(II))
               L2=ALPHA*ALPHA*(ONE-MUO(II)*MUO(II))
               Q2=BETA*BETA-(A*A-ALPHA*ALPHA)*MUO(II)*MUO(II)
C-  Now calculate the emission radius:
C   (If beta>0, there is a mu turning point; otherwise not)
               IF(BETA.GT.ZERO) THEN
                 TPM=.TRUE.
                 SM=ONE
                ELSE
                 TPM=.FALSE.
                 SM=-ONE
               ENDIF
C  The sign of the U integral starts out positive:
               SU=ONE
               IF(ABS(Q2).LT.1.D-16) Q2=ZERO
C  Calculate the emission radius (UE=1/RE)
               CALL GEOR3(ZERO,UE,MUO(II),ZERO,A,L,L2,Q2,TPM,TPR,SU,SM)
C  ....otherwise Newtonian version of the calculation:
               ELSE
                 UE=ONE/SQRT(ALPHA*ALPHA+BETA*BETA/MUO(II)/MUO(II))
               ENDIF
C  Check to see if the radius is within range:
               IF(UE.GT.ONE/ROUT.AND.UE.LE.ONE/RMS) THEN
                 RE0=ONE/UE
                 AMIN = MIN(AMIN,ALPHA)
                 BMIN = MIN(BMIN,BETA)
                 AMAX = MAX(AMAX,ALPHA)
                 BMAX = MAX(BMAX,BETA)
                ELSE
                 IF(NROTYPE.EQ.3) THEN
                   DO K=1,NKO
                      WRITE(800+K,200) ALPHA,BETA,0
                   ENDDO
                 ENDIF
                 GOTO 70
               ENDIF
C  Now, calculate the metric and relativisitic factors used in
C  computing the redshift and angle of emission:
               IF(GR.EQ.1) THEN
               DEL=RE0*RE0+A*A-TWO*RE0
               AEQ=(RE0*RE0+A*A)**TWO-A*A*DEL
               EPSINU=AEQ*UE*UE/SSQRT(DEL)
               ENU=SSQRT(DEL/AEQ)*RE0
C  The rotational frequency, velocity, and Lorentz factor of the gas:
               OMEGAG=ONE/(RE0*SSQRT(RE0)+A)
               VG=(RE0*RE0-TWO*A*SSQRT(RE0)+A*A)/SSQRT(DEL)*OMEGAG
               GAM=ONE/SSQRT(ONE-VG*VG)
               OMEGA=TWO*A*RE0/AEQ
C  Redshift of the emission:
               GE=ENU/GAM/(ONE-OMEGAG*L)
C  Angle of emission (cosine):
               MUE=SSQRT(Q2)/RE0*ENU/GAM/(ONE-OMEGAG*L)
C ....otherwise, compute Newtonian disk:
               ELSE
                 MUE=MUO(II)
                 GE=ONE
               ENDIF
               CALL LOCATE(MU,NW+2,MUE,MM,NW+2)
               U=(MUE-MU(MM))/(MU(MM+1)-MU(MM))
C  Find where the emission radius lies on the grid of calculated radii:
              IF(RE0.GT.MIN(RE(1),RE(NRE)).AND.
     &          RE0.LT.MAX(RE(1),RE(NRE))) THEN
                CALL LOCATE(RE,NRE,RE0,JJ,NRE)
                V=(RE0-RE(JJ))/(RE(JJ+1)-RE(JJ))
               ELSE
                JJ = 0
              ENDIF
C-  Calculate the emission+rotation angle PSI:
              IF(GR.EQ.1) THEN
              PSIP=PANGLE(RE0,TPR,L,SSQRT(Q2),BETA,ALPHA,A,ONE,MUO(II))
              PSIM=PANGLE(RE0,TPR,L,SSQRT(Q2),BETA,ALPHA,A,-ONE,MUO(II))
              ELSE
                PSIP=0.d0
                PSIM=0.d0
              ENDIF
C-  Loop over observed frequency points:
              DO 68 K=1,NKO
C-  Calculation of the emitted frequency:
                 NUE = NUO(K) / GE
C-  Calculation of the emitted specific intensity Ispec and polarization:
C   (I have a minimum temperature cutoff since my atmospheres don't converge
C    for low temperatures)
                 IF(JJ.GT.0.AND.TEFF(JJ).GT.TMIN) THEN
                   IF(NUE.LT.MIN(NU(1),NU(NK)).OR.
     &               NUE.GT.MAX(NU(1),NU(NK))) THEN
C*************************************************************** 
C First, check if the frequency is too low (NOTE: I'm assuming that NU(1)
C contains the lowest frequency and NU(NK) the highest, so beware if sorted 
C by increasing wavelength):
                       IF(NUE.LT.NU(1)) THEN
C Find the intensity at the lowest frequency:
                       ISPEC=(ONE-V)*((ONE-U)*INU(JJ,1,MM)+
     &                        U*INU(JJ,1,MM+1)) +
     &                        V*((ONE-U)*INU(JJ+1,1,MM)+
     &                        U*INU(JJ+1,1,MM+1))
                       IF(ISPEC.GT.0.D0) THEN
C Compute temperature of blackbody with same frequency and intensity:
                        TEX=LOG(HC/ISPEC*NU(1)**3+1.D0)
                        IF(TEX.GT.0.D0) THEN
                          TEX=HK*NU(1)/TEX
                         ELSE
                          TEX=HK/NU(1)**2*ISPEC/HC
C Extrapolate assuming a blackbody holds at lower frequencies:
                        ENDIF
                        ISPEC=FCTISP(NUE,MUE,DELTA,TEX)
                      ENDIF
C Assume the polarization is constant below NU(1):
                      DELTA=(ONE-V)*((ONE-U)*P(JJ,1,MM)+
     &                       U*P(JJ,1,MM+1)) +
     &                       V*((ONE-U)*P(JJ+1,1,MM)+
     &                       U*P(JJ+1,1,MM+1))
                     ELSE IF(NUE.GT.NU(NK)) THEN
C Find the intensity at the highest frequency:
                    ISPEC=(ONE-V)*((ONE-U)*INU(JJ, NK,MM)+
     &                    U*INU(JJ, NK,MM+1)) +
     &                    V*((ONE-U)*INU(JJ+1,NK,MM)+
     &                    U*INU(JJ+1,NK,MM+1))
                    IF(ISPEC.GT.0.D0) THEN
C Compute temperature of blackbody with same frequency and intensity:
                      TEX=LOG(HC/ISPEC*NU(NK)**3+1.D0)
                      IF(TEX.GT.0.D0) THEN
                        TEX=HK*NU(NK)/TEX
C Extrapolate assuming a blackbody holds at higher frequencies:
                      ELSE
                        TEX=HK/NU(NK)**2*ISPEC/HC
                      ENDIF
                      ISPEC=FCTISP(NUE,MUE,DELTA,TEX)
                    ENDIF   
C Assume the polarization is constant above NU(NK):
                    DELTA=(ONE-V)*((ONE-U)*P(JJ, NK,MM)+
     &                     U*P(JJ, NK,MM+1)) +
     &                     V*((ONE-U)*P(JJ+1,NK,MM)+
     &                     U*P(JJ+1,NK,MM+1))
                  ENDIF
C*************************************************************** 
                 ELSE
C-  Tri-linear interpolation in radius, frequency, and angle:
                  CALL LOCATE(NU,NK,NUE,KK,MK)
                  T=(NUE-NU(KK))/(NU(KK+1)-NU(KK))
                  ISPEC=(ONE-V)*((ONE-U)*((ONE-T)*INU(JJ,KK,MM)+T*
     %                    INU(JJ,KK+1,MM))+T*U*INU(JJ,KK+1,MM+1)+
     %                    (ONE-T)*U*INU(JJ,KK,MM+1))+
     %                    V*((ONE-U)*((ONE-T)*
     %                    INU(JJ+1,KK,MM)+T*INU(JJ+1,KK+1,MM))+
     %                    T*U*INU(JJ+1,KK+1,MM+1)+
     %                    (ONE-T)*U*INU(JJ+1,KK,MM+1))
                  DELTA=(ONE-V)*((ONE-U)*((ONE-T)*P(JJ,KK,MM)+T*
     %                    P(JJ,KK+1,MM))+T*U*P(JJ,KK+1,MM+1)+
     %                    (ONE-T)*U*P(JJ,KK,MM+1))+V*((ONE-U)*((ONE-T)*
     %                    P(JJ+1,KK,MM)+T*P(JJ+1,KK+1,MM))+
     %                    T*U*P(JJ+1,KK+1,MM+1)+
     %                    (ONE-T)*U*P(JJ+1,KK,MM+1))
                ENDIF
              ELSE
C- Just assume an isotropic blackbody and zero polarization for low Teff:
                BBB=ONE+A/RE0**1.5D0
                CCC=ONE-TRES/RE0+TWO*A/RE0**1.5D0
                X=SSQRT(RE0)
                QQQ=(ONE+A/X**TRES)/SSQRT(ONE-TRES/X**TWO+
     %              TWO*A/X**TRES)/X*
     %              (X-X0-1.5D0*A*LOG(X/X0)-TRES*(X1-A)**TWO/X1/
     %              (X1-X2)/(X1-X3)*LOG((X-X1)/(X0-X1))-
     %              TRES*(X2-A)**TWO/
     %              X2/(X2-X1)/(X2-X3)*LOG((X-X2)/(X0-X2))-TRES*(X3-A)
     %              **TWO/X3/(X3-X1)/(X3-X2)*LOG((X-X3)/(X0-X3)))
C-  Calculate flux and effective temperature:
                FLUX=TRES*GN*MASS*MDOT*QQQ/8.D0/PI/RE0**TRES/BBB/
     %               SSQRT(CCC)/RG**TRES
                TEFF0=(FLUX/SIGMAB)**FOURTH
		IF(RE0.LT.RE(1)) THEN
		  ISPEC = 0.D0
		ELSE
                  ISPEC= FCTISP(NUE,MUE,DELTA,TEFF0)
		ENDIF
              ENDIF
              IF(DELTA.GT.ZERO) PSI=PSIP
              IF(DELTA.LT.ZERO) PSI=PSIM
              SFAC(1)=ONE
              SFAC(2)=DELTA*COS(TWO*PSI)
              SFAC(3)=DELTA*SIN(TWO*PSI)
C-  Calculation of the integrand with respect to phi:
              IF(NROTYPE.EQ.1) THEN
                INTG = FOUR*PI*RO(I)*(RO(I)+RCUT)*ISPEC*GE**3*
     &                 WPHI*WRO(I)
              ELSE
c ******************************************************************
                INTG = FOUR*PI*ISPEC*GE**3*WAT*WBT
c ******************************************************************
              ENDIF
              DO M=1,3
                LSPEC(M,K) = LSPEC(M,K) + INTG*SFAC(M)*RG*RG
              END DO
              IF(NROTYPE.EQ.3) THEN
                WRITE(800+K,200) ALPHA,BETA,INTG
              ENDIF
 68         CONTINUE
 70         CONTINUE
 90      CONTINUE
c ******************************************************************
         A1T=0.5*A1T
         A2T=0.5*A2T
         B1T=0.5*B1T
         B2T=0.5*B2T
         WAT=0.5*WAT
         WBT=0.5*WBT
  100 CONTINUE
c ******************************************************************
C---------------------------------------------------------
C-------------- Output of the spec. luminosity -----------
c ***************************************************************************
c I think Eric's LSPEC(1,K) is 4*pi*d^2*F_nu, where F_nu is the flux (from
c one side of the disk) that an observer would measure at the particular
c viewing angle theta, and d is the distance.
c ***************************************************************************
      DO K=1,NKO
         WRITE(6,200) NUO(K),LSPEC(1,K),LSPEC(2,K),LSPEC(3,K),
     &                NMISS(K,II)
         DO M=1,3
            LTOT(M)=LTOT(M)+WT(II)*WK(K)*LSPEC(M,K)
         END DO
      END DO
      CALL FLUSH(6)
C---------------------------------------------------------
 200  FORMAT(10(1x,1pe13.5))
 201  FORMAT(1x,1pe13.5)
      END
C
C
C
C========================================================================
C
C
      double precision function FCTISP(NUE,MU,DELTA,T)
      implicit none
      double precision a,c,delta,h,k,mu,nue,t
      h=6.6262d-27
      k=1.3807d-16
      c=2.9979d10
      a=0.6d0
      if(t.eq.0.d0) then
        fctisp = 1.d0
        return
      endif
      if(h*nue/k/t.gt.200.d0) then
        fctisp=0.d0
      else
        fctisp=2.d0*h*nue**3/c/c/(exp(h*nue/k/t)-1.d0)
      endif
      return
      end

C
C
C========================================================================
C
C
       
       double precision FUNCTION ellf(tanphi,m1)
       implicit none
       double precision ef,ak,phi,m1,tanphi
C       USES rf
       double precision s,rf
       phi=atan2(abs(tanphi),sign(1.d0,tanphi))
       ak=sqrt(1.d0-m1)
       s=sin(phi)
       ef=s*rf(1.d0-s*s,(1.d0-s*ak)*(1.d0+s*ak),1.d0)
       if(tanphi.lt.0.d0) then
         ellf=2.d0*rf(0.d0,m1,1.d0)-ef
       else
         ellf=ef
       endif
       return
       END

C
C
C========================================================================
C
C
       double precision FUNCTION rf(x,y,z)
       implicit none
       double precision x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
       PARAMETER(ERRTOL=0.0025d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *              C1=1.d0/24.d0,C2=0.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
       double precision alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,
     *              sqrtz,ssqrt,xt,yt,zt
       if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.
     *              max(x,y,z).gt.BIG)pause 'invalid arguments in rf'
       xt=x
       yt=y
       zt=z
1       continue
              sqrtx=ssqrt(xt)
              sqrty=ssqrt(yt)
              sqrtz=ssqrt(zt)
              alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
              xt=0.25d0*(xt+alamb)
              yt=0.25d0*(yt+alamb)
              zt=0.25d0*(zt+alamb)
              ave=THIRD*(xt+yt+zt)
              if(ave.eq.0.d0) then
                delx=0.d0
                dely=0.d0
                delz=0.d0
              else
                delx=(ave-xt)/ave
                dely=(ave-yt)/ave
                delz=(ave-zt)/ave
              endif
       if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)go to 1
       e2=delx*dely-delz*delz
       e3=delx*dely*delz
       rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/ssqrt(ave)
       return
       END
C
C
C========================================================================
C
C

      SUBROUTINE SNCNDN(UU,EMMC,SN,CN,DN)
       implicit none
       integer i,ii,l
       double precision a,b,c,ca,cn,d,dn,emc,emmc,sn,ssqrt,u,uu
      PARAMETER (CA=3.d-8)
      LOGICAL BO
      double precision EM(13),EN(13)
      EMC=EMMC
      U=UU
      IF(EMC.NE.0.d0)THEN
        BO=(EMC.LT.0.d0)
        IF(BO)THEN
          D=1.d0-EMC
          EMC=-EMC/D
          D=SSQRT(D)
          U=D*U
        ENDIF
        A=1.d0
        DN=1.d0
        DO 11 I=1,13
          L=I
          EM(I)=A
          EMC=SSQRT(EMC)
          EN(I)=EMC
          C=0.5d0*(A+EMC)
          IF(ABS(A-EMC).LE.CA*A)GO TO 1
          EMC=A*EMC
          A=C
11      CONTINUE
1       U=C*U
        SN=DSIN(U)
        CN=DCOS(U)
        IF(SN.EQ.0.)GO TO 2
        A=CN/SN
        C=A*C
        DO 12 II=L,1,-1
          B=EM(II)
          A=C*A
          C=DN*C
          DN=(EN(II)+A)/(B+A)
          A=C/B
12      CONTINUE
        A=1.d0/SSQRT(C*C+1.d0)
        IF(SN.LT.0.)THEN
          SN=-A
        ELSE
          SN=A
        ENDIF
        CN=C*SN
2       IF(BO)THEN
          A=DN
          DN=CN
          CN=A
          SN=SN/D
        ENDIF
      ELSE
        CN=1.d0/DCOSH(U)
        DN=CN
        SN=DTANH(U)
      ENDIF
      RETURN
      END
       
C
C
***********************************************************************
      SUBROUTINE LAGUER(A,M,X,ITS)
***********************************************************************
*     PURPOSE:  Find one root of a polynomial.
*     ARGUMENTS:
*     ROUTINES CALLED:
*     ALGORITHM:  
*     ACCURACY:
*     REMARKS:  I don't have the documentation for this routine!
*     AUTHOR:  Numerical Recipes (new edition).
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
       IMPLICIT NONE
      INTEGER ITS,M
      COMPLEX*16 A(M+1),X
*       
      INTEGER ITER,J,MAXIT,MR,MT
      PARAMETER (MR=8,MT=10,MAXIT=MT*MR)
       DOUBLE PRECISION ABX,ABP,ABM,ERR,EPSS,FRAC(MR) 
       PARAMETER (EPSS=1.d-15)
      COMPLEX*16 DX,X1,B,D,F,G,H,SQ,GP,GM,G2
      DATA FRAC /0.5d0,0.25d0,0.75d0,0.13d0,0.38d0,0.62d0,0.88d0,1.d0/
*       Loop over iterations up to allowed maximum.
      DO 20 ITER=1,MAXIT
*         
        ITS=ITER
        B=A(M+1)
        ERR=ABS(B)
        D=DCMPLX(0.d0,0.d0)
        F=DCMPLX(0.d0,0.d0)
        ABX=ABS(X)
        DO 10 J=M,1,-1
*           Efficient computation of the polynomial and its first TWO
*           derivatives.
          F=X*F+D
          D=X*D+B
          B=X*B+A(J)
          ERR=ABS(B)+ABX*ERR
   10   CONTINUE
        ERR=EPSS*ERR
*         
        IF(ABS(B).LE.ERR) THEN
*           Special case: we are on the root.
          RETURN
        ELSE
*           The generic case; use Laguerre's formula.
          G=D/B
          G2=G*G
          H=G2-2.d0*F/B
          SQ=SQRT((M-1)*(M*H-G2))
          GP=G+SQ
          GM=G-SQ
          ABP=ABS(GP)
          ABM=ABS(GM)
          IF(ABP.LT.ABM) GP=GM
          IF (MAX(ABP,ABM).GT.0.d0) THEN
            DX=M/GP
          ELSE
            DX=CEXP(CMPLX(LOG(1.d0+ABX),DBLE(ITER)))
          ENDIF
        ENDIF
        X1=X-DX
*         Check if we've converged.
        IF(X.EQ.X1)RETURN
        IF (MOD(ITER,MT).NE.0) THEN
          X=X1
        ELSE
*           
          X=X-DX*FRAC(ITER/MT)
        ENDIF
   20 CONTINUE
      PAUSE 'Too many iterations'
      RETURN
      END
C
C
***********************************************************************
      SUBROUTINE ZROOTS(A,M,ROOTS,POLISH)
***********************************************************************
*     PURPOSE:  Find all roots of a polynomial.
*     ARGUMENTS:  Given the degree M and the M+1 complex coefficients
*       A of the polynomial (with A(0) being the constant term), this
*       routine returns all M roots in the complex array ROOTS.  The
*       logical variable POLISH should be input as .TRUE. if polishing
*       (by Laguerre's method) is desired, .FALSE. if the roots will be
*       subsequently polished by other means.
*     ROUTINES CALLED:  LAGUER.
*     ALGORITHM: Laguerre's method.
*     ACCURACY:  The parameter EPS sets the desired accuracy.
*     REMARKS:  
*     AUTHOR:  Numerical Recipes (new edition).
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
      INTEGER M
      COMPLEX*16 A(5),ROOTS(4)
      LOGICAL POLISH
*       
      INTEGER I,J,JJ,ITS,MAXM,ncc
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.d-6,MAXM=5)
      COMPLEX*16 AD(5),X,B,C,XSUM,XSTART
*       
      IF(M.GT.MAXM-1) THEN
        PAUSE 'M too large in ZROOTS'
      ENDIF
*       Copy of coefficients for successive deflation.
      DO 10 J=1,M+1
        AD(J)=A(J)
   10 CONTINUE
*       Loop over each root to be found.
      XSUM=0.d0
       XSTART=DCMPLX(0.D0,0.D0)
* If ncc=1, the previous root is a complex conjugate of another root
       ncc=0
      DO 20 J=M,1,-1
* Start at zero to favour convergence to smallest remaining
* root, or if the previous root was complex, start at its complex
* conjugate (since the coefficients are real):
        X=XSTART
         if(J.lt.M.and.dimag(ROOTS(J+1)).ne.0.d0.and.ncc.eq.0) then
           XSTART=DCMPLX(dble(roots(J+1)),-dble(roots(J+1)))
* Since we have chosen the second root to start at the complex conjugate,
* we don't want to use its complex conjugate again as a starting root:
           ncc=1
         else
           XSTART=DCMPLX(0.D0,0.D0)
           ncc=0
         endif
*        if(J.NE.1.or.a(M+1).eq.0.d0) then
* Find the root.
          CALL LAGUER(AD,J,X,ITS)
*           XSUM=XSUM+X
*        else
*          X=-a(M)/a(M+1)-XSUM
*        endif
        IF(ABS(DIMAG(X)).LE.2.d0*EPS*EPS*ABS(DBLE(X))) 
     &       X=DCMPLX(DBLE(X),0.d0)
        ROOTS(J)=X
        B=AD(J+1)
*         Forward deflation.
        DO 15 JJ=J,1,-1
          C=AD(JJ)
          AD(JJ)=B
          B=X*B+C
   15   CONTINUE
   20 CONTINUE
      IF(POLISH) THEN
*         Polish the roots using the undeflated coefficients.
        DO 30 J=1,M
          CALL LAGUER(A,M,ROOTS(J),ITS)
   30   CONTINUE
      ENDIF
      DO 40 J=2,M 
*         Sort roots by their real parts by straight insertion.
        X=ROOTS(J)
        DO 35 I=J-1,1,-1
          IF(DBLE(ROOTS(I)).LE.DBLE(X)) GO TO 37
          ROOTS(I+1)=ROOTS(I)
   35   CONTINUE
        I=0
   37   ROOTS(I+1)=X
   40 CONTINUE
      RETURN
      END


C
C
C========================================================================
      double precision function 
     & PANGLE(RE,turp,L,Q,BETA,ALPHA,A,delta,mu0)
      implicit none
      logical turp
      double precision a,aa,aeq,alpha,bb,beta,del,delta,
     &       enu,epsinu,f(4),gam,gamc,k(4),kt(4),L,mu0,
     &       omega,omegag,one,pi,q,ra2,re,ssqrt,vg
      one=1.d0
      pi=acos(-one)
      ra2=re*re+a*a
      del=ra2-2.d0*re
      aeq=ra2**2-a*a*del
      epsinu=aeq/re/re/ssqrt(del)
      enu=ssqrt(del/aeq)*re
      omegag=one/(re*ssqrt(re)+a)
      omega=2.d0*a*re/aeq
      vg=(omegag-omega)*epsinu
      gam=one/ssqrt(one-vg*vg)
C The following are the contravariant components of the photon
C momentum in the locally non-rotating frame in the equatorial plane:
      k(1)=(ra2/del*(ra2-L*a)-a*(a-L))/re**2
      k(2)=ssqrt((ra2-L*a)**2-del*((L-a)**2+Q*Q))/re**2
      if(turp) k(2)=-k(2)
      k(3)=-Q/re**2
C      if(theta0.gt.pi) k(3)=-k(3)
      k(4)=(-a+L+a/del*(ra2-L*a))/re**2
C The following are the same components, except in the tetrad
C notation (from Chandrasekhar, p. 347- note mistake in eq. 188):
      kt(1)=re*ssqrt(del/aeq)*k(1)
      kt(2)=re/ssqrt(del)*k(2)
      kt(3)=re*k(3)
      kt(4)=ssqrt(aeq)/re*(k(4)-omega*k(1))
      if(delta.ge.0.d0) then
C The following is the (unnormalized, but mod k) contravariant
C polarization vector in the LNRF, assuming it is in the plane
C of the disk in the emitting frame:
        f(1)=0.d0
        f(2)=ssqrt(del)/re*(vg*(kt(1)-kt(2)**2/kt(1))-kt(4))
        f(3)=-vg/re*kt(2)*kt(3)/kt(1)
        f(4)=re/ssqrt(aeq)*kt(2)*(one-vg*kt(4)/kt(1))
      else
        f(1)=0.d0
        f(2)=ssqrt(del)/re*kt(3)*kt(2)*(-one+vg*kt(4)/kt(1))
        f(3)=(kt(2)**2+(one+vg*vg)*kt(4)**2-2.d0*vg*kt(4)*kt(1)
     &       +vg*kt(3)**2*kt(4)/kt(1))/re
        f(4)=re/ssqrt(aeq)*kt(3)*
     &       (-(one+vg*vg)*kt(4)+vg*kt(1)+vg*kt(4)**2/kt(1))
      endif
C Evaluate Chandrasekhar's constants A & B in the LNRF (since
C we're in the equatorial plane, sin(theta)=1, cos(theta)=0):
C NOTE: Since mu0=0, k1=re*B, k2=re*A
      AA=k(1)*f(2)-k(2)*f(1)+a*(k(2)*f(4)-k(4)*f(2))
      BB=ra2*(k(4)*f(3)-k(3)*f(4))-a*(k(1)*f(3)-k(3)*f(1))
      gamc=-alpha-a*sqrt(one-mu0*mu0)
C Calculate the angle at infinity:
C     pangle=atan((beta*BB+gamc*AA)/(beta*AA-gamc*BB))
      pangle=atan2((beta*BB+gamc*AA),(beta*AA-gamc*BB))
C     if((beta*AA+gamc*BB).lt.0.d0) pangle=pangle+pi
      return
      end

C
C
C========================================================================
C
C
      SUBROUTINE LOCATE(XX,N,X,J,NS)
      INTEGER J,N,NS
      REAL*8 X,XX(NS)
      INTEGER JL,JM,JU
      JL=0
      JU=N+1
 10   IF(JU-JL.GT.1) THEN
        JM=INT((JU+JL)/2)
        IF((XX(N).GT.XX(1)).EQV.(X.GE.XX(JM))) THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
        GOTO 10
      ENDIF
      J=JL
      RETURN
      END

C
C
C========================================================================
C
C
 
C  Subroutine to calculate the abscissas and weights for the Gauss-Legendre-
C  Quadrature. The routine is based on the NUMERICAL RECIPES and uses an
C  algorithem of G.B. Rybicki.
C  Input: x1 ,x2: range of integration.
C          n: order of the orthogonal polynomials and the quadrature formula.
C  Output: x = x(n): array of the abscissas.
C          w = w(n): array of the weights.
 
 
      SUBROUTINE GAULEG(X1,X2,X,W,N,NS)
 
      INTEGER N,M,I,J
      REAL*8 X1,X2,X(NS),W(NS)
      REAL*8 PI,XM,XL,Z,P1,P2,P3,Z1,PP,EPS
      PARAMETER (PI = 3.14159265358979323846D0)
      PARAMETER (EPS=3.D-14)
 
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
         Z=COS(PI*(I-.25D0)/(N+.5D0))
 1       CONTINUE
         P1=1.D0
         P2=0.D0
         DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
 11      CONTINUE
         PP=N*(Z*P1-P2)/(Z*Z-1.D0)
         Z1=Z
         Z=Z1-P1/PP
         IF(ABS(Z-Z1).GT.EPS) GOTO 1
         X(I)=XM-XL*Z
         X(N+1-I)=XM+XL*Z
         W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
         W(N+1-I)=W(I)
 12   CONTINUE
      RETURN
      END
C
C
C=====================================================================
C The following subroutine is written to have safe square roots:
C
C
       function ssqrt(x)
       double precision x,ssqrt
       if(abs(x).lt.1.d-14) then
         ssqrt=0.d0
       elseif(x.lt.0.d0) then
         write(6,*) 'Negative value under square root sign = ', x
         ssqrt=sqrt(-x)
         pause
       else
         ssqrt=sqrt(x)
       endif
       return
       end
C
C
C========================================================================
C
C
      subroutine geor3(u0,uf,mu0,muf,a,l,l2,q2,tpm,tpr,su,sm)
      implicit none
       logical tpm,tpr
C Given a starting position (u0,mu0) where u0=1/r0, mu0=cos(theta0),
C and the final angle, muf, this subroutine calculates the final
C radius, uf, (radius of emission) for a geodesic with parameters (l,q)
C (lambda,q in Speith's notation), with tpm as an input (tpm is just true
C or false since there are at most 1 turning points since we're assuming
C the equatorial plane is opaque).  sr and sm are the sign of the u and
C mu integrals at the starting point.  The formulae are based
C on the appendix of Rauch and Blandford (1994), ApJ 421:46, tables 1-6,
C except for the cases where there are turning points.
C The original routine was modified so that the number of R turning points
C would not have to be taken into account.
      double precision a,aa,a2,al2,asech0,asechf,atanh,
     *       bb,b1,b2,b3,b4,c1,c2,c20,c2f2,c3,c4,c5,cn,cn2,
     *       dis,dn,ef,eff,ef0,ellf,f,fac,half,ib22,
     *      imu,jmu,km1,l,l2,m1,mneg,mpos,mu0,mtemp,muf,muneg,muplus,
     *       one,pi2,q2,ql2,qq,r,rb2,rf,s1,six,sixth,sm,su,
     *      sn,sn2,ssqrt,theta,third,tn2,TWO,u0,uf,uplus,x,ycheck,yy
      integer i,nreal
       complex*16 c(5),croot,root(4)
       tpr=.false.
      one=1.d0
       uplus=one/(one+sqrt(one-a*a))
      half=0.5d0
      third=one/3.d0
      TWO=2.d0
      six=6.d0
      sixth=one/six
      pi2=TWO*acos(-one)
      ql2=q2+l2
      if(a.eq.0.d0) then
        if(l2.eq.0.d0.and.q2.eq.0.d0) then
          muplus=one
          imu=su*sm*(muf-mu0)
        else
          mpos=q2/ql2
           muplus=ssqrt(mpos)
C          write(6,*) mu0
C If there is a turning point in mu, then we need to include  this.
C (We're assuming that muf<mu0, so that if we start out with increasing
C  mu, we have to have a turning point to decrease to muf):
           if(tpm) then
            imu=su*sm*(acos(abs(mu0/muplus))+acos(abs(muf/muplus)))/
     &          ssqrt(ql2)
C This is the second possibility:
           else
            imu=su*sm*(acos(abs(mu0/muplus))-acos(abs(muf/muplus)))/
     &          ssqrt(ql2)
           endif
        endif
C First few cases - if (l^2+q^2)>27, then we have 3 real roots 
C (Numerical Recipes, sec. 5.6).
        if(ql2.gt.27.d0) then
          theta=acos(54.d0/ql2-one)
C The roots are arranged so that b1>b2>b3
          b1=sixth*(one-TWO*cos((theta+pi2)*third))
          b2=sixth*(one-TWO*cos((theta-pi2)*third))
          b3=sixth*(one-TWO*cos(theta*third))
          m1=(b1-b2)/(b1-b3)
              c1=ssqrt((b1-b3)*ql2*half)
C Case 1- Unbound orbit:
           if(abs(u0-b2).lt.1.d-15) u0=b2
          if(u0.le.b2) then
            ef=ellf(ssqrt((b2-u0)/(u0-b3)/m1),m1)
C The maximum value of the integral is when u0=0, uf=b2:
            ef0=ellf(ssqrt(-b2/b3/m1),m1)
            jmu=imu*c1-ef
              if(abs(jmu).gt.ef0) then
                uf=-1.d0
                return
             endif
            call sncndn(jmu,m1,sn,cn,dn)
            cn2=cn**2
C Voila the final inverse radius (it doesn't depend on # r turning points):
              uf=(b1*(b2-b3)*cn2+b3*(b1-b2))/(b1-b2+(b2-b3)*cn2)
C Find if there was a radial turning point:
              if(jmu.gt.0.d0) tpr=.true.
C Case 2 - Bound orbit:
          elseif(u0.ge.b1) then
              if(u0.ne.b1) then
              ef=ellf(ssqrt((u0-b1)/(b1-b2)),m1)
              else
              ef=rf(0.d0,m1,one)
            endif
            jmu=imu*c1+ef
            call sncndn(jmu,m1,sn,cn,dn)
C Voila the final inverse radius:
              uf=(b1-b2*sn*sn)/cn/cn
              if(jmu.lt.0.d0) tpr=.true.
              if(uf.lt.b1.or.uf.ge.uplus) uf=-one
          else
            write(6,*) 'Unphysical a=0, q^+l^2>27'
            uf=-1.d0
            return
          endif
C Cases 3 & 4 - if (l^2+q^2)=27, the roots are 1/3,1/3,-1/6:
C Case 3 - Unbound orbit
C (This case orbits black hole, but never reaches r turning point).
        elseif(ql2.eq.27.d0.and.u0.lt.third) then
          f=ssqrt(TWO*u0+third)
          atanh=half*log(one+TWO*f/(one-f))
          x=ssqrt(27.d0)*half*imu+atanh
          uf=half*(((exp(TWO*x)-one)/(exp(TWO*x)+one))**2-third)
           if(uf.lt.0.d0.or.uf.gt.third) uf=-one 
C Case 4 - Bound orbit.
        elseif(ql2.eq.27.d0.and.u0.ge.third) then
           f=one/ssqrt(TWO*u0+third)
          atanh=half*log(one+TWO*f/(one-f))
          x=-ssqrt(27.d0)*half*imu+atanh
           uf=half*(((exp(TWO*x)-one)/(exp(TWO*x)+one))**-2-third) 
           if(uf.lt.third.or.uf.gt.uplus) uf=-one
C Handle case when ql2=0:
        elseif(ql2.eq.0.d0) then
          uf=u0+imu
           if(uf.lt.0.d0.or.uf.gt.uplus) uf=-one
C Case 5 - if (l^2+q2)<27, there are TWO complex roots.
        else
          aa=-(54.d0/ql2*(one+ssqrt(one-ql2/27.d0))-one)**third*sixth
          bb=one/36.d0/aa
          b3=aa+bb+sixth
          c2=ssqrt(b3*(3.d0*b3-one))
          c1=ssqrt(TWO*c2*ql2)
          m1=half+(six*b3-one)/8.d0/c2
          km1=rf(0.d0,m1,one)
          ef=ellf(TWO*ssqrt(c2*(u0-b3))/(c2-u0+b3),m1)
          jmu=imu*c1+ef
          if(abs(jmu).le.km1) then
C (c2f-uf+b3)>0:
            call sncndn(jmu,m1,sn,cn,dn)
            tn2=(sn/cn)**2
              uf=c2+b3+TWO*c2*(one+sqrt(one+tn2))
          else
C (c2f-uf+b3)<0:
            jmu=TWO*km1-abs(jmu)
            call sncndn(jmu,m1,sn,cn,dn)
            tn2=(sn/cn)**2
             uf=c2+b3+TWO*c2/tn2*(one+sqrt(one+tn2))
          endif
           if(uf.lt.0.d0.or.uf.gt.uplus) uf=-one
        endif
C Next, consider a=/=0, but q^2=0:
      elseif(q2.eq.0.d0) then
         al2=a*a-l2
         s1=sign(one,mu0)
         if(mu0.ne.0.d0) then
           x=ssqrt(one-(l/a)**2)/abs(mu0)
           asech0=log(x+ssqrt(x*x-one))
         else
           asech0=0.d0
         endif
         if(muf.ne.0.d0) then
           x=ssqrt(one-(l/a)**2)/abs(muf)
           asechf=log(x+ssqrt(x*x-one))
         else
           asechf=0.d0
         endif
         if(tpm) then
           imu=su*sm*s1*ssqrt(a*a-l2)*(asech0+asechf)
         else
           imu=su*sm*s1/ssqrt(a*a-l2)*(asech0-asechf)
         endif
        r=0.25d0/(a-l)**2+((a+l)/(a-l))**3/216.d0
        qq=((a+l)/(a-l))**2/36.d0
        dis=r*r-qq**3
C Case 6 is special - q^2=0, l=a:
        if(a.eq.l) then
C In this case, M=-a^2*mu^4, so M+-=0 => orbit lies in equatorial
C plane.  So, it will be absorbed by the accretion disk:
           if(mu0.eq.0.d0) then
            uf=u0-imu
           else
             uf=-1.d0
           endif
           return
C Case 7 is also special - q^2=0, l=-a:
C In this case, mu=0 for all r, so the orbit lies in the equatorial
C plane.  Then, it will be absorbed by the disk, so we don't have
C to worry about its propagation: 
        elseif(a.eq.-l) then
           if(u0.eq.0.d0.and.muf.eq.0.d0) then
             uf=-1.d0
             return
           endif
          m1=(TWO-sqrt(3.d0))*0.25d0
           km1=rf(0.d0,m1,one)
          c2=TWO*a**(2.d0/3.d0)
          c1=3.d0**0.25d0*c2
          if((sqrt(3.d0)-one-c2*u0).ne.0.d0) then
            ef0=ellf(3.d0**0.25d0*TWO*ssqrt(one+c2*u0)/
     &        (ssqrt(3.d0)-one-c2*u0),m1)
          else
            ef0=km1
          endif
           jmu=imu*c1+ef0
           if(abs(jmu).ge.km1) then
              jmu=TWO*km1-abs(jmu)
            call sncndn(jmu,m1,sn,cn,dn)
             tn2=(sn/cn)**2
             uf=(sqrt(3.d0)*(one+TWO/tn2*(one+sqrt(one+tn2)))-one)/c2
           else
            call sncndn(jmu,m1,sn,cn,dn)
             tn2=(sn/cn)**2
              uf=(sqrt(3.d0)*(one+TWO*(one+sqrt(one+tn2)))-one)/c2
           endif
C Next cases - three real roots (q^2=0, l=/=+-a):
C (Check to see if r^2-q^3<0)
        elseif(dis.le.0.d0) then
C          theta=acos(al2/abs(al2)+54.d0*abs(a-l)/abs((a+l)**3))
          theta=acos(r/qq**1.5d0)
C   The roots are arranged so that b1>b2>b3:
          fac=-sixth/(a-l)**2
          b1=fac*(al2+TWO*abs(al2)*cos((theta+pi2)*third))
          b2=fac*(al2+TWO*abs(al2)*cos((theta-pi2)*third))
          b3=fac*(al2+TWO*abs(al2)*cos(theta*third))
C   Next cases - 3 non-equal real roots:
          if(theta.ne.0.d0) then
            m1=(b1-b2)/(b1-b3)
              km1=rf(0.d0,m1,one)
              c1=ssqrt((b1-b3)*half)*abs(a-l)
C Case 8 - u<b2 (unbound orbit):
            if(u0.le.b2) then
              if(b2.ne.u0) then
                ef=ellf(ssqrt((b2-u0)/(u0-b3)/m1),m1)
              else
                ef=km1
              endif
              jmu=imu*c1-ef
              call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
C Voila the final inverse radius (it doesn't depend on # r turning points):
              uf=(b2+b3*tn2*m1)/(tn2*m1+one)
                if(jmu.gt.0.d0) tpr=.true.
                if(uf.lt.0.d0.or.uf.gt.b2) uf=-one
C Case 9 - u>b1 (bound orbit):
            elseif(u0.ge.b1) then
                if(u0.ne.b1) then
                ef=ellf(ssqrt((u0-b1)/(b1-b2)),m1)
                else
                ef=km1
              endif
              jmu=ef+imu*c1
              call sncndn(jmu,m1,sn,cn,dn)
C Voila the final inverse radius:
                uf=(b1-b2*sn*sn)/cn/cn
                if(jmu.lt.0.d0) tpr=.true.
                if(uf.lt.b1.or.uf.ge.uplus) uf=-one
            else
      write(6,*) 'wrong range - q^2=0, a=/=|l|, 3 distinct real roots'
              uf=-1.d0
              return
            endif
C   Next cases - 2 equal (b1=b2), 3 real roots:
          else
C  Case 10 - u<b2 (unbound orbit):
            if(u0.lt.b2) then
              c1=ssqrt((b2-b3)*half)*abs(a-l)
              f=ssqrt((u0-b3)/(b2-b3))
              atanh=half*log(one+TWO*f/(one-f))
              x=imu*c1+atanh
              uf=b3+(b2-b3)*((exp(TWO*x)-one)/(exp(TWO*x)+one))**2
                if(uf.lt.0.d0.or.uf.ge.b2) uf=-one
C  Case 11 - u>b1 (bound orbit):
            elseif(u0.gt.b1) then
              c1=-ssqrt((b2-b3)*half)*abs(a-l)
              f=ssqrt((b2-b3)/(u0-b3))
              atanh=half*log(one+TWO*f/(one-f))
              x=imu*c1+atanh          
              uf=b3+(b2-b3)*((exp(TWO*x)+one)/(exp(TWO*x)-one))**2
                if(uf.le.b1.or.uf.ge.uplus) uf=-one
            else
      write(6,*) 'wrong range-q^2=0, a=/=|l|, 2 equal, 3 real roots'
              uf=-one
              return
            endif
          endif
C  Case 12 - q^2=0, l=/=+-a, 2 complex roots:
C  These orbits terminate in singularity without passing
C  through the equatorial plane:
        else
          s1=sign(1.d0,r)
          aa=-s1*(abs(r)+ssqrt(dis))**third
          if(aa.ne.0.d0) then
            bb=qq/aa
          else
            bb=0.d0
           endif
          c3=(a+l)/(a-l)
           b3=aa+bb-c3*sixth
          c2=ssqrt(b3*(3.d0*b3+c3))
          c1=ssqrt(TWO*c2)*abs(a-l)
          m1=half+(six*b3+c3)/8.d0/c2
           km1=rf(0.d0,m1,one)
          if((c2-u0+b3).ne.0.d0) then
            ef0=ellf(TWO*ssqrt(c2*(u0-b3))/(c2-u0+b3),m1)
          else
            ef0=km1
          endif
           jmu=imu*c1+ef0
           if(abs(jmu).gt.km1) then
              jmu=TWO*km1-abs(jmu)
              call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
             uf=c2*(one+TWO/tn2*(one+sqrt(one+tn2)))+b3
           else
              call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
              uf=c2*(one+TWO/tn2*(one-sqrt(one+tn2)))+b3
           endif
           if(uf.lt.0.d0.or.uf.gt.uplus) uf=-one
        endif
C Now consider cases with q^2=/=0:
      else
        c(1)=dcmplx(one,0.d0)
        c(2)=dcmplx(0.d0,0.d0)
        c(3)=dcmplx(a*a-q2-l2,0.d0)
        c(4)=dcmplx(TWO*((a-l)*(a-l)+q2),0.d0)
        c(5)=dcmplx(-a*a*q2,0.d0)
        call zroots(c,4,root,.true.)
        nreal=0
        do i=1,4
          if(dimag(root(i)).eq.0.d0) nreal=nreal+1
        enddo
        if(q2.lt.0.d0) then
           a2=a*a
           yy=-0.5d0*(a2-ql2+sign(one,a2-ql2)*ssqrt((a2-ql2)**2+
     &         4.d0*q2*a2))
           if((a2-ql2).lt.0.d0) then
              mneg=-yy/a2
              mpos=q2/yy
           else
              mneg=q2/yy
              mpos=-yy/a2
           endif
           if(mneg.gt.mpos) then
              mtemp=mneg
              mneg=mpos
              mpos=mtemp
           endif
          muplus=ssqrt(mpos)
           muneg=ssqrt(mneg)
          s1=sign(one,mu0)
           if(muf.gt.s1*muplus.or.muf.lt.s1*muneg) then
             uf=-one
             return
                endif
          m1=mneg/mpos
           if(muf.eq.0.d0) then
              eff=rf(0.d0,m1,one)
           else
              eff=ellf(ssqrt((mpos-muf*muf)/(muf*muf-mneg)),m1)
           endif
           if(mu0.eq.0.d0) then
            ef0=rf(0.d0,m1,one)
           elseif(mu0.eq.1.d0) then
              ef0=0.d0
           else
              ef0=ellf(ssqrt((mpos-mu0*mu0)/(mu0*mu0-mneg)),m1)
           endif
            if(tpm) then
             imu=sm*su*s1/abs(a)/muplus*(ef0+eff)
           else
            imu=sm*su*s1/abs(a)/muplus*(ef0-eff)
               endif
C Case 13 - a =/= 0, q^2<0, 2 real roots, 2 complex roots:
          if(nreal.eq.2) then
            croot=0.d0
            b3=0.d0
            b4=0.d0
            do i=1,4
              if(dimag(root(i)).ne.0.d0.and.croot.eq.0.d0) croot=root(i)
              if(dimag(root(i)).eq.0.d0.and.
     &          dble(root(i)).lt.min(b3,b4)) then
                b3=b4
                b4=dble(root(i))
              endif
              if(dimag(root(i)).eq.0.d0.and.dble(root(i)).lt.b3.and.
     &            dble(root(i)).gt.b4) b3=dble(root(i))
            enddo
              rb2=dble(croot)
              ib22=dimag(croot)
            c4=ssqrt((rb2-b3)**2+ib22)
            c5=ssqrt((rb2-b4)**2+ib22)
            c20=ssqrt(c5*(u0-b3)/c4/(u0-b4))
            c1=ssqrt(-a*a*q2*c4*c5)
              m1=rb2*rb2-rb2*(b3+b4)+b3*b4
              m1=-(m1+ib22)/two/c4/c5+0.5d0
c            m1=((b3-b4)**2-(c4-c5)**2)*0.25d0/c4/c5
              km1=rf(0.d0,m1,one)
            if(one.ne.c20*c20) then
              ef0=ellf(TWO*c20/(one-c20*c20),m1)
            else
              ef0=km1
            endif
              jmu=imu*c1+ef0
            if(abs(jmu).le.km1) then
C c2f^2<1:
              call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
                yy=one+TWO/tn2
               c2f2=one/(yy+ssqrt(yy*yy-one))
C c2f^2>1:
              else
                jmu=TWO*km1-abs(jmu)
               call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
              yy=one+TWO/tn2
                c2f2=yy+ssqrt(yy*yy-one)
              endif
              uf=(c4*c2f2*b4-b3*c5)/(c2f2*c4-c5)
              if(uf.lt.0.d0.or.uf.gt.uplus) uf=-one
C Case 14 - four complex roots
          elseif(nreal.eq.0) then
            c5=ssqrt((dble(root(1))-dble(root(3)))**2+
     &           (dimag(root(1))-dimag(root(3)))**2)
            c4=ssqrt((dble(root(1))-dble(root(3)))**2+
     &           (dimag(root(1))+dimag(root(3)))**2)
            c2=ssqrt((4.d0*dimag(root(1))**2-(c4-c5)**2)/
     &            ((c4+c5)**2-4.d0*dimag(root(1))**2))
            c3=dble(root(1))+c2*dimag(root(1))
            c1=half*(c4+c5)*ssqrt(-a*a*q2)
C The following formula prevents some roundoff error:
              m1=16.d0*(dimag(root(1))*dimag(root(3)))**2/(c4+c5)**4
            if((dimag(root(1))*(one+c2*c2)+c2*(u0-c3)).ne.0.d0) then
              ef0=ellf((u0-c3)/(dimag(root(1))*(one+c2*c2)+
     &            c2*(u0-c3)),m1)
            else
              ef0=rf(0.d0,m1,one)
            endif
              jmu=imu*c1+ef0
            call sncndn(jmu,m1,sn,cn,dn)
              uf=c3+(sn/cn*dimag(root(1))*(one+c2*c2))/(one-c2*sn/cn)
              if(uf.lt.0.d0.or.uf.gt.uplus) uf=-one
          else
            write(6,*) 'Odd number of real roots for a=/=0, q^2<0 case'
          endif
C Final cases - q^2>0:
C (These do reach equatorial plane!)
        else
           a2=a*a
           yy=-0.5d0*(a2-ql2+sign(one,a2-ql2)*ssqrt((a2-ql2)**2+
     &         4.d0*q2*a2))
           if((a2-ql2).lt.0.d0) then
              mneg=-yy/a2
              mpos=q2/yy
           else
              mneg=q2/yy
              mpos=-yy/a2
           endif
           if(mneg.gt.mpos) then
              mtemp=mneg
              mneg=mpos
              mpos=mtemp
           endif
          muplus=ssqrt(mpos)
          s1=sign(one,mu0)
          m1=-mneg/(mpos-mneg)
           if(muf.eq.0.d0) then
              eff=rf(0.d0,m1,one)
           else
              eff=ellf(ssqrt(mpos/muf**2-one),m1)
           endif
           if(mu0.eq.0.d0) then
            ef0=rf(0.d0,m1,one)
           elseif(mu0.eq.1.d0) then
              ef0=0.d0
           else
              ef0=ellf(ssqrt((muplus/mu0)**2-one),m1)
           endif
            if(tpm) then
             imu=sm*su*s1/abs(a)/ssqrt(mpos-mneg)*(ef0+eff)
           else
            imu=sm*su*s1/abs(a)/ssqrt(mpos-mneg)*(ef0-eff)
               endif
          if(nreal.eq.4) then
            b1=dble(root(4))
            b2=dble(root(3))
            b3=dble(root(2))
            b4=dble(root(1))
            if(b2.ne.b3)then
C Case 15 - 4 distinct, real roots, u<b2:
                c1=half*ssqrt(a*a*q2*(b1-b3)*(b2-b4))
                if(abs(u0-b3).lt.1.d-12) u0=b3
                if(abs(u0-b2).lt.1.d-12) u0=b2
              if(u0.le.b3) then
                m1=(b1-b4)*(b2-b3)/(b1-b3)/(b2-b4)
                  km1=rf(0.d0,m1,one)
                if(u0.ne.b3) then
C The following is the u integral from u0 to b3:
                  ef=ellf(ssqrt((b3-u0)*(b2-b4)/(b2-b3)/(u0-b4)),m1)
                else
                  ef=km1
                endif
                jmu=imu*c1-ef
                call sncndn(jmu,m1,sn,cn,dn)
                  sn2=sn*sn
C Voila the final inverse radius:
                uf=((b2-b4)*b3-(b3-b4)*b2*sn2)/(b2-b4-(b3-b4)*sn2)
                  if(jmu.gt.0.d0) tpr=.true.
C Case 16 - 4 distinct, real roots, u>b2:
              elseif(u0.ge.b2) then
                m1=(b1-b2)*(b3-b4)/(b1-b3)/(b2-b4)
                  km1=rf(0.d0,m1,one)
                ef=ellf(ssqrt((b1-b3)*(u0-b2)/(b2-b3)/(b1-u0)),m1)
                  jmu=imu*c1+ef
                  call sncndn(jmu,m1,sn,cn,dn)
                  sn2=sn**2
                  uf=(b2*(b1-b3)-b3*(b1-b2)*sn2)/(b1-b3-(b1-b2)*sn2)
                  if(jmu.lt.0.d0) tpr=.true.
                  if (uf.lt.b2.or.uf.ge.b1) uf=-one
              else
      write(6,*) 'Range is unphysical for a=/=0, q^2>0,4 real, distinct'
                uf=-1.d0
                return
              endif
C Cases 17 & 18 - 4 real, TWO equal roots, u<b3:
            else
              c3=(b1-b3)*(b3-b4)
              c2=(TWO*b3-b1-b4)/(b1-b4)
                if(u0.lt.b3) c1=ssqrt(a*a*q2*c3)
                if(u0.lt.b1.and.u0.gt.b2) c1=-ssqrt(a*a*q2*c3)
              jmu=exp(imu*c1)*abs(c2+TWO*(c3+ssqrt(c3*
     &             (b1-u0)*(u0-b4)))/(b1-b4)/(b3-u0))
                uf=b3*(b1*(one-jmu)**2-b4*(one+jmu)**2+
     &             4.d0*jmu*b1*b4/b3)/
     &             (b1*(one+jmu)**2-b4*(one-jmu)**2-4.d0*jmu*b3)
                if(uf.gt.b3) then
      write(6,*) 'q2>0, a=/=0, 4 real, TWO = roots, uf out of range'
                endif
                ycheck=c2+TWO*(c3+ssqrt(c3*(b1-uf)*(uf-b4)))/
     &                 (b1-b4)/(b3-uf)
                if(ycheck.lt.0.d0) then
      write(6,*) 'q2>0, a=/=0, 4 real, TWO = roots, wrong sign'
                  pause
                endif
            endif
C Case 19 - 2 real, 2 complex roots:
          elseif(nreal.eq.2) then
            do i=1,4
              if(dimag(root(i)).eq.0.d0) then
                if(dble(root(i)).ge.0.d0) b1=dble(root(i))
                if(dble(root(i)).lt.0.d0) b4=dble(root(i))
              else
                rb2=dble(root(i))
                ib22=dimag(root(i))**2
              endif
            enddo
            c4=ssqrt((rb2-b1)**2+ib22)
            c5=ssqrt((rb2-b4)**2+ib22)
              m1=rb2*rb2-rb2*(b1+b4)+b1*b4
              m1=(m1+ib22)/two/c4/c5+0.5d0
              c1=ssqrt(a*a*q2*c4*c5)
            km1=rf(0.d0,m1,one)
            if(abs(b1-u0).gt.1.d-14) then
              c20=ssqrt(c4*(u0-b4)/c5/(b1-u0))
            else
              c20=0.d0
            endif
            if(c20.ne.one) then
              ef=ellf(TWO*c20/(one-c20*c20),m1)
            else
              ef=km1
            endif
            jmu=imu*c1+ef-TWO*km1
            if(jmu.gt.0.d0.and.b1.ge.uplus) then
              uf=-1.d0
              return
            elseif(jmu.gt.0.d0) then
              tpr=.true.
            endif
            jmu=TWO*km1-abs(jmu)
            if(abs(jmu).le.km1) then
C c2f^2<1:
              call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
                yy=one+TWO/tn2
               c2f2=one/(yy+ssqrt(yy*yy-one))
C c2f^2>1:
              else
                jmu=TWO*km1-abs(jmu)
               call sncndn(jmu,m1,sn,cn,dn)
              tn2=(sn/cn)**2
              yy=one+TWO/tn2
                c2f2=yy+ssqrt(yy*yy-one)
              endif
              uf=(b1*c2f2*c5+b4*c4)/(c4+c2f2*c5)
          else
            write(6,*) 'Odd number of real roots for a=/=0, q^2>0 case'
          endif
        endif
      endif
      return
      end



