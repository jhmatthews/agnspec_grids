      program dispar3
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (VELC=2.997925E10)
      PARAMETER (GRCON = 6.668D-8)
      PARAMETER (SIG4P = 4.5114062D-6)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0,
     *           HALF=ONE/TWO,THIRD=ONE/THREE)
      parameter (mre=100)
      DIMENSION R(MRE)
C
C ----------------------
C Basic input parameters
C ----------------------
C
C     MSTAR      - M(star), either in M(Sun), or in grams;
C     MDOT       - M(dot), either in M(Sun)/year; or in g/s
C     AA         - angular momentum of black hole
C     ALPHA      - Shakura-Sunyaev viscosity parameter
C     NRE        - number of radial annuli
C
      READ(5,*) XMSTAR,XMDOT,AA,ALPHA
      READ(5,*) ROUT,TMIN,DELTAR,NRE
C
      IF(ALPHA.EQ.0.) ALPHA=0.01
C
      xmstar=abs(xmstar)
      RSTAR=XMSTAR*1.989D33*GRCON/VELC/VELC
      IF(XMSTAR.GT.1.D16) XMSTAR=XMSTAR/1.989D33
      IF(XMDOT.GT.1.D3) XMDOT=XMDOT/6.3029D25
      IF(RSTAR.GT.1.D3) RSTAR=RSTAR/6.9598D10
C
C     Minimum radius for last stable circular orbit per unit mass, X0
C
      AA2=AA*AA
      Z1 = ONE + (ONE-AA2)**THIRD * ((ONE+AA)**THIRD + (ONE-AA)**THIRD)
      Z2 = SQRT(THREE*AA2 + Z1*Z1)
      RMS = THREE + Z2 - SQRT((THREE-Z1)*(THREE+Z1+TWO*Z2))
C
C     set up radii to represent a disk  
C
c     first, the outer cutoff
C
      IF(TMIN.GT.0.) THEN
         RE=100.
         R0=RSTAR*ABS(RE)
         QGRAV=5.9D0*GRCON*XMSTAR/R0**3
	 CALL GRCOR(AA,RE,QCOR,TCOR,ARH,BRH,CRH,DRH)
         TEFF0=(1.79049311D-1*QGRAV*3.34379D24*XMDOT/SIG4P)**0.25
	 TEFF=TEFF0*TCOR
         ROUT=100.*(TMIN/TEFF)**(-1.333)
      END IF
C
C     if DELTAR is set to a non-zero value, logarithmically
C     equidistant R between 1.05 *RMS and ROUT
C
      IF(DELTAR.GT.0.) THEN
         R1=LOG10(RMS*1.05)
         R2=LOG10(ROUT)
         NRE=(R2-R1)/DELTAR+1
         DO I=1,NRE
            R(I)=R1+(I-1)*DELTAR
            R(I)=10.**R(I)
         END DO
       ELSE
C
C     otherwise, a Gaussian integration between RMS and ROUT
C
         R1 = ZERO
         R2 = LOG(ROUT/RMS)
         CALL GAULEG(R1,R2,R,WR,NRE,MRE)
         DO I = 1, NRE
            R(I) = RMS * EXP( R(I) )
         END DO  
      END IF
C
C     basic parameters for the individual annuli
C
      DO I = 1, NRE
         R0=RSTAR*R(I)
         RE=R(I)
         QGRAV=5.9D0*GRCON*XMSTAR/R0**3
	 CALL GRCOR(AA,RE,QCOR,TCOR,ARH,BRH,CRH,DRH)
         TEFF0=(1.79049311D-1*QGRAV*3.34379D24*XMDOT/SIG4P)**0.25
	 TEFF=TEFF0*TCOR
	 QGRAV=QGRAV*QCOR
C
C Compute the Keplerian rotation frequency (omega=c/r_g/x^1.5);
C and Relativistic factors, using Krolik (1998) notation:
C
         OMEGA=VELC/RSTAR/6.9698D10/RE**1.5D0
         RELT=DRH/ARH
         RELR=DRH/BRH
         RL2=RE*(1.D0-2.D0*AA/RE**1.5D0+AA*AA/RE/RE)**2/BRH
         EINF=(1.D0-2.D0/RE+AA/RE**1.5D0)/SQRT(BRH)
         RELZ=(RL2-AA*AA*(EINF-1.D0))/RE
C
C Compute the surface mass density (assuming pure electron scattering, 
C pure hydrogen composition, and that T_rphi = \alpha P_tot):
C
         XMD=XMDOT*6.3029D25
         DMTOT=0.5*SIGMAR(ALPHA,XMD,TEFF,OMEGA,RELR,RELT,RELZ)
C
         TEFFL=LOG10(TEFF)
         QGRAL=LOG10(QGRAV)
         DMTOL=LOG10(DMTOT)
         WRITE(6,601) R(I),TEFFL,DMTOL,QGRAL,TEFF,DMTOT,QGRAV
      END DO
  601 FORMAT(4f10.3,f10.1,1p2e13.4)
      END
C
C
C
C     ****************************************************************
C
	SUBROUTINE GRCOR(AA,RR,QCOR,TCOR,ARH,BRH,CRH,DRH)
C	=================================================
C
C	Procedure for computing general-relativistic correction
C	factors to gravitational factor (QGRAV) and effective 
C	temperature (TEFF)
C       Also calculates all frour quantities in the Riffer-Herlod (RH)
C       notation - ARH, BRH, CRH, DRH
C
C       Input:
C             AA   - angular momentum (0.98 maximum)
C             RR   - R/R_g = r/(GM/c^2)
C       Outout:
C             QCOR - g-correction  = C/B   in RH notation
C             TCOR - T-correction  = (D/B)^(1/4)   in RH notation
C             ARH  - A  in RH notation
C             BRH  - B  in RH notation
C             CRH  - C  in RH notation
C             DRH  - D  in RH notation
C
        IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (THIRD=1.D0/3.D0, PI3=1.0471976)
C
C  ----------------
C  Imput parameters
C  ----------------
C
C	AA	- specific angular momentum/mass
C		  of the Kerr black hole
C	RR	- distance/mass of the Kerr black hole
C
C  ---------------------------------
C  Set correcion factors A through G  (see Novikov & Thorne,'73, eq.5.4.1a-g)
C  ---------------------------------
C
        rror=rr
        rr=abs(rr)
	AA2=AA*AA
	RR1=1/RR
	RR12=SQRT(RR1)
	RR2=RR1*RR1
	A2R2=AA2*RR2
	A4R4=A2R2*A2R2
	A2R3=AA2*RR2*RR1
	AR32=SQRT(A2R3)
C
	A = 1 +   A2R2 + 2*A2R3
	B = 1 +   AR32
	C = 1 - 3*RR1  + 2*AR32
	D = 1 - 2*RR1  +   A2R2
	E = 1 + 4*A2R2 - 4*A2R3 + 3*A4R4
C
C  -------------------------------
C  Set correction factor for QGRAV  (see Novikov & Thorne,'73, eq.5.7.2)
C  -------------------------------
C
	if(rror.lt.0) QCOR = B*B*D*E/(A*A*C)
c
c  correction - after Riffert and Harold 
c
        if(rror.gt.0) QCOR = (1. - 4.*AR32 + 3.*A2R2)/C
C
C  -----------------------
C  Set correction factor Q  (see Page & Thorne,'73, eq.35)
C  -----------------------
C
C	Minimum radius for last stable circular orbit per unit mass, X0
C
	Z1 = 1 + (1-AA2)**THIRD * ((1+AA)**THIRD + (1-AA)**THIRD)
	Z2 = SQRT(3*AA2 + Z1*Z1)
	R0 = 3 + Z2 - SQRT((3-Z1)*(3+Z1+2*Z2))
	X0 = SQRT(R0)
c
        if(rr.lt.r0) then
           arh=1.
           brh=1
           crh=1
           drh=1
           tcor=1.
           qcor=1.
           return
        end if
C
C	Roots of x^3 - 3x + 2a = 0
C
	CA3 = THIRD * ACOS(AA)
	X1 = 2*COS(CA3-PI3)
	X2 = 2*COS(CA3+PI3)	
	X3 = -2*COS(CA3)
C
C	FB = '[]' term in eq. (35) of Page&Thorne '73
C
	X = SQRT(RR)
	C1 = 3*(X1-AA)*(X1-AA)/(X1*(X1-X2)*(X1-X3))
	C2 = 3*(X2-AA)*(X2-AA)/(X2*(X2-X1)*(X2-X3))
	C3 = 3*(X3-AA)*(X3-AA)/(X3*(X3-X1)*(X3-X2))
	AL0 = 1.5*AA*log(X/X0)
	AL1 = log((X-X1)/(X0-X1))
	AL2 = log((X-X2)/(X0-X2))
	AL3 = log((X-X3)/(X0-X3))
	FB = (X-X0 - AL0 - C1*AL1 - C2*AL2 - C3*AL3)
	Q = FB*(1+AR32)*RR12/SQRT(1-3*RR1+2*AR32)

C  ------------------------------
C  Set correction factor for TEFF  (see Novikov & Thorne,'73, eq.5.5.14b)
C  ------------------------------
C
	TCOR = (Q/B/SQRT(C))**0.25
C
C  ------------------------------
C  RH quantities
C  ------------------------------
C 
        ARH = D
        BRH = C
        CRH = 1. - 4.*AR32 + 3.*A2R2
        DRH = Q/B*SQRT(C)
C 
	RETURN
	END
C
C
C     ****************************************************************
C
C
      FUNCTION SIGMAR(ALPHA,XMDOT,TEF,OMEGA,RELR,RELT,RELZ)
C     =====================================================
c
C--------------------------------------------------------------------
C The following function takes as inpute various disk parameters computed
C at a certain radius in cgs units, and outputs the surface mass density, 
C assuming that the opacity is electron-scattering dominated (or that kappa 
C is independent of the density).  The equations were derived assuming a 1-zone
C model (i.e. rho, T_g, mu are constant with height), assuming t_r,phi =
C -alpha P, and that the dissipation is constant per unit optical depth. 
C See chapter 7 of Krolik for information on the notation used here 
C (especially relativistic factors).
C     
C--------------------------------------------------------------------
C Uses Numerical Recipes subroutine LAGUER
C--------------------------------------------------------------------
C ALPHA - viscosity parameter
C XMDOT - accretion rate (in g/s)
C TEF   - temperature in Kelvins
C OMEGA - Keplerian frequency (in Hz)
C RELR,RELT, RELZ  - relativistic factors 
C KAPPA - opacity (cm^2/g)
C MU    - mean atomic mass (g) = rho/N
C--------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (HALF=0.5D0)
      REAL KAPPA,MU
      COMPLEX*16 COEFF(11),XGUESS
C
C We should check that the physical constants used here agree with those in
C disk195g:
C
      parameter (ZERO=0.D0,ONE=1.D0,TRES=3.D0,
     *           FOUR=4.D0,THIRD=ONE/TRES,FOURTH=0.25D0,EPS=1.D-5)
      parameter (C=2.9979D10,SIGMAB=5.6703D-5)
      parameter (BK=1.3807D-16)
      PI=ACOS(-ONE)
C
C We'll assume fully ionized, pure hydrogen:
C
      KAPPA=0.40D0
      MU=0.5D0*1.6726D-24
      FAC1=RELZ*(HALF*TRES*C*OMEGA/ALPHA/KAPPA/SIGMAB/TEF**4)**2
      FAC2=(HALF*KAPPA)**FOURTH*BK*TEF/MU
      FAC3=XMDOT*OMEGA*RELT/PI
C
C Coefficients of the equation for x^4=Sigma:
C
      COEFF(1)=DCMPLX(FAC1*(HALF*FAC3)**2,ZERO)
      COEFF(2)=ZERO
      COEFF(3)=ZERO
      COEFF(4)=ZERO
      COEFF(5)=DCMPLX(-(TRES*FAC3)/(8.D0*ALPHA),ZERO)
      COEFF(6)=DCMPLX(-FAC1*FAC3*ALPHA*FAC2,ZERO)
      COEFF(7)=ZERO
      COEFF(8)=ZERO
      COEFF(9)=ZERO
      COEFF(10)=DCMPLX(FOURTH*FAC2,ZERO)
      COEFF(11)=DCMPLX(FAC1*(ALPHA*FAC2)**2,ZERO)
C
C At small radii, we'll approximate P_rad >> P_gas
C First, compute sigma assuming radiation pressure dominates:
C
      SIGRAD=FOUR*OMEGA*C*C*RELT*RELZ/ALPHA/KAPPA**2/SIGMAB/TEF**4/RELR
C
C Next, compute sigma assuming gas pressure dominates:
C
      SIGGAS=((MU*XMDOT*OMEGA*RELT/PI/ALPHA/BK/TEF)**4/8./KAPPA)**0.2D0
C
C Use a starting guess which has the correct value for P_gas >> P_rad
C or P_rad >> P_gas.
C
      XGUESS=DCMPLX(ONE/(ONE/SIGRAD+ONE/SIGGAS)**FOURTH,ZERO)
C
C Look for root of the 10th order equation for x:
C
      CALL LAGUER(COEFF,10,XGUESS,ITS)
C
C Make sure that we haven't landed a wrong root:
C
      IF(ABS(DIMAG(XGUESS)).LT.EPS.AND.DBLE(XGUESS).GT.ZERO) THEN
        SIGMAR=DBLE(XGUESS)**4
      ELSE
        SIGMAR=ONE/(ONE/SIGRAD+ONE/SIGGAS)
        WRITE(6,*) 'Surface density approximated'
      ENDIF
c      WRITE(6,2000) TEF,SIGRAD,SIGGAS,SIGMAR
      RETURN
 2000 FORMAT(20(2x,1pe12.5))
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE laguer(a,m,x,its)
C     ============================
C
C     Routine from Numerical Recipees
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      REAL frac(MR)
      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
C
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.d0*f/b
          sq=sqrt(dble(m-1)*(dble(m)*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.D0) then
            dx=dble(m)/gp
          else
            dx=exp(dcmplx(log(1.d0+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END
C
C
C     ****************************************************************
C
      SUBROUTINE GAULEG(X1,X2,X,W,N,NS)
C     =================================
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(NS),W(NS)
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
C



