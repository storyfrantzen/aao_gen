c
c*******************Valery Model****************************************
c
      subroutine DVMPW(COST,W,Q2, Phi_g,E,heli,MESONMASS,
     *                   sigma0,S_T,S_L,S_TT,S_LT,S_LTP)
C
c  dsigma/dQ2 dphi dCosTheta* dW for ep-->ep pi0
C
C  exc pi0 x-section
c
c input:
c COST = Cos(Theta*), Theta* is the angle (pi0-gamma* or  p-p')in the CM system
c W is W
c xb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c
c del2=t (negative GeV^2)           NEGATIVE !!!!
C
C MESONMASS is the mass of the pi0 or eta.
c

      IMPLICIT NONE
      REAL    COST,W,del2,xb,Q2, Phi_g,E,mesonmass
      INTEGER heli
      REAL    Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     sigma0, S_T, S_L, S_LT, S_TT,S_LTP
      REAL     SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT,SIGMA_LTP
      REAL     EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,SIGMA_T,SIGMA_L, SIGMA_LT, SIGMA_TT
      EXTERNAL SIGMA_LTP
      REAL     JACG,JACR
      LOGICAL  CHECK_KINE
      REAL S_TOT
      DOUBLE PRECISION eps
C      INCLUDE 'pi0eta.par'
c
      CALL XSINIT(MESONMASS)
c     
      sigma0 = 0.0
      S_T = 0.0
      S_L = 0.0
      S_LT = 0.0
      S_TT = 0.0
      S_LTP = 0.0
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)) RETURN 
      EPS=EPSILON(XB,q2,e)

      

      S_T  =  0.001/(2*pi)*SIGMA_T     (COST,W,Q2,E)
      S_L  =  0.001/(2*pi)*SIGMA_L     (COST,W,Q2,E)
      S_LT =  0.001/(2*pi)*2*SIGMA_LT  (COST,W,Q2,E)
      S_TT =  0.001/(2*pi)*SIGMA_TT    (COST,W,Q2,E)
      S_LTP = 0.001/(2*pi)*2*SIGMA_LTP (COST,W,Q2,E)

      sigma0 = S_T   + EPS*S_L +         
     *   EPS                  *S_TT  * COS(2*PHI_G)   + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
       IF(sigma0 .LT. 0.0) sigma0=0.0     
c
c       S_TOT =FLUXW(xb,Q2,E)/(2.*PI)*(
c     *                       S_T   + EPS*S_L          + 
c     * EPS                  *S_TT  * COS(2*PHI_G)       + 
c     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
c     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
c     * )
      
c       IF(DVMPW.LT.0.0) DVMPW=0.
c       PRINT *,'DVMPW: S_tot ',S_TOT,COST,W,Q2,-del2,xb,phi_g,mesonmass,
c     *  eps,FLUXW(xb,Q2,E)

c       print *,'S_T ',  S_T,S_L,S_T   + EPS*S_L 
c       print *,'S_TT  ',  S_TT, EPS     *S_TT  * COS(2*PHI_G) 
c       print *, 'S_LT ',S_LT,SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G) 

      RETURN
      END

C====================================================================================c

      REAL FUNCTION SIGMA_T(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_T
      REAL          XSIGMA_T
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_T=JACR*XSIGMA_T(del2,x,Q2,E)
      ELSE
        SIGMA_T=0.0
      ENDIF
      RETURN
      END
C=============================================================      
      REAL FUNCTION SIGMA_L(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_L
      REAL          XSIGMA_L
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_L=JACR*XSIGMA_L(del2,x,Q2,E)
      ELSE
        SIGMA_L=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

      REAL FUNCTION SIGMA_U(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,Q2,E
      REAL          del2,x,JACG,JACR
      EXTERNAL      SIGMA_T,SIGMA_L,EPSILON
      REAL          SIGMA_T,SIGMA_L,EPSILON
      LOGICAL       CHECK_KINE

      SIGMA_U=0.
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) RETURN
      SIGMA_U=          SIGMA_T(COST,W,Q2,E)+
     * EPSILON(X,q2,e)*SIGMA_L(COST,W,Q2,E)

      RETURN
      END
C=============================================================      
      REAL FUNCTION SIGMA_TT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_TT
      REAL          XSIGMA_TT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_TT=JACR*XSIGMA_TT(del2,x,Q2,E)
      ELSE
        SIGMA_TT=0.0
      ENDIF
      RETURN
      END
      
C=============================================================      
      REAL FUNCTION SIGMA_LT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_LT
      REAL          XSIGMA_LT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_LT=JACR*XSIGMA_LT(del2,x,Q2,E)
      ELSE
        SIGMA_LT=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

       REAL FUNCTION SIGMA_LTP(COST,W,Q2,E)
       IMPLICIT NONE
       REAL                    COST,W,Q2,E

       SIGMA_LTP=0.

       RETURN
       END

C====================================================================================c

      LOGICAL FUNCTION CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      REAL COST,W,del2,xb,Q2,E,JACG,JACR
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
            

      real nu,W2
      real  xmin1,xmax1
      double precision y, e1, epsilon
      double precision tforward,tbackward
      logical valid_limits

      REAL AM(2)
      DATA AM/0.134976,0.547300/
      INTEGER K
      COMMON/PIETA/K

      MPI0=AM(K)

      CHECK_KINE=.FALSE.
      IF(.NOT.IEEE_IS_FINITE(COST).OR.
     *   .NOT.IEEE_IS_FINITE(W).OR.
     *   .NOT.IEEE_IS_FINITE(Q2).OR.
     *   .NOT.IEEE_IS_FINITE(E)) RETURN
      IF(Q2.LE.0..OR.E.LE.0..OR.W.LE.Mp+Mpi0) RETURN
      W2=W*W
      XB=Q2/(W2+Q2-Mp*Mp)
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      IF(W         .LE. Mp+Mpi0)              RETURN
      IF(W2        .GE. -Q2+2*Mp*E+Mp*Mp)     RETURN
      IF(Q2/(4*E*(E-NU)) .GE. 1.0)            RETURN
      IF(XB        .LE. Q2/(2*Mp*E))          RETURN
      IF(ABS(COST) .GE. 1.0)                  RETURN
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0

      CALL DVMP_T_LIMITS(Q2,xb,tforward,tbackward,valid_limits)
      IF(.NOT.valid_limits) RETURN
c     Use the same endpoints as XCHECK_KINE/tminq, then round only once.
      del2=-REAL(tforward+(tbackward-tforward)*
     *           (1D0-DBLE(COST))/2D0)
      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)

      IF(EPSILON.LE.0. .OR.EPSILON .GE. 1.0)         RETURN

C      jacobian
C      JAC=4*P1cm*P2CM*W*XB/(W2+Q2-Mp*Mp)
      JACG=2*W*XB/(W2+Q2-Mp*Mp)
      JACR=REAL((tbackward-tforward)/2D0)

      CHECK_KINE=.TRUE.
c      PRINT *,'Check_kine: ', W,xB,nu,qmod,E1cm,E2cm,del2,y,e1
      RETURN
      END


c*****************The end of Valery Model******************************
c
