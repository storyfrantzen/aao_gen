      REAL FUNCTION SIGMA_U_XBTPHI(XB_IN,Q2_IN,T_IN,PHI_IN)
c
c     Direct MC cross section interface:
c       sigma_u = SIGMA_U_XBTPHI(xB,Q2,t,phi)
c
c     T_IN is the Mandelstam t in GeV**2.  If a positive value is
c     supplied, it is interpreted as -t for convenience.
c     PHI_IN is in degrees, matching the existing MAID/AAO routines.
c
c     The four-argument function uses a configurable beam energy.  The
c     default is 10.6 GeV; call SET_SIGMA_U_XBTPHI_CONFIG first to change
c     the beam energy or the same physics options used by aao_norad.
c
      IMPLICIT NONE

      REAL XB_IN,Q2_IN,T_IN,PHI_IN
      REAL SIGMA_U_XBTPHI_E
      REAL CFG_EBEAM
      INTEGER CFG_THEORY,CFG_CHANNEL,CFG_RESONANCE

      COMMON /SIGU_XBTPHI_CFG/ CFG_EBEAM,CFG_THEORY,CFG_CHANNEL,
     1                         CFG_RESONANCE

      SIGMA_U_XBTPHI = SIGMA_U_XBTPHI_E(XB_IN,Q2_IN,T_IN,PHI_IN,
     1                                  CFG_EBEAM)

      RETURN
      END


      REAL FUNCTION SIGMA_U_XBTPHI_E(XB_IN,Q2_IN,T_IN,PHI_IN,
     1                               EBEAM_IN)
c
c     Same as SIGMA_U_XBTPHI, but with the beam energy supplied directly.
c     It converts xB,Q2,t to the variables used by aao_norad and calls
c     dsigma, so the returned SIGMA0 follows the MC event-generation
c     cross-section path.
c
      IMPLICIT NONE

      REAL XB_IN,Q2_IN,T_IN,PHI_IN,EBEAM_IN
      REAL MP,MPI0,MPIP
      REAL W2,W,NU,T_MANDEL,EP,SIN2
      REAL*8 THETA_E
      REAL E_PI_CM,P_PI_CM,QV_MAG_CM,NU_CM,COSTH
      REAL SIGMA0,SIGU,SIGT,SIGL,SIGI,SIGIP,ASYM_P
      REAL M_PI
      INTEGER EHEL
      REAL CFG_EBEAM
      INTEGER CFG_THEORY,CFG_CHANNEL,CFG_RESONANCE

      PARAMETER (MP=0.9382799)
      PARAMETER (MPI0=0.1349764)
      PARAMETER (MPIP=0.1395700)

      COMMON /SIGU_XBTPHI_CFG/ CFG_EBEAM,CFG_THEORY,CFG_CHANNEL,
     1                         CFG_RESONANCE

      SIGMA_U_XBTPHI_E = 0.0
      IF (XB_IN.LE.0.0 .OR. XB_IN.GE.1.0) RETURN
      IF (Q2_IN.LE.0.0 .OR. EBEAM_IN.LE.0.0) RETURN

      M_PI = MPI0
      IF (CFG_CHANNEL.NE.1) M_PI = MPIP

      W2 = MP**2 + Q2_IN*(1.0/XB_IN - 1.0)
      IF (W2.LE.(MP+M_PI)**2) RETURN
      W = SQRT(W2)

      NU = Q2_IN/(2.0*MP*XB_IN)
      IF (NU.GE.EBEAM_IN) RETURN

      EP = EBEAM_IN - NU
      IF (EP.LE.0.0) RETURN
      SIN2 = Q2_IN/(4.0*EBEAM_IN*EP)
      IF (SIN2.LE.0.0 .OR. SIN2.GE.1.0) RETURN
      THETA_E = 2.0D0*DASIN(DSQRT(DBLE(SIN2)))

      E_PI_CM = 0.5*(W2+M_PI**2-MP**2)/W
      P_PI_CM = E_PI_CM**2 - M_PI**2
      IF (P_PI_CM.LE.0.0) RETURN
      P_PI_CM = SQRT(P_PI_CM)

      QV_MAG_CM = ((W2+Q2_IN+MP**2)/(2.0*W))**2 - MP**2
      IF (QV_MAG_CM.LE.0.0) RETURN
      QV_MAG_CM = SQRT(QV_MAG_CM)
      NU_CM = (W2-MP**2-Q2_IN)/(2.0*W)

      T_MANDEL = T_IN
      IF (T_MANDEL.GT.0.0) T_MANDEL = -T_MANDEL
      COSTH = (T_MANDEL + Q2_IN - M_PI**2
     1        + 2.0*NU_CM*E_PI_CM)/(2.0*QV_MAG_CM*P_PI_CM)
      IF (ABS(COSTH).GE.1.0) RETURN

      EHEL = 0
      CALL DSIGMA(THETA_E,Q2_IN,W,COSTH,PHI_IN,CFG_THEORY,
     1            CFG_CHANNEL,CFG_RESONANCE,SIGMA0,SIGU,SIGT,SIGL,
     1            SIGI,SIGIP,ASYM_P,EHEL)

      SIGMA_U_XBTPHI_E = SIGMA0

      RETURN
      END


      SUBROUTINE SET_SIGMA_U_XBTPHI_CONFIG(EBEAM_IN,THEORY_OPT,
     1                                     CHANNEL_OPT,RESONANCE_OPT)
c
c     Configure SIGMA_U_XBTPHI defaults.
c       EBEAM_IN      : beam energy in GeV
c       THEORY_OPT    : same as aao_norad phys, usually 5
c       CHANNEL_OPT   : same as aao_norad epirea, 1=pi0, 3=pi+
c       RESONANCE_OPT : option passed through dsigma to maid_lee
c
      IMPLICIT NONE

      REAL EBEAM_IN
      INTEGER THEORY_OPT,CHANNEL_OPT,RESONANCE_OPT
      REAL CFG_EBEAM
      INTEGER CFG_THEORY,CFG_CHANNEL,CFG_RESONANCE

      COMMON /SIGU_XBTPHI_CFG/ CFG_EBEAM,CFG_THEORY,CFG_CHANNEL,
     1                         CFG_RESONANCE

      CFG_EBEAM = EBEAM_IN
      CFG_THEORY = THEORY_OPT
      CFG_CHANNEL = CHANNEL_OPT
      CFG_RESONANCE = RESONANCE_OPT

      RETURN
      END


      BLOCK DATA SIGU_XBTPHI_DEFAULTS
      IMPLICIT NONE

      REAL CFG_EBEAM
      INTEGER CFG_THEORY,CFG_CHANNEL,CFG_RESONANCE

      COMMON /SIGU_XBTPHI_CFG/ CFG_EBEAM,CFG_THEORY,CFG_CHANNEL,
     1                         CFG_RESONANCE

      DATA CFG_EBEAM /10.6/
      DATA CFG_THEORY /5/
      DATA CFG_CHANNEL /1/
      DATA CFG_RESONANCE /0/

      END
