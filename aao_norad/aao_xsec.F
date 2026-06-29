      PROGRAM AAO_XSEC
c
c     Command line interface:
c       ./aao_xsec -xB 0.3 -Q2 2.0 -t -0.2 -phi 30 -BeamEnergy 10.6
c     Returns the same dsigma SIGMA0 used by aao_norad before flux/Jacobian.
c
      IMPLICIT NONE

      INTEGER NARG,I
      CHARACTER*64 ARG
      CHARACTER*256 VAL
      CHARACTER*256 TABLE_PATH
      REAL XB,Q2,T,PHI,EBEAM,SIGU
      REAL SIGMA_U_XBTPHI_E
      INTEGER HAVE_XB,HAVE_Q2,HAVE_T,HAVE_PHI,HAVE_E
      INTEGER THEORY_OPT,CHANNEL_OPT,RESONANCE_OPT

      XB = 0.0
      Q2 = 0.0
      T = 0.0
      PHI = 0.0
      EBEAM = 10.6
      TABLE_PATH = ' '
      HAVE_XB = 0
      HAVE_Q2 = 0
      HAVE_T = 0
      HAVE_PHI = 0
      HAVE_E = 0
      THEORY_OPT = 5
      CHANNEL_OPT = 1
      RESONANCE_OPT = 0

      NARG = COMMAND_ARGUMENT_COUNT()
      IF (NARG.EQ.0) THEN
        CALL PRINT_AAO_XSEC_USAGE()
        STOP 1
      ENDIF

      I = 1
 10   IF (I.GT.NARG) GOTO 100
        CALL GET_COMMAND_ARGUMENT(I,ARG)

        IF (ARG.EQ.'-h' .OR. ARG.EQ.'--help') THEN
          CALL PRINT_AAO_XSEC_USAGE()
          STOP
        ENDIF

        IF (I.EQ.NARG) THEN
          WRITE(6,*) 'Missing value after ',ARG
          CALL PRINT_AAO_XSEC_USAGE()
          STOP 1
        ENDIF

        CALL GET_COMMAND_ARGUMENT(I+1,VAL)

        IF (ARG.EQ.'-xB' .OR. ARG.EQ.'--xB' .OR.
     1      ARG.EQ.'-xb' .OR. ARG.EQ.'--xb') THEN
          READ(VAL,*,ERR=900) XB
          HAVE_XB = 1
        ELSEIF (ARG.EQ.'-Q2' .OR. ARG.EQ.'--Q2' .OR.
     1          ARG.EQ.'-q2' .OR. ARG.EQ.'--q2') THEN
          READ(VAL,*,ERR=900) Q2
          HAVE_Q2 = 1
        ELSEIF (ARG.EQ.'-t' .OR. ARG.EQ.'--t') THEN
          READ(VAL,*,ERR=900) T
          HAVE_T = 1
        ELSEIF (ARG.EQ.'-phi' .OR. ARG.EQ.'--phi') THEN
          READ(VAL,*,ERR=900) PHI
          HAVE_PHI = 1
        ELSEIF (ARG.EQ.'-BeamEnergy' .OR.
     1          ARG.EQ.'--BeamEnergy' .OR.
     1          ARG.EQ.'-beam' .OR. ARG.EQ.'--beam') THEN
          READ(VAL,*,ERR=900) EBEAM
          HAVE_E = 1
        ELSEIF (ARG.EQ.'-theory' .OR. ARG.EQ.'--theory' .OR.
     1          ARG.EQ.'-phys' .OR. ARG.EQ.'--phys') THEN
          READ(VAL,*,ERR=900) THEORY_OPT
        ELSEIF (ARG.EQ.'-channel' .OR. ARG.EQ.'--channel' .OR.
     1          ARG.EQ.'-epirea' .OR. ARG.EQ.'--epirea') THEN
          READ(VAL,*,ERR=900) CHANNEL_OPT
        ELSEIF (ARG.EQ.'-resonance' .OR. ARG.EQ.'--resonance') THEN
          READ(VAL,*,ERR=900) RESONANCE_OPT
        ELSEIF (ARG.EQ.'-table' .OR. ARG.EQ.'--table') THEN
          TABLE_PATH = VAL
        ELSE
          WRITE(6,*) 'Unknown option ',ARG
          CALL PRINT_AAO_XSEC_USAGE()
          STOP 1
        ENDIF

        I = I + 2
        GOTO 10

 100  CONTINUE
      IF (HAVE_XB.EQ.0 .OR. HAVE_Q2.EQ.0 .OR. HAVE_T.EQ.0 .OR.
     1    HAVE_PHI.EQ.0 .OR. HAVE_E.EQ.0) THEN
        WRITE(6,*) 'Required options: -xB -Q2 -t -phi -BeamEnergy'
        CALL PRINT_AAO_XSEC_USAGE()
        STOP 1
      ENDIF

      CALL SET_SIGMA_U_XBTPHI_CONFIG(EBEAM,THEORY_OPT,CHANNEL_OPT,
     1                               RESONANCE_OPT)
      CALL SET_AAO_XSEC_QUIET(1)
      IF (TABLE_PATH.NE.' ') CALL SET_MAID07_TBL_PATH(TABLE_PATH)
      SIGU = SIGMA_U_XBTPHI_E(XB,Q2,T,PHI,EBEAM)

      WRITE(6,'(1PE16.8)') SIGU
      STOP

 900  CONTINUE
      WRITE(6,*) 'Could not read numeric value for ',ARG,' = ',VAL
      CALL PRINT_AAO_XSEC_USAGE()
      STOP 1

      END


      SUBROUTINE PRINT_AAO_XSEC_USAGE()
      IMPLICIT NONE

      WRITE(6,*) 'Usage:'
      WRITE(6,*) '  ./aao_xsec -xB X -Q2 Q2 -t T -phi PHI',
     1           ' -BeamEnergy E'
      WRITE(6,*) 'Optional: -theory/phys N -channel/epirea N',
     1           ' -resonance N'
      WRITE(6,*) 'Optional: -table /path/to/maid07-PPpi.tbl'
      WRITE(6,*) '  phi is in degrees; t may be negative t or',
     1           ' positive -t.'

      RETURN
      END


      SUBROUTINE SET_AAO_XSEC_QUIET(QUIET)
      IMPLICIT NONE

      INTEGER QUIET
      INTEGER AAO_XSEC_QUIET

      COMMON /AAO_XSEC_RUN_CFG/ AAO_XSEC_QUIET

      AAO_XSEC_QUIET = QUIET

      RETURN
      END


      SUBROUTINE SET_MAID07_TBL_PATH(TABLE_PATH)
      IMPLICIT NONE

      CHARACTER*(*) TABLE_PATH
      CHARACTER*256 MAID07_TBL_OVERRIDE
      INTEGER MAID07_TBL_SET

      COMMON /MAID07_TBL_CFG/ MAID07_TBL_OVERRIDE,MAID07_TBL_SET

      MAID07_TBL_OVERRIDE = TABLE_PATH
      MAID07_TBL_SET = 1

      RETURN
      END
