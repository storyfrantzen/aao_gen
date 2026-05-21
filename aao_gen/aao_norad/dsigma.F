      subroutine dsigma(the,q2,w,cscm,phicm,opt1,opt2,opt3
     *,sig0,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)

      implicit none

      real ki_mag,q2,w,cscm,phicm,kf_mag,s2
      real*8 the
      real fkt, w_thr, w2_thr, dvmpsig, maidsig
      real sig0,sigu,sigt,sigl,sigi,sigip,asym_p
      real nu,eps,eps1
      real q2_dummy, W_dummy, hcs_dummy, W2_dummy
      real Mp, mpi0, amp2, amu2, mesonmass, meta, mpip
      real xb_dummy, E_pi_cm, ppi_mag_cm, qv_mag_cm
      real nu_cm, t_dummy
      integer ehel
      real cthe
      integer opt1,opt2,opt3
      real pi
      parameter (Mp=0.93827)
      PARAMETER (MPI0=0.134976,MPIP=0.13957018,META=0.547853)
      PARAMETER (PI=3.1415926536)
      logical test1,test2,test3

      if (opt2.eq.1) then
        mesonmass = mpi0
      elseif (opt2.eq.3) then
        mesonmass = mpip
      elseif (opt2.eq.5) then
        mesonmass = meta
      endif

      nu 	= 0.5*(w**2 + q2 - 0.9382799**2)/0.9382799
      s2	= sin(0.5*the)**2
      ki_mag	= (nu+sqrt(q2/s2+nu**2))*0.5

      kf_mag 	= ki_mag-nu
      eps 	= 1. /(1+2.0*(1+nu*nu/q2)*tan(0.5*the)**2)

      test1	= ki_mag.lt.0.1.or.kf_mag.lt.0.01
      test3	= opt1.ge.4.and.(w.lt.1.1)

      if (test1.or.test3) then
c        print *, 'ABORT',ki_mag,kf_mag,q2,w
        sig0 = 0.
        sigu = 0.
        sigt = 0.
        sigl = 0.
        sigi = 0.
        sigip = 0.
        asym_p = 0.
        return
      endif

      if(opt1.eq.1) call aao(q2,w,eps,cscm,phicm,1,sig0,
     1             sigu,sigt,sigl,sigi)
      if(opt1.eq.2)
     1   call daresbury(q2,w,eps,cscm,phicm,1,sig0)
      if(opt1.eq.4.or.opt1.eq.5) then

c         if(w.gt.1.7.and.opt2.eq.1) then
CVPK         if(w.gt.2.and.opt2.eq.1) then
         if(w.gt.1.8 .and. opt2.eq.1) then
           amp2 = 0.93827**2
           amu2 = mesonmass**2
           q2_dummy = q2
           W_dummy  = w
           hcs_dummy = cscm
           W2_dummy = W_dummy**2
           xb_dummy = q2_dummy/(W2_dummy - amp2 + q2_dummy)
           E_pi_cm    = 0.5*(w2_dummy+amu2-amp2)/W_dummy
           ppi_mag_cm = E_pi_cm**2 - amu2
           ppi_mag_cm = sqrt(ppi_mag_cm)
           qv_mag_cm  = ((W2_dummy+q2_dummy+amp2)/2.0/W_dummy)**2
     1                   -amp2
           qv_mag_cm  = sqrt(qv_mag_cm)
           nu_cm  = (W2_dummy-amp2-q2_dummy)/(2*W_dummy)
           t_dummy =  -(-q2_dummy) - amu2 + 2*nu_cm*E_pi_cm
     1             - hcs_dummy * (2*qv_mag_cm*ppi_mag_cm)
          call dvmpw(hcs_dummy,w_dummy,q2_dummy,phicm*pi/180.,ki_mag,ehel,mesonmass,
     1           sig0,sigu,sigl,sigt,sigi,sigip)

c  divided individual structure functions by fkt to be consistent with structure functions from MAID
            fkt  = 2*W_dummy*ppi_mag_cm/(W2_dummy - amp2)

            sigu = sigu/fkt
            sigl = sigl/fkt
            sigt = sigt/fkt
            sigi = sigi/fkt
            sigip = sigip/fkt

c  the factor fkt does not need to be applied for dvmp parameterization
c           fkt  = 2*W_dummy*ppi_mag_cm/(W2_dummy - amp2)
c           sig0 = fkt*sig0

         else
           amp2 = 0.93827**2
           amu2 = mesonmass**2
           q2_dummy = q2
           W_thr  = 1.8
           hcs_dummy = cscm
           W2_thr = W_thr**2
           xb_dummy = q2_dummy/(W2_thr - amp2 + q2_dummy)
           E_pi_cm    = 0.5*(w2_thr+amu2-amp2)/W_thr
           ppi_mag_cm = E_pi_cm**2 - amu2
           ppi_mag_cm = sqrt(ppi_mag_cm)
           qv_mag_cm  = ((W2_thr+q2_dummy+amp2)/2.0/W_thr)**2
     1                   -amp2
           qv_mag_cm  = sqrt(qv_mag_cm)
           nu_cm  = (W2_thr-amp2-q2_dummy)/(2*W_thr)
           t_dummy =  -(-q2_dummy) - amu2 + 2*nu_cm*E_pi_cm
     1             - hcs_dummy * (2*qv_mag_cm*ppi_mag_cm)
           call dvmpw(hcs_dummy,w_thr,q2_dummy,phicm*pi/180.,ki_mag,ehel,mesonmass,
     1           dvmpsig,sigu,sigl,sigt,sigi,sigip)
           call maid_lee(q2,w_thr,eps,cscm,phicm,opt1,opt2,opt3,
     1           maidsig,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
           call maid_lee(q2,w,eps,cscm,phicm,opt1,opt2,opt3,
     1           sig0,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
           sig0 = sig0*dvmpsig/maidsig
         endif
      endif
c      print *, 'DSIGMA: ',sig0,sigu,sigl,sigt,sigi,sigip,asym_p,ehel
      end
