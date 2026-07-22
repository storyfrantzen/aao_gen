      program aao_norad
c
c     This program makes an n-tuple that can be used with Paw and GSIM to
c     make distributions of energies, angles, resonance mass resulting 
c     from non-radiative pion production on a proton.
c
c     The n-tuple contains the true hadronic invariant mass (W), the components 
c     of the proton momentum (PPX, PPY, PPZ),
c     the proton energy (EP), the pion momentum (PPIX, PPIY, PPIZ) and
c     pion energy (EPI), the  angles for the hadronic
c     decay in the hadronic frame (CSTCM, PHICM), the missing mass (MM).
c
c     Resolution of detectors is not folded into the results.  If this is 
c     desired it
c     should be done with a second program that can operate on the n-tuple
c     and make a new version.

      implicit none
      
#include "bcs.inc"
#include "mc.inc"
#include "names.inc"

      COMMON/ALPHA/ ALPHA,PI,MP,MN,MPI,MEL,WG,EPIREA
      common /random/idum
      
      real*8 alpha,pi,mp,mn,mpi,mel,wg,T0,th0,theta
      real*8 es,ep,ps,pp,rs,rp,nu
      real*8 cst0,snt0
      real*8 sig_tot,sig_sum
      real*8 ephot
      real lundmass(4)
     
      real pi0costh, pi0phi
      real g1px, g1py, g1pz, g1p
      real g2px, g2py, g2pz, g2p
      real pi0bx, pi0by, pi0bz
      real pi0bp, pi0b2, pi0gam, pi0gam2
      real csthcm, csthcm_max, cstk
      real epi, epmax, eprng, eprot, epw
      real ebeam,ep_min,ep_max,ep_test, evirt
      real*8 events
      real fmcall
      real gam_flux
      real itime1, itime2
      real jacob
      real analysis_volume,xb,xb_min,xb_max,xbrng
      real minus_t,mt_min,mt_max,mtrng
      real phi_min,phi_max,phirng,q2rng_analysis
      real q0cm,qcm,epicm,picm,tcm,cmden
      real myran,mm2, mpi0,mpip,mpi_s
      real ppx,ppy,ppz,px,py
      real ppix,ppiy,ppiz
      real phicm,phicm_max,phik,phir
      real q0,q2,q2_min,q2_max,q2max
      real qsq,qvecx,qvecz
      real ran,rotc,rots,rtest
      real s,sigi,sigl,sigt,sigu,sigma0,signr,sig_ratio,sigr_max,stest
      real sigip, asym_p
      real t_elapse
      real*8 tries
      real*8 uq2, uq2_min,uq2_max,uq2rng
      real w2,wreal,w_max
      real t_targ
      real*8 sig_int
      real vertex_x, vertex_y, vertex_z


      integer seed_source
      integer sampling_mode,input_status
      integer env_status, env_len
      integer boso
      integer phys,epirea
      integer i,iext1,iext2,iquest
      integer j,jj
      integer mcall,mcall_max
      integer ncycle
      integer nprint,ntell,ntold
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer get_spin
      integer flag_ehel
      integer ehel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*8 idum, count_rate, count_max
      integer*8 job_seed, attempt_seed
      integer*8 ntries
      
      real*4 ntp(24)

clcs     common /pawc/h(1000000)
clcs     common/quest/iquest(100)
      common/ntup/ntp
      
      integer h,nevent,nmax,istat,icycle
      
      character*132 filerz,file_out,file_sum,file_bos,file_lund,file_norm
      character*132 file_kin
      character*8 recname
      character*32 env_value

      DATA PI   /3.1415926/
      DATA MPIP /.1395/
      DATA MPI0 /.1349/
      DATA MP   /.93828/
      DATA MN   /.93957/
      DATA MEL  /.511E-3/
      
      data filerz	/'aao_norad.rz'/
      data file_out	/'aao_norad.out'/
      data file_sum	/'aao_norad.sum'/
c      data file_bos	/'aao_norad.evt'/
      data file_lund     /'aao_norad.lund'/
      data file_norm     /'aao_norad.norm'/
      data file_kin      /'aao_norad.kin'/

      alpha=1/137.
      
c     phys:     Physics model (1=AO, 4=MAID98, 5=MAID2000)
c     flag_ehel:  0=no polarized electron, 1= polarized electron  
c     npart: 	number of particles in BOS banks - 2(e-,h+) 3(e-,h+,h0)
c     epirea: 	final state hadron (1-pi0, 3-pi+)
c     ebeam: 	incident electron energy (GeV)
c     q2_*:	q**2 limits (GeV**2)
c     ep_*:	scattering electron limits (GeV)
c     nmax:	number of output events
c     fmcall:   scale factor for sigm_max
c     boso:	1=bos output, 0=no bos output
c     seed_source: 0= use unix timestamp from machine time to generate seed, otherwise use given value as seed

      read(5,*) phys,flag_ehel,npart,epirea,ebeam,q2_min,q2_max,ep_min,
     1ep_max,nmax,fmcall,boso,seed_source

c     Optional sampling extension.  EOF preserves the historical input format.
c       0 = legacy coordinates and full legacy domain
c       1 = uniform Q2, xB, -t, phi analysis coordinates
c       2 = legacy coordinates restricted to the same analysis bounds
c       3 = uniform 1/Q2, xB, -t, phi analysis coordinates
      sampling_mode = 0
      xb_min = 0.
      xb_max = 1.
      mt_min = 0.
      mt_max = 100.
      phi_min = 0.
      phi_max = 360.
      read(5,*,iostat=input_status) sampling_mode
      if (input_status .eq. 0 .and. sampling_mode .ne. 0) then
         read(5,*,iostat=input_status) xb_min,xb_max
         if (input_status .ne. 0) stop 'Missing xB sampling bounds'
         read(5,*,iostat=input_status) mt_min,mt_max
         if (input_status .ne. 0) stop 'Missing -t sampling bounds'
         read(5,*,iostat=input_status) phi_min,phi_max
         if (input_status .ne. 0) stop 'Missing phi sampling bounds'
      endif
      if (sampling_mode .lt. 0 .or. sampling_mode .gt. 3) then
         stop 'sampling_mode must be 0, 1, 2, or 3'
      endif
      if (sampling_mode .ne. 0) then
         if (xb_min .le. 0. .or. xb_max .le. xb_min) then
            stop 'Invalid xB sampling bounds'
         endif
         if (mt_min .lt. 0. .or. mt_max .le. mt_min) then
            stop 'Invalid -t sampling bounds'
         endif
         if (phi_min .lt. 0. .or. phi_max .gt. 360. .or.
     *       phi_max .le. phi_min) then
            stop 'Invalid phi sampling bounds'
         endif
      endif
    
      if (npart .ne. 3) npart=2

      id(1)		= 3		!Geant ID, e-
      pdgid(1)  	= 11		!PDG ID, e-
      q(1) 		= -1

      if (epirea.eq.1) then
         mpi		= mpi0
         id(2)	   	= 14		!Geant ID, proton
         pdgid(2) 	= 2212		!PDG ID, proton
         q(2) 		= 1
         if (npart.eq.3) then
            id(3)     	=  7            !Geant ID, pi-zero
            pdgid(3)  	= 111		!PDG ID, pi-zero
            q(3) 	= 0
         endif
      elseif (epirea.eq.2) then
         mpi		= mpip
         id(2)     	= 8             !Geant ID, pi-minus
         pdgid(2)  	= 212		!PDG ID, pi-minus
         q(2) 		= -1
         if (npart.eq.3) then
            id(3)	= 14		!Geant ID, proton
            pdgid(3)	= 2212		!PDG ID, proton
            q(3) 	= 1
         endif
      elseif (epirea.eq.3) then
         mpi		= mpip
         id(2)     	= 8             !Geant ID, pi-plus
         pdgid(2)  	= 211		!PDG ID, pi-plus
         q(2) 		= 1
         if (npart.eq.3) then
            id(3)	= 13		!Geant ID, neutron
            pdgid(3)	= 2112		!PDG ID, neutron
            q(3) 	= 0
         endif
      else
         stop
      endif

      mpi_s		= mpi		!Single precision pion mass
      wg		= mp+mpi+.0005	!minimum hadron mass for pion production

      es		= ebeam
      ps		= sqrt(es**2-mel**2)
      rs		= ps/es
      
      s			= .25		!cut off q2 at 60 deg sin(60)=0.25
      q2max		= 4.*ebeam**2*s/(1.+2.*ebeam*s/mp)

      iext1		= 48
      iext2		= 48

      if (q2_max.gt.q2max) q2_max=q2max
      
      uq2_min		= 1/q2_max
      uq2_max		= 1/q2_min
      uq2rng		= uq2_max-uq2_min
      q2rng_analysis     = q2_max-q2_min
      xbrng              = xb_max-xb_min
      mtrng              = mt_max-mt_min
      phirng             = (phi_max-phi_min)*pi/180.
      if (sampling_mode .eq. 3) then
         analysis_volume = 2.*pi*uq2rng*xbrng*mtrng*phirng
      else
         analysis_volume = 2.*pi*q2rng_analysis*xbrng*mtrng*phirng
      endif

      epmax		= es-(wg**2+q2_min-mp**2)/2./mp    
      if (ep_max.lt.epmax) epmax=ep_max     
      eprng		= epmax-ep_min

      nprint		= nmax/25
      ncycle		= 0
 1    mcall_max		= 0
      ntold		= 0
	     
      call timex(itime1)		!Initial seed from internal clock
      call system_clock(idum, count_rate, count_max)

      if(seed_source.eq.0) then
      job_seed = 0
      attempt_seed = 0
      call get_environment_variable('SWIF_JOB_ID', env_value,
     *     length=env_len, status=env_status)
      if (env_status .eq. 0 .and. env_len .gt. 0) then
         read(env_value(1:env_len),*,err=131) job_seed
      endif
 131  continue
      call get_environment_variable('SWIF_JOB_ATTEMPT_ID', env_value,
     *     length=env_len, status=env_status)
      if (env_status .eq. 0 .and. env_len .gt. 0) then
         read(env_value(1:env_len),*,err=132) attempt_seed
      endif
 132  continue
      idum = idum + 104729_8*job_seed + 1009_8*attempt_seed
      idum	= -idum
      else
      idum = seed_source
      endif

      Print *, 'Seed idum is', idum, ' '
    
      cstk		= myran(idum)
      write(*,*) cstk
      itime2		= itime1 
      nevent		= 0
      t_elapse		= 0.
      ntries		= 0
      sig_int		= 0.
      sig_tot		= 0.

      open(unit=12,file=file_out)
      open(unit=14,file=file_sum)
      open(unit=13,file=file_lund)
      open(unit=16,file=file_kin)
      write(16,*)'# event Q2 xB minus_t phi_deg signr jacobian'

clcs      call hlimit(1000000)
clcs      iquest(10)	= 65000
clcs      call hropen(1,'AONORAD',filerz,'NQ',8191,istat)
clcs      call hbset('BSIZE',8176,istat)
clcs      call hbnt(10,'AONORAD',' ')
clcs      call hbname(10,' ',0,'$clear')
clcs      call hbname(10,'AONORAD',ntp(1),'ES:R,EP:R,THETE:R,
clcs     *W:R,PPX:R,PPY:R,PPZ:R,EPROT:R,PPIX:R,PPIY:R,PPIZ:R,
clcs     *EPI:R,CSTCM:R,PHICM:R,MM:R,QX:R,QZ:R,Q0:R,CSTHE:R,
clcs     *HEL:R,ASYM:R,Q2:R,JAC:R,GAMF:R')

c      bosout = file_bos
c      recname = 'MCEVENT'
c      if (boso.eq.1) call bos_ini(recname)

      write(12,*)' AO Calculation of Single Pion Production'
      write(12,*)' Seed:',idum
      write(12,*)' Epirea (1 for pi0, 2 for pi-,  3 for pi+) =',epirea
      write(12,*)' Target thickness =',t_targ,' (r.l.)'
      write(12,*)' Incident electron energy =',ebeam,' GeV'
      write(12,*)' Electron Q**2 limits:',q2_min,q2_max
      write(12,*)' Lower and upper limit for scattered electron',
     * ' (GeV):',ep_min,epmax
      write(12,*)' Sampling mode:',sampling_mode
      if (sampling_mode .ne. 0) then
         write(12,*)' Analysis xB bounds:',xb_min,xb_max
         write(12,*)' Analysis -t bounds:',mt_min,mt_max
         write(12,*)' Analysis phi bounds (deg):',phi_min,phi_max
      endif

c     Do a preliminary calculation to estimate the maximum value
c     of the integrand

 10   q2		= q2_min
      sigr_max		= 0.     
      ehel              = 0  

      if(flag_ehel.eq.1) ehel = 1
       
      do i=1,10000
        if (sampling_mode .eq. 1 .or. sampling_mode .eq. 3) then
          if (sampling_mode .eq. 3) then
            uq2 = uq2_min+uq2rng*myran(idum)
            q2 = 1./uq2
          else
            q2 = q2_min+q2rng_analysis*myran(idum)
          endif
          xb = xb_min+xbrng*myran(idum)
          q0 = q2/(2.*mp*xb)
          ep = es-q0
        else
          uq2 = uq2_min+uq2rng*myran(idum)
          q2 = 1./uq2
          ep = epmax-eprng*myran(idum)
          q0 = es-ep
          if (q0 .le. 0.) cycle
          xb = q2/(2.*mp*q0)
        endif
        if (ep .lt. ep_min .or. ep .gt. epmax) cycle
        s = q2/(4.*es*ep)
        if (s .le. 0. .or. s .gt. .5) cycle
        th0 = 2.*asin(sqrt(s))
        theta = th0*180./pi
        t0 = th0
        snt0 = sin(th0)
        cst0 = cos(th0)
        w2 = mp**2+2.*mp*q0-q2
        if (w2 .le. (wg+0.002)**2) cycle
        epw = sqrt(w2)
        q0cm = (w2-mp**2-q2)/(2.*epw)
        if (q0cm**2+q2 .le. 0.) cycle
        qcm = sqrt(q0cm**2+q2)
        epicm = (w2+mpi**2-mp**2)/(2.*epw)
        if (epicm**2-mpi**2 .le. 0.) cycle
        picm = sqrt(epicm**2-mpi**2)
        cmden = 2.*qcm*picm
        if (cmden .le. 0.) cycle
        if (sampling_mode .eq. 1 .or. sampling_mode .eq. 3) then
          minus_t = mt_min+mtrng*myran(idum)
          tcm = -minus_t
          csthcm = (tcm-mpi**2+q2+2.*q0cm*epicm)/cmden
          phicm = phi_min+(phi_max-phi_min)*myran(idum)
          if (abs(csthcm) .gt. 1.) cycle
          if (sampling_mode .eq. 3) then
            jacob = q2**3/(8.*mp*xb**2*es*ep*qcm*picm)
          else
            jacob = q2/(8.*mp*xb**2*es*ep*qcm*picm)
          endif
        else
          csthcm = 2.*(myran(idum)-0.5)
          phicm = 360.*myran(idum)
          tcm = mpi**2-q2-2.*q0cm*epicm+cmden*csthcm
          minus_t = -tcm
          if (sampling_mode .eq. 2) then
            if (xb .lt. xb_min .or. xb .ge. xb_max) cycle
            if (minus_t .lt. mt_min .or. minus_t .ge. mt_max) cycle
            if (phicm .lt. phi_min .or. phicm .ge. phi_max) cycle
          endif
          jacob = q2**2/(2.*es*ep)
        endif
        if (jacob .le. 0.) cycle
        evirt = (epw**2-mp**2)/2./mp
        qsq = q2
        gam_flux = alpha/4/pi**2*ep/es*evirt/qsq*
     +    (2. + 1./tan(th0/2)**2/(1.+q0**2/qsq))
        call dsigma(th0,qsq,epw,csthcm,phicm,phys,epirea,0,sigma0
     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
        signr		= gam_flux*sigma0*jacob
        if (signr.gt.sigr_max) then
          sigr_max	= signr
          csthcm_max 	= csthcm
          phicm_max	= phicm
          w_max         = epw
        endif
      enddo

      write(6,*)'sigr_max,csthcm_max,phicm_max,w_max',
     *           sigr_max,csthcm_max,phicm_max,w_max
      write(12,*)'sigr_max,csthcm_max,phicm_max,w_max',
     *            sigr_max,csthcm_max,phicm_max,w_max

      if (sigr_max .le. 0.) then
         stop 'Maximum scan found no physical positive-cross-section point'
      endif
      sigr_max		= sigr_max*fmcall
      write(6,*) 'sigr_max changed to',sigr_max
      write(12,*)'sigr_max changed to',sigr_max

c   %%%%%%%%%%%%%%%%%%% Main Calculation  %%%%%%%%%%%%%%%%%%%%%%%
      ehel = 0
 20   continue
      ntries		= ntries+1
      if(flag_ehel.eq.1) ehel = get_spin(idum)
      es		= ebeam
      if (sampling_mode .eq. 1 .or. sampling_mode .eq. 3) then
         if (sampling_mode .eq. 3) then
            uq2 = uq2_min+uq2rng*myran(idum)
            q2 = 1./uq2
         else
            q2 = q2_min+q2rng_analysis*myran(idum)
         endif
         xb = xb_min+xbrng*myran(idum)
         q0 = q2/(2.*mp*xb)
         ep = es-q0
      else
         uq2 = uq2_min+uq2rng*myran(idum)
         q2 = 1./uq2
         ep = epmax-eprng*myran(idum)
         q0 = es-ep
         if (q0 .le. 0.) go to 20
         xb = q2/(2.*mp*q0)
      endif
      
      if (ep .lt. ep_min .or. ep .gt. epmax) go to 20

      s			= q2/4/es/ep

      if (s .le. 0. .or. s .gt. .5) go to 20
      
      th0		= 2.*asin(sqrt(s))
      theta		= th0*180./pi
      T0		= th0
      snt0		= sin(th0)
      cst0		= cos(th0)
      pp		= sqrt(ep**2-mel**2)
      qsq		= q2
      
      if (qsq .le. 0.)then
         write(6,*)' Unphysical Q**2 =',qsq
         go to 20
      endif

      qvecx             =   -pp*snt0
      qvecz             = ps-pp*cst0

      w2		= mp**2+2*mp*q0-q2
      
      if (w2 .lt. mp**2) go to 20
      
      epw		= sqrt(w2)
      
      if (epw .lt. wg+0.002)go to 20

      q0cm = (w2-mp**2-q2)/(2.*epw)
      if (q0cm**2+q2 .le. 0.) go to 20
      qcm = sqrt(q0cm**2+q2)
      epicm = (w2+mpi**2-mp**2)/(2.*epw)
      if (epicm**2-mpi**2 .le. 0.) go to 20
      picm = sqrt(epicm**2-mpi**2)
      cmden = 2.*qcm*picm
      if (cmden .le. 0.) go to 20

      if (sampling_mode .eq. 1 .or. sampling_mode .eq. 3) then
         minus_t = mt_min+mtrng*myran(idum)
         tcm = -minus_t
         csthcm = (tcm-mpi**2+q2+2.*q0cm*epicm)/cmden
         phicm = phi_min+(phi_max-phi_min)*myran(idum)
         if (abs(csthcm) .gt. 1.) go to 20
      else
         csthcm = -1.+2.*myran(idum)
         phicm = 360.*myran(idum)
         tcm = mpi**2-q2-2.*q0cm*epicm+cmden*csthcm
         minus_t = -tcm
         if (sampling_mode .eq. 2) then
            if (xb .lt. xb_min .or. xb .ge. xb_max) go to 20
            if (minus_t .lt. mt_min .or. minus_t .ge. mt_max) go to 20
            if (phicm .lt. phi_min .or. phicm .ge. phi_max) go to 20
         endif
      endif
      call dsigma(th0,qsq,epw,csthcm,phicm,phys,epirea,0,sigma0
     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)

      if (sigma0.le.0) go to 20

      evirt = (epw**2-mp**2)/2./mp
      gam_flux	= alpha/4/pi**2*ep/es*evirt/qsq*
     +  (2. + 1./tan(th0/2)**2/(1.+q0**2/qsq))

      signr = sigma0*gam_flux
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xsec      = sigma0
      asym      = asym_p
      hel       = ehel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
28    if (sampling_mode .eq. 1) then
         jacob = q2/(8.*mp*xb**2*es*ep*qcm*picm)
      else if (sampling_mode .eq. 3) then
         jacob = q2**3/(8.*mp*xb**2*es*ep*qcm*picm)
      else
         jacob = q2**2/(2.*es*ep)
      endif
      signr = signr*jacob
      sig_ratio = signr/sigr_max
      
c     Choose the number of times, mcall, to call the routine used
c     to calculate kinematic quantities for the n-tuple.

      rtest = myran(idum)
      mcall = sig_ratio
      stest = sig_ratio-mcall
      
      if (stest .gt. rtest) mcall = mcall + 1
      if (mcall .gt. mcall_max) mcall_max=mcall

c     Increment ntries by the value of mcall. It is also incremented by
c     1 at the beginning of the loop
c     if (mcall .gt. 0)ntries=ntries+mcall-1
c     if mcall .gt. 0 generate mcall n-tuple events.

      if (mcall.ge.0) sig_tot = sig_tot+signr
      if (mcall .le. 0) go to 30
      
      do j=1,mcall
        ephot = 0.
        cstk = 1.
        phik = 0.

        call missm(ebeam,es,ep,th0,ephot,cstk,phik,mpi_s,ppx,
     *  ppy,ppz,eprot,ppix,ppiy,ppiz,epi,csthcm,phicm,wreal,mm2)

        if(mm2 .eq. 0. .or. ep .lt. ep_min) go to 24

        w2 = mp**2+2*mp*(ebeam-ep)-2*es*ep*(1-cos(th0))+2.*mel**2
        epw = sqrt(w2)
        
        ntp(1) = es
        ntp(2) = ep
        ntp(3) = theta
        ntp(4) = epw
        ntp(5) = ppx
        ntp(6) = ppy
        ntp(7) = ppz
        ntp(8) = eprot
        ntp(9) = ppix
        ntp(10) = ppiy
        ntp(11) = ppiz
        ntp(12) = epi
        ntp(13) = csthcm
        ntp(14) = phicm
        ntp(15) = mm2
        ntp(16) = qvecx
        ntp(17) = qvecz
        ntp(18) = q0
        ntp(19) = cst0
        ntp(20) = ehel
        ntp(21) = asym_p
        ntp(22) = q2
        ntp(23) = jacob
        ntp(24) = gam_flux

clcs        call hfnt(10)
        
        nevent = nevent+1
        write(16,*)nevent,q2,xb,minus_t,phicm,signr,jacob
        
        do jj = 1,npart
          v(jj,1) = 0.
          v(jj,2) = 0.
          v(jj,3) = 0.
        enddo
         
        phir = 2.*pi*myran(idum) !Rotate all momentum around beam line randomly
        rotc = cos(phir)
        rots = sin(phir)

        px = ep*sin(theta*pi/180.) !Scattered electron
        py = 0.
        p(1,1) = px*rotc+py*rots
        p(1,2) = py*rotc-px*rots
        p(1,3) = ep*cos(theta*pi/180.)
        p(1,4) = ep

        if (epirea .eq. 1) then !Charged hadron
          p(2,1) = ppx*rotc+ppy*rots
          p(2,2) = ppy*rotc-ppx*rots
          p(2,3) = ppz
          p(2,4) = eprot
          if (npart .eq. 3)then !Scattered pion
            p(3,1) = ppix*rotc+ppiy*rots
            p(3,2) = ppiy*rotc-ppix*rots
            p(3,3) = ppiz
            p(3,4) = epi
          endif
        else
          p(2,1) = ppix*rotc+ppiy*rots !Scattered pion
          p(2,2) = ppiy*rotc-ppix*rots
          p(2,3) = ppiz
          p(2,4) = epi
          if (npart .eq. 3)then
            p(3,1) = ppx*rotc+ppy*rots
            p(3,2) = ppy*rotc-ppx*rots
            p(3,3) = ppz
            p(3,4) = eprot
          endif
        endif


         vertex_x = 0
         vertex_y = 0
         vertex_z = 5.0*(myran(idum)-0.5) - 3.0
         lundmass(1) = MEL
         lundmass(4) = 0

         if ((epirea .eq. 1).or.(epirea.eq.5)) then
            lundmass(2) = MP
            lundmass(3) = MPI

            pi0costh = 2.*(myran(idum)-0.5)
            pi0phi = 2.*pi*myran(idum)

            g1p = MPI/2.0
            g1px = g1p*sqrt(1-pi0costh*pi0costh)*cos(pi0phi)
            g1py = g1p*sqrt(1-pi0costh*pi0costh)*sin(pi0phi)
            g1pz = g1p*pi0costh
            g2p = g1p
            g2px = -g1px
            g2py = -g1py
            g2pz = -g1pz

            pi0bx = p(3,1)/p(3,4)
            pi0by = p(3,2)/p(3,4)
            pi0bz = p(3,3)/p(3,4)

            pi0b2 = pi0bx*pi0bx + pi0by*pi0by + pi0bz*pi0bz
            pi0gam = 1.0/(sqrt(1.0 - pi0b2))
            if(pi0b2.gt.0) then
              pi0gam2 = (pi0gam-1.0)/pi0b2
            else
              pi0gam2 = 0.0
            endif

            pi0bp = pi0bx*g1px + pi0by*g1py + pi0bz*g1pz
            g1px = g1px + pi0gam2*pi0bp*pi0bx + pi0gam*pi0bx*g1p
            g1py = g1py + pi0gam2*pi0bp*pi0by + pi0gam*pi0by*g1p
            g1pz = g1pz + pi0gam2*pi0bp*pi0bz + pi0gam*pi0bz*g1p
            g1p = pi0gam*(g1p + pi0bp)

            pi0bp = pi0bx*g2px + pi0by*g2py + pi0bz*g2pz
            g2px = g2px + pi0gam2*pi0bp*pi0bx + pi0gam*pi0bx*g2p
            g2py = g2py + pi0gam2*pi0bp*pi0by + pi0gam*pi0by*g2p
            g2pz = g2pz + pi0gam2*pi0bp*pi0bz + pi0gam*pi0bz*g2p
            g2p = pi0gam*(g2p + pi0bp)
            
            write(13,*) "4 1 1 0 0 ",q2/(2.0*MP*(ebeam-p(1,4))), 0, epw,
     *               q2, (ebeam-p(1,4))
c     *               phicm,sigma0,sigu,sigl,sigt,sigi,sigip
            write(13,*) 1, q(1), 1, pdgid(1), 0, 0,
     *               p(1,1), p(1,2), p(1,3), p(1,4),
     *               lundmass(1), vertex_x, vertex_y, vertex_z
            write(13,*) 2, q(2), 1, pdgid(2), 0, 0,
     *               p(2,1), p(2,2), p(2,3), p(2,4),
     *               lundmass(2), vertex_x, vertex_y, vertex_z
            write(13,*) 3, 0, 1, 22, 0, 0,
     *               g1px, g1py, g1pz, g1p, 0,
     *               vertex_x, vertex_y, vertex_z
            write(13,*) 4, 0, 1, 22, 0, 0,
     *               g2px, g2py, g2pz, g2p, 0,
     *               vertex_x, vertex_y, vertex_z

         else
           
            lundmass(2) = MPIP                                                                                               
            if(epirea.eq.2) lundmass(3) = MP
            if(epirea.eq.3) lundmass(3) = MN

            write(13,*) "3 1 1 0 0 ",q2/(2.0*MP*(ebeam-p(1,4))), 0, epw,
     *               q2, (ebeam-p(1,4))
c     *               phicm,sigma0,sigu,sigl,sigt,sigi,sigip
            write(13,*) 1, q(1), 1, pdgid(1), 0, 0,
     *               p(1,1), p(1,2), p(1,3), p(1,4),
     *               lundmass(1), vertex_x, vertex_y, vertex_z
            write(13,*) 2, q(2), 1, pdgid(2), 0, 0,
     *               p(2,1), p(2,2), p(2,3), p(2,4),
     *               lundmass(2), vertex_x, vertex_y, vertex_z
            write(13,*) 3, q(3), 1, pdgid(3), 0, 0,
     *               p(3,1), p(3,2), p(3,3), p(3,4),
     *               lundmass(3), vertex_x, vertex_y, vertex_z
         endif

c        if (boso.eq.1) call bos_out	! Pack the BOS banks and write out to file
 24     continue
      enddo
      
      if (nprint .eq. 0) then
        ntell = 0
      else 
        ntell		= nevent/nprint-ntold
      endif      

      if (ntell .gt. 0)then
        write(6,*)' ntries, nevent, mcall_max: '
     *             ,ntries,nevent,mcall_max
        write(12,*)' ntries, nevent, mcall_max: '
     *              ,ntries,nevent,mcall_max
        events		= nevent
        tries		= ntries
        sig_int		= events/tries
          
        if (sampling_mode .eq. 1 .or. sampling_mode .eq. 3) then
          sig_int = sig_int*sigr_max*analysis_volume
          sig_sum = sig_tot*analysis_volume/tries
        else
c         hadron phase space = 4*pi
c         electron phase space =2*pi*uq2rng*eprng
          sig_int = sig_int*sigr_max*(4.*pi)*
     *      (2.*pi*uq2rng*eprng)
          sig_sum = sig_tot*(4.*pi)*
     *      (2.*pi*uq2rng*eprng)/tries
        endif

        write(6,*)' Integrated cross section =',sig_int,sig_sum,
     *    ' micro-barns'
        write(6,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
        write(12,*)' Integrated cross section =',sig_int,sig_sum
     *    ,' micro-barns'
        write(12,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
        call timex(itime2)
        t_elapse	= t_elapse+INT(itime2-itime1)
        itime1		= itime2
        write(6,*) ' Elapsed CPU time = ',t_elapse/60,' minutes'
        write(12,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
        ntold		= ntold+1
      endif

      if (nevent .ge. nmax) go to 50
 30   go to 20

 50   continue
      events = nevent
      tries = ntries
      sig_int = events/tries
      if (sampling_mode .eq. 1 .or. sampling_mode .eq. 3) then
         sig_int = sig_int*sigr_max*analysis_volume
         sig_sum = sig_tot*analysis_volume/tries
      else
         sig_int = sig_int*sigr_max*(4.*pi)*
     *      (2.*pi*uq2rng*eprng)
         sig_sum = sig_tot*(4.*pi)*
     *      (2.*pi*uq2rng*eprng)/tries
      endif
clcs      call hrout(0,icycle,' ')
clcs      call hrend('AONORAD')
clcs      close(unit=1)
c      if (boso.eq.1) call bos_end(recname)
      
      write(14,*)' AO Calculation of Single Pion Production',
     + ' without radiative corrections'
      write(14,*)' Epirea (1 for pi0, 3 for pi+) =',epirea
      write(14,*)' Target thickness =',t_targ,' (r.l.)'
      write(14,*)' Incident electron energy =',ebeam,' GeV'
      write(14,*)' Electron Q**2 limits:',q2_min,q2_max
      write(14,*)' Lower and upper limit for scattered electron',
     * ' energy(GeV):',ep_min,epmax
      write(14,*)'sigr_max,csthcm_max,phicm_max',
     * sigr_max,csthcm_max,phicm_max
      write(14,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max

      write(14,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
      write(14,*)' Integrated cross section =',sig_int,sig_sum
     *    ,' micro-barns'
      write(14,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
      open(unit=15,file=file_norm)
      write(15,*)'generator=aao_norad'
      write(15,*)'integrated_cross_section_units=micro-barns'
      write(15,*)'sig_int=',sig_int
      write(15,*)'sig_sum=',sig_sum
      write(15,*)'events=',events
      write(15,*)'ntries=',ntries
      write(15,*)'nevent=',nevent
      write(15,*)'mcall_max=',mcall_max
      write(15,*)'sigr_max=',sigr_max
      write(15,*)'ebeam=',ebeam
      write(15,*)'q2_min=',q2_min
      write(15,*)'q2_max=',q2_max
      write(15,*)'ep_min=',ep_min
      write(15,*)'ep_max_effective=',epmax
      write(15,*)'epirea=',epirea
      write(15,*)'phys=',phys
      write(15,*)'fmcall=',fmcall
      write(15,*)'sampling_mode=',sampling_mode
      if (sampling_mode .ne. 0) then
         write(15,*)'xb_min=',xb_min
         write(15,*)'xb_max=',xb_max
         write(15,*)'minus_t_min=',mt_min
         write(15,*)'minus_t_max=',mt_max
         write(15,*)'phi_min_deg=',phi_min
         write(15,*)'phi_max_deg=',phi_max
      endif
      close(unit=15)
      close(unit=12)
      close(unit=14)
      close(unit=16)

 99   continue
 
      stop
      end

C======================================================================
      FUNCTION MYRAN(IDUM)
C----------------------------------------------------------------------
C-
C- Purpose and Methods : Flat Random Numbers (from Numerical Recipes:"RAN2").
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
      SAVE
C
C Function parameters
      INTEGER*8 IDUM
      REAL    MYRAN
C
C Local parameters
      INTEGER M, IA, IC
      REAL RM
      PARAMETER (M=714025, IA=1366, IC=150889, RM=1.4005112E-6)
C
C Local variables
      INTEGER J, IY, IFF, IR(97)
      DATA IFF /0/
C
C
      IF(IDUM.LT.0 .OR. IFF.EQ.0) THEN
        IFF = 1
        IDUM = MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM = MOD(IA*IDUM+IC,M)
          IR(J) = IDUM
 11     CONTINUE
        IDUM = MOD(IA*IDUM+IC,M)
        IY = IDUM
      ENDIF
C
      J = 1 + (97*IY)/M
      IF(J.GT.97 .OR. J.LT.1) STOP
      IY   = IR(J)
      MYRAN = IY*RM
      IDUM = MOD(IA*IDUM+IC,M)
      IR(J)= IDUM
      
C
      RETURN
      END

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine missm(ebeam,es,ep,th0,ephot,cstk,phik,mpi,ppx,
     * ppy,ppz,eprot,ppix,ppiy,ppiz,epi,csthcm,phicm,wreal,mm2)


c     Input:
c         ebeam = incident electron beam energy
c         es  = incident electron energy at interaction point
c         ep = scattered electron energy
c         th0 = electron scattering angle
c         ephot = energy of radiated photon
c         cstk = cosine of the photon angle (relative to the q vector)
c         phik = azimuthal angle of photon
c         csthcm = cosine of hadronic decay angle in the hadronic frame
c         phicm = phi angle in the hadronic frame
c     Output:
c         ppx, ppy, ppz = proton momentum components
c         eprot = proton energy
c         ppix,ppiy,ppiz = pion momentum components
c         epi = pion energy
c         wreal = true hadronic invariant mass
c         mm2 = experimental missing mass


c     Choose the hadronic decay angles randomly and calculate the
c     missing mass, the proton momenta and pion momenta.

      implicit none
      common /random/idum

      real*8 es,ep,ps,pp,ephot,th0,csthe,snthe,qx,qz


      real beta
      real cstk
      real cstq
      real csthcm
      real cxx,cxy,cxz,cyx,cyy,cyz,czx,czy,czz
      real csphk
      real csphi
      real ebeam
      real ephotx,ephoty,ephotz,eph_dot_q
      real epcm
      real epi
      real epicm
      real eprot
      real ewreal
      real gamma
      real mel
      real mp
      real mpi
      real mm2
      real nu
      real pfac
      real phicm
      real phik
      real pi
      real ppix,ppiy,ppiz
      real ppx,ppy,ppz
      real ppiwx,ppiwy,ppiwz
      real ppwx,ppwy,ppwz
      real pstar
      real pwrx,pwry,pwrz,pwr
      real q2,qvec
      real q_dot_pp
      real snphk
      real snphi
      real snthcm
      real sntk
      real sntq
      real w2
      real wmin
      real wreal

      integer*8 idum
c

      mp	= .938
      mel	= 0.511E-3
      pi	= 3.14159
      wmin	= mp+mpi
      csthe	= cos(th0)
      snthe	= sin(th0)
      nu	= es-ep
      ps	= abs(es**2-mel**2)
      pp	= abs(ep**2-mel**2)
      ps	= sqrt(ps)
      pp	= sqrt(pp)
      q2	= 2.*es*ep-2.*ps*pp*csthe-2.*mel**2
      w2	= mp**2-q2+2.*mp*nu
      qx	= -pp*snthe
      qz	= ps-pp*csthe
      qvec	= sqrt(qx**2+qz**2)
      if (abs(cstk) .gt. 1.)then
         write(6,*)' missm-1: cstk =',cstk
         cstk	= cstk/abs(cstk)
      endif
      sntk	= sqrt(1.-cstk**2)
      csphk	= cos(phik)
      snphk	= sin(phik)

      cstq	= qz/qvec
      sntq	= sqrt(1.-cstq**2)
      
c      ephotx=ephot*(sntk*csphk*cstq-cstk*sntq)
c      ephoty=ephot*sntk*snphk
c      ephotz=ephot*(cstk*cstq+sntk*csphk*sntq)

      ephot	= 0.
      ephotx	= 0.
      ephoty	= 0.
      ephotz	= 0.
      
c     Calculate the dot product of the photon vector and the q-vector
c     eph_dot_q=ephotx*qx+ephotz*qz

      eph_dot_q	= 0.

c     Calculate the mass of the actual hadronic system for the two
c     photon directions.

      wreal	= (w2-2.*ephot*(nu+mp)+2.*eph_dot_q)
      if (wreal .le. wmin**2) go to 30	!Hadronic mass below pion threshold
      wreal	= sqrt(wreal)
      
c     Calculate the energy of the actual hadronic system, with radiation

      ewreal	= nu+mp-ephot

c     Calculate the momentum components of the nucleon and pion in the lab

      snthcm	= (1.-csthcm**2)
      if (snthcm .le. 0.)then
         write(6,*)' missm-2: snthcm =',snthcm
         snthcm = 0.0000001
      endif
      snthcm	= sqrt(snthcm)

c     Calculate cos and sin of phicm

      csphi	= cos(phicm*pi/180.)
      snphi	= sin(phicm*pi/180.)

c     Calculate laboratory components of the resonance momentum vector

      pwrz	= qz-ephotz
      pwrx	= qx-ephotx
      pwry	= -ephoty
      pwr	= sqrt(pwrz**2+pwrx**2+pwry**2)
      
c     Calculate the relativistic factors, gamma and beta, for the resonance

      beta	= pwr/ewreal
      gamma	= ewreal/wreal

c     Define angle cosines for the laboratory system with the resonance
c     as the z axis.  Choose the y axis of the resonance frame perpendicular
c     to the laboratory x axis.

      pfac	= sqrt(pwry**2+pwrz**2)
      cxx	= pfac/pwr
      cxy	= -pwrx*pwry/pfac/pwr
      cxz	= -pwrx*pwrz/pfac/pwr
      cyx	= 0
      cyy	= pwrz/pfac
      cyz	= -pwry/pfac
      czx	= pwrx/pwr
      czy	= pwry/pwr
      czz	= pwrz/pwr

c     Calculate the momentum of the pion and proton in the resonance frame

      pstar	= ((wreal**2-mp**2-mpi**2)**2/4.-(mp*mpi)**2 )/wreal**2
      if (pstar .le. 0.) then
        write(6,*)'missm-3: pstar=',pstar
        pstar	= 0.000001
      endif
      pstar	= sqrt(pstar)
      
c     Calculate the energy of the proton and pion in the resonance center
c     of mass frame

      epcm  	= sqrt(pstar**2+mp**2)
      epicm 	= sqrt(pstar**2+mpi**2)

c     Calculate the pion momentum components and energy in the lab frame
c     where the z axis is the direction of the momentum of the resonance.
c     The center of mass angles are pion angles.

      ppiwx	= pstar*snthcm*csphi
      ppiwy	= pstar*snthcm*snphi
      ppiwz	= gamma*(pstar*csthcm+beta*epicm)
      epi	= gamma*(epicm+beta*pstar*csthcm)
      
c     Calculate momentum components of the pion in
c     the lab frame where the z axis is along the incident beam

      ppix	= ppiwx*cxx+ppiwy*cyx+ppiwz*czx
      ppiy	= ppiwx*cxy+ppiwy*cyy+ppiwz*czy
      ppiz	= ppiwx*cxz+ppiwy*cyz+ppiwz*czz

c     Calculate the proton energy and momentum components
c     in the lab resonance frame
c     The center of mass angles are pion angles

      ppwx	= -ppiwx
      ppwy	= -ppiwy

c     ppwz=-ppiwz+gamma*beta*(epicm+epcm)
c     epicm+epcm=wreal

      ppwz	= gamma*beta*wreal-ppiwz
      eprot	= gamma*wreal-epi

c     Rotate the lab momentum components of the proton into the frame
c     where the z axis is along the incident electron beam

      ppx	= ppwx*cxx+ppwy*cyx+ppwz*czx
      ppy	= ppwx*cxy+ppwy*cyy+ppwz*czy
      ppz	= ppwx*cxz+ppwy*cyz+ppwz*czz

c     Calculate the square of the missing mass, associated with the
c     momentum components of the charged hadron

      if (mpi .lt. .137) then
        q_dot_pp	= qx*ppx+qz*ppz
        mm2		= -q2+2*mp**2+2*mp*(nu-eprot)-2*nu*eprot+2*q_dot_pp
      else
        q_dot_pp	= qx*ppix+qz*ppiz
        mm2		= -q2+mp**2+mpi**2+2*mp*(nu-epi)-2*nu*epi+2*q_dot_pp
      endif

      return
      
 30   mm2=0.
 
      return
      end
C======================================================================
      function get_spin(iseed)
C----------------------------------------------------------------------
C-
C-   Purpose : Get spin (1 or -1)
C-
C-   Inputs  : random seed
C-
C-   Outputs : get_spin
C----------------------------------------------------------------------
      implicit none
      integer get_spin
      integer*8 iseed
      real    random, myran 

      random = myran(iseed)
      random = 0.5 - myran(iseed)
      get_spin = 1
      if(random.lt.0) get_spin = -1

      return
      end
C======================================================================
