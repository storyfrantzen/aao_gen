
      program aaorad_gen

c     This program makes an n-tuple that can be used with Paw to
c     make distributions of energies, angles, resonance
c     mass resulting from internal bremmstrahlung associated with pion
c     production on a proton. The exact integration formula of Mo and Tsai
c     is used.
c     The n-tuple contains the photon energy(EG), the true hadronic invariant
c     mass (W), the components of the proton momentum (PPX, PPY, PPZ),
c     the proton energy (EP), the pion momentum (PPIX, PPIY, PPIZ) and
c     pion energy (EPI), the  angles for the hadronic
c     decay in the hadronic frame (CSTCM, PHICM), the missing mass (MM),
c     and the photon angles relative to the q vector, (CSTHK, PHIK).
c
c     This program forces the monte carlo to concentrate on the regions
c     of photon emission along the directions of the incident and
c     scattered electrons.
c
c     The electrons are radiated as they pass through the target. Resolution
c     of detectors is not folded into the results.  If this is desired it
c     should be done with a second program that can operate on the n-tuple
c     and make a new version.

      implicit none
      
#include "bcs.inc"
#include "mc.inc"
#include "names.inc"

      COMMON/ALPHA/ ALPHA,PI,MP,MN,MPI,MEL,WG,EPIREA,TH_OPT,RES_OPT
      common /radcal/T0,es,ep,ps,pp,rs,rp,u0,pu,uu,cst0,snt0,csths,csthp
     * ,snths,snthp,pdotk,sdotk
      common /random/idum

      real lundmass(4)
      integer timevals(8)
      integer sseed(64)
      integer ssize
      integer env_status, env_len
      integer*8 job_seed, attempt_seed, seed_value
      character*32 env_value
      real*8 ek,Tk,delta
      real*8 alpha,pi,mp,mn,mpi,mel,wg,T0
      real*8 es,ep,ps,pp,rs,rp,u0,pu,uu,pdotk,sdotk
      real*8 cst0,snt0,csths,csthp,snths,snthp

      real csran,csrng,csrnge,csrngb,delphi
      real csthcm
      real csthcm_max
      real cstk
      real cstk1,cstk2
      real cstmp
      real csdotk,cpdotk,cqdotk
      real delinf
      real deltar
      real ek_max
      real ekmax
      real ekx,eky,ekz
      real epeps
      real epi
      real epmax
      real eprng
      real eprot
      real epw
      real ep_min,ep_max,ep_test
      real ep_sav
      real ep_out
      real*8 events
      real f
      real fkt
      real fmcall
      real g
      real jacob
      real kexp
      real kfac
      real mcfac
      real mm2
      real mm_exp,mm_cut
      real mpfac
      real mpi0
      real mpip
      real meta
      real mpi_s
      real myran
      real nu
      real phik
      real ppx,ppy,ppz
      real ppix,ppiy,ppiz
      real phicm
      real phicm_max
      real phir
      real q0
      real q2
      real q2_min,q2_max,q2max
      real qsq
      real qvecx
      real qvecz
      real ran
      real reg1,reg2,reg3,reg4
      real rn1,rn2
      real rtest
      real s
      real s1
      real s2
      real sigi
      real sigl
      real sigma
      real sigma0
      real signr
      real sigr,sig_ratio
      real sigr_observed
      real sigr_max
      real sigr1
      real*8 sig_tot,sig_sum
      real sigt
      real sigu
      real sigip, asym_p
      real sntk
      real sp
      real spence
      real stest
      real t_elapse
      real th0
      real theta
      real tk_max
      real tp
      real*8 tries
      real itime1, itime2
      real ts
      real uek
      real uq2, uq2_min,uq2_max,uq2rng
      real w2
      real wreal, w_min_input
      real wreal_h
      real ppx_h,ppy_h,ppz_h,eprot_h
      real q2_lep,xb_lep,w_lep
      real q2_hard,xb_hard,minus_t_hard
      real q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs
      real r_u,r_ep
      real vertex_fraction,incoming_loss_base
      real incoming_loss_accept_test
      real photon_side_base,photon_cosine_base,photon_phi_base
      real hadron_cosine_base,hadron_phi_base
      real outgoing_loss_base,outgoing_loss_accept_test
      real incoming_loss_fraction,outgoing_loss_fraction
      real x1
      real x2
      real targs,targp,xs,eloss,gxs,xtest,ebeam,t_targ,bfac,r_targ,temp
      real*8 sig_int
      real*8 analysis_volume
      real*8 analysis_mp
      real*8 survey_sig_internal_tot
      real*8 survey_internal_contribution
      real*8 survey_observed_contribution
      real hydrogen_rad, vertex_x, vertex_y, vertex_z, vz
      integer max_survey_bins
      parameter (max_survey_bins=64)
      real survey_q2_edges(max_survey_bins+1)
      real survey_xb_edges(max_survey_bins+1)
      real survey_t_edges(max_survey_bins+1)
      real survey_phi_edges(max_survey_bins+1)
      real survey_legacy_fraction,survey_balanced_fraction
      real proposal_selector,proposal_density_ratio
      real proposal_xb,proposal_minus_t,proposal_t_jacobian
      real proposal_t_at_cosine_zero
      real proposal_balanced_to_legacy
      real mode4_core_fraction,mode4_legacy_fraction
      real mode4_base_padding,mode4_core_volume
      real mode4_q2_low,mode4_q2_high
      real mode4_xb_low,mode4_xb_high
      real mode4_t_low,mode4_t_high
      real mode4_phi_low,mode4_phi_high
      real mode4_w_minimum,mode4_y_maximum
      real mode4_ru_low,mode4_ru_high
      real mode4_rep_low,mode4_rep_high
      real mode4_ugamma_low,mode4_ugamma_high
      real mode4_hcos_low,mode4_hcos_high
      real mode4_hphi_origin,mode4_hphi_relative_low
      real mode4_hphi_relative_high,mode4_hphi_relative
      real mode4_density_denominator
      real mode4_component_core_fraction
      real mode4_component_weight,mode4_trial_integrand

      integer dismc(6,100)
      integer intreg

      integer epirea, th_opt, res_opt
      integer i
      integer ir1
c      integer iext
      integer j
      integer jj
      integer mcall
      integer mcall_max
      integer sampling_mode,input_status
      integer survey_proposal_mode
      integer survey_nq2_bins,survey_nxb_bins
      integer survey_nt_bins,survey_nphi_bins
      integer proposal_component,proposal_q2_bin,proposal_xb_bin
      integer proposal_t_bin,proposal_phi_bin,proposal_mapping_status
      integer candidate_status,coordinate_status,final_valid
      integer mode4_apply_y_max,mode4_stratum_status
      integer mode4_operation,mode4_inside_core
      integer mode4_calibration_proposal
      integer stratum_flat,stratum_iq2,stratum_ixb,stratum_it
      integer stratum_iphi
      integer nprint,ntell,ntold
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer get_spin
      integer flag_ehel
      integer ehel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 idum
      integer*8 ntries
      integer*8 survey_ntrials,survey_seed,survey_replica
      integer*8 survey_rows,survey_final_valid
      integer*8 survey_legacy_trials,survey_balanced_trials
      integer*8 mode4_seed,mode4_replica
      integer*8 mode4_core_trials,mode4_legacy_trials
      integer*8 mode4_internal_rows,mode4_final_candidates
      integer*8 mode4_target_candidates
      integer*8 mode4_core_events,mode4_legacy_events
      integer*8 mode4_core_targets,mode4_legacy_targets
      integer*8 mode4_trial_limit,mode4_heartbeat_interval
      integer*8 mode4_last_heartbeat
      
      real cfac,asig

c     Parameters for the n-tuple, which is named func1 and contains
c     15 elements per event.

clcs      common /pawc/h(5000000)
      integer h,n,nevent,nmax,lrecl,istat,icycle
      parameter (n=32)
      real*4 ntp(n)
      
c     tag is the an array of names for the variables in the n-tuple.

c     character*1 ich1
      character*3 month
      character*2 day
      character*2 year
      character*5 tag(n)
      character*13 filerz
      character*13 file_out
      character*13 lund_out
      character*13 file_sum
      character*13 file_norm
      character*13 file_bos
      character*18 file_survey
      character*18 file_mode4
      character*32 file_mode4_heartbeat
      character*32 file_mode4_calibration
      character*8 recname
      character*28  ctime

      data tag /' ES  ',' EP  ','THETE','  W  ','WREAL',
     * ' PPX ',' PPY ',' PPZ ','EPROT','PPIX ',
     *  'PPIY ','PPIZ ',' EPI ','CSTCM','PHICM',' MM  ','  EG ',
     * 'CSTHK','PHIK ', ' QX  ',' QZ  ',' Q0  ','CSTHE',' EGX ',
     * ' EGY ',' EGZ ', ' VX  ', ' VY  ', ' VZ  ', ' Q2  ',
     * ' HEL ','ASYM '/


      DATA PI   /3.1415926/
      DATA MPIP /.1395/
      DATA MPI0 /.1349/
      DATA META /.547853/
      DATA MP   /.93828/
      DATA MN   /.93957/
      DATA MEL  /.511E-3/
      DATA analysis_mp /0.9382720813D0/

      data filerz   /'aao_rad.hbook'/
      data file_out /'aao_rad.out'/
      data lund_out /'aao_rad.lund'/
      data file_sum /'aao_rad.sum'/
      data file_norm /'aao_rad.norm'/
      data file_bos /'aao_rad.bos'/
      data file_survey /'aao_rad.survey.csv'/
      data file_mode4 /'aao_rad.mode4.csv'/
      data file_mode4_heartbeat /'aao_rad.mode4.heartbeat.csv'/
      data file_mode4_calibration /'aao_rad.mode4.calibration.csv'/
      data ctime    /'                            '/

      do j=1,6
        do i=1,100
         dismc(j,i)=0
        enddo
      enddo

c     set up parameters for breaking the monte-carlo integration region
c     over csthk into 5 parts:

      csrng=.04

c      Region sizes suggested for 4 GeV: reg1=.23, reg2=.14, reg3=.11, reg4=.10

 14   print *, 'Enter theory_opt,resonance_opt: (1=AO,4=MAID98,5=MAID2000)'
      read(5,*) th_opt
c
      print *, th_opt

      print *, 'Enter 1 for polarized electron, 0 for unpolarized electron'
      read(5,*) flag_ehel
c      
      print *, flag_ehel
      res_opt = 0
      
      print *, 'Input the sizes of the integration regions'
      read(5,*)reg1,reg2,reg3,reg4
      
      reg2 = reg1 + reg2
      reg3 = reg2 + reg3
      reg4 = reg3 + reg4
      
      if (reg4 .gt. .95) then
        write(6,*)' The sum of the region sizes must be less than .95'
        go to 14
      endif
      
      alpha = 1/137.

c     set up parameters for bos bank input to GSIM

      print *, ' Input 2 for two charged particles in the bos bank'
      print *, ' Input 4 to get the neutral hadron and photon as well'              
      read(5,*)npart
      print *, ' Input=',npart
      
      if (npart .ne. 4) npart=2

      q(1) 	= -1
      id(1)	= 3	!Geant ID, e-
      pdgid(1)  = 11	!PDG ID, e-
      
      if (npart .eq. 4)then
         q(4)	   = 0
         id(4)     = 1	!Geant ID, photon
         pdgid(4)  = 22	!PDG ID, photon
      endif
c   
c     Choose whether a neutral or charged pion is made in the reaction

 5    print *, ' Input epirea (1 for pi0, 2 for pi-, 3 for pi+, 5 for eta)'
      read(5,*)epirea
      print *, 'epiarea=', epirea
      
      print *, ' Input a limit on the error in (mm)**2'
      read(5,*)mm_cut
      print *, 'mm_cut=', mm_cut
      
      IF(epirea.eq.1)then
         MPI	= MPI0
         mm_exp = MPI0**2
         id(2)	= 14		!Geant ID, proton
         q(2) 	= 1
         pdgid(2)  = 2212		!PDG ID, proton
         if (npart .eq. 4)then
            id(3)     =  7              !Geant ID, pi-zero
            pdgid(3)  = 111		!PDG ID, pi-zero
            q(3)      = 0
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
      elseif(epirea.eq.3)then
         MPI	= MPIP
         mm_exp	= mp**2
         id(2)     = 8                !Geant ID, pi-plus
         pdgid(2)  = 211		!PDG ID, pi-plus
         q(2) 	   = 1
         if (npart .eq. 4)then
            id(3)	= 13		!Geant ID, neutron
            pdgid(3)    = 2112		!PDG ID, neutron
            q(3) 	= 0
         endif
      elseif(epirea.eq.5)then
         MPI	= META
         mm_exp = META**2
         id(2)	= 14		!Geant ID, proton
         q(2) 	= 1
         pdgid(2)  = 2212		!PDG ID, proton
         if (npart .eq. 4)then
            id(3)     = 17              !Geant ID, eta-zero
            pdgid(3)  = 221		!PDG ID, pi-zero
            q(3)      = 0
         endif
      else
         go to 5
      endif
      
      if (epirea .eq. 3)filerz(13:13)='p'

c     Set single precision version of pion mass
      mpi_s=mpi

c     Calculate the minimum hadronic mass for pion production:
      wg	= mp + mpi + .0005

      write(6,*)' Input the target thickness (cm)'
      read(5,*)t_targ            
      print *, 't_target=', t_targ

      write(6,*)' Input the target radius, (cm)'
      read(5,*)r_targ
      print *, 'r_target=', r_targ

      write(6,*)' Input the x-coordinate of the beam position, (cm)'
      read(5,*)vertex_x
      write(6,*)' Input the y-coordinate of the beam position, (cm)'
      read(5,*)vertex_y
      write(6,*)' Input the z-coordinate of the beam position, (cm)'
      read(5,*)vz
      print *, 'VX,VY,VZ=', vertex_x,vertex_y,vertex_z

      bfac		= 4./3.
      hydrogen_rad 	= 865 		! hydrogen radiation length (cm)
      t_targ 		= bfac * t_targ / hydrogen_rad 

      write(6,*)' Input the incident electron energy (GeV)'
      read(5,*)ebeam
      print *, 'E beam=', ebeam
      
c     calculate the incident momentum

      es	= ebeam
      ps	= sqrt(es**2-mel**2)
      rs	= ps/es
      
c     cut off q2 at the value for 90 degree elastic scattering

      s		= 0.5
      q2max	= 4.*ebeam**2*s/(1.+2.*ebeam*s/mp)
      
c     Choose two limits for Q**2

      write(6,*)'Input lower and upper limit for Q**2'
      read(5,*)q2_min,q2_max
      print *, ' Q2 min max=', q2_min,q2_max
      
      if (q2_max .gt. q2max) q2_max = q2max
      uq2_min	= 1/q2_max
      uq2_max	= 1/q2_min
      uq2rng	= uq2_max-uq2_min

c     Set the limits on the range of scattered electron energies

      write(6,*)'Input lower and upper limit for scat e- energy(GeV).'
      read(5,*)ep_min,ep_max
      print *, 'Eelectron min max', ep_min,ep_max
      
      epmax	= es - (wg**2 + q2_min - mp**2)/2./mp
      
      if (ep_max .lt. epmax)epmax=ep_max
      eprng	= epmax - ep_min

c     Choose a maximum value for the range of photon energies

      write(6,*)' Input minimum photon energy for integration'
      read(5,*)delta
      print *, ' Min photon energy=', delta

c     Select the number of events desired in the rz file

      write(6,*)' Input the desired number of events in the ntuple'
      read(5,*)nmax
      print *, 'Number of events in the ntuple=', nmax
      
      nprint	= max(1,nmax/25)

      write(6,*)' Input a multiplication factor for sigr_max'
      read(5,*)fmcall
      print *, ' Factor=', fmcall
      
      if (fmcall .eq. 0.)then
         write(6,*)' Input sigr_max'
         read(5,*)sigr_max
         print *, 'sigr_max=', sigr_max
      endif
c
      if (th_opt .gt. 10) then
         write(6,*)'W mininum input'
         read(5,*) w_min_input
         print *, 'W min=', w_min_input
      endif
c
c     Optional backward-compatible sampling trailer.  Legacy input files end
c     above and therefore retain sampling_mode=0.  Fixed-trial survey mode 1
c     requires an exact proposal count, explicit nonzero seed, and replica ID.
c     Mode 4 generates one final-LUND analysis stratum from an exact mixture
c     of a learned continuous core and the unrestricted legacy proposal.

      sampling_mode=0
      survey_ntrials=0
      survey_seed=0
      survey_replica=0
      survey_proposal_mode=0
      survey_legacy_fraction=1.
      survey_balanced_fraction=0.
      survey_nq2_bins=0
      survey_nxb_bins=0
      survey_nt_bins=0
      survey_nphi_bins=0
      mode4_seed=0
      mode4_replica=0
      mode4_core_fraction=0.
      mode4_legacy_fraction=0.
      mode4_base_padding=0.
      mode4_core_volume=0.
      mode4_apply_y_max=0
      mode4_y_maximum=1.
      mode4_operation=0
      mode4_calibration_proposal=0
      mode4_trial_limit=0
      mode4_heartbeat_interval=100000
      mode4_component_core_fraction=0.
      stratum_flat=-1
      stratum_iq2=-1
      stratum_ixb=-1
      stratum_it=-1
      stratum_iphi=-1
      read(5,*,iostat=input_status)sampling_mode
      if (input_status .lt. 0) sampling_mode=0
      if (input_status .gt. 0) then
         write(6,*)'Invalid optional sampling mode input'
         stop 2
      endif
      if (sampling_mode .ne. 0 .and. sampling_mode .ne. 1 .and.
     *    sampling_mode .ne. 4) then
         write(6,*)'sampling_mode must be 0 (legacy), 1 (survey),'
     *      //' or 4 (bin conditional)'
         stop 2
      endif
      if (sampling_mode .eq. 1) then
         read(5,*,iostat=input_status)survey_ntrials
         if (input_status .ne. 0 .or. survey_ntrials .le. 0) then
            write(6,*)'Survey mode requires a positive fixed trial count'
            stop 2
         endif
         read(5,*,iostat=input_status)survey_seed
         if (input_status .ne. 0 .or. survey_seed .eq. 0) then
            write(6,*)'Survey mode requires an explicit nonzero seed'
            stop 2
         endif
         read(5,*,iostat=input_status)survey_replica
         if (input_status .ne. 0 .or. survey_replica .lt. 0) then
            write(6,*)'Survey mode requires a nonnegative replica ID'
            stop 2
         endif
         read(5,*,iostat=input_status)survey_proposal_mode
         if (input_status .lt. 0) survey_proposal_mode=0
         if (input_status .gt. 0) then
            write(6,*)'Invalid optional survey proposal mode'
            stop 2
         endif
         if (survey_proposal_mode .ne. 0 .and.
     *       survey_proposal_mode .ne. 1) then
            write(6,*)'Survey proposal mode must be 0 or 1'
            stop 2
         endif
         if (survey_proposal_mode .eq. 1) then
            read(5,*,iostat=input_status)survey_legacy_fraction
            if (input_status .ne. 0 .or.
     *          survey_legacy_fraction .le. 0. .or.
     *          survey_legacy_fraction .ge. 1.) then
               write(6,*)'Balanced survey legacy fraction must be in (0,1)'
               stop 2
            endif
            survey_balanced_fraction=1.-survey_legacy_fraction
            read(5,*,iostat=input_status)survey_nq2_bins
            if (input_status .ne. 0 .or. survey_nq2_bins .le. 0 .or.
     *          survey_nq2_bins .gt. max_survey_bins) then
               write(6,*)'Invalid balanced survey Q2 bin count'
               stop 2
            endif
            read(5,*,iostat=input_status)
     *         (survey_q2_edges(i),i=1,survey_nq2_bins+1)
            if (input_status .ne. 0) then
               write(6,*)'Cannot read balanced survey Q2 edges'
               stop 2
            endif
            read(5,*,iostat=input_status)survey_nxb_bins
            if (input_status .ne. 0 .or. survey_nxb_bins .le. 0 .or.
     *          survey_nxb_bins .gt. max_survey_bins) then
               write(6,*)'Invalid balanced survey xB bin count'
               stop 2
            endif
            read(5,*,iostat=input_status)
     *         (survey_xb_edges(i),i=1,survey_nxb_bins+1)
            if (input_status .ne. 0) then
               write(6,*)'Cannot read balanced survey xB edges'
               stop 2
            endif
            read(5,*,iostat=input_status)survey_nt_bins
            if (input_status .ne. 0 .or. survey_nt_bins .le. 0 .or.
     *          survey_nt_bins .gt. max_survey_bins) then
               write(6,*)'Invalid balanced survey -t bin count'
               stop 2
            endif
            read(5,*,iostat=input_status)
     *         (survey_t_edges(i),i=1,survey_nt_bins+1)
            if (input_status .ne. 0) then
               write(6,*)'Cannot read balanced survey -t edges'
               stop 2
            endif
            read(5,*,iostat=input_status)survey_nphi_bins
            if (input_status .ne. 0 .or.
     *          survey_nphi_bins .le. 0 .or.
     *          survey_nphi_bins .gt. max_survey_bins) then
               write(6,*)'Invalid balanced survey phi bin count'
               stop 2
            endif
            read(5,*,iostat=input_status)
     *         (survey_phi_edges(i),i=1,survey_nphi_bins+1)
            if (input_status .ne. 0) then
               write(6,*)'Cannot read balanced survey phi edges'
               stop 2
            endif
            do i=1,survey_nq2_bins
               if (survey_q2_edges(i+1) .le.
     *             survey_q2_edges(i)) then
                  write(6,*)'Balanced survey Q2 edges not increasing'
                  stop 2
               endif
            enddo
            do i=1,survey_nxb_bins
               if (survey_xb_edges(i+1) .le.
     *             survey_xb_edges(i)) then
                  write(6,*)'Balanced survey xB edges not increasing'
                  stop 2
               endif
            enddo
            do i=1,survey_nt_bins
               if (survey_t_edges(i+1) .le. survey_t_edges(i)) then
                  write(6,*)'Balanced survey -t edges not increasing'
                  stop 2
               endif
            enddo
            do i=1,survey_nphi_bins
               if (survey_phi_edges(i+1) .le.
     *             survey_phi_edges(i)) then
                  write(6,*)'Balanced survey phi edges not increasing'
                  stop 2
               endif
            enddo
            if (survey_q2_edges(1) .lt. q2_min .or.
     *          survey_q2_edges(survey_nq2_bins+1) .gt. q2_max) then
               write(6,*)'Balanced survey Q2 edges exceed input range'
               stop 2
            endif
            if (survey_xb_edges(1) .le. 0.) then
               write(6,*)'Balanced survey xB edges must be positive'
               stop 2
            endif
            if (survey_t_edges(1) .lt. 0.) then
               write(6,*)'Balanced survey -t edges must be nonnegative'
               stop 2
            endif
            if (abs(survey_phi_edges(1)) .gt. 1.e-5 .or.
     *          abs(survey_phi_edges(survey_nphi_bins+1)-360.)
     *          .gt. 1.e-4) then
               write(6,*)'Balanced survey phi edges must span 0 to 360'
               stop 2
            endif
         endif
         if (epirea .ne. 1 .or. npart .ne. 4) then
            write(6,*)'Survey mode currently requires pi0 and npart=4'
            stop 2
         endif
         write(6,*)'Fixed-trial radiative survey mode:',survey_ntrials
         write(6,*)'Explicit survey seed:',survey_seed
         write(6,*)'Survey replica ID:',survey_replica
         write(6,*)'Survey proposal mode:',survey_proposal_mode
         write(6,*)'Survey legacy fraction:',survey_legacy_fraction
      endif
      if (sampling_mode .eq. 4) then
         read(5,*,iostat=input_status)mode4_seed
         if (input_status .ne. 0 .or. mode4_seed .eq. 0) then
            write(6,*)'Mode 4 requires an explicit nonzero seed'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_replica
         if (input_status .ne. 0 .or. mode4_replica .lt. 0) then
            write(6,*)'Mode 4 requires a nonnegative replica ID'
            stop 2
         endif
         read(5,*,iostat=input_status)stratum_flat,stratum_iq2,
     *      stratum_ixb,stratum_it,stratum_iphi
         if (input_status .ne. 0 .or. stratum_flat .lt. 0 .or.
     *       stratum_iq2 .lt. 0 .or. stratum_ixb .lt. 0 .or.
     *       stratum_it .lt. 0 .or. stratum_iphi .lt. 0) then
            write(6,*)'Mode 4 requires nonnegative stratum indices'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_q2_low,mode4_q2_high
         if (input_status .ne. 0) stop 2
         read(5,*,iostat=input_status)mode4_xb_low,mode4_xb_high
         if (input_status .ne. 0) stop 2
         read(5,*,iostat=input_status)mode4_t_low,mode4_t_high
         if (input_status .ne. 0) stop 2
         read(5,*,iostat=input_status)mode4_phi_low,mode4_phi_high
         if (input_status .ne. 0) stop 2
         if (mode4_q2_high .le. mode4_q2_low .or.
     *       mode4_xb_high .le. mode4_xb_low .or.
     *       mode4_t_high .le. mode4_t_low .or.
     *       mode4_phi_high .le. mode4_phi_low) then
            write(6,*)'Mode 4 stratum bounds must increase'
            stop 2
         endif
         if (mode4_q2_low .lt. q2_min .or.
     *       mode4_q2_high .gt. q2_max .or.
     *       mode4_xb_low .le. 0. .or. mode4_t_low .lt. 0. .or.
     *       mode4_phi_low .lt. 0. .or.
     *       mode4_phi_high .gt. 360.) then
            write(6,*)'Mode 4 stratum bounds exceed generator domain'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_w_minimum
         if (input_status .ne. 0 .or. mode4_w_minimum .le. 0.) then
            write(6,*)'Mode 4 requires a positive W minimum'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_apply_y_max,
     *      mode4_y_maximum
         if (input_status .ne. 0 .or.
     *       (mode4_apply_y_max .ne. 0 .and.
     *        mode4_apply_y_max .ne. 1) .or.
     *       mode4_y_maximum .le. 0. .or.
     *       mode4_y_maximum .gt. 1.) then
            write(6,*)'Invalid optional mode-4 y selection'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_core_fraction
         if (input_status .ne. 0 .or. mode4_core_fraction .le. 0.
     *       .or. mode4_core_fraction .ge. 1.) then
            write(6,*)'Mode 4 core fraction must lie in (0,1)'
            stop 2
         endif
         mode4_legacy_fraction=1.-mode4_core_fraction
         read(5,*,iostat=input_status)mode4_base_padding
         if (input_status .ne. 0 .or. mode4_base_padding .lt. 0.) then
            write(6,*)'Mode 4 base padding must be nonnegative'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_core_volume
         if (input_status .ne. 0 .or. mode4_core_volume .le. 0.
     *       .or. mode4_core_volume .gt. 1.) then
            write(6,*)'Mode 4 core volume must lie in (0,1]'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_ru_low,mode4_ru_high
         if (input_status .ne. 0) stop 2
         read(5,*,iostat=input_status)mode4_rep_low,mode4_rep_high
         if (input_status .ne. 0) stop 2
         read(5,*,iostat=input_status)mode4_ugamma_low,
     *      mode4_ugamma_high
         if (input_status .ne. 0) stop 2
         read(5,*,iostat=input_status)mode4_hcos_low,
     *      mode4_hcos_high
         if (input_status .ne. 0) stop 2
         if (mode4_ru_low .lt. 0. .or. mode4_ru_high .gt. 1. .or.
     *       mode4_ru_high .le. mode4_ru_low .or.
     *       mode4_rep_low .lt. 0. .or. mode4_rep_high .gt. 1. .or.
     *       mode4_rep_high .le. mode4_rep_low .or.
     *       mode4_ugamma_low .lt. 0. .or.
     *       mode4_ugamma_high .gt. 1. .or.
     *       mode4_ugamma_high .le. mode4_ugamma_low .or.
     *       mode4_hcos_low .lt. 0. .or. mode4_hcos_high .gt. 1.
     *       .or. mode4_hcos_high .le. mode4_hcos_low) then
            write(6,*)'Invalid nonperiodic mode-4 guard bounds'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_hphi_origin,
     *      mode4_hphi_relative_low,mode4_hphi_relative_high
         if (input_status .ne. 0 .or. mode4_hphi_origin .lt. 0.
     *       .or. mode4_hphi_origin .ge. 1. .or.
     *       mode4_hphi_relative_low .lt. -.5 .or.
     *       mode4_hphi_relative_high .gt. .5 .or.
     *       mode4_hphi_relative_high .le.
     *          mode4_hphi_relative_low) then
            write(6,*)'Invalid periodic mode-4 guard bounds'
            stop 2
         endif
         mode4_density_denominator=
     *      (mode4_ru_high-mode4_ru_low)*
     *      (mode4_rep_high-mode4_rep_low)*
     *      (mode4_ugamma_high-mode4_ugamma_low)*
     *      (mode4_hcos_high-mode4_hcos_low)*
     *      (mode4_hphi_relative_high-mode4_hphi_relative_low)
         if (abs(mode4_density_denominator-mode4_core_volume)
     *       .gt. 2.e-5*mode4_core_volume) then
            write(6,*)'Mode 4 core volume disagrees with guard bounds'
            stop 2
         endif
         read(5,*,iostat=input_status)mode4_operation
         if (input_status .lt. 0) then
            mode4_operation=0
            mode4_trial_limit=0
            mode4_heartbeat_interval=100000
            mode4_component_core_fraction=mode4_core_fraction
         elseif (input_status .ne. 0 .or.
     *           (mode4_operation .ne. 0 .and.
     *            mode4_operation .ne. 1)) then
            write(6,*)'Mode 4 operation must be 0 or 1'
            stop 2
         else
            read(5,*,iostat=input_status)mode4_trial_limit
            if (input_status .ne. 0 .or. mode4_trial_limit .lt. 0) then
               write(6,*)'Invalid mode-4 fixed trial limit'
               stop 2
            endif
            read(5,*,iostat=input_status)mode4_heartbeat_interval
            if (input_status .ne. 0 .or.
     *          mode4_heartbeat_interval .le. 0) then
               write(6,*)'Mode-4 heartbeat interval must be positive'
               stop 2
            endif
            read(5,*,iostat=input_status)
     *         mode4_component_core_fraction
            if (input_status .ne. 0 .or.
     *          mode4_component_core_fraction .le. 0. .or.
     *          mode4_component_core_fraction .ge. 1.) then
               write(6,*)'Mode-4 component fraction must be in (0,1)'
               stop 2
            endif
            read(5,*,iostat=input_status)mode4_calibration_proposal
            if (input_status .lt. 0) mode4_calibration_proposal=0
            if (input_status .gt. 0 .or.
     *          (mode4_calibration_proposal .ne. 0 .and.
     *           mode4_calibration_proposal .ne. 1)) then
               write(6,*)'Mode-4 calibration proposal must be 0 or 1'
               stop 2
            endif
         endif
         if (mode4_operation .eq. 0) then
            mode4_trial_limit=0
            mode4_component_core_fraction=mode4_core_fraction
            mode4_calibration_proposal=0
         elseif (mode4_trial_limit .le. 0) then
            write(6,*)'Mode-4 calibration requires positive trials'
            stop 2
         endif
         if (mode4_operation .eq. 1 .and.
     *       mode4_calibration_proposal .eq. 1 .and.
     *       mode4_core_volume .ge. 1.-1.e-7) then
            write(6,*)'Mode-4 guard complement has zero volume'
            stop 2
         endif
         if (fmcall .ne. 0. .or. sigr_max .le. 0.) then
            write(6,*)'Mode 4 requires fmcall=0 and positive sigr_max'
            stop 2
         endif
         if (epirea .ne. 1 .or. npart .ne. 4) then
            write(6,*)'Mode 4 currently requires pi0 and npart=4'
            stop 2
         endif
         write(6,*)'Radiative bin-conditional mode 4'
         write(6,*)'Mode 4 seed/replica:',mode4_seed,mode4_replica
         write(6,*)'Mode 4 stratum:',stratum_flat,stratum_iq2,
     *      stratum_ixb,stratum_it,stratum_iphi
         write(6,*)'Mode 4 core/tail fractions:',
     *      mode4_core_fraction,mode4_legacy_fraction
         write(6,*)'Mode 4 core volume:',mode4_core_volume
         write(6,*)'Mode 4 operation/trial limit:',
     *      mode4_operation,mode4_trial_limit
         write(6,*)'Mode 4 heartbeat/component fraction:',
     *      mode4_heartbeat_interval,mode4_component_core_fraction
         write(6,*)'Mode 4 calibration proposal:',
     *      mode4_calibration_proposal
      endif
c
 1    mcall_max	= 0
      ntold	= 0
      events	= 0
      survey_rows=0
      survey_final_valid=0
      survey_sig_internal_tot=0.
      survey_legacy_trials=0
      survey_balanced_trials=0
      mode4_core_trials=0
      mode4_legacy_trials=0
      mode4_internal_rows=0
      mode4_final_candidates=0
      mode4_target_candidates=0
      mode4_core_events=0
      mode4_legacy_events=0
      mode4_core_targets=0
      mode4_legacy_targets=0
      mode4_last_heartbeat=-1

c     Use the internal clock to initialize the random number generator

      call timex(itime1)
      call getunixtime(idum)
      call getasciitime(idum,ctime)

      idum	= -idum
      month 	= ctime(5:7)
      day   	= ctime(9:10)
      year  	= ctime(23:24)
      
      if (day(1:1).eq. ' ')then
        ir1=48
        day(1:1)=char(ir1)
      endif

      call date_and_time(values=timevals)
      call random_seed(size=ssize)
      write(6,*)'seed:',ssize
      job_seed = 0
      attempt_seed = 0
      call get_environment_variable('SWIF_JOB_ID', env_value,
     *     length=env_len, status=env_status)
      if (env_status .eq. 0 .and. env_len .gt. 0) then
         read(env_value(1:env_len),*,err=431) job_seed
      endif
 431  continue
      call get_environment_variable('SWIF_JOB_ATTEMPT_ID', env_value,
     *     length=env_len, status=env_status)
      if (env_status .eq. 0 .and. env_len .gt. 0) then
         read(env_value(1:env_len),*,err=432) attempt_seed
      endif
 432  continue
      do i=1,ssize
         if (sampling_mode .eq. 1) then
            seed_value = abs(survey_seed) + 7919_8*i
         elseif (sampling_mode .eq. 4) then
            seed_value = abs(mode4_seed) + 7919_8*i
         else
            seed_value = timevals(8) + 1000_8*timevals(7)
     *         + 100000_8*timevals(6) + 10000000_8*timevals(5)
     *         + 104729_8*job_seed + 1009_8*attempt_seed
     *         + 7919_8*i
         endif
         sseed(i) = int(mod(abs(seed_value), 2147483647_8))
         if (sseed(i) .eq. 0) sseed(i) = i + 1
      enddo
      write(6,*)'seed mix:',job_seed,attempt_seed,sseed(1)
      call random_seed(put=sseed)

c      call random_number(rnd)
c      write(6,*)'seed:',rnd
c      call random_number(rnd)
c      write(6,*)'seed:',rnd
c      call random_number(rnd)
c      write(6,*)'seed:',rnd
c      write(6,*)'seed:',idum,' from start time ',ctime


      nevent	= 0
      t_elapse	= 0.
      itime2	= itime1
      ntries	= 0
      sig_int	= 0.
      sig_tot	= 0.

c     Name the output rz file according to meson type and beam energy.
c     filerz=aaoradgen-pi0-1.6-0811.rz.0, for example.

      open(unit=12,file=file_out)
      open(unit=13,file=lund_out)
      if (sampling_mode .eq. 1) then
         open(unit=16,file=file_survey,status='replace')
         write(16,'(A)')'# schema=aao-rad-survey-v2'
         write(16,'(A)')'replica,trial,final_valid,candidate_status,'
     *   //'intreg,helicity,proposal_component,proposal_q2_bin,'
     *   //'proposal_xb_bin,proposal_t_bin,proposal_phi_bin,'
     *   //'proposal_density_ratio,r_u,r_ep,u_gamma,channel_selector,'
     *   //'vertex_fraction,incoming_loss_base,'
     *   //'incoming_loss_accept_test,photon_side_base,'
     *   //'photon_cosine_base,photon_phi_base,'
     *   //'hadron_cosine_base,hadron_phi_base,'
     *   //'outgoing_loss_base,outgoing_loss_accept_test,'
     *   //'cos_theta_gamma,'
     *   //'phi_gamma_deg,cos_theta_cm,phi_cm_deg,energy_in_vertex,'
     *   //'energy_e_pre_external,energy_e_lund,energy_gamma,'
     *   //'final_e_px,final_e_py,final_e_pz,final_e_energy,'
     *   //'final_p_px,final_p_py,final_p_pz,final_p_energy,'
     *   //'q2_leptonic,xb_leptonic,w_leptonic,q2_hard,xb_hard,'
     *   //'w_hard,minus_t_hard,q2_observed,xb_observed,'
     *   //'minus_t_observed,phi_observed_deg,w_observed,y_observed,'
     *   //'integrand_internal,integrand_observed,'
     *   //'trial_xsec_internal_microbarn,'
     *   //'trial_xsec_observed_microbarn,rotation_deg,'
     *   //'incoming_loss_fraction,outgoing_loss_fraction'
      endif
      if (sampling_mode .eq. 4 .and. mode4_operation .eq. 0) then
         open(unit=17,file=file_mode4,status='replace')
         write(17,'(A)')'# schema=aao-rad-mode4-events-v1'
         write(17,'(A)')'event,proposal_component,'
     *      //'proposal_density_ratio,r_u,r_ep,u_gamma,'
     *      //'hadron_cosine_base,hadron_phi_base,'
     *      //'q2_observed,xb_observed,'
     *      //'minus_t_observed,phi_observed_deg,w_observed,'
     *      //'y_observed,integrand_corrected'
         flush(17)
      endif
      if (sampling_mode .eq. 4 .and. mode4_operation .eq. 1) then
         open(unit=19,file=file_mode4_calibration,status='replace')
         write(19,'(A)')'# schema=aao-rad-mode4-calibration-events-v2'
         write(19,'(A)')'trial,proposal_component,inside_core,'
     *      //'proposal_density_ratio,integrand_corrected,'
     *      //'component_importance_weight,r_u,r_ep,u_gamma,'
     *      //'hadron_cosine_base,hadron_phi_base,q2_observed,'
     *      //'xb_observed,minus_t_observed,phi_observed_deg,'
     *      //'w_observed,y_observed'
         flush(19)
      endif
      if (sampling_mode .eq. 4) then
         open(unit=18,file=file_mode4_heartbeat,status='replace')
         write(18,'(A)')'# schema=aao-rad-mode4-heartbeat-v2'
         write(18,'(A)')'proposals,events,core_trials,noncore_trials,'
     *      //'internal_rows,final_candidates,target_candidates,'
     *      //'core_targets,noncore_targets,core_events,noncore_events,'
     *      //'mcall_max'
         write(18,9010)ntries,nevent,mode4_core_trials,
     *      mode4_legacy_trials,mode4_internal_rows,
     *      mode4_final_candidates,mode4_target_candidates,
     *      mode4_core_targets,mode4_legacy_targets,
     *      mode4_core_events,mode4_legacy_events,mcall_max
         flush(18)
         mode4_last_heartbeat=ntries
      endif

c     Initialize BOS
     
      bosout 	= file_bos
      recname 	= 'MCEVENT'
     
Cvpk      call bos_ini(recname)
  
c     set up the ntuple file
      print *,'     set up the ntuple file'

      lrecl	= 1024
      
Clcs      call hlimit(5000000)
Clcs      call hropen(1,'aaoradgen',filerz,'n',lrecl,istat)
Clcs      call hbookn(10,'func1',n,'aaoradgen',1000,tag)

      open(unit=12,file=file_out)

      write(12,*)' AO Calculation of Single Pion Production'
      write(12,*)' Starting time:', ctime
      write(12,*)' Epirea (1 for pi0, 2 for pi-, 3 for pi+, 5 for eta) =',epirea
      write(12,*)' Target thickness =',t_targ*3./4.,' (r.l.)'
      write(12,*)' Incident electron energy =',ebeam,' GeV'

      write(12,*)' Electron Q**2 limits:',q2_min,q2_max
      write(12,*)' Lower and upper limit for scattered electron',
     * ' energy(GeV):',ep_min,epmax
      write(12,*)' Minimum photon energy for integration (delta):',delta

c     Use a new variable in place of ek. Let uek=exp(-kek*ek)
c     ek=-(1/kexp)alog(uek).  The factor of 5. was chosen empirically for kexp
c     by looking at the ek spectrum for E0=1.6 GeV.
c     Let uek range from 0 to 1. Then ek will range from 0 and infinity.
c     This requires a jacobian. Jacobian=1./(kexp*uek)=(1./kexp)exp(kexp*ek)

      kexp	= 5.
      analysis_volume=(4.*pi)**2*(2.*pi*uq2rng*eprng)

      if (sampling_mode .eq. 1)then
         sigr_max=0.
         go to 20
      endif
      if (sampling_mode .eq. 4)go to 20
      if (fmcall .eq. 0.)then
         write(6,*)' sigr_max from input data =',sigr_max
         write(12,*)' sigr_max from input data =',sigr_max
         go to 20
      endif

c     Do a preliminary calculation to estimate the maximum value
c     of the integrand

c     calculate the energy and momentum of the scattered electron,
c     and calculate Q**2 at the delta mass, 1.232 GeV

 10   q2	= q2_min
      q0	= (1.232**2 - mp**2 + q2)/2./mp
      if (th_opt.gt.10) then
        q0	= (w_min_input**2 - mp**2 + q2)/2./mp   
      endif
      ep	= es-q0
      pp	= sqrt(ep**2 - mel**2)
      rp	= pp/ep
      s		= q2/4/es/ep
      th0	= 2.*asin(sqrt(s))
      theta	= th0*180./pi
      T0	= th0
      snt0	= sin(th0)
      cst0	= cos(th0)

c     calculate kinematic quantities needed for the Mo and Tsai calculation

      u0	= es - ep + mp
      pu	= sqrt(ps**2 + pp**2 - 2*ps*pp*cst0)
      uu	= u0**2 - pu**2
      csths	= (ps - pp*cst0)/pu
      csthp	= (ps*cst0 - pp)/pu
      snths	= sqrt(1. - csths**2)
      snthp	= sqrt(1. - csthp**2)
      ts	= acos(csths)
      tp	= acos(csthp)
      qsq	= q2
      sp	= es*ep - ps*pp*cst0

      sigr_max	= 0.
      cstk1	= (es - ep*cst0)/(sqrt(es**2 + ep**2 - 2*es*ep*cst0))
      cstk2	= (es*cst0 - ep)/(sqrt(es**2 + ep**2 - 2*es*ep*cst0))
      csrnge	= csrng
      
      if ((1.-cstk1) .lt. csrnge)       csrnge = 1.-cstk1
      if ((cstk1-cstk2) .lt. 2.*csrnge) csrnge = 0.5*(cstk1-cstk2)
      
      csrngb	= csrng/40.
      
      if (csrngb .gt. csrnge/5.) csrngb = csrnge/5.
      
      delphi	= pi/9.
      phik	= (myran()-0.5)*delphi
      mpfac	= delphi/2./pi

      ek	= delta
      
      jacob	= exp(kexp*ek)/kexp*(q2**2/(2*es*ep))
c
      ehel = 0   
      if(flag_ehel.eq.1) ehel = 1    
c    
      do i=1,10000
        csthcm	= 2.*(myran()-0.5)
        phicm	= 360.*myran()

        csran	= myran()
        
        if (csran .gt. .3) then
            cstk = 2.*csrngb*(myran()-0.5) + cstk1
        else
            cstk = 2.*csrngb*(myran()-0.5) + cstk2
        endif
        
        mcfac	= csrngb /reg1
        phik	= (myran()-0.5)*delphi
        mpfac	= delphi/2./pi

        Tk	= acos(cstk)
        sntk	= sin(Tk)
        sdotk	= es*ek - ps*ek*cstk*csths - ps*ek*sntk*snths*cos(phik)
        pdotk	= ep*ek - pp*ek*cstk*csthp - pp*ek*sntk*snthp*cos(phik)

        sigr	= sigma(ek,Tk,csthcm,phicm,ehel)
        sigr	= sigr*mcfac*mpfac
        sigr	= sigr*jacob      
        if (sigr .gt. sigr_max) then
           sigr_max	= sigr
           ek_max	= ek
           tk_max	= Tk
           csthcm_max	= csthcm
           phicm_max	= phicm
        endif
        
      enddo

      write(6,*)'sigr_max,ek_max,tk_max,csthcm_max,phicm_max',
     * sigr_max,ek_max,tk_max,csthcm_max,phicm_max
      write(12,*)'sigr_max,ek_max,tk_max,csthcm_max,phicm_max',
     * sigr_max,ek_max,tk_max,csthcm_max,phicm_max
     
      sigr_max	= sigr_max*fmcall
      
 25   write(6,*) ' sigr_max changed to ',sigr_max
      write(12,*)' sigr_max changed to ',sigr_max

c   %%%%%%%%%%%%%%%%%%% Main Calculation  %%%%%%%%%%%%%%%%%%%%%%%
c     Use a Monte-Carlo to calculate a distribution of nmax events
c     distributed according to the Mo-Tsai integrand.
      ehel = 0      
 20   continue
      if (sampling_mode .eq. 1 .and.
     *    ntries .ge. survey_ntrials) go to 50
      if (sampling_mode .eq. 4 .and. ntries .gt. 0 .and.
     *    ntries .ne. mode4_last_heartbeat .and.
     *    mod(ntries,mode4_heartbeat_interval) .eq. 0) then
         write(18,9010)ntries,nevent,mode4_core_trials,
     *      mode4_legacy_trials,mode4_internal_rows,
     *      mode4_final_candidates,mode4_target_candidates,
     *      mode4_core_targets,mode4_legacy_targets,
     *      mode4_core_events,mode4_legacy_events,mcall_max
         flush(18)
         flush(12)
         flush(13)
         if (mode4_operation .eq. 0)flush(17)
         if (mode4_operation .eq. 1)flush(19)
         mode4_last_heartbeat=ntries
      endif
      if (sampling_mode .eq. 4 .and. mode4_operation .eq. 1
     *    .and. ntries .ge. mode4_trial_limit) go to 50
      ntries	= ntries+1
      proposal_component=0
      proposal_q2_bin=0
      proposal_xb_bin=0
      proposal_t_bin=0
      proposal_phi_bin=0
      proposal_density_ratio=1.
      proposal_balanced_to_legacy=0.
      proposal_xb=0.
      proposal_minus_t=0.
      proposal_t_jacobian=0.
      proposal_t_at_cosine_zero=0.
      if (sampling_mode .eq. 1 .and.
     *    survey_proposal_mode .eq. 1) then
         proposal_selector=myran()
         if (proposal_selector .ge. survey_legacy_fraction) then
            proposal_component=1
            survey_balanced_trials=survey_balanced_trials+1
         else
            survey_legacy_trials=survey_legacy_trials+1
         endif
      elseif (sampling_mode .eq. 1) then
         survey_legacy_trials=survey_legacy_trials+1
      elseif (sampling_mode .eq. 4) then
         proposal_selector=myran()
         if (proposal_selector .lt.
     *       mode4_component_core_fraction) then
            proposal_component=1
            mode4_core_trials=mode4_core_trials+1
         else
            mode4_legacy_trials=mode4_legacy_trials+1
            if (mode4_operation .eq. 1 .and.
     *          mode4_calibration_proposal .eq. 1) then
               call mode4_sample_complement(r_u,r_ep,uek,
     *            hadron_cosine_base,hadron_phi_base,
     *            mode4_ru_low,mode4_ru_high,
     *            mode4_rep_low,mode4_rep_high,
     *            mode4_ugamma_low,mode4_ugamma_high,
     *            mode4_hcos_low,mode4_hcos_high,
     *            mode4_hphi_origin,mode4_hphi_relative_low,
     *            mode4_hphi_relative_high)
            endif
         endif
      endif
      if(flag_ehel.eq.1) ehel = get_spin(idum)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hel       = ehel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c     calculate the energy of the electron at the scattering point
c     after making its way through the target.  First, randomly
c     choose the interaction point.
      
c     Change from units of r.l. to cm
      t_targ 	=  hydrogen_rad * t_targ / bfac 
      vertex_fraction=myran()
      targs	= t_targ*vertex_fraction

c     Change into proper coordinate system
      vertex_z	= vz + targs - t_targ / 2.0
 
c     calculates the distance that the electron stays in the target
      targp	= r_targ/sin(T0)
      temp	= (t_targ - targs)/cos(T0)
      
      if (temp .lt. targp) targp = temp
      
c     Back to units of r.l.      
      t_targ 	=  bfac * t_targ / hydrogen_rad 
      targs  	=  bfac * targs / hydrogen_rad
      targp  	=  bfac * targp / hydrogen_rad      

c     Now calculate the radiation loss

 22   xs	= myran()
      eloss	= xs**(1./targs)
      gxs	= 1.-eloss
      xtest	= myran()
      
      if (xtest.gt.gxs) go to 22
      incoming_loss_base=xs
      incoming_loss_accept_test=xtest
      
      es	= ebeam*(1.-eloss)
      incoming_loss_fraction=1.-es/ebeam
      
c     Cut off the incident energy at e_s=ebeam/4.

      if (es .lt. ebeam/4.) go to 20
      
      ps	= sqrt(es**2-mel**2)
      rs	= ps/es

      if (sampling_mode .eq. 4 .and.
     *    proposal_component .eq. 1) then
         r_u=mode4_ru_low+
     *      (mode4_ru_high-mode4_ru_low)*myran()
         uq2=uq2_min+uq2rng*r_u
         q2=1./uq2
         r_ep=mode4_rep_low+
     *      (mode4_rep_high-mode4_rep_low)*myran()
         ep=epmax-eprng*r_ep
      elseif (sampling_mode .eq. 4 .and.
     *        mode4_operation .eq. 1 .and.
     *        mode4_calibration_proposal .eq. 1) then
         uq2=uq2_min+uq2rng*r_u
         q2=1./uq2
         ep=epmax-eprng*r_ep
      elseif (proposal_component .eq. 1) then
         proposal_q2_bin=min(int(survey_nq2_bins*myran())+1,
     *      survey_nq2_bins)
         q2=survey_q2_edges(proposal_q2_bin)+
     *      (survey_q2_edges(proposal_q2_bin+1)-
     *       survey_q2_edges(proposal_q2_bin))*myran()
         proposal_xb_bin=min(int(survey_nxb_bins*myran())+1,
     *      survey_nxb_bins)
         proposal_xb=survey_xb_edges(proposal_xb_bin)+
     *      (survey_xb_edges(proposal_xb_bin+1)-
     *       survey_xb_edges(proposal_xb_bin))*myran()
         ep=es-q2/(2.*real(analysis_mp)*proposal_xb)
         uq2=1./q2
         r_u=(uq2-uq2_min)/uq2rng
         r_ep=(epmax-ep)/eprng
      else
         uq2=uq2_min+uq2rng*myran()
         q2=1./uq2
         r_u=(uq2-uq2_min)/uq2rng

c        calculate the energy and momentum of the scattered electron,
c        and calculate Q**2

         ep=epmax-eprng*myran()
         r_ep=(epmax-ep)/eprng
      endif
      
c     check to see if the scattered electron energy is below the
c     detector threshold

      if (ep .lt. ep_min .or. ep .gt. epmax)go to 20

      q0	= es-ep
      s		= q2/4/es/ep
      
c     cut off scattering at 90 degree

      if (s .gt. .5)go to 20
      
      th0	= 2.*asin(sqrt(s))
      theta	= th0*180./pi
      T0	= th0
      snt0	= sin(th0)
      cst0	= cos(th0)
      
c     check to see if the scattered electron energy is above
c     the pion threshold for this angle.

      ep_test	= (mp**2 + 2*mp*es - wg**2)/2./(mp + 2.*es*s)
      
      if (ep .gt. ep_test)go to 20
      
      pp	= sqrt(ep**2-mel**2)
      rp	= pp/ep
      qsq	= q2
      
      if (qsq .le. 0.)then
         write(6,*)' Main-1:, qsq =',qsq
         go to 20
      endif
      
      qvecx	=     -pp*sin(th0)
      qvecz	= ps - pp*cos(th0)

      w2	= mp**2 + 2*mp*q0 - q2
      
      if (w2 .lt. mp**2)go to 20
      
      epw	= sqrt(w2)
      
      if (th_opt.gt.10) then
        if (epw .lt. W_min_input) go to 20
      else
        if (epw .lt. wg+0.002)go to 20
      endif
      
c     calculate kinematic quantities needed for the Mo and Tsai calculation

      u0	= es - ep + mp
      pu	= ps**2 + pp**2 - 2*ps*pp*cst0
      
      if (pu .le. 0.)then
         write(6,*)' Main-2, pu**2 =',pu
         go to 20
      endif
      
      pu	= sqrt(pu)
      uu	= u0**2 - pu**2
      csths	= (ps - pp*cst0)/pu
      csthp	= (ps*cst0 - pp)/pu
      snths	= 1.-csths**2
      
      if (snths**2 .le. 0.)then
         write(6,*)' Main-3: snths =',snths
         go to 20
      endif
      
      snths	= sqrt(snths)
      snthp	= 1.-csthp**2
      
      if (snthp .le. 0.)then
         write(6,*)' Main-4: snthp**2 =',snthp
         go to 20
      endif
      
      snthp	= sqrt(snthp)
      ts	= acos(csths)
      tp	= acos(csthp)
      sp	= es*ep - ps*pp*cst0

      cstk1	= csths
      cstk2	= csthp
      csrnge 	= csrng
      
      if (cstk1 .lt. cstk2)then
         cstmp=cstk1
         cstk1=cstk2
         cstk2=cstmp
      endif
      
      if ((1.-cstk1) .lt. csrnge)       csrnge=1.-cstk1
      if ((cstk1-cstk2) .lt. 2.*csrnge) csrnge=0.5*(cstk1-cstk2)
      
      csrngb	= csrng/40.
      
      if (csrngb .gt. csrnge/5.) csrngb=csrnge/5.
      
      csran	= myran()
      rn1	= myran()
      photon_side_base=rn1
      
      if (rn1 .gt. 0.5)then
        rn1=1.
      else
        rn1=-1.
      endif
      
      rn2	= myran()
      photon_cosine_base=rn2
      delphi	= ps*cstk1*csrngb/pp/sqrt(1.-cst0**2)/sqrt(1.-cstk1**2)
      
      if (delphi .lt. pi/9.) delphi=pi/9.
      if (delphi .gt. 2.*pi) delphi=2.*pi
      if (cstk1 .gt. .995)   delphi=2.*pi
      
      delphi	= pi/9.
       
      if (csran .lt. reg1)then
         intreg	= 1
         cstk	= cstk1+(2.*rn2-1.)*csrngb
         mcfac	= csrngb /reg1
         photon_phi_base=myran()
         phik	= (photon_phi_base-0.5)*delphi
         mpfac	= delphi/2./pi
      elseif(csran .lt. reg2)then
         intreg	= 2
         cstk	= cstk2+(2.*rn2-1.)*csrngb
         mcfac	= csrngb /(reg2-reg1)
         photon_phi_base=myran()
         phik	= (photon_phi_base-0.5)*delphi
         mpfac	= delphi/2./pi
      elseif(csran .lt. reg3)then
         intreg	= 3
         cstk	= cstk1+rn1*(csrngb+rn2*(csrnge-csrngb))
         mcfac	= (csrnge-csrngb) /(reg3-reg2)
         photon_phi_base=myran()
         phik	= (photon_phi_base-0.5)*delphi
         mpfac	= delphi/2./pi
      elseif(csran .lt. reg4)then
         intreg	= 4
         cstk	= cstk2+rn1*(csrngb+rn2*(csrnge-csrngb))
         mcfac	= (csrnge-csrngb) /(reg4-reg3)
         photon_phi_base=myran()
         phik	= (photon_phi_base-0.5)*delphi
         mpfac	= delphi/2./pi
      else
         intreg	= 5
 45      cstk	= 2.*rn2-1.
         photon_phi_base=myran()
         phik	= 2.*pi*(photon_phi_base-0.5)
         if (abs(cstk-cstk1) .lt. csrnge .or.
     &       abs(cstk-cstk2) .lt. csrnge) then
         if(abs(phik).lt. delphi/2.) go to 45
        endif
           
c     combine mcfac and mpfac into one factor and set mpfac=1

        mcfac=(1.-csrnge*delphi/pi)/(1.-reg4)
        mpfac=1.
      endif

      Tk	= acos(cstk)
      sntk	= sin(Tk)
      
c     change the following on Jan. 19, 1999
c      phran=myran()
c      if (phran .lt. .2)then
c      else
c 48      phik=2*pi*(myran()-0.5)
c         if(abs(phik).lt. pi/180.)go to 48
c         mpfac=1.25
c      endif
c     end of jan 19, 1999 correction

      ekmax	= 0.5*(uu - wg**2)/(u0 - pu*cstk)
      
      if (ekmax .gt. ebeam)then
        write(6,*)' Main-5: ekmax =',ekmax
        ekmax=ebeam
      endif
      
c     choose ek by making a change of variables

 78   if (sampling_mode .eq. 4 .and.
     *    proposal_component .eq. 1) then
         uek=mode4_ugamma_low+
     *      (mode4_ugamma_high-mode4_ugamma_low)*myran()
      elseif (sampling_mode .eq. 4 .and.
     *        mode4_operation .eq. 1 .and.
     *        mode4_calibration_proposal .eq. 1) then
c        The guard-complement component was sampled once at trial start.
         uek=uek
      else
         uek=myran()
      endif
 
      if (uek .lt. 0.1E-20)then
        write(6,*)' Main-6: uek =',uek
        go to 20
      endif
      
      ek	= -alog(uek)/kexp
      
      if (ek .gt. ekmax)go to 20

      if (sampling_mode .eq. 1 .and.
     *    survey_proposal_mode .eq. 1) then
         call radiative_t_mapping(es,ep,th0,ek,cstk,mpi,0.938D0,
     *      analysis_mp,proposal_t_at_cosine_zero,
     *      proposal_t_jacobian,proposal_mapping_status)
         if (proposal_mapping_status .ne. 0)go to 20
      endif

      if (sampling_mode .eq. 4 .and.
     *    proposal_component .eq. 1) then
         hadron_cosine_base=mode4_hcos_low+
     *      (mode4_hcos_high-mode4_hcos_low)*myran()
         csthcm=-1.+2.*hadron_cosine_base
         mode4_hphi_relative=mode4_hphi_relative_low+
     *      (mode4_hphi_relative_high-mode4_hphi_relative_low)*
     *      myran()
         hadron_phi_base=modulo(mode4_hphi_origin+
     *      mode4_hphi_relative,1.)
         phicm=360.*hadron_phi_base
      elseif (sampling_mode .eq. 4 .and.
     *        mode4_operation .eq. 1 .and.
     *        mode4_calibration_proposal .eq. 1) then
         csthcm=-1.+2.*hadron_cosine_base
         phicm=360.*hadron_phi_base
      elseif (proposal_component .eq. 1) then
         proposal_t_bin=min(int(survey_nt_bins*myran())+1,
     *      survey_nt_bins)
         proposal_minus_t=survey_t_edges(proposal_t_bin)+
     *      (survey_t_edges(proposal_t_bin+1)-
     *       survey_t_edges(proposal_t_bin))*myran()
         csthcm=(proposal_t_at_cosine_zero-proposal_minus_t)/
     *      proposal_t_jacobian
         if (csthcm .lt. -1. .or. csthcm .gt. 1.)go to 20
         hadron_cosine_base=(csthcm+1.)/2.
         proposal_phi_bin=min(int(survey_nphi_bins*myran())+1,
     *      survey_nphi_bins)
         phicm=survey_phi_edges(proposal_phi_bin)+
     *      (survey_phi_edges(proposal_phi_bin+1)-
     *       survey_phi_edges(proposal_phi_bin))*myran()
         hadron_phi_base=phicm/360.
      else
         hadron_cosine_base=myran()
         hadron_phi_base=myran()
         csthcm=-1.+2.*hadron_cosine_base
         phicm=360.*hadron_phi_base
         if (sampling_mode .eq. 1 .and.
     *       survey_proposal_mode .eq. 1) then
            proposal_minus_t=proposal_t_at_cosine_zero-
     *         proposal_t_jacobian*csthcm
         endif
      endif
      
c********************************** Ek < Delta ****************************************      

      if (ek .lt. delta)then
        intreg=6
        
c      print *, th0,qsq,epw,csthcm,phicm
c     calculate the non-radiative cross section

        call dsigma(th0,qsq,epw,csthcm,phicm,th_opt,epirea,res_opt,sigma0
     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
     
        if (sigma0.le.0.) go to 20
c       PRINT *,'DSIGMA=',
c     *  th0,qsq,epw,csthcm,phicm,th_opt,epirea,res_opt,sigma0
c     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel
 
c       print *, 'AAO_RAD s0,u,t,l,i:',sigma0,sigu,sigt,sigl,sigi
     
c      if (abs(epw-1.232).lt.0.010.and.abs(qsq-0.35).lt.0.1) then 
c      print *, qsq,epw,csthcm,phicm
c      print *, 'AO:',sigma0,sigu,sigl,sigt,sigi
c      call dsigma(th0,qsq,epw,csthcm,phicm,4,epirea,res_opt,sigma0
c     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
c      print *, 'MAID:',sigma0,sigu,sigl,sigt,sigi
c      print *, ''
c     endif

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     The AO program calculates sigu,sigl,sigt and sigi with
c     the kinematic factors included. To convert them to the
c     form needed for the above problem we must use

        epeps	= 1. + 2.*(1 + q0**2/qsq)*s/(1-s)
      
        if (epeps .le. 1.)then
           write(6,*)' Main-6: epeps-inv =',epeps
           go to 20
        endif
        
        epeps	= 1./epeps
        kfac	= (w2 - mp**2)/2./mp
        
c        cfac=1/Gamma_T
c        cfac=2*pi**2*qsq*(es/ep)*(1-epeps)/kfac/alpha
c        sigu=sigu*cfac
c        sigl=sigl*cfac
c        sigt=sigt*cfac
c        sigi=sigi*cfac

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Mo and Tsai give a formula for the cross section in terms
c     of amplitudes, f and g.  Here we try to extract f and g
c     from the AO cross section.  I may not have this part right.

c     sig0=4.*pi**2*alpha/mp**2

c     The following would work if there were no polarization terms

c     f=2.*kfac/(mp**3)*qsq/(qsq+q0**2)*(sigu+sigl)/sig0
c     g=(2.*kfac/mp)*(sigu)/sig0

c     I don\'t know what to do with the polarization terms: these give
c     the dependence on the center of mass phi decay angle.  I try
c     the following:

        nu	= (w2 + qsq - mp**2)/2/mp
        
        f = (1./(2*pi**2*alpha*mp))*(kfac/(1.+nu**2/qsq))*(sigu+sigl)
        g = mp/(2*pi**2*alpha) * kfac * sigu
        
        fkt = (w2-mp**2+mpi**2)/2./epw
        fkt = sqrt(fkt**2-mpi**2)*2.*epw/(w2-mp**2)
        f   = f*fkt
        g   = g*fkt

c     Calculate the non-radiative cross section

        signr = 2*(alpha*ep/qsq)**2
     &   * (f*mp*cos(th0/2)**2+2*g*sin(th0/2)**2/mp)
     
c        asig = fkt*(sigu+epeps*sigl)/cfac
c        print *, 'Compare x-sections',signr,asig 
     
c     Modulate the cross section with the polarization and interference terms

        signr = signr * (1.+(epeps*sigt*cos(phicm*pi/90.)
     &  + sqrt(epeps*(1.+epeps)/2)*sigi*cos(phicm*pi/180.)
     &  + ehel*sqrt(epeps*(1.-epeps)/2)*sigip*sin(phicm*pi/180.))
     &  /(sigu+epeps*sigl))

        if (signr .le. 0.)then
           write(6,*)'signr,sigu,sigl,sigt,sigi,epeps',
     &     signr,sigu,sigl,sigt,sigi,epeps
           write(6,*)'qsq,epw,th0,csthcm,phicm,epirea',
     &      qsq,epw,th0,csthcm,phicm,epirea
        endif

c     Calculate the radiative correction factor deltar for the cross
c      section.  This includes vertex corrections and the integration
c     up to photon energies of delta.

        x1	= (ep-es)/ep
        x2	= (es-ep)/es
        s1	= spence(x1)
        s2	= spence(x2)
        
        deltar	= -(alpha/pi)*(28./9. -(13./6.)*dlog(2.*sp/mel**2)
     &   -s1-s2)
        delinf	= -(alpha/pi)*dlog(es*ep/delta**2)*(dlog(2.*sp/mel**2)-1.)
      
        sigr1	= signr*(1.+deltar)*exp(delinf)
        
c	calculate average differential cross section in the region
c	from ek=0 to delta, and from cos(thetak)=-1 to 1.

	sigr	= sigr1/delta/4./pi
	
        if (sigr .gt. 0.) then
          go to 28
        else
          go to 20
        endif
        
c       end of section for calculation with ek < delta.
c       Normally, go to statement 28

      endif

      sdotk = es*ek - ps*ek*cstk*csths - ps*ek*sntk*snths*cos(phik)
      pdotk = ep*ek - pp*ek*cstk*csthp - pp*ek*sntk*snthp*cos(phik)

      sigr	= sigma(ek,Tk,csthcm,phicm,ehel)
      
      if (sigr .le. 0.) go to 20
      
28    jacob	= exp(kexp*ek)/kexp/(2*es*ep)*q2**2
      sigr	= sigr*jacob
      
      call missm(epirea,ebeam,es,ep,th0,ek,cstk,phik,mpi_s,
     *  ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     *  csthcm,phicm,wreal,mm2)
      
      if (abs(mm2-mm_exp) .gt. mm_cut)go to 20

c     Compare sigr to the sigr_max to determine whether to generate
c     an event.

      sigr	= mcfac*mpfac*sigr

c     Fixed-trial mode constructs exactly one post-external-loss candidate
c     for every internally valid proposal.  It records the original proposal
c     integrand and gives invalid final candidates zero observed contribution.

      if (sampling_mode .eq. 1 .and.
     *    survey_proposal_mode .eq. 1) then
         proposal_xb=0.
         if (es-ep .gt. 0.) then
            proposal_xb=q2/(2.*real(analysis_mp)*real(es-ep))
         endif
         proposal_minus_t=ppx**2+ppy**2+ppz**2
     *      -(real(analysis_mp)-eprot)**2
         call balanced_proposal_ratio(q2,proposal_xb,
     *      proposal_minus_t,phicm,proposal_t_jacobian,
     *      uq2rng,eprng,analysis_mp,survey_legacy_fraction,
     *      survey_nq2_bins,survey_q2_edges,survey_nxb_bins,
     *      survey_xb_edges,survey_nt_bins,survey_t_edges,
     *      survey_nphi_bins,survey_phi_edges,
     *      proposal_balanced_to_legacy,proposal_density_ratio)
         sigr=sigr*proposal_density_ratio
      endif

c     Mode 4 samples a full-support mixture in the five learned native
c     coordinates.  Correct every point by q_legacy/q_mixture, including
c     legacy-component points that happen to lie inside the learned core.

      if (sampling_mode .eq. 4) then
         call mode4_proposal_ratio(r_u,r_ep,uek,
     *      hadron_cosine_base,hadron_phi_base,
     *      mode4_core_fraction,mode4_core_volume,
     *      mode4_ru_low,mode4_ru_high,
     *      mode4_rep_low,mode4_rep_high,
     *      mode4_ugamma_low,mode4_ugamma_high,
     *      mode4_hcos_low,mode4_hcos_high,
     *      mode4_hphi_origin,mode4_hphi_relative_low,
     *      mode4_hphi_relative_high,proposal_density_ratio,
     *      mode4_inside_core)
         sigr=sigr*proposal_density_ratio
         mode4_internal_rows=mode4_internal_rows+1
      endif

c     A radiative stratum is defined from the final LUND electron and proton.
c     Construct one complete post-external-loss candidate before applying the
c     stratum indicator and acceptance-rejection.  This makes outgoing loss
c     part of the sampled probability space rather than an afterthought.

      if (sampling_mode .eq. 4) then
         ep_sav=ep
         call build_final_candidate(ebeam,es,ep_sav,ep,ep_min,th0,
     *      theta,ek,cstk,phik,mpi_s,csthcm,phicm,targp,delta,
     *      ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     *      wreal,mm2,epw,phir,outgoing_loss_base,
     *      outgoing_loss_accept_test,candidate_status)
         if (candidate_status .ne. 0)go to 30
         mode4_final_candidates=mode4_final_candidates+1
         call observed_analysis_coordinates(ebeam,analysis_mp,
     *      q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs,
     *      coordinate_status)
         if (coordinate_status .ne. 0)go to 30
         call mode4_stratum_membership(q2_obs,xb_obs,minus_t_obs,
     *      phi_obs,w_obs,y_obs,mode4_q2_low,mode4_q2_high,
     *      mode4_xb_low,mode4_xb_high,mode4_t_low,mode4_t_high,
     *      mode4_phi_low,mode4_phi_high,mode4_w_minimum,
     *      mode4_apply_y_max,mode4_y_maximum,
     *      mode4_stratum_status)
         if (mode4_stratum_status .ne. 0)go to 30
         mode4_target_candidates=mode4_target_candidates+1
         if (proposal_component .eq. 1) then
            mode4_core_targets=mode4_core_targets+1
         else
            mode4_legacy_targets=mode4_legacy_targets+1
         endif
         if (mode4_operation .eq. 1) then
            if (mode4_calibration_proposal .eq. 1) then
               if (proposal_component .eq. 1) then
                  mode4_component_weight=(mode4_core_fraction+
     *               mode4_legacy_fraction*mode4_core_volume)/
     *               mode4_component_core_fraction
               else
                  mode4_component_weight=(mode4_legacy_fraction*
     *               (1.-mode4_core_volume))/
     *               (1.-mode4_component_core_fraction)
               endif
            elseif (proposal_component .eq. 1) then
               mode4_component_weight=mode4_core_fraction/
     *               mode4_component_core_fraction
            else
               mode4_component_weight=mode4_legacy_fraction/
     *            (1.-mode4_component_core_fraction)
            endif
            mode4_trial_integrand=sigr*mode4_component_weight
            sig_tot=sig_tot+dble(mode4_trial_integrand)
            write(19,9003)ntries,proposal_component,
     *         mode4_inside_core,proposal_density_ratio,sigr,
     *         mode4_component_weight,r_u,r_ep,uek,
     *         hadron_cosine_base,hadron_phi_base,
     *         q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs
            go to 30
         endif
         sig_tot=sig_tot+dble(sigr)
         sig_ratio=sigr/sigr_max
         rtest=myran()
         mcall=sig_ratio
         stest=sig_ratio-mcall
         if (stest .gt. rtest)mcall=mcall+1
         if (mcall .gt. mcall_max)mcall_max=mcall
         if (mcall .eq. 0)go to 30
         if (mcall .lt. 100)then
            dismc(intreg,mcall)=dismc(intreg,mcall)+1
         else
            dismc(intreg,100)=dismc(intreg,100)+1
         endif
         lundmass(1)=MEL
         lundmass(4)=0.
         lundmass(2)=MP
         lundmass(3)=MPI
         do j=1,mcall
            write(13,*) "4 1 1 0 0 ",
     *         q2/(2.0*MP*(ebeam-p(1,4))),0,epw,
     *         q2,(ebeam-p(1,4))
            write(13,*)1,q(1),1,pdgid(1),0,0,
     *         p(1,1),p(1,2),p(1,3),p(1,4),
     *         lundmass(1),vertex_x,vertex_y,vertex_z
            write(13,*)2,q(2),1,pdgid(2),0,0,
     *         p(2,1),p(2,2),p(2,3),p(2,4),
     *         lundmass(2),vertex_x,vertex_y,vertex_z
            write(13,*)3,q(3),1,pdgid(3),0,0,
     *         p(3,1),p(3,2),p(3,3),p(3,4),
     *         lundmass(3),vertex_x,vertex_y,vertex_z
            write(13,*)4,q(4),1,pdgid(4),0,0,
     *         p(4,1),p(4,2),p(4,3),p(4,4),
     *         lundmass(4),vertex_x,vertex_y,vertex_z
            nevent=nevent+1
            if (proposal_component .eq. 1) then
               mode4_core_events=mode4_core_events+1
            else
               mode4_legacy_events=mode4_legacy_events+1
            endif
            write(17,9002)nevent,proposal_component,
     *         proposal_density_ratio,r_u,r_ep,uek,
     *         hadron_cosine_base,hadron_phi_base,
     *         q2_obs,xb_obs,minus_t_obs,
     *         phi_obs,w_obs,y_obs,sigr
         enddo
         go to 124
      endif

      if (sampling_mode .eq. 1) then
         ep_sav=ep
         q2_lep=q2
         w_lep=epw
         xb_lep=0.
         if (es-ep_sav .gt. 0.) then
            xb_lep=q2_lep/(2.*mp*(es-ep_sav))
         endif
         wreal_h=wreal
         ppx_h=ppx
         ppy_h=ppy
         ppz_h=ppz
         eprot_h=eprot
         q2_hard=mp**2+2.*mp*(es-ep_sav-ek)-wreal_h**2
         xb_hard=0.
         if (es-ep_sav-ek .gt. 0.) then
            xb_hard=q2_hard/(2.*mp*(es-ep_sav-ek))
         endif
         minus_t_hard=ppx_h**2+ppy_h**2+ppz_h**2
     *      -(real(analysis_mp)-eprot_h)**2

         survey_sig_internal_tot=survey_sig_internal_tot+dble(sigr)
         survey_internal_contribution=analysis_volume*dble(sigr)
         survey_observed_contribution=0.
         sigr_observed=0.
         final_valid=0
         candidate_status=0
         coordinate_status=1
         q2_obs=0.
         xb_obs=0.
         minus_t_obs=0.
         phi_obs=0.
         w_obs=0.
         y_obs=0.
         phir=0.
         do jj=1,4
            p(1,jj)=0.
            p(2,jj)=0.
         enddo

         call build_final_candidate(ebeam,es,ep_sav,ep,ep_min,th0,
     *      theta,ek,cstk,phik,mpi_s,csthcm,phicm,targp,delta,
     *      ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     *      wreal,mm2,epw,phir,outgoing_loss_base,
     *      outgoing_loss_accept_test,candidate_status)

         ep_out=ep
         outgoing_loss_fraction=1.-ep_out/ep_sav
         if (candidate_status .eq. 0) then
            call observed_analysis_coordinates(ebeam,analysis_mp,
     *         q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs,
     *         coordinate_status)
            if (coordinate_status .eq. 0) then
               final_valid=1
               survey_final_valid=survey_final_valid+1
               sigr_observed=sigr
               sig_tot=sig_tot+dble(sigr_observed)
               survey_observed_contribution=
     *            analysis_volume*dble(sigr_observed)
            else
               candidate_status=3
            endif
         endif

         survey_rows=survey_rows+1
         write(16,9001)survey_replica,ntries,final_valid,
     *      candidate_status,intreg,ehel,
     *      proposal_component,proposal_q2_bin-1,proposal_xb_bin-1,
     *      proposal_t_bin-1,proposal_phi_bin-1,
     *      proposal_density_ratio,
     *      r_u,r_ep,uek,csran,
     *      vertex_fraction,incoming_loss_base,
     *      incoming_loss_accept_test,photon_side_base,
     *      photon_cosine_base,photon_phi_base,
     *      hadron_cosine_base,hadron_phi_base,
     *      outgoing_loss_base,outgoing_loss_accept_test,
     *      cstk,phik*180./pi,csthcm,phicm,
     *      es,ep_sav,ep_out,real(ek),
     *      p(1,1),p(1,2),p(1,3),p(1,4),
     *      p(2,1),p(2,2),p(2,3),p(2,4),
     *      q2_lep,xb_lep,w_lep,
     *      q2_hard,xb_hard,wreal_h,minus_t_hard,
     *      q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs,
     *      sigr,sigr_observed,survey_internal_contribution,
     *      survey_observed_contribution,phir*180./pi,
     *      incoming_loss_fraction,outgoing_loss_fraction
         go to 30
      endif

      sig_ratio	= sigr/sigr_max
      sig_tot	= sig_tot + sigr
      
c     Choose the number of times, mcall, to call the routine used
c     to calculate kinematic quantities for the n-tuple.

      rtest	= myran()
      mcall	= sig_ratio
      stest	= sig_ratio - mcall
      
      if (stest .gt. rtest)mcall = mcall + 1
       
      if (mcall .gt. mcall_max) mcall_max=mcall
      if (mcall .gt. 10)then
c        write(6,*)' mcall =',mcall,' intreg=',intreg
c        write(6,*)es,ep,th0,ek
c        write(6,*)'cstk1,cstk2,cstk,phik',cstk1,cstk2,cstk,phik
c        write(6,*)' csrnge=',csrnge
c        write(6,*) ppx,ppy,ppz,eprot,ekx,eky,ekz
c        write(6,*)csthcm,phicm,wreal,mm2

         csdotk	= cstk*csths+sntk*snths*cos(phik)
         cpdotk	= cstk*csthp+sntk*snthp*cos(phik)
         cqdotk	= (ps*csdotk-pp*cpdotk)/sqrt(ps**2+pp**2-2*ps*pp*cos(th0))
         
c         write(6,*)'es,ep,ek',es,ep,ek
c         write(6,*)'cstk,phik,csdotk,cpdotk'
c     +   ,cstk,phik,csdotk,cpdotk
c         write(6,*)'cqdotk',cqdotk

         write(12,*)' mcall =',mcall,' intreg=',intreg
         write(12,*)'es,ep,ek',es,ep,ek
         write(12,*)'cstk,phik,csdotk,cpdotk'
     +   ,cstk,phik,csdotk,cpdotk
         write(12,*)'cqdotk',cqdotk
      endif
      
c     If mcall .gt. 0 generate mcall n-tuple events

      if (mcall .eq. 0)go to 30

      ep_sav	= ep

      if (mcall .lt. 100)then
         dismc(intreg,mcall)=dismc(intreg,mcall)+1
      else
         dismc(intreg,100)=dismc(intreg,100)+1
      endif

      do j=1,mcall

c     Calculate the radiation loss for the electron leaving the target

         call build_final_candidate(ebeam,es,ep_sav,ep,ep_min,th0,
     *      theta,ek,cstk,phik,mpi_s,csthcm,phicm,targp,delta,
     *      ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     *      wreal,mm2,epw,phir,outgoing_loss_base,
     *      outgoing_loss_accept_test,candidate_status)

         if (candidate_status .eq. 1)then
            sig_tot = sig_tot-sigr_max
         endif
         if (candidate_status .ne. 0)go to 24

c     Calculate the members of the n-tuple and ouput it to the rz file.

         ntp(1)	= es
         ntp(2)	= ep
         ntp(3) = theta
         ntp(4) = epw
         ntp(5) = wreal
         ntp(6) = ppx
         ntp(7) = ppy
         ntp(8) = ppz
         ntp(9) = eprot
         ntp(10)= ppix
         ntp(11)= ppiy
         ntp(12)= ppiz
         ntp(13)= epi
         ntp(14)= csthcm
         ntp(15)= phicm
         ntp(16)= mm2
         ntp(17)= ek
         ntp(18)= cstk
         ntp(19)= phik*180./pi
         ntp(20)= qvecx
         ntp(21)= qvecz
         ntp(22)= q0
         ntp(23)= cst0
         ntp(24)= ekx
         ntp(25)= eky
         ntp(26)= ekz
         ntp(27)= vertex_x
         ntp(28)= vertex_y
         ntp(29)= vertex_z
         ntp(30)= qsq                  
         ntp(31)= ehel
         ntp(32)= asym_p

Clcs         call hfn(10,ntp)
c         nevent	= nevent+1
        
         do jj = 1,npart
           v(jj,1) = vertex_x
           v(jj,2) = vertex_y
           v(jj,3) = vertex_z
         enddo
*
*   Cut of theta>16 deg applied. For e1-dvcs eta meson production.
*   12/12/2011, Ivan Bedlinskiy
*         
c         if (theta.lt.16) goto 24

         lundmass(1) = MEL
         lundmass(4) = 0
         if ((epirea .eq. 1).or.(epirea.eq.5)) then
            lundmass(2) = MP
            lundmass(3) = MPI
         else 
            lundmass(2) = MPIP
            if(epirea.eq.2) lundmass(3) = MP
            if(epirea.eq.3) lundmass(3) = MN
         endif

         write(13,*) "4 1 1 0 0 ",q2/(2.0*MP*(ebeam-p(1,4))), 0, epw, 
     *   q2, (ebeam-p(1,4))
         write(13,*) 1, q(1), 1, pdgid(1), 0, 0,
     *               p(1,1), p(1,2), p(1,3), p(1,4),
     *               lundmass(1), vertex_x, vertex_y, vertex_z
         write(13,*) 2, q(2), 1, pdgid(2), 0, 0,
     *               p(2,1), p(2,2), p(2,3), p(2,4),
     *               lundmass(2), vertex_x, vertex_y, vertex_z
         write(13,*) 3, q(3), 1, pdgid(3), 0, 0,
     *               p(3,1), p(3,2), p(3,3), p(3,4),
     *               lundmass(3), vertex_x, vertex_y, vertex_z
         write(13,*) 4, q(4), 1, pdgid(4), 0, 0,
     *               p(4,1), p(4,2), p(4,3), p(4,4),
     *               lundmass(4), vertex_x, vertex_y, vertex_z

C         print *,'Event to be written'
Cvpk         call bos_out	! Pack the BOS banks and write out to file   
         nevent	= nevent+1
C         print *,'Events=',nevent,nprint
         goto 124
 24      continue
       enddo

c     Talk to the user every now and then
 124    continue
       ntell	= nevent/nprint-ntold
       
       if (ntell .gt. 0)then
          write(6,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max
          write(12,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max
     
         events	= nevent
         tries	= ntries
         sig_int= events/tries
         
c     photon phase space = 4*pi*delta-omega=4*pi (after change of variables)
c     hadron phase space = 4*pi
c     electron phase space =2*pi*ucrng*eprng

         sig_int = sig_int*sigr_max*(4.*pi)**2*(2.*pi*uq2rng*eprng)
         sig_sum = sig_tot*(4.*pi)**2*(2.*pi*uq2rng*eprng)/tries

         write(6,*)' Integrated cross section (MC, numerical) ='
     +   ,sig_int,sig_sum, ' mu-barns'
         write(6,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
         write(12,*)' Integrated cross section =',sig_int,sig_sum
     *    ,' micro-barns'
         write(12,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'

          call timex(itime2)
          
          WRITE(12,*)' t_elapse,itime1,iteme2=',
     *                t_elapse,itime1,itime2          
          WRITE(6,*)' t_elapse,itime1,iteme2=',
     *                t_elapse,itime1,itime2          
C          t_elapse	= t_elapse+INT(itime2-itime1)
          t_elapse	= t_elapse+itime2-itime1
          itime1	= itime2
          write(6,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
          write(12,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
          
          ntold	= ntold+1
          
       endif

c     Do we have enough events in the n-tuple?

       if (sampling_mode .eq. 4 .and. nevent .ge. nmax)go to 50
       if (sampling_mode .ne. 4 .and. nevent .gt. nmax)go to 50
       
 30    go to 20

c     Close out the n-tuple file

 50   continue
Clcs      call hrout(0,icycle,' ')
Clcs      call hrend('aaoradgen')

      if (sampling_mode .eq. 4 .and.
     *    ntries .ne. mode4_last_heartbeat) then
         write(18,9010)ntries,nevent,mode4_core_trials,
     *      mode4_legacy_trials,mode4_internal_rows,
     *      mode4_final_candidates,mode4_target_candidates,
     *      mode4_core_targets,mode4_legacy_targets,
     *      mode4_core_events,mode4_legacy_events,mcall_max
         flush(18)
         flush(12)
         flush(13)
         if (mode4_operation .eq. 0)flush(17)
         if (mode4_operation .eq. 1)flush(19)
         mode4_last_heartbeat=ntries
      endif

      if (sampling_mode .eq. 1)then
         tries=ntries
         events=0.
         sig_int=0.
         sig_sum=sig_tot*analysis_volume/tries
         write(6,*)' Fixed-trial survey proposals =',ntries
         write(6,*)' Internally valid survey rows =',survey_rows
         write(6,*)' Valid final survey candidates =',
     *      survey_final_valid
         write(6,*)' Survey legacy/balanced trials =',
     *      survey_legacy_trials,survey_balanced_trials
         write(6,*)' Survey internal cross section =',
     *      survey_sig_internal_tot*analysis_volume/tries,
     *      ' micro-barns'
         write(6,*)' Survey observed cross section =',sig_sum,
     *      ' micro-barns'
         write(12,*)' Fixed-trial survey proposals =',ntries
         write(12,*)' Internally valid survey rows =',survey_rows
         write(12,*)' Valid final survey candidates =',
     *      survey_final_valid
         write(12,*)' Survey legacy/balanced trials =',
     *      survey_legacy_trials,survey_balanced_trials
         write(12,*)' Survey internal cross section =',
     *      survey_sig_internal_tot*analysis_volume/tries,
     *      ' micro-barns'
         write(12,*)' Survey observed cross section =',sig_sum,
     *      ' micro-barns'
      endif
      if (sampling_mode .eq. 4)then
         tries=ntries
         events=nevent
         sig_sum=sig_tot*analysis_volume/tries
         if (mode4_operation .eq. 1) then
            sig_int=0.
         else
            sig_int=events/tries*sigr_max*analysis_volume
         endif
         write(6,*)' Mode-4 proposals/events =',ntries,nevent
         write(6,*)' Mode-4 core/legacy trials =',
     *      mode4_core_trials,mode4_legacy_trials
         write(6,*)' Mode-4 final/target candidates =',
     *      mode4_final_candidates,mode4_target_candidates
         write(6,*)' Mode-4 core/legacy events =',
     *      mode4_core_events,mode4_legacy_events
         write(6,*)' Mode-4 integrated cross section =',
     *      sig_int,sig_sum,' micro-barns'
         write(12,*)' Mode-4 proposals/events =',ntries,nevent
         write(12,*)' Mode-4 core/legacy trials =',
     *      mode4_core_trials,mode4_legacy_trials
         write(12,*)' Mode-4 final/target candidates =',
     *      mode4_final_candidates,mode4_target_candidates
         write(12,*)' Mode-4 core/legacy events =',
     *      mode4_core_events,mode4_legacy_events
         write(12,*)' Mode-4 integrated cross section =',
     *      sig_int,sig_sum,' micro-barns'
      endif

      close(12)
      close(13)
      if (sampling_mode .eq. 1)close(16)
      if (sampling_mode .eq. 4 .and. mode4_operation .eq. 0)close(17)
      if (sampling_mode .eq. 4 .and. mode4_operation .eq. 1)close(19)
      if (sampling_mode .eq. 4)close(18)

      open(unit=14,file=file_sum)
      
      write(14,*)' AO Calculation of Single Pion Production'
      write(14,*)' AO Calculation of Single Pion Production'

      write(14,*)' Starting time: ',ctime
      write(14,*)' Epirea (1 for pi0, 3 for pi+) =',epirea
      write(14,*)'Target thickness =',t_targ*3./4.,' (r.l.)'
      write(14,*)' Incident electron energy =',ebeam,' GeV'

      write(14,*)'Electron Q**2 limits:',q2_min,q2_max
      write(14,*)'lower and upper limit for scattered electron',
     * ' energy(GeV):',ep_min,epmax
      write(14,*)' Minimum photon energy for integration (delta):',delta
      write(14,*)'sigr_max,ek_max,tk_max,csthcm_max,phicm_max',
     * sigr_max,ek_max,tk_max,csthcm_max,phicm_max
      write(14,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max
      write(14,*)' Missing-mass squared cut at:',mm_cut
      write(14,*)' Integrated cross section =',sig_int,sig_sum
     *    ,' micro-barns'
      if (sampling_mode .eq. 0)then
         write(14,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *      ,' seconds'
      elseif (sampling_mode .eq. 1)then
         write(14,*)' Survey rows, final valid =',
     *      survey_rows,survey_final_valid
         write(14,*)' Survey internal cross section =',
     *      survey_sig_internal_tot*analysis_volume/tries
      else
         write(14,*)' Mode-4 core/legacy trials =',
     *      mode4_core_trials,mode4_legacy_trials
         write(14,*)' Mode-4 final/target candidates =',
     *      mode4_final_candidates,mode4_target_candidates
         write(14,*)' Mode-4 core/legacy events =',
     *      mode4_core_events,mode4_legacy_events
         write(14,*)' Mode-4 operation/trial limit =',
     *      mode4_operation,mode4_trial_limit
      endif
      write(14,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
      if ((sampling_mode .eq. 0 .or. sampling_mode .eq. 4)
     *    .and. nevent .gt. 0)then
         write(14,*)' CPU time/event = ', t_elapse/nevent,' sec'
      endif

      open(unit=15,file=file_norm)
      write(15,*)'generator=aao_rad'
      write(15,*)'sampling_mode=',sampling_mode
      write(15,*)'fixed_trial_survey=',merge(1,0,sampling_mode.eq.1)
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
      write(15,*)'delta=',delta
      write(15,*)'epirea=',epirea
      write(15,*)'th_opt=',th_opt
      write(15,*)'res_opt=',res_opt
      write(15,*)'fmcall=',fmcall
      if (sampling_mode .eq. 1)then
         write(15,*)'survey_schema=aao-rad-survey-v2'
         write(15,*)'survey_file=',file_survey
         write(15,*)'survey_seed=',survey_seed
         write(15,*)'survey_replica=',survey_replica
         write(15,*)'survey_ntrials_requested=',survey_ntrials
         write(15,*)'survey_proposal_mode=',survey_proposal_mode
         if (survey_proposal_mode .eq. 1)then
            write(15,*)'survey_proposal=analysis_balanced_mixture'
         else
            write(15,*)'survey_proposal=legacy'
         endif
         write(15,*)'survey_legacy_fraction=',survey_legacy_fraction
         write(15,*)'survey_balanced_fraction=',
     *      survey_balanced_fraction
         write(15,*)'survey_legacy_trials=',survey_legacy_trials
         write(15,*)'survey_balanced_trials=',survey_balanced_trials
         write(15,*)'survey_balance_q2_bins=',survey_nq2_bins
         write(15,*)'survey_balance_xb_bins=',survey_nxb_bins
         write(15,*)'survey_balance_t_bins=',survey_nt_bins
         write(15,*)'survey_balance_phi_bins=',survey_nphi_bins
         write(15,*)'survey_rows=',survey_rows
         write(15,*)'survey_final_valid=',survey_final_valid
         write(15,*)'survey_phase_volume=',analysis_volume
         write(15,*)'survey_internal_sig_sum=',
     *      survey_sig_internal_tot*analysis_volume/tries
         write(15,*)'survey_observed_sig_sum=',sig_sum
         write(15,*)'survey_emits_lund=0'
      endif
      if (sampling_mode .eq. 4)then
         write(15,*)'mode4_schema=aao-rad-mode4-v3'
         write(15,*)'mode4_operation=',mode4_operation
         write(15,*)'mode4_calibration_proposal=',
     *      mode4_calibration_proposal
         write(15,*)'mode4_trial_limit=',mode4_trial_limit
         write(15,*)'mode4_heartbeat_interval=',
     *      mode4_heartbeat_interval
         write(15,*)'mode4_component_core_fraction=',
     *      mode4_component_core_fraction
         write(15,*)'mode4_heartbeat_file=',file_mode4_heartbeat
         if (mode4_operation .eq. 0) then
            write(15,*)'mode4_event_file=',file_mode4
            write(15,*)'mode4_calibration_file=none'
         else
            write(15,*)'mode4_event_file=none'
            write(15,*)'mode4_calibration_file=',
     *         file_mode4_calibration
         endif
         write(15,*)'mode4_seed=',mode4_seed
         write(15,*)'mode4_replica=',mode4_replica
         write(15,*)'stratum_flat_index=',stratum_flat
         write(15,*)'stratum_iq2=',stratum_iq2
         write(15,*)'stratum_ixb=',stratum_ixb
         write(15,*)'stratum_it=',stratum_it
         write(15,*)'stratum_iphi=',stratum_iphi
         write(15,*)'mode4_q2_low=',mode4_q2_low
         write(15,*)'mode4_q2_high=',mode4_q2_high
         write(15,*)'mode4_xb_low=',mode4_xb_low
         write(15,*)'mode4_xb_high=',mode4_xb_high
         write(15,*)'mode4_t_low=',mode4_t_low
         write(15,*)'mode4_t_high=',mode4_t_high
         write(15,*)'mode4_phi_low=',mode4_phi_low
         write(15,*)'mode4_phi_high=',mode4_phi_high
         write(15,*)'mode4_w_minimum=',mode4_w_minimum
         write(15,*)'mode4_apply_y_max=',mode4_apply_y_max
         write(15,*)'mode4_y_maximum=',mode4_y_maximum
         write(15,*)'mode4_core_fraction=',mode4_core_fraction
         write(15,*)'mode4_legacy_fraction=',mode4_legacy_fraction
         write(15,*)'mode4_base_padding=',mode4_base_padding
         write(15,*)'mode4_core_volume=',mode4_core_volume
         write(15,*)'mode4_core_trials=',mode4_core_trials
         write(15,*)'mode4_legacy_trials=',mode4_legacy_trials
         write(15,*)'mode4_internal_rows=',mode4_internal_rows
         write(15,*)'mode4_final_candidates=',mode4_final_candidates
         write(15,*)'mode4_target_candidates=',mode4_target_candidates
         write(15,*)'mode4_core_targets=',mode4_core_targets
         write(15,*)'mode4_legacy_targets=',mode4_legacy_targets
         write(15,*)'mode4_core_events=',mode4_core_events
         write(15,*)'mode4_legacy_events=',mode4_legacy_events
         write(15,*)'mode4_phase_volume=',analysis_volume
         if (mode4_operation .eq. 0 .and. events .gt. 0.) then
            write(15,*)'mode4_event_weight_microbarn=',sig_sum/events
            write(15,*)'mode4_lund_structure=legacy-four-particle'
            write(15,*)'mode4_emits_lund=1'
         else
            write(15,*)'mode4_event_weight_microbarn=0'
            write(15,*)'mode4_lund_structure=none-calibration'
            write(15,*)'mode4_emits_lund=0'
         endif
      endif
      close(unit=15)

      write(14,1425)reg1,reg2-reg1,reg3-reg2, reg4-reg3,1.-reg4
      write(14,*)'csrng =',csrng
 1425 format(' size of cosine regions:',5f5.2)
      do intreg=1,6
        write(14,1410)intreg,(dismc(intreg,j),j=1,100)
      enddo

 1410 format(' Distribution of mcall values, region', i2 /10(1x,10i7/))
 9001 format(I0,10(',',I0),51(',',ES24.16E3))
 9002 format(I0,',',I0,13(',',ES24.16E3))
 9003 format(I0,2(',',I0),14(',',ES24.16E3))
 9010 format(I0,11(',',I0))

      close(14)
Cvpk      CALL bos_end(recname)

 99   continue

       write(6,*)' END OF LOOP'
       stop
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine mode4_sample_complement(r_u,r_ep,u_gamma,hcos,hphi,
     * ru_low,ru_high,rep_low,rep_high,ugamma_low,ugamma_high,
     * hcos_low,hcos_high,hphi_origin,hphi_relative_low,
     * hphi_relative_high)

c     Draw uniformly from the complement of the five-dimensional guard.
c     Rejection from the unit hypercube is exact.  This operation is used
c     only by fixed-trial calibration; production retains the unrestricted
c     legacy component and its full-support mixture correction.

      implicit none

      real r_u,r_ep,u_gamma,hcos,hphi
      real ru_low,ru_high,rep_low,rep_high
      real ugamma_low,ugamma_high,hcos_low,hcos_high
      real hphi_origin,hphi_relative_low,hphi_relative_high
      real relative_phi,myran
      integer inside

 10   r_u=myran()
      r_ep=myran()
      u_gamma=myran()
      hcos=myran()
      hphi=myran()
      inside=1
      if (r_u .lt. ru_low .or. r_u .gt. ru_high)inside=0
      if (r_ep .lt. rep_low .or. r_ep .gt. rep_high)inside=0
      if (u_gamma .lt. ugamma_low .or.
     *    u_gamma .gt. ugamma_high)inside=0
      if (hcos .lt. hcos_low .or. hcos .gt. hcos_high)inside=0
      relative_phi=modulo(hphi-hphi_origin+.5,1.)-.5
      if (relative_phi .lt. hphi_relative_low .or.
     *    relative_phi .gt. hphi_relative_high)inside=0
      if (inside .eq. 1)go to 10
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine mode4_proposal_ratio(r_u,r_ep,u_gamma,hcos,hphi,
     * core_fraction,core_volume,ru_low,ru_high,rep_low,rep_high,
     * ugamma_low,ugamma_high,hcos_low,hcos_high,hphi_origin,
     * hphi_relative_low,hphi_relative_high,density_ratio,inside)

c     Return q_legacy/q_mixture for
c
c       q_mixture = core_fraction*q_core
c                 + (1-core_fraction)*q_legacy.
c
c     q_core is uniform in one five-dimensional native-coordinate box and
c     zero outside it.  q_legacy is uniform on the complete unit hypercube.
c     The periodic hadron-phi coordinate is compared relative to the frozen
c     recipe origin.

      implicit none

      real r_u,r_ep,u_gamma,hcos,hphi
      real core_fraction,core_volume
      real ru_low,ru_high,rep_low,rep_high
      real ugamma_low,ugamma_high,hcos_low,hcos_high
      real hphi_origin,hphi_relative_low,hphi_relative_high
      real density_ratio,legacy_fraction,relative_phi,denominator
      integer inside

      legacy_fraction=1.-core_fraction
      density_ratio=1./legacy_fraction
      inside=1
      if (r_u .lt. ru_low .or. r_u .gt. ru_high)inside=0
      if (r_ep .lt. rep_low .or. r_ep .gt. rep_high)inside=0
      if (u_gamma .lt. ugamma_low .or.
     *    u_gamma .gt. ugamma_high)inside=0
      if (hcos .lt. hcos_low .or. hcos .gt. hcos_high)inside=0
      relative_phi=modulo(hphi-hphi_origin+.5,1.)-.5
      if (relative_phi .lt. hphi_relative_low .or.
     *    relative_phi .gt. hphi_relative_high)inside=0
      if (inside .eq. 1)then
         denominator=legacy_fraction+core_fraction/core_volume
         density_ratio=1./denominator
      endif
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine mode4_stratum_membership(q2,xb,minus_t,phi,w,y,
     * q2_low,q2_high,xb_low,xb_high,t_low,t_high,phi_low,phi_high,
     * w_minimum,apply_y_max,y_maximum,status)

c     Classify the exact final-LUND analysis coordinates used by mode 4.
c     Bounds are inclusive within roundoff.  The unrestricted proposal and
c     physical checks remain outside this routine.

      implicit none

      real q2,xb,minus_t,phi,w,y
      real q2_low,q2_high,xb_low,xb_high,t_low,t_high
      real phi_low,phi_high,w_minimum,y_maximum,phi_wrapped
      integer apply_y_max,status

      status=1
      phi_wrapped=modulo(phi,360.)
      if (q2 .lt. q2_low .or. q2 .gt. q2_high)return
      if (xb .lt. xb_low .or. xb .gt. xb_high)return
      if (minus_t .lt. t_low .or. minus_t .gt. t_high)return
      if (phi_wrapped .lt. phi_low .or.
     *    phi_wrapped .gt. phi_high)return
      if (w .lt. w_minimum)return
      if (apply_y_max .eq. 1 .and. y .gt. y_maximum)return
      status=0
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine radiative_t_mapping(es,ep,th0,ephot,cstk,mpi,mp,
     * analysis_mp,t_at_cosine_zero,t_jacobian,status)

c     At fixed lepton and internal-photon kinematics, -t is affine in the
c     pion center-of-mass cosine:
c
c       -t(cos(theta*)) = t_at_cosine_zero - t_jacobian*cos(theta*).
c
c     This lets the balanced survey sample an analysis -t bin directly while
c     retaining the exact legacy cosine density in the proposal correction.

      implicit none

      real*8 es,ep,ephot,mpi,mp,analysis_mp
      real th0,cstk,t_at_cosine_zero,t_jacobian
      real*8 mel,ps,pp,nu,q2,qvec,pwr2,pwr
      real*8 ewreal,wreal2,wreal,beta,gamma,pstar2,pstar,epcm
      integer status

      status=1
      t_at_cosine_zero=0.
      t_jacobian=0.
      mel=0.511D-3
      if (es .le. mel .or. ep .le. mel) return
      ps=sqrt(es**2-mel**2)
      pp=sqrt(ep**2-mel**2)
      nu=es-ep
      q2=2.D0*es*ep-2.D0*ps*pp*dcos(dble(th0))-2.D0*mel**2
      if (q2 .le. 0.D0) return
      qvec=sqrt(q2+nu**2)
      ewreal=nu+mp-ephot
      pwr2=qvec**2+ephot**2-2.D0*qvec*ephot*dble(cstk)
      if (pwr2 .le. 0.D0 .or. ewreal .le. 0.D0) return
      pwr=sqrt(pwr2)
      wreal2=ewreal**2-pwr2
      if (wreal2 .le. (mp+mpi)**2) return
      wreal=sqrt(wreal2)
      beta=pwr/ewreal
      gamma=ewreal/wreal
      if (beta .le. 0.D0 .or. beta .ge. 1.D0) return
      pstar2=((wreal**2-mp**2-mpi**2)**2/4.D0-
     *   (mp*mpi)**2)/wreal**2
      if (pstar2 .le. 0.D0) return
      pstar=sqrt(pstar2)
      epcm=sqrt(pstar**2+mp**2)
      t_at_cosine_zero=real(2.D0*analysis_mp*gamma*epcm-
     *   analysis_mp**2-mp**2)
      t_jacobian=real(2.D0*analysis_mp*gamma*beta*pstar)
      if (t_jacobian .le. 0.) return
      status=0
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine balanced_proposal_ratio(q2,xb,minus_t,phicm,
     * t_jacobian,uq2rng,eprng,analysis_mp,legacy_fraction,
     * nq2,q2_edges,nxb,xb_edges,nt,t_edges,nphi,phi_edges,
     * balanced_to_legacy,density_ratio)

c     Return q_legacy/q_mixture for
c
c       q_mixture = legacy_fraction*q_legacy
c                 + (1-legacy_fraction)*q_balanced.
c
c     q_balanced chooses analysis Q2, xB, -t, and phi bins uniformly,
c     then samples uniformly inside the chosen bin.  The xB-to-E' and
c     -t-to-cos(theta*) Jacobians are included explicitly.  A point outside
c     the configured balanced support has q_balanced=0 and therefore receives
c     density_ratio=1/legacy_fraction.

      implicit none

      integer nq2,nxb,nt,nphi
      integer iq2,ixb,it,iphi,i
      real q2,xb,minus_t,phicm,t_jacobian,uq2rng,eprng
      real legacy_fraction,q2_edges(*),xb_edges(*),t_edges(*)
      real phi_edges(*),balanced_to_legacy,density_ratio
      real dq,dx,dt,dphi,phi_wrapped
      real*8 analysis_mp

      iq2=0
      ixb=0
      it=0
      iphi=0
      balanced_to_legacy=0.
      density_ratio=1./legacy_fraction
      if (q2 .le. 0. .or. xb .le. 0. .or.
     *    t_jacobian .le. 0.) return

      do i=1,nq2
         if (q2 .ge. q2_edges(i) .and.
     *       q2 .lt. q2_edges(i+1)) iq2=i
      enddo
      do i=1,nxb
         if (xb .ge. xb_edges(i) .and.
     *       xb .lt. xb_edges(i+1)) ixb=i
      enddo
      do i=1,nt
         if (minus_t .ge. t_edges(i) .and.
     *       minus_t .lt. t_edges(i+1)) it=i
      enddo
      phi_wrapped=modulo(phicm,360.)
      do i=1,nphi
         if (phi_wrapped .ge. phi_edges(i) .and.
     *       phi_wrapped .lt. phi_edges(i+1)) iphi=i
      enddo
      if (iq2 .eq. 0 .or. ixb .eq. 0 .or.
     *    it .eq. 0 .or. iphi .eq. 0) return

      dq=q2_edges(iq2+1)-q2_edges(iq2)
      dx=xb_edges(ixb+1)-xb_edges(ixb)
      dt=t_edges(it+1)-t_edges(it)
      dphi=phi_edges(iphi+1)-phi_edges(iphi)
      balanced_to_legacy=
     *   uq2rng*q2**2/(real(nq2)*dq)*
     *   eprng*2.*real(analysis_mp)*xb**2/
     *      (q2*real(nxb)*dx)*
     *   2.*t_jacobian/(real(nt)*dt)*
     *   360./(real(nphi)*dphi)
      density_ratio=1./(legacy_fraction+
     *   (1.-legacy_fraction)*balanced_to_legacy)
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine build_final_candidate(ebeam,es,ep_pre,ep_out,ep_min,
     * th0,theta,ephot,cstk,phik,mpi_s,csthcm,phicm,targp,delta,
     * ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     * wreal,mm2,epw,phir,outgoing_loss_base,
     * outgoing_loss_accept_test,status)

c     Construct the exact post-outgoing-external-loss candidate used for
c     legacy LUND output.  Survey mode calls the same routine before observed
c     coordinate classification so the candidate definitions cannot drift.
c
c     status=0: valid final candidate
c     status=1: final electron below the configured energy threshold
c     status=2: no physical final hadronic candidate

      implicit none

#include "mc.inc"

      COMMON/ALPHA/ ALPHA,PI,MP,MN,MPI,MEL,WG,EPIREA,TH_OPT,RES_OPT

      real*8 alpha,pi,mp,mn,mpi,mel,wg
      real*8 es,ep_out,ephot,delta
      real ebeam,ep_pre,ep_min,th0,theta,cstk,phik,mpi_s
      real csthcm,phicm,targp
      real ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi
      real ekx,eky,ekz,wreal,mm2,epw,phir
      real xs,eloss,gxs,xtest
      real outgoing_loss_base,outgoing_loss_accept_test
      real w2,rotc,rots,px,py
      real myran
      integer epirea,th_opt,res_opt,status

      status=0

 222  xs=myran()
      eloss=xs**(1./targp)
      gxs=1.-eloss
      xtest=myran()
      if (xtest .gt. gxs)go to 222
      outgoing_loss_base=xs
      outgoing_loss_accept_test=xtest

      ep_out=ep_pre*(1.-eloss)
      if (ep_out .lt. ep_min)then
         status=1
         return
      endif

      call missm(epirea,ebeam,es,ep_out,th0,ephot,cstk,phik,mpi_s,
     * ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     * csthcm,phicm,wreal,mm2)

      if (mm2 .eq. 0. .or. ep_out .lt. ep_min)then
         status=2
         return
      endif

      w2=mp**2+2*mp*(ebeam-ep_out)
     *   -2*es*ep_out*(1-cos(th0))+2.*mel**2
      epw=sqrt(w2)

c     Rotate all momenta by one common random angle around the beam line.

      phir=2.*pi*myran()
      rotc=cos(phir)
      rots=sin(phir)

      px=ep_out*sin(theta*pi/180.)
      py=0.
      p(1,1)=px*rotc+py*rots
      p(1,2)=py*rotc-px*rots
      p(1,3)=ep_out*cos(theta*pi/180.)
      p(1,4)=ep_out

      if ((epirea .eq. 1).or.(epirea.eq.5)) then
         p(2,1)=ppx*rotc+ppy*rots
         p(2,2)=ppy*rotc-ppx*rots
         p(2,3)=ppz
         p(2,4)=eprot

         if (npart .eq. 4)then
            p(3,1)=ppix*rotc+ppiy*rots
            p(3,2)=ppiy*rotc-ppix*rots
            p(3,3)=ppiz
            p(3,4)=epi
         endif
      else
         p(2,1)=ppix*rotc+ppiy*rots
         p(2,2)=ppiy*rotc-ppix*rots
         p(2,3)=ppiz
         p(2,4)=epi

         if (npart .eq. 4)then
            p(3,1)=ppx*rotc+ppy*rots
            p(3,2)=ppy*rotc-ppx*rots
            p(3,3)=ppz
            p(3,4)=eprot
         endif
      endif

      if (npart .eq. 4)then
         if (ephot .le. delta)then
            p(4,1)=0.
            p(4,2)=0.
            p(4,3)=1.e-5
            p(4,4)=1.e-5
         else
            p(4,1)=ekx*rotc+eky*rots
            p(4,2)=eky*rotc-ekx*rots
            p(4,3)=ekz
            p(4,4)=ephot
         endif
      endif

      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine observed_analysis_coordinates(ebeam,analysis_mp,
     * q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs,status)

c     Reproduce the analysis coordinate convention from final generator
c     particles: nominal beam plus the post-external-loss LUND electron and
c     proton.  The dedicated analysis proton mass is passed explicitly.

      implicit none

#include "mc.inc"

      real ebeam,q2_obs,xb_obs,minus_t_obs,phi_obs,w_obs,y_obs
      real*8 analysis_mp
      real*8 ex,ey,ez,ee,pxp,pyp,pzp,epr
      real*8 qx,qy,qz,q0,qnorm2,qnorm
      real*8 lnx,lny,lnz,hnx,hny,hnz,lnorm,hnorm
      real*8 cx,cy,cz,cosphi,sinphi,w2_obs,pi_local
      integer status

      status=1
      q2_obs=0.
      xb_obs=0.
      minus_t_obs=0.
      phi_obs=0.
      w_obs=0.
      y_obs=0.

      ex=dble(p(1,1))
      ey=dble(p(1,2))
      ez=dble(p(1,3))
      ee=dble(p(1,4))
      pxp=dble(p(2,1))
      pyp=dble(p(2,2))
      pzp=dble(p(2,3))
      epr=dble(p(2,4))

      q0=dble(ebeam)-ee
      qx=-ex
      qy=-ey
      qz=dble(ebeam)-ez
      qnorm2=qx*qx+qy*qy+qz*qz
      q2_obs=real(qnorm2-q0*q0)
      if (q0 .le. 0.D0 .or. q2_obs .le. 0.)return

      xb_obs=real(dble(q2_obs)/(2.D0*analysis_mp*q0))
      y_obs=real(q0/dble(ebeam))
      w2_obs=analysis_mp**2+2.D0*analysis_mp*q0-dble(q2_obs)
      if (w2_obs .le. 0.D0)return
      w_obs=real(sqrt(w2_obs))

      minus_t_obs=real(pxp*pxp+pyp*pyp+pzp*pzp
     *   -(analysis_mp-epr)**2)

c     Lepton normal = beam x electron.

      lnx=-dble(ebeam)*ey
      lny=dble(ebeam)*ex
      lnz=0.D0

c     Hadron normal = proton x q.

      hnx=pyp*qz-pzp*qy
      hny=pzp*qx-pxp*qz
      hnz=pxp*qy-pyp*qx
      lnorm=sqrt(lnx*lnx+lny*lny+lnz*lnz)
      hnorm=sqrt(hnx*hnx+hny*hny+hnz*hnz)
      qnorm=sqrt(qnorm2)
      if (lnorm .le. 0.D0 .or. hnorm .le. 0.D0 .or.
     *    qnorm .le. 0.D0)return

      cosphi=(lnx*hnx+lny*hny+lnz*hnz)/(lnorm*hnorm)
      cx=lny*hnz-lnz*hny
      cy=lnz*hnx-lnx*hnz
      cz=lnx*hny-lny*hnx
      sinphi=(qx*cx+qy*cy+qz*cz)/(qnorm*lnorm*hnorm)
      pi_local=acos(-1.D0)
      phi_obs=real(atan2(sinphi,cosphi)*180.D0/pi_local)
      if (phi_obs .lt. 0.)phi_obs=phi_obs+360.

      status=0
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real function sigma(ek,Tk,epcos,epphi,ehel)

c     Calculate the single pion electroproduction cross section with
c      radiative tail, according to the prescription in Mo and Tsai.
c     Mo and Tsai  calculate dsigma/d_omega dp for (omega>delta)
c     omega is the energy of the radiated photon.

c     The Mo and Tsai cross section for the 3-3 resonance is replaced with
c     the AO cross section for single pion production from the proton.


      implicit none
      COMMON/ALPHA/ ALPHA,PI,MP,MN,MPI,MEL,WG,EPIREA,TH_OPT,RES_OPT
      common /radcal/T0,es,ep,ps,pp,rs,rp,u0,pu,uu,cst0,snt0,csths,csthp
     * ,snths,snthp,pdotk,sdotk

      real*8 ek,Tk
      real*8 alpha,pi,mp,mn,mpi,mel,wg
      real*8 T0,es,ep,ps,pp,rs,rp,u0,pu,uu,pdotk,sdotk
      real*8 cst0,snt0,csths,csthp,snths,snthp
      real*8 csthk,snthk,qq,mf2
      real*8 sp
      real*4 ffac1,ffac2,ffac3,ffac4,ffac5,ffac6,ffac
      real*4 gfac1,gfac2,gfac3,gfac4,gfac

      real*4 f,g,fkt
      real*4 sig_r,sigf
      real nu
c     arguments for aaosub_1:
      real*4 qsq,epw,th0,sigma0,sigu,sigt,sigl,sigi
      real*4 sigip, asym_p
      real epcos,epphi
      real*4 q0,kfac
      real s2,epeps
      integer epirea, th_opt, res_opt
      integer ehel
c
      th0=T0
      csthk=cos(Tk)
      snthk=sin(Tk)
      qq=2*mel**2-2*es*ep+2*ps*pp*cst0-2*ek*(es-ep)+2*ek*pu*csthk
      mf2=uu-2*ek*(u0-pu*csthk)
      if (mf2 .lt. wg**2 .or. qq .ge. 0.)then
         sigma=0.
         return
      endif

      epw=sqrt(mf2)

      sp=es*ep-ps*pp*cst0

      ffac1=-(mel/pdotk)**2*(2.*es*(ep+ek)+qq/2)
      ffac2=-(mel/sdotk)**2*(2.*ep*(es-ek)+qq/2)
      ffac3=-2.
      ffac4=2/sdotk/pdotk*(mel**2*(sp-ek**2)+sp*(2*es*ep-sp+ek*(es-ep)))       
      ffac5=(2*(es*ep+es*ek+ep*ep)+qq/2-sp-mel**2)/pdotk
      ffac6=-(2*(es*ep-ep*ek+es*es)+qq/2-sp-mel**2)/sdotk
      ffac=ffac1+ffac2+ffac3+ffac4+ffac5+ffac6
      if (ffac .le. 0.)then
         sigma=0.1e-30
         return
      endif

      gfac1=mel**2*(2*mel**2+qq)* (1./(pdotk**2)+1./(sdotk**2) )
      gfac2=4.
      gfac3=4.*sp*(sp-2*mel**2)/pdotk/sdotk
      gfac4=(2*sp+2*mel**2-qq)*(1./pdotk-1./sdotk)

      gfac=gfac1+gfac2+gfac3+gfac4
      if (gfac .le. 0.)then
         sigma=0.1e-30
         return
      endif
      qsq=-qq
      q0=es-ep

c     Use the AO cross section to fake up the values of f and g
c     needed for the Mo and Tsai integration.
c     This should be compared to the formulas given in Mo and Tsai
c     for the 3-3 resonance.  We ought to agree pretty well in the
c     3-3 resonance region.

     
      call dsigma(th0,qsq,epw,epcos,epphi,th_opt,epirea,res_opt,sigma0
     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
c      print *, 'AO:',sigma0,sigu,sigt,sigl,sigi
c      call dsigma(th0,qsq,epw,epcos,epphi,4,epirea,res_opt,sigma0
c     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
c      print *, 'MAID:',sigma0,sigu,sigt,sigl,sigi
c      print *, ''
     
      if (sigma0 .le. 0.)then
         sigma = 0.0
         return
      endif

c     According to my calculations
c     F=(2K/Mp**3)[q2/(q2+nu**2)][(sigt+sigl)/sig0]
c     G=(2K/Mp)(sigt/sig0)
c     where K=(W**2-Mp**2)/(2Mp)
c     sig0=4*pi**2*alpha/mp**2=127 microbarns
C     In terms of the sigt and sigl, the cross section,
c     differential in omega_e and E_e, is
c     sigma=[alpha/(4 pi**2][2 * K * E_s]/[q2 * E_s * (1-epsilon)]
c     *  (sigt+epsilon*sigl)
c     where
c     1/epsilon = 1+2tan(theta_e/2)**2(1+nu**2/q2)

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     The AO program calculates sigu,sigl,sigt and sigi with
c     the kinematic factors included. To convert them to the
c     form needed for the above problem we must use
      s2=(1.-cst0)/2.
      epeps=1+2*(1+q0**2/qsq)*s2/(1-s2)

      if (epeps .gt. 1.)go to 120
        write(6,*)' sigma-1: epeps-inv =',epeps
        sigma=0.1e-30
        return
 120  epeps=1./epeps

      Kfac=(mf2-mp**2)/2./mp
c      cfac=2*pi**2*qsq*es*(1-epeps)/kfac/ep/alpha
c      sigu=sigu*cfac
c      sigl=sigl*cfac
c      sigt=sigt*cfac
c      sigi=sigi*cfac
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c      sig0=4.*pi**2*alpha/mp**2
c     The following would work if there were no polarization terms
c      f=2.*kfac/(mp**3)*qsq/(qsq+q0**2)*(sigu+sigl)/sig0
c      g=(2.*kfac/mp)*(sigu)/sig0
c     I don't know what to do with the polarization terms: these give
c     the dependence on the center of mass phi decay angle.  I try
c     the following and modulate the cross section later:
      nu=(mf2-mp**2+qsq)/2./mp
      if (nu .gt. 0.)go to 121
       sigma=0.1e-30
       write(6,*)' sigma-2: nu = ',nu
       return

 121  f=kfac/(2.*pi**2*alpha*mp)/(1.+nu**2/qsq)*(sigu+sigl)
      g=mp/(2.*pi**2*alpha)*kfac*sigu
        fkt=(epw**2-mp**2+mpi**2)/2./epw
        fkt=sqrt(fkt**2-mpi**2)*2.*epw/(epw**2-mp**2)
        f=f*fkt
        g=g*fkt

c      f_old=2.*kfac/(mp**3)*qsq/(qsq+q0**2)*(sigu+sigl)/sig0
c      g_old=(2.*kfac/mp)*(sigu)/sig0

c     ??????????????????????????
c     The following formula is the same as B.8 in Mo and Tsai,
c     except I have divided by 2pi.  It seems to me that the
c     result doesn't make sense otherwise.
c     ?????????????????????????


      sig_r=((alpha**3/(2*pi*qq)**2)/mp)*(ep/es)*ek
      sigf=mp**2*f*ffac+g*gfac
      if (sigf .gt. 0.)go to 122
c        print *, 'SIGMA:',th0,qsq,epw,epcos,epphi,theory_opt
c        write(6,*)' sigma-3: ffac, gfac, f, g,sigu,sigl, sigf ='
c     *  ,ffac,gfac,f,g,sigu,sigl,sigf       
        sigma=0.1e-30
        return

c     Modulate the cross section with the polarization and interference
c     terms
 122   if (sigu+epeps*sigl .gt. 0.)go to 123
         write(6,*)' sigma-4:, sigu+epeps*sigl =',sigu+epeps*sigl
         sigma=0.1e-30
         return

 123   sigf=sigf * (1.+(epeps*sigt*cos(epphi*pi/90.)
     &  +sqrt(epeps*(1.+epeps)/2)*sigi*cos(epphi*pi/180.)
     &  +ehel*sqrt(epeps*(1.-epeps)/2)*sigip*sin(epphi*pi/180.))
c     &  /(sigu+epeps*sigl))
     &  /sigf)

       if (sigf .gt. 0.)go to 124
         write(6,*)' sigma-5: sigf =',sigf
         
         sigma=0.1e-30
         return

 124   sigma=sig_r*sigf
       
      return
      end


      function myran()
c     Random number generator used because I can't find one in the
c     library.

      implicit none
      real myran

      call random_number(myran)

      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine missm(type,ebeam,es,ep,th0,ephot,cstk,phik,mpi,ppx,
     * ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ephotx,ephoty,ephotz,
     * csthcm,phicm,wreal,mm2)

c
c     Input:
c         type= 1- pi0 , 3 - p+, 5 -eta
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

      real*8 es,ep,ps,pp,ephot

      real beta
      real cstk
      real csthe
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
      real qx,qz
      real q_dot_pp
      real snphk
      real snphi
      real snthcm
      real snthe
      real sntk
      real sntq
      real th0
      real w2
      real wmin
      real wreal
      integer type
      integer*4 idum
c

      mp=.938
      mel=0.511E-3
      pi=3.14159
      wmin=mp+mpi
      csthe=cos(th0)
      snthe=sin(th0)
      nu=es-ep
      ps=abs(es**2-mel**2)
      pp=abs(ep**2-mel**2)
      ps=sqrt(ps)
      pp=sqrt(pp)
      q2=2.*es*ep-2.*ps*pp*csthe-2.*mel**2
      w2=mp**2-q2+2.*mp*nu
c     get components of the q vector
      qx=-pp*snthe
      qz=ps-pp*csthe
      qvec=sqrt(qx**2+qz**2)
      q2=2.*es*ep-2.*ps*pp*csthe-2.*mel**2
      w2=mp**2-q2+2.*mp*nu
c     get components of the q vector
      qx=-pp*snthe
      qz=ps-pp*csthe
      qvec=sqrt(qx**2+qz**2)
c     get components of the photon vector
      if (abs(cstk) .gt. 1.)then
         write(6,*)' missm-1: cstk =',cstk
         cstk=cstk/abs(cstk)
      endif
      sntk=sqrt(1.-cstk**2)
      csphk=cos(phik)
      snphk=sin(phik)

      cstq=qz/qvec
      sntq=sqrt(1.-cstq**2)
      ephotx=ephot*(sntk*csphk*cstq-cstk*sntq)
      ephoty=ephot*sntk*snphk
      ephotz=ephot*(cstk*cstq+sntk*csphk*sntq)

c     calculate the dot product of the photon vector and the q-vector
      eph_dot_q=ephotx*qx+ephotz*qz
c     calculate the mass of the actual hadronic system for the two
c     photon directions.

      wreal=(w2-2.*ephot*(nu+mp)+2.*eph_dot_q)
      if (wreal .le. wmin**2)go to 30
      wreal=sqrt(wreal)
c     go to the end of the loop if the hadronic mass is below the pion
c     threshold.

c     calculate the energy of the actual hadronic system, with radiation
      ewreal=nu+mp-ephot


c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Calculate the momentum components of the nucleon and pion in the
c     Lab frame

      snthcm=(1.-csthcm**2)
      if (snthcm .le. 0.)then
         write(6,*)' missm-2: snthcm =',snthcm
         snthcm=0.0000001
      endif
      snthcm=sqrt(snthcm)
      

c     Calculate cos and sin of phicm
      csphi=cos(phicm*pi/180.)
      snphi=sin(phicm*pi/180.)

c     calculate laboratory components of the resonance momentum vector
         pwrz=qz-ephotz
         pwrx=qx-ephotx
         pwry=-ephoty
         pwr=sqrt(pwrz**2+pwrx**2+pwry**2)
c     calculate the relativistic factors, gamma and beta, for the resonance
         beta=pwr/ewreal
         gamma=ewreal/wreal

c     define angle cosines for the laboratory system with the resonance
c     as the z axis.  Choose the y axis of the resonance frame perpendicular
c     to the laboratory x axis.
         pfac=sqrt(pwry**2+pwrz**2)
         cxx=pfac/pwr
         cxy=-pwrx*pwry/pfac/pwr
         cxz=-pwrx*pwrz/pfac/pwr
         cyx=0
         cyy=pwrz/pfac
         cyz=-pwry/pfac
         czx=pwrx/pwr
         czy=pwry/pwr
         czz=pwrz/pwr

c     calculate the momentum of the pion and proton in the resonance frame
         pstar=( (wreal**2-mp**2-mpi**2)**2/4.-(mp*mpi)**2 )/wreal**2
         if (pstar .le. 0.)then
            write(6,*)'missm-3: pstar=',pstar
            pstar=0.000001
         endif
         pstar=sqrt(pstar)
c     Calculate the energy of the proton and pion in the resonance center
c     of mass frame.
         epcm  = sqrt(pstar**2+mp**2)
         epicm = sqrt(pstar**2+mpi**2)

c     Calculate the pion momentum components and energy in the lab frame
c     where the z axis is the direction of the momentum of the resonance.
c     The center of mass angles are pion angles.
         ppiwx=pstar*snthcm*csphi
         ppiwy=pstar*snthcm*snphi
         ppiwz=gamma*(pstar*csthcm+beta*epicm)
         epi=gamma*(epicm+beta*pstar*csthcm)

c         print *,'M(eta)=',sqrt(epi**2-ppiwz**2-ppiwy**2-ppiwx**2)
c     Calculate momentum components of the pion in
c     the lab frame where the z axis is along the incident beam
         ppix=ppiwx*cxx+ppiwy*cyx+ppiwz*czx
         ppiy=ppiwx*cxy+ppiwy*cyy+ppiwz*czy
         ppiz=ppiwx*cxz+ppiwy*cyz+ppiwz*czz

c     Calculate the proton energy and momentum components
c     in the lab resonance frame.
c     The center of mass angles are pion angles.
         ppwx=-ppiwx
         ppwy=-ppiwy

c         ppwz=-ppiwz+gamma*beta*(epicm+epcm)
c     epicm+epcm=wreal
         ppwz=gamma*beta*wreal-ppiwz

         eprot=gamma*wreal-epi

c     Rotate the lab momentum components of the proton into the frame
c     where the z axis is along the incident electron beam.
         ppx=ppwx*cxx+ppwy*cyx+ppwz*czx
         ppy=ppwx*cxy+ppwy*cyy+ppwz*czy
         ppz=ppwx*cxz+ppwy*cyz+ppwz*czz


c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c     Calculate the square of the missing mass, associated with the
c     momentum components of the charged hadron:
      q2=2.*ebeam*ep*(1.-csthe)
      qz=ebeam-pp*csthe
      qvec=sqrt(qx**2+qz**2)
      nu=ebeam-ep
         if ((type .eq. 1).or.(type.eq.5)) then
         q_dot_pp=qx*ppx+qz*ppz
         mm2=-q2+2*mp**2+2*mp*(nu-eprot)-2*nu*eprot+2*q_dot_pp
         else       
         q_dot_pp=qx*ppix+qz*ppiz
         mm2=-q2+mp**2+mpi**2+2*mp*(nu-epi)-2*nu*epi+2*q_dot_pp
         endif

         return
 30      mm2=0.
         return
         end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real function spence(x)
c     Calculate the Spence function needed for the Mo and Tsai formula.

      implicit none
      real x
      real pi
      real sintp,sintn
      pi=3.14159

      if (abs(x) .lt. 0.1)then
         spence=x+x**2/4.
c         write(6,*)' spence: abs(x) .lt. 0.1'
         return
      endif

      if (x .gt. 0.99 .and. x .lt. 1.01)then
         spence=pi**2/6.
c         write(6,*)' spence: x=1.'
         return
      endif

      if (x .gt. -1.01 .and. x .lt. -0.99)then
         spence=-pi**2/12.
c         write(6,*)' spence: x= -1.'
         return
      endif

      if (x .gt. 0.)then
         spence=.1025+sintp(x)
c         write(6,*)' x .gt. 0.'
         return
      endif
      spence=-0.0975+sintn(x)
c      write(6,*)' spence: x .lt. 0.'
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real function sintp(x)
      implicit none
      real x
      real xstep,sum,y,arg
      integer i

      xstep=(x-.1)/100.
      sum=0.
      y=.1-xstep/2.
      do i=1,100
        y=y+xstep
        arg=abs(1.-y)
        sum=sum-alog(arg)/y
      enddo
      sintp=sum*xstep
      return
      end

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      real function sintn(x)
      implicit none
      real x,xa,ystep,y,sum
      integer i

      xa=abs(x)
      ystep=(xa-0.1)/100.
      sum=0.
      y=.1-ystep/2.
      do i=1,100
        y=y+ystep
        sum=sum-alog(1.+y)/y
      enddo
      sintn=sum*ystep
      return
      end

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine gauss(x,y,sigma_x,sigma_y)

c     calculate two random numbers, x, y,  for gaussian distributions
c     with s.d. of sigma_x and sigma_y.

      implicit none

      common /random/idum

      real x,y,sigma_x,sigma_y
      real r1,r2,pi
      real myran

      integer*4 idum

         pi=3.14159
         r1=myran()
         r2=myran()
         r1=sqrt(-2.*alog(r1))
         r2=2.*pi*r2
         x=sigma_x*r1*cos(r2)
         y=sigma_y*r1*sin(r2)
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
	  integer iseed
      real    random, myran 

	  random = myran()
      random = 0.5 - myran()
      get_spin = 1
      if(random.lt.0) get_spin = -1

	  return
	  end
C======================================================================
