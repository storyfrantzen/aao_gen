      subroutine xsection

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'

      real vl, vt, vlt, vtt, vltp
      real phicm_rad, ekin
      real sigma_p, sigma_u
      real ratio

      call multipole_amps ! calc. multipole amplitudes
      call helicity_amps  ! calc. helicity amplitudes

      phicm_rad = phicm*pi/180.0
      
      vt   = 1.0
      vl   = epsilon
      vtt  = epsilon
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     vlt and vltp are divided by 2 compared with spp_int_e1
c     to match AO formalism
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      vlt  = sqrt(epsilon*(1+epsilon)/2)
      vltp = sqrt(epsilon*(1-epsilon)/2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ekin = sqrt(Q2) / nu_cm
      fkt  = 2*W*ppi_mag_cm/(W**2 - m_p**2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5 Response fucntions     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      sigma_t  = (cabs(hh1)**2+cabs(hh2)**2
     1           +cabs(hh3)**2+cabs(hh4)**2)/2.
      sigma_l  =  cabs(hh5)**2+cabs(hh6)**2
      sigma_tt =  real(hh3*conjg(hh2)-hh4*conjg(hh1))
      sigma_lt =  sqrt(2.0)*real(conjg(hh5)*(hh1-hh4)  + 
     1                           conjg(hh6)*(hh2+hh3))     
      sigma_ltp = sqrt(2.0)*aimag(conjg(hh5)*(hh4-hh1) - 
     1                  conjg(hh6)*(hh2+hh3))
c
      sigma_l  = sigma_l*ekin**2
      sigma_lt = sigma_lt*ekin
      sigma_ltp= sigma_ltp*ekin
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    sigma_u: unpolarized cross section
c    simga_h: polarized cross section
c    asym_ltp: single spin beam asymmetry
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      sigma_u  = fkt*(vt*sigma_t   
     1        +  vl*sigma_l  
     1        + vtt*sigma_tt*cos(2*phicm_rad)
     1        + vlt*sigma_lt*cos(phicm_rad))
c
c      if (q2.ge.1) sigma_u = sigma_u*(0.0918*q2**2.1485)
c      if (q2.ge.1) sigma_p = sigma_p*(0.0918*q2**2.1485)
       sigma_u=sigma_u * ratio(q2)
       sigma_p=sigma_p * ratio(q2)

      sigma_p  = fkt*vltp*sigma_ltp*sin(phicm_rad)

      sigma_0  = sigma_u + e_hel*sigma_p
c
      asym_ltp = sigma_p/sigma_u
c      print *, 'xsection: ',sigma_0,sigma_u, sigma_p,asym_ltp, e_hel
c
      end
C
C------     RATIO(Q2) =  DVMP(W=1.8)/MAID(W=1.8) as a function of Q2
C
      REAL FUNCTION RATIO(Q2)
      implicit none
      real Q2,x(19),y(19)
      INTEGER I,N
      data x/ 0.,0.1387, 0.243, 0.3446, 0.4455,
     *    0.669, 1.193, 1.706, 2.,   2.21, 2.71, 3.4,
     *    4.,    4.41,  5.45,  6.46, 7.46, 8.47, 9.45/


      data y/7., 4.25, 1.6538, 0.990, 0.8494,
     *       0.493151, 0.368132, 0.444444, 0.454545,
     *       0.514286, 0.692308, 1.21429, 1.81416,
     *       2.55738, 4.0, 4.15385, 4.84444, 5.75, 7.04545/
      N=19
      do I=1,N-1
        if(Q2.ge.x(i) .and. Q2.le.x(i+1)) then
          ratio=y(i)+(y(i+1)-y(i))/(x(i+1)-x(i))*(Q2-x(i))
        endif
      enddo
      if(Q2.gt.x(N)) ratio=Y(N) + (Q2-X(N)) * 1.0538
      RETURN
      END

