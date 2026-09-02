! Shared-model regression probe, linked separately against Born and rad sources.
program test_dvmp_boundary
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  real, external :: tminq, xsigma_t, xsigma_tt, xsigma_lt
  logical, external :: xcheck_kine, check_kine
  real :: q2, xb, t, beam, w, cost, jacg, jacr, sigma, st, sl, stt, slt, sltp
  real :: values(6), mass, nan
  real, parameter :: mp = 0.93827
  real, parameter :: qs(5) = [1.0, 2.1358285, 4.64465284, 8.0, 10.5]
  real, parameter :: xs(5) = [0.05, 0.2, 0.399552494, 0.586575747, 0.70]
  real, parameter :: beams(2) = [6.535, 10.604]
  real, parameter :: masses(2) = [0.134976, 0.547300]
  real, parameter :: costs(7) = [-0.99999994, -0.9, 0.0, 0.9, 0.99, 0.9999, 0.99999994]
  double precision :: lo, hi, lo_ref, hi_ref, w2, wd, e1, e3, p1, p3
  integer :: im, iq, ix, ib, ic, checked
  logical :: valid, invalid_flag

  call xsinit(masses(1))
  ! These exact REAL values previously made sqrt(T-tminq) invalid.
  q2=4.64465284; xb=0.586575747; t=-0.609785080; beam=6.535
  call ieee_set_flag(ieee_invalid, .false.)
  slt=xsigma_lt(t,xb,q2,beam)
  call require(ieee_is_finite(slt) .and. slt > 0.0, 'RGK failing point')
  q2=2.13582850; xb=0.399552494; t=-0.204734638; beam=10.604
  slt=xsigma_lt(t,xb,q2,beam)
  call require(ieee_is_finite(slt), 'RGA failing point')
  call ieee_get_flag(ieee_invalid, invalid_flag)
  call require(.not.invalid_flag, 'no invalid arithmetic at reported points')

  checked=0
  do im=1,size(masses)
    mass=masses(im)
    call xsinit(mass)
    do iq=1,size(qs)
      q2=qs(iq)
      do ix=1,size(xs)
        xb=xs(ix)
        call dvmp_t_limits(q2,xb,lo,hi,valid)
        w2=dble(q2)*(1d0/dble(xb)-1d0)+dble(mp)**2
        if (w2 <= (dble(mp)+dble(mass))**2) then
          call require(.not.valid, 'reject below hadronic threshold')
          cycle
        endif
        call require(valid .and. lo >= 0d0 .and. hi > lo, 'ordered endpoints')
        ! Independent double-precision CM four-vector expression.
        wd=sqrt(w2)
        e1=(w2+dble(q2)+dble(mp)**2)/(2d0*wd)
        e3=(w2-dble(mass)**2+dble(mp)**2)/(2d0*wd)
        p1=sqrt(e1**2-dble(mp)**2); p3=sqrt(e3**2-dble(mp)**2)
        lo_ref=2d0*(e1*e3-dble(mp)**2-p1*p3)
        hi_ref=2d0*(e1*e3-dble(mp)**2+p1*p3)
        call require(abs(lo-lo_ref) < 1d-10*max(1d0,lo), 'analytic forward limit')
        call require(abs(hi-hi_ref) < 1d-10*max(1d0,hi), 'analytic backward limit')
        call require(tminq(q2,xb) == real(lo), 'shared rounded endpoint')
        do ib=1,size(beams)
          beam=beams(ib)
          t=-real(lo)
          if (.not.xcheck_kine(t,xb,q2,beam)) cycle
          call require(xsigma_lt(t,xb,q2,beam) == 0.0, 'LT vanishes at forward endpoint')
          call require(xsigma_tt(t,xb,q2,beam) == 0.0, 'TT vanishes at forward endpoint')
          t=nearest(t,1.0)
          call require(.not.xcheck_kine(t,xb,q2,beam), 'reject one REAL step below -t minimum')
          call require(xsigma_lt(t,xb,q2,beam) == 0.0, 'outside point has zero cross section')
          t=nearest(-real(hi),-1.0)
          call require(.not.xcheck_kine(t,xb,q2,beam), 'reject one REAL step above -t maximum')
          w=real(wd)
          do ic=1,size(costs)
            cost=costs(ic)
            if (.not.check_kine(cost,w,q2,beam,t,xb,jacg,jacr)) cycle
            call require(xcheck_kine(t,xb,q2,beam), 'angular/direct domain agreement')
            call require(-t >= tminq(q2,xb), 'nonnegative LT radicand')
            call require(jacg > 0.0 .and. jacr > 0.0, 'positive angular Jacobians')
            call ieee_set_flag(ieee_invalid, .false.)
            call dvmpw(cost,w,q2,1.2,beam,1,mass,sigma,st,sl,stt,slt,sltp)
            values=[sigma,st,sl,stt,slt,sltp]
            call require(all(ieee_is_finite(values)), 'finite angular structure functions')
            call require(ieee_is_finite(xsigma_t(t,xb,q2,beam)), 'finite direct transverse term')
            call require(ieee_is_finite(xsigma_tt(t,xb,q2,beam)), 'finite direct TT term')
            call require(ieee_is_finite(xsigma_lt(t,xb,q2,beam)), 'finite direct LT term')
            call ieee_get_flag(ieee_invalid, invalid_flag)
            call require(.not.invalid_flag, 'no invalid arithmetic across grid')
            checked=checked+1
          enddo
          ! CHECK_KINE returns xb by reference; restore the grid value.
          xb=xs(ix)
        enddo
      enddo
    enddo
  enddo
  call require(checked >= 100, 'sufficient valid grid coverage')
  nan=ieee_value(1.0,ieee_quiet_nan)
  call dvmp_t_limits(nan,0.3,lo,hi,valid)
  call require(.not.valid, 'reject nonfinite Q2')
  call dvmp_t_limits(2.0,0.0,lo,hi,valid)
  call require(.not.valid, 'reject zero xB before division')
  call dvmp_t_limits(-2.0,0.3,lo,hi,valid)
  call require(.not.valid, 'reject negative Q2')
  call require(.not.xcheck_kine(nan,0.3,2.0,6.535), 'reject nonfinite t')
  print *, 'PASS boundary grid points:', checked
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ', label
      stop 1
    endif
  end subroutine
end program
