module io_units
! define unit numbers for standard input and output
! Time-stamp: "2001-03-20 15:00:02 cjn"
  implicit none
  private
  public fi, fo
  integer, parameter      :: fi = 5    ! standard input
  integer, parameter      :: fo = 6    ! standard output
end module io_units
    
module precisn
! define working precision kind value
! Time-stamp: "2001-03-20 09:46:41 cjn"
  implicit none
  public
  integer, parameter    :: wp = selected_real_kind(12)
end module precisn
   
module weniger
! Time-stamp: "2001-03-27 14:58:36 cjn"
  use precisn, only: wp
  implicit none
  private
  public levin_weniger

  integer, save               :: minn = 30 ! 60 30
  integer, save               :: maxn = 50 ! 80 50
  !$OMP THREADPRIVATE(minn, maxn)

contains
  subroutine levin_weniger (a, f, n, err)
! wrapper for levin_wen
    real(wp), intent(in)     :: a(0:maxn+1) ! elements of series
    real(wp), intent(out)    :: f        ! final value returned
    integer, intent(out)     :: n        ! # of terms giving 'best' result
    real(wp), intent(out)    :: err      ! estimated uncertainty of result
    integer                  :: l, m
    l = minn
    m = maxn
    call levin_wen (a, l, m, f, n, err)
  end subroutine levin_weniger

  subroutine levin_wen (a, l, m, f, n, err)
! computes the accelerated sum of a series or limit of a sequence.
    integer, intent(in)      :: m        ! most # of terms to be used
    real(wp), intent(in)     :: a(0:m+1) ! elements of series
    integer, intent(in)      :: L        ! least # of terms to be used
    real(wp), intent(out)    :: f        ! final value returned
    integer, intent(out)     :: n        ! # of terms giving 'best' result
    real(wp), intent(out)    :: err      ! estimated uncertainty of result
    logical                  :: improved, start_testing
    logical                  :: previous_iteration
    integer                  :: i
    real(wp)                 :: result(0:m+1), trunc(m+1)
    real(wp)                 :: w1(0:m+1), w2(0:m+1)
    real(wp)                 :: s, test
    real(wp), parameter      :: small = 0.01_wp
    real(wp), parameter      :: huge = 1.0e75_wp

    call wen0 (a, L, w1, w2, result(L), s) ! no test of initial terms
    trunc(L) = 0.0_wp
    improved = .false.
    start_testing = .false.
    previous_iteration = .false.
    test = huge
    n = L
    i_loop: do i = L+1, M
       call wenn (a, i, w1, w2, result(i), s) ! accelerate i-th term
       trunc(i) = ABS(result(i)-result(i-1))  ! truncation error
       improved = (trunc(i) < trunc(i-1)) .or. trunc(i) < small * &
            ABS(result(i))
       start_testing = start_testing .or. (improved .and. previous_iteration)
       previous_iteration = improved
       if (.NOT.start_testing) cycle
       if (trunc(i) < test) then ! find term giving best result
          n = i
          test = trunc(i)
       else
          if (n /= i-1) test = (test + trunc(i)) / 2.0_wp
       end if
    end do i_loop
    f = result(n)
    err = trunc(n)
  end subroutine levin_wen

  subroutine wen0 (a, n, qnum, qden, value, s)
! Weniger algorithm for accelerating a series
    integer, intent(in)      :: n         ! # of elements in a
    real(wp), intent(in)     :: a(0:n+1)  ! terms of series
    real(wp), intent(out)    :: qnum(0:n) ! backward diag of num array
    real(wp), intent(out)    :: qden(0:n) ! backward diag of den array
    real(wp), intent(out)    :: value     ! accelerated value of sum
    real(wp), intent(out)    :: s         ! partial sum of series
    integer                  :: next, j, nmj
    real(wp)                 :: term, termp, c, num

    term = a(0)
    termp = 1.0_wp / a(1)
    s = a(0)
    qden(0) = termp
    qnum(0) = s * termp

    if (n > 0.0_wp) then
       term = a(1)
       termp = 1.0_wp / a(2)
       s = a(1) + s
       qden(1) = termp
       qnum(1) = s * termp
       qden(0) = qden(1) - qden(0)
       qnum(0) = qnum(1) -  qnum(0)
    end if

    do next = 2, n
       term = a(next)
       termp = 1.0_wp / a(next+1)
       s = a(next) + s
       qden(next) = termp
       qnum(next) = s * termp

       num = REAL(next*(next-1),wp)
       j_loop: do j = 1, next
          nmj = next - j
          c = num / REAL((next + j - 1) * (next + j - 2),wp)
          qden(nmj) = qden(nmj+1) - c * qden(nmj)
          qnum(nmj) = qnum(nmj+1) - c *  qnum(nmj)
       end do j_loop
    end do

    value = qnum(0) / qden(0)
  end subroutine wen0

  subroutine wenn (a, n, qnum, qden, value, s)
! Weniger algorithm for accelerating a series
! get n-terms result from result of n-1 terms
    integer, intent(in)      :: n         ! # of elements in a
    real(wp), intent(in)     :: a(0:n+1)  ! elements of series
    real(wp), intent(out)    :: qnum(0:n) ! backward diag of num array
    real(wp), intent(out)    :: qden(0:n) ! backward diag of den array
    real(wp), intent(out)    :: value     ! accelerated value of sum
    real(wp), intent(inout)  :: s         ! partial sum of series
    integer                  :: next, j, nmj
    real(wp)                 :: term, termp, c, num

    next = n
    term = a(next)
    termp = 1.0_wp / a(next+1)
    s = a(next) + s
    qden(next) = termp
    qnum(next) = s * termp

    if (next > 1) then
       num = REAL(next * (next - 1), wp)
       j_loop: do j = 1, next
          nmj = next - j
          c = num / REAL((next + j - 1) * (next + j - 2), wp)
          qden(nmj) = qden(nmj+1) - c * qden(nmj)
          qnum(nmj) = qnum(nmj+1) - c * qnum(nmj)
       end do j_loop
    else if (next == 1) then    ! unlikely case
       qden(0) = qden(1) - qden(0)
       qnum(0) = qnum(1) - qnum(0)
    end if

    value = qnum(0) / qden(0)
  end subroutine wenn
end module weniger

 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
module whittaker_w
! Negative energy decaying Coulomb function
! Whittaker function W(-eta, l+0.5, 2*rho)
! Time-stamp: "2001-05-15 08:30:23 cjn"

! module whittaker_w uses modules:
!    (1) precisn  (definition of working precision kind value)
!    (2) io_units (definition of unit number for error messages)
!    (3) weniger  (sequence transformation)
  use io_units, only: fo
  use precisn, only: wp
  implicit none

  real(wp), save        :: a, c, x ! confluent hypergeometric parameters
  real(wp), save        :: c_true

  integer, save         :: kc                   ! Miller recursion point
  real(wp), save        :: eps = 1.0e-14_wp     ! default error value
  real(wp), parameter   :: small = 1.0e-30_wp   ! scale factor
  integer, parameter    :: mag = 200
  real(wp), parameter   :: vbig = 1.0e+200_wp
  real(wp), parameter   :: vsmall = 1.0e-200_wp

  integer, save         :: method                ! method switch
  logical, save         :: method_reset = .true. ! method override
  logical, save         :: c_recur = .false.     ! c-recurrence switch
  logical, save         :: recur_reset = .true.  ! recurrence override
  
  !$OMP THREADPRIVATE(a, c, x, c_true, kc, eps, method, method_reset, c_recur, recur_reset)

  private
  public   coulomb_whittaker, cw_force_method, cw_force_recur, scale_real

  interface
!  -- LAPACK routine (version 3.0) --
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
     SUBROUTINE DSTEQR (COMPZ, N, D, E, Z, LDZ, WORK, INFO)
       use precisn, only: wp
       CHARACTER(1), intent(in)  :: COMPZ
       INTEGER, intent(in)       :: LDZ, N
       INTEGER, intent(out)      :: INFO
       REAL(wp), intent(inout)   :: D(*), E(*)
       REAL(wp), intent(out)     :: WORK(*)
       REAL(wp), intent(inout)   :: Z(LDZ,*)
     end SUBROUTINE DSTEQR
  end interface

contains

  subroutine cw_force_method (meth)
! override default method selection procedure
    integer, intent(in)     :: meth    ! method to be used
    method_reset = .false.
    method = meth
  end subroutine cw_force_method

  subroutine cw_force_recur (recur)
! override default c-recurrence procedure
    logical, intent(in)     :: recur    ! recurrence flag
    recur_reset = .false.
    c_recur = recur
  end subroutine cw_force_recur

  subroutine coulomb_whittaker (eta, l, rho, w, wd, sf, err)
! Decaying negative energy Coulomb function and derivative wrt rho
    real(wp), intent(in)           :: eta   ! Coulomb parameter
    integer, intent(in)            :: l     ! orbital angular momentum
    real(wp), intent(in)           :: rho   ! radial variable
    real(wp), intent(out)          :: w     ! Whittaker function
    real(wp), intent(out)          :: wd    ! Whittaker derivative
    integer, intent(out)           :: sf    ! scale factor
    real(wp), intent(in), optional :: err   ! error threshold
    real(wp)                       :: kappa, mu, z

    if (PRESENT(err)) eps = err

    kappa = - eta
    mu = REAL(l,wp) + 0.5_wp
    z = 2.0_wp * rho
    call whittaker (kappa, mu, z, w, wd, sf)
    wd = 2.0_wp * wd                   ! return derivative wrt rho
  end subroutine coulomb_whittaker

  subroutine whittaker (kappa, mu, z, w, wd, sf)
! Whittaker W function, W(kappa, mu, z) and derivative
    real(wp), intent(in)            :: kappa
    real(wp), intent(in)            :: mu
    real(wp), intent(in)            :: z
    real(wp), intent(out)           :: w     ! Whittaker function
    real(wp), intent(out), optional :: wd    ! Derivative
    integer, intent(out)            :: sf    ! Scale factor
    real(wp)                        :: y, yn, wt, a_t, c_t, c_r
    integer                         :: n, sfd, nrecur, sft, np, sfr

! set global confluent hypergeometric parameters U(a, c, x):
    a_t = 0.5_wp + mu - kappa
    c_t = 1.0_wp + 2.0_wp * mu
    a = a_t
    c = c_t
    x = z
    
    y = -0.5_wp * z + (0.5_wp + mu) * LOG(z)
    call scale_exponent (y, yn, n)

    if (method_reset) method = select_method()

    c_true = c
    if (recur_reset) then                     ! default choice
       if (c > 12.0_wp .and. method > 2) then ! c-recursion
          nrecur = NINT(c) - 2
          c_r = c - REAL(nrecur,wp)
          c = c_r
       else                                   ! no c-recursion
          nrecur = 0
          c_r = c_t
       end if
    else if (c_recur) then                    ! force recursion
       if (c >= 3.0_wp) then
          nrecur = NINT(c) - 2
          c_r = c - REAL(nrecur,wp)
          c = c_r
       else
          write (fo,'(a,e16.8)') 'whittaker: c ( < 3 ) = ', c
          write (fo,'(a)') 'c-recurrence error'
          stop
       end if
    else                                      ! force no recursion
       nrecur = 0
       c_r = c_t
    end if

    w = psi(sf)
    if (nrecur > 0) then
       a = a_t
       c = c_r + 1.0_wp
       wt = psi(sft)
       if (sft /= sf) wt = wt * 10.0_wp**(sft-sf)
       call c_recursion (w, wt, nrecur, sfr)
       sf = sf + sfr
    end if
    y = yn * w
    call scale_real (y, w, np)
    sf = sf + np

    if (present(wd)) then ! calculate derivative
       a = a_t + 1.0_wp
       if (method_reset) method = select_method()
       if (nrecur > 0) then  ! recur on c
          c = c_r + 1.0_wp
          wd = psi(sfd)
          c = c_r + 2.0_wp
          wt = psi(sft)
          if (sft /= sfd) wt = wt * 10.0_wp**(sft-sfd)
          call c_recursion (wd, wt, nrecur, sfr)
          sfd = sfd + sfr
          wd = - a_t * wd
          y = wd
          call scale_real (y, wd, np)
          np = np + sfd
          if (np /= sf) wd = wd * 10.0_wp**(np-sf)
       else
          c = c_t + 1.0_wp
          wd = - a_t * psi(sfd)
          if (sfd /= sf) wd = wd * 10.0_wp**(sfd-sf)
       end if
       wd = ((mu + 0.5_wp) / z - 0.5_wp) * w  + wd * yn
    end if
    sf = sf + n
  end subroutine whittaker

  function select_method ()
! select method for calculation of function U(a,c,x)
! method determined from globally defined parameters, a, c, x.
! NB: assume x > 0; c >= 2 (integral)
    integer                  :: select_method
    real(wp)                 :: b, d, eta
    integer                  :: meth, na, nd

    na = NINT(a)
    d = c - a - 1.0_wp
    nd = NINT(d)
    eta = a - 0.5_wp * c
    if (ABS(a - REAL(na,wp)) <= small .and. na <= 0) then
       meth = 1      ! Laguerre polynomial case
    else if (ABS(d - REAL(nd,wp)) <= small .and. nd >= 0) then
       meth = 2      ! finite polynomial case

    else if (x <= 0.8_wp) then
       meth = 3      ! accelerated power series expansion

    else if (eta > 0.0_wp) then
       if (x <= 12.0_wp .and. eta > 3.0_wp) then
          meth = 4   ! quadrature
       else
          meth = 5   ! Miller recursion
       end if

    else
       if (x >= 2.0_wp * (0.00250479_wp * eta**2 - &
            0.340344_wp * eta + 8.07783_wp)) then
          meth = 5   ! Miller recursion
       else
          meth = 6   ! recursion from computed psi values
       end if
    end if
    select_method = meth
  end function select_method

  subroutine c_recursion (u1, u2, n, sf)
! upwards recursion on parameter c
! note parameters a, c, x are defined globally
    real(wp), intent(inout)    :: u1
    real(wp), intent(inout)    :: u2
    integer, intent(in)        :: n
    integer, intent(out)       :: sf
    real(wp)                   :: rx, up, u, um, ct
    integer                    :: k

    sf = 0
    u = u2
    um = u1
    rx = 1.0_wp / x
    do k = 1, n-1
       ct = c + REAL(k,wp)
       up = rx * ((2.0_wp + a - ct) * um + (ct + x - 2.0_wp) * u)
       um = u
       u = up
       if (ABS(up) > vbig) then ! scale to avoid overflow/underflow
          sf = sf + mag
          up = vsmall * up
          um = vsmall * um
          u = vsmall * u
       else if (ABS(up) < vsmall) then
          sf = sf - mag
          up = vbig * up
          um = vbig * um
          u = vbig * u
       end if
    end do
    u1= up
    u2 = um
  end subroutine c_recursion

  function psi (sf)
! confluent hypergeometric function psi(a,c;x)
! parameters a, c, x defined by global variables
    real(wp)               :: psi
    integer, intent(out)   :: sf
    real(wp)               :: lgn, ex, ex1, d, psi_t
    integer                :: n1, i, sf1, nd

    select case (method)
    case (1)    ! Special case: Laguerre polynomial
       n1 = - NINT(a)
       psi_t = laguerre_polynomial (n1, c-1.0_wp, x)
       if (MOD(n1,2) == 1) psi_t = - psi_t
       call scale_real (psi_t, psi, sf)
       lgn = 0.0_wp
       do i = 2, n1
          lgn = lgn + LOG(REAL(i,wp))
       end do
       call scale_exponent(lgn, ex1, sf1)
       psi = psi * ex1
       sf = sf + sf1
    case (2)   ! Special case: finite polynomial
       d = c - a - 1.0_wp
       nd = NINT(d)
       psi_t = finite (nd, a, x, sf)
       call scale_real (psi_t, psi, sf1)
       sf = sf + sf1
       ex = d * LOG(x)
       call scale_exponent(ex, ex1, sf1)
       psi = psi * ex1
       sf = sf + sf1
    case (3)   ! sequence transformed power series expansion
       psi = psi_ps (sf)
    case (4)   ! integral representation by numerical quadrature
       psi = psi_qd (sf)
    case (5)   ! Gautschi/Miller recursion
       psi = wimp (sf)
    case (6)   ! downward recursion from positive eta values
       psi = psi_rc (sf)
    case (7)   ! sequence transformed asyptotic series
       psi = psi_as (sf)
    end select
  end function psi

  function psi_rc (sf)
! downward recursion from psi computed at positve a values
    real(wp)                :: psi_rc
    integer, intent(out)    :: sf
    real(wp)                :: wn, wn1, rk, cp, tn, s0, sk
    real(wp)                :: a_true, a_plus_n, eta
    integer                 :: ia_true, na, k, sf1
    integer                 :: sf_recur
    
    a_true = a               ! required a-value
    eta = a - 0.5_wp * c_true
    if (eta >= 0.0_wp) then
       write (fo,'(a,e16.8,a)') 'psi_rc: called with eta = ', eta, &
            ' >= 0'
       stop
    end if
    ia_true = FLOOR(eta)
    na = - ia_true
    a_plus_n = a_true + REAL(na,wp)  

    a = a_plus_n + 1.0_wp
    method = 5               ! force Gautschi method
    wn1 = psi(sf1)
    a = a_plus_n
    method = 5
    wn = psi(sf)

    if (sf1 /= sf) wn1 = wn1 * 10.0_wp**(sf1-sf)
    method = 6               ! restore method switch

    cp = a_true + 1.0_wp - c
    tn = 1.0_wp              ! coefficient calculation
    do k = 1, na+1
       rk = REAL(k,wp) - 1.0_wp
       tn = (a_true + rk) * (cp + rk) * tn / (rk + 1.0_wp)
       
       if (ABS(tn) > vbig) then
          tn = tn * vsmall
          sf1 = sf1 + mag
       else if (ABS(tn) < vsmall .and. tn /= 0.0_wp) then
          tn = tn * vbig
          sf1 = sf1 - mag
       end if
       
    end do
    tn = tn * wn1
    call scale_real(tn, sk, sf1)
    sf = sf + sf1

    rk = a_plus_n * (cp + REAL(na,wp)) * wn1 / (REAL(na+1,wp) * wn)
    a = a_true
    call r_recur (na, rk, sk, s0, sf_recur)
    call scale_real(s0, psi_rc, sf1)
    sf = sf1 + sf + sf_recur
  end function psi_rc

  subroutine r_recur (np, rk, sk, s0, sf)
! downward recursion of r(n)
    integer, intent(in)   :: np ! # iterations
    real(wp), intent(in)  :: rk ! r(n) = w(n+1)/w(n) ratio
    real(wp), intent(in)  :: sk ! w(n+1) value
    real(wp), intent(out) :: s0 ! w(0)
    integer, intent(out)  :: sf
    real(wp)              :: rn, sn, cp, dp, r_k, ak, bk
    integer               :: k

    sf = 0
    rn = rk
    sn = rn / sk
    cp = a - c
    dp = a + cp + x
    recur_down: do k = np, 1, -1
       r_k= REAL(k,wp)
       ak = -(2.0_wp * r_k + dp)
       bk = r_k + 1.0_wp
       rn = - ((r_k + a - 1.0_wp) * (r_k + cp)) / (r_k * &
            (ak + bk * rn))
       sn = rn * sn
       
       if (ABS(sn) > vbig) then
          sn = sn * vsmall
          sf = sf + mag 
       else if (ABS(sn) < vsmall .and. sn /= 0.0_wp) then
          sn = sn * vbig
          sf = sf - mag
       end if
       
    end do recur_down
    s0 = 1.0_wp / sn
    sf = - sf
  end subroutine r_recur

  function wimp (sf) result(psi)
! Gautschi version of Miller algorthm
    real(wp)               :: psi
    integer, intent(out)   :: sf
    real(wp)               :: rk_est, sk_est, s0, psi_t, ex, ex1, chk
    integer                :: sf1, nchk

! if a(n) can vanish, reset kc:
    kc = 1
    chk = 0.5_wp * (c - x) - a - 1.0_wp
    nchk = NINT(chk)
    if (ABS(chk - REAL(nchk,wp)) < small) kc = MAX(kc, nchk+1)

    rk_est = rk() ! evaluate continued fraction, r(kc)
    sk_est = sk() ! evaluate normalization sum, s_n(kc)
    call r_iterate (rk_est, sk_est, s0) ! downward recursion
    psi_t = 1.0_wp / (1.0_wp + s0)
    call scale_real (psi_t, psi, sf)
    ex = - a * LOG(x)
    call scale_exponent(ex, ex1, sf1)
    psi = psi * ex1
    sf = sf + sf1
  end function wimp

  function a_fun(n) RESULT(g)
! psi recursion a(n) coefficient
    real(wp)               :: g
    integer, intent(in)    :: n
    real(wp)               :: rn, rn1
    rn = REAL(n,wp)
    rn1 = REAL(n+1,wp)
    g = - rn1 * (2.0_wp * (rn1 + a) + x - c) / ((rn + a) * &
         (rn1 + a - c))
  end function a_fun

  function b_fun(n) RESULT(h)
! psi recusion b(n) coefficient
    real(wp)               :: h
    integer, intent(in)    :: n
    real(wp)               :: rn, rn1
    rn = REAL(n,wp)
    rn1 = REAL(n+1,wp)
    h = REAL((n+1)*(n+2),wp) / ((rn + a) * (rn1 + a - c))
  end function b_fun

  function rk()
! continued fraction, r(kc) by Lenz-Thompson algorithm
    real(wp)              :: rk
    real(wp)              :: g, h, an, bn, f, df
    integer               :: k
    logical               :: converged
    integer, parameter    :: kmax = 5000

    f = small
    g = f
    h = 0.0_wp
    converged = .false.
    do k = 1, kmax
       if (k == 1) then
          an = - 1.0_wp
       else
          an = - b_fun(kc + k - 2)
       end if
       bn = a_fun(kc + k - 1)

       g = bn + an / g
       if (ABS(g) < small) g = small
       h = bn + an * h
       if (ABS(h) < small) h = small
       h = 1.0_wp / h
       df = h * g
       f  = f * df
       if (ABS(df-1.0_wp) < eps) then
          converged = .true.
          exit
       end if
    end do
    if (.NOT.converged) then
       write (fo,'(a,1p,5e12.3)') 'rk: f, df, eps = ', f, df, eps
    end if
    rk = f
  end function rk

  function sk()
! normalization sum s_N(kc) for N = kc
    real(wp)              :: sk
    real(wp)              :: y_p, y, y_m
    real(wp)              :: rho, rho_m, s_p, s, s_m, af, bf
    integer               :: n
    logical               :: converged, include_y
    integer, parameter    :: kmax = 15000

    s = 0.0_wp
    s_m = 0.0_wp
    y = 1.0_wp
    y_m = 0.0_wp
    rho_m = 0.0_wp
    converged = .false.
    include_y = .true.
    do n = kc+1, kc+kmax
       af = a_fun(n-1)
       bf = b_fun(n-2)
       rho = - 1.0_wp / (af + bf * rho_m)
       if (include_y) then
          y_p = - af * y - bf * y_m
          if (ABS(y_p) <= vsmall) then
             write (fo,'(a,e16.8)') 'sk: y_p too small, y_p = ', y_p
             stop
          end if
          if (ABS(y_p) >= vbig) then
             include_y = .false.
             s_p = - rho * (af * s + rho_m * bf * s_m)
          else
             s_p = - rho * (af * s + rho_m * bf * s_m) + 1.0_wp / y_p
             y_m = y
             y = y_p
          end if
       else
          s_p = - rho * (af * s + rho_m * bf * s_m)
       end if
       if (ABS(s_p - s) <= eps*ABS(s_p)) then
          converged = .true.
          exit
       end if
       s_m = s
       s = s_p
       rho_m = rho
    end do
    if (.NOT.converged) then
       write (fo,'(a,2e24.15)') 'sk: convergence failure, ', s, s_m
    end if
    sk = s_p
  end function sk

  subroutine r_iterate (rk, sk, s0)
! downward recursion of r(n), s(n)
    real(wp), intent(in)  :: rk
    real(wp), intent(in)  :: sk
    real(wp), intent(out) :: s0 ! normalization sum s(0)
    real(wp)              :: rp, rn, sp, sn, cp, dp, r_k, ak, bk
    integer               :: k

    rp = rk
    sp = sk
    cp = a + 1.0_wp - c
    dp = a + cp + 1.0_wp + x   ! dp = 2(a + 1) + x -c
    iteration_down: do k = kc-1, 0, -1
       r_k= REAL(k,wp)
       ak = -(2.0_wp * r_k + dp)
       bk = r_k + 2.0_wp
       rn = - ((r_k + a) * (r_k + cp)) / ((r_k + 1.0_wp) * &
            (ak + bk * rp))
       sn = rn * (1.0_wp + sp)
       rp = rn
       sp = sn
    end do iteration_down
    s0 = sn
  end subroutine r_iterate

  function psi_ps (sf)  result(fx)
! confluent hypergeometric function by seqence transforation of the
! power series about the origin
    use weniger, only: levin_weniger ! sequence transformation
    real(wp)              :: fx
    integer, intent(out)  :: sf
    integer, parameter    :: rmax = 51
    real(wp)              :: u(0:rmax)
    real(wp), parameter   :: psi1 = -0.5772156649015328606512_wp ! psi(1)
    real(wp)              :: f1, f2, sm, rr, aa, bb, cc, s, t, v
    real(wp)              :: est_err, psiq, z
    integer               :: r, nu, m, n

    z = x
    m = NINT(c)     ! assume that c is integral
    n = m - 1
    f1 = -1.0_wp
    sm = 0.0_wp
    do r = 1, n
       rr = REAL(r,wp)
       f1 = - f1 / rr
       sm = sm + 1.0_wp / rr
    end do

    f2 = 1.0_wp / REAL(n,wp)
    do r = 1, n
       rr = REAL(r,wp)
       f2 = f2 * rr / ((a - rr) * z)
    end do

    aa = a - 1.0_wp
    bb = REAL(n,wp)
    cc = a - REAL(n+1,wp)
    t = 1.0_wp
    v = 1.0_wp
    call DIGAMM (a, psiq)
    s = LOG(z) + psiq - 2.0_wp * psi1 - sm
    u(0) = f1 * s + f2
    do r = 1, n-1
       rr = REAL(r,wp)
       t = (aa + rr) * z * t / ((bb + rr) * rr)
       s = s + 1.0_wp / (aa + rr) - 1.0_wp / rr - 1.0_wp / (bb + rr)
       v = (cc + rr) * z * v / ((rr - bb) * rr)
       u(r) = f1 * t * s + f2 * v
    end do
    do r = n, rmax
       rr = REAL(r,wp)
       t = (aa + rr) * z * t / ((bb + rr) * rr)
       s = s + 1.0_wp / (aa + rr) - 1.0_wp / rr - 1.0_wp / (bb + rr)
       u(r) = f1 * t * s
    end do
    call levin_weniger (u, fx, nu, est_err)
    fx = fx / gamma1(a - bb, sf)
    sf = - sf
  end function psi_ps

  function psi_as (sf)
! confluent hypergeometric function by seqence transformation of the
! asymptotic series
    use weniger, only: levin_weniger ! sequence transformation
    real(wp)                 :: psi_as
    integer, intent(out)     :: sf     ! scale factor
    integer                  :: i, nused, sf1
    integer, parameter       :: maxt = 52  ! 52 82
    real(wp)                 :: term(maxt)
    real(wp)                 :: est_err, f2, f3, f4, an, q, fx, fx1
    real(wp)                 :: ex, ex1

    q = 1.0_wp
    term(1) = q
    f2 = a - 1.0_wp
    f3 = a - c
    f4 = - 1.0_wp / x
    do i = 1, maxt-1
       an = REAL(i,wp)
       q = (f2 + an) * (f3 + an) * f4 * q / an
       term(i+1) = q
    end do
    call levin_weniger (term, fx, nused, est_err)
    ex = - a * LOG(x)
    call scale_exponent (ex, ex1, sf1)
    call scale_real (fx, fx1, sf)
    sf = sf + sf1
    psi_as = fx1 * ex1
  end function psi_as

  function laguerre_polynomial (n, k, z)
! laguerroe polynomial for integer k by recursion
    real(wp)                :: laguerre_polynomial
    integer, intent(in)     :: n  ! lower index
    real(wp), intent(in)    :: k  ! upper index
    real(wp), intent(in)    :: z
    real(wp)                :: lnk, lnm1, lnm2
    integer                 :: i

    select case (n)
    case (0)
       lnk = 1.0_wp
    case (1)
       lnk= k + 1.0_wp - z
    case default
       lnm2 = 1.0_wp
       lnm1 = k + 1.0_wp - z
       do i = 1, n-1
          lnk = ((REAL(2*i+1,wp) + k - z) * lnm1 - &
               (k + REAL(i,wp)) * lnm2) / REAL(i+1,wp)
          lnm2 = lnm1
          lnm1 = lnk
       end do
    end select
    laguerre_polynomial = lnk
  end function laguerre_polynomial

  function finite (n, a, z, sf)
! finite series case, c - a -1 = n = 0, 1, 2
    real(wp)              :: finite
    integer, intent(in)   :: n
    real(wp), intent(in)  :: a
    real(wp), intent(in)  :: z
    integer, intent(out)  :: sf
    real(wp)              :: s, t, zi, rk
    integer               :: k

    t = 1.0_wp
    s = t
    zi = 1.0_wp / z
    sf = 0
    do k = 1, n
       rk = REAL(k,wp)
       t = t * REAL(n-k+1,wp) * zi * (a + REAL(k-1,wp)) / REAL(k,wp)
       s = s + t
       if (ABS(t) > vbig .or. ABS(s) > vbig) then ! scale to avoid overflow/underflow
          sf = sf + mag
          t = vsmall * t
          s = vsmall * s
       else if (ABS(t) < vsmall .or. ABS(s) < vsmall) then
          sf = sf - mag
          t = vbig * t
          s = vbig * s
       end if
    end do
    finite = s
  end function finite

  function psi_qd (sf)
! confluent hypergoemetric function psi(a,c,x) by Gaussian quadrature
! confluent hypergeometric parameters a, c, x are global
    real(wp)                    :: psi_qd
    integer, intent(out)        :: sf
    integer                     :: n, status, sf1
    real(wp)                    :: xi, eta_gauss, e1, e2, f
    logical, save               :: gauss = .false.
    integer, parameter          :: n_gauss = 72
    real(wp), save, allocatable :: xg(:), wg(:)

    xi = 1.0_wp / x
    !$OMP CRITICAL(whittaker_gauss_init)
    if (.not.gauss) then ! obtain Gauss quadrature
       n = n_gauss
       allocate (xg(n), wg(n), stat=status)
       if (status /= 0) then
          write (fo,'(a,i6)') 'psi_qd: allocation error = ', status
          stop
       end if
       eta_gauss = 0.0_wp
       call gaulag (n, xg, wg, eta_gauss)
       gauss = .true.
    end if
    !$OMP END CRITICAL(whittaker_gauss_init)

    e1 = a - 1.0_wp
    e2 = c - a - 1.0_wp
    if (a > 100.0_wp) then
       f = SUM(wg(:) * (1.0_wp / (1.0_wp + 1.0_wp / (xi * xg(:))))**e1 &
            * (1.0_wp + xi * xg(:))**(c - 2.0_wp))
    else
       f = SUM(wg(:) * (xi * xg(:))**e1 * (1.0_wp + xi * xg(:))**e2)
    end if
    f = xi * f / gamma1(a, sf1)
    call scale_real (f, psi_qd, sf)
    sf = sf - sf1
  end function psi_qd

  subroutine digamm (p, psi)
! digamm = d(ln gamma(p))/dp; Chebyshev polynomial approximation
! Y.L.Luke, The Special Functions and their Applications Vol.2 (1969)
! p301ff (Academic Press New York and London)
    real(wp), intent(in)       :: p
    real(wp), intent(out)      :: psi
    real(wp), parameter        :: c(18) = (/                      &
         1.09632653123280158001e0_wp, 1.6629180127147181741e-1_wp,     &
         -6.85272202005329673e-3_wp,  3.7339634237860654e-4_wp,        &
         -2.270956946803564e-5_wp,    1.46238871457953e-6_wp,          &
         -9.741775395695e-8_wp,       6.63197083199e-9_wp,             &
         -4.581337814e-10_wp,         3.197121417e-11_wp,              &
         -2.24737976e-12_wp,          1.588114e-13_wp,                 &
         -1.126604e-14_wp,            8.0152e-16_wp,                   &
         -5.715e-17_wp,               4.08e-18_wp,                     &
         -2.9e-19_wp,                 2.0e-20_wp/)
    real(wp)  :: ps, y, t0, t1, t2, co, pi
    integer   :: is, i, np

    ps = p
    if (ps <= 0.0_wp) ps = 1.0_wp - p
    if (ps >= 3.0_wp .and. ps <= 4.0_wp) then
       is = 1
       y = ps - 3.0_wp
    else if (ps >= 4.0_wp) then
       is = 2
       do i = 1, 100000
          y = ps - REAL(3 + i,wp)
          np = i
          if (y <= 1.0_wp) exit
       end do
    else if (ps >= 2.0_wp) then
       is = 3
       y = ps - 2.0_wp
    else if (ps >= 1.0_wp) then
       is = 4
       y = ps - 1.0_wp
    else
       y = ps
       is = 5
    end if

    t0 = 1.0_wp
    t1 = 2.0_wp * y - 1.0_wp
    psi = c(1) * t0 + c(2) * t1
    co = 2.0_wp * t1
    do i = 3, 18
       t2 = co * t1 - t0
       psi = psi + c(i) * t2
       t0 = t1
       t1 = t2
    end do
    select case(is)
    case (2)
       do i = 1, np
          psi = psi + 1.0_wp / (y + REAL(2+i,wp))
       end do
    case (3)
       psi = psi - 1.0_wp / (y + 2.0_wp)
    case (4)
       psi = psi - 1.0_wp / (y + 1.0_wp) - 1.0_wp / (y + 2.0_wp)
    case (5)
       psi = psi - 1.0_wp / y - 1.0_wp / (y + 1.0_wp) - 1.0_wp / &
            (y + 2.0_wp)
    end select
    if (p <= 0.0_wp) then
       pi = 3.141592653589793238463e0_wp
       psi = psi + pi * cos(pi * ps) / sin(pi * ps)
    end if
  end subroutine digamm

  function gamma1 (X, sf)
! GAMMA function for all X except X=0 and X a negative integer
    real(wp)             :: gamma1
    real(wp), intent(in) :: x              ! gamma argument
    integer, intent(out) :: sf             ! scale factor
    real(wp)             :: c, y, b, d, g
    integer              :: cs
    
    
    integer              :: sf_c
    real(wp), parameter  :: vbig_l = 1.0e50_wp
    real(wp), parameter  :: vsmall_l = 1.0e-50_wp
    integer, parameter   :: mag_l = 50

    sf = 0
    c = 1.0_wp
    y = x
    if (y > 10.0_wp) then
       cs = 2
    else if (y >= 2.0_wp) then
       do while (y > 3.0_wp)
          y = y - 1.0_wp
          c = c * y
          if (ABS(c) > vbig_l) then
             c = c * vsmall_l
             sf = sf + mag_l
          end if
       end do
       cs = 1
    else
       sf_c = 0
       do while (y < 2.0_wp)
          c = y * c
          y = y + 1.0_wp
          if (ABS(c) > vbig_l) then
             c = c * vsmall_l
             sf_c = sf_c + mag_l
          end if
       end do
       c = 1.0_wp / c
       cs = 1
       sf = sf - sf_c
    end if
    select case (cs)
    case (1)
       b = y - 2.0_wp
       g = ((((((((((((((-.5113262726698e-6_wp * B +                   &
            .51063592072582e-5_wp) * B -.248410053848712e-4_wp) * B   &
            +.815530498066373e-4_wp) * B -.2064476319159326e-3_wp) *   &
            B + .4677678114964956e-3_wp) * B -                        &
            .9083465574200521e-3_wp) * B + .002099759035077036_wp) *  &
            B -.002851501243034649_wp) * B + .0111538196719067_wp) * B &
            -.2669510287555266e-3_wp) * B + .07424900794340127_wp) *  &
            B + .08157691940138868_wp) * B + .4118403304219815_wp) *  &
            B + .4227843350985181_wp) * B + .9999999999999999_wp
       gamma1 = g * c
    case (2)
       d = 1.0_wp / y
       c = d * d
       g = (((((((((-1.392432216905901_wp * c + .1796443723688306_wp) * &
            c - .02955065359477124_wp) * c + .00641025641025641_wp) *  &
            c - .0019175269175269175_wp) * c + .8417508417508418e-3_wp)&
            * c - .5952380952380952e-3_wp) * c +                        &
            .79365079365079365e-3_wp) * c - .002777777777777778_wp) * c &
            + .08333333333333333_wp) *  d + .9189385332046727_wp +     &
            (y - .5_wp) * LOG(Y) - y
       if (g <= LOG(vbig_l) .and. g >= LOG(vsmall_l)) then
          gamma1 = EXP(g)
       else
          call scale_exponent (g, gamma1, sf)
       end if
    end select
  end function gamma1

  subroutine gaulag (n, x, w, alpha)
! Gauss-Laguerre integration: int(0,inf)  exp(-x) x**alpha f(x) dx
    integer, intent(in)      :: n      ! order of rule
    real(wp), intent(out)    :: x(n)  ! quadrature nodes
    real(wp), intent(out)    :: w(n)  ! quadrature weights
    real(wp), intent(in)     :: alpha  ! quadrature exponent
    integer                  :: i, err, sf
    real(wp)                 :: a(n), b(n), z(n,n)
    real(wp)                 :: work(MAX(1,2*n-2))     ! automatic
    character(1)             :: job

! Laguerre polynomial recurrence coefficients
    do i = 1, n
       a(i) = REAL(2 * i - 1, wp) + alpha
       b(i) = REAL(i * i, wp) + REAL(i, wp) * alpha
    end do

    x = a                                ! diagonal elements
    w(:n-1) = - SQRT(b(1:n-1))           ! sub-diagonal elements

! eigenvalues and eigenvectors of symmetric tridiagonal matrix
    z = 0.0_wp
    do i = 1, n
       z(i,i) = 1.0_wp
    end do
    job = 'i'
    call dsteqr (job, n, x, w(:n-1), z, n, work, err) ! diagonalize
    if (err /= 0) then
       write (fo,'(a,i6)') 'GAULAG: error in dsteqr =', err
       stop
    end if

! calculate weights
    w = gamma1(alpha+1.0_wp, sf) * z(1,:)**2
  end subroutine gaulag

  subroutine scale_exponent (y, x, n)
! return mantissa and exponent of an exponential factor
! EXP(y) = x * 10**n
    real(wp), intent(in)   :: y ! expontial argument
    real(wp), intent(out)  :: x ! mantissa
    integer, intent(out)   :: n ! power of 10
    real(wp)               :: z

    z = y * LOG10(EXP(1.0_wp))
    n = INT(z)
    x = 10.0_wp**(z - REAL(n,wp))
    if (y < 0.0_wp) then
       x = 10.0_wp * x
       n = n - 1
    end if
  end subroutine scale_exponent

  subroutine scale_real (y, x, n)
! return mantissa and exponent of an real number: y = x * 10**n
    real(wp), intent(in)   :: y
    real(wp), intent(out)  :: x
    integer, intent(out)   :: n
    real(wp)               :: z, sgn

    if (y /= 0.0_wp) then
       sgn = SIGN(1.0_wp, y)
       z = LOG10(ABS(y))
       n = INT(z)
       if (ABS(y) > 1.0_wp) n = n + 1
       x = sgn * 10.0_wp**(z - REAL(n,wp))
    else
       x = 0.0_wp
       n = 0
    end if
  end subroutine scale_real

end module whittaker_w
