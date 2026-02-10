!-----------------------------------------------------------------------
!  Copyright (c) 1983, 2007, 2012 A. R Barnett and I. J. Thompson
!  Permission is hereby granted, free of charge,
!  to any person obtaining a copy of this software and associated documentation
!  files (the "Software"), to deal in the Software without restriction, including without
!  limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
!  and/or sell copies of the Software, and to permit persons to whom the Software is
!  furnished to do so, subject to the following conditions:
!
!  The above copyright notice and this permission notice shall be included in all copies or
!  substantial portions of the Software.
!
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
!  INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
!  PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
!  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
!  OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!  OTHER DEALINGS IN THE SOFTWARE.
! ----------------------------------------------------------------------
MODULE KINDS_M
!-----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  INTEGER,PARAMETER :: HPI = SELECTED_INT_KIND(4)        ! 2-BYTE INTEGER
  INTEGER,PARAMETER :: SPI = SELECTED_INT_KIND(9)        ! 4-BYTE INTEGER
  INTEGER,PARAMETER :: DPI = SELECTED_INT_KIND(18)       ! 8-BYTE INTEGER
! ----------------------------------------------------------------------
  INTEGER,PARAMETER :: SPF = SELECTED_REAL_KIND(6,37)    ! 4-BYTE REAL (FLOAT)
  INTEGER,PARAMETER :: DPF = SELECTED_REAL_KIND(15,307)  ! 8-BYTE REAL (FLOAT)
! ----------------------------------------------------------------------
  INTEGER,PARAMETER :: QPI = SELECTED_INT_KIND(9)        !=> 8 BYTE INT
  INTEGER,PARAMETER :: QPF = SELECTED_REAL_KIND(15,307)  !=> 8 BYTE REAL
! ----------------------------------------------------------------------
END MODULE KINDS_M
!-----------------------------------------------------------------------
MODULE STDIO_M
!-----------------------------------------------------------------------
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INPUT_UNIT,OUTPUT_UNIT,ERROR_UNIT
!-----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  INTEGER, PARAMETER :: STDIN =INPUT_UNIT
  INTEGER, PARAMETER :: STDOUT=OUTPUT_UNIT
  INTEGER, PARAMETER :: STDERR=ERROR_UNIT
! ----------------------------------------------------------------------
  CHARACTER(len=72), PARAMETER :: BREAK=REPEAT("*",72)
! ----------------------------------------------------------------------
END MODULE STDIO_M
!-----------------------------------------------------------------------
MODULE LOGAM_M
!-----------------------------------------------------------------------
  USE KINDS_M,ONLY: SPI,DPF
  USE STDIO_M,ONLY: STDOUT
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
  INTEGER(spi),PARAMETER, PRIVATE :: NB=15
  REAL(dpf),DIMENSION(NB),PRIVATE :: B
  REAL(dpf),              PRIVATE :: ACCUR
  REAL(dpf),              PRIVATE :: FPLMIN
  INTEGER(spi),           PRIVATE :: NX0,NT
  !$OMP THREADPRIVATE(B, ACCUR, FPLMIN, NX0, NT)
!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------
    SUBROUTINE LOGAM_INIT(ACC,FPLMIN_IN)
!-----------------------------------------------------------------------
      REAL(dpf),INTENT(IN) :: ACC
      REAL(dpf),INTENT(IN) :: FPLMIN_IN
!-----------------------------------------------------------------------
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: TWO=2._dpf
      INTEGER(spi) :: K
      REAL(dpf) :: X0,F21,ERROR
      LOGICAL :: ACCUR_REACHED
!-----------------------------------------------------------------------
      REAL(dpf),DIMENSION(NB),PARAMETER :: BN= [             +1._dpf &
                                               &,            -1._dpf &
                                               &,            +1._dpf &
                                               &,            -1._dpf &
                                               &,            +5._dpf &
                                               &,          -691._dpf &
                                               &,          +  7._dpf &
                                               &,         -3617._dpf &
                                               &,         43867._dpf &
                                               &,       -174611._dpf &
                                               &,        854513._dpf &
                                               &,    -236364091._dpf &
                                               &,     + 8553103._dpf &
                                               &,  -23749461029._dpf &
                                               &, 8615841276005._dpf]
      REAL(dpf),DIMENSION(NB),PARAMETER :: BD= [              6._dpf &
                                               &,            30._dpf &
                                               &,            42._dpf &
                                               &,            30._dpf &
                                               &,            66._dpf &
                                               &,          2730._dpf &
                                               &,             6._dpf &
                                               &,           510._dpf &
                                               &,           798._dpf &
                                               &,           330._dpf &
                                               &,           138._dpf &
                                               &,          2730._dpf &
                                               &,             6._dpf &
                                               &,           870._dpf &
                                               &,         14322._dpf]
!-----------------------------------------------------------------------
      ACCUR=ACC
      FPLMIN=FPLMIN_IN
!-----------------------------------------------------------------------
      NX0 = 6
      X0  = NX0 + ONE
!-----------------------------------------------------------------------
      ACCUR_REACHED=.FALSE.
!-----------------------------------------------------------------------
      LOOP_120: DO K=1,NB
        F21 = K*2 - ONE
        B(K) = BN(K) / (BD(K) * K*TWO * F21)
        ERROR = ABS(B(K)) * K*TWO / X0**F21
        IF (ERROR.LT.ACCUR) THEN
          ACCUR_REACHED=.TRUE.
          EXIT LOOP_120
        END IF
      END DO LOOP_120
!-----------------------------------------------------------------------
      IF (.NOT.ACCUR_REACHED) THEN
        NX0 = INT((ERROR/ACCUR)**(ONE/F21) * X0)
        K = NB
      END IF
      NT = K
!-----------------------------------------------------------------------
    END SUBROUTINE LOGAM_INIT
!-----------------------------------------------------------------------
!   this routine computes the logarithm of the gamma function gamma(z)
!   for any complex argument 'Z' to any accuracy preset by CALL LOGAM
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CLOGAM(Z)
!-----------------------------------------------------------------------
      COMPLEX(dpf),INTENT(IN) :: Z
!-----------------------------------------------------------------------
      REAL(dpf), PARAMETER :: ZERO=0._dpf
      REAL(dpf), PARAMETER :: HALF=0.5_dpf
      REAL(dpf), PARAMETER :: QUART=0.25_dpf
      REAL(dpf), PARAMETER :: ONE=1._dpf
      REAL(dpf), PARAMETER :: TWO=2._dpf
      REAL(dpf), PARAMETER :: FOUR=4._dpf
      REAL(dpf), PARAMETER :: PI=FOUR*ATAN(ONE)
      REAL(dpf), PARAMETER :: ALPI = LOG(PI)
      REAL(dpf), PARAMETER :: HL2P = LOG(TWO*PI)*HALF
      COMPLEX(dpf) :: V,H,R,SER
      INTEGER(spi) :: N,J,K,MX,I
      REAL(dpf) :: X,A,C,D,E,F,T
!-----------------------------------------------------------------------
 1000 FORMAT(1X,A6,' ... ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)
!-----------------------------------------------------------------------
      X=Z%RE
      T=Z%IM
      MX=INT(REAL(ACCUR*100._dpf-X,KIND=dpf),KIND=spi)
      IF (ABS(ABS(X)-MX)+ABS(T).LT.ACCUR*50._dpf) THEN
        WRITE(STDOUT,1000) 'CLOGAM',X
        CLOGAM = ZERO
      ELSE
        F=ABS(T)
        V=CMPLX(X,F,KIND=dpf)
        IF (X.LT.ZERO) V=ONE-V
        H=ZERO
        C=V%RE
        N=NX0-INT(C,KIND=spi)
        IF (N.GE.0_spi) THEN
          H=V
          D=V%IM
          A=ATAN2(D,C)
          IF (N.NE.0) THEN
            DO I = 1,N
              C=C+ONE
              V=CMPLX(C,D,KIND=dpf)
              H=H*V
              A=A+ATAN2(D,C)
            END DO
          END IF
          H=CMPLX(HALF*LOG(H%RE**2+H%IM**2),A,KIND=dpf)
          V=V+ONE
        END IF
        R=ONE/V**2
        SER = B(NT)
        DO J=2,NT
          K = NT+1 - J
          SER = B(K) + R*SER
        END DO
!-----------------------------------------------------------------------
        CLOGAM = HL2P+(V-HALF)*LOG(V)-V + SER/V - H
!-----------------------------------------------------------------------
        IF (X.LT.ZERO) THEN
!-----------------------------------------------------------------------
          A=INT(X,KIND=spi)-ONE
          C=PI*(X-A)
          D=PI*F
          E=ZERO
          F=-TWO*D
          IF (F.GT.FPLMIN) E = EXP(F)
          F=SIN(C)
          E= D + HALF*LOG(E*F**2+QUART*(ONE-E)**2)
          F=ATAN2(COS(C)*TANH(D),F)-A*PI
          CLOGAM=ALPI-CMPLX(E,F,KIND=dpf)-CLOGAM
        END IF
!-----------------------------------------------------------------------
        IF (SIGN(ONE,T).LT.-HALF) CLOGAM=CONJG(CLOGAM)
      END IF
!-----------------------------------------------------------------------
    END FUNCTION CLOGAM
!-----------------------------------------------------------------------
!   this routine computes the logarithmic derivative of the gamma
!   function  psi(Z) = digamma(Z) = d (ln gamma(Z))/dZ  for any
!   complex argument Z, to any accuracy preset by CALL LOGAM(ACC)
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CDIGAM(Z)
!-----------------------------------------------------------------------
      COMPLEX(dpf),INTENT(IN) :: Z
!-----------------------------------------------------------------------
      REAL(dpf), PARAMETER :: ZERO=0._dpf
      REAL(dpf), PARAMETER :: HALF=0.5_dpf
      REAL(dpf), PARAMETER :: ONE=1._dpf
      REAL(dpf), PARAMETER :: FOUR=4._dpf
      REAL(dpf), PARAMETER :: PI=FOUR*ATAN(ONE)
      COMPLEX(dpf) :: U,V,H,R,SER
      INTEGER(spi) :: N,J,K,I
      REAL(dpf) :: X,A
!-----------------------------------------------------------------------
 1000 FORMAT(1X,A6,' ... ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)
!-----------------------------------------------------------------------
      U=Z
      X=U%RE
      A=ABS(X)
      IF (ABS(U%IM) + ABS(A+INT(X,KIND=spi)).LT.ACCUR) THEN
!-----------------------------------------------------------------------
        WRITE(STDOUT,1000) 'CDIGAM',X
        CDIGAM=ZERO
!-----------------------------------------------------------------------
      ELSE
        IF (X.LT.ZERO) U=-U
        V=U
        H=ZERO
        N=NX0-INT(A,KIND=spi)
        IF (N.GE.0_spi) THEN
          H=ONE/V
          IF (N.NE.0) THEN
            DO I = 1,N
              V=V+ONE
              H=H+ONE/V
            END DO
          END IF
          V=V+ONE
        END IF
        R=ONE/V**2
        SER = B(NT) * (2*NT-1)
        DO J=2,NT
          K = NT+1 - J
          SER = B(K)*(2*K-1) + R*SER
        END DO
!-----------------------------------------------------------------------
        CDIGAM = LOG(V) - HALF/V - R*SER - H
!-----------------------------------------------------------------------
        IF (X.LT.ZERO) THEN
          H=PI*U
          CDIGAM = CDIGAM + ONE/U + PI*COS(H)/SIN(H)
        END IF
!-----------------------------------------------------------------------
      END IF
!-----------------------------------------------------------------------
    END FUNCTION CDIGAM
!-----------------------------------------------------------------------
END MODULE LOGAM_M
!-----------------------------------------------------------------------
MODULE RCF_M
!-----------------------------------------------------------------------
  USE KINDS_M,ONLY: SPI,DPF
  USE STDIO_M,ONLY: STDOUT
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
  LOGICAL     , PRIVATE :: EVEN
  COMPLEX(dpf), PRIVATE :: X1
  INTEGER(spi), PRIVATE :: M2M1,MP12,M
  !$OMP THREADPRIVATE(EVEN, X1, M2M1, MP12, M)
!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------
!  RCF converts polynomial A to the corresponding continued
!         fraction, in 'normal'  form with coefficients B
!         by the 'P algorithmn' of Patry & Gupta
!-----------------------------------------------------------------------
!   A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)
!-----------------------------------------------------------------------
!   B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)
!-----------------------------------------------------------------------
!  data:
!   A     vector A(k), k=1,INUM         input
!   B     vector B(k), k=IBEG,INUM      output
!   IBEG  order of first coef. calc.    input
!   INUM  order of A, even or odd       input
!   XX    auxiliary vector of length .ge. length of vector B
!         caller provides space for A,B,XX
!     Note that neither of the first two terms A(1) A(2) should be zero
!             & the user can start the calculation with any value of
!                IBEG provided the c.f. coefs have been already
!                calculated up to INUM = IBEG-1
!             & the method breaks down as soon as the absolute value
!                of a c.f. coef. is less than EPS.    At the time of the
!                break up XX(1) has been replaced by 1E-50, and INUM has
!                been replaced by minus times the number of this coef.
!   algorithm: J.Patry & S.Gupta,
!              EIR-bericht nr. 247,
!              Eidg. Institut fur Reaktorforschung Wuerenlingen
!              Wueringlingen, Schweiz.
!              November 1973
!   see also:  Haenggi,Roesel & Trautmann,
!              Jnl. Computational Physics, vol 137, pp242-258 (1980)
!   note:      restart procedure modified by I.J.Thompson
!   note:      ierror added by A.R. FLORES
!-----------------------------------------------------------------------
    SUBROUTINE RCF(A,B,IBEG,INUM,XX,EPS,IERROR)
!-----------------------------------------------------------------------
      INTEGER(spi),                 INTENT(IN)    :: IBEG
      INTEGER(spi),                 INTENT(IN)    :: INUM
      COMPLEX(dpf),DIMENSION(100),  INTENT(IN)    :: A
      COMPLEX(dpf),DIMENSION(100),  INTENT(INOUT) :: B
      COMPLEX(dpf),DIMENSION(2,100),INTENT(INOUT) :: XX
      REAL(dpf),                    INTENT(IN)    :: EPS
      INTEGER(spi),                 INTENT(OUT)   :: IERROR
!-----------------------------------------------------------------------
      INTEGER(spi) :: K,IBN,II
      COMPLEX(dpf) :: X0
!-----------------------------------------------------------------------
 1000 FORMAT('0RCF: LAST CALL SET M =',I4,', BUT RESTART REQUIRES',I4)
!-----------------------------------------------------------------------
      IERROR=0_spi
!-----------------------------------------------------------------------
      X0 = CMPLX(0._dpf,0._dpf,KIND=dpf)
!-----------------------------------------------------------------------
      IF (IBEG.GT.4 .AND. M .NE. IBEG-1) THEN
        WRITE(STDOUT,1000) M,IBEG-1
        STOP ("RCF HAS FAILED")
      END IF
!-----------------------------------------------------------------------
!     B(IBN) is last value set on this call
!     B(M) is last value set in previous call
!-----------------------------------------------------------------------
      IBN = INUM
!-----------------------------------------------------------------------
      IF (IBEG.GT.4) THEN
!-----------------------------------------------------------------------
        DO K=M2M1,2,-1
          XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)
        END DO
        XX(2,1) = XX(1,1) + B(M)
        DO K=1,M2M1
          X0 = XX(2,K)
          XX(2,K) = XX(1,K)
          XX(1,K) = X0
        END DO
        X0 = X1
        XX(1,M2M1+1) = 0.
        M = M+1
        EVEN = .NOT.EVEN
!-----------------------------------------------------------------------
      ELSE !IBEG.LE.4
!-----------------------------------------------------------------------
        IF (IBEG.LT.4) THEN
          B(1) = A(1)
          IF (IBN.GE.2) B(2) = - A(2)/A(1)
          IF (IBN.LT.3) RETURN
          X0 = A(3) / A(2)
          XX(2,1) = B(2)
          XX(1,1) = - X0
          XX(1,2) = 0.
          B(3) = -X0 - B(2)
          X0 = -B(3) * A(2)
          M = 3
          MP12 = 2
          EVEN = .TRUE.
          IF (IBN.LE.3) RETURN
        END IF
!-----------------------------------------------------------------------
        IF (ABS(B(3)) .LT. EPS*ABS(X0)) THEN
          IERROR = -M
          RETURN
        END IF
!-----------------------------------------------------------------------
        M = 4
!-----------------------------------------------------------------------
      END IF
!-----------------------------------------------------------------------
      LOOP_RCF: DO II=1,IBN-M+1
        X1 = A(M)
        M2M1 = MP12
        MP12 = M2M1 + 1
        IF (EVEN) MP12 = M2M1
        DO K=2,MP12
          X1 = X1 + A(M-K+1) * XX(1,K-1)
        END DO
        B(M) = - X1/X0
!-----------------------------------------------------------------------
        IF (M.GE.IBN) EXIT LOOP_RCF
!-----------------------------------------------------------------------
        IF (ABS(B(M)).LT.EPS*ABS(X0)) THEN
          IERROR = -M
          EXIT LOOP_RCF
        END IF
!-----------------------------------------------------------------------
        DO K=M2M1,2,-1
          XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)
        END DO
        XX(2,1) = XX(1,1) + B(M)
        DO K=1,M2M1
          X0 = XX(2,K)
          XX(2,K) = XX(1,K)
          XX(1,K) = X0
        END DO
        X0 = X1
        XX(1,M2M1+1) = 0.
        M = M+1
        EVEN = .NOT.EVEN
!-----------------------------------------------------------------------
      END DO LOOP_RCF
!-----------------------------------------------------------------------
    END SUBROUTINE RCF
!-----------------------------------------------------------------------
END MODULE RCF_M
!-----------------------------------------------------------------------
MODULE CSTEED_M
!-----------------------------------------------------------------------
  USE KINDS_M,ONLY: SPI,DPF
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
!***  This module is for information & storage only.
!     (it is not essential to working of the code)
!-----------------------------------------------------------------------
    REAL(dpf)    :: RERR
    INTEGER(spi) :: NFP,N11,N20
    INTEGER(spi),DIMENSION(2) :: NPQ,KAS
    !$OMP THREADPRIVATE(RERR, NFP, N11, N20, NPQ, KAS)
END MODULE CSTEED_M
!-----------------------------------------------------------------------
MODULE COULCC_M
!-----------------------------------------------------------------------
  USE KINDS_M,ONLY: SPI,DPF,QPF
  USE STDIO_M,ONLY: STDOUT
!-----------------------------------------------------------------------
  USE LOGAM_M
  USE RCF_M
  USE CSTEED_M
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------
!  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD
!-----------------------------------------------------------------------
!  A CONTINUED-FRACTION ALGORITHM FOR COULOMB
!  FUNCTIONS OF COMPLEX ORDER WITH COMPLEX ARGUMENTS
!-----------------------------------------------------------------------
!  REF. IN COMP. PHYS. COMMUN. 36 (1985) 363
!-----------------------------------------------------------------------
!  A. R. Barnett           Manchester  March   1981
!  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions
!  modified A.R. FLORES    WashU, April 2025 for Modern Compilers
!-----------------------------------------------------------------------
!  original program  RCWFN       in    CPC  8 (1974) 377-395
!                 +  RCWFF       in    CPC 11 (1976) 141-142
!                 +  COULFG      in    CPC 27 (1982) 147-166
!  description of real algorithm in    CPC 21 (1981) 297-314
!  description of complex algorithm    JCP XX (1985) YYY-ZZZ
!  this version written up       in    CPC XX (1985) YYY-ZZZ
!-----------------------------------------------------------------------
!  COULCC returns F,G,G',G',SIG for complex XX, ETA1, and ZLMIN,
!   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL-1 inclusive,
!   thus giving  complex-energy solutions to the Coulomb Schrodinger
!   equation,to the Klein-Gordon equation and to suitable forms of
!   the Dirac equation ,also spherical & cylindrical Bessel equations
!
!  if /MODE1/= 1  get F,G,F',G'   for integer-spaced lambda values
!            = 2      F,G      unused arrays must be dimensioned in
!            = 3      F,  F'          call to at least length (1)
!            = 4      F
!            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )
!            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in
!            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC
!            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )
!
!     if MODE1<0 then the values returned are scaled by an exponential
!                factor (dependent only on XX) to bring nearer unity
!                the functions for large /XX/, small ETA & /ZL/ < /XX/
!        Define SCALE = (  0        if MODE1 > 0
!                       (  IMAG(XX) if MODE1 < 0  &  KFN < 3
!                       (  REAL(XX) if MODE1 < 0  &  KFN = 3
!        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)
!         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )
!               or EXP(SCALE)       * ( H+, H(1), or K)
!               or EXP(-SCALE)      * ( H- or H(2) )
!
!  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G
!           =  1   spherical Bessel      "      "     "       j & y
!           =  2 cylindrical Bessel      "      "     "       J & Y
!           =  3 modified cyl. Bessel    "      "     "       I & K
!
!          and where Coulomb phase shifts put in SIG if KFN=0 (not -1)
!
!  The use of MODE and KFN is independent
!    (except that for KFN=3,  H(1) & H(2) are not given)
!
!  With negative orders lambda, COULCC can still be used but with
!    reduced accuracy as CF1 becomes unstable. The user is thus
!    strongly advised to use reflection formulae based on
!    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi)
!
!  Precision:  results to within 2-3 decimals of 'machine accuracy',
!               but if CF1A fails because X too small or ETA too large
!               the F solution  is less accurate if it decreases with
!               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)
!              RERR in CSTEED module traces the main roundoff errors.
!
!   COULCC is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14
!          with a section of doubled REAL*16 for less roundoff errors.
!          (If no doubled precision available, increase JMAX to eg 100)
!   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32
!   For single precision CDC (48 bits) reassign REAL*8=REAL etc.
!
!   IFAIL  on input   = 0 : no printing of error messages
!                    ne 0 : print error messages on file 6
!   IFAIL  in output = -2 : argument out of range
!                    = -1 : one of the continued fractions failed,
!                           or arithmetic check before final recursion
!                    =  0 : All Calculations satisfactory
!                    ge 0 : results available for orders up to & at
!                             position NL-IFAIL in the output arrays.
!                    = -3 : values at ZLMIN not found as over/underflow
!                    = -4 : roundoff errors make results meaningless
!-----------------------------------------------------------------------
!     ROUTINES CALLED :       LOGAM/CLOGAM/CDIGAM,
!                             F20, CF1A, RCF, CF1C, CF2, F11, CF1R
!     Intrinsic functions :   MIN, MAX, SQRT, REAL, IMAG, ABS, LOG, EXP
!      (Generic names)        NINT, MOD, ATAN, ATAN2, COS, SIN, CMPLX,
!                             SIGN, CONJG, INT, TANH
!-----------------------------------------------------------------------
!     Machine dependent Parameters
!-----------------------------------------------------------------------
!     ACCUR  ::  target bound on relative error (except near 0 crossings)
!                (ACCUR should be at least 100 * ACC8)
!     ACC8   ::  smallest number with 1+ACC8 .ne.1 in REAL*8  arithmetic
!     ACC16  ::  smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic
!     FPMAX  ::  magnitude of largest floating point number * ACC8
!     FPMIN  ::  magnitude of smallest floating point number / ACC8
!     FPLMAX ::  LOG(FPMAX)
!     FPLMIN ::  LOG(FPMIN)
!-----------------------------------------------------------------------
!     Parameters determining region of calculations
!-----------------------------------------------------------------------
!     R20   ::   estimate of (2F0 iterations)/(CF2 iterations)
!     ASYM  ::   minimum X/(ETA**2+L) for CF1A to converge easily
!     XNEAR ::   minimum ABS(X) for CF2 to converge accurately
!     LIMIT ::   maximum no. iterations for CF1, CF2, and 1F1 series
!     JMAX  ::   size of work arrays for Pade accelerations
!     NDROP ::   number of successive decrements to define instability
!-----------------------------------------------------------------------
    SUBROUTINE COULCC(XX,ETA1,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE1,KFN,IFAIL)
!-----------------------------------------------------------------------
      COMPLEX(dpf),              INTENT(IN)    :: XX,ETA1,ZLMIN
      INTEGER(spi),              INTENT(IN)    :: NL,MODE1,KFN
      COMPLEX(dpf),DIMENSION(NL),INTENT(OUT)   :: FC,GC,FCP,GCP,SIG
      INTEGER(spi),              INTENT(INOUT) :: IFAIL
!-----------------------------------------------------------------------
      INTEGER(spi),PARAMETER :: JMAX=50_spi
      INTEGER(spi),PARAMETER :: LIMIT=20000_spi
      REAL(dpf),   PARAMETER :: R20=3._dpf
      REAL(dpf),   PARAMETER :: ASYM=3._dpf
      REAL(dpf),   PARAMETER :: XNEAR=0.5_dpf
      INTEGER(spi),PARAMETER :: NDROP=5_spi
!-----------------------------------------------------------------------
      REAL(dpf),   PARAMETER :: ZERO=0._dpf
      REAL(dpf),   PARAMETER :: HALF=0.5_dpf
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: TWO=2._dpf
      COMPLEX(dpf),PARAMETER :: CI=CMPLX(0._dpf,1._dpf,KIND=dpf)
      COMPLEX(dpf),PARAMETER :: CZERO=CMPLX(ZERO,ZERO,KIND=dpf)
      REAL(dpf),   PARAMETER :: HPI=TWO*ATAN(ONE)
      REAL(dpf),   PARAMETER :: TLOG=LOG(TWO)
!-----------------------------------------------------------------------
      REAL(dpf),   PARAMETER :: ACC8  =EPSILON(1._dpf)
      REAL(dpf),   PARAMETER :: ACC16 =EPSILON(1._qpf)
      REAL(dpf),   PARAMETER :: FPMAX =HUGE(1._dpf) * ACC8
      REAL(dpf),   PARAMETER :: FPMIN =TINY(1._dpf) / ACC8
      REAL(dpf),   PARAMETER :: FPLMAX=LOG(FPMAX)
      REAL(dpf),   PARAMETER :: FPLMIN=LOG(FPMIN)
!-----------------------------------------------------------------------
      COMPLEX(dpf),DIMENSION(JMAX,4) :: XRCF
!-----------------------------------------------------------------------
      LOGICAL      :: PR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2
      REAL(dpf)    :: ACCUR,ERROR,ACCT,ACCH,ACCB,PACCQ,EPS,OFF,SCALE,SF
      REAL(dpf)    :: SFSH,TA,RK,OMEGA,ABSX
      INTEGER(spi) :: ID,I,IH,KASE,L,L1,LAST,LF,LH,M1,MODE,MONO,N
      COMPLEX(dpf) :: X,ETA,ETAI,DELL,ZM1,AA,AB,ALPHA,BB,BETA,CHI,CIK
      COMPLEX(dpf) :: CLGAA,CLGAB,CLGBB,CLL,DF,DSIG,EK,ETAP,F,F11V,F20V
      COMPLEX(dpf) :: FCL,FCL1,FCM,FESL,FEST,FIRST,FPL,GAM,HCL,HCL1,P,P11
      COMPLEX(dpf) :: PK,PL,PM,PQ1,PQ2,Q,RL,SIGMA,SL,THETA,THETAM,TPK1,W
      COMPLEX(dpf) :: XI,XLOG,Z11,ZID,ZL,ZLL,ZLM,ZLOG,HPL,UC
!-----------------------------------------------------------------------
      NOCF2=.TRUE.
      KASE = 0_spi
      PACCQ = ZERO
      SFSH = ZERO
      OMEGA = ZERO
      F20V = CZERO
      CLL = CZERO
      PQ1 = CZERO
      HPL = CZERO
      THETAM = CZERO
      Q = CZERO
!-----------------------------------------------------------------------
      MODE = MOD(ABS(MODE1),10)
      IFCP = MOD(MODE,2).EQ.1
      PR = IFAIL.NE.0 !PRINT REPORT
      IFAIL = -2
      N11   = 0
      NFP   = 0
      KAS(:)= 0
      NPQ(:)= 0
      N20 = 0
      ACCUR = MAX((1._dpf)*(10._dpf)**(-14_spi), 50*ACC8)
      ACCT = ACCUR * .5
!-----------------------------------------------------------------------
      CALL LOGAM_INIT(ACC8,FPLMIN)
!-----------------------------------------------------------------------
      ACCH = SQRT(ACCUR)
      ACCB = SQRT(ACCH)
      RERR = ACCT
!-----------------------------------------------------------------------
      CIK = ONE
      IF (KFN.GE.3) CIK = CI * SIGN(ONE,ACC8-XX%IM)
      X     = XX * CIK
      ETA   = ETA1
      IF (KFN .GT. 0) ETA = ZERO
      ETANE0  = ABSC(ETA).GT.ACC8
      ETAI = ETA*CI
      DELL  = ZERO
      IF (KFN .GE. 2)  DELL = HALF
      ZM1   = ZLMIN - DELL
      SCALE = ZERO
      IF (MODE1.LT.0) SCALE = X%IM
!-----------------------------------------------------------------------
      M1 = 1
      L1  = M1 + NL - 1
      RLEL = ABS(ETA%IM) + ABS(ZM1%IM) .LT. ACC8
      ABSX = ABS(X)
      AXIAL = RLEL .AND. ABS(X%IM) .LT. ACC8 * ABSX
      IF (MODE.LE.2 .AND. ABSX.LT.FPMIN) GOTO 310
      XI  = ONE/X
      XLOG = LOG(X)
!-----------------------------------------------------------------------
! log with cut along the negative real axis see also OMEGA
!-----------------------------------------------------------------------
      ID = 1
      DONEM = .FALSE.
      UNSTAB = .FALSE.
      LF = M1
      IFAIL = -1
  10  ZLM = ZM1 + LF - M1
      ZLL = ZM1 + L1 - M1
!-----------------------------------------------------------------------
! *** ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels
!-----------------------------------------------------------------------
      Z11 = ZLL
      IF (ID.LT.0) Z11 = ZLM
      P11 = CI*SIGN(ONE,ACC8-ETA%IM)
      LAST = L1
!-----------------------------------------------------------------------
! *** Find phase shifts and Gamow factor at lambda = ZLL
!-----------------------------------------------------------------------
      PK = ZLL + ONE
      AA = PK - ETAI
      AB = PK + ETAI
      BB = TWO*PK
      ZLNEG = NPINT(BB,ACCB)
      CLGAA = CLOGAM(AA)
      CLGAB = CLGAA
      IF (ETANE0.AND..NOT.RLEL)  CLGAB = CLOGAM(AB)
      IF (ETANE0.AND.     RLEL)  CLGAB = CONJG(CLGAA)
      SIGMA = (CLGAA - CLGAB) * CI*HALF
      IF (KFN.EQ.0) SIG(L1) = SIGMA
      IF(.NOT.ZLNEG) CLL = ZLL*TLOG- HPI*ETA - CLOGAM(BB) + (CLGAA+CLGAB)*HALF
      THETA  = X - ETA*(XLOG+TLOG) - ZLL*HPI + SIGMA
!-----------------------------------------------------------------------
      TA = (AA%IM**2+AB%IM**2+ABS(AA%RE)+ABS(AB%RE))*HALF
      IF (ID.GT.0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GOTO 20
!-----------------------------------------------------------------------
! ***  use CF1 instead of CF1A, if predicted to converge faster,
!          (otherwise using CF1A as it treats negative lambda &
!           recurrence-unstable cases properly)
!-----------------------------------------------------------------------
      RK = SIGN(ONE, X%RE + ACC8)
      P =  THETA
      IF (RK.LT.0) P = -X + ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA
!-----------------------------------------------------------------------
! FLORES -- I.T. HAS jmax/2 --  seems jmax/2 fixes axial problem for second test case
      !F = RK * CF1A(X*RK,ETA*RK,ZLL,P,ACCT,JMAX,NFP,FEST,ERROR,FPMAX,XRCF,XRCF(1,3), XRCF(1,4))
!-----------------------------------------------------------------------
      F = RK * CF1A(X*RK,ETA*RK,ZLL,P,ACCT,JMAX/2,NFP,FEST,ERROR,FPMAX,XRCF,XRCF(1,3), XRCF(1,4))
      FESL = LOG(FEST) + ABS(X%IM)
      NFP = - NFP
      IF (NFP.LT.0   .OR.(UNSTAB.AND.ERROR.LT.ACCB)) GOTO 40
      IF(.NOT.ZLNEG .OR. UNSTAB.AND.ERROR.GT.ACCB)  GOTO 20
      IF (PR) WRITE(STDOUT,1060) '-L',ERROR
      IF (ERROR.GT.ACCB) GOTO 280
      GOTO 40
!-----------------------------------------------------------------------
! ***  evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)
!-----------------------------------------------------------------------
  20  IF (AXIAL) THEN
!-----------------------------------------------------------------------
!  REAL VERSION
!-----------------------------------------------------------------------
        F = CF1R(X%RE,ETA%RE,ZLL%RE,ACC8,SF,RK,ETANE0,LIMIT,ERROR,NFP,ACCH,FPMIN,FPMAX,PR,'COULCC')
        FCL = SF
        TPK1= RK
      ELSE
!-----------------------------------------------------------------------
!  COMPLEX VERSION
!-----------------------------------------------------------------------
        F = CF1C(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERROR,NFP,ACCH,FPMIN,FPMAX,PR,'COULCC')
      END IF
!-----------------------------------------------------------------------
      IF (ERROR.GT.ONE) GOTO 390
!-----------------------------------------------------------------------
! ***  Make a simple check for CF1 being badly unstable:
!-----------------------------------------------------------------------
      IF (ID.LT.0) GOTO 30
      UC = (CMPLX(ONE,0._dpf,KIND=dpf)-ETA*XI)*CI*THETA%IM/F
      UNSTAB = UC%RE.GT.ZERO     &
        & .AND. .NOT.AXIAL &
        & .AND. ABS(THETA%IM).GT.-LOG(ACC8)*HALF &
        & .AND. ABSC(ETA)+ABSC(ZLL).LT.ABSC(X)
      IF (UNSTAB) GOTO 60
!-----------------------------------------------------------------------
! *** compare accumulated phase FCL with asymptotic phase for G(k+1) :
!     to determine estimate of F(ZLL) (with correct sign) to start recur
!-----------------------------------------------------------------------
! FLORES -- I.T. HAS -- 
! this flips the sign of test: 0X=20.,-0.0010, ETA=0.,0., ZLMIN=0.,0. NL=11 MODE=2 KFN=1 LENTZ: J,Y SPH.BESSE
!   30 W   = LOG(ONE+TWO*ETA*X/(TPK1-ONE)+X*X/TPK1*(TWO*ETA*ETA/(TPK1-ONE)-HALF))
!      FESL= (ZLL+ONE) * XLOG + CLL + W - LOG(FCL)
!-----------------------------------------------------------------------
  30  W    =  X*X  *(HALF/TPK1 + ONE/TPK1**2) + ETA*(ETA-TWO*X)/TPK1
      FESL = (ZLL+ONE) * XLOG + CLL - W - LOG(FCL)
!-----------------------------------------------------------------------
  40  FESL = FESL - ABS(SCALE)
      RK   = MAX(FESL%RE,FPLMIN*HALF)
      FESL = CMPLX(MIN(RK,FPLMAX*HALF),FESL%IM,KIND=dpf)
      FEST = EXP(FESL)
!-----------------------------------------------------------------------
      RERR = MAX(RERR,ERROR,ACC8*ABS(THETA%RE))
!-----------------------------------------------------------------------
      FCL = FEST
      FPL = FCL*F
      IF (IFCP) FCP(L1) = FPL
      FC (L1) = FCL
!-----------------------------------------------------------------------
! *** downward recurrence to lambda = ZLM. array GC,if present,stores RL
!-----------------------------------------------------------------------
      I  = MAX(-ID, 0)
      ZL  = ZLL + I
      MONO = 0
      OFF = ABS(FCL)
      TA = ABSC(SIGMA)
      LOOP_70: DO L = L1-ID,LF,-ID
        IF (ETANE0) THEN
          IF (RLEL) THEN
            DSIG = ATAN2(ETA%RE,ZL%RE)
            RL = SQRT(ZL%RE**2 + ETA%RE**2)
          ELSE
            AA = ZL - ETAI
            BB = ZL + ETAI
            IF (ABSC(AA).LT.ACCH.OR.ABSC(BB).LT.ACCH) GOTO 50
            DSIG = (LOG(AA) - LOG(BB)) * CI*HALF
            RL = AA * EXP(CI*DSIG)
          END IF
          IF (ABSC(SIGMA).LT.TA*HALF) THEN
!-----------------------------------------------------------------------
! re-calculate SIGMA because of accumulating roundoffs:
!-----------------------------------------------------------------------
            SL =(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CI*HALF
            RL = (ZL - ETAI) * EXP(CI*ID*(SIGMA - SL))
            SIGMA = SL
            TA = ZERO
          ELSE
            SIGMA = SIGMA - DSIG * ID
          END IF
          TA = MAX(TA, ABSC(SIGMA))
          SL =  ETA  + ZL*ZL*XI
          PL = ZERO
          IF (ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZL
          FCL1  = (FCL *SL + ID*ZL*FPL)/RL
          SF = ABS(FCL1)
          IF (SF.GT.FPMAX) GOTO 350
          FPL   = (FPL *SL + ID*PL*FCL)/RL
          IF (MODE .LE. 1) GCP(L+ID)= PL * ID
        ELSE
!-----------------------------------------------------------------------
!  ETA = 0, including Bessels.  NB RL==SL
!-----------------------------------------------------------------------
          RL = ZL* XI
          FCL1 = FCL * RL + FPL*ID
          SF = ABS(FCL1)
          IF (SF.GT.FPMAX) GOTO 350
          FPL  =(FCL1* RL - FCL) * ID
        END IF
        IF (SF.LT.OFF) THEN
          MONO = MONO + 1
        ELSE
          MONO = 0
        END IF
        FCL =  FCL1
        OFF = SF
        FC(L) =  FCL
        IF (IFCP) FCP(L)  = FPL
        IF (KFN.EQ.0) SIG(L) = SIGMA
        IF (MODE .LE. 2) GC(L+ID) = RL
        ZL = ZL - ID
        IF (MONO.LT.NDROP) CYCLE LOOP_70
        IF (AXIAL .OR. ZLM%RE*ID.GT.-NDROP.AND..NOT.ETANE0) CYCLE LOOP_70
        UNSTAB = .TRUE.
!-----------------------------------------------------------------------
! ***    take action if cannot or should not recur below this ZL:
!-----------------------------------------------------------------------
  50    ZLM = ZL
        LF = L
        IF (ID.LT.0) GOTO 380
        IF(.NOT.UNSTAB) LF = L + 1
        IF (L+MONO.LT.L1-2 .OR. ID.LT.0 .OR. .NOT.UNSTAB) GOTO 80
!-----------------------------------------------------------------------
! otherwise, all L values (for stability) should be done
!            in the reverse direction:
!-----------------------------------------------------------------------
        GOTO 60
      END DO LOOP_70
      GOTO 80
  60  ID = -1
      LF = L1
      L1 = M1
      RERR = ACCT
      GOTO 10
  80  IF (FCL .EQ. ZERO) FCL = + ACC8
      F = FPL/FCL
!-----------------------------------------------------------------------
! *** Check, if second time around, that the 'f' values agree]
!-----------------------------------------------------------------------
      IF (ID.GT.0) FIRST = F
      IF (DONEM) RERR = MAX(RERR, ABSC(F-FIRST)/ABSC(F))
      IF (.NOT.DONEM) THEN 
!-----------------------------------------------------------------------
        NOCF2 = .FALSE.
        THETAM  = X - ETA*(XLOG+TLOG) - ZLM*HPI + SIGMA
!-----------------------------------------------------------------------
! *** on left x-plane, determine OMEGA by requiring cut on -x axis
!     on right x-plane, choose OMEGA (using estimate based on THETAM)
!       so H(omega) is smaller and recurs upwards accurately.
!     (x-plane boundary is shifted to give CF2(LH) a chance to converge)
!-----------------------------------------------------------------------
        OMEGA = SIGN(ONE,X%IM+ACC8)
        IF (X%RE.GE.XNEAR) OMEGA = SIGN(ONE,THETAM%IM+ACC8)
        IF (AXIAL) OMEGA = ONE
!-----------------------------------------------------------------------
        SFSH = EXP(OMEGA*SCALE - ABS(SCALE))
        OFF = EXP(MIN(TWO*MAX(ABS(X%IM),ABS(THETAM%IM),ABS(ZLM%IM)*3),FPLMAX))
        EPS = MAX(ACC8 , ACCT * HALF / OFF)
      END IF
!-----------------------------------------------------------------------
! ***    Try first estimated omega, then its opposite,
!        to find the H(omega) linearly independent of F
!        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega)
!-----------------------------------------------------------------------
      DO L=1,2
        LH = 1
        IF (OMEGA.LT.ZERO) LH = 2
        PM = CI*OMEGA
        ETAP = ETA * PM
        IF (DONEM) GOTO 130
        PQ1 = ZERO
        PACCQ = ONE
        KASE = 0
!-----------------------------------------------------------------------
! ***   Check for small X, i.e. whether to avoid CF2 :
!-----------------------------------------------------------------------
! FLORES -- I.T. has:
! no change is test output -- keep original
!       IF(NOCF2 .OR. ABSX.LT.XNEAR .AND. &
!          ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM).LT.4) THEN
!         KASE = 5
!         GO TO 120
!       END IF
        IF (MODE.GE.3 .AND. ABSX.LT.ONE ) GOTO 190
        IF (MODE.LT.3 .AND. (NOCF2 .OR. ABSX.LT.XNEAR .AND.  &
          &  ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM).LT.4)) THEN
          KASE = 5
          GOTO 120
        END IF
!-----------------------------------------------------------------------
! ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM
!-----------------------------------------------------------------------
        PQ1 = CF2(X,ETA,ZLM,PM,EPS,LIMIT,ERROR,NPQ(LH),ACC8,ACCH,PR,ACCUR,DELL,'COULCC')
!-----------------------------------------------------------------------
        ERROR = ERROR * MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8))
        IF (ERROR.LT.ACCH) GOTO 110
!-----------------------------------------------------------------------
! *** check if impossible to get F-PQ accurately because of cancellation
!-----------------------------------------------------------------------
        NOCF2 = X%RE.LT.XNEAR .AND. ABS(X%IM).LT.-LOG(ACC8)
!-----------------------------------------------------------------------
!     original guess for OMEGA (based on THETAM) was wrong
!     Use KASE 5 or 6 if necessary if Re(X) < XNEAR
!-----------------------------------------------------------------------
        OMEGA = - OMEGA
      END DO
!-----------------------------------------------------------------------
      IF (UNSTAB) GOTO 360
      IF (X%RE.LT.-XNEAR .AND. PR) WRITE(STDOUT,1060) '-X',ERROR
  110 RERR = MAX(RERR,ERROR)
!-----------------------------------------------------------------------
! ***  establish case of calculation required for irregular solution
!-----------------------------------------------------------------------
  120 IF (KASE.GE.5) GOTO 130
!-----------------------------------------------------------------------
!  estimate errors if KASE 2 or 3 were to be used:
!-----------------------------------------------------------------------
! FLORES I.T. REMOVES THE FOLLOWING IF STATEMENT -- ALWAYS SETS PACCQ
      IF (X%RE .GT. XNEAR) THEN
        PACCQ = EPS * OFF * ABSC(PQ1) / MAX(ABS(PQ1%IM),ACC8)
      END IF
      IF (PACCQ .LT. ACCUR) THEN
        KASE = 2
        IF (AXIAL) KASE = 3
      ELSE
        KASE = 1
! FLORES I.T. HAS -- keep modification -- fixes maxed out iters in F20
        IF (NPQ(1).GT.0 .AND. NPQ(1) * R20 .LT. JMAX) KASE = 4
        !IF (NPQ(1) * R20 .LT. JMAX) KASE = 4
!-----------------------------------------------------------------------
!  i.e. change to kase=4 if the 2F0 predicted to converge
!-----------------------------------------------------------------------
      END IF
  130 SELECT CASE(ABS(KASE))
        CASE(1)
          GOTO 190
        CASE(2)
          GOTO 140
        CASE(3)
          GOTO 150
        CASE(4)
          GOTO 170
        CASE(5)
          GOTO 190
        CASE(6)
          GOTO 190
      END SELECT
  140 IF(.NOT.DONEM) PQ2=CF2(X,ETA,ZLM,-PM,EPS,LIMIT,ERROR,NPQ(3-LH),ACC8,ACCH,PR,ACCUR,DELL,'COULCC')
!-----------------------------------------------------------------------
! ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2)
!-----------------------------------------------------------------------
      P = (PQ2 + PQ1) * HALF
      Q = (PQ2 - PQ1) * HALF*PM
!-----------------------------------------------------------------------
      GOTO 160
!-----------------------------------------------------------------------
  150 P = PQ1%RE
      Q = PQ1%IM
!-----------------------------------------------------------------------
! ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*
!-----------------------------------------------------------------------
      PQ2 = CONJG(PQ1)
!-----------------------------------------------------------------------
! *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL
!-----------------------------------------------------------------------
  160 W   = (PQ1 - F) * (PQ2 - F)
      SF  = EXP(-ABS(SCALE))
      FCM = SQRT(Q / W) * SF
!-----------------------------------------------------------------------
!  any SQRT given here is corrected by
!  using sign for FCM nearest to phase of FCL
!-----------------------------------------------------------------------
      IF (REAL(FCM/FCL,KIND=dpf).LT.ZERO) FCM  = - FCM
      GAM = (F - P)/Q
      TA = ABSC(GAM + PM)
      PACCQ= EPS * MAX(TA,ONE/TA)
      HCL = FCM * (GAM + PM) * (SFSH/(SF*SF))
!-----------------------------------------------------------------------
      IF (PACCQ.GT.ACCUR .AND. KASE.GT.0) THEN
!-----------------------------------------------------------------------
!  Consider a KASE = 1 Calculation
!-----------------------------------------------------------------------
        F11V= F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERROR,N11,FPMAX,ACC8,ACC16)
        IF (ERROR.LT.PACCQ) GOTO 200
      END IF
      RERR=MAX(RERR,PACCQ)
!-----------------------------------------------------------------------
      GOTO 230
!-----------------------------------------------------------------------
! *** Arrive here if KASE = 4
!     to evaluate the exponentially decreasing H(LH) directly.
!-----------------------------------------------------------------------
  170 IF (DONEM) GOTO 180
      AA = ETAP - ZLM
      BB = ETAP + ZLM + ONE
      F20V = F20(AA,BB,-HALF*PM*XI, ACCT,JMAX,ERROR,FPMAX,N20,XRCF)
      IF (N20.LE.0) GOTO 190
      RERR = MAX(RERR,ERROR)
      HCL = FPMIN
      IF (ABS(REAL(PM*THETAM,KIND=dpf)+OMEGA*SCALE).GT.FPLMAX) GOTO 330
  180 HCL = F20V * EXP(PM * THETAM + OMEGA*SCALE)
      FCM = SFSH / ((F - PQ1) * HCL )
!-----------------------------------------------------------------------
      GOTO 230
!-----------------------------------------------------------------------
! *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)
!-----------------------------------------------------------------------
!  for small values of X, calculate F(X,SL) directly from 1F1
!      using REAL*16 arithmetic if possible.
!  where Z11 = ZLL if ID>0, or = ZLM if ID<0
!-----------------------------------------------------------------------
  190 F11V = F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERROR,N11,FPMAX,ACC8,ACC16)
!-----------------------------------------------------------------------
  200 IF (N11.LT.0) THEN
!-----------------------------------------------------------------------
!    F11 failed from BB = negative integer
!-----------------------------------------------------------------------
        WRITE(STDOUT,1060) '-L',ONE
        GOTO 390
      END IF
!-----------------------------------------------------------------------
      IF (ERROR.GT.PACCQ .AND. PACCQ.LT.ACCB) THEN
!-----------------------------------------------------------------------
!  Consider a KASE 2 or 3 calculation
!-----------------------------------------------------------------------
        KASE = -2
        IF (AXIAL) KASE = -3
        GOTO 130
      END IF
      RERR = MAX(RERR, ERROR)
      IF (ERROR.GT.FPMAX) GOTO 370
      IF (ID.LT.0) CLL = Z11*TLOG-HPI*ETA-CLOGAM(BB)+CLOGAM(Z11+ONE+P11*ETA)-P11*SIGMA
      EK = (Z11+ONE)*XLOG - P11*X + CLL  - ABS(SCALE)
      IF (ID.GT.0) EK = EK - FESL + LOG(FCL)
      IF (EK%RE.GT.FPLMAX) GOTO 350
      IF (EK%RE.LT.FPLMIN) GOTO 340
      FCM = F11V * EXP(EK)
!-----------------------------------------------------------------------
      IF (KASE.GE.5) THEN
        IF (ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)).LT.ACCH) KASE = 6
!-----------------------------------------------------------------------
! ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so
! ***      use an expansion for irregular soln from origin :
!-----------------------------------------------------------------------
        SL = ZLM
        ZLNEG = ZLM%RE .LT. -ONE + ACCB
        IF (KASE.EQ.5 .OR. ZLNEG) SL = - ZLM - ONE
        PK = SL + ONE
        AA = PK - ETAP
        AB = PK + ETAP
        BB = TWO*PK
        CLGAA = CLOGAM(AA)
        CLGAB = CLGAA
        IF (ETANE0) CLGAB = CLOGAM(AB)
        CLGBB = CLOGAM(BB)
        IF (KASE.EQ.6 .AND. .NOT.ZLNEG) THEN
          IF (NPINT(AA,ACCUR)) CLGAA = CLGAB - TWO*PM*SIGMA
          IF (NPINT(AB,ACCUR)) CLGAB = CLGAA + TWO*PM*SIGMA
        END IF
        CLL = SL*TLOG - HPI*ETA - CLGBB + (CLGAA + CLGAB) * HALF
        DSIG = (CLGAA - CLGAB) * PM*HALF
        IF (KASE.EQ.6) P11 = - PM
        EK = PK * XLOG - P11*X + CLL  - ABS(SCALE)
        SF = EXP(-ABS(SCALE))
        CHI = ZERO
        IF(.NOT.( KASE.EQ.5 .OR. ZLNEG ) ) GOTO 210
!-----------------------------------------------------------------------
! *** Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)
!-----------------------------------------------------------------------
!      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2
!-----------------------------------------------------------------------
        CHI = SIGMA - DSIG - (ZLM-SL) * HPI
        F11V=F11(X,ETA,SL,P11,ACCT,LIMIT,0,ERROR,NPQ(1),FPMAX,ACC8,ACC16)
        RERR = MAX(RERR,ERROR)
        IF (KASE.EQ.6) GOTO 210
        FESL = F11V * EXP(EK)
        FCL1 = EXP(PM*CHI) * FCM
        HCL = FCL1 - FESL
        RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL))
        HCL = HCL / SIN(CHI) * (SFSH/(SF*SF))
!-----------------------------------------------------------------------
        GOTO 220
!-----------------------------------------------------------------------
! *** Use the logarithmic expansion for the irregular solution (KASE 6)
!        for the case that BB is integral so sin(CHI) would be zero.
!-----------------------------------------------------------------------
  210   RL = BB - ONE
        N  = NINTC(RL)
        ZLOG = XLOG + TLOG - PM*HPI
        CHI = CHI + PM * THETAM + OMEGA * SCALE + AB * ZLOG
        AA  = ONE - AA
        IF (NPINT(AA,ACCUR)) THEN
          HCL = ZERO
        ELSE
          IF (ID.GT.0 .AND. .NOT.ZLNEG) F11V = FCM * EXP(-EK)
          HCL = EXP(CHI-CLGBB-CLOGAM(AA))*(-1)**(N+1)*(F11V*ZLOG &
            & +F11(X,ETA,SL,-PM,ACCT,LIMIT,2,ERROR,NPQ(2),FPMAX,ACC8,ACC16))
          RERR = MAX(RERR,ERROR)
        END IF
        IF (N.GT.0) THEN
           EK  = CHI + CLOGAM(RL) - CLGAB - RL*ZLOG
           DF  = F11(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERROR,L,FPMAX,ACC8,ACC16)
           HCL = HCL + EXP(EK) * DF
        END IF
!-----------------------------------------------------------------------
  220   PQ1 = F - SFSH/(FCM * HCL)
!-----------------------------------------------------------------------
      ELSE !KASE<=4
        IF (MODE.LE.2) HCL = SFSH/((F - PQ1) * FCM)
        KASE = 1
      END IF
!-----------------------------------------------------------------------
! ***  Now have absolute normalisations for Coulomb Functions
!          FCM & HCL  at lambda = ZLM
!      so determine linear transformations for Functions required :
!-----------------------------------------------------------------------
  230 IH = ABS(MODE1) / 10
      IF (KFN.EQ.3) IH = INT((3._dpf-CIK%IM)/2._dpf  + HALF,KIND=spi)
      P11 = ONE
      IF (IH.EQ.1) P11 = CI
      IF (IH.EQ.2) P11 = -CI
      DF = - PM
      IF (IH.GE.1) DF = - PM + P11
      IF (ABSC(DF).LT.ACCH) DF = ZERO
!-----------------------------------------------------------------------
! *** Normalisations for spherical or cylindrical Bessel functions
!-----------------------------------------------------------------------
      ALPHA = ZERO
      IF (KFN .EQ. 1) ALPHA = XI
      IF (KFN .GE. 2) ALPHA = XI*HALF
      BETA  = ONE
      IF (KFN .EQ. 1) BETA  = XI
      IF (KFN .GE. 2) BETA  = SQRT(XI/HPI)
      IF (KFN .GE. 2 .AND. BETA%RE.LT.ZERO) BETA  = - BETA
!-----------------------------------------------------------------------
      AA = ONE
      IF (KFN.GT.0) AA = -P11 * BETA
      IF (KFN.GE.3) THEN
!-----------------------------------------------------------------------
! Calculate rescaling factors for I & K output
!-----------------------------------------------------------------------
        P = EXP((ZLM+DELL) * HPI * CIK)
        AA= BETA * HPI * P
        BETA = BETA / P
        Q = CIK * ID
      END IF
!-----------------------------------------------------------------------
! Calculate rescaling factors for GC output
!-----------------------------------------------------------------------
      IF (IH.EQ.0) THEN
        TA = ABS(SCALE) + PM%IM*SCALE
        RK = ZERO
        IF (TA.LT.FPLMAX) RK = EXP(-TA)
      ELSE
        TA = ABS(SCALE) + P11%IM*SCALE
!-----------------------------------------------------------------------
        IF (ABSC(DF).GT.ACCH .AND. TA.GT.FPLMAX) GOTO 320
        IF (ABSC(DF).GT.ACCH) DF = DF * EXP(TA)
        SF = TWO * (LH-IH) * SCALE
        RK = ZERO
        IF (SF.GT.FPLMAX) GOTO 320
        IF (SF.GT.FPLMIN) RK = EXP(SF)
      END IF
!-----------------------------------------------------------------------
      KAS((3-ID)/2) = KASE
      W = FCM / FCL
      IF (LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GOTO 340
      IF (MODE.LT.3) THEN
        IF (ABSC(F-PQ1) .LT. ACCH*ABSC(F) .AND. PR) WRITE(STDOUT,1020) LH,ZLM+DELL
        HPL = HCL * PQ1
        IF (ABSC(HPL).LT.FPMIN.OR.ABSC(HCL).LT.FPMIN) GOTO 330
      END IF
!-----------------------------------------------------------------------
! *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)
! *** renormalise FC,FCP at each lambda
! ***    ZL   = ZLM - MIN(ID,0) here
!-----------------------------------------------------------------------
     DO L = LF,L1,ID
!-----------------------------------------------------------------------
        FCL = W* FC(L)
        IF (ABSC(FCL).LT.FPMIN) GOTO 340
        IF (IFCP) FPL = W*FCP(L)
        FC(L)  = BETA * FCL
        IF (IFCP) FCP(L) = BETA * (FPL - ALPHA * FCL) * CIK
        FC(L)  = TIDY(FC(L),ACCUR)
        IF (IFCP) FCP(L) = TIDY(FCP(L),ACCUR)
        IF (MODE .GE. 3) GOTO 260
        IF (L.EQ.LF)  GOTO 250
        ZL = ZL + ID
        ZID= ZL * ID
        RL = GC(L)
        IF (ETANE0) THEN
          SL = ETA + ZL*ZL*XI
          IF (MODE.EQ.1) THEN
            PL = GCP(L)
          ELSE
            PL = ZERO
            IF (ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZID
          END IF
          HCL1= (SL*HCL - ZID*HPL) / RL
          HPL = (SL*HPL - PL *HCL) / RL
        ELSE
          HCL1= RL * HCL - HPL * ID
          HPL = (HCL - RL * HCL1) * ID
        END IF
        HCL = HCL1
        IF (ABSC(HCL).GT.FPMAX) GOTO 320
  250   GC(L) = AA * (RK * HCL + DF * FCL)
        IF (MODE.EQ.1) GCP(L) = (AA*(RK*HPL+DF*FPL)-ALPHA*GC(L))*CIK
        GC(L) = TIDY(GC(L),ACCUR)
        IF (MODE.EQ.1) GCP(L) = TIDY(GCP(L),ACCUR)
        IF (KFN.GE.3) AA = AA * Q
  260   IF (KFN.GE.3) BETA = - BETA * Q
!-----------------------------------------------------------------------
        LAST = MIN(LAST,(L1-L)*ID)
      END DO
!-----------------------------------------------------------------------
! *** Come here after all soft errors to determine how many L values ok
!-----------------------------------------------------------------------
  280 IF (ID.GT.0 .OR.  LAST.EQ.0) IFAIL = LAST
      IF (ID.LT.0 .AND. LAST.NE.0) IFAIL = -3
!-----------------------------------------------------------------------
! *** Come here after ALL errors for this L range (ZLM,ZLL)
!-----------------------------------------------------------------------
  290 IF (ID.GT.0 .AND. LF.NE.M1) GOTO 300
      IF (IFAIL.LT.0) RETURN
      IF (RERR.GT.ACCB) WRITE(STDOUT,1070) RERR
      IF (RERR.GT.0.1) IFAIL = -4
!-----------------------------------------------------------------------
      RETURN
!-----------------------------------------------------------------------
! *** so on first block, 'F' started decreasing monotonically,
!                        or hit bound states for low ZL.
!     thus redo M1 to LF-1 in reverse direction
!      i.e. do CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)
!-----------------------------------------------------------------------
  300 ID = -1
      IF(.NOT.UNSTAB) LF = LF - 1
      DONEM = UNSTAB
      LF = MIN(LF,L1)
      L1 = M1
!-----------------------------------------------------------------------
      GOTO 10
!-----------------------------------------------------------------------
! ***    error messages
!-----------------------------------------------------------------------
 1000 FORMAT(/' COULCC: CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',1P,2D10.2,', AS ABS(X) IS TOO SMALL'/)
 1010 FORMAT(' COULCC: AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P,2E10.1,') WILL BE ',A4,' THAN',E10.1)
 1020 FORMAT(' COULCC WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''H(',I1,')'' &
             & IS LOST AT ZL =',2F7.2,' (EG. COULOMB EIGENSTATE, OR CF1 UNSTABLE)'/)
 1030 FORMAT(' COULCC: (ETA&L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE ATL =',2F8.2)
 1040 FORMAT(' COULCC: OVERFLOW IN 1F1 SERIES AT ZL =',2F8.3,' AT TERM',I5)
 1050 FORMAT(' COULCC: BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR'&
        &  ,', OR MULTIPLE INSTABILITIES PRESENT.'                         &
        & ,/,' TRY CALLING TWICE,  FIRST FOR ZL FROM',2F8.3,' TO',2F8.3,   &
        &  ' (INCL.)',/,20X,     'SECOND FOR ZL FROM',2F8.3,' TO',2F8.3)
 1060 FORMAT(' COULCC WARNING: AS ''',A2,''' REFLECTION RULES NOT USED,ERRORS CAN BE UP TO',1P,D12.2/)
 1070 FORMAT(' COULCC WARNING: OVERALL ROUNDOFF ERROR APPROX.',1P,E11.1)
!-----------------------------------------------------------------------
  310 IF (PR) WRITE(STDOUT,1000) XX
      RETURN
  320 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'IR',HCL,'MORE',FPMAX
      GOTO 280
  330 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'IR',HCL,'LESS',FPMIN
      GOTO 280
  340 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'  ',FCL,'LESS',FPMIN
      GOTO 280
  350 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'  ',FCL,'MORE',FPMAX
      GOTO 280
  360 IF (PR) WRITE(STDOUT,1030) ZLL+DELL
      GOTO 280
  370 IF (PR) WRITE(STDOUT,1040) Z11,I
      GOTO 390
  380 IF (PR) WRITE(STDOUT,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL-ONE
  390 IFAIL = -1
      GOTO 290
!-----------------------------------------------------------------------
    END SUBROUTINE COULCC
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CF1C(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERROR,NFP &
                            &, ACCH,FPMIN,FPMAX,PR,CALLER)
!-----------------------------------------------------------------------
! *** Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
!     using complex arithmetic
!-----------------------------------------------------------------------
      COMPLEX(dpf),    INTENT(IN)  :: X,ETA,ZL
      COMPLEX(dpf),    INTENT(OUT) :: TPK1,FCL
      INTEGER(spi),    INTENT(IN)  :: LIMIT
      INTEGER(spi),    INTENT(OUT) :: NFP
      REAL(dpf),       INTENT(OUT) :: ERROR
      REAL(dpf),       INTENT(IN)  :: EPS,ACCH,FPMIN,FPMAX
      LOGICAL,         INTENT(IN)  :: PR,ETANE0
      CHARACTER(len=6),INTENT(IN)  :: CALLER
!-----------------------------------------------------------------------
      REAL(dpf), PARAMETER :: ONE=1._dpf
      REAL(dpf), PARAMETER :: TWO=2._dpf
!-----------------------------------------------------------------------
      REAL(dpf) :: RK,PX,SMALL
      COMPLEX(dpf) :: XI,PK,EK,RK2,F,PK1,TK,SL,D,DF
      LOGICAL :: CONVERGED
!-----------------------------------------------------------------------
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',/1X,1P,13D9.2/)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATIONS AS ABS(X) =',F15.0)
!-----------------------------------------------------------------------
      SL = CMPLX(0._dpf,0._dpf,KIND=dpf)
!-----------------------------------------------------------------------
      FCL = ONE
      XI  = ONE/X
      PK  = ZL+ONE
      PX  = PK%RE+LIMIT
!-----------------------------------------------------------------------
! ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.
!-----------------------------------------------------------------------
      IF (ETANE0) THEN
        LOOP_10: DO
          EK  = ETA / PK
          RK2 = ONE + EK*EK
          F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
          PK1 =  PK + ONE
          TPK1 = PK + PK1
          TK  = TPK1*(XI + EK/PK1)
          IF (ABSC(TK) .GT. ACCH) EXIT LOOP_10
          FCL  = RK2/(ONE + (ETA/PK1)**2)
          SL   = TPK1*XI * (TPK1+TWO)*XI
          PK   =  TWO + PK
        END DO LOOP_10
      ELSE
        EK  = ETA / PK
        RK2 = ONE + EK*EK
        F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
        PK1 =  PK + ONE
        TPK1 = PK + PK1
        TK  = TPK1*(XI + EK/PK1)
      END IF
!-----------------------------------------------------------------------
      D  =  ONE/TK
      DF = -FCL*RK2*D
      IF (PK%RE.GT.ZL%RE+TWO) FCL = - RK2 * SL
      FCL = FCL * D * TPK1 * XI
      F   =  F  + DF
!-----------------------------------------------------------------------
! *** begin CF1 loop on PK = k = lambda + 1
!-----------------------------------------------------------------------
      RK    = ONE
      SMALL = SQRT(FPMIN)
      CONVERGED=.FALSE.
!-----------------------------------------------------------------------
      LOOP_30: DO
        PK    = PK1
        PK1 = PK1 + ONE
        TPK1 = PK + PK1
        IF (ETANE0) THEN
          EK  = ETA / PK
          RK2 = ONE + EK*EK
        END IF
        TK = TPK1*(XI + EK/PK1)
        D  =  TK - D*RK2
        IF (ABSC(D) .LE. ACCH) THEN
          IF (PR) WRITE (STDOUT,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X
          RK= RK +   ONE
          IF( RK .GT. TWO) EXIT LOOP_30
        END IF
        D = ONE/D
        FCL = FCL * D * TPK1*XI
        IF (ABSC(FCL).LT.SMALL) FCL = FCL / SMALL
        IF (ABSC(FCL).GT.FPMAX) FCL = FCL / FPMAX
        DF = DF*(D*TK - ONE)
        F  = F  + DF
        IF (PK%RE .GT. PX) EXIT LOOP_30
        IF (ABSC(DF) .LT. ABSC(F)*EPS) THEN
          CONVERGED=.TRUE.
          EXIT LOOP_30
        END IF
      END DO LOOP_30
!-----------------------------------------------------------------------
      IF (CONVERGED) THEN
        NFP = INT(PK%RE - ZL%RE,KIND=spi) - 1_spi
        ERROR = EPS * SQRT(REAL(NFP,KIND=dpf))
        CF1C = F
      ELSE
        IF (PR) WRITE (STDOUT,1010) CALLER,LIMIT,ABS(X)
        ERROR = TWO
        CF1C = CMPLX(0._dpf,0._dpf,KIND=dpf)
      END IF
!-----------------------------------------------------------------------
    END FUNCTION CF1C
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CF2(X,ETA,ZL,PM,EPS,LIMIT,ERROR,NPQ,ACC8,ACCH  &
                            &, PR,ACCUR,DELL,CALLER)
!-----------------------------------------------------------------------
!                                    (omega)        (omega)
! *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)
!                                    ZL             ZL
!     where PM = omega.i
!-----------------------------------------------------------------------
      COMPLEX(dpf),    INTENT(IN)  :: X,ETA,ZL,PM,DELL
      INTEGER(spi),    INTENT(IN)  :: LIMIT
      INTEGER(spi),    INTENT(OUT) :: NPQ
      REAL(dpf),       INTENT(OUT) :: ERROR
      REAL(dpf),       INTENT(IN)  :: EPS,ACC8,ACCH,ACCUR
      LOGICAL,         INTENT(IN)  :: PR
      CHARACTER(len=6),INTENT(IN)  :: CALLER
!-----------------------------------------------------------------------
      REAL(dpf), PARAMETER :: ZERO=0._dpf
      REAL(dpf), PARAMETER :: HALF=0.5_dpf
      REAL(dpf), PARAMETER :: ONE=1._dpf
      REAL(dpf), PARAMETER :: TWO=2._dpf
!-----------------------------------------------------------------------
      REAL(dpf) :: TA,RK
      COMPLEX(dpf) :: E2MM1,ETAP,XI,WI,PQ,AA,BB,RL,DD,DL
!-----------------------------------------------------------------------
      TA = TWO*LIMIT
      E2MM1 = ETA*ETA + ZL*ZL + ZL
      ETAP = ETA * PM
      XI = ONE/X
      WI = TWO*ETAP
      RK = ZERO
      PQ = (ONE - ETA*XI) * PM
      AA = -E2MM1 + ETAP
      BB = TWO*(X - ETA + PM)
      RL = XI * PM
      IF (ABSC(BB).LT.ACCH) THEN
        RL = RL * AA / (AA + RK + WI)
        PQ = PQ + RL * (BB + TWO*PM)
        AA = AA + TWO*(RK+ONE+WI)
        BB = BB + (TWO+TWO)*PM
        RK = RK + (TWO+TWO)
      END IF
      DD = ONE/BB
      DL = AA*DD* RL
      ERROR = HUGE(1._dpf)
      DO WHILE (ERROR.GE.MAX(EPS,ACC8*RK*HALF) .AND. RK.LE.TA)
        PQ = PQ + DL
        RK = RK + TWO
        AA = AA + RK + WI
        BB = BB + TWO*PM
        DD = ONE/(AA*DD + BB)
        DL = DL*(BB*DD - ONE)
        ERROR = ABSC(DL)/ABSC(PQ)
      END DO
!-----------------------------------------------------------------------
      NPQ = INT(RK/TWO,KIND=spi)
      PQ  = PQ + DL
      CF2 = PQ
!-----------------------------------------------------------------------
 1000 FORMAT(' ',A6,': CF2(',I2,') NOT CONVERGED FULLY IN ',I7,         &
     & ' ITERATIONS, SO ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL &
     & =', 0P,2F8.3)
!-----------------------------------------------------------------------
      IF (PR.AND.NPQ.GE.LIMIT-1 .AND. ERROR.GT.ACCUR) WRITE(STDOUT,1000) CALLER,INT(PM%IM),NPQ,ERROR,ZL+DELL
!-----------------------------------------------------------------------
    END FUNCTION CF2
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION F11(X,ETA,ZL,P,EPS,LIMIT,KIND,ERROR,NITS,FPMAX,ACC8,ACC16)
!-----------------------------------------------------------------------
! *** evaluate the HYPERGEOMETRIC FUNCTION 1F1
!                                        i
!            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i] )
!           1 1              i       i            i
!
!     to accuracy EPS with at most LIMIT terms.
!  If KIND = 0 : using extended precision but real arithmetic only,
!            1 : using normal precision in complex arithmetic,
!   or       2 : using normal complex arithmetic, but with CDIGAM factor
!-----------------------------------------------------------------------
      COMPLEX(dpf),INTENT(IN)  :: X,ETA,ZL,P
      REAL(dpf),   INTENT(IN)  :: EPS,FPMAX,ACC8,ACC16
      REAL(dpf),   INTENT(OUT) :: ERROR
      INTEGER(spi),INTENT(OUT) :: NITS
      INTEGER(spi),INTENT(IN)  :: LIMIT,KIND
!-----------------------------------------------------------------------
      COMPLEX(dpf) :: AA,BB,Z,DD,G,F,AI,BI,T
      REAL(dpf)    :: R,RK,TA
      INTEGER(spi) :: I,EXIT_COND
      LOGICAL      :: ZLLIN,ZLLIN2
!-----------------------------------------------------------------------
      REAL(qpf) :: AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN !128 bit real
!-----------------------------------------------------------------------
      REAL(dpf),   PARAMETER :: ZERO=0._dpf
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: TWO=2._dpf
      COMPLEX(dpf),PARAMETER :: CI=CMPLX(0._dpf,1._dpf,KIND=dpf)
!-----------------------------------------------------------------------
      DR = 0._qpf
      DI = 0._qpf
!-----------------------------------------------------------------------
!  FUNCTION INPUTS
!-----------------------------------------------------------------------
      AA = ZL+ONE - ETA*P
      BB = TWO*(ZL+ONE)
      Z  = TWO*P*X
!-----------------------------------------------------------------------
      EXIT_COND=0_spi
!-----------------------------------------------------------------------
      ZLLIN = BB%RE.LE.ZERO .AND. ABS(BB-CMPLX(NINTC(BB),0._dpf,KIND=dpf)).LT.ACC8**0.25_dpf
      ZLLIN2= BB%RE+LIMIT.LT.1.5_dpf
      IF(.NOT.(.NOT.ZLLIN.OR.ZLLIN2)) EXIT_COND=-1
!-----------------------------------------------------------------------
      IF (LIMIT.LE.0) EXIT_COND=1
!-----------------------------------------------------------------------
      IF (.NOT.(EXIT_COND.EQ.-1 .OR. EXIT_COND.EQ.1)) THEN
!-----------------------------------------------------------------------
        TA = ONE
        RK = ONE
        IF (KIND.LE.0 .AND. ABSC(Z)*ABSC(AA).GT.ABSC(BB)*1.0) THEN
          DR = ONE
          DI = ZERO
          GR = ONE
          GI = ZERO
          AR = AA%RE
          BR = BB%RE
          FI = ZERO
          LOOP_20: DO I=2,LIMIT
            FI1 = FI + ONE
            TR = BR * FI1
            TI = BB%IM * FI1
            DEN= ONE / (TR*TR + TI*TI)
            UR = (AR*TR + AA%IM*TI) * DEN
            UI = (AA%IM*TR - AR*TI) * DEN
            TR = UR*GR - UI*GI
            TI = UR*GI + UI*GR
            GR = Z%RE * TR - Z%IM*TI
            GI = Z%RE * TI + Z%IM*TR
            DR = DR + GR
            DI = DI + GI
            ERROR = REAL(ABS(GR) + ABS(GI),KIND=dpf)
            IF (ERROR.GT.FPMAX) THEN
              EXIT_COND=-99
              EXIT LOOP_20
            END IF
            RK  = REAL(ABS(DR) + ABS(DI),KIND=dpf)
            TA = MAX(TA,RK)
            IF (ERROR.LT.RK*EPS .OR. I.GE.4.AND.ERROR.LT.ACC16) THEN
              EXIT_COND=2
              EXIT LOOP_20
            END IF
            FI = FI1
            AR = AR + ONE
            BR = BR + ONE
          END DO LOOP_20
!-----------------------------------------------------------------------
          IF (EXIT_COND.NE.-99) EXIT_COND=2
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
!* -------------------alternative code----------------------------------
!-----------------------------------------------------------------------
!*    If REAL*16 arithmetic is not available, (or already using it]),
!*    then use KIND > 0
!-----------------------------------------------------------------------
          G = ONE
          F = ONE
          IF (KIND.GE.2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)
          DD = F
          LOOP_40: DO I=2,LIMIT
            AI = AA + (I-2)
            BI = BB + (I-2)
            R  = I-ONE
            G = G * Z * AI / (BI * R)
!-----------------------------------------------------------------------
!  multiply by (psi(a+r)-psi(b+r)-psi(1+r))
!-----------------------------------------------------------------------
            IF (KIND.GE.2) F = F + ONE/AI - ONE/BI - ONE/R
            T  = G * F
            DD = DD + T
            ERROR = ABSC(T)
            IF (ERROR.GT.FPMAX) THEN
              EXIT_COND=-99
              EXIT LOOP_40
            END IF
            RK = ABSC(DD)
            TA = MAX(TA,RK)
            IF (ERROR.LT.RK*EPS.OR.ERROR.LT.ACC8.AND.I.GE.4) THEN
              EXIT_COND=3
              EXIT LOOP_40
            END IF
          END DO LOOP_40
          IF (EXIT_COND.NE.-99) EXIT_COND=3
!* --------------------------------------------- end of alternative code
        END IF
      END IF
!-----------------------------------------------------------------------
      SELECT CASE(EXIT_COND)
        CASE(-1)
          F11 = ZERO
          NITS = -1
        CASE(1)
          F11 = ZERO
          ERROR = ZERO
          NITS= 1
        CASE(2)
          F11 = CMPLX(DR + CI * DI,KIND=dpf)
          ERROR = ACC16 * TA / RK
          NITS = I
        CASE(3)
          ERROR = ACC8 * TA / RK
          F11 = DD
          NITS = I
        CASE(-99)
          F11 = ZERO
          NITS = I
      END SELECT
!-----------------------------------------------------------------------
    END FUNCTION F11
!-----------------------------------------------------------------------
    REAL(dpf) FUNCTION CF1R(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERROR,NFP &
                         &,  ACCH,FPMIN,FPMAX,PR,CALLER)
!-----------------------------------------------------------------------
! ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
!        using real arithmetic
!-----------------------------------------------------------------------
      REAL(dpf),       INTENT(IN)     :: X,ETA,ZL,EPS
      REAL(dpf),       INTENT(INOUT)  :: FCL,TPK1,ERROR
      REAL(dpf),       INTENT(IN)     :: ACCH,FPMIN,FPMAX
      INTEGER(spi),    INTENT(IN)     :: LIMIT
      INTEGER(spi),    INTENT(INOUT)  :: NFP
      LOGICAL,         INTENT(IN)     :: PR,ETANE0
      CHARACTER(len=6),INTENT(IN)     :: CALLER
!-----------------------------------------------------------------------
      REAL(dpf),PARAMETER :: ONE=1._dpf
      REAL(dpf),PARAMETER :: TWO=2._dpf
!-----------------------------------------------------------------------
      REAL(dpf) :: XI,PK,PX,EK,RK2,F,PK1,TK,D,DF,RK,SMALL,SL
      LOGICAL :: CONVERGED
!-----------------------------------------------------------------------
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',/1X,1P,7D9.2/)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATIONS AS ABS(X) =',F15.0)
!-----------------------------------------------------------------------
      SL=0._dpf
!-----------------------------------------------------------------------
      FCL = ONE
      XI  = ONE/X
      PK  = ZL + ONE
      PX  = PK  + LIMIT
!-----------------------------------------------------------------------
! ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.
!-----------------------------------------------------------------------
      IF (ETANE0) THEN
!-----------------------------------------------------------------------
        LOOP_10: DO
          EK  = ETA / PK
          RK2 = ONE + EK*EK
          F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
          PK1 =  PK + ONE
          TPK1 = PK + PK1
          TK  = TPK1*(XI + EK/PK1)
!-----------------------------------------------------------------------
          IF (ABS(TK) .GT. ACCH) EXIT LOOP_10
!-----------------------------------------------------------------------
          FCL  = RK2/(ONE + (ETA/PK1)**2)
          SL   = TPK1*XI * (TPK1+TWO)*XI
          PK   =  TWO + PK
        END DO LOOP_10
!-----------------------------------------------------------------------
      ELSE
        EK  = ETA / PK
        RK2 = ONE + EK*EK
        F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
        PK1 =  PK + ONE
        TPK1 = PK + PK1
        TK  = TPK1*(XI + EK/PK1)
      END IF
!-----------------------------------------------------------------------
      D =  ONE/TK
      DF= -FCL*RK2*D
      IF (PK.GT.ZL+TWO) FCL = - RK2 * SL
      FCL=FCL * D * TPK1 * XI
      F  = F  + DF
!-----------------------------------------------------------------------
! ***   begin CF1 loop on PK = k = lambda + 1
!-----------------------------------------------------------------------
      RK    = ONE
      SMALL = SQRT(FPMIN)
!-----------------------------------------------------------------------
      CONVERGED=.FALSE.
!-----------------------------------------------------------------------
      LOOP_30: DO
        PK    = PK1
        PK1   = PK1 + ONE
        TPK1 = PK + PK1
        IF (ETANE0) THEN
          EK  = ETA / PK
          RK2 = ONE + EK*EK
        END IF
        TK=TPK1*(XI + EK/PK1)
        D = TK - D*RK2
        IF (ABS(D) .LE. ACCH) THEN
          IF (PR) WRITE(STDOUT,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X
          RK= RK + ONE
          IF (RK .GT. TWO) EXIT LOOP_30
        END IF
        D = ONE/D
        FCL = FCL * D * TPK1*XI
        IF (ABS(FCL).LT.SMALL) FCL = FCL / SMALL
        IF (ABS(FCL).GT.FPMAX) FCL = FCL / FPMAX
        DF  = DF*(D*TK - ONE)
        F   = F  + DF
        IF (PK .GT. PX) EXIT LOOP_30
        IF (ABS(DF) .LT. ABS(F)*EPS) THEN
          CONVERGED=.TRUE.
          EXIT LOOP_30
        END IF
      END DO LOOP_30
!-----------------------------------------------------------------------
      IF (CONVERGED) THEN
        NFP = INT(PK - ZL,KIND=spi) - 1_spi
        ERROR = EPS * SQRT(REAL(NFP,KIND=dpf))
        CF1R = F
      ELSE
        IF (PR) WRITE(STDOUT,1010) CALLER,LIMIT,ABS(X)
        ERROR = TWO
        CF1R = 0._dpf
      END IF
!-----------------------------------------------------------------------
    END FUNCTION CF1R
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION F20(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)
!-----------------------------------------------------------------------
!     evaluate the HYPERGEOMETRIC FUNCTION 2F0
!-----------------------------------------------------------------------
!            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i]
!           2 0              i       i     i
!-----------------------------------------------------------------------
!     to accuracy EPS with at most JMAX terms.
!-----------------------------------------------------------------------
!     if the terms start diverging,
!     the corresponding continued fraction is found by RCF
!     & evaluated progressively by Steed's method to obtain convergence.
!-----------------------------------------------------------------------
!      useful number also input:  FPMAX = near-largest f.p. number
!-----------------------------------------------------------------------
      COMPLEX(dpf),                  INTENT(IN)    :: AA,BB,Z
      INTEGER(spi),                  INTENT(IN)    :: JMAX
      COMPLEX(dpf),DIMENSION(JMAX,4),INTENT(INOUT) :: X
      INTEGER(spi),                  INTENT(OUT)   :: N
      REAL(dpf),                     INTENT(INOUT) :: RE
      REAL(dpf),                     INTENT(IN)    :: FPMAX
!-----------------------------------------------------------------------
      REAL(dpf), PARAMETER :: ZERO=0._dpf
      REAL(dpf), PARAMETER :: ONE=1._dpf
!-----------------------------------------------------------------------
      COMPLEX(dpf) :: SUM,F,D,DF
      REAL(dpf)    :: EP,EPS,AT,ATL
      INTEGER(spi) :: I,J,K,MA,MB,IMAX,IERROR,EXIT_COND
      LOGICAL      :: FINITE
!-----------------------------------------------------------------------
      EXIT_COND=20 !DEFAULT
!-----------------------------------------------------------------------
      RE = 0.0
      X(1,1) = ONE
      SUM = X(1,1)
      AT=ZERO
      ATL = ABSC(X(1,1))
      F = SUM
      D = ONE
      DF= SUM
      J = 0
      EP = EPS * JMAX*10._dpf
      MA = - NINTC(AA)
      MB = - NINTC(BB)
      FINITE = ABS(ABS(AA%RE)-MA).LT.EP .AND. ABS(AA%IM).LT.EP &
     &    .OR. ABS(ABS(BB%RE)-MB).LT.EP .AND. ABS(BB%IM).LT.EP
      IMAX = JMAX
      IF (FINITE.AND.MA.GE.0) IMAX = MIN(MA+1,IMAX)
      IF (FINITE.AND.MB.GE.0) IMAX = MIN(MB+1,IMAX)
      LOOP_10: DO I=2,IMAX
        X(I,1) = X(I-1,1) * Z * (AA+I-2) * (BB+I-2) / (I-1)
        IF (ABSC(X(I,1)).GT.FPMAX) THEN
          EXIT_COND=-1
          EXIT LOOP_10
        END IF
        AT = ABSC(X(I,1))
        IF (J.EQ.0) THEN
          SUM = SUM + X(I,1)
          IF (AT.LT.ABSC(SUM)*EPS) EXIT LOOP_10
        END IF
        IF (FINITE) THEN
          ATL = AT
          CYCLE LOOP_10
        END IF
        IF (J.GT.0 .OR. AT.GT.ATL .OR. I.GE.JMAX-2) J = J + 1
        IF (J.EQ.0) THEN
          ATL = AT
          CYCLE LOOP_10
        END IF
        CALL RCF(X(1,1),X(1,2),J,I,X(1,3),EPS,IERROR)
        IF (IERROR.LT.0) THEN
          EXIT_COND=-1
          EXIT LOOP_10
        END IF
        DO K=MAX(J,2),I
          D = ONE/(D*X(K,2) + ONE)
          DF = DF*(D - ONE)
          F = F + DF
          IF (ABSC(DF) .LT. ABSC(F)*EPS) THEN
            EXIT_COND=0
            EXIT LOOP_10
          END IF
          IF (DF.EQ.ZERO.AND.F.EQ.ZERO.AND.I.GE.4) THEN
            EXIT_COND=0
            EXIT LOOP_10
          END IF
        END DO
        J = I
        ATL = AT
      END DO LOOP_10
!-----------------------------------------------------------------------
      IF(.NOT.FINITE.AND.I.EQ.IMAX+1) I = -JMAX
!-----------------------------------------------------------------------
      SELECT CASE(EXIT_COND)
        CASE(0)
          F20 = F
          RE = ABSC(DF) / ABSC(F)
          N = K
        CASE(-1) !BAD_EXIT
          N = 0
          F20 = SUM
          IF(.NOT.FINITE) RE  = AT / ABSC(SUM)
        CASE DEFAULT
          N = I
          F20 = SUM
          IF(.NOT.FINITE) RE  = AT / ABSC(SUM)
      END SELECT
!-----------------------------------------------------------------------
    END FUNCTION F20
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CF1A(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
!-----------------------------------------------------------------------
!     evaluate the ASYMPTOTIC EXPANSION for the
!            LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION
!-----------------------------------------------------------------------
! ***        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)
!-----------------------------------------------------------------------
!      that is valid for REAL(RHO)>0, and best for RHO >> ETA**2, XL,
!      and is derived from the 2F0 expansions for H+ and H-
!      e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)
!      Some lines of this subprogram are for convenience copied from
!           Takemasa, Tamura & Wolter CPC 17 (1979) 351.
!-----------------------------------------------------------------------
!     Evaluate to accuracy EPS with at most NMAX terms.
!-----------------------------------------------------------------------
!     If the terms start diverging,
!     the corresponding continued fraction is found by RCF
!     & evaluated progressively by Steed's method to obtain convergence.
!-----------------------------------------------------------------------
!      useful number also input:  FPMAX = near-largest f.p. number
!-----------------------------------------------------------------------
      COMPLEX(dpf),                   INTENT(IN)    :: RHO,ETA,XL,PSI
      COMPLEX(dpf),                   INTENT(INOUT) :: FCL
      INTEGER(spi),                   INTENT(IN)    :: NMAX
      INTEGER(spi),                   INTENT(OUT)   :: NUSED
      COMPLEX(dpf), DIMENSION(2,NMAX),INTENT(INOUT) :: XX
      COMPLEX(dpf), DIMENSION(NMAX)  ,INTENT(INOUT) :: G,C
      REAL(dpf),                      INTENT(OUT)   :: RE
      REAL(dpf),                      INTENT(IN)    :: EPS,FPMAX
!-----------------------------------------------------------------------
      REAL(dpf),    PARAMETER :: ZERO=0._dpf
      REAL(dpf),    PARAMETER :: ONE=1._dpf
      REAL(dpf),    PARAMETER :: TWO=2._dpf
!-----------------------------------------------------------------------
      COMPLEX(dpf) :: GLAST,GSUM,XLL1,SL1,SL2,SL,SC1,SC,TL1,TL,TC1,TC
      COMPLEX(dpf) :: SC2,TC2,TL2,F,D,DF,COSL,TANL,C1,C2,DENOM,ETASQ
      REAL(dpf)    :: T1,T2,T3,AT,ATL
      INTEGER(spi) :: K,N,J,IERROR
      LOGICAL :: BAD_EXIT
!-----------------------------------------------------------------------
      BAD_EXIT=.FALSE.
!-----------------------------------------------------------------------
      K=0_spi
!-----------------------------------------------------------------------
      T1 = SIN(PSI%RE)
      T2 = COS(PSI%RE)
      ATL= TANH(PSI%IM)
!-----------------------------------------------------------------------
! GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGN
!-----------------------------------------------------------------------
      COSL = CMPLX(T2,-T1*ATL,KIND=dpf)
      TANL = CMPLX(T1, T2*ATL,KIND=dpf) / COSL
!-----------------------------------------------------------------------
      RE = ZERO
      XLL1= XL*(XL+ONE)
      ETASQ = ETA*ETA
      SL1=ONE
      SL=SL1
      SC1=ZERO
      SC=SC1
      TL1=SC
      TL=TL1
      TC1=ONE-ETA/RHO
      TC=TC1
      FCL  = TL + SL*TANL
      G(1) = (TC + SC*TANL) / FCL
      GLAST = G(1)
      ATL = ABSC(GLAST)
      F   = GLAST
      D   = ONE
      DF  = GLAST
      J = 0
!-----------------------------------------------------------------------
      LOOP_10: DO N=2,NMAX
        T1=N-1
        T2=TWO*T1-ONE
        T3=T1*(T1-ONE)
        DENOM=TWO*RHO*T1
        C1=(ETA*T2)/DENOM
        C2=(ETASQ+XLL1-T3)/DENOM
        SL2=C1*SL1-C2*TL1
        TL2=C1*TL1+C2*SL1
        SC2=C1*SC1-C2*TC1-SL2/RHO
        TC2=C1*TC1+C2*SC1-TL2/RHO
        SL=SL+SL2
        TL=TL+TL2
        SC=SC+SC2
        TC=TC+TC2
        SL1=SL2
        TL1=TL2
        SC1=SC2
        TC1=TC2
        FCL  =  TL + SL*TANL
        IF (ABSC(FCL).GT.FPMAX .OR. ABSC(FCL).LT.1./FPMAX) THEN
          BAD_EXIT=.TRUE.
          EXIT LOOP_10
        END IF
        GSUM = (TC + SC*TANL) / FCL
        G(N) = GSUM - GLAST
        GLAST = GSUM
        AT = ABSC(G(N))
        IF (AT.LT.ABSC(GSUM)*EPS) THEN
          FCL = FCL * COSL
          CF1A = GSUM
          RE = AT / ABSC(GSUM)
          NUSED = N
          RETURN
        END IF
        IF (J.GT.0 .OR. AT.GT.ATL .OR. N.GE.NMAX-2) J = J + 1
        IF (J.EQ.0) THEN
          ATL=AT
          CYCLE LOOP_10
        END IF
        CALL RCF(G,C,J,N,XX,EPS,IERROR)
        IF (IERROR.LT.0) THEN
          BAD_EXIT=.TRUE.
          EXIT LOOP_10
        END IF
        DO K=MAX(J,2),N
          D = ONE/(D*C(K) + ONE)
          DF = DF*(D - ONE)
          F = F + DF
          IF (ABSC(DF) .LT. ABSC(F)*EPS) EXIT LOOP_10
          IF (DF.EQ.ZERO.AND.F.EQ.ZERO.AND.N.GE.4) EXIT LOOP_10
        END DO
        J = N
        ATL = AT
      END DO LOOP_10
!-----------------------------------------------------------------------
      IF (N.EQ.NMAX+1) K = -NMAX
!-----------------------------------------------------------------------
      IF (BAD_EXIT) THEN
        CF1A = G(1)
        FCL = 1.0
        RE = 1.0
        NUSED = 0
      ELSE
        CF1A = F
        FCL = FCL * COSL
        RE = ABSC(DF) / ABSC(F)
        NUSED = K
      END IF
!-----------------------------------------------------------------------
    END FUNCTION CF1A
!-----------------------------------------------------------------------
    PURE LOGICAL FUNCTION NPINT(Z,ACCB)
      COMPLEX(dpf),INTENT(IN) :: Z
      REAL(dpf),   INTENT(IN) :: ACCB
      REAL(dpf), PARAMETER :: HALF=0.5_dpf
      NPINT= ABSC(CMPLX(NINTC(Z),0._dpf,KIND=dpf)-Z).LT.ACCB .AND. Z%RE.LT.HALF
    END FUNCTION NPINT
!-----------------------------------------------------------------------
    PURE INTEGER(spi) FUNCTION NINTC(Z)
      !!integer nearest to a complex no.
      COMPLEX(dpf),INTENT(IN) :: Z
      NINTC = NINT(Z%RE)
    END FUNCTION NINTC
!-----------------------------------------------------------------------
    PURE REAL(dpf) FUNCTION ABSC(Z)
      COMPLEX(dpf),INTENT(IN) :: Z
      ABSC = ABS(Z%RE) + ABS(Z%IM)
    END FUNCTION ABSC
!-----------------------------------------------------------------------
    PURE COMPLEX(dpf) FUNCTION TIDY(Z,ACC)
      !!tidy a complex number
      COMPLEX(dpf),INTENT(IN) :: Z
      REAL(dpf),   INTENT(IN) :: ACC
!-----------------------------------------------------------------------
      REAL(dpf) :: X,Y,AZ
      REAL(dpf),PARAMETER :: ZERO=0._dpf
!-----------------------------------------------------------------------
      X = Z%RE
      Y = Z%IM
      AZ= (ABS(X) + ABS(Y)) * ACC * 5
      IF (ABS(X) .LT. AZ) X = ZERO
      IF (ABS(Y) .LT. AZ) Y = ZERO
      TIDY = CMPLX(X,Y,KIND=dpf)
!-----------------------------------------------------------------------
    END FUNCTION TIDY
!-----------------------------------------------------------------------
END MODULE COULCC_M