!-------------------------------------< Green >-------------------------------------!
!                                                                                   !    
! Contains: Green function for the equation [E-T_N-V]|psi> = |rightside>            !                                                                                     
!                                                                                   !
! Last revision:    04/02/2026                                                      !
!                                                                                   !
!-----------------------------------------------------------------------------------!

    
MODULE XXGreen
    USE Math
    USE Parameters
    USE COULCC_M, only: COULCC
    USE whittaker_w, only: coulomb_whittaker
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: apply_green, apply_green_coulomb_bound, apply_green_analytic, general_asymptotic, analytic_asymptotic, general_free_solution

    
    INTERFACE apply_green
        MODULE PROCEDURE apply_green_real_vector
        MODULE PROCEDURE apply_green_real_matrix
        MODULE PROCEDURE apply_green_complex_vector
        MODULE PROCEDURE apply_green_complex_matrix
    END INTERFACE apply_green

    INTERFACE apply_green_coulomb_bound
        MODULE PROCEDURE apply_green_coulomb_bound_real_vector
        MODULE PROCEDURE apply_green_coulomb_bound_real_matrix
    END INTERFACE apply_green_coulomb_bound

    INTERFACE apply_green_analytic
        MODULE PROCEDURE apply_green_analytic_real_vector
        MODULE PROCEDURE apply_green_analytic_real_matrix
        MODULE PROCEDURE apply_green_analytic_complex_vector
        MODULE PROCEDURE apply_green_analytic_complex_matrix
    END INTERFACE apply_green_analytic


    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_green_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: x
        END FUNCTION potential_green_interface
    END INTERFACE

    CONTAINS

    ! Function to compute a generalized exp(i*k*x) asymptotics for higher partial waves and/or Coulomb potential. For E>0 it always returns H+ (outgoing wave), for E<0 it returns "real" decay function (Whittaker W for Coulomb, rotated H+ for free case).
    FUNCTION general_asymptotic(x, E, m, Z, l) RESULT(val)
        REAL(KIND=idk), INTENT(IN) :: x, E, m, Z
        INTEGER, INTENT(IN) :: l
        COMPLEX(KIND=idk) :: val
        
        REAL(KIND=idk) :: k, eta_real, kappa
        COMPLEX(KIND=idk) :: zrho, zeta, zl_min
        COMPLEX(KIND=idk) :: fc(1), gc(1), fcp(1), gcp(1), sig(1)
        REAL(KIND=idk) :: w_whit, wd_whit
        INTEGER :: ifail, sf_whit
        INTEGER :: kfn, mode
        
        zl_min = CMPLX(l, 0.0d0, KIND=idk)
        mode = 12
        
        IF (E >= 0.0d0) THEN
            k = SQRT(2.0d0 * m * E) / hbar
            eta_real = -Z * m / (hbar**2 * k)         
            zrho = CMPLX(k * x, 0.0d0, KIND=idk)
            zeta = CMPLX(eta_real, 0.0d0, KIND=idk)
            kfn = 0
            ifail = 0        
            CALL COULCC(zrho, zeta, zl_min, 1, fc, gc, fcp, gcp, sig, mode, kfn, ifail)
            val = gc(1)
            
        ELSE
            kappa = SQRT(2.0d0 * m * ABS(E)) / hbar    
            IF (ABS(Z) > 1.0d-10) THEN
                eta_real = -Z * m / (hbar**2 * kappa) 
                CALL coulomb_whittaker(eta_real, l, kappa * x, w_whit, wd_whit, sf_whit)
                val = CMPLX(w_whit * (10.0d0**sf_whit), 0.0d0, KIND=idk)
            ELSE
                zrho = CMPLX(0.0d0, kappa * x, KIND=idk)
                zeta = (0.0d0, 0.0d0)
                kfn = 0
                ifail = 0
                CALL COULCC(zrho, zeta, zl_min, 1, fc, gc, fcp, gcp, sig, mode, kfn, ifail)
                val = gc(1) * CMPLX(0.0d0, 1.0d0, KIND=idk)**(-l)              
            END IF
        END IF
        
    END FUNCTION general_asymptotic


    FUNCTION analytic_asymptotic(x, E, m, Z, l) RESULT(val)
        REAL(KIND=idk), INTENT(IN) :: x, E, m, Z
        INTEGER, INTENT(IN) :: l
        COMPLEX(KIND=idk) :: val
        
        REAL(KIND=idk) :: k, eta_real, kappa
        COMPLEX(KIND=idk) :: zrho, zeta, zl_min
        COMPLEX(KIND=idk) :: fc(1), gc(1), fcp(1), gcp(1), sig(1)
        INTEGER :: ifail
        INTEGER :: kfn, mode
        
        zl_min = CMPLX(l, 0.0d0, KIND=idk)
        
        
        IF (E >= 0.0d0) THEN
            mode = 12
            k = SQRT(2.0d0 * m * E) / hbar
            eta_real = -Z * m / (hbar**2 * k)         
            zrho = CMPLX(k * x, 0.0d0, KIND=idk)
            zeta = CMPLX(eta_real, 0.0d0, KIND=idk)
            kfn = 0
            ifail = 0        
            CALL COULCC(zrho, zeta, zl_min, 1, fc, gc, fcp, gcp, sig, mode, kfn, ifail)
            val = gc(1)
            
        ELSE
            mode = 22
            kappa = SQRT(2.0d0 * m * ABS(E)) / hbar    
            eta_real = -Z * m / (hbar**2 * kappa) 
            zrho = CMPLX(0.0d0, kappa * x, KIND=idk)
            zeta = CMPLX(0.0d0, -eta_real, KIND=idk)
            kfn = 0
            ifail = 0
            CALL COULCC(zrho, zeta, zl_min, 1, fc, gc, fcp, gcp, sig, mode, kfn, ifail)
            val = gc(1)           
        END IF
        
    END FUNCTION analytic_asymptotic


    FUNCTION general_free_solution(x, E, m, Z, l) RESULT(val)
        REAL(KIND=idk), INTENT(IN) :: x, E, m, Z
        INTEGER, INTENT(IN) :: l
        REAL(KIND=idk) :: val

        REAL(KIND=idk) :: k, eta, norm_factor
        COMPLEX(KIND=idk) :: zrho, zeta, zl_min
        COMPLEX(KIND=idk) :: fc(1), gc(1), fcp(1), gcp(1), sig(1)
        INTEGER :: ifail, kfn, mode

        IF (E < 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Energy must be positive for free solution')
            val = 0.0d0
            RETURN
        END IF

        k = SQRT(2.0d0 * m * E) / hbar

        IF (ABS(Z) > 1.0d-10) THEN
            eta = -Z * m / (hbar**2 * k) 
        ELSE
            eta = 0.0d0
        END IF

        zrho = CMPLX(k * x, 0.0d0, KIND=idk)
        zeta = CMPLX(eta, 0.0d0, KIND=idk)
        zl_min = CMPLX(l, 0.0d0, KIND=idk)
        mode = 12
        kfn = 0
        ifail = 0

        CALL COULCC(zrho, zeta, zl_min, 1, fc, gc, fcp, gcp, sig, mode, kfn, ifail)

        norm_factor = SQRT(2.0_idk * m / (pi * k * hbar**2))

        val = REAL(fc(1), KIND=idk) * norm_factor

    END FUNCTION general_free_solution






    ! Functions to compute wronskians
    FUNCTION wronskian_c(f, g, dx)
        IMPLICIT NONE
        COMPLEX(KIND=qdk), INTENT(IN) :: f(:), g(:)
        REAL(KIND=qdk), INTENT(IN) :: dx
        COMPLEX(KIND=qdk) :: wronskian_c
        INTEGER :: Nc
        COMPLEX(KIND=qdk) :: df, dg

        Nc = CEILING(mp * R0/ABS(xmax-xmin))
        df = (-f(Nc+2) + 8*f(Nc+1) - 8*f(Nc-1) + f(Nc-2)) / (12*dx)
        dg = (-g(Nc+2) + 8*g(Nc+1) - 8*g(Nc-1) + g(Nc-2)) / (12*dx)
        wronskian_c = dg*f(Nc) - df*g(Nc)
        
        IF (wronskian_c == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Wronskian equal to zero')
        END IF
        IF (ISNAN(ABS(df))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF
        IF (ISNAN(ABS(dg))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF

    END FUNCTION wronskian_c

    FUNCTION wronskian_r(f, g, dx)
        IMPLICIT NONE
        REAL(KIND=qdk), INTENT(IN) :: f(:), g(:)
        REAL(KIND=qdk), INTENT(IN) :: dx
        REAL(KIND=qdk) :: wronskian_r
        INTEGER :: Nc
        REAL(KIND=qdk) :: df, dg

        Nc = CEILING(mp * R0/ABS(xmax-xmin))
        df = (-f(Nc+2) + 8*f(Nc+1) - 8*f(Nc-1) + f(Nc-2)) / (12*dx)
        dg = (-g(Nc+2) + 8*g(Nc+1) - 8*g(Nc-1) + g(Nc-2)) / (12*dx)
        wronskian_r = dg*f(Nc) - df*g(Nc)
        
        IF (wronskian_r == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Wronskian equal to zero')
        END IF
        IF (ISNAN(ABS(df))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF
        IF (ISNAN(ABS(dg))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF
    
    END FUNCTION wronskian_r







    ! Subroutines to compute regular and irregular solutions of the homogeneous equation [E-T_N-V]|psi> = 0
    SUBROUTINE compute_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=qdk), INTENT(OUT) :: psi_R(:), psi_I(:)

        INTEGER :: N, i, j
        REAL(KIND=qdk) :: dx
        REAL(KIND=qdk), ALLOCATABLE :: k2(:)
        COMPLEX(KIND=qdk) :: wronski

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        psi_R(1) = 0.0_qdk
        psi_R(2) = dx ** (l+1)
        j = 0
        DO 
            IF (j>N-3) EXIT
            psi_I(N-j) = CMPLX(general_asymptotic(x(N-j), E, m, Z, l), KIND=qdk)
            psi_I(N-j-1) = CMPLX(general_asymptotic(x(N-j-1), E, m, Z, l), KIND=qdk)
            IF (ABS(psi_I(N-j))>0) EXIT
            j = j + 2
        END DO

        ALLOCATE(k2(N))
        DO i = 1, N
            k2(i) = 2.0_qdk * REAL(m/hbar**2, KIND=qdk) * (REAL(E,kind=qdk) - potential_green(x(i)))
        END DO

        DO i = 2, N-1
            psi_R(i+1) = (2.0_qdk*(1.0_qdk-5.0_qdk*dx**2*k2(i)/12.0_qdk)*psi_R(i) - (1.0_qdk+dx**2*k2(i-1)/12.0_qdk)*psi_R(i-1)) / (1.0_qdk+dx**2*k2(i+1)/12.0_qdk)
        END DO

        IF (ISNAN(ABS(psi_R(N)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = N-j-1, 2, -1
            psi_I(i-1) = (2.0_qdk*(1.0_qdk-5.0_qdk*dx**2*k2(i)/12.0_qdk)*psi_I(i) - (1.0_qdk+dx**2*k2(i+1)/12.0_qdk)*psi_I(i+1)) / (1.0_qdk+dx**2*k2(i-1)/12.0_qdk)
        END DO

        IF ((ABS(psi_I(1))) == 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF

        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF

        DEALLOCATE(k2)

        wronski = wronskian_c(psi_R, psi_I, dx)
        IF (wronski == 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF

        psi_I = psi_I / wronski * 2.0_qdk * REAL(m/hbar**2, KIND = qdk)

    END SUBROUTINE compute_homogeneous_solutions


    SUBROUTINE compute_analytic_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=qdk), INTENT(OUT) :: psi_R(:), psi_I(:)

        INTEGER :: N, i, j
        REAL(KIND=qdk) :: dx
        REAL(KIND=qdk), ALLOCATABLE :: k2(:)
        COMPLEX(KIND=qdk) :: wronski

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        psi_R(1) = 0.0_qdk
        psi_R(2) = dx ** (l+1)
        j = 0
        DO 
            IF (j>N-3) EXIT
            psi_I(N-j) = CMPLX(analytic_asymptotic(x(N-j), E, m, Z, l), KIND=qdk)
            psi_I(N-j-1) = CMPLX(analytic_asymptotic(x(N-j-1), E, m, Z, l), KIND=qdk)
            IF (ABS(psi_I(N-j))>0) EXIT
            j = j + 2
        END DO

        ALLOCATE(k2(N))
        DO i = 1, N
            k2(i) = 2.0_qdk * REAL(m/hbar**2, KIND=qdk) * (REAL(E,kind=qdk) - potential_green(x(i)))
        END DO

        DO i = 2, N-1
            psi_R(i+1) = (2.0_qdk*(1.0_qdk-5.0_qdk*dx**2*k2(i)/12.0_qdk)*psi_R(i) - (1.0_qdk+dx**2*k2(i-1)/12.0_qdk)*psi_R(i-1)) / (1.0_qdk+dx**2*k2(i+1)/12.0_qdk)
        END DO

        IF (ISNAN(ABS(psi_R(N)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = N-j-1, 2, -1
            psi_I(i-1) = (2.0_qdk*(1.0_qdk-5.0_qdk*dx**2*k2(i)/12.0_qdk)*psi_I(i) - (1.0_qdk+dx**2*k2(i+1)/12.0_qdk)*psi_I(i+1)) / (1.0_qdk+dx**2*k2(i-1)/12.0_qdk)
        END DO

        IF ((ABS(psi_I(1))) == 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF

        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF

        DEALLOCATE(k2)

        wronski = wronskian_c(psi_R, psi_I, dx)
        IF (wronski == 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF

        psi_I = psi_I / wronski * 2.0_qdk * REAL(m/hbar**2, KIND = qdk)

    END SUBROUTINE compute_analytic_homogeneous_solutions
    

    SUBROUTINE compute_bound_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        REAL(KIND=qdk), INTENT(OUT) :: psi_R(:), psi_I(:)

        INTEGER :: N, i, j
        REAL(KIND=qdk) :: dx
        REAL(KIND=qdk), ALLOCATABLE :: k2(:)
        REAL(KIND=qdk) :: wronski

        IF (E >= 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: Energy must be negative for bound state solutions')
            psi_R = 0.0_qdk
            psi_I = 0.0_qdk
            RETURN
        END IF

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        psi_R(1) = 0.0_qdk
        psi_R(2) = dx ** (l+1)
        j = 0
        DO 
            IF (j>N-3) EXIT
            psi_I(N-j) = REAL(general_asymptotic(x(N-j), E, m, Z, l), KIND=qdk)
            psi_I(N-j-1) = REAL(general_asymptotic(x(N-j-1), E, m, Z, l), KIND=qdk)
            IF (ABS(psi_I(N-j))>0) EXIT
            j = j + 2
        END DO

        ALLOCATE(k2(N))
        DO i = 1, N
            k2(i) = 2.0_qdk * REAL(m/hbar**2, KIND=qdk) * (REAL(E,kind=qdk) - potential_green(x(i)))
        END DO

        DO i = 2, N-1
            psi_R(i+1) = (2.0_qdk*(1.0_qdk-5.0_qdk*dx**2*k2(i)/12.0_qdk)*psi_R(i) - (1.0_qdk+dx**2*k2(i-1)/12.0_qdk)*psi_R(i-1)) / (1.0_qdk+dx**2*k2(i+1)/12.0_qdk)
        END DO

        IF (ISNAN(ABS(psi_R(N)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = N-j-1, 2, -1
            psi_I(i-1) = (2.0_qdk*(1.0_qdk-5.0_qdk*dx**2*k2(i)/12.0_qdk)*psi_I(i) - (1.0_qdk+dx**2*k2(i+1)/12.0_qdk)*psi_I(i+1)) / (1.0_qdk+dx**2*k2(i-1)/12.0_qdk)
        END DO

        IF ((ABS(psi_I(1))) == 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF

        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF

        DEALLOCATE(k2)

        wronski = wronskian_r(psi_R, psi_I, dx)
        IF (wronski == 0.0_qdk) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF

        psi_I = psi_I / wronski * 2.0_qdk * REAL(m/hbar**2, KIND = qdk)

    END SUBROUTINE compute_bound_homogeneous_solutions


    ! Subroutine to apply the Green's function method to solve the previously mentioned differential equation.
    ! |rightside> can be either a real/complex vector or a rectangular matrix (j,:), then the Green's function is applied on each j-th vector
    SUBROUTINE apply_green_real_vector(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z, rightside(:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)

        INTEGER :: N, i
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        yR = psi_R * REAL(rightside, KIND=qdk) 
        yI = psi_I * REAL(rightside, KIND=qdk) 
        CALL REVERSE(yI)
        CALL primitive(yR, dx, intR)
        CALL primitive(yI, dx, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)

        DEALLOCATE(intR,intI,psi_I,psi_R)

    END SUBROUTINE apply_green_real_vector


    SUBROUTINE apply_green_complex_vector(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)

        INTEGER :: N, i
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        yR = psi_R * CMPLX(rightside, KIND=qdk) 
        yI = psi_I * CMPLX(rightside, KIND=qdk) 
        CALL REVERSE(yI)
        CALL primitive(yR, dx, intR)
        CALL primitive(yI, dx, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)

        DEALLOCATE(intR,intI,psi_I,psi_R)

    END SUBROUTINE apply_green_complex_vector


    SUBROUTINE apply_green_real_matrix(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z, rightside(:,:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)

        INTEGER :: N, i, Nrightside
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        Nrightside = SIZE(rightside,1)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        DO i = 1, Nrightside
            yR = psi_R * REAL(rightside(i,:), KIND=qdk) 
            yI = psi_I * REAL(rightside(i,:), KIND=qdk) 
            CALL REVERSE(yI)
            CALL primitive(yR, dx, intR)
            CALL primitive(yI, dx, intI)
            CALL REVERSE(intI)
            psi(i,:) = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)
        END DO

        DEALLOCATE(intR,intI,psi_I,psi_R,yR,yI)

    END SUBROUTINE apply_green_real_matrix


    SUBROUTINE apply_green_complex_matrix(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:,:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)

        INTEGER :: N, i, Nrightside
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        Nrightside = SIZE(rightside,1)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        DO i = 1, Nrightside
            yR = psi_R * CMPLX(rightside(i,:), KIND=qdk) 
            yI = psi_I * CMPLX(rightside(i,:), KIND=qdk) 
            CALL REVERSE(yI)
            CALL primitive(yR, dx, intR)
            CALL primitive(yI, dx, intI)
            CALL REVERSE(intI)
            psi(i,:) = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)
        END DO

        DEALLOCATE(intR,intI,psi_I,psi_R,yR,yI)

    END SUBROUTINE apply_green_complex_matrix




    ! Subroutine to apply the Green's function method to solve the previously mentioned differential equation with ANALYTIC continuation of the irregular solution.
    ! |rightside> can be either a real/complex vector or a rectangular matrix (j,:), then the Green's function is applied on each j-th vector
    SUBROUTINE apply_green_analytic_real_vector(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z, rightside(:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)

        INTEGER :: N, i
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_analytic_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        yR = psi_R * REAL(rightside, KIND=qdk) 
        yI = psi_I * REAL(rightside, KIND=qdk) 
        CALL REVERSE(yI)
        CALL primitive(yR, dx, intR)
        CALL primitive(yI, dx, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)

        DEALLOCATE(intR,intI,psi_I,psi_R)

    END SUBROUTINE apply_green_analytic_real_vector


    SUBROUTINE apply_green_analytic_complex_vector(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)

        INTEGER :: N, i
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_analytic_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        yR = psi_R * CMPLX(rightside, KIND=qdk) 
        yI = psi_I * CMPLX(rightside, KIND=qdk) 
        CALL REVERSE(yI)
        CALL primitive(yR, dx, intR)
        CALL primitive(yI, dx, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)

        DEALLOCATE(intR,intI,psi_I,psi_R)

    END SUBROUTINE apply_green_analytic_complex_vector


    SUBROUTINE apply_green_analytic_real_matrix(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z, rightside(:,:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)

        INTEGER :: N, i, Nrightside
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        Nrightside = SIZE(rightside,1)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_analytic_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        DO i = 1, Nrightside
            yR = psi_R * REAL(rightside(i,:), KIND=qdk) 
            yI = psi_I * REAL(rightside(i,:), KIND=qdk) 
            CALL REVERSE(yI)
            CALL primitive(yR, dx, intR)
            CALL primitive(yI, dx, intI)
            CALL REVERSE(intI)
            psi(i,:) = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)
        END DO

        DEALLOCATE(intR,intI,psi_I,psi_R,yR,yI)

    END SUBROUTINE apply_green_analytic_real_matrix


    SUBROUTINE apply_green_analytic_complex_matrix(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:,:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)

        INTEGER :: N, i, Nrightside
        REAL(KIND=qdk) :: dx
        COMPLEX(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        N = SIZE(x)
        Nrightside = SIZE(rightside,1)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_analytic_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        DO i = 1, Nrightside
            yR = psi_R * CMPLX(rightside(i,:), KIND=qdk) 
            yI = psi_I * CMPLX(rightside(i,:), KIND=qdk) 
            CALL REVERSE(yI)
            CALL primitive(yR, dx, intR)
            CALL primitive(yI, dx, intI)
            CALL REVERSE(intI)
            psi(i,:) = CMPLX((psi_I*intR + psi_R*intI), KIND=idk)
        END DO

        DEALLOCATE(intR,intI,psi_I,psi_R,yR,yI)

    END SUBROUTINE apply_green_analytic_complex_matrix




    ! Subroutine to apply the Green's function method to solve the previously mentioned differential equation for E<0 with only REAL solutions.
    ! |rightside> can be either a real/complex vector or a rectangular matrix (j,:), then the Green's function is applied on each j-th vector
    SUBROUTINE apply_green_coulomb_bound_real_vector(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z, rightside(:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        REAL(KIND=idk), INTENT(OUT) :: psi(:)

        INTEGER :: N, i
        REAL(KIND=qdk) :: dx
        REAL(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        IF (E >= 0.0_idk) THEN
            CALL CONSOLE('[ERROR]: Energy must be negative for bound state solutions')
            psi = 0.0_idk
            RETURN
        END IF

        N = SIZE(x)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_bound_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        yR = psi_R * REAL(rightside, KIND=qdk) 
        yI = psi_I * REAL(rightside, KIND=qdk) 
        CALL REVERSE(yI)
        CALL primitive(yR, dx, intR)
        CALL primitive(yI, dx, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = REAL((psi_I*intR + psi_R*intI), KIND=idk)

        DEALLOCATE(intR,intI,psi_I,psi_R)

    END SUBROUTINE apply_green_coulomb_bound_real_vector



    SUBROUTINE apply_green_coulomb_bound_real_matrix(x, E, m, Z, l, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z, rightside(:,:)
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        REAL(KIND=idk), INTENT(OUT) :: psi(:,:)

        INTEGER :: N, i, Nrightside
        REAL(KIND=qdk) :: dx
        REAL(KIND=qdk), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI

        IF (E >= 0.0_idk) THEN
            CALL CONSOLE('[ERROR]: Energy must be negative for bound state solutions')
            psi = 0.0_idk
            RETURN
        END IF

        N = SIZE(x)
        Nrightside = SIZE(rightside,1)
        dx = REAL(ABS(x(2) - x(1)), KIND=qdk)

        ALLOCATE(psi_R(N), psi_I(N))

        call compute_bound_homogeneous_solutions(x, E, m, Z, l, potential_green, psi_R, psi_I)

        ALLOCATE(yR(N), yI(N), intR(N), intI(N))

        DO i = 1, Nrightside
            yR = psi_R * REAL(rightside(i,:), KIND=qdk) 
            yI = psi_I * REAL(rightside(i,:), KIND=qdk) 
            CALL REVERSE(yI)
            CALL primitive(yR, dx, intR)
            CALL primitive(yI, dx, intI)
            CALL REVERSE(intI)
            psi(i,:) = REAL((psi_I*intR + psi_R*intI), KIND=idk)
        END DO

        DEALLOCATE(intR,intI,psi_I,psi_R,yR,yI)

    END SUBROUTINE apply_green_coulomb_bound_real_matrix

    
END MODULE XXGreen