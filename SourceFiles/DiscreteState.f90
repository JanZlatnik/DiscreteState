!--------------------------------< Discrete State >--------------------------------!
!                                                                                  !    
! Contains: Computation of Delta and Gamma for a given potential and DS            !                                                   
!                                                                                  !
! Last revision:    25/09/2025                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!
    

MODULE DiscreteState
    USE Math
    USE Green
    USE Parameters
    USE COULCC_M, only: COULCC
    USE whittaker_w, only: coulomb_whittaker
    USE OMP_LIB
    USE, INTRINSIC :: ieee_arithmetic
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: compute_DS, unwrap_phaseshift, compute_DSenergy

    INTEGER, PARAMETER :: xpick = 20, xn = 10

    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_interface(R)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: R
        END FUNCTION potential_interface
    END INTERFACE
    
    
    CONTAINS
    
    
    
    ! Subroutine to compute discrete state properties (Delta, Gamma, phase shift) for a given potential and discrete state wavefunction
    SUBROUTINE compute_DS(x, E, mass, Z, l, potential, dstate, Delta, Gamma2, DeltaA, Gamma2A, phaseshift, DS_phaseshift, Vde)
        USE OMP_LIB
        REAL(KIND = idk), INTENT(IN) :: x(:), E(:), mass, dstate(:), Z
        INTEGER, INTENT(IN) :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: Delta(:), Gamma2(:), DeltaA(:), Gamma2A(:), phaseshift(:), DS_phaseshift(:), Vde(:)
        
        INTEGER :: i, j, N
        COMPLEX(KIND = idk), ALLOCATABLE :: green_dstate(:), scattered_state(:), greened_mat(:,:)
        REAL(KIND = idk), ALLOCATABLE :: ham_dstate(:), freestate(:), mat(:,:)
        REAL(KIND = idk) :: k, dx, hamdot
        COMPLEX(KIND=idk) :: greendot, scattdot
        
        dx = ABS(x(2) - x(1)) 
        N = SIZE(x)
  
        ALLOCATE(Delta(SIZE(E)))
        ALLOCATE(DeltaA(SIZE(E)))
        ALLOCATE(Gamma2(SIZE(E)))
        ALLOCATE(Gamma2A(SIZE(E)))
        ALLOCATE(phaseshift(SIZE(E)))
        ALLOCATE(DS_phaseshift(SIZE(E)))
        ALLOCATE(Vde(SIZE(E)))
        Delta = 0.0d0
        DeltaA = 0.0d0
        Gamma2 = 0.0d0
        Gamma2A = 0.0d0
        phaseshift = 0.0d0
        DS_phaseshift = 0.0d0
        Vde = 0.0d0
        
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(j, i, hamdot, greendot, scattdot, &
        !$OMP     ham_dstate, green_dstate, mat, greened_mat, freestate, scattered_state) &
        !$OMP SHARED(x, dstate, E, mass, Z, dx, Delta, Gamma2, DeltaA, Gamma2A, phaseshift, DS_phaseshift, Vde, N) &
        !$OMP SCHEDULE(DYNAMIC)
        DO j = 1, SIZE(E)

            greendot = (0.0d0, 0.0d0)
            scattdot = (0.0d0, 0.0d0)
            
            call compute_DSenergy(x, mass, potential, dstate, hamdot)
            
            
            if (e(j) >= 0.0d0) then
                allocate(mat(2,N))
                allocate(freestate(N))
                do i = 1, N
                    freestate(i) = general_free_solution(x(i), e(j), mass, Z, l)
                    mat(2,i) = freestate(i) * (potential(x(i)) + Z /(x(i)) - REAL(l*(l+1),KIND=idk)/(2.0d0*mass*x(i)**2)* hbar**2 )
                end do
                allocate(greened_mat(2,N))
                mat(1,:) = dstate(:)
                call apply_green(x, e(j), mass, z, l, potential, mat, greened_mat)
                deallocate(mat)

                allocate(green_dstate(N))
                do i = 1, N
                    green_dstate(i) = greened_mat(1,i) * dstate(i)
                end do

                allocate(scattered_state(N))
                scattered_state = freestate + greened_mat(2,:)
                deallocate(freestate)

                deallocate(greened_mat)
                call definite_integral(green_dstate, dx, greendot)
                greendot = 1.0d0 / greendot
                deallocate(green_dstate)

                delta(j) = e(j) - hamdot - real(greendot)
                gamma2(j) = 2.0d0 * aimag(greendot)

                call compute_phaseshift(scattered_state, x, e(j), Z, mass, l, phaseshift(j))
                ds_phaseshift(j) = - atan2(aimag(greendot),real(greendot))

                do i = 1, N
                    scattered_state(i) = scattered_state(i) * dstate(i)
                end do
                call definite_integral(scattered_state, dx, scattdot)
                deallocate(scattered_state)
                vde(j) = abs(scattdot * greendot)

            else

                allocate(green_dstate(N))
                call apply_green(x, e(j), mass, Z, l, potential, dstate, green_dstate)
                do i = 1, N
                    green_dstate(i) = green_dstate(i) * dstate(i)
                end do
                call definite_integral(green_dstate, dx, greendot)
                greendot = 1.0d0 / greendot
                delta(j) = e(j) - hamdot - real(greendot)
                gamma2(j) = 2.0d0 * aimag(greendot)
                vde(j) = sqrt(abs(gamma2(j)) / (2.0d0 * pi))

                if (abs(Z) > 1.0d-10 .AND. e(j)>cutoff_energy(l)) then
                    call compute_analytic_defect(x, e(j), mass, Z, l, potential, phaseshift(j))
                    call apply_green_analytic(x, e(j), mass, Z, l, potential, dstate, green_dstate)
                    do i = 1, N
                        green_dstate(i) = green_dstate(i) * dstate(i)
                    end do
                    call definite_integral(green_dstate, dx, greendot)
                    deallocate(green_dstate)
                    greendot = 1.0d0 / greendot
                    ds_phaseshift(j) = atan2(aimag(greendot),real(greendot))
                    deltaA(j) = e(j) - hamdot - real(greendot)
                    gamma2A(j) = -2.0d0 * aimag(greendot)

                else
                    deallocate(green_dstate)
                    phaseshift(j) = ieee_value(0.0d0, ieee_quiet_nan)
                    ds_phaseshift(j) = ieee_value(0.0d0, ieee_quiet_nan)
                    deltaA(j) = ieee_value(0.0d0, ieee_quiet_nan)
                    gamma2A(j) = ieee_value(0.0d0, ieee_quiet_nan)
                end if
                

            end if
            
        END DO
        !$OMP END PARALLEL DO
    
    END SUBROUTINE compute_DS
    
    
    
    
    
    ! Subroutine to extract phase shift from scattered state in a general potential
    SUBROUTINE compute_phaseshift(psi, x, e, Z, m, l, delta)
      IMPLICIT NONE
      COMPLEX(KIND = idk), INTENT(IN)  :: psi(:)
      REAL(KIND = idk),    INTENT(IN)  :: x(:)
      REAL(KIND = idk),    INTENT(IN)  :: e
      REAL(KIND = idk),    INTENT(IN)  :: Z
      REAL(KIND = idk),    INTENT(IN)  :: m
      INTEGER,    INTENT(IN)           :: l
      REAL(KIND = idk),    INTENT(OUT) :: delta

      COMPLEX(KIND = idk) :: M11, M12, M21, M22, R1, R2
      COMPLEX(KIND = idk) :: A_in, A_out, b_minus, b_plus, S
      INTEGER :: N, i0, i1, i, sf

      N  = SIZE(psi)
      i0 = MAX(1, N - xpick - xn)
      i1 = MIN(N, N - xpick + xn)

      M11 = (0.0d0,0.0d0); M12 = (0.0d0,0.0d0)
      M21 = (0.0d0,0.0d0); M22 = (0.0d0,0.0d0)
      R1  = (0.0d0,0.0d0); R2  = (0.0d0,0.0d0)

      DO i = i0, i1
         b_plus  = general_asymptotic(x(i), e, m, Z, l, sf)
         b_minus = CONJG(b_plus)

         M11 = M11 + CONJG(b_minus)*b_minus
         M12 = M12 + CONJG(b_minus)*b_plus
         M21 = M21 + CONJG(b_plus)*b_minus
         M22 = M22 + CONJG(b_plus)*b_plus
         R1  = R1  + CONJG(b_minus)*psi(i)
         R2  = R2  + CONJG(b_plus) *psi(i)
      END DO

      A_in  = ( M22*R1 - M12*R2 ) 
      A_out = (-M21*R1 + M11*R2 ) 
      S     = -A_out / A_in 

      delta = 0.5d0 * ATAN2(AIMAG(S), REAL(S))

    END SUBROUTINE compute_phaseshift



    ! Subroutine to compute analytical continuation of quantum defect
    SUBROUTINE compute_analytic_defect(x, e, Z, m, l, potential, defect_phase)
        USE OMP_LIB
        IMPLICIT NONE
        
        REAL(KIND = idk), INTENT(IN)  :: x(:)
        REAL(KIND = idk), INTENT(IN)  :: e, m, Z
        INTEGER, INTENT(IN)           :: l
        REAL(KIND = idk), INTENT(OUT) :: defect_phase
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential

        INTEGER :: N, i, i0, i1, i_SR
        REAL(KIND = idk) :: dx, kappa, nu
        REAL(KIND = idk), ALLOCATABLE :: psi(:), k2(:)
        COMPLEX(KIND=idk) :: H_plus, H_minus
        COMPLEX(KIND=idk) :: M11, M12, M21, M22, R1, R2
        COMPLEX(KIND=idk) :: A_coef, B_coef, S_ratio
        
        IF (e >= 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Energy must be negative for analytic defect computation')
            defect_phase = ieee_value(0.0d0, ieee_quiet_nan)
            RETURN
        END IF

        N = SIZE(x)
        dx = ABS(x(2) - x(1))
        i_SR = NINT((xmax_SR-x(1))/dx) + 1 + xpick
        i_SR = MIN(MAX(i_SR,1),N)
        i0 = MAX(1, i_SR - xpick - xn)
        i1 = MIN(N, i_SR - xpick + xn)
        
        ALLOCATE(psi(i1))
        ALLOCATE(k2(i1))
        DO i = 1, i1
            k2(i) = 2.0d0 * m / (hbar**2) * (e - potential(x(i)))
        END DO
        
        psi(1) = 0.0d0
        psi(2) = dx**(l+1) 
        DO i = 2, i1-1
            psi(i+1) = (2.0d0*(1.0d0-5.0d0*dx**2*k2(i)/12.0d0)*psi(i) - (1.0d0+dx**2*k2(i-1)/12.0d0)*psi(i-1)) / (1.0d0+dx**2*k2(i+1)/12.0d0)
        END DO
        DEALLOCATE(k2)
        
        M11 = (0.0d0, 0.0d0); M12 = (0.0d0, 0.0d0)
        M21 = (0.0d0, 0.0d0); M22 = (0.0d0, 0.0d0)
        R1  = (0.0d0, 0.0d0); R2  = (0.0d0, 0.0d0)
        
        DO i = i0, i1
            H_plus  = analytic_asymptotic(x(i), e, m, Z, l, 1, .FALSE.)
            H_minus = analytic_asymptotic(x(i), e, m, Z, l, -1, .FALSE.)
            
            M11 = M11 + CONJG(H_plus)  * H_plus
            M12 = M12 + CONJG(H_plus)  * H_minus
            M21 = M21 + CONJG(H_minus) * H_plus  
            M22 = M22 + CONJG(H_minus) * H_minus
            
            R1  = R1  + CONJG(H_plus)  * CMPLX(psi(i), 0.0d0, KIND=idk)
            R2  = R2  + CONJG(H_minus) * CMPLX(psi(i), 0.0d0, KIND=idk)
        END DO
        DEALLOCATE(psi)


        A_coef = (R1 * M22 - R2 * M12)
        B_coef = (R2 * M11 - R1 * M21) 

        nu = Z * m / (hbar * SQRT(2.0d0 * m * ABS(e)) )
        
        S_ratio = - A_coef / B_coef * EXP(CMPLX(0.0d0, pi * REAL(FLOOR(nu),KIND=idk), KIND=idk))
        
        defect_phase = ATAN2(AIMAG(S_ratio), REAL(S_ratio)) 
        
    END SUBROUTINE compute_analytic_defect

    
    
    
    ! Subroutine to unwrap (make continuous) phase shifts over energy grid
    SUBROUTINE unwrap_phaseshift(phaseshift)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(INOUT) :: phaseshift(:)

        INTEGER :: i, N
        REAL(KIND = idk) :: raw_diff, k_pi, current_offset
        
        REAL(KIND = idk), PARAMETER :: tolerance = 0.1d0 

        N = SIZE(phaseshift)
        current_offset = 0.0_idk

        DO i = 2, N
            phaseshift(i) = phaseshift(i) + current_offset
        
            raw_diff = phaseshift(i) - phaseshift(i-1)

            k_pi = ANINT(raw_diff / pi)
            
            IF (k_pi /= 0.0_idk .AND. ABS((raw_diff / pi) - k_pi) < tolerance) THEN
                current_offset = current_offset - (k_pi * pi)
                phaseshift(i) = phaseshift(i) - (k_pi * pi)
            END IF
        END DO

    END SUBROUTINE unwrap_phaseshift



    
    ! Subroutine to compute discrete state energy from wavefunction and potential
    SUBROUTINE compute_DSenergy(x, mass, potential, dstate, DSenergy)
        REAL(KIND = idk), INTENT(IN)                        :: x(:), mass, dstate(:)
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), INTENT(OUT)                       :: DSenergy
        
        
        INTEGER :: i, N
        REAL(KIND = idk) :: prefactor
        REAL(KIND = idk), ALLOCATABLE :: ham_dstate(:), dx
        
        N = SIZE(x)
        dx = ABS(x(2) - x(1))
        prefactor = - (hbar**2 / (2.0d0 * mass)) / (12.0d0 * dx**2)
        
        ALLOCATE(ham_dstate(N))

        ham_dstate(1) = prefactor * (35.0d0*dstate(1) - 104.0d0*dstate(2) + 114.0d0*dstate(3) - 56.0d0*dstate(4) + 11.0d0*dstate(5)) + potential(x(1)) * dstate(1)
        ham_dstate(1) = ham_dstate(1) * dstate(1)

        ham_dstate(2) = prefactor * (11.0d0*dstate(1) - 20.0d0*dstate(2) + 6.0d0*dstate(3) + 4.0d0*dstate(4) - dstate(5)) + potential(x(2)) * dstate(2)
        ham_dstate(2) = ham_dstate(2) * dstate(2)

        DO i = 3, N-2
            ham_dstate(i) = prefactor * (-dstate(i+2) + 16.0d0*dstate(i+1) - 30.0d0*dstate(i) + 16.0d0*dstate(i-1) - dstate(i-2)) + potential(x(i)) * dstate(i)
            ham_dstate(i) = ham_dstate(i) * dstate(i)
        END DO

        ham_dstate(N-1) = prefactor * (11.0d0*dstate(N) - 20.0d0*dstate(N-1) + 6.0d0*dstate(N-2) + 4.0d0*dstate(N-3) - dstate(N-4)) + potential(x(N-1)) * dstate(N-1)
        ham_dstate(N-1) = ham_dstate(N-1) * dstate(N-1)

        ham_dstate(N) = prefactor * (35.0d0*dstate(N) - 104.0d0*dstate(N-1) + 114.0d0*dstate(N-2) - 56.0d0*dstate(N-3) + 11.0d0*dstate(N-4)) + potential(x(N)) * dstate(N)
        ham_dstate(N) = ham_dstate(N) * dstate(N)
        
        CALL definite_integral(ham_dstate, dx, DSenergy)
        DEALLOCATE(ham_dstate)  
    
    END SUBROUTINE compute_DSenergy




END MODULE DiscreteState