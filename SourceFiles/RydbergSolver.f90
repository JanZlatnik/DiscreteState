!------------------------------------< Rydberg >-----------------------------------!
!                                                                                  !    
! Contains: Computation of Rydberg states and PHP states in Coulomb-like potential !                                                   
!                                                                                  !
! Last revision:    14/02/2026                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!




MODULE RydbergSolver
    USE Parameters
    USE Math
    USE Green
    USE whittaker_w, ONLY: coulomb_whittaker
    USE, INTRINSIC :: ieee_arithmetic
    USE OMP_LIB
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ComputeRydbergSystem, ComputeVdn, rydberg_grid, compute_dstate

    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_interface(R)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: R
        END FUNCTION potential_interface
    END INTERFACE   


    CONTAINS

    SUBROUTINE calc_H_wavefunction(x, E, m, l, potential, psi)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), E, m
        INTEGER, INTENT(IN)                                     :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)              :: psi(:)
        
        INTEGER :: N, i
        REAL(KIND = idk) :: h
        REAL(KIND = idk), ALLOCATABLE :: k2(:)
        
        IF(ALLOCATED(psi)) DEALLOCATE(psi)
        N = SIZE(x)
        ALLOCATE(k2(N), psi(N))
        
        h = ABS(x(2)-x(1))
        
        DO i = 1, N
            k2(i) = 2.0d0*m/hbar**2 * (E - potential(x(i)))
        END DO
        
        psi(1) = 0.0d0
        psi(2) = h**(l+1)
        DO i = 2, N-1
            psi(i+1) = (2.0d0*(1.0d0-5.0d0*h**2*k2(i)/12.0d0)*psi(i) - (1.0d0+h**2*k2(i-1)/12.0d0)*psi(i-1)) / (1.0d0+h**2*k2(i+1)/12.0d0)
        END DO
        
        DEALLOCATE(k2)
    
    END SUBROUTINE calc_H_wavefunction




    SUBROUTINE calc_H_logder_diff(x, E, m, l, potential, Z, logderdiff)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), E, m, Z
        INTEGER, INTENT(IN)                                     :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: potential
        REAL(KIND = idk), INTENT(OUT)                           :: logderdiff
        
        INTEGER :: N, sf
        REAL(KIND = idk) :: k, eta, w, wd
        REAL(KIND = idk) :: h
        REAL(KIND = idk), ALLOCATABLE :: psi(:)
        
        N = SIZE(x)
        k = SQRT(2.0d0*m*ABS(E))/hbar
        h = ABS(x(2)-x(1))
        eta = -Z*m/(hbar**2*k)

        CALL calc_H_wavefunction(x, E, m, l, potential, psi)
        CALL coulomb_whittaker(eta, l, k*x(N-1), w, wd, sf)

        logderdiff = ((psi(N) - psi(N-2)) / (2.0d0*h) * w - ( k * wd) * psi(N-1)) *10.0d0 ** sf

        DEALLOCATE(psi)
    
    
    END SUBROUTINE calc_H_logder_diff




    SUBROUTINE calc_PHP_overlap(x, Ndstate, dstate, E, m, l, potential, Z, overlap)
        REAL(KIND = idk), INTENT(IN)                        :: x(:), dstate(:), E, m, Z
        INTEGER, INTENT(IN)                                 :: l
        INTEGER, INTENT(IN)                                 :: Ndstate
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), INTENT(OUT)                       :: overlap
    
        
        REAL(KIND = idk), ALLOCATABLE       :: gdstate(:)
        REAL(KIND = idk), ALLOCATABLE       :: overlapstate(:)
        REAL(KIND = idk)                    :: dx
        INTEGER :: i, Nmin
        CHARACTER(LEN=256) :: message
        
        dx = ABS(x(2)-x(1))
        
        ALLOCATE(gdstate(SIZE(x)))
        CALL apply_green_coulomb_bound(x, E, m, Z, l, potential, dstate, gdstate)
        
        Nmin = MIN(Ndstate,SIZE(x))
        ALLOCATE(overlapstate(Nmin))
        DO i = 1, Nmin
            overlapstate(i) = gdstate(i) * dstate(i)
        END DO
        DEALLOCATE(gdstate) 
        CALL definite_integral(overlapstate, dx, overlap)
        DEALLOCATE(overlapstate)
        
        overlap = TANH(overlap)

        IF (ieee_is_nan(overlap)) THEN
            !$OMP CRITICAL
            WRITE(message, '(A, G0.6)') '[WARNING]: calc_PHP_overlap resulted in NaN at E = ', E
            CALL CONSOLE(message)
            !$OMP END CRITICAL
        END IF

        
        
    END SUBROUTINE calc_PHP_overlap




    SUBROUTINE calc_PHP_normstates(x, dstate, E, mass, l, potential, Z, psi)
        REAL(KIND = idk), INTENT(IN)                        :: x(:), dstate(:), E, mass, Z
        INTEGER, INTENT(IN)                                 :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)          :: psi(:)
        
        REAL(KIND = idk), ALLOCATABLE       :: gdstate(:)
        REAL(KIND = idk), ALLOCATABLE       :: normstate(:)
        REAL(KIND = idk)                    :: dx, norm, norm_tail
        INTEGER :: i
        CHARACTER(LEN=256) :: message
        
        REAL(KIND = idk) :: k_num, eta, R_match, w_match, w_deriv
        INTEGER          :: sf_match, sf_curr, Ngrid
        REAL(KIND = idk) :: A_scale_sq, r_curr, w_curr, term, term_prev, dr_adaptive
        REAL(KIND = idk) :: R_turn, E_kin_local, V_coulomb, local_k
        LOGICAL          :: tail_converged
        
        REAL(KIND = idk), PARAMETER         :: eps = 1.0d-30
        
        dx = ABS(x(2)-x(1))
        Ngrid = SIZE(x)
        
        ALLOCATE(gdstate(Ngrid))
        CALL apply_green_coulomb_bound(x, E, mass, Z, l, potential, dstate, gdstate)
        
        ALLOCATE(normstate(Ngrid))
        DO i = 1, Ngrid
            normstate(i) = gdstate(i) ** 2
        END DO
        CALL definite_integral(normstate, dx, norm)
        DEALLOCATE(normstate)

        k_num = SQRT(2.0d0*m*ABS(E)) / hbar
        eta = - Z * m / (hbar**2 * k_num)
        R_match = x(Ngrid)
        CALL coulomb_whittaker(eta, l, k_num*R_match, w_match, w_deriv, sf_match)
        A_scale_sq = (gdstate(Ngrid) / w_match) ** 2
        
        R_turn = Z / ABS(E)
        norm_tail = 0.0d0
        r_curr = R_match
        tail_converged = .FALSE.
        term_prev = w_match ** 2
        dr_adaptive = dx
        
        DO WHILE (.NOT. tail_converged)

            V_coulomb = - Z / r_curr
            E_kin_local = E - V_coulomb

            IF (E_kin_local > 0.0d0) THEN
                local_k = SQRT(2.0d0 * mass * E_kin_local) / hbar
                IF(local_k > 1.0d-10) THEN
                    dr_adaptive = 2.0d0 * PI / local_k / 1000.0d0
                ELSE
                    dr_adaptive = dr_adaptive * 1.1d0
                END IF
            ELSE
                dr_adaptive = dr_adaptive * 1.05d0
            END IF
                
            dr_adaptive = MAX(dr_adaptive, dx)
            r_curr = r_curr + dr_adaptive

            CALL coulomb_whittaker(eta, l, k_num*r_curr, w_curr, w_deriv, sf_curr)
            term = w_curr ** 2 * (10.0d0) ** (2*(sf_curr-sf_match))
            norm_tail = norm_tail + 0.5d0 * (term + term_prev) * dr_adaptive
            term_prev = term
            
            IF (r_curr > MAX(R_turn, R_match) ) THEN 
                IF (term < eps * norm_tail) THEN
                    tail_converged = .TRUE.
                END IF
            END IF
            
            IF (r_curr > 100.0d0 * MAX(R_turn, R_match)) THEN
                WRITE(message, '(A,F0.8,A)') '[WARNING]: Normalization failed to converge for PHP state at En = ', E*phys_h0, ' eV.'
                CALL CONSOLE(message)
                EXIT
            END IF
                
        END DO
        
        norm_tail = norm_tail * A_scale_sq
        norm = norm + norm_tail
        norm = SQRT(ABS(norm))
        
        IF (ALLOCATED(psi)) DEALLOCATE(psi)
        ALLOCATE(psi(SIZE(x)))
        psi = 0.0d0
        DO i = 1, SIZE(x)
            psi(i) = gdstate(i) / norm
        END DO
        
        DEALLOCATE(gdstate) 
        
    END SUBROUTINE calc_PHP_normstates




    SUBROUTINE ComputeVdn(x, dstate, mass, potential, phin, Vdn)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), dstate(:), mass, phin(:,:)
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: potential
        REAL(KIND = idk), INTENT(OUT)                           :: Vdn(:)
        
        INTEGER :: i, j, nstates, Ncut, N
        REAL(KIND = idk), ALLOCATABLE :: ham_dstate(:), psi(:)
        REAL(KIND = idk) :: dx


        
        N = SIZE(x)
        nstates = SIZE(phin, 2)
        dx = ABS(x(2)-x(1))
        
        ALLOCATE(ham_dstate(N))
        ham_dstate = 0.0d0
        ham_dstate(1) = - (hbar**2 / (2.0d0 * mass)) * (dstate(2) - 2.0d0*dstate(1)) / (dx**2) + potential(x(1)) * dstate(1)
        ham_dstate(1) = ham_dstate(1) * dstate(1)
        ham_dstate(N) = - (hbar**2 / (2.0d0 * mass)) * (dstate(N) - 2.0d0*dstate(N) + dstate(N-1)) / (dx**2) + potential(x(N)) * dstate(N)
        ham_dstate(N) = ham_dstate(N)
        DO i = 2, N-1
            ham_dstate(i) = - (hbar**2 / (2.0d0 * mass)) * (dstate(i+1) - 2.0d0*dstate(i) + dstate(i-1)) / (dx**2) + potential(x(i)) * dstate(i)
        END DO
        
        Ncut = N
        DO i = N, 1, -1
            IF (ABS(ham_dstate(i)) > 1.0d-200) THEN
                Ncut = i
                EXIT
            END IF
        END DO
        
        !$OMP PARALLEL DO &
        !$OMP&   PRIVATE(j, i, psi) &
        !$OMP&   SHARED(nstates, Ncut, phin, ham_dstate, dx, Vdn) &
        !$OMP&   DEFAULT(NONE) &
        !$OMP&   SCHEDULE(DYNAMIC)
        DO j = 1, nstates
            ALLOCATE(psi(Ncut))
            psi = 0.0d0
            DO i = 1, Ncut
                psi(i) = phin(i,j) * ham_dstate(i)
            END DO
            CALL definite_integral(psi, dx, Vdn(j))
            DEALLOCATE(psi)
        END DO
        !$OMP END PARALLEL DO
        
        
        DEALLOCATE(ham_dstate)
    
    
    END SUBROUTINE ComputeVdn




    SUBROUTINE rydberg_grid(N, p, t, Emin, egrid)
        INTEGER, INTENT(IN)                         :: N, p
        REAL(KIND = idk), INTENT(IN)                :: Emin, t
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)  :: egrid(:)
        
        INTEGER :: jmin, jmax, j, i, M
        
        jmin = p
        jmax = CEILING(p * N * (1.0d0 + t))
        M = jmax - jmin + 1
        IF (ALLOCATED(egrid)) DEALLOCATE(egrid)
        ALLOCATE(egrid(M))
        
        i = 1
        DO j = jmin, jmax
            egrid(i) = Emin * p / j * p / j 
            i = i + 1
        END DO
    
    END SUBROUTINE rydberg_grid





    SUBROUTINE ComputeRydbergSystem(x, dstate, e_test, m, Z, l, potential, Ntotal, Nstartgrid, Nfitpoints, logdergrid, eigE_H, eigE_PHP, eigFunc_PHP, defects_H, defects_PHP)
        REAL(KIND=idk), INTENT(IN)                                          :: x(:), dstate(:), e_test(:), m, Z
        INTEGER, INTENT(IN)                                                 :: l, Ntotal, Nstartgrid, Nfitpoints
        PROCEDURE(potential_interface), POINTER, INTENT(IN)                 :: potential
        REAL(KIND=idk), ALLOCATABLE, INTENT(OUT)                            :: logdergrid(:), eigE_H(:), eigE_PHP(:), eigFunc_PHP(:,:), defects_H(:), defects_PHP(:)

        INTEGER :: i, j, k, Ne, n_H_found_grid
        LOGICAL, ALLOCATABLE :: mask_H(:)
        INTEGER, ALLOCATABLE :: idx_H(:)

        REAL(KIND=idk), ALLOCATABLE :: H_E_work(:), H_defects_work(:)
        INTEGER :: N_work_H
        REAL(KIND=idk) :: EH_L, EH_R, EH_mid, fL, fM
        REAL(KIND=idk) :: dist_left, dist_right

        REAL(KIND=idk) :: S_x, S_y, S_xy, S_xx, x_bar, y_bar, B_fit, A_fit
        REAL(KIND=idk) :: mu_pred, En_pred, mu_prev, E_prev, E_next, mu_next
        REAL(KIND=idk), ALLOCATABLE :: H_brackets(:,:)

        REAL(KIND=idk) :: width, eps, E_PHP_R, E_PHP_L, E_PHP_M, f_PHP_M
        REAL(KIND=idk) :: f_check_L, f_check_R
        REAL(KIND=idk), ALLOCATABLE :: psi_tmp(:)

        CHARACTER(256) :: message

        REAL(KIND=idk), PARAMETER :: machine_eps_factor = 1.0d9
        REAL(KIND=idk), PARAMETER :: width_eps_factor = 1.0d-6


        ! -------------------------------------------------------------------------
        ! 1. ALOCATION AND INITIALIZATION
        ! -------------------------------------------------------------------------

        IF (ALLOCATED(logdergrid)) DEALLOCATE(logdergrid)
        IF (ALLOCATED(eigE_H)) DEALLOCATE(eigE_H)
        IF (ALLOCATED(eigE_PHP)) DEALLOCATE(eigE_PHP)
        IF (ALLOCATED(eigFunc_PHP)) DEALLOCATE(eigFunc_PHP)
        IF (ALLOCATED(defects_H)) DEALLOCATE(defects_H)
        IF (ALLOCATED(defects_PHP)) DEALLOCATE(defects_PHP)

        N_work_H = Ntotal + 1
        ALLOCATE(H_E_work(N_work_H), H_defects_work(N_work_H))
        H_E_work = 0.0d0
        H_defects_work = 0.0d0

        ALLOCATE(eigE_H(Ntotal), eigE_PHP(Ntotal), eigFunc_PHP(SIZE(x), Ntotal), defects_H(Ntotal), defects_PHP(Ntotal))
        eigE_PHP = ieee_value(0.0d0, ieee_quiet_nan)
        defects_PHP = ieee_value(0.0d0, ieee_quiet_nan)
        eigFunc_PHP = 0.0d0
        
        CALL CONSOLE('=================================================')
        CALL CONSOLE('      RYDBERG SOLVER (Unified H & PHP)           ')
        CALL CONSOLE('=================================================')

        ! -------------------------------------------------------------------------
        ! 2. GRID SEARCH FOR H-STATE EIGENVALUES
        ! -------------------------------------------------------------------------
        CALL CONSOLE('Scanning Grid for H-states...')

        Ne = SIZE(e_test)
        ALLOCATE(logdergrid(Ne))

        !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
        DO i = 1, Ne
            CALL calc_H_logder_diff(x, e_test(i), m, l, potential, Z, logdergrid(i))
        END DO
        !$OMP END PARALLEL DO

        ALLOCATE(mask_H(Ne-1))
        DO i = 1, Ne-1
            mask_H(i) = (logdergrid(i) * logdergrid(i+1) <= 0.0d0)
        END DO

        idx_H = PACK([(i, i=1, Ne-1)], mask_H)
        DEALLOCATE(mask_H)

        n_H_found_grid = MIN(SIZE(idx_H), Nstartgrid)
        WRITE(message, '(A,I0,A)') 'Found ', SIZE(idx_H), ' potential roots.'
        CALL CONSOLE(message)

        DO k = 1, n_H_found_grid
            i = idx_H(k)
            WRITE(message, '(A,G0.6,A)') '  -> Bound H state candidate found at E ~ ', (e_test(i)+e_test(i+1))*phys_h0/2.0d0, ' eV'
            CALL CONSOLE(message)
        END DO

        WRITE(message, '(A,I0,A)') 'Starting parallel bisection for ', n_H_found_grid, ' states...'
        CALL CONSOLE(message)
        

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, i, EH_L, EH_R, fL, fM, j, EH_mid)
        DO k = 1, n_H_found_grid
            i = idx_H(k)
            EH_L = e_test(i)
            EH_R = e_test(i+1)

            CALL calc_H_logder_diff(x, EH_L, m, l, potential, Z, fL)
            DO j = 1, max_iter
                EH_mid = 0.5d0 * (EH_L + EH_R)
                CALL calc_H_logder_diff(x, EH_mid, m, l, potential, Z, fM)
                IF (fL * fM <= 0.0d0) THEN
                    EH_R = EH_mid
                ELSE
                    EH_L = EH_mid
                    fL = fM
                END IF
            END DO
            H_E_work(k) = 0.5d0 * (EH_L + EH_R)
            H_defects_work(k) = REAL(k, idk) - SQRT(0.5d0/ABS(H_E_work(k)))
        END DO
        !$OMP END PARALLEL DO
        DEALLOCATE(idx_H)


        ! -------------------------------------------------------------------------
        ! 3. QDT FIT & BRACKETING FOR H-STATES
        ! -------------------------------------------------------------------------
        WRITE(message, '(A,I0,A,I0)')  'Performing QDT fit on n = ', n_H_found_grid-Nfitpoints+1, ' to ', n_H_found_grid
        CALL CONSOLE(message)

        S_x = 0.0d0; S_y = 0.0d0; S_xy = 0.0d0; S_xx = 0.0d0
        DO k = 1, N_fit_points
            j = (n_H_found_grid - N_fit_points) + k
            S_x = S_x + 1.0d0/(REAL(j, idk)**2)
            S_y = S_y + H_defects_work(j)
            S_xy = S_xy + (1.0d0/(REAL(j, idk)**2)) * H_defects_work(j)
            S_xx = S_xx + (1.0d0/(REAL(j, idk)**2))**2
        END DO

        x_bar = S_x / REAL(N_fit_points, idk)
        y_bar = S_y / REAL(N_fit_points, idk)
        B_fit = (S_xy - REAL(N_fit_points,idk)*x_bar*y_bar) / (S_xx - REAL(N_fit_points,idk)*x_bar*x_bar)
        A_fit = y_bar - B_fit * x_bar

        WRITE(message, '(A,F0.8,A,F0.8,A)') 'Fit: mu(n) = A + B/n^2. A (mu_inf) = ', A_fit, ', B = ', B_fit
        CALL CONSOLE(message)

        ALLOCATE(H_brackets(2, N_work_H))

        DO k = n_H_found_grid + 1, N_work_H
            mu_pred = A_fit + B_fit / REAL(k, idk)**2
            En_pred = - 0.5d0 / (REAL(k, idk)-mu_pred)**2

            IF (k == n_H_found_grid + 1) THEN
                E_prev = H_E_work(n_H_found_grid)
            ELSE 
                mu_prev = A_fit + B_fit / REAL(k-1, idk)**2
                E_prev = - 0.5d0 / (REAL(k-1, idk)-mu_prev)**2
            END IF

            mu_next = A_fit + B_fit / REAL(k+1, idk)**2
            E_next = - 0.5d0 / (REAL(k+1, idk)-mu_next)**2

            dist_left = ABS(En_pred - E_prev)
            dist_right = ABS(En_pred - E_next)

            H_brackets(1,k) = En_pred - dist_left / 3.0d0
            H_brackets(2,k) = En_pred + dist_right / 3.0d0
        END DO


        ! -------------------------------------------------------------------------
        ! 4. PARALLEL BISECTION FOR H-STATES IN BRACKETS
        ! -------------------------------------------------------------------------
        CALL CONSOLE('Starting parallel bisection for H states in brackets...')

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, EH_L, EH_R, fL, fM, j, EH_mid) SCHEDULE(DYNAMIC)
        DO k = n_H_found_grid + 1, N_work_H
            EH_L = H_brackets(1,k)
            EH_R = H_brackets(2,k)

            CALL calc_H_logder_diff(x, EH_L, m, l, potential, Z, fL)
            DO j = 1, max_iter2
                EH_mid = 0.5d0 * (EH_L + EH_R)
                CALL calc_H_logder_diff(x, EH_mid, m, l, potential, Z, fM)
                IF (fL * fM <= 0.0d0) THEN
                    EH_R = EH_mid
                ELSE
                    EH_L = EH_mid
                    fL = fM
                END IF
            END DO
            H_E_work(k) = 0.5d0 * (EH_L + EH_R)
            H_defects_work(k) = REAL(k, idk) - SQRT(0.5d0/ABS(H_E_work(k)))
        END DO
        !$OMP END PARALLEL DO
        DEALLOCATE(H_brackets)

        DO k = 1, N_work_H
            WRITE(message, '(A,G0.6,A)') '  -> Bound H state found at E ~ ', H_E_work(k)*phys_h0, ' eV'
            CALL CONSOLE(message)
        END DO


        ! -------------------------------------------------------------------------
        ! 5. FINDING PHP STATES
        ! -------------------------------------------------------------------------
        CALL CONSOLE('Searching for PHP states...')

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, i, EH_L, EH_R, width, eps, E_PHP_L, E_PHP_R, E_PHP_M, f_check_L, f_check_R, f_PHP_M, j, psi_tmp, message) SCHEDULE(DYNAMIC, 1)  
        DO k = 1, Ntotal
            EH_L = H_E_work(k)
            EH_R = H_E_work(k+1)
            width = ABS(EH_R - EH_L)

            eps = MIN(machine_eps_factor * EPSILON(EH_L) * ABS(EH_L), width_eps_factor * width)
            DO i = 1,200
                E_PHP_L = EH_L + eps
                CALL calc_PHP_overlap(x, SIZE(x), dstate, E_PHP_L, m, l, potential, Z, f_check_L)

                IF (f_check_L > 0.0d0) THEN
                    EXIT
                ELSE
                    eps = eps * 0.7d0
                END IF

                IF (eps < 1.0d0 * EPSILON(EH_L) * ABS(EH_L)) THEN
                    f_check_L = -1.0d0 
                    EXIT
                END IF
            END DO

            IF (f_check_L <= 0.0d0) THEN
                !$OMP CRITICAL
                WRITE(message, '(A,I0,A,I0,A)') '[INFO]: PHP state n=', k, ' coincides with H-pole n=', k, '.'
                CALL CONSOLE(message)
                !$OMP END CRITICAL

                eigE_PHP(k) = EH_L
                defects_PHP(k) = H_defects_work(k)
                CALL calc_PHP_normstates(x, dstate, eigE_PHP(k), m, l, potential, Z, psi_tmp)
                eigFunc_PHP(:,k) = psi_tmp(:)
                DEALLOCATE(psi_tmp)
                CYCLE 
            END IF

            eps = MIN(machine_eps_factor * EPSILON(EH_R) * ABS(EH_R), width_eps_factor * width)
            DO i = 1,200
                E_PHP_R = EH_R - eps
                CALL calc_PHP_overlap(x, SIZE(x), dstate, E_PHP_R, m, l, potential, Z, f_check_R)

                IF (f_check_R < 0.0d0) THEN
                    EXIT
                ELSE
                    eps = eps * 0.7d0
                END IF

                IF (eps < 1.0d0 * EPSILON(EH_R) * ABS(EH_R)) THEN
                    f_check_R = 1.0d0 
                    EXIT
                END IF
            END DO

            IF (f_check_R >= 0.0d0) THEN
                !$OMP CRITICAL
                WRITE(message, '(A,I0,A,I0,A)') '[INFO]: PHP state n=', k, ' coincides with H-pole n=', k+1, '.'
                CALL CONSOLE(message)
                !$OMP END CRITICAL
                
                eigE_PHP(k) = EH_R
                defects_PHP(k) = H_defects_work(k+1) - 1.0d0
                CALL calc_PHP_normstates(x, dstate, eigE_PHP(k), m, l, potential, Z, psi_tmp)
                eigFunc_PHP(:,k) = psi_tmp(:)
                DEALLOCATE(psi_tmp)
                CYCLE
            END IF


            DO j = 1, max_iter
                E_PHP_M = 0.5d0 * (E_PHP_L + E_PHP_R)
                CALL calc_PHP_overlap(x, SIZE(x), dstate, E_PHP_M, m, l, potential, Z, f_PHP_M)

                IF (f_PHP_M < 0.0d0) THEN
                    E_PHP_R = E_PHP_M
                ELSE
                    E_PHP_L = E_PHP_M
                END IF
            END DO
            eigE_PHP(k) = 0.5d0 * (E_PHP_L + E_PHP_R)
            defects_PHP(k) = REAL(k, idk) - SQRT(0.5d0/ABS(eigE_PHP(k)))

            CALL calc_PHP_normstates(x, dstate, eigE_PHP(k), m, l, potential, Z, psi_tmp)
            eigFunc_PHP(:,k) = psi_tmp(:)
            DEALLOCATE(psi_tmp)

        END DO
        !$OMP END PARALLEL DO


        ! -------------------------------------------------------------------------
        ! 6. EXPORTING RESULTS
        ! -------------------------------------------------------------------------
        eigE_H = H_E_work(1:Ntotal)
        defects_H = H_defects_work(1:Ntotal)

        DEALLOCATE(H_E_work, H_defects_work)

        CALL CONSOLE('Rydberg System Computation Successful.')




    END SUBROUTINE ComputeRydbergSystem





    SUBROUTINE compute_dstate(x, m, l, Z, V_asymptotic, Emin_guess, dstate, e_asympt, info_dstate)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), m, Z, Emin_guess
        INTEGER, INTENT(IN)                                     :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: V_asymptotic
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)              :: dstate(:)
        REAL(KIND = idk), INTENT(OUT)                           :: e_asympt
        INTEGER, INTENT(OUT)                                    :: info_dstate
        
        REAL(KIND = idk) :: E_curr, E_next, f_curr, f_next, dE
        REAL(KIND = idk) :: E_L, E_R, E_mid, f_L, f_mid, norm, norm_tail, dx
        INTEGER          :: i, N
        REAL(KIND = idk), ALLOCATABLE :: psi_tmp(:), normstate(:)
        CHARACTER(LEN=256) :: message
        
        REAL(KIND = idk) :: k_num, eta, R_match, w_match, w_deriv
        INTEGER          :: sf_match, sf_curr
        REAL(KIND = idk) :: A_scale_sq, r_curr, w_curr, term, term_prev, dr_adaptive
        REAL(KIND = idk) :: R_turn, E_kin_local, V_coul, local_k
        LOGICAL          :: tail_converged
        REAL(KIND = idk), PARAMETER :: eps = 1.0d-30

        N = SIZE(x)
        dx = ABS(x(2) - x(1))
        dE = ABS(Emin_guess) * 5.0d-3 / m 
        E_curr = Emin_guess
        
        CALL calc_H_logder_diff(x, E_curr, m, l, V_asymptotic, Z, f_curr)

        DO
            E_next = E_curr + dE
            IF (E_next >= 0.0d0) THEN
                CALL CONSOLE('[ERROR]: No bound state found for in dstate initialization!')
                info_dstate = -1
                e_asympt = ieee_value(0.0d0, ieee_quiet_nan)
                RETURN
            END IF
            
            CALL calc_H_logder_diff(x, E_next, m, l, V_asymptotic, Z, f_next)
            IF (f_curr * f_next <= 0.0d0) THEN
                E_L = E_curr; E_R = E_next; f_L = f_curr
                EXIT
            END IF
            E_curr = E_next
            f_curr = f_next
        END DO

        DO i = 1, max_iter * 2
            E_mid = 0.5d0 * (E_L + E_R)
            CALL calc_H_logder_diff(x, E_mid, m, l, V_asymptotic, Z, f_mid)
            IF (f_L * f_mid <= 0.0d0) THEN
                E_R = E_mid
            ELSE
                E_L = E_mid; f_L = f_mid
            END IF
        END DO
        e_asympt = 0.5d0 * (E_L + E_R)

        CALL calc_H_wavefunction(x, e_asympt, m, l, V_asymptotic, psi_tmp)
        
        ALLOCATE(normstate(N))
        DO i = 1, N
            normstate(i) = psi_tmp(i)**2
        END DO
        CALL definite_integral(normstate, dx, norm)
        DEALLOCATE(normstate)
        
        k_num = SQRT(2.0d0 * m * ABS(e_asympt)) / hbar
        eta = - Z * m / (hbar**2 * k_num)
        R_match = x(N)
        
        CALL coulomb_whittaker(eta, l, k_num * R_match, w_match, w_deriv, sf_match)
        A_scale_sq = (psi_tmp(N) / w_match) ** 2
        
        R_turn = Z / ABS(e_asympt)
        norm_tail = 0.0d0
        r_curr = R_match
        tail_converged = .FALSE.
        term_prev = w_match ** 2
        dr_adaptive = dx
        
        DO WHILE (.NOT. tail_converged)
            V_coul = - Z / r_curr
            E_kin_local = e_asympt - V_coul

            IF (E_kin_local > 0.0d0) THEN
                local_k = SQRT(2.0d0 * m * E_kin_local) / hbar
                IF(local_k > 1.0d-10) THEN
                    dr_adaptive = 2.0d0 * PI / local_k / 1000.0d0
                ELSE
                    dr_adaptive = dr_adaptive * 1.1d0
                END IF
            ELSE
                dr_adaptive = dr_adaptive * 1.05d0
            END IF
                
            dr_adaptive = MAX(dr_adaptive, dx)
            r_curr = r_curr + dr_adaptive

            CALL coulomb_whittaker(eta, l, k_num * r_curr, w_curr, w_deriv, sf_curr)
            term = w_curr ** 2 * (10.0d0) ** (2*(sf_curr - sf_match))
            norm_tail = norm_tail + 0.5d0 * (term + term_prev) * dr_adaptive
            term_prev = term
            
            IF (r_curr > MAX(R_turn, R_match) ) THEN 
                IF (term < eps * norm_tail) tail_converged = .TRUE.
            END IF
            
            IF (r_curr > 100.0d0 * MAX(R_turn, R_match)) THEN
                WRITE(message, '(A)') '[ERROR]: Normalization of dstate tail did not converge!'
                CALL CONSOLE(message)
                info_dstate = -1
                RETURN
            END IF
        END DO
        
        norm_tail = norm_tail * A_scale_sq
        norm = SQRT(ABS(norm + norm_tail))
        
        IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
        ALLOCATE(dstate(N))
        DO i = 1, N
            dstate(i) = psi_tmp(i) / norm
        END DO
        DEALLOCATE(psi_tmp)

        WRITE(message, '(A,G0.8,A)') ' -> dstate found at energy E = ', e_asympt*phys_h0, ' eV'
        info_dstate = 0
        CALL CONSOLE(message)

    END SUBROUTINE compute_dstate
        


END MODULE RydbergSolver