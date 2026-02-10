!-------------------------------< Radial Shooting >--------------------------------!
!                                                                                  !    
! Contains: Computation of bound states in Coulomb-like potential                  !                                                   
!                                                                                  !
! Last revision:    02/10/2025                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!
    
    
MODULE RadialShooting   
    USE Parameters
    USE Math
    USE whittaker_w, ONLY: coulomb_whittaker
    USE, INTRINSIC :: ieee_arithmetic
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ComputeBoundStates, rydberg_grid, compute_defect
    
    
    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_interface(R)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: R
        END FUNCTION potential_interface
    END INTERFACE
    
    
    CONTAINS
    
    SUBROUTINE compute_wavefunction(x, E, m, potential, psi)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), E, m
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
            k2(i) = 2*m/hbar**2 * (E - potential(x(i)))
        END DO
        
        psi(1) = 0.0d0
        psi(2) = h
        DO i = 2, N-1
            psi(i+1) = (2.0d0*(1.0d0-5.0d0*h**2*k2(i)/12.0d0)*psi(i) - (1.0d0+h**2*k2(i-1)/12.0d0)*psi(i-1)) / (1.0d0+h**2*k2(i+1)/12.0d0)
        END DO
        
        DEALLOCATE(k2)
    
    END SUBROUTINE compute_wavefunction
    
    
    SUBROUTINE compute_logder_diff(x, E, m, potential, Z, logderdiff)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), E, m, Z
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: potential
        REAL(KIND = idk), INTENT(OUT)                           :: logderdiff
        
        INTEGER :: N, sf
        REAL(KIND = idk) :: k, eta, w, wd
        REAL(KIND = idk) :: h
        REAL(KIND = idk), ALLOCATABLE :: psi(:)
        
        N = SIZE(x)
        k = SQRT(2*m*ABS(E))/hbar
        h = ABS(x(2)-x(1))
        eta = -Z*m/(hbar**2*k)

        CALL compute_wavefunction(x, E, m, potential, psi)
        CALL coulomb_whittaker(eta, 0, k*x(N-1), w, wd, sf)

        logderdiff = ((psi(N) - psi(N-2)) / (2.0d0*h) * w - ( k * wd) * psi(N-1)) *10.0d0 ** sf

        DEALLOCATE(psi)
    
    
    END SUBROUTINE compute_logder_diff
    
    
    
    SUBROUTINE compute_normalized_wavefunction(x, E, mass, potential, Z, psi)
        REAL(KIND = idk), INTENT(IN)                        :: x(:), E, mass, Z
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)          :: psi(:)
        
        REAL(KIND = idk), ALLOCATABLE       :: normstate(:)
        REAL(KIND = idk)                    :: dx, norm, norm_tail
        INTEGER :: i
        CHARACTER(LEN=256) :: message
        
        REAL(KIND = idk) :: k_num, eta, R_match, w_match, w_deriv
        INTEGER          :: sf_match, sf_curr, Ngrid
        REAL(KIND = idk) :: A_scale_sq, r_curr, w_curr, term, term_prev
        REAL(KIND = idk) :: R_turn
        LOGICAL          :: tail_converged
        
        REAL(KIND = idk), PARAMETER         :: eps = 1.0d-50
        
        dx = ABS(x(2)-x(1))
        Ngrid = SIZE(x)
        
        CALL compute_wavefunction(x, E, mass, potential, psi)
        
        ALLOCATE(normstate(Ngrid))
        DO i = 1, Ngrid
            normstate(i) = psi(i) ** 2
        END DO
        CALL definite_integral(normstate, dx, norm)
        DEALLOCATE(normstate)

        k_num = SQRT(2.0d0*m*ABS(E)) / hbar
        eta = - Z * m / (hbar**2 * k_num)
        R_match = x(Ngrid)
        CALL coulomb_whittaker(eta, 0, k_num*R_match, w_match, w_deriv, sf_match)
        A_scale_sq = (psi(Ngrid) / w_match) ** 2
        
        R_turn = Z / ABS(E)
        norm_tail = 0.5d0 * w_match ** 2
        r_curr = R_match
        tail_converged = .FALSE.
        
        DO WHILE (.NOT. tail_converged)
            
            r_curr = r_curr + dx
            CALL coulomb_whittaker(eta, 0, k_num*r_curr, w_curr, w_deriv, sf_curr)
            term = w_curr ** 2 * (10.0d0) ** (2*(sf_curr-sf_match))
            norm_tail = norm_tail + term
            
            IF (r_curr > 2.0d0 * R_turn .AND. term < eps ) THEN 
                tail_converged = .TRUE.
            END IF
            
            IF (r_curr > 20.0d0 * MAX(R_turn, R_match)) THEN
                WRITE(message, '(A,F0.8,A)') '[WARNING]: Normalization failed to converge for bound state at En = ', E*phys_h0, ' eV.'
                CALL CONSOLE(message)
                EXIT
            END IF
                
        END DO
        
        norm_tail = norm_tail * A_scale_sq * dx
        norm = norm + norm_tail
        norm = SQRT(ABS(norm))
    
        psi = psi / norm
        
    END SUBROUTINE compute_normalized_wavefunction
        

     
    SUBROUTINE ComputeBoundStates(x, e_test, m, potential, Z, Ntotal, Nstartgrid, eigE, eigFunc, Nfound, logderdiffs, defects)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), e_test(:), m, Z
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: potential
        INTEGER, INTENT(IN)                                     :: Ntotal, Nstartgrid
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)              :: eigE(:), eigFunc(:,:), logderdiffs(:), defects(:)
        INTEGER, INTENT(OUT)                                    :: Nfound
    
        INTEGER                                                 :: max_log_candidates
        
        REAL(KIND = idk) :: dx
        INTEGER :: i, j, N, Ne, Nadd, Nfound1
        CHARACTER(LEN=256) :: message

        INTEGER, ALLOCATABLE :: root_indices(:) 
        INTEGER :: N_potential_roots, k
        LOGICAL, ALLOCATABLE :: root_mask(:)  
        
        REAL(KIND = idk) :: Eleft, Eright, Emid, fleft, fright, fmid, Etotal
        INTEGER :: ilast
        REAL(KIND = idk), ALLOCATABLE :: psi_local(:)
        
        REAL(KIND = idk) :: mu_diff, mu_low, mu_high, mu_conv
        REAL(KIND = idk) :: E1, E2
        
        REAL(KIND = idk) :: S_x, S_y, S_xy, S_xx, y_j, x_j, x_bar, y_bar, B_fit, A_fit
        INTEGER          :: k_fit, j_fit

        
        max_log_candidates = Nstartgrid
        Nadd = Ntotal - Nstartgrid
        IF (Nadd < 0) THEN
            CALL CONSOLE('[ERROR]: Ntotal must be >= Nstart.')
            RETURN
        END IF
        
        Nfound = 0
        N = SIZE(x)
        dx = ABS(x(2)-x(1))
        Ne = SIZE(e_test)
        
        CALL CONSOLE('Searching for bound states...')
        
        IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
        IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
        IF (ALLOCATED(logderdiffs)) DEALLOCATE(logderdiffs)
        IF (ALLOCATED(defects)) DEALLOCATE(defects)
        ALLOCATE(eigE(Ntotal))
        ALLOCATE(eigFunc(N, Ntotal)) 
        ALLOCATE(logderdiffs(Ne))
        ALLOCATE(defects(Ntotal))
        
        eigE =  ieee_value(0.0d0, ieee_quiet_nan)
        eigFunc = ieee_value(0.0d0, ieee_quiet_nan) 
        defects = ieee_value(0.0d0, ieee_quiet_nan)
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(DYNAMIC, 1) 
        DO i = 1, Ne
            CALL compute_logder_diff(x, e_test(i), m, potential, Z, logderdiffs(i))
        END DO
        !$OMP END PARALLEL DO
        
        
        CALL CONSOLE('Locating potential root intervals...')
        ALLOCATE(root_mask(Ne-1))
        DO i = 1, Ne-1
            root_mask(i) = (logderdiffs(i)*logderdiffs(i+1) < 0.0d0 .OR. logderdiffs(i) == 0.0d0)
        END DO
        root_indices = PACK([(i, i=1, Ne-1)], mask=root_mask)
        DEALLOCATE(root_mask)
        N_potential_roots = SIZE(root_indices)
        
        WRITE(message, '(A,I0,A)') 'Found ', N_potential_roots, ' potential roots.'
        CALL CONSOLE(message)

        IF (N_potential_roots == 0) THEN
             DEALLOCATE(root_indices)
             WRITE(message, '(A,I0,A)') 'Computation of bound states finished. ', Nfound, ' bound states found.'
             CALL CONSOLE(message)
             RETURN
        END IF

        Nfound = MIN(N_potential_roots, Nstartgrid)


        DO k = 1, Nfound
            i = root_indices(k)
            WRITE(message, '(A,F0.4,A)') '  -> Bound state candidate found at E ~ ', (e_test(i)+e_test(i+1))*phys_h0/2.0d0, ' eV'
            CALL CONSOLE(message)
        END DO
        

        WRITE(message, '(A,I0,A)') 'Starting parallel bisection for ', Nfound, ' states...'
        CALL CONSOLE(message)
        
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP&   PRIVATE(k, i, Eleft, Eright, Emid, fleft, fright, fmid, j, psi_local) &
        !$OMP&   SCHEDULE(DYNAMIC, 1)
        DO k = 1, Nfound
            i = root_indices(k) 
            
            Eleft = e_test(i)
            Eright = e_test(i+1)
            fleft = logderdiffs(i)
            fright = logderdiffs(i+1)
            
            DO j = 1, max_iter
                Emid = 0.5d0 * (Eleft + Eright)
                CALL compute_logder_diff(x, Emid, m, potential, Z, fmid)
                
                IF (fleft * fmid <= 0.0d0) THEN
                    Eright = Emid
                    fright = fmid
                ELSE
                    Eleft = Emid
                    fleft = fmid
                END IF
            END DO
            
            eigE(k) = 0.5d0 * (Eleft + Eright)
            defects(k) = REAL(k, idk) - SQRT(0.5d0/ABS(eigE(k)))
        
            CALL compute_normalized_wavefunction(x, eigE(k), m, potential, Z, psi_local)
            
            eigFunc(:, k) = psi_local(:) 
            DEALLOCATE(psi_local)
            
        END DO
        !$OMP END PARALLEL DO
        Nfound1 = Nfound
        
        ilast = root_indices(Nfound1)
        DEALLOCATE(root_indices)
        
        WRITE(message, '(A,I0,A,F0.8,A,F0.8)') 'Last bound state with n = ', Nfound1, ' E = ', eigE(Nfound1) * phys_h0, ' found with mu = ', defects(Nfound1)
        CALL CONSOLE(message)
        
        
        DO i = Nfound1 + 1, Ntotal

            IF (i-1>N_fit_points+N_fit_min_n) THEN
                IF (ABS(defects(i-1) / defects(i-2) - 1.0d0) < mueps) THEN
                    WRITE(message, '(A,I0,A,I0,A,I0)') 'Defect converged at n = ', i-1, '. Performing QDT fit on n = ', i-N_fit_points, ' to ', i-1
                    CALL CONSOLE(message)
                    
                    S_x = 0.0d0
                    S_y = 0.0d0
                    S_xy = 0.0d0
                    S_xx = 0.0d0
                    
                    DO k_fit = 1, N_fit_points
                        j_fit = (i-N_fit_points) + k_fit - 1
                        
                        y_j = defects(j_fit)
                        x_j = 1.0d0 / (REAL(j_fit,idk)*REAL(j_fit,idk))
                        
                        S_x = S_x + x_j
                        S_y = S_y + y_j
                        S_xy = S_xy + x_j*y_j
                        S_xx = S_xx + x_j*x_j
                    END DO
                    
                    x_bar = S_x / REAL(N_fit_points,idk)
                    y_bar = S_y / REAL(N_fit_points,idk)
                    
                    B_fit = (S_xy - REAL(N_fit_points,idk)*x_bar*y_bar) / (S_xx - REAL(N_fit_points, idk) * x_bar * x_bar)
                    A_fit = y_bar - B_fit * x_bar
                    
                    WRITE(message, '(A,F0.8,A,F0.8,A)') 'Fit: mu(n) = A + B/n^2. A (mu_inf) = ', A_fit, ', B = ', B_fit
                    CALL CONSOLE(message)
                    
                
                    !$OMP PARALLEL DO DEFAULT(SHARED) &
                    !$OMP&   PRIVATE(j, E1, psi_local) &
                    !$OMP&   SCHEDULE(DYNAMIC)
                    DO j = i, Ntotal
                        mu_conv = A_fit + B_fit / (REAL(j, idk) * REAL(j, idk))
                        E1 = -0.5d0 / (REAL(j, idk) - mu_conv)**2
                        eigE(j) = E1
                        defects(j) = mu_conv
                    
                        CALL compute_normalized_wavefunction(x, E1, m, potential, Z, psi_local)
                        eigFunc(:, j) = psi_local(:)
                        DEALLOCATE(psi_local)
                    END DO
                    !$OMP END PARALLEL DO

                    Nfound = Ntotal
                    EXIT 
                END IF
            END IF

            mu_diff = ABS(defects(i-1)-defects(i-2))
            mu_low = defects(i-1) - 2.0d0 * mu_diff
            mu_high = defects(i-1) + 2.0d0 * mu_diff
        
            E1 = -0.5d0 / (REAL(i, idk) - mu_low)**2
            E2 = -0.5d0 / (REAL(i, idk) - mu_high)**2
            Eleft = MIN(E1, E2)
            Eright = MAX(E1, E2)
        
            CALL compute_logder_diff(x, Eleft, m, potential, Z, fleft)
            CALL compute_logder_diff(x, Eright, m, potential, Z, fright)
        
            IF (fleft * fright > 0.0d0) THEN
                WRITE(message, '(A,I0,A)') '[WARNING]: 200% mu_diff bracket failed for n = ', i, '. Widening to 500%.'
                CALL CONSOLE(message)
                mu_low = defects(i-1) - 5.0d0 * mu_diff
                mu_high = defects(i-1) + 5.0d0 * mu_diff
                E1 = -0.5d0 / (REAL(i, idk) - mu_low)**2; E2 = -0.5d0 / (REAL(i, idk) - mu_high)**2
                Eleft = MIN(E1, E2); Eright = MAX(E1, E2)
            
                CALL compute_logder_diff(x, Eleft, m, potential, Z, fleft)
                CALL compute_logder_diff(x, Eright, m, potential, Z, fright)

                IF (fleft * fright > 0.0d0) THEN
                    WRITE(message, '(A,I0,A)') '[WARNING]: 500% mu_diff bracket failed for n = ', i, '. Widening to 1000%.'
                    CALL CONSOLE(message)
                    mu_low = defects(i-1) - 10.0d0 * mu_diff
                    mu_high = defects(i-1) + 10.0d0 * mu_diff
                    E1 = -0.5d0 / (REAL(i, idk) - mu_low)**2; E2 = -0.5d0 / (REAL(i, idk) - mu_high)**2
                    Eleft = MIN(E1, E2); Eright = MAX(E1, E2)
                
                    CALL compute_logder_diff(x, Eleft, m, potential, Z, fleft)
                    CALL compute_logder_diff(x, Eright, m, potential, Z, fright)

                    IF (fleft * fright > 0.0d0) THEN
                        WRITE(message, '(A,I0,A)') 'Error: Cannot bracket root for n = ', i, '. Stopping.'
                        CALL CONSOLE(message)
                        Nfound = i-1
                        EXIT 
                    END IF
                END IF
            END IF
        
            DO j = 1, max_iter2
                Emid = 0.5d0 * (Eleft + Eright)
                CALL compute_logder_diff(x, Emid, m, potential, Z, fmid)
                IF (fleft * fmid <= 0.0d0) THEN
                    Eright = Emid; fright = fmid
                ELSE
                    Eleft = Emid; fleft = fmid
                END IF
            END DO
            eigE(i) = 0.5d0 * (Eleft + Eright)
        
            defects(i) = REAL(i, idk) - SQRT(0.5d0/ABS(eigE(i)))
            CALL compute_normalized_wavefunction(x, eigE(i), m, potential, Z, psi_local)
            eigFunc(:, i) = psi_local(:)
            DEALLOCATE(psi_local)
        
            Nfound = i 
            
            WRITE(message, '(A,I0,A,F0.8)') 'Bound state with n = ', Nfound, ' found with mu = ', defects(i)
            CALL CONSOLE(message)
        
        END DO 
        
        WRITE(message, '(A,I0,A)') 'Computation of bound states finished. ', Nfound, ' bound states found.'
        CALL CONSOLE(message)
    
    END SUBROUTINE ComputeBoundStates
    
    
    
    
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
    
    
    
    
    SUBROUTINE compute_defect(boundE, defect)
        REAL(KIND = idk), INTENT(IN)                :: boundE(:)
        REAL(KIND = idk), INTENT(OUT)               :: defect(:)
        
        INTEGER :: N, i
        
        N = SIZE(boundE)
        
        DO i = 1, N
            defect(i) = i - SQRT(0.5d0/ABS(boundE(i)))
        END DO
    
    END SUBROUTINE compute_defect
    
    
    
END MODULE RadialShooting    
    
    
    
