!----------------------------------< PHP Bound >-----------------------------------!
!                                                                                  !    
! Contains: Computation of bound states of PHP in Coulomb-like potential           !                                                   
!                                                                                  !
! Last revision:    19/10/2025                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!
    
    
MODULE PHPBound
    USE Parameters
    USE Math
    USE Green
    USE OMP_LIB
    USE whittaker_w, only: coulomb_whittaker
    USE, INTRINSIC :: ieee_arithmetic
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ComputePHPBoundStates, ComputeVdn
  
    
    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_interface(R)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: R
        END FUNCTION potential_interface
    END INTERFACE
    
    
    CONTAINS
    
    SUBROUTINE ComputePHPBoundStates(x, dstate, e_test, m, potential, Z, l, Ntotal, Nstart, eigE, eigFunc, Nfound, overlaps, defects)
        REAL(KIND = idk), INTENT(IN)                            :: x(:), dstate(:), e_test(:), m, Z
        INTEGER, INTENT(IN)                                     :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN)     :: potential
        INTEGER, INTENT(IN)                                     :: Ntotal, Nstart
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)              :: eigE(:), eigFunc(:,:), overlaps(:), defects(:)
        INTEGER, INTENT(OUT)                                    :: Nfound
        
        INTEGER, PARAMETER                                      :: npts_check = 15
        INTEGER                                                 :: max_log_candidates
        
        REAL(KIND = idk) :: dx
        INTEGER :: i, j, N, Ne, Ndstate, Nadd, Nfound1
        CHARACTER(LEN=256) :: message

        INTEGER, ALLOCATABLE :: root_indices(:) 
        INTEGER :: N_potential_roots, k, N_good_logged
        INTEGER, ALLOCATABLE :: root_status(:)  
        REAL(KIND = idk), ALLOCATABLE :: approx_E(:)  
        LOGICAL, ALLOCATABLE :: root_mask(:)  
        INTEGER, ALLOCATABLE :: good_root_indices(:) 
        INTEGER :: Nfound_good 
        
        REAL(KIND = idk) :: Eleft, Eright, Emid, fleft, fright, fmid, Etotal
        REAL(KIND = idk) :: EHfull, defectHfull
        INTEGER :: ilast, iHfull, kH
        REAL(KIND = idk), ALLOCATABLE :: floc(:), Eloc(:), psi_local(:)
        REAL(KIND = idk) :: mean_abs, edge_abs, ratio
        
        REAL(KIND = idk) :: mu_diff, mu_low, mu_high, mu_conv
        REAL(KIND = idk) :: E1, E2
        
        REAL(KIND = idk) :: S_x, S_y, S_xy, S_xx, y_j, x_j, x_bar, y_bar, B_fit, A_fit
        INTEGER          :: k_fit, j_fit

        
        max_log_candidates = Nstart
        Nadd = Ntotal - Nstart
        IF (Nadd < 0) THEN
            CALL CONSOLE('[ERROR]: Ntotal must be >= Nstart.')
            RETURN
        END IF
        
        Nfound = 0
        N = SIZE(x)
        dx = ABS(x(2)-x(1))
        Ne = SIZE(e_test)

        Ndstate = N
        DO i = N, 1, -1
            IF (dstate(i) /= 0.0d0) THEN
                Ndstate = i
                EXIT
            END IF
        END DO
        
        CALL CONSOLE('Searching for bound PHP states...')
        
        IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
        IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
        IF (ALLOCATED(overlaps)) DEALLOCATE(overlaps)
        IF (ALLOCATED(defects)) DEALLOCATE(defects)
        ALLOCATE(eigE(Ntotal))
        ALLOCATE(eigFunc(N, Ntotal)) 
        ALLOCATE(overlaps(Ne))
        ALLOCATE(defects(Ntotal))
        
        eigE =  ieee_value(0.0d0, ieee_quiet_nan)
        eigFunc = ieee_value(0.0d0, ieee_quiet_nan) 
        defects = ieee_value(0.0d0, ieee_quiet_nan)
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(DYNAMIC, 1) 
        DO i = 1, Ne
            CALL compute_overlap(x, Ndstate, dstate, e_test(i), m, l, potential, Z, overlaps(i))
        END DO
        !$OMP END PARALLEL DO
        
        
        CALL CONSOLE('Locating potential root intervals...')
        ALLOCATE(root_mask(Ne-1))
        DO i = 1, Ne-1
            root_mask(i) = (overlaps(i)*overlaps(i+1) < 0.0d0 .OR. overlaps(i) == 0.0d0)
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

        ALLOCATE(root_status(N_potential_roots))
        ALLOCATE(approx_E(N_potential_roots))
        
        root_status = -1
        approx_E = 0.0d0

        CALL CONSOLE('Filtering all potential roots...')
        
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP&   PRIVATE(k, i, Eleft, Eright, j, floc, Eloc, mean_abs, edge_abs, ratio) &
        !$OMP&   SCHEDULE(DYNAMIC, 1)
        DO k = 1, N_potential_roots
            i = root_indices(k)
            Eleft = e_test(i)
            Eright = e_test(i+1)
            approx_E(k) = 0.5d0 * (Eleft + Eright) 

            ALLOCATE(floc(npts_check), Eloc(npts_check))
            DO j = 1, npts_check
                Eloc(j) = Eleft + (Eright - Eleft) * (REAL(j-1, idk) / REAL(npts_check-1, idk))
                CALL compute_overlap(x, Ndstate, dstate, Eloc(j), m, l, potential, Z, floc(j))
            END DO
            
            mean_abs = SUM(ABS(floc)) / REAL(npts_check, idk)
            edge_abs = MAX(ABS(floc(1)), ABS(floc(npts_check)))
            ratio = 0.0d0
            IF (edge_abs > 1.0d-100) ratio = mean_abs / edge_abs 
            DEALLOCATE(floc, Eloc)

            IF (ratio > 0.9d0) THEN
                root_status(k) = 2  
            ELSE
                root_status(k) = 0
            END IF
        END DO
        !$OMP END PARALLEL DO
        
        
        good_root_indices = PACK(root_indices, mask=(root_status == 0))
        Nfound_good = SIZE(good_root_indices)

        WRITE(message, '(A,I0,A)') 'Found ', Nfound_good, ' PHP state candidates.'
        CALL CONSOLE(message)

        N_good_logged = 0
        DO k = 1, N_potential_roots
            SELECT CASE (root_status(k))
                CASE (0) 
                    IF (N_good_logged < max_log_candidates) THEN
                        WRITE(message, '(A,F0.4,A)') '  -> PHP state candidate found at E ~ ', approx_E(k)*phys_h0, ' eV'
                        CALL CONSOLE(message)
                        N_good_logged = N_good_logged + 1
                    ELSE 
                        EXIT
                    END IF
                CASE (2) 
                    WRITE(message, '(A,F0.4,A)') '  -> H bound state skipped at E ~ ', approx_E(k)*phys_h0, ' eV'
                    CALL CONSOLE(message)
            END SELECT
        END DO
        
        Nfound = MIN(Nfound_good, Nstart)
        
        IF (Nfound == 0) THEN
             DEALLOCATE(good_root_indices)
             WRITE(message, '(A,I0,A)') 'Computation of bound states finished. ', Nfound, ' bound states found.'
             CALL CONSOLE(message)
             RETURN
        END IF
        

        WRITE(message, '(A,I0,A)') 'Starting parallel bisection for ', Nfound, ' states...'
        CALL CONSOLE(message)
        
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP&   PRIVATE(k, i, Eleft, Eright, Emid, fleft, fright, fmid, j, psi_local) &
        !$OMP&   SCHEDULE(DYNAMIC, 1)
        DO k = 1, Nfound
            i = good_root_indices(k) 
            
            Eleft = e_test(i)
            Eright = e_test(i+1)
            fleft = overlaps(i)
            fright = overlaps(i+1)
            
            DO j = 1, max_iter
                Emid = 0.5d0 * (Eleft + Eright)
                CALL compute_overlap(x, Ndstate, dstate, Emid, m, l, potential, Z, fmid)
                
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
        
            CALL compute_normalized_states(x, dstate, eigE(k), m, l, potential, Z, psi_local)
            
            eigFunc(:, k) = psi_local(:) 
            DEALLOCATE(psi_local)
            
        END DO
        !$OMP END PARALLEL DO
        Nfound1 = Nfound
        
        ilast = good_root_indices(Nfound1)
        DEALLOCATE(good_root_indices)
        iHfull = -1
        DO kH = 1, N_potential_roots
            IF (root_indices(kH) > ilast .AND. root_status(kH)==2) THEN
                iHfull = root_indices(kH)
                EXIT
            END IF
        END DO
        DEALLOCATE(root_indices)
        
        IF (iHfull == -1) THEN
            CALL CONSOLE('[ERROR]: Could not find bracketing H-state. e_test grid too small?')
            DEALLOCATE(root_status, approx_E, defects)
            RETURN
        END IF
        
        WRITE(message, '(A,F0.4,A)') 'Found bracketing H-state near E ~ ', approx_E(kH)*phys_h0, ' eV. Converging...'
        CALL CONSOLE(message)
        DEALLOCATE(root_status, approx_E)
        
        Eleft = e_test(iHfull)
        Eright = e_test(iHfull + 1)
        fleft = overlaps(iHfull)
        fright = overlaps(iHfull + 1)
        
        DO j = 1, max_iter
            Emid = 0.5d0 * (Eleft + Eright)
            CALL compute_overlap(x, Ndstate, dstate, Emid, m, l, potential, Z, fmid)
            IF (fleft * fmid <= 0.0d0) THEN
                Eright = Emid; fright = fmid
            ELSE
                Eleft = Emid; fleft = fmid
            END IF
        END DO
        EHfull = 0.5d0 * (Eleft + Eright)
        
        defectHfull = REAL(Nfound1 + 1, idk) - SQRT(0.5d0/ABS(EHfull))
        WRITE(message, '(A,F0.4,A,F0.6)') ' H-state E = ', EHfull*phys_h0, ' eV, with defect mu = ', defectHfull
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
                    
                        CALL compute_normalized_states(x, dstate, E1, m, l, potential, Z, psi_local)
                        eigFunc(:, j) = psi_local(:)
                        DEALLOCATE(psi_local)
                    END DO
                    !$OMP END PARALLEL DO

                    Nfound = Ntotal
                    EXIT 
                END IF
            END IF

        
            mu_diff = ABS(defects(i-1) - defectHfull)
            mu_low = defects(i-1) - 0.04d0 * mu_diff
            mu_high = defects(i-1) + 0.06d0 * mu_diff
        
            E1 = -0.5d0 / (REAL(i, idk) - mu_low)**2
            E2 = -0.5d0 / (REAL(i, idk) - mu_high)**2
            Eleft = MIN(E1, E2)
            Eright = MAX(E1, E2)
        
            CALL compute_overlap(x, Ndstate, dstate, Eleft, m, l, potential, Z, fleft)
            CALL compute_overlap(x, Ndstate, dstate, Eright, m, l, potential, Z, fright)
        
            IF (fleft * fright > 0.0d0) THEN
                WRITE(message, '(A,I0,A)') '[WARNING]: 5% defect bracket failed for n = ', i, '. Widening to 50%.'
                CALL CONSOLE(message)
                mu_low = defects(i-1) - 0.5d0 * mu_diff
                mu_high = defects(i-1) + 0.5d0 * mu_diff
                E1 = -0.5d0 / (REAL(i, idk) - mu_low)**2; E2 = -0.5d0 / (REAL(i, idk) - mu_high)**2
                Eleft = MIN(E1, E2); Eright = MAX(E1, E2)
            
                CALL compute_overlap(x, Ndstate, dstate, Eleft, m, l, potential, Z, fleft)
                CALL compute_overlap(x, Ndstate, dstate, Eright, m, l, potential, Z, fright)

                IF (fleft * fright > 0.0d0) THEN
                     WRITE(message, '(A,I0,A)') 'Error: Cannot bracket root for n = ', i, '. Stopping.'
                     CALL CONSOLE(message)
                     Nfound = i-1
                     EXIT 
                END IF
            END IF
        
            DO j = 1, max_iter2
                Emid = 0.5d0 * (Eleft + Eright)
                CALL compute_overlap(x, Ndstate, dstate, Emid, m, l, potential, Z, fmid)
                IF (fleft * fmid <= 0.0d0) THEN
                    Eright = Emid; fright = fmid
                ELSE
                    Eleft = Emid; fleft = fmid
                END IF
            END DO
            eigE(i) = 0.5d0 * (Eleft + Eright)
        
            defects(i) = REAL(i, idk) - SQRT(0.5d0/ABS(eigE(i)))
            CALL compute_normalized_states(x, dstate, eigE(i), m, l, potential, Z, psi_local)
            eigFunc(:, i) = psi_local(:)
            DEALLOCATE(psi_local)
        
            Nfound = i 
            
            WRITE(message, '(A,I0,A,F0.8)') 'PHP state with n = ', Nfound, ' found with mu = ', defects(i)
            CALL CONSOLE(message)
        
        END DO 
        
        WRITE(message, '(A,I0,A)') 'Computation of bound states finished. ', Nfound, ' bound states found.'
        CALL CONSOLE(message)
    
    END SUBROUTINE ComputePHPBoundStates
    
    
    SUBROUTINE compute_overlap(x, Ndstate, dstate, E, m, l, potential, Z, overlap)
        REAL(KIND = idk), INTENT(IN)                        :: x(:), dstate(:), E, m, Z
        INTEGER, INTENT(IN)                                 :: l
        INTEGER, INTENT(IN)                                 :: Ndstate
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), INTENT(OUT)                       :: overlap
    
        
        REAL(KIND = idk), ALLOCATABLE       :: gdstate(:)
        REAL(KIND = idk), ALLOCATABLE       :: overlapstate(:)
        REAL(KIND = idk)                    :: dx
        INTEGER :: i, Nmin
        
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
        
    END SUBROUTINE compute_overlap
    
    
    
    SUBROUTINE compute_normalized_states(x, dstate, E, mass, l, potential, Z, psi)
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
        
        REAL(KIND = idk), PARAMETER         :: eps = 1.0d-25
        
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
                    dr_adaptive = 2.0d0 * PI / local_k / 25.0d0
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
        
        norm_tail = norm_tail * A_scale_sq * dx
        norm = norm + norm_tail
        norm = SQRT(ABS(norm))
        
        IF (ALLOCATED(psi)) DEALLOCATE(psi)
        ALLOCATE(psi(SIZE(x)))
        psi = 0.0d0
        DO i = 1, SIZE(x)
            psi(i) = gdstate(i) / norm
        END DO
        
        DEALLOCATE(gdstate) 
        
    END SUBROUTINE compute_normalized_states
    
    
    
    
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
  
  
  
  
  
  
END MODULE PHPBound