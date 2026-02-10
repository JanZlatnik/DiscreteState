!-----------------------------------< GD Trans >-----------------------------------!
!                                                                                  !    
! Contains: Computation of Hilbert transform from Gamma to Delta                   !                                                   
!                                                                                  !
! Last revision:    26/10/2025                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!

INCLUDE 'mkl_dfti.f90'
    
MODULE GDTrans
    USE Parameters
    USE Math
    USE MKL_DFTI
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ComputeDeltaContinuous, ComputeDeltaRydberg
    
    CONTAINS
    
    
    SUBROUTINE ComputeDeltaContinuous(Vde_in, N_in, dE, Delta_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN)             :: N_in
        REAL(KIND = idk), INTENT(IN)    :: dE
        REAL(KIND = idk), INTENT(IN)    :: Vde_in(N_in)
        REAL(KIND = idk), INTENT(OUT)   :: Delta_out(N_in)
        
        INTEGER :: N, M, k_size, i, status
        INTEGER :: MKL_SUCCES = 0
        
        REAL(KIND = idk), ALLOCATABLE       :: f_temp(:), f_padded(:), delta_padded(:), W_k(:)
        REAL(KIND = idk)                    :: scale_factor
        COMPLEX(KIND =  idk)                :: ipineg 
        COMPLEX(KIND = idk), ALLOCATABLE    :: F_k(:), D_k(:)
        TYPE(DFTI_DESCRIPTOR), POINTER      :: desc_handle
        CHARACTER(LEN=256)                  :: message
        
        
        N = N_in
        M = 2 * N_in
        
        ALLOCATE(f_temp(N))
        ALLOCATE(f_padded(M))
        
        DO i = 1,N
            f_temp(i) = Vde_in(i) ** 2
        END DO
        
        f_padded(1:N) = f_temp(1:N)
        f_padded(N+1:M) = 0.0d0
        DEALLOCATE(f_temp)
        
        k_size = M / 2 + 1
        ALLOCATE(F_k(k_size))
        
        status = DftiCreateDescriptor(desc_handle, DFTI_DOUBLE, DFTI_REAL, 1, M)
        IF (status /= MKL_SUCCES) THEN
            WRITE(message, '(A,I0)') '[ERROR]: DFTI_CREATE_DESCRIPTOR failed with status ', status
            CALL CONSOLE(message)
        END IF
        
        status = DftiSetValue( desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        IF (status /= MKL_SUCCES) THEN
            WRITE(message, '(A,I0)') '[ERROR]: DFTI_SET_VALUE failed with status ', status
            CALL CONSOLE(message)
        END IF
        
        status = DftiCommitDescriptor(desc_handle)
        IF (status /= MKL_SUCCES) THEN
            WRITE(message, '(A,I0)') '[ERROR]: DFTI_COMMIT_DESCRIPTOR failed with status ', status
            CALL CONSOLE(message)
        END IF
        
        status = DftiComputeForward(desc_handle, f_padded, F_k)
        IF (status /= MKL_SUCCES) THEN
            WRITE(message, '(A,I0)') '[ERROR]: DFTI_COMPUTE_FORWARD failed with status ', status
            CALL CONSOLE(message)
        END IF
        
        DEALLOCATE(f_padded)
        
        ALLOCATE(D_k(k_size))
        ipineg = CMPLX(0.0d0, -pi, KIND=idk)
        D_k(1) = CMPLX(0.0d0, 0.0d0, KIND=idk)
        D_k(k_size) = CMPLX(0.0d0, 0.0d0, KIND=idk)
        ALLOCATE(W_k(k_size))
        DO i = 1,k_size
            W_k(i) = 0.5d0 * (1.0d0 + COS(pi * REAL(i-1)/REAL(M/2)))
        END DO
        
        DO i = 2,k_size-1
            D_k(i) =  F_k(i) * ipineg  * W_k(i)
        END DO
        
        DEALLOCATE(F_k)
        
        ALLOCATE(delta_padded(M))
        
        status = DftiComputeBackward(desc_handle, D_k, delta_padded)
        IF (status /= MKL_SUCCES) THEN
            WRITE(message, '(A,I0)') '[ERROR]: DFTI_COMPUTE_BACKWARD failed with status ', status
            CALL CONSOLE(message)
        END IF
        
        DEALLOCATE(D_k)
        scale_factor = 1/REAL(M,idk)
        delta_padded = delta_padded * scale_factor
        Delta_out(1:N) = delta_padded(1:N)
        DEALLOCATE(delta_padded)
        
        status = DftiFreeDescriptor(desc_handle)
        IF (status /= MKL_SUCCES) THEN
            WRITE(message, '(A,I0)') '[ERROR]: DFTI_FREE_DESCRIPTOR failed with status ', status
            CALL CONSOLE(message)
        END IF
    
    
    END SUBROUTINE ComputeDeltaContinuous
    
    
    
    
    
    
    SUBROUTINE ComputeDeltaRydberg(egrid, Vdn, En, Delta)
        REAL(KIND = idk), INTENT(IN)                :: egrid(:), Vdn(:), En(:)    
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)  :: Delta(:)
    
        INTEGER, PARAMETER                  :: Nfit = 20
        INTEGER, PARAMETER                  :: Ncorrectionmult = 10
        
        INTEGER :: i, j, k, N, Nmax, Nkorekce
        REAL(KIND = idk) :: ee, Delta_exact, Delta_tail_near, Delta_tail_far
    
        REAL(KIND = idk) :: mu0, B, A_avg, n_star, En_fit, Vdn_fit, mu_fit_k
        REAL(KIND = idk) :: S_x, S_y, S_xy, S_xx, y_j, x_j, x_bar, y_bar
        REAL(KIND = idk) :: n_star_start_far, E_cutoff_far
        REAL(KIND = idk), ALLOCATABLE :: mu_fit(:), A_fit(:)
        
        CHARACTER(LEN=256) :: message
    
        IF(ALLOCATED(Delta)) DEALLOCATE(Delta)
        Nmax = SIZE(Vdn)
        N = SIZE(egrid)
        ALLOCATE(Delta(N))
        Nkorekce = Nmax * Ncorrectionmult
        
        WRITE(message, '(A,I0,A,I0)') 'Performing QDT fit on n=', Nmax - Nfit + 1, ' to ', Nmax
        CALL CONSOLE(message)
        
        S_x = 0.0d0
        S_y = 0.0d0
        S_xy = 0.0d0
        S_xx = 0.0d0        
    
        ALLOCATE(mu_fit(Nfit), A_fit(Nfit))
        DO k = 1, Nfit
            j = Nmax - Nfit + k
            mu_fit(k) = j - SQRT(-0.5d0 / En(j))
            A_fit(k)  = Vdn(j)*Vdn(j) * (REAL(j, idk)**3)
            
            y_j = mu_fit(k)
            x_j = 1.0d0 / (REAL(j,idk)**2)
                        
            S_x = S_x + x_j
            S_y = S_y + y_j
            S_xy = S_xy + x_j*y_j
            S_xx = S_xx + x_j*x_j
        END DO
                    
        x_bar = S_x / REAL(Nfit,idk)
        y_bar = S_y / REAL(Nfit,idk)
                    
        B = (S_xy - REAL(Nfit,idk)*x_bar*y_bar) / (S_xx - REAL(Nfit, idk) * x_bar * x_bar)
        mu0 = y_bar - B * x_bar
        A_avg  = SUM(A_fit) / Nfit
        DEALLOCATE(mu_fit, A_fit)
        
        WRITE(message, '(A,F0.8,A,F0.8,A)') 'Fit: mu(n) = A + B/n^2. A (mu_inf) = ', mu0, ', B = ', B
        CALL CONSOLE(message)
        WRITE(message, '(A,F0.8)') 'Fit: Vdn^2(n) = A/n^3. A = ', A_avg * phys_h0 * phys_h0
        CALL CONSOLE(message)
    
        
        n_star_start_far = (Nmax + Nkorekce) - mu0
        E_cutoff_far = -0.5d0 / (n_star_start_far * n_star_start_far)
    
    
        !$OMP PARALLEL DO PRIVATE(j, k, ee, Delta_exact, Delta_tail_near, Delta_tail_far, &
        !$OMP&                    n_star, En_fit, Vdn_fit, mu_fit_k)
        DO i = 1,N
            ee = egrid(i)

            Delta_exact = 0.0d0
            DO j = 1,Nmax
                Delta_exact = Delta_exact + Vdn(j)*Vdn(j)/(ee-En(j))
            END DO
    
            Delta_tail_near = 0.0d0
            DO k = Nmax + 1, Nmax + Nkorekce
        
                mu_fit_k = mu0 + B / (REAL(k, idk)**2) 
                n_star   = REAL(k, idk) - mu_fit_k
                En_fit   = -0.5d0 / (n_star * n_star)
                Vdn_fit  = A_avg / (REAL(k, idk)**3)
        
                Delta_tail_near = Delta_tail_near + Vdn_fit / (ee - En_fit)
            END DO
    
            Delta_tail_far = A_avg * LOG((ee - E_cutoff_far) / ee)
            
            Delta(i) = Delta_exact + Delta_tail_near + Delta_tail_far
    
        END DO
        !$OMP END PARALLEL DO
    
    END SUBROUTINE
    
    
    
    
    
END MODULE GDTrans