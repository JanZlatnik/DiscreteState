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
    PUBLIC ::compute_DS, unwrap_phaseshift, compute_coulomb_DS, compute_DSenergy

    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_interface(R)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: R
        END FUNCTION potential_interface
    END INTERFACE
    
    
    CONTAINS
    
    ! Subroutine to compute discrete state properties
    SUBROUTINE compute_DS(x, E, mass, potential, dstate, Delta, Gamma2, phaseshift, DS_phaseshift, Vde)
        REAL(KIND = idk), INTENT(IN) :: x(:), E(:), mass, dstate(:)
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: Delta(:), Gamma2(:), phaseshift(:), DS_phaseshift(:), Vde(:)
        
        INTEGER :: i, j, N
        COMPLEX(KIND = idk), ALLOCATABLE :: green_dstate(:), green_mat(:,:), scattered_state(:)
        COMPLEX(KIND = idk) :: greendot, scattdot
        REAL(KIND = idk), ALLOCATABLE :: ham_dstate(:), mat(:,:), freestate(:)
        REAL(KIND = idk) :: tail, k, dx, hamdot
        CHARACTER(LEN=50) :: message
        
        dx = ABS(x(2) - x(1)) 
        N = SIZE(x)
        tail = 2.0d0
        
        ALLOCATE(Delta(SIZE(E)))
        ALLOCATE(Gamma2(SIZE(E)))
        ALLOCATE(phaseshift(SIZE(E)))
        ALLOCATE(DS_phaseshift(SIZE(E)))
        ALLOCATE(Vde(SIZE(E)))
        Delta = 0.0d0
        Gamma2 = 0.0d0
        phaseshift = 0.0d0
        DS_phaseshift = 0.0d0
        Vde = 0.0d0
        
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(j, i, k, hamdot, greendot, scattdot, tail, &
        !$OMP         ham_dstate, green_dstate, mat, green_mat, freestate, scattered_state, message) &
        !$OMP SHARED(x, E, mass, N, dstate, dx, Delta, Gamma2, phaseshift, DS_phaseshift, Vde) &
        !$OMP SCHEDULE(DYNAMIC)
        DO j = 1, SIZE(E)
            
            hamdot = 0.0d0
            greendot = (0.0d0, 0.0d0)
            scattdot = (0.0d0, 0.0d0)
            
            ALLOCATE(ham_dstate(N))
            ham_dstate = 0.0d0
            ham_dstate(1) = - (hbar**2 / (2.0d0 * mass)) * (dstate(2) - 2.0d0*dstate(1)) / (dx**2) + potential(x(1)) * dstate(1)
            ham_dstate(1) = ham_dstate(1) * dstate(1)
            ham_dstate(N) = - (hbar**2 / (2.0d0 * mass)) * (dstate(N) - 2.0d0*dstate(N) + dstate(N-1)) / (dx**2) + potential(x(N)) * dstate(N)
            ham_dstate(N) = ham_dstate(N) * dstate(N)
            DO i = 2, N-1
                ham_dstate(i) = - (hbar**2 / (2.0d0 * mass)) * (dstate(i+1) - 2.0d0*dstate(i) + dstate(i-1)) / (dx**2) + potential(x(i)) * dstate(i)
                ham_dstate(i) = ham_dstate(i) * dstate(i)
            END DO
            CALL definite_integral(ham_dstate, dx, hamdot)
            DEALLOCATE(ham_dstate)  
    
            ALLOCATE(mat(2,N))
            mat = 0.0d0
            
            mat(1,:) = dstate
            
            IF (E(j) >= 0.0d0) THEN
                ALLOCATE(freestate(N))
                freestate = 0.0d0
                k = SQRT(2.0d0*mass*E(j))/hbar
                DO i = 1, N
                    freestate(i) = SIN(k*x(i)) * SQRT(2*mass/(pi*k*hbar**2))
                    mat(2,i) = freestate(i) * potential(x(i))
                END DO
            END IF
            
                 
            tail = 1.0d0
            DO i = 1, 100
                IF (ABS(potential(x(N)*tail)/E(j)) < 1.0d-2) THEN
                    EXIT
                END IF
                tail = tail + 1.0d0
            END DO
            IF (i >= 100000) THEN
                !$OMP CRITICAL
                WRITE(message, '(A,F0.4,A)') '[WARNING]: Tail not found for E = ', E(j)*phys_h0, ' eV'
                CALL CONSOLE(message)
                !$OMP END CRITICAL
            END IF
            
            ALLOCATE(green_mat(2,N))
            CALL apply_green(x, E(j), mass, tail, potential, mat, green_mat)
            
            DEALLOCATE(mat)
            ALLOCATE(green_dstate(N))
            IF (E(j) >= 0.0d0) THEN
                ALLOCATE(scattered_state(N))
                scattered_state = freestate + green_mat(2,:)
                DEALLOCATE(freestate)
            END IF
            green_dstate = green_mat(1,:)
            DEALLOCATE(green_mat)
            
            DO i = 1, N
                green_dstate(i) = green_dstate(i) * dstate(i)
            END DO
            CALL definite_integral(green_dstate, dx, greendot)
            greendot = 1.0d0 / greendot
            DEALLOCATE(green_dstate)
            
            Delta(j) = E(j) - hamdot - REAL(greendot)
            Gamma2(j) = 2.0d0 * AIMAG(greendot)
            
            IF (E(j) >= 0.0d0) THEN
                CALL compute_phaseshift(scattered_state, x, k, phaseshift(j))
                DS_phaseshift(j) = - ATAN2(AIMAG(greendot),REAL(greendot))
                DO i = 1, N
                    scattered_state(i) = scattered_state(i) * dstate(i)
                END DO
                CALL definite_integral(scattered_state, dx, scattdot)
                DEALLOCATE(scattered_state)
                Vde(j) = ABS(scattdot * greendot)
            ELSE
                phaseshift(j) = ieee_value(0.0d0, ieee_quiet_nan)
                ds_phaseshift(j) = ieee_value(0.0d0, ieee_quiet_nan)
                Vde(j) = 0.0d0
            END IF
            
            !$OMP CRITICAL
            WRITE(message, '(A,F0.4,A)') 'Computation finsished for E = ', E(j)*phys_h0, ' eV'
            CALL CONSOLE(message)
            !$OMP END CRITICAL
            
        END DO
        !$OMP END PARALLEL DO
    
    END SUBROUTINE compute_DS
    
    
    
    
    ! Subroutine to compute discrete state properties in Coulomb potential
    SUBROUTINE compute_coulomb_DS(x, E, mass, potential, Z, dstate, Delta, Gamma2, phaseshift, DS_phaseshift, Vde)
        USE OMP_LIB
        REAL(KIND = idk), INTENT(IN) :: x(:), E(:), mass, dstate(:), Z
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: Delta(:), Gamma2(:), phaseshift(:), DS_phaseshift(:), Vde(:)
        
        INTEGER :: i, j, N
        COMPLEX(KIND = idk), ALLOCATABLE :: green_dstate(:), scattered_state(:), greened_mat(:,:)
        REAL(KIND = idk), ALLOCATABLE :: ham_dstate(:), freestate(:), mat(:,:), auxstate(:)
        REAL(KIND = idk) :: k, dx, hamdot
        COMPLEX(KIND=idk) :: XX, ETA1, greendot, scattdot
        COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
        INTEGER :: IFAIL
        
        dx = ABS(x(2) - x(1)) 
        N = SIZE(x)
  
        ALLOCATE(Delta(SIZE(E)))
        ALLOCATE(Gamma2(SIZE(E)))
        ALLOCATE(phaseshift(SIZE(E)))
        ALLOCATE(DS_phaseshift(SIZE(E)))
        ALLOCATE(Vde(SIZE(E)))
        Delta = 0.0d0
        Gamma2 = 0.0d0
        phaseshift = 0.0d0
        DS_phaseshift = 0.0d0
        Vde = 0.0d0
        
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(j, i, hamdot, greendot, scattdot, &
        !$OMP     ham_dstate, green_dstate, mat, greened_mat, freestate, &
        !$OMP     scattered_state, k, XX, ETA1, FC, GC, FCP, GCP, SIG, IFAIL) &
        !$OMP SHARED(x, dstate, E, mass, Z, dx, Delta, Gamma2, phaseshift, DS_phaseshift, Vde, N) &
        !$OMP SCHEDULE(DYNAMIC)
        DO j = 1, SIZE(E)

            hamdot = 0.0d0
            greendot = (0.0d0, 0.0d0)
            scattdot = (0.0d0, 0.0d0)
            
            ALLOCATE(ham_dstate(N))
            ham_dstate = 0.0d0
            ham_dstate(1) = - (hbar**2 / (2.0d0 * mass)) * (dstate(2) - 2.0d0*dstate(1)) / (dx**2) + potential(x(1)) * dstate(1)
            ham_dstate(1) = ham_dstate(1) * dstate(1)
            ham_dstate(N) = - (hbar**2 / (2.0d0 * mass)) * (dstate(N) - 2.0d0*dstate(N) + dstate(N-1)) / (dx**2) + potential(x(N)) * dstate(N)
            ham_dstate(N) = ham_dstate(N) * dstate(N)
            DO i = 2, N-1
                ham_dstate(i) = - (hbar**2 / (2.0d0 * mass)) * (dstate(i+1) - 2.0d0*dstate(i) + dstate(i-1)) / (dx**2) + potential(x(i)) * dstate(i)
                ham_dstate(i) = ham_dstate(i) * dstate(i)
            END DO
            CALL definite_integral(ham_dstate, dx, hamdot)
            DEALLOCATE(ham_dstate)  
            
            
            
            if (e(j) >= 0.0d0) then
                allocate(mat(2,N))
                mat = 0.0d0
                k = sqrt(2.0d0*mass*e(j))/hbar
                eta1 = cmplx(- z * m / (hbar**2 * k), 0.0d0, kind=idk)
                allocate(freestate(N))
                freestate = 0.0d0
                do i = 1, N
                    !freestate(i) = sin(k*x(i)) * sqrt(2*mass/(pi*k*hbar**2))
                    !mat(2,i) = freestate(i) * potential(x(i))
                    ifail = 0
                    xx = cmplx(k * x(i), 0.0d0, kind=idk)
                    call coulcc(xx, eta1, (0.0d0,0.0d0), 1, fc, gc, fcp, gcp, sig, 12, 0, ifail)
                    freestate(i) = REAL(fc(1)) * sqrt(2*mass/(pi*k*hbar**2))
                    mat(2,i) = freestate(i) * (potential(x(i)) + z /(x(i)+creg) )
                end do
                allocate(greened_mat(2,N))
            else
                allocate(mat(1,N))
                mat = 0.0d0
                allocate(greened_mat(1,N))
            end if
            
            mat(1,:) = dstate(:)
            
            call apply_green_coulomb(x, e(j), mass, z, potential, mat, greened_mat)
            deallocate(mat)
            
            allocate(green_dstate(N))
            do i = 1, N
                green_dstate(i) = greened_mat(1,i) * dstate(i)
            end do
            if (e(j) >= 0.0d0) then
                allocate(scattered_state(N))
                scattered_state = freestate + greened_mat(2,:)
                deallocate(freestate)
            end if
            deallocate(greened_mat)
            call definite_integral(green_dstate, dx, greendot)
            greendot = 1.0d0 / greendot
            deallocate(green_dstate)
        
            
            delta(j) = e(j) - hamdot - real(greendot)
            gamma2(j) = 2.0d0 * aimag(greendot)
            
            if (e(j) >= 0.0d0) then
                call compute_coulomb_phaseshift(scattered_state, x, k, z, mass, phaseshift(j))
                ds_phaseshift(j) = - atan2(aimag(greendot),real(greendot))
                do i = 1, N
                    scattered_state(i) = scattered_state(i) * dstate(i)
                end do
                call definite_integral(scattered_state, dx, scattdot)
                deallocate(scattered_state)
                vde(j) = abs(scattdot * greendot)
            else
                call compute_analytic_defect(x, e(j), mass, potential, Z, phaseshift(j))
                allocate(green_dstate(N))
                call apply_green_analytic_continuation(x, e(j), mass, Z, potential, dstate, green_dstate)
                do i = 1, N
                    green_dstate(i) = green_dstate(i) * dstate(i)
                end do
                call definite_integral(green_dstate, dx, greendot)
                deallocate(green_dstate)
                greendot = 1.0d0 / greendot
                ds_phaseshift(j) = atan2(aimag(greendot),real(greendot))
                vde(j) = 0.0d0
            end if
            
        END DO
        !$OMP END PARALLEL DO
    
    END SUBROUTINE compute_coulomb_DS
    
    
    
    
    
    
    
    
    ! Subroutine to extract phase shift from scattered state
    SUBROUTINE compute_phaseshift(psi, x, k, delta)
      IMPLICIT NONE
      COMPLEX(KIND = idk), INTENT(IN)  :: psi(:)     ! scattered wavefunction
      REAL(KIND = idk),    INTENT(IN)  :: x(:)       ! grid positions
      REAL(KIND = idk),    INTENT(IN)  :: k          ! momentum
      REAL(KIND = idk),    INTENT(OUT) :: delta      ! phase shift

      INTEGER, PARAMETER :: xpick = 20   ! grid offset
      INTEGER, PARAMETER :: xn    = 10   ! halfwidth
      COMPLEX(KIND = idk) :: M11, M12, M21, M22, R1, R2
      COMPLEX(KIND = idk) :: A_in, A_out, b1, b2, S
      INTEGER :: N, i0, i1, i

      N  = SIZE(psi)
      i0 = MAX(1, N - xpick - xn)
      i1 = MIN(N, N - xpick + xn)

      M11 = (0.0d0,0.0d0); M12 = (0.0d0,0.0d0)
      M21 = (0.0d0,0.0d0); M22 = (0.0d0,0.0d0)
      R1  = (0.0d0,0.0d0); R2  = (0.0d0,0.0d0)

      DO i = i0, i1
         b1 = EXP(CMPLX(0.0d0,-1.0d0, KIND=idk) * k * x(i))  
         b2 = EXP(CMPLX(0.0d0, 1.0d0, KIND=idk) * k * x(i)) 
         M11 = M11 + CONJG(b1)*b1
         M12 = M12 + CONJG(b1)*b2
         M21 = M21 + CONJG(b2)*b1
         M22 = M22 + CONJG(b2)*b2
         R1  = R1  + CONJG(b1)*psi(i)
         R2  = R2  + CONJG(b2)*psi(i)
      END DO


      A_in  = ( M22*R1 - M12*R2 )
      A_out = (-M21*R1 + M11*R2 )
      S     = -A_out / A_in

      delta = 0.5d0 * ATAN2(AIMAG(S), REAL(S))

    END SUBROUTINE compute_phaseshift
    
    ! Subroutine to extract phase shift from scattered state in Coulomb potential
    SUBROUTINE compute_coulomb_phaseshift(psi, x, k, Z, m, delta)
      IMPLICIT NONE
      COMPLEX(KIND = idk), INTENT(IN)  :: psi(:)
      REAL(KIND = idk),    INTENT(IN)  :: x(:)
      REAL(KIND = idk),    INTENT(IN)  :: k
      REAL(KIND = idk),    INTENT(IN)  :: Z
      REAL(KIND = idk),    INTENT(IN)  :: m
      REAL(KIND = idk),    INTENT(OUT) :: delta

      INTEGER, PARAMETER :: xpick = 20, xn = 10
      COMPLEX(KIND = idk) :: M11, M12, M21, M22, R1, R2
      COMPLEX(KIND = idk) :: A_in, A_out, b_minus, b_plus, S
      INTEGER :: N, i0, i1, i
      COMPLEX(KIND=idk) :: XX, ETA1
      COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
      INTEGER :: IFAIL
      REAL(KIND=idk) :: eta

      N  = SIZE(psi)
      i0 = MAX(1, N - xpick - xn)
      i1 = MIN(N, N - xpick + xn)

      M11 = (0.0d0,0.0d0); M12 = (0.0d0,0.0d0)
      M21 = (0.0d0,0.0d0); M22 = (0.0d0,0.0d0)
      R1  = (0.0d0,0.0d0); R2  = (0.0d0,0.0d0)

      eta = - Z * m / (hbar**2 * k)
      ETA1 = CMPLX(eta, 0.0d0, KIND=idk)
      IFAIL = 1

      DO i = i0, i1
         XX = CMPLX(k * x(i), 0.0d0, KIND=idk)
         CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
         b_plus  = GC(1)
         CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 22, 0, IFAIL)
         b_minus = GC(1)

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

    END SUBROUTINE compute_coulomb_phaseshift



! Subroutine to compute analytical continuation of quantum defect
    SUBROUTINE compute_analytic_defect(x, E_val, mass, potential, Z, defect_phase)
        USE OMP_LIB
        IMPLICIT NONE
        
        REAL(KIND = idk), INTENT(IN)  :: x(:)
        REAL(KIND = idk), INTENT(IN)  :: E_val, mass, Z
        REAL(KIND = idk), INTENT(OUT) :: defect_phase
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential

        INTEGER, PARAMETER :: xpick = 20, xn = 10
        INTEGER :: N, i, ifail, i0, i1
        REAL(KIND = idk) :: dx, kappa, nu
        REAL(KIND = idk), ALLOCATABLE :: psi(:), k2(:)
        REAL(KIND = idk), PARAMETER :: pi_const = 3.14159265358979_idk
        
        COMPLEX(KIND=idk) :: XX, ETA_C, i_unit
        COMPLEX(KIND=idk) :: F_val, G_val
        COMPLEX(KIND=idk) :: H_plus, H_minus
        COMPLEX(KIND=idk) :: M11, M12, M21, M22, R1, R2
        COMPLEX(KIND=idk) :: A_coef, B_coef, S_ratio
        COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
        
        N = SIZE(x)
        dx = ABS(x(2) - x(1))
        
        ! 1. Numerov propagace
        ALLOCATE(psi(N)); ALLOCATE(k2(N))
        DO i = 1, N
            k2(i) = 2.0d0 * mass / (hbar**2) * (E_val - potential(x(i)))
        END DO
        
        psi(1) = 0.0d0
        psi(2) = dx**2 
        DO i = 2, n-1
            psi(i+1) = (2.0d0*(1.0d0-5.0d0*dx**2*k2(i)/12.0d0)*psi(i) - (1.0d0+dx**2*k2(i-1)/12.0d0)*psi(i-1)) / (1.0d0+dx**2*k2(i+1)/12.0d0)
        END DO
        DEALLOCATE(k2)
        
        ! 2. Příprava parametrů pro COULCC
        kappa = SQRT(2.0d0 * mass * ABS(E_val)) / hbar
        nu = Z * mass / (hbar**2 * kappa)
        i_unit = CMPLX(0.0d0, 1.0d0, KIND=idk)
        
        ! 3. Nastavení fitovacího okna (analogicky k phaseshift rutině)
        i0 = MAX(1, N - xpick - xn)
        i1 = MIN(N, N - xpick + xn)
        
        M11 = (0.0d0, 0.0d0); M12 = (0.0d0, 0.0d0)
        M21 = (0.0d0, 0.0d0); M22 = (0.0d0, 0.0d0)
        R1  = (0.0d0, 0.0d0); R2  = (0.0d0, 0.0d0)
        
        ! 4. Smyčka fitování (Least Squares)
        DO i = i0, i1
            XX = CMPLX(0.0d0, kappa * x(i), KIND=idk) 
            ETA_C = CMPLX(0.0d0, nu, KIND=idk)
            
            ifail = 1
            CALL COULCC(XX, ETA_C, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, ifail)
            H_plus  = GC(1)
            CALL COULCC(XX, ETA_C, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 22, 0, ifail)
            H_minus = GC(1)
            
            ! Plnění matice soustavy
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


    
    
    
    
    
    
    SUBROUTINE compute_DSenergy(x, mass, potential, dstate, DSenergy)
        REAL(KIND = idk), INTENT(IN)                        :: x(:), mass, dstate(:)
        PROCEDURE(potential_interface), POINTER, INTENT(IN) :: potential
        REAL(KIND = idk), INTENT(OUT)                       :: DSenergy
        
        
        INTEGER :: i, N
        REAL(KIND = idk), ALLOCATABLE :: ham_dstate(:), dx
        
        N = SIZE(x)
        dx = ABS(x(2) - x(1))
        
        ALLOCATE(ham_dstate(N))
        ham_dstate = 0.0d0
        ham_dstate(1) = - (hbar**2 / (2.0d0 * mass)) * (dstate(2) - 2.0d0*dstate(1)) / (dx**2) + potential(x(1)) * dstate(1)
        ham_dstate(1) = ham_dstate(1) * dstate(1)
        ham_dstate(N) = - (hbar**2 / (2.0d0 * mass)) * (dstate(N) - 2.0d0*dstate(N) + dstate(N-1)) / (dx**2) + potential(x(N)) * dstate(N)
        ham_dstate(N) = ham_dstate(N) * dstate(N)
        DO i = 2, N-1
            ham_dstate(i) = - (hbar**2 / (2.0d0 * mass)) * (dstate(i+1) - 2.0d0*dstate(i) + dstate(i-1)) / (dx**2) + potential(x(i)) * dstate(i)
            ham_dstate(i) = ham_dstate(i) * dstate(i)
        END DO
        CALL definite_integral(ham_dstate, dx, DSenergy)
        DEALLOCATE(ham_dstate)  
    
    END SUBROUTINE compute_DSenergy





END MODULE DiscreteState