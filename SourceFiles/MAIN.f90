!-------------------------------------< MAIN >-------------------------------------!
!                                                                                  !    
! Contains: Main program                                                           !
!                                                                                  !
! Last revision:    25/09/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!

PROGRAM MAIN
    USE Parameters
    USE DiscreteState
    USE COULCC_M, only: COULCC
    USE whittaker_w, only: coulomb_whittaker
    USE Green
    USE RydbergSolver
    USE GDTrans
    USE, INTRINSIC :: ieee_arithmetic

    IMPLICIT NONE
    REAL(KIND = idk), ALLOCATABLE :: x(:), e(:), V_params(:)
    REAL(KIND = idk), ALLOCATABLE :: Delta(:), Gamma2(:), DeltaA(:), Gamma2A(:), Phaseshift(:), DS_phaseshift(:), Vde(:)
    REAL(KIND = idk), ALLOCATABLE :: DeltaFull(:,:), Gamma2Full(:,:), DeltaAFull(:,:), Gamma2AFull(:,:), PhaseshiftFull(:,:), DS_phaseshiftFull(:,:), VdeFull(:,:)
    REAL(KIND = idk), ALLOCATABLE :: DeltaContinuous(:), DeltaContinuousFull(:,:), DeltaRydbergFull(:,:), DeltaRydberg(:)
    REAL(KIND = idk), ALLOCATABLE :: dstate_bound(:), dstate(:,:), freestate(:), moddedfreestate(:), realscatstate(:), xext(:), psi_R(:), psi_I(:), k2(:), wvalues(:), wdvalues(:), DSenergies(:)
    REAL(KIND = idk), ALLOCATABLE :: boundEgrid(:), eigEFull_H(:,:), defects_H(:,:), defect_H(:), logdevsFull(:,:), logdevs(:), eigE_H(:)
    REAL(KIND = idk), ALLOCATABLE :: eigEFull_PHP(:,:), defects_PHP(:,:), eigE_PHP(:), eigFunc(:,:), defect_PHP(:)
    REAL(KIND = idk), ALLOCATABLE :: VdnFull(:,:), Vdn(:)
    COMPLEX(KIND = idk), ALLOCATABLE :: cmplxscatstate(:)
    INTEGER :: status, iounit, i, j, Negrid, jj, Nfull, l
    REAL(KIND = idk) :: dx, dE, dV, norm, k, x_cut, wronski, h, Ascale, E_dstate_asymptotic
    CHARACTER(LEN=256) :: message
    REAL(KIND=idk) :: xxx, eta, w, wd
    COMPLEX(KIND=idk) :: XX, ETA1, ZLMIN
    COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
    INTEGER :: NL, MODE1, KFN, IFAIL, sf, sfscale, status_dstate
    CHARACTER(LEN=200) :: filename
    INTEGER :: istat

    REAL(KIND = idk)                    :: V_A = 1.0d0
    
    ! Interface for R->R function
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.0d0)) FUNCTION real_function_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.0d0)), INTENT(IN) :: x
        END FUNCTION real_function_interface
    END INTERFACE 

    PROCEDURE(real_function_interface), POINTER :: V_ptr, V_asymptotic_ptr, dstate_ptr



    V_ptr => V_2D_2
    V_asymptotic_ptr => V_2D_2_asymptotic
    dstate_ptr => dstate1
    


    
    print*, '=============================================================='
    print*, '           DISCRETE STATE IN CONTINUUM COMPUTATION            '
    print*, '=============================================================='
    
    
    ! Parameters output
    CALL parameters_o()
    

    ! Computational coordinate grid
    CALL make_mesh(xmin, xmax, mp, x, dx)
    CALL CONSOLE('Coordinate computational grid created successfully.')

    ! Energy computational grid
    CALL make_mesh(Emin, Emax, ep, e, dE)
    CALL CONSOLE('Energy computational grid created successfully.')
    
    ! Potential parameter grid
    !CALL make_mesh(Vmin, Vmax, nv, V_params, dV)
    CALL make_mesh_ndyn(Vmin, Vmax, 1.5d0, 4.0d0, nv, V_params)
    CALL CONSOLE('Potential parameters grid created successfully.')
    
    
    
    ! Create directories
    INQUIRE(DIRECTORY="DATA", EXIST=status)
    IF (.NOT. status) THEN
        CALL SYSTEM("mkdir -p DATA")
        CALL CONSOLE('Directory "DATA" created successfully.')
    ELSE IF (status == -1) THEN
        CALL CONSOLE('Directory "DATA" already exists.')
    END IF


    ! Creating initial discrete state
    IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
    ALLOCATE(dstate(nv,SIZE(x)))
    IF (dstate_boundstate) THEN
        CALL compute_dstate(x, m, l_ang, Z, V_asymptotic_ptr, 2.0d0*Ebound_min, dstate_bound, E_dstate_asymptotic, status_dstate)
        IF (status_dstate < 0) THEN
            CALL CONSOLE('[ERROR]: Discrete state computation failed. Check the parameters and potential settings.') 
            STOP
        ELSE
            DO j = 1, nv
                dstate(j,:) = dstate_bound(:)
            END DO
            IF (ALLOCATED(dstate_bound)) DEALLOCATE(dstate_bound)
            CALL CONSOLE('Initial discrete state created successfully.')
        END IF
    ELSE
        DO j = 1, nv
            V_A = V_params(j)
            DO i = 1, SIZE(x)
                dstate(j,i) = dstate_ptr(x(i))
            END DO
        norm = SUM( ABS(dstate(j,:))**2 ) * dx
        dstate(j,:) = dstate(j,:) / SQRT(norm)
        END DO
        CALL CONSOLE('Initial discrete state created successfully.')
    END IF
    
    
    !=================================================================================
    !                     Rydberg states in Coulomb-like potential
    !=================================================================================
    IF (Rydbergstates) THEN
    IF (ABS(Z) > 1.0d-10) THEN
        INQUIRE(DIRECTORY="DATA/RydbergStates", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/RydbergStates")
            CALL CONSOLE('Subdirectory "RydbergStates" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "RydbergStates" already exists.')
        END IF
        
        INQUIRE(DIRECTORY="DATA/RydbergStates/Eigenfunctions", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/RydbergStates/Eigenfunctions")
            CALL CONSOLE('Subdirectory "Eigenfunctions" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "Eigenfunctions" already exists.')
        END IF
        
        ! Printing V(x) to file
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/V.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','V(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            DO j = 1, nv
                V_A = V_params(j)
                WRITE(iounit, '(1E20.12)', advance='no') V_ptr(x(i)) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/V.txt')

        ! Printing V_A parameters to file
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/Rvalues.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17)') '#', 'V_A_params [au]'
        DO j = 1, nv
            WRITE(iounit, '(1E20.12)') V_params(j)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/Rvalues.txt')

    
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') dstate(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/dstate.txt')
        
        
        
        CALL rydberg_grid(Nbound, NperN, gridtail, Ebound_min, boundEgrid)
        Negrid = SIZE(boundEgrid)
        
        IF (ALLOCATED(eigEFull_H)) DEALLOCATE(eigEFull_H)
        IF (ALLOCATED(defects_H)) DEALLOCATE(defects_H)
        IF (ALLOCATED(logdevsFull)) DEALLOCATE(logdevsFull)
        IF (ALLOCATED(eigEFull_PHP)) DEALLOCATE(eigEFull_PHP)
        IF (ALLOCATED(defects_PHP)) DEALLOCATE(defects_PHP)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)

        
        ALLOCATE(eigEFull_H(nv,Nbound))
        ALLOCATE(defects_H(nv,Nbound))
        ALLOCATE(logdevsFull(nv,Negrid))
        ALLOCATE(eigEFull_PHP(nv,Nbound))
        ALLOCATE(defects_PHP(nv,Nbound))
        ALLOCATE(VdnFull(nv,Nbound))
        
        DO j = 1, nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing in short-range potential with V_A = ', V_A, ' ...'
            CALL CONSOLE(message)
            
            IF(ALLOCATED(eigE_H)) DEALLOCATE(eigE_H)
            IF(ALLOCATED(eigE_PHP)) DEALLOCATE(eigE_PHP)
            IF(ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
            IF(ALLOCATED(logdevs)) DEALLOCATE(logdevs)
            IF(ALLOCATED(defect_H)) DEALLOCATE(defect_H)
            IF(ALLOCATED(defect_PHP)) DEALLOCATE(defect_PHP)
            
            CALL ComputeRydbergSystem(x, dstate(j,:), boundEgrid, m, Z, l_ang, V_ptr, Nbound, Nstart, N_fit_points, logdevs, eigE_H, eigE_PHP, eigFunc, defect_H, defect_PHP)   
            
            eigEFull_H(j,:) = eigE_H
            defects_H(j,:) = defect_H
            eigEFull_PHP(j,:) = eigE_PHP
            defects_PHP(j,:) = defect_PHP
            logdevsFull(j,:) = logdevs
            
            ! Printing eigenfunctions
            WRITE(filename, '(F0.3)') V_A
            filename = 'DATA/RydbergStates/Eigenfunctions/eigFunc_VA_' // TRIM(ADJUSTL(filename)) // '.txt'
            OPEN(NEWUNIT=iounit, FILE=filename, STATUS='replace', ACTION='write')
            WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'R [au]','|psi_n>'
            WRITE(iounit,*)
            DO jj = 1, SIZE(x)
                WRITE(iounit, '(1E20.12)', advance='no') x(jj)
                DO i = 1, MIN(Nbound,Nprint)
                    WRITE(iounit, '(1E20.12)', advance='no') eigFunc(jj,i)
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            
            IF(ALLOCATED(eigE_H)) DEALLOCATE(eigE_H)
            IF(ALLOCATED(eigE_PHP)) DEALLOCATE(eigE_PHP)
            IF(ALLOCATED(logdevs)) DEALLOCATE(logdevs)
            IF(ALLOCATED(defect_H)) DEALLOCATE(defect_H)
            IF(ALLOCATED(defect_PHP)) DEALLOCATE(defect_PHP)
            
            ALLOCATE(Vdn(Nbound))
            CALL ComputeVdn(x, dstate(j,:), m, V_ptr, eigFunc, Vdn)
            VdnFull(j,:) = Vdn
    
            IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
            IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
            
        END DO

        
        
        
        ! Printing log-derivatives on test grid
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/logderivs.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','d/dr log(u) [1/Bohr]'
        WRITE(iounit,*)
        DO i = 1, Negrid
            WRITE(iounit, '(I15)', advance='no') i
            WRITE(iounit, '(1E20.12)', advance='no') boundEgrid(i) * phys_h0
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') TANH(logdevsFull(j,i))
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/logderivs.txt')
        
        ! Printing eigenenergies
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/eigE_H.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') eigEFull_H(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/eigE_H.txt')

        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/eigE_PHP.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') eigEFull_PHP(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/eigE_PHP.txt')
        
        
        ! Printing defects
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/defects_H.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','mu_n'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') defects_H(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/defects_H.txt')
        
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/defects_PHP.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','mu_n'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') defects_PHP(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/defects_PHP.txt')    

        ! Printing Vdn
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/Vdn.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n', 'Vdn [eV]'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') VdnFull(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/RydbergStates/Vdn.txt')

        ! Saving Vdn and eigE to bin files for later use
        CALL CONSOLE('Caching V_params to binary file...')
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/V_params.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) V_params
        CLOSE(iounit)
        CALL CONSOLE('Caching eigEFull_H to binary files...')
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/eigE_H.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) eigEFull_H
        CLOSE(iounit)
        CALL CONSOLE('Caching eigEFull_PHP and VdnFull to binary files...')
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/eigE_PHP.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) eigEFull_PHP
        CLOSE(iounit)
        OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/Vdn.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) VdnFull
        CLOSE(iounit)
        CALL CONSOLE('Binary cache saved successfully.')

        
        IF (ALLOCATED(eigEFull_H)) DEALLOCATE(eigEFull_H)
        IF (ALLOCATED(defects_H)) DEALLOCATE(defects_H)
        IF (ALLOCATED(logdevsFull)) DEALLOCATE(logdevsFull)
        IF (ALLOCATED(eigEFull_PHP)) DEALLOCATE(eigEFull_PHP)
        IF (ALLOCATED(defects_PHP)) DEALLOCATE(defects_PHP)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        
        
    ELSE
        CALL CONSOLE('[ERROR]: Z is zero or very close to zero. Skipping Rydberg state computation.') 
    END IF   
    END IF
    
    
    
    !=================================================================================
    !                           Discrete state computation
    !=================================================================================
    IF (dstate_computation) THEN
        INQUIRE(DIRECTORY="DATA/DSState", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/DSState")
            CALL CONSOLE('Subdirectory "DSState" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "DSState" already exists.')
        END IF
    

       ! Printing V(x) to file
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/V.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','V(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            DO j = 1, nv
                V_A = V_params(j)
                WRITE(iounit, '(1E20.12)', advance='no') V_ptr(x(i)) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/V.txt')

        ! Printing V_A parameters to file
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/Rvalues.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17)') '#', 'V_A_params [au]'
        DO j = 1, nv
            WRITE(iounit, '(1E20.12)') V_params(j)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/Rvalues.txt')
        
        ! Printing V Short Range b to file
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/V_SR.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','V_SR(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            DO j = 1, nv
                V_A = V_params(j)
                WRITE(iounit, '(1E20.12)', advance='no') (V_ptr(x(i))+Z/x(i)-REAL(l_ang*(l_ang+1),KIND=idk)/(2.0d0*m*x(i)**2)* hbar**2) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/V_SR.txt')

    
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') dstate(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/dstate.txt')
    
    
        ! Computation of Gamma & Delta
        IF (ALLOCATED(DeltaFull)) DEALLOCATE(DeltaFull)
        IF (ALLOCATED(Gamma2Full)) DEALLOCATE(Gamma2Full)
        IF (ALLOCATED(DeltaAFull)) DEALLOCATE(DeltaAFull)
        IF (ALLOCATED(Gamma2AFull)) DEALLOCATE(Gamma2AFull)
        IF (ALLOCATED(PhaseshiftFull)) DEALLOCATE(PhaseshiftFull)
        IF (ALLOCATED(DS_PhaseshiftFull)) DEALLOCATE(DS_PhaseshiftFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        IF (ALLOCATED(DSenergies)) DEALLOCATE(DSenergies)
        ALLOCATE(DeltaFull(nv,ep), Gamma2Full(nv,ep), DeltaAFull(nv,ep), Gamma2AFull(nv,ep), PhaseshiftFull(nv,ep), DS_PhaseshiftFull(nv,ep), VdeFull(nv,ep), DSenergies(nv))
    
        DO j = 1, nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing in Coulomb-like potential with V_A = ', V_A, ' ...'
            CALL CONSOLE(message)
            CALL compute_DSenergy(x, m, V_ptr, dstate(j,:), DSenergies(j))
            CALL compute_DS(x, e, m, Z, l_ang, V_ptr, dstate(j,:), Delta, Gamma2, DeltaA, Gamma2A, Phaseshift, DS_phaseshift, Vde)
            DeltaFull(j,:) = Delta
            Gamma2Full(j,:) = Gamma2
            DeltaAFull(j,:) = DeltaA
            Gamma2AFull(j,:) = Gamma2A
            IF (unwrap) THEN
                CALL unwrap_phaseshift(Phaseshift)
                CALL unwrap_phaseshift(DS_phaseshift)
            END IF
            PhaseshiftFull(j,:) = Phaseshift
            DS_PhaseshiftFull(j,:) = DS_phaseshift
            VdeFull(j,:) = Vde
            IF (ALLOCATED(Delta)) DEALLOCATE(Delta)
            IF (ALLOCATED(Gamma2)) DEALLOCATE(Gamma2)
            IF (ALLOCATED(DeltaA)) DEALLOCATE(DeltaA)
            IF (ALLOCATED(Gamma2A)) DEALLOCATE(Gamma2A)
            IF (ALLOCATED(Phaseshift)) DEALLOCATE(Phaseshift)
            IF (ALLOCATED(DS_phaseshift)) DEALLOCATE(DS_phaseshift)
            IF (ALLOCATED(Vde)) DEALLOCATE(Vde)
        END DO
        
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/DSenergies.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A3, A21)', advance='no') '#', 'V_A','E_DS [eV]'
        WRITE(iounit,*)
        DO j = 1, nv
            WRITE(iounit, '(F0.3, 1E20.12)') V_params(j), DSenergies(j) * phys_h0
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/DSenergies.txt')
        
    
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/delta.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Delta [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') DeltaFull(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/delta.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/gamma.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Gamma2 [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') Gamma2Full(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/gamma.txt')

        IF (ABS(Z) > 1.0d-10) THEN
            OPEN(NEWUNIT=iounit, FILE='DATA/DSState/deltaA.txt', STATUS='replace', ACTION='write')
            WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Delta [eV]'
            WRITE(iounit,*)
            DO i = 1, SIZE(e)
                IF (e(i) < 0.0d0 .AND. e(i) > cutoff_energy(l_ang)) THEN
                    WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                    DO j = 1, nv
                        WRITE(iounit, '(1E20.12)', advance='no') DeltaAFull(j,i) * phys_h0
                    END DO
                    WRITE(iounit,*)
                END IF
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Data successfully written to DATA/DSState/deltaA.txt')
        
            OPEN(NEWUNIT=iounit, FILE='DATA/DSState/gammaA.txt', STATUS='replace', ACTION='write')
            WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Gamma2 [eV]'
            WRITE(iounit,*)
            DO i = 1, SIZE(e)
                IF (e(i) < 0.0d0 .AND. e(i) > cutoff_energy(l_ang)) THEN
                    WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                    DO j = 1, nv
                        WRITE(iounit, '(1E20.12)', advance='no') Gamma2AFull(j,i) * phys_h0
                    END DO
                    WRITE(iounit,*)
                END IF
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Data successfully written to DATA/DSState/gammaA.txt')
        END IF

    
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/phaseshift.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','phaseshift'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) >= 0.0d0) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') PhaseshiftFull(j,i)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/phaseshift.txt')

        IF (ABS(Z) > 1.0d-10) THEN
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/defect.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','defect_phase'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) < 0.0d0 .AND. e(i) > cutoff_energy(l_ang)) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') PhaseshiftFull(j,i)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/defect.txt')
        ELSE
            CALL CONSOLE('Z is zero or very close to zero. Skipping defect phase output.')
        END IF
    
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/DS_phaseshift.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','DS_phaseshift'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) >= 0.0d0) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') DS_phaseshiftFull(j,i)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/DS_phaseshift.txt')

        IF (ABS(Z) > 1.0d-10) THEN
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/DS_defect.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','DS_defect_phase'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) < 0.0d0 .AND. e(i) > cutoff_energy(l_ang)) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') DS_phaseshiftFull(j,i)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/DS_defect.txt')
        ELSE
            CALL CONSOLE('Z is zero or very close to zero. Skipping DS defect phase output.')
        END IF
    
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/Vde.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Vde [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') VdeFull(j,i) * SQRT(phys_h0)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/DSState/Vde.txt')
        
        ! Saving Vde to bin files for later use
        CALL CONSOLE('Caching VdeFull to binary files...')
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/Vde.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) VdeFull
        CLOSE(iounit)
        CALL CONSOLE('Binary cache saved successfully.')
        
        
        IF (ALLOCATED(DeltaFull)) DEALLOCATE(DeltaFull)
        IF (ALLOCATED(Gamma2Full)) DEALLOCATE(Gamma2Full)
        IF (ALLOCATED(PhaseshiftFull)) DEALLOCATE(PhaseshiftFull)
        IF (ALLOCATED(DS_PhaseshiftFull)) DEALLOCATE(DS_PhaseshiftFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
    
    END IF
    

    
    
    
    !=================================================================================
    !                       Hilbert transform Vde -> Delta 
    !=================================================================================
    IF (hilbert_computation) THEN
        
        INQUIRE(DIRECTORY="DATA/Hilbert", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/Hilbert")
            CALL CONSOLE('Subdirectory "Hilbert" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "Hilbert" already exists.')
        END IF
        
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        IF (ALLOCATED(DeltaContinuousFull)) DEALLOCATE(DeltaContinuousFull)
        IF (ALLOCATED(eigEFull_PHP)) DEALLOCATE(eigEFull_PHP)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        IF (ALLOCATED(DeltaRydbergFull)) DEALLOCATE(DeltaRydbergFull)
        
        CALL CONSOLE('Loading eigEFull, VdnFull & VdeFull from binary files...')
        
        IF (ABS(Z) > 1.0d-10) THEN

            ALLOCATE(eigEFull_PHP(nv,Nbound))
            OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/eigE_PHP.bin', FORM='unformatted', STATUS='old', ACTION='read', IOSTAT=istat)
            IF (istat /= 0) THEN
                CALL CONSOLE('[ERROR]: File eigE_PHP.bin cannot be found or openned.')
                STOP 'Error reading eigE_PHP.bin'
            END IF
            READ(iounit) eigEFull_PHP
            CLOSE(iounit)

            ALLOCATE(VdnFull(nv,Nbound))
            OPEN(NEWUNIT=iounit, FILE='DATA/RydbergStates/Vdn.bin', &
                FORM='unformatted', STATUS='old', ACTION='read', IOSTAT=istat)
            IF (istat /= 0) THEN
                CALL CONSOLE('[ERROR]: File Vdn.bin cannot be found or openned.')
                STOP 'Error reading Vdn.bin'
            END IF
            READ(iounit) VdnFull
            CLOSE(iounit)

        END IF
        
        ALLOCATE(VdeFull(nv,SIZE(e)))
        OPEN(NEWUNIT=iounit, FILE='DATA/DSState/Vde.bin', &
             FORM='unformatted', STATUS='old', ACTION='read', IOSTAT=istat)
        IF (istat /= 0) THEN
            CALL CONSOLE('[ERROR]: File Vde.bin cannot be found or openned.')
            STOP 'Error reading Vde.bin'
        END IF
        READ(iounit) VdeFull
        CLOSE(iounit)

        CALL CONSOLE('Binary cache loaded successfully.')
        
        
        IF (ABS(Z) > 1.0d-10) THEN
            ! Printing eigenenergies
            OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/En.txt', STATUS='replace', ACTION='write')
            WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
            WRITE(iounit,*)
            DO i = 1, Nbound
                WRITE(iounit, '(I3)', advance='no') i
                DO j = 1, nv 
                    WRITE(iounit, '(1E20.12)', advance='no') eigEFull_PHP(j,i) * phys_h0
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Data successfully written to DATA/Hilbert/En.txt')
            
            ! Printing Vdn
            OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/Vdn.txt', STATUS='replace', ACTION='write')
            WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n', 'Vdn [SQRT(eV)]'
            WRITE(iounit,*)
            DO i = 1, Nbound
                WRITE(iounit, '(I3)', advance='no') i
                DO j = 1, nv 
                    WRITE(iounit, '(1E20.12)', advance='no') VdnFull(j,i) * phys_h0
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Data successfully written to DATA/Hilbert/Vdn.txt')
        END IF
        
        ! Printing Vde
        OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/Vde.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Vde [SQRT(eV)]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') VdeFull(j,i) * SQRT(phys_h0)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Hilbert/Vde.txt')
        
        
        ! Computining Delta continuous by Hilbert transform
        ALLOCATE(DeltaContinuousFull(nv,SIZE(e)))
        ALLOCATE(DeltaRydbergFull(nv,SIZE(e)))
        DO j = 1,nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing Hilbert transform for V_A = ', V_A, '...'
            CALL CONSOLE(message)
            IF (ALLOCATED(DeltaContinuous)) DEALLOCATE(DeltaContinuous)
            IF (ALLOCATED(Vde)) DEALLOCATE(Vde)
            ALLOCATE(DeltaContinuous(SIZE(e)))
            ALLOCATE(Vde(SIZE(e)))
            DeltaContinuous(:) = DeltaContinuousFull(j,:)
            Vde(:) = VdeFull(j,:)
            CALL ComputeDeltaContinuous(Vde, SIZE(e), ABS(e(2)-e(1)), DeltaContinuous)
            DeltaContinuousFull(j,:) = DeltaContinuous
            IF (ALLOCATED(DeltaContinuous)) DEALLOCATE(DeltaContinuous)
            IF (ALLOCATED(Vde)) DEALLOCATE(Vde)
            WRITE(message, '(A,F0.3,A)') 'Computation of Hilbert transform for V_A = ', V_A, ' completed.'
            CALL CONSOLE(message)
            IF (ABS(Z) > 1.0d-10) THEN
                WRITE(message, '(A,F0.3,A)') 'Computing Rydberg-part of Delta for V_A = ', V_A, '...'
                CALL CONSOLE(message)
                IF (ALLOCATED(DeltaRydberg)) DEALLOCATE(DeltaRydberg)
                IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
                IF (ALLOCATED(eigE_PHP)) DEALLOCATE(eigE_PHP)
                ALLOCATE(DeltaRydberg(SIZE(e)))
                ALLOCATE(Vdn(Nbound))
                ALLOCATE(eigE_PHP(Nbound))
                Vdn(:) = VdnFull(j,:)
                eigE_PHP(:) = eigEFull_PHP(j,:)
                CALL ComputeDeltaRydberg(e, Vdn, eigE_PHP, DeltaRydberg)
                DeltaRydbergFull(j,:) = DeltaRydberg
                IF (ALLOCATED(DeltaRydberg)) DEALLOCATE(DeltaRydberg)
                IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
                IF (ALLOCATED(eigE_PHP)) DEALLOCATE(eigE_PHP)      
                WRITE(message, '(A,F0.3,A)') 'Computation of Rydberg-part of Delta for V_A = ', V_A, ' completed.'
                CALL CONSOLE(message)
            ELSE
                DeltaRydbergFull(j,:) = 0.0d0
            END IF
        END DO
        
        IF (ABS(Z) > 1.0d-10) THEN
            ! Printing Delta Continuous
            OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/DeltaContinuous.txt', STATUS='replace', ACTION='write')
            WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','DeltaC [eV]'
            WRITE(iounit,*)
            DO i = 1, SIZE(e)
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') DeltaContinuousFull(j,i) * phys_h0
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Data successfully written to DATA/Hilbert/DeltaContinuous.txt')
            
            
            ! Printing Delta Rydberg
            OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/DeltaRydberg.txt', STATUS='replace', ACTION='write')
            WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','DeltaC [eV]'
            WRITE(iounit,*)
            DO i = 1, SIZE(e)
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') DeltaRydbergFull(j,i) * phys_h0
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Data successfully written to DATA/Hilbert/DeltaRydberg.txt')
        END IF
        
        
        ! Printing Delta Full
        OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/DeltaFull.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','DeltaC [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
            DO j = 1, nv
                WRITE(iounit, '(1E20.12)', advance='no') (DeltaRydbergFull(j,i)+DeltaContinuousFull(j,i))* phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Hilbert/DeltaFull.txt')
        
        
        
        IF (ALLOCATED(eigEFull_PHP)) DEALLOCATE(eigEFull_PHP)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        IF (ALLOCATED(DeltaContinuousFull)) DEALLOCATE(DeltaContinuousFull)
        IF (ALLOCATED(DeltaRydbergFull)) DEALLOCATE(DeltaRydbergFull)
        
        
    END IF
    
    
    
    

    ! Release allocated resources
    IF (ALLOCATED(x)) DEALLOCATE(x)
    IF (ALLOCATED(e)) DEALLOCATE(e)
    IF (ALLOCATED(V_params)) DEALLOCATE(V_params)
    IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
    
    CONTAINS
    
    ! Potential V(R)
    REAL(KIND = idk) FUNCTION V(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        V = ( 0.5d0 * R**2 - 1.5d0 + V_A ) / ( 1.0d0 + EXP( (4.0d0-2.0d0*V_A) * (R-2.0d0-2.0d0*V_A) ) )
    END FUNCTION V
    
    
    
    ! Potential V_2(R)
    REAL(KIND = idk) FUNCTION V_coulomb(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk) :: sR
        sR = 1.0d0 / ( 1.0d0 + EXP( (4.0d0-2.0d0*V_A) * (R- 2.25d0-1.35d0*V_A) ) )
        V_coulomb = ( 0.5d0 * R**2 - 1.5d0 + V_A ) * sR - (1.0 - sR) / (R + 2.0d0)
    END FUNCTION V_coulomb
    
    
    ! Potential V_II(R)
    REAL(KIND = idk) FUNCTION V_coulomb2(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk) :: sR
        sR = 1.0d0 / ( 1.0d0 + EXP( (4.0d0-2.0d0*V_A) * (R- 2.35d0-1.15d0*V_A) ) )
        V_coulomb2 = ( 0.5d0 * R**2 - 1.5d0 + V_A ) * sR - (1.0 - sR) * Z / R
    END FUNCTION V_coulomb2






    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_2D_2(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk), PARAMETER :: a = 0.4d0
        REAL(KIND = idk), PARAMETER :: b = 2.0d0
        REAL(KIND = idk), PARAMETER :: c = 1.5d0
        REAL(KIND = idk), PARAMETER :: d = 0.72d0
        REAL(KIND = idk), PARAMETER :: rc = 5.0d0
        REAL(KIND = idk), PARAMETER :: rrc = 1.5d0
        V_2D_2 = - 1.0d0 / R + 0.5d0 * REAL(l_ang*(l_ang+1),KIND=idk) / R**2 + a * EXP(-(R-rc)**2/b**2) - d * EXP(-R**2/4.0d0) * TANH( (V_A-rrc)/c )
    END FUNCTION V_2D_2

    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_2D_2_asymptotic(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk), PARAMETER :: a = 0.4d0
        REAL(KIND = idk), PARAMETER :: b = 2.0d0
        REAL(KIND = idk), PARAMETER :: c = 1.5d0
        REAL(KIND = idk), PARAMETER :: d = 0.72d0
        REAL(KIND = idk), PARAMETER :: rc = 5.0d0
        V_2D_2_asymptotic = - 1.0d0 / R + 0.5d0 * REAL(l_ang*(l_ang+1),KIND=idk) / R**2 + a * EXP(-(R-rc)**2/b**2) - d * EXP(-R**2/4.0d0)
    END FUNCTION V_2D_2_asymptotic




    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_2D_1(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk), PARAMETER :: a = 0.5d0
        REAL(KIND = idk), PARAMETER :: b = 1.2d0
        REAL(KIND = idk), PARAMETER :: c = 1.8d0
        REAL(KIND = idk), PARAMETER :: d = 0.6519d0
        REAL(KIND = idk), PARAMETER :: rc = 4.5d0
        REAL(KIND = idk), PARAMETER :: rrc = 1.5d0
        V_2D_1 = - 1.0d0 / R + 0.5d0 * REAL(l_ang*(l_ang+1),KIND=idk) / R**2 + a * EXP(-(R-rc)**2/b**2) - d * EXP(-R**2/4.0d0) * TANH( (V_A-rrc)/c )
    END FUNCTION V_2D_1

    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_2D_1_asymptotic(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk), PARAMETER :: a = 0.5d0
        REAL(KIND = idk), PARAMETER :: b = 1.2d0
        REAL(KIND = idk), PARAMETER :: c = 1.8d0
        REAL(KIND = idk), PARAMETER :: d = 0.6519d0
        REAL(KIND = idk), PARAMETER :: rc = 4.5d0
        V_2D_1_asymptotic = - 1.0d0 / R + 0.5d0 * REAL(l_ang*(l_ang+1),KIND=idk) / R**2 + a * EXP(-(R-rc)**2/b**2) - d * EXP(-R**2/4.0d0)
    END FUNCTION V_2D_1_asymptotic
    



    
    ! Pure Coulombic potential for testing 
    REAL(KIND = idk) FUNCTION V_pure_coulomb(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        V_pure_coulomb = -1.0d0 / (R + 0.000001d0)
    END FUNCTION V_pure_coulomb



    ! Dstate pure harmonic 1
    REAL(KIND = idk) FUNCTION dstate1(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        dstate1 = R * EXP(-R**2 / 2.0d0)
    END FUNCTION dstate1

    
END PROGRAM MAIN
  