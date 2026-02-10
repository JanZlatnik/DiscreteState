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
    USE RadialShooting
    USE PHPBound
    USE GDTrans
    USE, INTRINSIC :: ieee_arithmetic

    IMPLICIT NONE
    REAL(KIND = idk), ALLOCATABLE :: x(:), e(:), V_params(:)
    REAL(KIND = idk), ALLOCATABLE :: Delta(:), Gamma2(:), Phaseshift(:), DS_phaseshift(:), Vde(:)
    REAL(KIND = idk), ALLOCATABLE :: DeltaFull(:,:), Gamma2Full(:,:), PhaseshiftFull(:,:), DS_phaseshiftFull(:,:), VdeFull(:,:)
    REAL(KIND = idk), ALLOCATABLE :: DeltaContinuous(:), DeltaContinuousFull(:,:), DeltaRydbergFull(:,:), DeltaRydberg(:)
    REAL(KIND = idk), ALLOCATABLE :: dstate(:), freestate(:), moddedfreestate(:), realscatstate(:), xext(:), psi_R(:), psi_I(:), k2(:), wvalues(:), wdvalues(:)
    REAL(KIND = idk), ALLOCATABLE :: boundEgrid(:), eigEFull(:,:), defects(:,:), logdevsFull(:,:), overlapsFull(:,:), eigE(:), eigFunc(:,:), logdevs(:), DSenergies(:), overlaps(:), VdnFull(:,:), Vdn(:), defect(:)
    COMPLEX(KIND = idk), ALLOCATABLE :: cmplxscatstate(:)
    INTEGER :: status, iounit, i, j, Negrid, Nfound, jj, Nfull, l
    REAL(KIND = idk) :: dx, dE, dV, norm, k, x_cut, wronski, h, Ascale
    CHARACTER(LEN=256) :: message
    REAL(KIND=idk) :: xxx, eta, w, wd
    COMPLEX(KIND=idk) :: XX, ETA1, ZLMIN
    COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
    INTEGER :: NL, MODE1, KFN, IFAIL, sf, sfscale
    CHARACTER(LEN=200) :: filename
    INTEGER :: istat

    REAL(KIND = idk)                    :: V_A = 1.0d0!-0.15d0
    
    ! Interface for R->R function
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.0d0)) FUNCTION real_function_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.0d0)), INTENT(IN) :: x
        END FUNCTION real_function_interface
    END INTERFACE 

    PROCEDURE(real_function_interface), POINTER :: V_ptr, dstateptr
    
    
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
    CALL make_mesh(Vmin, Vmax, nv, V_params, dV)
    CALL CONSOLE('Potential parameters grid created successfully.')
    
    
    
    ! Create directories
    INQUIRE(DIRECTORY="DATA", EXIST=status)
    IF (.NOT. status) THEN
        CALL SYSTEM("mkdir -p DATA")
        CALL CONSOLE('Directory "DATA" created successfully.')
    ELSE IF (status == -1) THEN
        CALL CONSOLE('Directory "DATA" already exists.')
    END IF
    
    !=================================================================================
    !                                    Tests
    !=================================================================================
    IF (test) THEN
        INQUIRE(DIRECTORY="DATA/Test", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/Test")
            CALL CONSOLE('Subdirectory "Test" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "Test" already exists.')
        END IF
        CALL CONSOLE('Testing Coulomb functions ...')
        
        V_ptr => V_coulomb2
        dstateptr => dstate1
        
        ETA1 = cmplx(-1.0d0 / SQRT(2.0d0*e_test), 0.0d0, KIND=idk)  
        eta = REAL(ETA1)
        ZLMIN = cmplx(0.0d0, 0.0d0, KIND=idk)
        NL    = 1
        MODE1 = 12    ! 11 for computing derivatives
        KFN   = 0
        IFAIL = 0
        
        
        
        !! === 
        !CALL coulomb_whittaker(eta, 0, SQRT(2.0d0*e_test)*x(SIZE(x)), w, wd, sf)
        !PRINT*, w
        !PRINT*, sf
        !CALL CONSOLE('Press ENTER to exit...')
        !READ(*,*)
        !STOP
        !
        !
        !! ===
        !
        
        
        
        x_cut = MAX((turnplus+1.0d0)*Z/e_test, xmax)
        Nfull = CEILING( (x_cut - x(1)) / dx ) + 1
        
        ALLOCATE(xext(Nfull))
        ALLOCATE(dstate(Nfull))
        
        xext(1) = x(1)
        DO i = 1, Nfull-1
            xext(i+1) = xext(i) + dx
        END DO

        
        DO i = 1, Nfull
            dstate(i) = dstateptr(xext(i))
        END DO
        norm = SUM( ABS(dstate)**2 ) * dx
        dstate = dstate / SQRT(norm)

        OPEN(NEWUNIT=iounit, FILE='DATA/Test/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, Nfull
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') dstate(i)
            WRITE(iounit, '(1E20.12)', advance='no') (V_ptr(xext(i)) + Z /(xext(i)+Creg) ) * phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/dstate.txt')
        
        CALL CONSOLE('Testing Coulomb H+ function...')

        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/coulomb_Hplus.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21, A21, A21)') '#', 'R [au]','Re(H+)','Im(H+)'
        
        DO i = 1, Nfull
            XX = cmplx(SQRT(2.0d0*e_test)*xext(i), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, ZLMIN, NL, FC, GC, FCP, GCP, SIG, MODE1, KFN, IFAIL)
            WRITE(iounit,'(3E20.12)') xext(i), REAL(GC(1)), AIMAG(GC(1))
        END DO
        
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/coulomb_Hplus.txt')
        
        CALL CONSOLE('Testing Whittaker W function...')
        
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/coulomb_whittakerW.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21, A21, A21)') '#', 'R [au]','W'
        
        ALLOCATE(wvalues(Nfull))
        ALLOCATE(wdvalues(Nfull))
        DO i = 1, Nfull
            CALL coulomb_whittaker(eta, 0, SQRT(2.0d0*e_test)*xext(i), w, wd, sf)
            wvalues(i) = w * 10.0d0**sf
            wdvalues(i) = wd * 10.0d0**sf
        END DO

        DO i = 1, Nfull
            WRITE(iounit,'(1E20.12)', advance='no') xext(i)
            WRITE(iounit,'(1E20.12)', advance='no') wvalues(i)
            WRITE(iounit,'(1E20.12)', advance='no') SQRT(2.0d0*e_test) * wdvalues(i)
            IF (i>1 .AND. i<Nfull) THEN
                WRITE(iounit,'(1E20.12)', advance='no') (wvalues(i+1) - wvalues(i-1)) / (2.0d0*dx)
            ELSE
                WRITE(iounit,'(1E20.12)', advance='no') ieee_value(0.0d0, ieee_quiet_nan)
            END IF
            WRITE(iounit,*)
        END DO

        DEALLOCATE(wvalues, wdvalues)
        
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/coulomb_whittakerW.txt')
        
        CALL CONSOLE('Testing Coulomb Green function for positive energies...')
        
        
        k = SQRT(2.0d0*m*e_test)/hbar
        ETA1 = CMPLX(- Z * m / (hbar**2 * k), 0.0d0, KIND=idk)      
        

        ALLOCATE(freestate(Nfull))
        ALLOCATE(moddedfreestate(Nfull))
        freestate = 0.0d0
        
        DO i = 1, Nfull
            IFAIL = 0
            XX = CMPLX(k * xext(i), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            freestate(i) = FC(1) * SQRT(2*m/(pi*k*hbar**2))
            moddedfreestate(i) = freestate(i) * (V_ptr(xext(i)) + Z /(xext(i)+Creg) )
        END DO
        
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/freestate+.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, Nfull
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') freestate(i)
            WRITE(iounit, '(1E20.12)', advance='no') moddedfreestate(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/freestate+.txt')
        
        
        ALLOCATE(cmplxscatstate(Nfull))
        CALL apply_green_coulomb(xext, e_test, m, Z, V_ptr, moddedfreestate, cmplxscatstate)
        cmplxscatstate = cmplxscatstate + freestate
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/scatstate+full.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, Nfull
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') REAL(cmplxscatstate(i))
            WRITE(iounit, '(1E20.12)', advance='no') AIMAG(cmplxscatstate(i))
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/scatstate+full.txt')   
        
        
        DEALLOCATE(cmplxscatstate)
        ALLOCATE(cmplxscatstate(mp))
        CALL apply_green_coulomb(xext(1:mp), e_test, m, Z, V_ptr, moddedfreestate(1:mp), cmplxscatstate)
        cmplxscatstate = cmplxscatstate + freestate(1:mp)
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/scatstate+short.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, mp
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') REAL(cmplxscatstate(i))
            WRITE(iounit, '(1E20.12)', advance='no') AIMAG(cmplxscatstate(i))
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/scatstate+short.txt')   
        
        
        DEALLOCATE(cmplxscatstate,freestate,moddedfreestate)
        
        CALL CONSOLE('Testing Coulomb Green function for negative energies...')
        
        
        ALLOCATE(realscatstate(Nfull))
        CALL apply_green_coulomb_bound(xext, -e_test, m, Z, V_ptr, dstate, realscatstate)
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/gdstatefull.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, Nfull
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') realscatstate(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/gdstatefull.txt')   
        
        
        DEALLOCATE(realscatstate)
        ALLOCATE(realscatstate(mp))
        CALL apply_green_coulomb_bound(xext(1:mp), -e_test, m, Z, V_ptr, dstate(1:mp), realscatstate)
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/gdstateshort.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, mp
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') realscatstate(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/gdstateshort.txt')   
        
        CALL coulomb_whittaker(eta, 0, SQRT(2.0d0*e_test)*xext(mp), w, wd, sfscale)
        Ascale = realscatstate(mp) / w
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/gdstateext.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = mp+1, Nfull
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            CALL coulomb_whittaker(eta, 0, SQRT(2.0d0*e_test)*xext(i), w, wd, sf)
            WRITE(iounit, '(1E20.12)', advance='no') w * Ascale * 10.0d0 ** (sf-sfscale)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/gdstateext.txt')
        
        
        DEALLOCATE(realscatstate)
        
        CALL CONSOLE('Testing regular & irregular solutions for negative energies...')
        
        ALLOCATE(k2(Nfull))
        DO i = 1, Nfull
            k2(i) = 2*m/hbar**2 * (-e_test - V_ptr(xext(i)))
        END DO
        
        ALLOCATE(psi_R(Nfull), psi_I(Nfull))
        h = dx
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        DO 
            IF (l>Nfull-3) EXIT
            CALL coulomb_whittaker(eta, 0, k*xext(Nfull-l), w, wd, sf) 
            psi_I(Nfull-l) = w 
            CALL coulomb_whittaker(eta, 0, k*xext(Nfull-l-1), w, wd, sf)
            psi_I(Nfull-l-1) = w
            IF (ABS(psi_I(Nfull-l)) > 0.0d0) EXIT 
            l = l + 2
        END DO
        PRINT*, l
        
        DO i = 2, Nfull-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        DO i = Nfull-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
        wronski = wronskian_16_real(REAL(psi_R, KIND=16),REAL(psi_I, KIND=16),dx)
        psi_R = psi_R * 2*m/hbar**2 / wronski
        psi_I = psi_I * 2*m/hbar**2 / wronski
        
        
                
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/psifull.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, Nfull
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') psi_R(i)
            WRITE(iounit, '(1E20.12)', advance='no') psi_I(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/psufull.txt') 
        
        DEALLOCATE(psi_R,psi_I)
        
        ALLOCATE(psi_R(mp), psi_I(mp))
        h = dx
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        DO 
            IF (l>mp-3) EXIT
            CALL coulomb_whittaker(eta, 0, k*xext(mp-l), w, wd, sf) 
            psi_I(mp-l) = w 
            CALL coulomb_whittaker(eta, 0, k*xext(mp-l-1), w, wd, sf)
            psi_I(mp-l-1) = w 
            IF (ABS(psi_I(mp-l)) > 0.0d0) EXIT 
            l = l + 2
        END DO
        PRINT*, l
        
        DO i = 2, mp-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        DO i = mp-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
        wronski = wronskian_16_real(REAL(psi_R, KIND=16),REAL(psi_I, KIND=16),dx)
        psi_R = psi_R * 2*m/hbar**2 / wronski
        psi_I = psi_I * 2*m/hbar**2 / wronski
                
        OPEN(NEWUNIT=iounit, FILE='DATA/Test/psishort.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, mp
            WRITE(iounit, '(1E20.12)', advance='no') xext(i)
            WRITE(iounit, '(1E20.12)', advance='no') psi_R(i) 
            WRITE(iounit, '(1E20.12)', advance='no') psi_I(i) 
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Test/psishort.txt') 
        
        DEALLOCATE(psi_R,psi_I)
        
        
        
        DEALLOCATE(xext, dstate, k2)
        
    END IF
    
    
    !=================================================================================
    !                     Bound states in Coulomb-like potential
    !=================================================================================
    IF (boundstates) THEN
        INQUIRE(DIRECTORY="DATA/BoundStates", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/BoundStates")
            CALL CONSOLE('Subdirectory "BoundStates" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "BoundStates" already exists.')
        END IF
        
        INQUIRE(DIRECTORY="DATA/BoundStates/Eigenfunctions", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/BoundStates/Eigenfunctions")
            CALL CONSOLE('Subdirectory "Eigenfunctions" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "Eigenfunctions" already exists.')
        END IF
        
        V_ptr => V_2D
        
        ! Printing V(x) to file
        OPEN(NEWUNIT=iounit, FILE='DATA/BoundStates/V.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/BoundStates/V.txt')
        
        
        
        CALL rydberg_grid(Nbound, NperN, gridtail, Ebound_min, boundEgrid)
        Negrid = SIZE(boundEgrid)
        
        IF (ALLOCATED(eigEFull)) DEALLOCATE(eigEFull)
        IF (ALLOCATED(defects)) DEALLOCATE(defects)
        IF (ALLOCATED(logdevsFull)) DEALLOCATE(logdevsFull)
        
        ALLOCATE(eigEFull(nv,Nbound), defects(nv,Nbound), logdevsFull(nv,Negrid))
        
        DO j = 1, nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing in short-range potential with V_A = ', V_A, ' ...'
            CALL CONSOLE(message)
            
            IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
            IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
            IF (ALLOCATED(logdevs)) DEALLOCATE(logdevs)
            IF (ALLOCATED(defect)) DEALLOCATE(defect)
            
            CALL ComputeBoundStates(x, boundEgrid, m, V_ptr, Z, Nbound, Nstart, eigE, eigFunc, Nfound, logdevs, defect)    
            
            eigEFull(j,:) = eigE
            defects(j,:) = defect
            logdevsFull(j,:) = logdevs
            
            ! Printing eigenfunctions
            WRITE(filename, '(F0.3)') V_A
            filename = 'DATA/BoundStates/Eigenfunctions/eigFunc_VA_' // TRIM(ADJUSTL(filename)) // '.txt'
            OPEN(NEWUNIT=iounit, FILE=filename, STATUS='replace', ACTION='write')
            WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'R [au]','|psi_n>'
            WRITE(iounit,*)
            DO jj = 1, SIZE(x)
                WRITE(iounit, '(1E20.12)', advance='no') x(jj)
                DO i = 1, MIN(Nfound,Nprint)
                    WRITE(iounit, '(1E20.12)', advance='no') eigFunc(jj,i)
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            
            IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
            IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
            IF (ALLOCATED(logdevs)) DEALLOCATE(logdevs)
            IF (ALLOCATED(defect)) DEALLOCATE(defect)
            
        END DO

        
        
        
        ! Printing log-derivatives on test grid
        OPEN(NEWUNIT=iounit, FILE='DATA/BoundStates/logderivs.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/BoundStates/logderivs.txt')
        
        ! Printing eigenenergies
        OPEN(NEWUNIT=iounit, FILE='DATA/BoundStates/eigE.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nfound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') eigEFull(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/BoundStates/eigE.txt')
        
        
        ! Printing defects
        OPEN(NEWUNIT=iounit, FILE='DATA/BoundStates/defects.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nfound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') defects(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/BoundStates/defects.txt')        
        
        
        IF (ALLOCATED(eigEFull)) DEALLOCATE(eigEFull)
        IF (ALLOCATED(defects)) DEALLOCATE(defects)
        IF (ALLOCATED(logdevsFull)) DEALLOCATE(logdevsFull)
        
        
        
    END IF
    
    
    !=================================================================================
    !                PHP Bound states in Coulomb-like potential
    !=================================================================================
    IF (PHPstates) THEN
        INQUIRE(DIRECTORY="DATA/PHPStates", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/PHPStates")
            CALL CONSOLE('Subdirectory "PHPStates" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "PHPStates" already exists.')
        END IF
        
        INQUIRE(DIRECTORY="DATA/PHPStates/Eigenfunctions", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/PHPStates/Eigenfunctions")
            CALL CONSOLE('Subdirectory "Eigenfunctions" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "Eigenfunctions" already exists.')
        END IF
        
        V_ptr => V_2D
        dstateptr => dstate2D
        
        ! Creating initial discrete state
        IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
        ALLOCATE(dstate(SIZE(x)))
        DO i = 1, SIZE(x)
            dstate(i) = dstateptr(x(i))
        END DO
        norm = SUM( ABS(dstate)**2 ) * dx
        dstate = dstate / SQRT(norm)
        CALL CONSOLE('Initial discrete state created successfully.')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            WRITE(iounit, '(1E20.12)', advance='no') dstate(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/PHPStates/dstate.txt')
        
        ! Printing V(x) to file
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/V.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/PHPStates/V.txt')
        
        
        
        CALL rydberg_grid(Nstart, NperN, gridtail, Ebound_min, boundEgrid)
        Negrid = SIZE(boundEgrid)
        
        IF (ALLOCATED(eigEFull)) DEALLOCATE(eigEFull)
        IF (ALLOCATED(overlapsFull)) DEALLOCATE(overlapsFull)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        IF (ALLOCATED(defects)) DEALLOCATE(defects)
        
        ALLOCATE(eigEFull(nv,Nbound), overlapsFull(nv,Negrid), VdnFull(nv,Nbound), defects(nv,Nbound))
        
        DO j = 1, nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing in short-range potential with V_A = ', V_A, ' ...'
            CALL CONSOLE(message)
            
            IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
            IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
            IF (ALLOCATED(overlaps)) DEALLOCATE(overlaps)
            IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
            IF (ALLOCATED(defect)) DEALLOCATE(defect)
            
            CALL ComputePHPBoundStates(x, dstate, boundEgrid, m, V_ptr, Z, Nbound, Nstart, eigE, eigFunc, Nfound, overlaps, defect)
            
            eigEFull(j,:) = eigE
            overlapsFull(j,:) = overlaps
            defects(j,:) = defect
            
            IF (ALLOCATED(defect)) DEALLOCATE(defect)
            
            ! Printing eigenfunctions
            WRITE(filename, '(F0.3)') V_A
            filename = 'DATA/PHPStates/Eigenfunctions/eigFunc_VA_' // TRIM(ADJUSTL(filename)) // '.txt'
            OPEN(NEWUNIT=iounit, FILE=filename, STATUS='replace', ACTION='write')
            WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'R [au]','|psi_n>'
            WRITE(iounit,*)
            DO jj = 1, SIZE(x)
                WRITE(iounit, '(1E20.12)', advance='no') x(jj)
                DO i = 1, MIN(Nfound,Nprint)
                    WRITE(iounit, '(1E20.12)', advance='no') eigFunc(jj,i)
                END DO
                WRITE(iounit,*)
            END DO
            CLOSE(iounit)
            
            IF (ALLOCATED(overlaps)) DEALLOCATE(overlaps)
            IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
            
            ALLOCATE(Vdn(Nbound))
            CALL ComputeVdn(x, dstate, m, V_ptr, eigFunc, Vdn)
            VdnFull(j,:) = Vdn
    
            IF (ALLOCATED(eigFunc)) DEALLOCATE(eigFunc)
            IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
            
        END DO

        
        
        
        ! Printing log-derivatives on test grid
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/overlaps.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','d/dr log(u) [1/Bohr]'
        WRITE(iounit,*)
        DO i = 1, Negrid
            WRITE(iounit, '(I15)', advance='no') i
            WRITE(iounit, '(1E20.12)', advance='no') boundEgrid(i) * phys_h0
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') overlapsFull(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/PHPStates/overlaps.txt')
        
        ! Printing eigenenergies
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/eigE.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') eigEFull(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/PHPStates/eigE.txt')
        
        ! Printing Vdn
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/Vdn.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n', 'Vdn [sqrt(eV)]'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') VdnFull(j,i) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/PHPStates/Vdn.txt')
        
        ! Printing defects
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/defects.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nfound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') defects(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/PHPStates/defects.txt')  
        
        
        ! Saving Vdn and eigE to bin files for later use
        CALL CONSOLE('Caching eigEFull and VdnFull to binary files...')
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/eigE.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) eigEFull
        CLOSE(iounit)
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/Vdn.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) VdnFull
        CLOSE(iounit)
        CALL CONSOLE('Binary cache saved successfully.')
               
        
        
        IF (ALLOCATED(eigEFull)) DEALLOCATE(eigEFull)
        IF (ALLOCATED(logdevsFull)) DEALLOCATE(logdevsFull)
        IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        IF (ALLOCATED(defects)) DEALLOCATE(defects)
        
        
        
    END IF
    
    
    
    !=================================================================================
    !            Discrete state computation for Short Range Potential
    !=================================================================================
    IF (NO_coulomb_computation) THEN
        INQUIRE(DIRECTORY="DATA/ShortRange", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/ShortRange")
            CALL CONSOLE('Subdirectory "ShortRange" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "ShortRange" already exists.')
        END IF
        
        ! Assignment of pointers to corresponding functions
        V_ptr => V
        dstateptr => dstate1

       ! Printing V(x) to file
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/V.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/ShortRange/V.txt')

    
        ! Creating initial discrete state
        IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
        ALLOCATE(dstate(SIZE(x)))
        DO i = 1, SIZE(x)
            dstate(i) = dstateptr(x(i))
        END DO
        norm = SUM( ABS(dstate)**2 ) * dx
        dstate = dstate / SQRT(norm)
        CALL CONSOLE('Initial discrete state created successfully.')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            WRITE(iounit, '(1E20.12)', advance='no') dstate(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/ShortRange/dstate.txt')
    
    
        ! Computation of Gamma & Delta
        IF (ALLOCATED(DeltaFull)) DEALLOCATE(DeltaFull)
        IF (ALLOCATED(Gamma2Full)) DEALLOCATE(Gamma2Full)
        IF (ALLOCATED(PhaseshiftFull)) DEALLOCATE(PhaseshiftFull)
        IF (ALLOCATED(DS_PhaseshiftFull)) DEALLOCATE(DS_PhaseshiftFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        ALLOCATE(DeltaFull(nv,ep), Gamma2Full(nv,ep), PhaseshiftFull(nv,ep), DS_PhaseshiftFull(nv,ep), VdeFull(nv,ep))
    
        DO j = 1, nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing in short-range potential with V_A = ', V_A, ' ...'
            CALL CONSOLE(message)
            CALL compute_DS(x, e, m, V_ptr, dstate, Delta, Gamma2, Phaseshift, DS_phaseshift, Vde)
            DeltaFull(j,:) = Delta
            Gamma2Full(j,:) = Gamma2
            IF (unwrap) THEN
                CALL unwrap_phaseshift(Phaseshift)
                CALL unwrap_phaseshift(DS_phaseshift)
            END IF
            PhaseshiftFull(j,:) = Phaseshift
            DS_PhaseshiftFull(j,:) = DS_phaseshift
            VdeFull(j,:) = Vde
            IF (ALLOCATED(Delta)) DEALLOCATE(Delta)
            IF (ALLOCATED(Gamma2)) DEALLOCATE(Gamma2)
            IF (ALLOCATED(Phaseshift)) DEALLOCATE(Phaseshift)
            IF (ALLOCATED(DS_phaseshift)) DEALLOCATE(DS_phaseshift)
            IF (ALLOCATED(Vde)) DEALLOCATE(Vde)
        END DO
    
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/delta.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/ShortRange/delta.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/gamma.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/ShortRange/gamma.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/phaseshift.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/ShortRange/phaseshift.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/DS_phaseshift.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/ShortRange/DS_phaseshift.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/ShortRange/Vde.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Vde [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) >= 0.0d0) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') VdeFull(j,i) * SQRT(phys_h0)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/ShortRange/Vde.txt')
    
    END IF
    
    
    !=================================================================================
    !            Discrete state computation for Coulomb-like Potential
    !=================================================================================
    IF (coulomb_computation) THEN
        INQUIRE(DIRECTORY="DATA/Coulomb", EXIST=status)
        IF (.NOT. status) THEN
            CALL SYSTEM("mkdir -p DATA/Coulomb")
            CALL CONSOLE('Subdirectory "Coulomb" created successfully.')
        ELSE
            CALL CONSOLE('Subdirectory "Coulomb" already exists.')
        END IF
        
        ! Assignment of pointers to corresponding functions
        V_ptr => V_coulomb2 !V_2D
        dstateptr => dstate1 !dstate2D

       ! Printing V(x) to file
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/V.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/Coulomb/V.txt')
        
        ! Printing V-Coulomb to file
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/V_NoC.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','V(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            DO j = 1, nv
                V_A = V_params(j)
                WRITE(iounit, '(1E20.12)', advance='no') (V_ptr(x(i))+Z/(x(i)+Creg)) * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Coulomb/V_NoC.txt')

    
        ! Creating initial discrete state
        ALLOCATE(dstate(SIZE(x)))
        DO i = 1, SIZE(x)
            dstate(i) = dstateptr(x(i))
        END DO
        norm = SUM( ABS(dstate)**2 ) * dx
        dstate = dstate / SQRT(norm)
        CALL CONSOLE('Initial discrete state created successfully.')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'R [au]','dstate [a.u.]'
        WRITE(iounit,*)
        DO i = 1, SIZE(x)
            WRITE(iounit, '(1E20.12)', advance='no') x(i)
            WRITE(iounit, '(1E20.12)', advance='no') dstate(i)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Coulomb/dstate.txt')
    
    
        ! Computation of Gamma & Delta
        IF (ALLOCATED(DeltaFull)) DEALLOCATE(DeltaFull)
        IF (ALLOCATED(Gamma2Full)) DEALLOCATE(Gamma2Full)
        IF (ALLOCATED(PhaseshiftFull)) DEALLOCATE(PhaseshiftFull)
        IF (ALLOCATED(DS_PhaseshiftFull)) DEALLOCATE(DS_PhaseshiftFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        IF (ALLOCATED(DSenergies)) DEALLOCATE(DSenergies)
        ALLOCATE(DeltaFull(nv,ep), Gamma2Full(nv,ep), PhaseshiftFull(nv,ep), DS_PhaseshiftFull(nv,ep), VdeFull(nv,ep), DSenergies(nv))
    
        DO j = 1, nv
            V_A = V_params(j)
            WRITE(message, '(A,F0.3,A)') 'Computing in Coulomb-like potential with V_A = ', V_A, ' ...'
            CALL CONSOLE(message)
            CALL compute_DSenergy(x, m, V_ptr, dstate, DSenergies(j))
            CALL compute_coulomb_DS(x, e, m, V_ptr, Z, dstate, Delta, Gamma2, Phaseshift, DS_phaseshift, Vde)
            DeltaFull(j,:) = Delta
            Gamma2Full(j,:) = Gamma2
            IF (unwrap) THEN
                CALL unwrap_phaseshift(Phaseshift)
                CALL unwrap_phaseshift(DS_phaseshift)
            END IF
            PhaseshiftFull(j,:) = Phaseshift
            DS_PhaseshiftFull(j,:) = DS_phaseshift
            VdeFull(j,:) = Vde
            IF (ALLOCATED(Delta)) DEALLOCATE(Delta)
            IF (ALLOCATED(Gamma2)) DEALLOCATE(Gamma2)
            IF (ALLOCATED(Phaseshift)) DEALLOCATE(Phaseshift)
            IF (ALLOCATED(DS_phaseshift)) DEALLOCATE(DS_phaseshift)
            IF (ALLOCATED(Vde)) DEALLOCATE(Vde)
        END DO
        
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/DSenergies.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A3, A21)', advance='no') '#', 'V_A','E_DS [eV]'
        WRITE(iounit,*)
        DO j = 1, nv
            WRITE(iounit, '(F5.3, 1E20.12)') V_params(j), DSenergies(j) * phys_h0
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Coulomb/DSenergies.txt')
        
    
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/delta.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/Coulomb/delta.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/gamma.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/Coulomb/gamma.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/phaseshift.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/Coulomb/phaseshift.txt')

        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/defect.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','defect_phase'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) < 0.0d0) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') PhaseshiftFull(j,i)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Coulomb/phaseshift.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/DS_phaseshift.txt', STATUS='replace', ACTION='write')
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
        CALL CONSOLE('Data successfully written to DATA/Coulomb/DS_phaseshift.txt')

        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/DS_defect.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','DS_defect_phase'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) < 0.0d0) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') DS_phaseshiftFull(j,i)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Coulomb/DS_defect.txt')
    
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/Vde.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)', advance='no') '#', 'E [eV]','Vde [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(e)
            IF (e(i) >= 0.0d0) THEN
                WRITE(iounit, '(1E20.12)', advance='no') e(i) * phys_h0
                DO j = 1, nv
                    WRITE(iounit, '(1E20.12)', advance='no') VdeFull(j,i) * SQRT(phys_h0)
                END DO
                WRITE(iounit,*)
            END IF
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/Coulomb/Vde.txt')
        
        ! Saving Vde to bin files for later use
        CALL CONSOLE('Caching VdeFull to binary files...')
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/Vde.bin', &
             FORM='unformatted', STATUS='replace', ACTION='write')
        WRITE(iounit) VdeFull
        CLOSE(iounit)
        CALL CONSOLE('Binary cache saved successfully.')
        
        
        IF (ALLOCATED(DeltaFull)) DEALLOCATE(DeltaFull)
        IF (ALLOCATED(Gamma2Full)) DEALLOCATE(Gamma2Full)
        IF (ALLOCATED(PhaseshiftFull)) DEALLOCATE(PhaseshiftFull)
        IF (ALLOCATED(DS_PhaseshiftFull)) DEALLOCATE(DS_PhaseshiftFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        IF (ALLOCATED(dstate)) DEALLOCATE(dstate)
    
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
        
        IF (ALLOCATED(eigEFull)) DEALLOCATE(eigEFull)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        IF (ALLOCATED(DeltaContinuousFull)) DEALLOCATE(DeltaContinuousFull)
        IF (ALLOCATED(DeltaRydbergFull)) DEALLOCATE(DeltaRydbergFull)
        
        CALL CONSOLE('Loading eigEFull, VdnFull & VdeFull from binary files...')
        
        ALLOCATE(eigEFull(nv,Nbound))
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/eigE.bin', FORM='unformatted', STATUS='old', ACTION='read', IOSTAT=istat)
        IF (istat /= 0) THEN
            CALL CONSOLE('[ERROR]: File eigE.bin cannot be found or openned.')
            STOP 'Error reading eigE.bin'
        END IF
        READ(iounit) eigEFull
        CLOSE(iounit)

        ALLOCATE(VdnFull(nv,Nbound))
        OPEN(NEWUNIT=iounit, FILE='DATA/PHPStates/Vdn.bin', &
             FORM='unformatted', STATUS='old', ACTION='read', IOSTAT=istat)
        IF (istat /= 0) THEN
            CALL CONSOLE('[ERROR]: File Vdn.bin cannot be found or openned.')
            STOP 'Error reading Vdn.bin'
        END IF
        READ(iounit) VdnFull
        CLOSE(iounit)
        
        ALLOCATE(VdeFull(nv,SIZE(e)))
        OPEN(NEWUNIT=iounit, FILE='DATA/Coulomb/Vde.bin', &
             FORM='unformatted', STATUS='old', ACTION='read', IOSTAT=istat)
        IF (istat /= 0) THEN
            CALL CONSOLE('[ERROR]: File Vde.bin cannot be found or openned.')
            STOP 'Error reading Vde.bin'
        END IF
        READ(iounit) VdeFull
        CLOSE(iounit)

        CALL CONSOLE('Binary cache loaded successfully.')
        
        
        ! Printing eigenenergies
        OPEN(NEWUNIT=iounit, FILE='DATA/Hilbert/En.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A21)', advance='no') '#', 'n','E_n [eV]'
        WRITE(iounit,*)
        DO i = 1, Nbound
            WRITE(iounit, '(I3)', advance='no') i
            DO j = 1, nv 
                WRITE(iounit, '(1E20.12)', advance='no') eigEFull(j,i) * phys_h0
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
            WRITE(message, '(A,F0.3,A)') 'Computing Rydberg-part of Delta for V_A = ', V_A, '...'
            CALL CONSOLE(message)
            IF (ALLOCATED(DeltaRydberg)) DEALLOCATE(DeltaRydberg)
            IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
            IF (ALLOCATED(eigE)) DEALLOCATE(eigE)
            ALLOCATE(DeltaRydberg(SIZE(e)))
            ALLOCATE(Vdn(Nbound))
            ALLOCATE(eigE(Nbound))
            Vdn(:) = VdnFull(j,:)
            eigE(:) = eigEFull(j,:)
            CALL ComputeDeltaRydberg(e, Vdn, eigE, DeltaRydberg)
            DeltaRydbergFull(j,:) = DeltaRydberg
            IF (ALLOCATED(DeltaRydberg)) DEALLOCATE(DeltaRydberg)
            IF (ALLOCATED(Vdn)) DEALLOCATE(Vdn)
            IF (ALLOCATED(eigE)) DEALLOCATE(eigE)      
            WRITE(message, '(A,F0.3,A)') 'Computation of Rydberg-part of Delta for V_A = ', V_A, ' completed.'
            CALL CONSOLE(message)
        END DO
        
        
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
        
        
        
        IF (ALLOCATED(eigEFull)) DEALLOCATE(eigEFull)
        IF (ALLOCATED(VdnFull)) DEALLOCATE(VdnFull)
        IF (ALLOCATED(VdeFull)) DEALLOCATE(VdeFull)
        IF (ALLOCATED(DeltaContinuousFull)) DEALLOCATE(DeltaContinuousFull)
        IF (ALLOCATED(DeltaRydbergFull)) DEALLOCATE(DeltaRydbergFull)
        
        
    END IF
    
    
    
    

    ! Release allocated resources
    IF (ALLOCATED(x)) DEALLOCATE(x)
    IF (ALLOCATED(e)) DEALLOCATE(e)
    IF (ALLOCATED(V_params)) DEALLOCATE(V_params)
    
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
        V_coulomb2 = ( 0.5d0 * R**2 - 1.5d0 + V_A ) * sR - (1.0 - sR) * Z / (R + Creg)
    END FUNCTION V_coulomb2


    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_2D(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk), PARAMETER :: a = 0.5d0
        REAL(KIND = idk), PARAMETER :: b = 1.2d0
        REAL(KIND = idk), PARAMETER :: c = 1.8d0
        REAL(KIND = idk), PARAMETER :: d = 0.6519d0
        REAL(KIND = idk), PARAMETER :: rc = 4.5d0
        REAL(KIND = idk), PARAMETER :: rrc = 1.5d0
        V_2D = - 1.0d0 / R + 1.0d0 / R**2 + a * EXP(-(R-rc)**2/b**2) - d * EXP(-R**2/4.0d0) * TANH( (V_A-rrc)/c )
    END FUNCTION V_2D
    
    
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


    ! Dstate 2D model
    REAL(KIND = idk) FUNCTION dstate2D(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        dstate2D = R * EXP(-(R-1.825d0)**2 / 2.0d0 * 0.65d0)
    END FUNCTION dstate2D
    
END PROGRAM MAIN
  