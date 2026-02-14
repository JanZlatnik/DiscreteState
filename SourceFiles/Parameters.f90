
!------------------------< Parameters for Photodetachment >------------------------!
!                                                                                  !    
! Contains: Basic mathematical and physical constants, parameters, mesh creation   !                                                                             
!                                                                                  !
! Last revision:    13/08/2025                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!
    
MODULE Parameters
    IMPLICIT NONE
    INTEGER, PARAMETER :: idk = 8
    INTEGER, PARAMETER :: qdk = 16
    
    ! Mathematical constants:
    REAL(KIND = idk), PARAMETER :: pi = 3.141592653589793238462643383279502d0
    REAL(KIND = idk), PARAMETER :: euler = 2.718281828459045235d0
    
    ! Physical constants:
    REAL(KIND = idk), PARAMETER :: &
        phys_a0 = 0.52917706d0, &   ! Bohr radius [angstroms]
        phys_h0 = 27.2113961d0, &   ! 1 Hartree in [eV]
        phys_k = 8.617385d-5, &     ! Boltzman konstant in [eV/K]
        phys_me = 9.1093897d-31,&   ! electron mass [kg]
        phys_u = 1.6605402d-27, &   ! (unified) atomic mass unit [kg]
        hbar = 1.0d0, &             ! Planck constant [au]
        cmtoht = 4.556d-6           ! cm^-1 to hartree
    
    
    ! Computational settings
    LOGICAL, PARAMETER                  :: Rydbergstates = .TRUE.
    LOGICAL, PARAMETER                  :: PHPRydbergstates = .FALSE.
    LOGICAL, PARAMETER                  :: dstate_computation = .FALSE.
    LOGICAL, PARAMETER                  :: hilbert_computation = .FALSE.
    
        
    LOGICAL, PARAMETER                  :: unwrap = .TRUE.
    
    
    
    ! Potential settings
    REAL(KIND = idk), PARAMETER         :: Vmin = 0.7d0                 ! potential parameter V_A - minimum
    REAL(KIND = idk), PARAMETER         :: Vmax = 10.0d0                ! potential parameter V_A - maximum
    INTEGER, PARAMETER                  :: nv = 10                      ! number of diferent potential parameters calculated
    REAL(KIND = idk), PARAMETER         :: Z = 1.0d0                    ! strenght of Coulombic potential, i.e. Z/r

    ! Bound state parameters
    INTEGER, PARAMETER          :: max_iter = 60                        ! maximum number of bisections on a test grid
    INTEGER, PARAMETER          :: max_iter2 = 50                       ! maximum number of bisections in defect convergence process
    REAL(KIND = idk), PARAMETER :: mueps = 1.0d-6                       ! threshold for defect convergence before fitting muint + B/n^2 
    INTEGER, PARAMETER          :: N_fit_min_n = 100                    ! threshold for minimal number of states computed before fitting
    INTEGER, PARAMETER          :: N_fit_points = 15                    ! number of fitted points
    INTEGER, PARAMETER          :: NperN = 400                          ! number of test grid points per expected state
    REAL(KIND = idk), PARAMETER :: gridtail = 1.0d0                     ! tail size for test grid
    REAL(KIND = idk), PARAMETER :: Ebound_min = -0.18d0                  ! minimum energy in test grid
    INTEGER, PARAMETER          :: Nbound = 200                         ! number of bound states explicitely calculated
    INTEGER, PARAMETER          :: Nstart = 10                          ! number of bound states computed on test grid
    INTEGER, PARAMETER          :: Nprint = 50                          ! number of wavefucntions of bound states printed
    
    ! Mesh settings
    REAL(KIND = idk), PARAMETER :: xmin = 1.0d-10         ! in [au]
    REAL(KIND = idk), PARAMETER :: xmax = 14.0d0          ! in [au]
    INTEGER, PARAMETER          :: mp = 3000
    
    ! Energy mesh settings
    REAL(KIND = idk), PARAMETER :: Emin = -4000.0d0/phys_h0       ! [eV] to [au]
    REAL(KIND = idk), PARAMETER :: Emax = 4000.0d0/phys_h0        ! [eV] to [au]
    INTEGER, PARAMETER          :: ep = 10000000

    ! Other parameters
    REAL(KIND = idk), PARAMETER :: m = 1.0d0        ! mass in [au]
    INTEGER, PARAMETER          :: l_ang = 1        ! angular momentum quantum number of the detached electron (l=0 for s-wave detachment, l=1 for p-wave detachment, etc.)
    REAL(KIND = idk), PARAMETER :: R0 = 2.0d0       ! position of wronskian evaluation [au]
    
    CONTAINS
    
    SUBROUTINE make_mesh(xmin, xmax, mp, x, dx)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(OUT) :: dx
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: x(:)
        INTEGER, INTENT(IN) :: mp
        REAL(KIND = idk), INTENT(IN) :: xmin, xmax 
        INTEGER :: i
        
        IF (ALLOCATED(x)) DEALLOCATE(x)
        ALLOCATE(x(mp))
        dx = (xmax - xmin) / (mp-1.0d0)
        DO i = 1,mp
            x(i) = xmin + (i-1)*dx
        END DO
        
    END SUBROUTINE
    
    
    SUBROUTINE CONSOLE(message)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: message
        CHARACTER(LEN=10) :: timestr
        INTEGER :: h, m, s
        
        CALL DATE_AND_TIME(TIME=timestr)

        READ(timestr(1:2),*) h
        READ(timestr(3:4),*) m
        READ(timestr(5:6),*) s

        WRITE(*,'(A,I2.2,A,I2.2,A,I2.2,A,1X,A)') '[', h, ':', m, ':', s, ']:', TRIM(message)
    END SUBROUTINE CONSOLE
    
    FUNCTION STR(k) RESULT(str_out)
        INTEGER, INTENT(IN) :: k
        CHARACTER(LEN=20)   :: str_out
        WRITE(str_out, '(I0)') k 
        str_out = ADJUSTL(str_out)    
    END FUNCTION STR


    FUNCTION cutoff_energy(l) RESULT(Ecut)
        INTEGER, INTENT(IN) :: l
        REAL(KIND = idk) :: Ecut
        IF (l == 0) THEN
            Ecut = -1.0d0
        ELSE
            Ecut = -1.0d0 / (2.0d0 * REAL(l, KIND=idk)**2)
        END IF
    END FUNCTION cutoff_energy
    
    
    
    
    SUBROUTINE parameters_write(unit)
        INTEGER, INTENT(IN)  :: unit
        
        WRITE(unit,*) "=============================================================="
        WRITE(unit,*) "                  Compuational Settings                       "
        WRITE(unit,*) "=============================================================="
        
        WRITE(unit,*) ""
        WRITE(unit,'(A35,1X,L1)') "RydbergStates:", Rydbergstates 
        WRITE(unit,'(A35,1X,L1)') "PHPRydbergStates:", PHPRydbergstates
        WRITE(unit,'(A35,1X,L1)') "dstate_computation:", dstate_computation
        WRITE(unit,'(A35,1X,L1)') "hilbert_computation:", hilbert_computation
        WRITE(unit,*) ""

        
        WRITE(unit,*) "=============================================================="
        WRITE(unit,*) "                 Computation Parameters                       "
        WRITE(unit,*) "=============================================================="
        
        WRITE(unit,*) ""
        WRITE(unit,*) "[x-Grid]"
        WRITE(unit,'(A25,1X,E0.6)') "xmin:", xmin
        WRITE(unit,'(A25,1X,F0.6)') "xmax:", xmax
        WRITE(unit,'(A25,1X,I0)')   "mp:", mp 
        WRITE(unit,*) ""
        
        WRITE(unit,*) ""
        WRITE(unit,*) "[E-Grid]"
        WRITE(unit,'(A25,1X,F0.6)') "Emin:", Emin * phys_h0
        WRITE(unit,'(A25,1X,F0.6)') "Emax:", Emax * phys_h0
        WRITE(unit,'(A25,1X,I0)')   "ep:", ep 
        WRITE(unit,*) ""     
        
        WRITE(unit,*) ""
        WRITE(unit,*) "[Potential parameters]"
        WRITE(unit,'(A25,1X,F0.6)') "Vmin:", Vmin
        WRITE(unit,'(A25,1X,F0.6)') "Vmax:", Vmax
        WRITE(unit,'(A25,1X,F0.6)') "Z:", Z
        WRITE(unit,'(A25,1X,I0)') "nv:", nv
        WRITE(unit,*) ""   

        WRITE(unit,*) ""
        WRITE(unit,*) "[Bound states settings]"
        WRITE(unit,'(A25,1X,I0)') "max_iter:", max_iter
        WRITE(unit,'(A25,1X,I0)') "max_iter2:", max_iter2
        WRITE(unit,'(A25,1X,E0.6)') "mueps:", mueps
        WRITE(unit,'(A25,1X,I0)') "N_fit_min_n:", N_fit_min_n
        WRITE(unit,'(A25,1X,I0)') "N_fit_points:", N_fit_points
        WRITE(unit,'(A25,1X,I0)') "NperN:", NperN
        WRITE(unit,'(A25,1X,F0.6)') "gridtail:", gridtail
        WRITE(unit,'(A25,1X,F0.6)') "Ebound_min:", Ebound_min * phys_h0
        WRITE(unit,'(A25,1X,I0)') "Nbound:", Nbound
        WRITE(unit,'(A25,1X,I0)') "Nstart:", Nstart
        WRITE(unit,'(A25,1X,I0)') "Nprint:", Nprint
        WRITE(unit,*) ""   

        WRITE(unit,*) ""
        WRITE(unit,*) "[Other Parameters]"
        WRITE(unit,'(A25,1X,F0.6)') "m:", m
        WRITE(unit,'(A25,1X,I0)')   "l:", l_ang
        WRITE(unit,'(A25,1X,F0.6)') "R0:", R0
        WRITE(unit,*) ""   
        
        
        WRITE(unit,*) "=============================================================="
    
    END SUBROUTINE parameters_write
    
    
    SUBROUTINE parameters_o()
        INTEGER :: param_u
        INTEGER :: console_u = 6
        
        CALL parameters_write(console_u)
        
        OPEN(newunit = param_u, file = 'parameters.txt', status = 'replace', form = 'formatted')
        CALL parameters_write(param_u)
        CALL CONSOLE('Parameters have been succesfully written to parameters.txt.')
    
    
    END SUBROUTINE



    
    
END MODULE Parameters