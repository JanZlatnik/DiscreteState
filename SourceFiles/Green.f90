
!-------------------------------------< Green >-------------------------------------!
!                                                                                   !    
! Contains: Green function for the equation [E-T_N-V]|psi> = |rightside>            !                                                                             
!                                                                                   !
! Last revision:    25/09/2025                                                      !
!                                                                                   !
!-----------------------------------------------------------------------------------!

    
MODULE Green
    USE Math
    USE Parameters
    USE COULCC_M, only: COULCC
    USE whittaker_w, only: coulomb_whittaker
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: apply_green, apply_green_LCP, apply_green_coulomb, apply_green_coulomb_bound, wronskian_16_real, apply_green_analytic_continuation
    
    INTERFACE apply_green
        MODULE PROCEDURE apply_green_real_vector
        MODULE PROCEDURE apply_green_real_matrix
        MODULE PROCEDURE apply_green_complex_vector
        MODULE PROCEDURE apply_green_complex_matrix
    END INTERFACE apply_green
    
    INTERFACE apply_green_coulomb
        MODULE PROCEDURE apply_green_coulomb_real_matrix
        MODULE PROCEDURE apply_green_coulomb_complex_matrix
        MODULE PROCEDURE apply_green_coulomb_real_vector
    END INTERFACE apply_green_coulomb

    INTERFACE apply_green_coulomb_bound
        MODULE PROCEDURE apply_green_coulomb_bound_real_vector
    END INTERFACE apply_green_coulomb_bound
    
    INTERFACE apply_green_LCP
        MODULE PROCEDURE apply_green_real_vector_LCP
        MODULE PROCEDURE apply_green_real_matrix_LCP
        MODULE PROCEDURE apply_green_complex_vector_LCP
        MODULE PROCEDURE apply_green_complex_matrix_LCP
    END INTERFACE apply_green_LCP

    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_green_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: x
        END FUNCTION potential_green_interface
    END INTERFACE

    CONTAINS
    ! Subroutine to apply the Green's function method to solve the previously mentioned differential equation.
    ! |rightside> can be either a real/complex vector or a rectangular matrix (j,:), then the Green's function is applied on each j-th vector

    SUBROUTINE apply_green_real_vector(x, E, m, ext, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, rightside(:), ext
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:), x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(SIZE(x)), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x_ext(i)))
        END DO

        DO i = 2, SIZE(x)-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO

        IF (ISNAN(ABS(psi_R(SIZE(x)-1)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
        
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
        
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        DEALLOCATE(intR,intI,psi_I,psi_R)
        
    END SUBROUTINE apply_green_real_vector
    
    
    SUBROUTINE apply_green_complex_vector(x, E, m, ext, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, ext
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:), x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x_ext(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO

        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
        
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
        
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        DEALLOCATE(intR,intI,psi_I,psi_R)
        
    END SUBROUTINE apply_green_complex_vector
    
    
    SUBROUTINE apply_green_real_matrix(x, E, m, ext, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, rightside(:,:), ext
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)
        INTEGER :: n, i, j, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:), x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(SIZE(x)), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x_ext(i)))
        END DO

        DO i = 2, SIZE(x)-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO

        IF (ISNAN(ABS(psi_R(SIZE(x)-1)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
        IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
           
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
        
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        DO j = 1, SIZE(rightside,1)
            yR = psi_R(1:n) * rightside(j,:) * 2*m/hbar**2
            yI = psi_I(1:n) * rightside(j,:) * 2*m/hbar**2
            CALL REVERSE(yI)
            CALL primitive(yR, h, intR)
            CALL primitive(yI, h, intI)
            CALL REVERSE(intI)
            psi(j,:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        END DO
        DEALLOCATE(yR,yI,intR,intI)

    
    END SUBROUTINE apply_green_real_matrix
    
    
    SUBROUTINE apply_green_complex_matrix(x, E, m, ext, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, ext
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:,:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)
        INTEGER :: n, i, j, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:), x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x_ext(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
    
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        DO j = 1, SIZE(rightside,1)
            yR = psi_R(1:n) * rightside(j,:) * 2*m/hbar**2
            yI = psi_I(1:n) * rightside(j,:) * 2*m/hbar**2
            CALL REVERSE(yI)
            CALL primitive(yR, h, intR)
            CALL primitive(yI, h, intI)
            CALL REVERSE(intI)

            psi(j,:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        END DO
        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)
    
    END SUBROUTINE apply_green_complex_matrix
    
    SUBROUTINE apply_green_real_vector_LCP(x, E, m, ext, LCP, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, rightside(:), ext
        COMPLEX(KIND=idk), INTENT(IN) :: LCP(:)
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        COMPLEX(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk), ALLOCATABLE :: x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - LCP(i))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO

        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
        
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
        
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        DEALLOCATE(intR,intI,psi_I,psi_R)
        
    END SUBROUTINE apply_green_real_vector_LCP
    
    
    SUBROUTINE apply_green_complex_vector_LCP(x, E, m, ext, LCP, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, ext
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:)
        COMPLEX(KIND=idk), INTENT(IN) :: LCP(:)
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        COMPLEX(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk), ALLOCATABLE :: x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - LCP(i))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO

        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
        
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
        
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)
        DEALLOCATE(yR,yI)

        psi = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        DEALLOCATE(intR,intI,psi_I,psi_R)
        
    END SUBROUTINE apply_green_complex_vector_LCP
    
    
    SUBROUTINE apply_green_real_matrix_LCP(x, E, m, ext, LCP, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, rightside(:,:), ext
        COMPLEX(KIND=idk), INTENT(IN) :: LCP(:)
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)
        INTEGER :: n, i, j, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        COMPLEX(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk), ALLOCATABLE :: x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - LCP(i))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO

        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
           
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
        
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        DO j = 1, SIZE(rightside,1)
            yR = psi_R(1:n) * rightside(j,:) * 2*m/hbar**2
            yI = psi_I(1:n) * rightside(j,:) * 2*m/hbar**2
            CALL REVERSE(yI)
            CALL primitive(yR, h, intR)
            CALL primitive(yI, h, intI)
            CALL REVERSE(intI)
            psi(j,:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        END DO
        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)

    
    END SUBROUTINE apply_green_real_matrix_LCP
    
    
    SUBROUTINE apply_green_complex_matrix_LCP(x, E, m, ext, LCP, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, ext
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:,:)
        COMPLEX(KIND=idk), INTENT(IN) :: LCP(:)
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)
        INTEGER :: n, i, j, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        COMPLEX(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk), ALLOCATABLE :: x_ext(:)

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        IF (ext > 1.0d0) THEN
            n = CEILING(SIZE(x) * ext)
            ALLOCATE(x_ext(n))
            DO i = 1, SIZE(x)
                x_ext(i) = x(i)
            END DO
            DO i = SIZE(x), n - 1
                x_ext(i+1) = x_ext(i) + dx
            END DO
        ELSE
            ALLOCATE(x_ext(n))
            DO i = 1, n
                x_ext(i) = x(i)
            END DO
        END IF

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            psi_I(n) = EXP(DCMPLX(0, k*x_ext(n)))
            psi_I(n-1) = EXP(DCMPLX(0, k*x_ext(n-1)))
        ELSE IF (E<0) THEN
            psi_I(n) = EXP(-k*x_ext(n))
            psi_I(n-1) = EXP(-k*x_ext(n-1))
            DO 
                IF (REAL(psi_I(n-l))>0) EXIT
                l = l + 2
                IF (l>n-3) EXIT
                psi_I(n-l) = EXP(-k*x_ext(n-l))
                psi_I(n-l-1) = EXP(-k*x_ext(n-l-1))
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - LCP(i))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
        END IF
    
        DEALLOCATE(k2)
        DEALLOCATE(x_ext)
        
        n = SIZE(x)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        DO j = 1, SIZE(rightside,1)
            yR = psi_R(1:n) * rightside(j,:) * 2*m/hbar**2
            yI = psi_I(1:n) * rightside(j,:) * 2*m/hbar**2
            CALL REVERSE(yI)
            CALL primitive(yR, h, intR)
            CALL primitive(yI, h, intI)
            CALL REVERSE(intI)

            psi(j,:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        END DO
        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)
    
    END SUBROUTINE apply_green_complex_matrix_LCP



    ! Function to calculate the Wronskian of two complex functions.
    FUNCTION wronskian(f, g, dx)
        IMPLICIT NONE
        COMPLEX(KIND=idk), INTENT(IN) :: f(:), g(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        COMPLEX(KIND=idk) :: wronskian
        INTEGER :: Npul
        COMPLEX(KIND=idk) :: df, dg

        Npul = CEILING(mp * R0/ABS(xmax-xmin))
        df = (-f(Npul+2) + 8*f(Npul+1) - 8*f(Npul-1) + f(Npul-2)) / (12*dx)
        dg = (-g(Npul+2) + 8*g(Npul+1) - 8*g(Npul-1) + g(Npul-2)) / (12*dx)
        wronskian = dg*f(Npul) - df*g(Npul)
        
        IF (wronskian == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Wronskian equal to zero')
        END IF
        IF (ISNAN(ABS(df))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF
        IF (ISNAN(ABS(dg))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF
    
    END FUNCTION wronskian
    
    FUNCTION wronskian_16(f, g, dx)
        IMPLICIT NONE
        COMPLEX(KIND=16), INTENT(IN) :: f(:), g(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        COMPLEX(KIND=16) :: wronskian_16
        INTEGER :: Npul
        COMPLEX(KIND=16) :: df, dg

        Npul = CEILING(mp * R0/ABS(xmax-xmin))
        df = (-f(Npul+2) + 8*f(Npul+1) - 8*f(Npul-1) + f(Npul-2)) / (12*dx)
        dg = (-g(Npul+2) + 8*g(Npul+1) - 8*g(Npul-1) + g(Npul-2)) / (12*dx)
        wronskian_16 = dg*f(Npul) - df*g(Npul)
        
        IF (wronskian_16 == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Wronskian equal to zero')
        END IF
        IF (ISNAN(ABS(df))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF
        IF (ISNAN(ABS(dg))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF
    
    END FUNCTION wronskian_16


    FUNCTION wronskian_16_real(f, g, dx)
        IMPLICIT NONE
        REAL(KIND=16), INTENT(IN) :: f(:), g(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        REAL(KIND=16) :: wronskian_16_real
        INTEGER :: Npul
        REAL(KIND=16) :: df, dg

        Npul = CEILING(mp * R0/ABS(xmax-xmin))
        df = (-f(Npul+2) + 8*f(Npul+1) - 8*f(Npul-1) + f(Npul-2)) / (12*dx)
        dg = (-g(Npul+2) + 8*g(Npul+1) - 8*g(Npul-1) + g(Npul-2)) / (12*dx)
        wronskian_16_real = dg*f(Npul) - df*g(Npul)
        
        IF (wronskian_16_real == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: Wronskian equal to zero')
        END IF
        IF (ISNAN(ABS(df))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF
        IF (ISNAN(ABS(dg))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
        END IF
    
    END FUNCTION wronskian_16_real
    
    
    
    
    
    
    
! === Green in Coulomb ===
    
    SUBROUTINE apply_green_coulomb_real_vector(x, E, m, Z, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        REAL(KIND=idk), INTENT(IN) :: rightside(:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk) :: eta, w, wd
        COMPLEX(KIND=idk) :: XX, ETA1
        COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
        INTEGER :: IFAIL, sf1, sf2

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        eta = -Z*m/(hbar**2*k)
        ETA1 = CMPLX(eta,0.0d0,KIND=idk)
        IFAIL = 0 

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            XX = CMPLX(k*x(n), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n) = GC(1)
            XX = CMPLX(k*x(n-1), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n-1) = GC(1)
        ELSE IF (E<0) THEN
            DO 
                IF (l>n-3) EXIT
                CALL coulomb_whittaker(eta, 0, k*x(n-l), w, wd, sf1)
                psi_I(n-l) = w 
                CALL coulomb_whittaker(eta, 0, k*x(n-l-1), w, wd, sf2)
                psi_I(n-l-1) = w * 10.0d0 ** (sf2-sf1)
                IF (ABS(psi_I(n-l))>0) EXIT
                l = l + 2
            END DO
        END IF
     
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF
         
        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
            PRINT*, E
            PRINT*, psi_I(n-l)
            PRINT*, psi_I(n-l-1)
        END IF
    
        DEALLOCATE(k2)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)
        psi(:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)

        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)
    
    END SUBROUTINE apply_green_coulomb_real_vector
    
    
    
    
    SUBROUTINE apply_green_coulomb_real_matrix(x, E, m, Z, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        REAL(KIND=idk), INTENT(IN) :: rightside(:,:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)
        INTEGER :: n, i, j, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk) :: eta, w, wd
        COMPLEX(KIND=idk) :: XX, ETA1
        COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
        INTEGER :: IFAIL, sf1, sf2

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        eta = -Z*m/(hbar**2*k)
        ETA1 = CMPLX(eta,0.0d0,KIND=idk)
        IFAIL = 0

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            XX = CMPLX(k*x(n), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n) = GC(1)
            XX = CMPLX(k*x(n-1), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n-1) = GC(1)
        ELSE IF (E<0) THEN
            DO 
                IF (l>n-3) EXIT
                CALL coulomb_whittaker(eta, 0, k*x(n-l), w, wd, sf1)
                psi_I(n-l) = w 
                CALL coulomb_whittaker(eta, 0, k*x(n-l-1), w, wd, sf2)
                psi_I(n-l-1) = w * 10.0d0 ** (sf2-sf1)
                IF (ABS(psi_I(n-l))>0) EXIT
                l = l + 2
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF
         
        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
            PRINT*, E
            PRINT*, psi_I(n-l)
            PRINT*, psi_I(n-l-1)
        END IF
    
        DEALLOCATE(k2)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        DO j = 1, SIZE(rightside,1)
            yR = psi_R(1:n) * rightside(j,:) * 2*m/hbar**2
            yI = psi_I(1:n) * rightside(j,:) * 2*m/hbar**2
            CALL REVERSE(yI)
            CALL primitive(yR, h, intR)
            CALL primitive(yI, h, intI)
            CALL REVERSE(intI)

            psi(j,:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        END DO
        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)
    
    END SUBROUTINE apply_green_coulomb_real_matrix
    
    
    
    SUBROUTINE apply_green_coulomb_complex_matrix(x, E, m, Z, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        COMPLEX(KIND=idk), INTENT(IN) :: rightside(:,:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:,:)
        INTEGER :: n, i, j, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk) :: eta, w, wd
        COMPLEX(KIND=idk) :: XX, ETA1
        COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
        INTEGER :: IFAIL, sf1, sf2

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        eta = -Z*m/(hbar**2*k)
        ETA1 = CMPLX(eta,0.0d0,KIND=idk)
        IFAIL = 0

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            XX = CMPLX(k*x(n), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n) = GC(1)
            XX = CMPLX(k*x(n-1), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n-1) = GC(1)
        ELSE IF (E<0) THEN
            DO 
                IF (l>n-3) EXIT
                CALL coulomb_whittaker(eta, 0, k*x(n-l), w, wd, sf1)
                psi_I(n-l) = w 
                CALL coulomb_whittaker(eta, 0, k*x(n-l-1), w, wd, sf2)
                psi_I(n-l-1) = w * 10.0d0 ** (sf2-sf1)
                IF (ABS(psi_I(n-l))>0) EXIT
                l = l + 2
            END DO
        END IF
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF
         
        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
            PRINT*, E
            PRINT*, psi_I(n-l)
            PRINT*, psi_I(n-l-1)
        END IF
    
        DEALLOCATE(k2)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        DO j = 1, SIZE(rightside,1)
            yR = psi_R(1:n) * rightside(j,:) * 2*m/hbar**2
            yI = psi_I(1:n) * rightside(j,:) * 2*m/hbar**2
            CALL REVERSE(yI)
            CALL primitive(yR, h, intR)
            CALL primitive(yI, h, intI)
            CALL REVERSE(intI)

            psi(j,:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)
        END DO
        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)
    
    END SUBROUTINE apply_green_coulomb_complex_matrix



    SUBROUTINE apply_green_coulomb_bound_real_vector(x, E, m, Z, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        REAL(KIND=idk), INTENT(IN) :: rightside(:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        REAL(KIND=idk), INTENT(OUT) :: psi(:)
        
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        REAL(KIND=16) :: wronski
        REAL(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk) :: eta, w, wd
        INTEGER :: IFAIL, sf1, sf2

        IF (E >= 0.0d0) THEN
            CALL CONSOLE('[ERROR] apply_green_coulomb_bound_real_vector called with E >= 0')
            RETURN
        END IF

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        eta = -Z*m/(hbar**2*k)
        IFAIL = 0
        

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        DO 
            IF (l>n-3) EXIT
            CALL coulomb_whittaker(eta, 0, k*x(n-l), w, wd, sf1) 
            psi_I(n-l) = w 
            CALL coulomb_whittaker(eta, 0, k*x(n-l-1), w, wd, sf2)
            psi_I(n-l-1) = w * 10.0d0 ** (sf2-sf1)
            IF (ABS(psi_I(n-l)) > 0.0d0) EXIT 
            l = l + 2
        END DO 
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(psi_R(n))) THEN 
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF (ABS(psi_I(1)) == 0.0d0) THEN
             CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF
         
        IF (ISNAN(psi_I(1))) THEN 
            CALL CONSOLE('[ERROR]: psi_I is NaN')
            PRINT*, E
            PRINT*, psi_I(n-l)
            PRINT*, psi_I(n-l-1)
        END IF
    
        DEALLOCATE(k2)


        wronski = wronskian_16_real(psi_R, psi_I, h) 
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(wronski)) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)

        psi(:) = REAL( (psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk )

        DEALLOCATE(yR,yI,intR,intI, psi_R, psi_I)
    
    END SUBROUTINE apply_green_coulomb_bound_real_vector
    
    
    
    
    
    
    SUBROUTINE apply_green_analytic_continuation(x, E, m, Z, potential_green, rightside, psi)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: x(:), E, m, Z
        REAL(KIND=idk), INTENT(IN) :: rightside(:)
        PROCEDURE(potential_green_interface), POINTER, INTENT(IN) :: potential_green
        COMPLEX(KIND=idk), INTENT(OUT) :: psi(:)
        INTEGER :: n, i, l
        REAL(KIND=idk) :: dx, k, h
        COMPLEX(KIND=16) :: wronski
        COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: psi_R, psi_I, yR, yI, intR, intI
        REAL(KIND=idk), ALLOCATABLE :: k2(:)
        REAL(KIND=idk) :: eta, w, wd
        COMPLEX(KIND=idk) :: XX, ETA1
        COMPLEX(KIND=idk), DIMENSION(1) :: FC, GC, FCP, GCP, SIG
        INTEGER :: IFAIL, sf1, sf2

        n = SIZE(x)
        dx = ABS(x(2) - x(1))
        k = SQRT(2*m*ABS(E))/hbar
        h = dx
        
        eta = -Z*m/(hbar**2*k)
        IFAIL = 0 

        ALLOCATE(psi_R(n), psi_I(n))
        psi_R(1) = 0
        psi_R(2) = h
        l = 0
        IF (E >= 0) THEN
            ETA1 = CMPLX(eta,0.0d0,KIND=idk)
            XX = CMPLX(k*x(n), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n) = GC(1)
            XX = CMPLX(k*x(n-1), 0.0d0, KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 12, 0, IFAIL)
            psi_I(n-1) = GC(1)
        ELSE IF (E<0) THEN
            ETA1 = CMPLX(0,-eta,KIND=idk)
            XX = CMPLX(0.0d0, k*x(n), KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 22, 0, IFAIL)
            psi_I(n) = GC(1)
            XX = CMPLX(0.0d0, k*x(n-1), KIND=idk)
            CALL COULCC(XX, ETA1, (0.0d0,0.0d0), 1, FC, GC, FCP, GCP, SIG, 22, 0, IFAIL)
            psi_I(n-1) = GC(1)
        END IF
     
        
        ALLOCATE(k2(n))
        DO i = 1, n
            k2(i) = 2*m/hbar**2 * (E - potential_green(x(i)))
        END DO

        DO i = 2, n-1
            psi_R(i+1) = (2*(1-5*h**2*k2(i)/12)*psi_R(i) - (1+h**2*k2(i-1)/12)*psi_R(i-1)) / (1+h**2*k2(i+1)/12)
        END DO
        
        IF (ISNAN(ABS(psi_R(n)))) THEN
            CALL CONSOLE('[ERROR]: psi_R is NaN')
        END IF

        DO i = n-l-1, 2, -1
            psi_I(i-1) = (2*(1-5*h**2*k2(i)/12)*psi_I(i) - (1+h**2*k2(i+1)/12)*psi_I(i+1)) / (1+h**2*k2(i-1)/12)
        END DO
        
         IF ((ABS(psi_I(1))) == 0.0d0) THEN
            CALL CONSOLE('[ERROR]: psi_I is 0')
         END IF
         
        IF (ISNAN(ABS(psi_I(1)))) THEN
            CALL CONSOLE('[ERROR]: psi_I is NaN')
            PRINT*, E
            PRINT*, psi_I(n-l)
            PRINT*, psi_I(n-l-1)
        END IF
    
        DEALLOCATE(k2)

        wronski = wronskian_16(psi_R, psi_I, h)
        IF (wronski == 0) THEN
            CALL CONSOLE('[ERROR]: Wronski == 0')
        END IF
        IF (ISNAN(ABS(wronski))) THEN
            CALL CONSOLE('[ERROR]: Wronski is NaN')
        END IF
    
        ALLOCATE(yR(n), yI(n), intR(n), intI(n))
        yR = psi_R(1:n) * rightside * 2*m/hbar**2
        yI = psi_I(1:n) * rightside * 2*m/hbar**2
        CALL REVERSE(yI)
        CALL primitive(yR, h, intR)
        CALL primitive(yI, h, intI)
        CALL REVERSE(intI)
        psi(:) = CMPLX((psi_I(1:n)*intR + psi_R(1:n)*intI) / wronski, KIND=idk)

        DEALLOCATE(yR,yI,intR,intI,psi_I,psi_R)
    
    END SUBROUTINE apply_green_analytic_continuation
    
    

END MODULE Green



