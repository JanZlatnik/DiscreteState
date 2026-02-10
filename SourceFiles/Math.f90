!-------------------------------------< Math >-------------------------------------!
!                                                                                  !    
! Contains: Basic numerical operations like integration, primitive function,       !
!           reversing order of vectors, special functions, solve linear system     !
!           of equations                                                           !
!                                                                                  !
! Last revision:    23/03/2024                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!

MODULE Math
    USE Parameters
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: primitive, definite_integral, REVERSE, h1F1, gammaincup, solve_SLE
    
    ! Interface for function f being complex / real
    INTERFACE primitive
        MODULE PROCEDURE primitive_real
        MODULE PROCEDURE primitive_complex
        MODULE PROCEDURE primitive_real_qdk
        MODULE PROCEDURE primitive_complex_qdk
    END INTERFACE primitive

    INTERFACE definite_integral
        MODULE PROCEDURE definite_integral_real
        MODULE PROCEDURE definite_integral_complex
    END INTERFACE definite_integral
    
    ! Interface for vector being complex / real
    INTERFACE REVERSE
        MODULE PROCEDURE REVERSE_real
        MODULE PROCEDURE REVERSE_complex
        MODULE PROCEDURE REVERSE_real_qdk
        MODULE PROCEDURE REVERSE_complex_qdk
    END INTERFACE REVERSE
    
    ! Interface for Ax=b (solve_SLE) A,b being real / complex 
    INTERFACE solve_SLE
        MODULE PROCEDURE solve_SLE_r
        MODULE PROCEDURE solve_SLE_c
        MODULE PROCEDURE solve_SLE_r_m
        MODULE PROCEDURE solve_SLE_c_m
    END INTERFACE solve_SLE
    
    CONTAINS
    
    ! Primitive function for real function f using trapezoid rule
    SUBROUTINE primitive_real(f, dx, primitive_f)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: f(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        REAL(KIND=idk), ALLOCATABLE, INTENT(OUT) :: primitive_f(:)
        INTEGER :: i, n

        n = SIZE(f)
        IF (ALLOCATED(primitive_f)) DEALLOCATE(primitive_f)
        ALLOCATE(primitive_f(n))
        primitive_f(1) = 0.0d0
        primitive_f(2) = f(1) * 3.0d0/8.0d0 + f(2) * 19.0d0/24.0d0 - f(3) * 5.0d0/24.0d0 + f(4) * 1.0d0/24.0d0
        DO i = 3, n-1
            primitive_f(i) = primitive_f(i-1) + (f(i-1)+f(i)) * 13.0d0/24.0d0 - (f(i-2)+f(i+1)) / 24.0d0
        END DO
        primitive_f(n) = primitive_f(n-1) + f(n) * 3.0d0/8.0d0 + f(n-1) * 19.0d0/24.0d0 - f(n-3) * 5.0d0/24.0d0 + f(n-3) * 1.0d0/24.0d0
        primitive_f = primitive_f * dx
    END SUBROUTINE primitive_real
    
    
    ! Primitive function for complex function f using trapezoid rule
    SUBROUTINE primitive_complex(f, dx, primitive_f)
        IMPLICIT NONE
        COMPLEX(KIND=idk), INTENT(IN) :: f(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        COMPLEX(KIND=idk), ALLOCATABLE, INTENT(OUT) :: primitive_f(:)
        INTEGER :: i, n

        n = SIZE(f)
        IF (ALLOCATED(primitive_f)) DEALLOCATE(primitive_f)
        ALLOCATE(primitive_f(n))
        primitive_f(1) = CMPLX(0.0d0, 0.0d0, KIND=idk)
        primitive_f(2) = f(1) * 3.0d0/8.0d0 + f(2) * 19.0d0/24.0d0 - f(3) * 5.0d0/24.0d0 + f(4) * 1.0d0/24.0d0
        DO i = 3, n-1
            primitive_f(i) = primitive_f(i-1) + (f(i-1)+f(i)) * 13.0d0/24.0d0 - (f(i-2)+f(i+1)) / 24.0d0
        END DO
        primitive_f(n) = primitive_f(n-1) + f(n) * 3.0d0/8.0d0 + f(n-1) * 19.0d0/24.0d0 - f(n-3) * 5.0d0/24.0d0 + f(n-3) * 1.0d0/24.0d0
        primitive_f = primitive_f * dx
    END SUBROUTINE primitive_complex
    
    
    ! Primitive function for real function f using trapezoid rule
    SUBROUTINE primitive_real_qdk(f, dx, primitive_f)
        IMPLICIT NONE
        REAL(KIND=qdk), INTENT(IN) :: f(:)
        REAL(KIND=qdk), INTENT(IN) :: dx
        REAL(KIND=qdk), ALLOCATABLE, INTENT(OUT) :: primitive_f(:)
        INTEGER :: i, n

        n = SIZE(f)
        IF (ALLOCATED(primitive_f)) DEALLOCATE(primitive_f)
        ALLOCATE(primitive_f(n))
        primitive_f(1) = 0.0_qdk
        primitive_f(2) = f(1) * 3.0_qdk/8.0_qdk + f(2) * 19.0_qdk/24.0_qdk - f(3) * 5.0_qdk/24.0_qdk + f(4) * 1.0_qdk/24.0_qdk
        DO i = 3, n-1
            primitive_f(i) = primitive_f(i-1) + (f(i-1)+f(i)) * 13.0_qdk/24.0_qdk - (f(i-2)+f(i+1)) / 24.0_qdk
        END DO
        primitive_f(n) = primitive_f(n-1) + f(n) * 3.0_qdk/8.0_qdk + f(n-1) * 19.0_qdk/24.0_qdk - f(n-3) * 5.0_qdk/24.0_qdk + f(n-3) * 1.0_qdk/24.0_qdk
        primitive_f = primitive_f * dx
    END SUBROUTINE primitive_real_qdk
    
    
    ! Primitive function for complex function f using trapezoid rule
    SUBROUTINE primitive_complex_qdk(f, dx, primitive_f)
        IMPLICIT NONE
        COMPLEX(KIND=qdk), INTENT(IN) :: f(:)
        REAL(KIND=qdk), INTENT(IN) :: dx
        COMPLEX(KIND=qdk), ALLOCATABLE, INTENT(OUT) :: primitive_f(:)
        INTEGER :: i, n

        n = SIZE(f)
        IF (ALLOCATED(primitive_f)) DEALLOCATE(primitive_f)
        ALLOCATE(primitive_f(n))
        primitive_f(1) = CMPLX(0.0_qdk, 0.0_qdk, KIND=qdk)
        primitive_f(2) = f(1) * 3.0_qdk/8.0_qdk + f(2) * 19.0_qdk/24.0_qdk - f(3) * 5.0_qdk/24.0_qdk + f(4) * 1.0_qdk/24.0_qdk
        DO i = 3, n-1
            primitive_f(i) = primitive_f(i-1) + (f(i-1)+f(i)) * 13.0_qdk/24.0_qdk - (f(i-2)+f(i+1)) / 24.0_qdk
        END DO
        primitive_f(n) = primitive_f(n-1) + f(n) * 3.0_qdk/8.0_qdk + f(n-1) * 19.0_qdk/24.0_qdk - f(n-3) * 5.0_qdk/24.0_qdk + f(n-3) * 1.0_qdk/24.0_qdk
        primitive_f = primitive_f * dx
    END SUBROUTINE primitive_complex_qdk

    
    ! Definite integral for real function f using trapezoid rule
    SUBROUTINE definite_integral_real(f, dx, integral)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: f(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        REAL(KIND=idk), INTENT(OUT) :: integral
        INTEGER :: n

        n = SIZE(f)
        integral = (SUM(f)-(f(1) + f(n))*5.0d0/8.0d0 + (f(2)+f(n-1))/6.0d0 - (f(3)+f(n-2))/24.0d0)*dx
    END SUBROUTINE definite_integral_real

    ! Definite integral for complex function f using trapezoid rule
    SUBROUTINE definite_integral_complex(f, dx, integral)
        IMPLICIT NONE
        COMPLEX(KIND=idk), INTENT(IN) :: f(:)
        REAL(KIND=idk), INTENT(IN) :: dx
        COMPLEX(KIND=idk), INTENT(OUT) :: integral
        INTEGER :: n

        n = SIZE(f)
        integral = (SUM(f)-(f(1) + f(n))*5.0d0/8.0d0 + (f(2)+f(n-1))/6.0d0 - (f(3)+f(n-2))/24.0d0)*dx
    END SUBROUTINE definite_integral_complex
    
    
    ! Reversing vectors
    SUBROUTINE REVERSE_real(arr)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(INOUT) :: arr(:)
        INTEGER :: n, i
        REAL(KIND=idk) :: temp

        n = SIZE(arr)

        DO i = 1, n/2
            temp = arr(i)
            arr(i) = arr(n+1-i)
            arr(n+1-i) = temp
        END DO
    END SUBROUTINE REVERSE_real

    SUBROUTINE REVERSE_complex(arr)
        IMPLICIT NONE
        COMPLEX(KIND=idk), INTENT(INOUT) :: arr(:)
        INTEGER :: n, i
        COMPLEX(KIND=idk) :: temp

        n = SIZE(arr)

        DO i = 1, n/2
            temp = arr(i)
            arr(i) = arr(n+1-i)
            arr(n+1-i) = temp
        END DO
    END SUBROUTINE REVERSE_complex
    
    SUBROUTINE REVERSE_real_qdk(arr)
        IMPLICIT NONE
        REAL(KIND=qdk), INTENT(INOUT) :: arr(:)
        INTEGER :: n, i
        REAL(KIND=qdk) :: temp

        n = SIZE(arr)

        DO i = 1, n/2
            temp = arr(i)
            arr(i) = arr(n+1-i)
            arr(n+1-i) = temp
        END DO
    END SUBROUTINE REVERSE_real_qdk

    SUBROUTINE REVERSE_complex_qdk(arr)
        IMPLICIT NONE
        COMPLEX(KIND=qdk), INTENT(INOUT) :: arr(:)
        INTEGER :: n, i
        COMPLEX(KIND=qdk) :: temp

        n = SIZE(arr)

        DO i = 1, n/2
            temp = arr(i)
            arr(i) = arr(n+1-i)
            arr(n+1-i) = temp
        END DO
    END SUBROUTINE REVERSE_complex_qdk
    
    
    ! Confluent hypergeometric function 1F1
    FUNCTION h1F1(alpha,beta,x)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: alpha, beta, x
        INTEGER :: n, max_iter
        REAL(KIND=idk) :: term, h1F1, precision, a, b
        
        ! Precision setings
        max_iter = 10000
        precision = 1.0d-10
        
        ! Series
        a = alpha
        b = beta
        h1F1 = 1.0d0
        term = 1.0d0
        
        DO n = 1, max_iter
            term = term * a * x / (b * n)
            IF (ABS(term) < precision) EXIT
            a = a + 1
            b = b + 1
            h1F1 = h1F1 + term
        END DO
        IF (n - 1 == max_iter) PRINT*, "Confluent hypergeometric function 1F1 does NOT CONVERGE"
        RETURN
    END FUNCTION h1F1
    
    ! Upper incomplete gamma function using chain fraction
    FUNCTION gammaincup(a,x)
        IMPLICIT NONE
        REAL(KIND=idk), INTENT(IN) :: a, x
        INTEGER :: n, max_iter
        REAL(KIND=idk) :: frac, gammaincup
        
        ! Precision setings
        max_iter = 20
        
        ! Series
        gammaincup = EXP(-x) * x**a
        
        frac = x + (max_iter+1-a) / (2+max_iter)
        
        DO n = max_iter, 1, -1
            frac = x + (n-a)/(1+n/frac)
        END DO
        
        gammaincup = gammaincup / frac
        RETURN
    END FUNCTION gammaincup
    
    
    ! Solve system of linear equations Ax=b where A is a square matrix, A and b are real
    SUBROUTINE solve_SLE_r(coeffMatrix,rhsVec,solutionVec)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: coeffMatrix(:,:), rhsVec(:)
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: solutionVec(:)
        REAL(KIND = idk), ALLOCATABLE :: coeffMatrixCopy(:,:)
        INTEGER :: order_n, error_info
        INTEGER, ALLOCATABLE :: IPVT(:)
        
        EXTERNAL DGESV
        
        order_n = SIZE(rhsVec)
        IF (ALLOCATED(solutionVec)) DEALLOCATE(solutionVec)
        ALLOCATE(solutionVec(order_n),coeffMatrixCopy(order_n,order_n),IPVT(order_n))
        solutionVec = rhsVec
        coeffMatrixCopy = coeffMatrix
        
        CALL DGESV(order_n,1,coeffMatrixCopy,order_n,IPVT,solutionVec,order_n,error_info)
        
        IF (error_info /= 0) THEN
            PRINT *, 'DGESV failed with INFO:', error_info
        END IF
        
        DEALLOCATE(coeffMatrixCopy,IPVT)
    
    END SUBROUTINE solve_SLE_r
    
    ! Solve system of linear equations Ax=b where A is a square matrix, A and b are complex
    SUBROUTINE solve_SLE_c(coeffMatrix,rhsVec,solutionVec)
        IMPLICIT NONE
        COMPLEX(KIND = idk), INTENT(IN) :: coeffMatrix(:,:), rhsVec(:)
        COMPLEX(KIND = idk), ALLOCATABLE, INTENT(OUT) :: solutionVec(:)
        COMPLEX(KIND = idk), ALLOCATABLE :: coeffMatrixCopy(:,:)
        INTEGER :: order_n, error_info
        INTEGER, ALLOCATABLE :: IPVT(:)
        EXTERNAL ZGESV
        
        order_n = SIZE(rhsVec)
        IF (ALLOCATED(solutionVec)) DEALLOCATE(solutionVec)
        ALLOCATE(solutionVec(order_n),coeffMatrixCopy(order_n,order_n),IPVT(order_n))
        solutionVec = rhsVec
        coeffMatrixCopy = coeffMatrix
        
        CALL ZGESV(order_n,1,coeffMatrixCopy,order_n,IPVT,solutionVec,order_n,error_info)
        
        IF (error_info /= 0) THEN
            PRINT *, 'ZGESV failed with INFO:', error_info
        END IF
        
        DEALLOCATE(coeffMatrixCopy,IPVT)
    
    END SUBROUTINE solve_SLE_c
    
    ! Solve system of linear equations Ax=b where A is a square matrix for multiple rhs vectors b, A and b are real
    SUBROUTINE solve_SLE_r_m(coeffMatrix,rhsVecs,solution)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: coeffMatrix(:,:), rhsVecs(:,:)
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: solution(:,:)
        REAL(KIND = idk), ALLOCATABLE :: LU(:,:), solutionVecs(:,:)
        INTEGER :: order_n, num_rhs, error_info
        INTEGER, ALLOCATABLE :: IPVT(:)
        
        EXTERNAL DGETRF, DGETRS

        order_n = SIZE(coeffMatrix, 1)
        num_rhs = SIZE(rhsVecs, 1)
        IF (ALLOCATED(solution)) DEALLOCATE(solution)
        ALLOCATE(solutionVecs(order_n, num_rhs), LU(order_n, order_n), IPVT(order_n))
        LU = coeffMatrix
        solutionVecs = TRANSPOSE(rhsVecs)

        CALL DGETRF(order_n, order_n, LU, order_n, IPVT, error_info)
        IF (error_info /= 0) THEN
            PRINT *, 'DGETRF failed with INFO:', error_info
            RETURN
        END IF

        CALL DGETRS('N', order_n, num_rhs, LU, order_n, IPVT, solutionVecs, order_n, error_info)
        IF (error_info /= 0) THEN
            PRINT *, 'DGETRS failed with INFO:', error_info
        END IF

        ALLOCATE(solution(num_rhs,order_n))
        solution = TRANSPOSE(solutionVecs)
        DEALLOCATE(LU, IPVT)
        

    END SUBROUTINE solve_SLE_r_m


    
    ! Solve system of linear equations Ax=b where A is a square matrix for multiple rhs vectors b, A and b are complex
    SUBROUTINE solve_SLE_c_m(coeffMatrix,rhsVecs,solution)
        IMPLICIT NONE
        COMPLEX(KIND = idk), INTENT(IN) :: coeffMatrix(:,:), rhsVecs(:,:)
        COMPLEX(KIND = idk), ALLOCATABLE, INTENT(OUT) :: solution(:,:)
        COMPLEX(KIND = idk), ALLOCATABLE :: LU(:,:), solutionVecs(:,:)
        INTEGER :: order_n, num_rhs, error_info
        INTEGER, ALLOCATABLE :: IPVT(:)
        
        EXTERNAL ZGETRF, ZGETRS

        order_n = SIZE(coeffMatrix, 1)
        num_rhs = SIZE(rhsVecs, 1)
        IF (ALLOCATED(solution)) DEALLOCATE(solution)
        ALLOCATE(solutionVecs(order_n, num_rhs), LU(order_n, order_n), IPVT(order_n))
        LU = coeffMatrix
        solutionVecs = TRANSPOSE(rhsVecs)

        CALL ZGETRF(order_n, order_n, LU, order_n, IPVT, error_info)
        IF (error_info /= 0) THEN
            PRINT *, 'ZGETRF failed with INFO:', error_info
            RETURN
        END IF

        CALL ZGETRS('N', order_n, num_rhs, LU, order_n, IPVT, solutionVecs, order_n, error_info)
        IF (error_info /= 0) THEN
            PRINT *, 'ZGETRS failed with INFO:', error_info
        END IF
        
        ALLOCATE(solution(num_rhs,order_n))
        solution = TRANSPOSE(solutionVecs)
        DEALLOCATE(LU, IPVT)

    END SUBROUTINE solve_SLE_c_m





END MODULE Math
