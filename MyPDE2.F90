MODULE MyLinearSolver

  USE DefUtils
  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE my_matvec(x,r)
    REAL(KIND=dp) :: x(:), r(:)
    TYPE(Matrix_t), POINTER :: A
    A => CurrentModel % Solver % Matrix
    CALL CRS_MatrixVectorMultiply( A,x,r )
  END SUBROUTINE my_matvec
    
  SUBROUTINE my_rpcond(x,r)
    REAL(KIND=dp) :: x(:), r(:)
    TYPE(Matrix_t), POINTER :: A
    LOGICAL :: Found
    IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
        == 'diagonal' ) THEN
      A => CurrentModel % Solver % Matrix
      x = r / A % Values(A % Diag)
    ELSE
      x = r 
    END IF    
  END SUBROUTINE my_rpcond

  FUNCTION my_normfun(n, b) RESULT ( bnorm ) 
    REAL(KIND=dp) :: b(:)
    INTEGER :: n
    REAL(KIND=dp) :: bnorm
    bnorm = SQRT(SUM(b(1:n)**2))
  END FUNCTION my_normfun
          
  FUNCTION my_dotprodfun(n, t1, t2 ) RESULT ( beta ) 
    INTEGER :: n
    REAL(KIND=dp) :: t1(:), t2(:)        
    REAL(KIND=dp) :: beta    
    beta = SUM(t1(1:n)*t2(1:n))
  END FUNCTION my_dotprodfun
      
  ! Write something similar to this routine using the above functions.
  ! We can then embed this to the library rather easily. No need to change the variablenames
  ! but the they should be pretty much the same. Probably have some additional vector and logical
  ! for the upper/lower limit. 
  
  SUBROUTINE my_GCR( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
      Converged, Diverged, OutputInterval, m, MinIter) 
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: Rounds, MinIter
    REAL(KIND=dp) :: x(n),b(n)
    LOGICAL :: Converged, Diverged
    REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
    INTEGER :: n, OutputInterval, m
    REAL(KIND=dp) :: bnorm,rnorm
    REAL(KIND=dp), ALLOCATABLE :: R(:)    
    REAL(KIND=dp), ALLOCATABLE :: S(:,:), V(:,:), T1(:), T2(:)

!------------------------------------------------------------------------------
    INTEGER :: i,j,k
    REAL(KIND=dp) :: alpha, beta
!------------------------------------------------------------------------------
    INTEGER :: allocstat

    ALLOCATE( R(n), T1(n), T2(n), STAT=allocstat )
    IF( allocstat /= 0 ) THEN
      CALL Fatal('GCR','Failed to allocate memory of size: '//I2S(n))
    END IF

    IF ( m > 1 ) THEN
      ALLOCATE( S(n,m-1), V(n,m-1), STAT=allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal('GCR','Failed to allocate memory of size: '&
            //I2S(n)//' x '//I2S(m-1))
      END IF

      V(1:n,1:m-1) = 0.0d0	
      S(1:n,1:m-1) = 0.0d0
    END IF

!    CALL C_matvec( x, r, ipar, matvecsubr )
    CALL my_matvec( x, r )
    r(1:n) = b(1:n) - r(1:n)

!    bnorm = normfun(n, b, 1)
!    rnorm = normfun(n, r, 1)
    bnorm = my_normfun(n, b)
    rnorm = my_normfun(n, r)

    !IF (UseStopCFun) THEN
    !  Residual = stopcfun(x,b,r,ipar,dpar)
    !ELSE
      Residual = rnorm / bnorm
    !END IF

      Converged = (Residual < MinTolerance) .AND. ( MinIter <= 0 )
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
      IF( Converged .OR. Diverged) RETURN

    
    DO k=1,Rounds
      !----------------------------------------------
      ! Check for restarting
      !----------------------------------------------
      IF ( MOD(k,m)==0 ) THEN
        j = m
      ELSE
        j = MOD(k,m)
        !--------------------------------------------
        ! Compute the true residual when restarting:
        !--------------------------------------------
        IF ( (j==1) .AND. (k>1) ) THEN
!          CALL C_matvec( x, r, ipar, matvecsubr )
          CALL my_matvec( x, r )
          r(1:n) = b(1:n) - r(1:n)
        END IF
      END IF
      
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !---------------------------------------------------------------
      !CALL C_rpcond( T1, r, ipar, pcondrsubr )
      CALL my_rpcond( T1, r )
      !CALL C_matvec( T1, T2, ipar, matvecsubr )
      CALL my_matvec( T1, T2 )
      
      !--------------------------------------------------------------
      ! Perform the orthogonalization of the search directions....
      !--------------------------------------------------------------
      DO i=1,j-1
!        beta = dotprodfun(n, V(1:n,i), 1, T2(1:n), 1 )
        beta = my_dotprodfun(n, V(1:n,i), T2(1:n) )

        T1(1:n) = T1(1:n) - beta * S(1:n,i)
        T2(1:n) = T2(1:n) - beta * V(1:n,i)        
      END DO

!      alpha = normfun(n, T2(1:n), 1 )
      alpha = my_normfun(n, T2(1:n) )
      T1(1:n) = 1.0d0/alpha * T1(1:n)
      T2(1:n) = 1.0d0/alpha * T2(1:n)
      
      !-------------------------------------------------------------
      ! The update of the solution and save the search data...
      !------------------------------------------------------------- 
!      beta = dotprodfun(n, T2(1:n), 1, r(1:n), 1 )
      beta = my_dotprodfun(n, T2(1:n), r(1:n) )

      x(1:n) = x(1:n) + beta * T1(1:n)      
      r(1:n) = r(1:n) - beta * T2(1:n)

      IF ( j /= m ) THEN
        S(1:n,j) = T1(1:n)
        V(1:n,j) = T2(1:n)
      END IF

      !--------------------------------------------------------------
      ! Check whether the convergence criterion is met 
      !--------------------------------------------------------------
      !rnorm = normfun(n, r, 1)
      rnorm = my_normfun(n, r)

      !IF (UseStopCFun) THEN
        !Residual = stopcfun(x,b,r,ipar,dpar)
      !  IF( MOD(k,OutputInterval) == 0) THEN
      !    WRITE (*, '(A, I6, 2E12.4)') '   gcr:',k, rnorm / bnorm, residual
      !    CALL FLUSH(6)
      !  END IF
      !ELSE
        Residual = rnorm / bnorm
        IF( MOD(k,OutputInterval) == 0) THEN
          WRITE (*, '(A, I6, 2E12.4)') '   gcr:',k, residual, beta
          CALL FLUSH(6)
        END IF
      !END IF

      Converged = (Residual < MinTolerance) .AND. ( k >= MinIter )
      !-----------------------------------------------------------------
      ! Make an additional check that the true residual agrees with 
      ! the iterated residual:
      !-----------------------------------------------------------------
      IF (Converged ) THEN
        WRITE( Message,'(A,I0,A,ES12.3)') 'Iterated residual norm after ',k,' iters:', rnorm
        CALL Info('IterMethod_GCR', Message, Level=5)
        CALL Info('IterMethod_GCR','Total number of GCR iterations: '//I2S(k), Level=5)                     
      END IF
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)    
      IF( Converged .OR. Diverged) EXIT

    END DO

    DEALLOCATE( R, T1, T2 )
    IF ( m > 1 ) DEALLOCATE( S, V)

  END SUBROUTINE my_GCR

!------------------------------------------------------------------------------
  SUBROUTINE my_MPRGP(n, x, b, c, epsr, maxit, Gamma, precond, adapt, bound, &
                      ncg, ne, np, iters, converged, final_norm_gp)
    USE DefUtils
    IMPLICIT NONE

    ! ---------------------------
    ! Arguments
    ! ---------------------------
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), INTENT(INOUT) :: x(n)    ! initial guess in, final solution out
    REAL(KIND=dp), INTENT(IN) :: b(n), c(n) ! rhs and bound
    REAL(KIND=dp), INTENT(IN) :: epsr
    INTEGER, INTENT(IN) :: maxit
    REAL(KIND=dp), INTENT(IN) :: Gamma
    CHARACTER(*), INTENT(IN) :: precond    ! 'none' or 'jacobi'
    LOGICAL, INTENT(IN) :: adapt
    CHARACTER(*), INTENT(IN) :: bound      ! 'lower' or 'upper'

    INTEGER, INTENT(OUT) :: ncg, ne, np, iters
    LOGICAL, INTENT(OUT) :: converged
    REAL(KIND=dp), INTENT(OUT) :: final_norm_gp

    ! ---------------------------
    ! Local declarations (all here)
    ! ---------------------------
    REAL(KIND=dp), ALLOCATABLE :: g(:), gf(:), gc(:), gr(:), gp(:)
    REAL(KIND=dp), ALLOCATABLE :: z(:), p(:), Ap(:), yy(:), D(:)
    LOGICAL, ALLOCATABLE :: J(:)
    INTEGER :: bs
    REAL(KIND=dp) :: lAl, alpha, a_f
    INTEGER :: i, allocstat
    REAL(KIND=dp) :: rtp, pAp, acg, beta
    LOGICAL :: use_jacobi
    REAL(KIND=dp) :: tmp_norm
    REAL(KIND=dp) :: eps_local
    TYPE(Matrix_t), POINTER :: MatA

    ! ---------------------------
    ! Allocate scratch vectors
    ! ---------------------------
    ALLOCATE(g(n), gf(n), gc(n), gr(n), gp(n), z(n), p(n), Ap(n), yy(n), J(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('my_MPRGP','Allocation failed for scratch vectors')
    END IF

    ! ---------------------------
    ! Preliminaries
    ! ---------------------------
    eps_local = EPSILON(1.0_dp)

    ! decide jacobi preconditioner
    use_jacobi = .FALSE.
    IF (TRIM(ADJUSTL(precond)) == 'jacobi') THEN
      use_jacobi = .TRUE.
    END IF

    ! set bound sign
    IF (TRIM(ADJUSTL(bound)) == 'upper') THEN
      bs = -1
    ELSE
      bs = 1
    END IF

    ! set matrix pointer from CurrentModel like my_matvec does
    MatA => CurrentModel % Solver % Matrix

    ! Prepare diagonal D if jacobi requested
    IF (use_jacobi) THEN
      IF (ASSOCIATED(MatA)) THEN
        ALLOCATE(D(n), STAT=allocstat)
        IF (allocstat /= 0) THEN
          CALL Fatal('my_MPRGP','Failed to allocate D(n)')
        END IF
        ! Copy diagonal from matrix
        DO i = 1, n
          D(i) = MatA % Values(MatA % Diag(i))
          IF (ABS(D(i)) < EPSILON(1.0_dp)) D(i) = 1.0_dp
        END DO
      ELSE
        use_jacobi = .FALSE.
      END IF
    END IF

    ! ---------------------------
    ! Initialization (g = A*x - b)
    ! ---------------------------
    CALL my_matvec(x, g)       ! g = A*x
    DO i = 1, n
      g(i) = g(i) - b(i)
    END DO

    DO i = 1, n
      J(i) = (bs * x(i) > bs * c(i))
    END DO

    DO i = 1, n
      IF (J(i)) THEN
        gf(i) = g(i)
      ELSE
        gf(i) = 0.0_dp
      END IF
    END DO

    IF (bs == 1) THEN
      DO i = 1, n
        IF (.NOT. J(i)) THEN
          gc(i) = MIN(g(i), 0.0_dp)
        ELSE
          gc(i) = 0.0_dp
        END IF
      END DO
    ELSE
      DO i = 1, n
        IF (.NOT. J(i)) THEN
          gc(i) = MAX(g(i), 0.0_dp)
        ELSE
          gc(i) = 0.0_dp
        END IF
      END DO
    END IF

    ! estimate matrix norm ||A|| with a small power iteration
    CALL estimate_matrix_norm(n, 10, lAl, MatA)
    IF (lAl <= 0.0_dp) lAl = 1.0_dp
    alpha = 1.0_dp / lAl

    ! reduced free gradient gr
    IF (bs == 1) THEN
      DO i = 1, n
        IF (J(i)) THEN
          gr(i) = MIN(lAl * (x(i) - c(i)), gf(i))
        ELSE
          gr(i) = 0.0_dp
        END IF
      END DO
    ELSE
      DO i = 1, n
        IF (J(i)) THEN
          gr(i) = MAX(lAl * (x(i) - c(i)), gf(i))
        ELSE
          gr(i) = 0.0_dp
        END IF
      END DO
    END IF

    DO i = 1, n
      gp(i) = gf(i) + gc(i)
    END DO

    ! preconditioning: z = M^{-1} * g on free set
    IF (use_jacobi) THEN
      DO i = 1, n
        IF (J(i)) THEN
          z(i) = g(i) / D(i)
        ELSE
          z(i) = 0.0_dp
        END IF
      END DO
    ELSE
      DO i = 1, n
        IF (J(i)) THEN
          z(i) = g(i)
        ELSE
          z(i) = 0.0_dp
        END IF
      END DO
    END IF

    p = z

    ! counters
    ncg = 0
    ne = 0
    np = 0
    iters = 0
    converged = .FALSE.

    ! ---------------------------
    ! Main loop
    ! ---------------------------
    DO WHILE ( my_normfun(n, gp) > epsr .AND. iters < maxit )

      iters = iters + 1
      WRITE(*,*) "iteration", iters

      IF ( my_dotprodfun(n, gc, gc) <= (Gamma**2) * my_dotprodfun(n, gr, gf) ) THEN
        ! CG-like step
        CALL my_matvec(p, Ap)
        rtp = my_dotprodfun(n, z, g)
        pAp = my_dotprodfun(n, p, Ap)

        IF (ABS(pAp) < eps_local) THEN
          CALL Info('my_MPRGP','p''*A*p nearly zero, stopping',Level=5)
          EXIT
        END IF

        acg = rtp / pAp
        DO i = 1, n
          yy(i) = x(i) - acg * p(i)
        END DO

        IF (ALL(bs * yy(:) >= bs * c(:))) THEN
          ! accept full CG step
          x = yy
          DO i = 1, n
            g(i) = g(i) - acg * Ap(i)
            J(i) = (bs * x(i) > bs * c(i))
          END DO

          ! precondition
          IF (use_jacobi) THEN
            DO i = 1, n
              IF (J(i)) THEN
                z(i) = g(i) / D(i)
              ELSE
                z(i) = 0.0_dp
              END IF
            END DO
          ELSE
            DO i = 1, n
              IF (J(i)) THEN
                z(i) = g(i)
              ELSE
                z(i) = 0.0_dp
              END IF
            END DO
          END IF

          beta = my_dotprodfun(n, z, Ap) / pAp
          DO i = 1, n
            p(i) = z(i) - beta * p(i)
          END DO

          ! update gf,gc,gr,gp
          DO i = 1, n
            IF (J(i)) THEN
              gf(i) = g(i)
            ELSE
              gf(i) = 0.0_dp
            END IF
          END DO

          IF (bs == 1) THEN
            DO i = 1, n
              IF (.NOT. J(i)) THEN
                gc(i) = MIN(g(i), 0.0_dp)
              ELSE
                gc(i) = 0.0_dp
              END IF
            END DO
            DO i = 1, n
              IF (J(i)) THEN
                gr(i) = MIN(lAl * (x(i) - c(i)), gf(i))
              ELSE
                gr(i) = 0.0_dp
              END IF
            END DO
          ELSE
            DO i = 1, n
              IF (.NOT. J(i)) THEN
                gc(i) = MAX(g(i), 0.0_dp)
              ELSE
                gc(i) = 0.0_dp
              END IF
            END DO
            DO i = 1, n
              IF (J(i)) THEN
                gr(i) = MAX(lAl * (x(i) - c(i)), gf(i))
              ELSE
                gr(i) = 0.0_dp
              END IF
            END DO
          END IF

          DO i = 1, n
            gp(i) = gf(i) + gc(i)
          END DO
          ncg = ncg + 1

        ELSE
          ! expansion step (feasible step)
          a_f = -1.0_dp
          DO i = 1, n
            IF ( (bs * p(i) > 0.0_dp) .AND. J(i) ) THEN
              tmp_norm = (x(i) - c(i)) / p(i)
              IF (a_f < 0.0_dp) THEN
                a_f = tmp_norm
              ELSE
                a_f = MIN(a_f, tmp_norm)
              END IF
            END IF
          END DO
          IF (a_f < 0.0_dp) a_f = 0.0_dp

          IF (bs == 1) THEN
            DO i = 1, n
              x(i) = MAX(x(i) - a_f * p(i), c(i))
            END DO
          ELSE
            DO i = 1, n
              x(i) = MIN(x(i) - a_f * p(i), c(i))
            END DO
          END IF

          DO i = 1, n
            J(i) = (bs * x(i) > bs * c(i))
            g(i) = g(i) - a_f * Ap(i)
          END DO

          ! recompute preconditioned residual z
          IF (use_jacobi) THEN
            DO i = 1, n
              IF (J(i)) THEN
                z(i) = g(i) / D(i)
              ELSE
                z(i) = 0.0_dp
              END IF
            END DO
          ELSE
            DO i = 1, n
              IF (J(i)) THEN
                z(i) = g(i)
              ELSE
                z(i) = 0.0_dp
              END IF
            END DO
          END IF

          ! adaptive alpha if requested
          IF (adapt) THEN
            CALL my_matvec(gr, Ap)
            tmp_norm = my_dotprodfun(n, gr, g)
            pAp = my_dotprodfun(n, gr, Ap)
            IF (ABS(pAp) < eps_local) THEN
              alpha = 1.0_dp / lAl
            ELSE
              alpha = tmp_norm / pAp
              IF (alpha <= 0.0_dp .OR. alpha > 1.0_dp / lAl) alpha = 1.0_dp / lAl
            END IF
          END IF

          IF (bs == 1) THEN
            DO i = 1, n
              IF (J(i)) THEN
                x(i) = MAX(x(i) - alpha * g(i), c(i))
              END IF
            END DO
          ELSE
            DO i = 1, n
              IF (J(i)) THEN
                x(i) = MIN(x(i) - alpha * g(i), c(i))
              END IF
            END DO
          END IF

          CALL my_matvec(x, g)
          DO i = 1, n
            g(i) = g(i) - b(i)
          END DO

          ! recompute z and p
          IF (use_jacobi) THEN
            DO i = 1, n
              IF (J(i)) THEN
                z(i) = g(i) / D(i)
              ELSE
                z(i) = 0.0_dp
              END IF
            END DO
          ELSE
            DO i = 1, n
              IF (J(i)) THEN
                z(i) = g(i)
              ELSE
                z(i) = 0.0_dp
              END IF
            END DO
          END IF

          p = z

          ! update gf,gc,gr,gp
          DO i = 1, n
            IF (J(i)) THEN
              gf(i) = g(i)
            ELSE
              gf(i) = 0.0_dp
            END IF
          END DO

          IF (bs == 1) THEN
            DO i = 1, n
              IF (.NOT. J(i)) THEN
                gc(i) = MIN(g(i), 0.0_dp)
              ELSE
                gc(i) = 0.0_dp
              END IF
            END DO
            DO i = 1, n
              IF (J(i)) THEN
                gr(i) = MIN(lAl * (x(i) - c(i)), gf(i))
              ELSE
                gr(i) = 0.0_dp
              END IF
            END DO
          ELSE
            DO i = 1, n
              IF (.NOT. J(i)) THEN
                gc(i) = MAX(g(i), 0.0_dp)
              ELSE
                gc(i) = 0.0_dp
              END IF
            END DO
            DO i = 1, n
              IF (J(i)) THEN
                gr(i) = MAX(lAl * (x(i) - c(i)), gf(i))
              ELSE
                gr(i) = 0.0_dp
              END IF
            END DO
          END IF

          DO i = 1, n
            gp(i) = gf(i) + gc(i)
          END DO
          ne = ne + 1
        END IF

      ELSE
        ! proportioning step
        CALL my_matvec(gc, Ap)
        pAp = my_dotprodfun(n, gc, Ap)
        IF (ABS(pAp) < eps_local) THEN
          CALL Info('my_MPRGP','denominator in proportioning nearly zero, stopping',Level=5)
          EXIT
        END IF
        acg = my_dotprodfun(n, gc, g) / pAp
        DO i = 1, n
          x(i) = x(i) - acg * gc(i)
        END DO

        IF (bs == 1) THEN
          DO i = 1, n
            IF (x(i) < c(i)) x(i) = c(i)
          END DO
        ELSE
          DO i = 1, n
            IF (x(i) > c(i)) x(i) = c(i)
          END DO
        END IF

        DO i = 1, n
          J(i) = (bs * x(i) > bs * c(i))
          g(i) = g(i) - acg * Ap(i)
        END DO

        IF (use_jacobi) THEN
          DO i = 1, n
            IF (J(i)) THEN
              z(i) = g(i) / D(i)
            ELSE
              z(i) = 0.0_dp
            END IF
          END DO
        ELSE
          DO i = 1, n
            IF (J(i)) THEN
              z(i) = g(i)
            ELSE
              z(i) = 0.0_dp
            END IF
          END DO
        END IF

        p = z

        ! update gf,gc,gr,gp
        DO i = 1, n
          IF (J(i)) THEN
            gf(i) = g(i)
          ELSE
            gf(i) = 0.0_dp
          END IF
        END DO

        IF (bs == 1) THEN
          DO i = 1, n
            IF (.NOT. J(i)) THEN
              gc(i) = MIN(g(i), 0.0_dp)
            ELSE
              gc(i) = 0.0_dp
            END IF
          END DO
          DO i = 1, n
            IF (J(i)) THEN
              gr(i) = MIN(lAl * (x(i) - c(i)), gf(i))
            ELSE
              gr(i) = 0.0_dp
            END IF
          END DO
        ELSE
          DO i = 1, n
            IF (.NOT. J(i)) THEN
              gc(i) = MAX(g(i), 0.0_dp)
            ELSE
              gc(i) = 0.0_dp
            END IF
          END DO
          DO i = 1, n
            IF (J(i)) THEN
              gr(i) = MAX(lAl * (x(i) - c(i)), gf(i))
            ELSE
              gr(i) = 0.0_dp
            END IF
          END DO
        END IF

        DO i = 1, n
          gp(i) = gf(i) + gc(i)
        END DO
        np = np + 1
      END IF

    END DO  ! main loop

    final_norm_gp = my_normfun(n, gp)
    converged = (final_norm_gp <= epsr)

    ! cleanup
    IF (ALLOCATED(D)) DEALLOCATE(D)
    DEALLOCATE(g, gf, gc, gr, gp, z, p, Ap, yy, J)

    RETURN
  CONTAINS
    ! power-method norm estimator using A
    SUBROUTINE estimate_matrix_norm(nloc, niter, out_norm, Aptr)
      INTEGER, INTENT(IN) :: nloc, niter
      REAL(KIND=dp), INTENT(OUT) :: out_norm
      TYPE(Matrix_t), POINTER :: Aptr
      REAL(KIND=dp), ALLOCATABLE :: v(:), w(:)
      INTEGER :: itl, allocstat
      REAL(KIND=dp) :: normv

      ALLOCATE(v(nloc), w(nloc), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('estimate_matrix_norm','alloc fail')
      END IF

      DO itl = 1, nloc
        v(itl) = 1.0_dp / REAL(nloc, dp)
      END DO

      DO itl = 1, niter
        CALL CRS_MatrixVectorMultiply(Aptr, v, w)
        normv = my_normfun(nloc, w)
        IF (normv <= 0.0_dp) EXIT
        v = w / normv
      END DO

      CALL CRS_MatrixVectorMultiply(Aptr, v, w)
      out_norm = my_normfun(nloc, w)

      DEALLOCATE(v, w)
    END SUBROUTINE estimate_matrix_norm

  END SUBROUTINE my_MPRGP
  
END MODULE MyLinearSolver
  


!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE MyLinearSolver
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: i, n, nb, nd, t, active
  INTEGER :: iter, maxiter
  LOGICAL :: Found, UpperLim, LowerLim  
  TYPE(ValueList_t), POINTER :: Params  
  REAL(KIND=dp), ALLOCATABLE :: LimVal(:)
  
!------------------------------------------------------------------------------

  CALL DefaultStart()

  Params => GetSolverParams()
  maxiter = ListGetInteger( Params,'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  LowerLim = ListCheckPresentAnyBodyForce(Model,'Temp Lower Limit')
  UpperLim = ListCheckPresentAnyBodyForce(Model,'Temp Upper Limit')
  IF(LowerLim .AND. UpperLim) THEN
    CALL Warn('AdvDiffSolver','This code cannot have both upper and lower limit at same time!')
  END IF

  IF(UpperLim .OR. LowerLim) THEN
    n = SIZE(Solver % Variable % Values)
    ALLOCATE(LimVal(n))
    IF(UpperLim) THEN
      LimVal = HUGE(Norm)
    ELSE
      LimVal = -HUGE(Norm)
    END IF
  END IF
    
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()

1   Active = GetNOFActive()

    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()
      CALL LocalMatrix(  Element, n, nd+nb )
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        CALL LocalMatrixBC(  Element, n, nd )
      END IF
    END DO

    IF(DefaultCutFEM()) GOTO 1
    
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    
    ! And finally, solve:
    !--------------------
  IF(ListGetLogical(Solver % Values,'My linear solver',Found ) ) THEN
    ! Use custom solver - choose between GCR and MPRGP
    
    ! Use MPRGP for bound-constrained problems
    BLOCK
      INTEGER :: n_mprgp, ncg, ne, np, iters
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp), POINTER :: x(:), b(:)
      REAL(KIND=dp), ALLOCATABLE :: c(:)
      LOGICAL :: converged_mprgp
      REAL(KIND=dp) :: final_norm_gp
      CHARACTER(LEN=10) :: bound_type
      
      n_mprgp = SIZE(Solver % Variable % Values)
      A => Solver % Matrix
      b => Solver % Matrix % Rhs
      x => Solver % Variable % Values
      
      ! Prepare bound constraint vector
      ALLOCATE(c(n_mprgp))
      IF(UpperLim) THEN
        c = LimVal
        bound_type = 'upper'
      ELSE
        c = LimVal
        bound_type = 'lower'
      END IF
      
      CALL my_MPRGP(n_mprgp, x, b, c, 1.0e-8_dp, 100000, 1.0_dp, 'jacobi', .FALSE., &
            bound_type, ncg, ne, np, iters, converged_mprgp, final_norm_gp)
      
      CALL Info('AdvDiffSolver','MPRGP finished: iters='//I2S(iters)// &
                ', ncg='//I2S(ncg)//', ne='//I2S(ne)//', np='//I2S(np), Level=5)
      ! Save the x

      OPEN(1,FILE="x_mprgp.dat", STATUS='Unknown')
      DO i=1,SIZE(x)
        WRITE(1,*) x(i)    
      END DO
      CLOSE(1)

      DEALLOCATE(c)
    END BLOCK
    
  ELSE
    Norm = DefaultSolve()      
    IF( DefaultConverged() ) EXIT    
  END IF

  END DO

  ! Save the limit
  IF(UpperLim .OR. LowerLim) THEN
    OPEN(1,FILE="lim.dat", STATUS='Unknown')
    DO i=1,SIZE(LimVal)
      WRITE(1,*) LimVal(i)    
    END DO
    CLOSE(1)
  END IF
    
  CALL DefaultFinish()
  
CONTAINS

  

  
! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: diff_coeff(n), conv_coeff(n),react_coeff(n), &
                     time_coeff(n), D,C,R, rho,Velo(3,n),a(3), Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n), Lim(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce()
    
    IF ( ASSOCIATED(BodyForce) ) THEN
      Load(1:n) = GetReal( BodyForce,'field source', Found )
      Found = .FALSE.
      IF(UpperLim) THEN
        Lim(1:n) = GetReal(BodyForce,'Temp Upper Limit',Found) 
      ELSE IF(LowerLim) THEN
        Lim(1:n) = GetReal(BodyForce,'Temp Lower Limit',Found) 
      END IF
      IF(Found) THEN
        LimVal(Solver % Variable % Perm(Element % NodeIndexes)) = Lim(1:n)
      END IF
    END IF
    
       
    Material => GetMaterial()
    diff_coeff(1:n)=GetReal(Material,'diffusion coefficient',Found)
    react_coeff(1:n)=GetReal(Material,'reaction coefficient',Found)
    conv_coeff(1:n)=GetReal(Material,'convection coefficient',Found)
    time_coeff(1:n)=GetReal(Material,'time derivative coefficient',Found)

    
    Velo = 0._dp
    DO i=1,dim
      Velo(i,1:n)=GetReal(Material,&
          'convection velocity '//I2S(i),Found)
    END DO

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    IF( Element % ElementIndex == 1 ) THEN
      CALL Info('AdvDiffSolver','Integration points in 1st element: '//I2S(IP % n),Level=8)
    END IF


    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      rho = SUM(Basis(1:n)*time_coeff(1:n))
      a = MATMUL(Velo(:,1:n),Basis(1:n))
      D = SUM(Basis(1:n)*diff_coeff(1:n))
      C = SUM(Basis(1:n)*conv_coeff(1:n))
      R = SUM(Basis(1:n)*react_coeff(1:n))

      Weight = IP % s(t) * DetJ

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
             D * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      DO p=1,nd
        DO q=1,nd
          ! advection term (C*grad(u),v)
          ! -----------------------------------
          STIFF (p,q) = STIFF(p,q) + Weight * &
             C * SUM(a(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! reaction term (R*u,v)
          ! -----------------------------------
          STIFF(p,q) = STIFF(p,q) + Weight * R*Basis(q) * Basis(p)

          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * rho * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BC,'field flux', Found )
    Coeff(1:n) = GetReal( BC,'robin coefficient', Found )
    Ext_t(1:n) = GetReal( BC,'external field', Found )

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = SUM(Basis(1:n)*flux(1:n))

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
      C = SUM(Basis(1:n)*coeff(1:n))
      Ext = SUM(Basis(1:n)*ext_t(1:n))

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
    END DO
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------
