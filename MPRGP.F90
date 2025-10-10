!------------------------------------------------------------------------------
! Optimized MPRGP (Modified Projected Refined Gradient Projection) solver
! Uses Elmer's efficient built-in functions for matrix-vector products, 
! norms, and dot products for optimal performance.
!------------------------------------------------------------------------------
SUBROUTINE my_MPRGP(n, x, b, c, epsr, maxit, Gamma, adapt, bound, &
                    ncg, ne, np, iters, converged, final_norm_gp)
!------------------------------------------------------------------------------
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
    REAL(KIND=dp) :: tmp_norm
    REAL(KIND=dp) :: eps_local
    TYPE(Matrix_t), POINTER :: MatA
    LOGICAL :: Found

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

    ! set bound sign
    IF (TRIM(ADJUSTL(bound)) == 'upper') THEN
      bs = -1
    ELSE
      bs = 1
    END IF

    ! set matrix pointer from CurrentModel like my_matvec does
    MatA => CurrentModel % Solver % Matrix

    ! ---------------------------
    ! Initialization (g = A*x - b)
    ! ---------------------------
    CALL MatrixVectorMultiply(MatA, x, g)       ! g = A*x

    g = g - b
    J = (bs * x > bs * c)
    gf = MERGE(g, 0.0_dp, J)

    IF (bs == 1) THEN
      WHERE (.NOT. J)
        gc = MIN(g, 0.0_dp)
      ELSEWHERE
        gc = 0.0_dp
      END WHERE
    ELSE
      WHERE (.NOT. J)
        gc = MAX(g, 0.0_dp)
      ELSEWHERE
        gc = 0.0_dp
      END WHERE
    END IF

    ! estimate matrix norm ||A|| with a small power iteration
    CALL estimate_matrix_norm(n, 10, lAl, MatA)
    IF (lAl <= 0.0_dp) lAl = 1.0_dp
    alpha = 1.0_dp / lAl

    ! reduced free gradient gr
    IF (bs == 1) THEN
      gr = MERGE(MIN(lAl * (x - c), gf), 0.0_dp, J)
    ELSE
      gr = MERGE(MAX(lAl * (x - c), gf), 0.0_dp, J)
    END IF

    gp = gf + gc

    ! preconditioning: z = M^{-1} * g on free set
    ! Apply diagonal preconditioning with free set masking
    IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
        == 'diagonal' ) THEN
      WHERE (J)
        WHERE (ABS(MatA % Values(MatA % Diag)) > AEPS)
          z = g / MatA % Values(MatA % Diag)
        ELSEWHERE
          z = g
        END WHERE
      ELSEWHERE
        z = 0.0_dp
      END WHERE
    ELSE
      ! No preconditioning, but still mask by free set
      WHERE (J)
        z = g
      ELSEWHERE
        z = 0.0_dp
      END WHERE
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
    DO WHILE ( ComputeNorm(CurrentModel % Solver, n, gp) > epsr .AND. iters < maxit )
      iters = iters + 1
      IF ( DOT_PRODUCT(gc(1:n), gc(1:n)) <= (Gamma**2) * DOT_PRODUCT(gr(1:n), gf(1:n)) ) THEN
        ! CG-like step
        CALL MatrixVectorMultiply(MatA, p, Ap)
        rtp = DOT_PRODUCT(z(1:n), g(1:n))
        pAp = DOT_PRODUCT(p(1:n), Ap(1:n))

        IF (ABS(pAp) < eps_local) THEN
          CALL Info('my_MPRGP','p''*A*p nearly zero, stopping',Level=5)
          EXIT
        END IF

        acg = rtp / pAp
        yy = x - acg * p


        IF (ALL(bs * yy(:) >= bs * c(:))) THEN
          ! accept full CG step
          x = yy
          g = g - acg * Ap
          J = (bs * x > bs * c)

          ! precondition
          ! Apply diagonal preconditioning with free set masking
          IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
              == 'diagonal' ) THEN
            WHERE (J)
              WHERE (ABS(MatA % Values(MatA % Diag)) > AEPS)
                z = g / MatA % Values(MatA % Diag)
              ELSEWHERE
                z = g
              END WHERE
            ELSEWHERE
              z = 0.0_dp
            END WHERE
          ELSE
            ! No preconditioning, but still mask by free set
            WHERE (J)
              z = g
            ELSEWHERE
              z = 0.0_dp
            END WHERE
          END IF

          beta = DOT_PRODUCT(z(1:n), Ap(1:n)) / pAp
          p = z - beta * p

          ! update gf,gc,gr,gp
          gf = MERGE(g, 0.0_dp, J)

          IF (bs == 1) THEN
            gc = 0.0_dp
            WHERE (.NOT. J)
              gc = MIN(g, 0.0_dp)
            END WHERE
            gr = MERGE(MIN(lAl * (x - c), gf), 0.0_dp, J)
          ELSE
            gc = 0.0_dp
            WHERE (.NOT. J)
              gc = MAX(g, 0.0_dp)
            END WHERE
            gr = MERGE(MAX(lAl * (x - c), gf), 0.0_dp, J)
          END IF
          gp = gf + gc
          ncg = ncg + 1
        ELSE
          ! expansion step (feasible step)
          ! TODO revise if this is correct
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
              x = MAX(x - a_f * p, c)
          ELSE
              x = MIN(x - a_f * p, c)
          END IF

          J = (bs * x > bs * c)
          g = g - a_f * Ap

          ! recompute preconditioned residual z
          ! Apply diagonal preconditioning with free set masking
          IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
              == 'diagonal' ) THEN
            WHERE (J)
              WHERE (ABS(MatA % Values(MatA % Diag)) > AEPS)
                z = g / MatA % Values(MatA % Diag)
              ELSEWHERE
                z = g
              END WHERE
            ELSEWHERE
              z = 0.0_dp
            END WHERE
          ELSE
            ! No preconditioning, but still mask by free set
            WHERE (J)
              z = g
            ELSEWHERE
              z = 0.0_dp
            END WHERE
          END IF

          ! adaptive alpha if requested
          IF (adapt) THEN
            CALL MatrixVectorMultiply(MatA, gr, Ap)
            tmp_norm = DOT_PRODUCT(gr(1:n), g(1:n))
            pAp = DOT_PRODUCT(gr(1:n), Ap(1:n))
            IF (ABS(pAp) < eps_local) THEN
              alpha = 1.0_dp / lAl
            ELSE
              alpha = tmp_norm / pAp
              IF (alpha <= 0.0_dp .OR. alpha > 1.0_dp / lAl) alpha = 1.0_dp / lAl
            END IF
          END IF

          IF (bs == 1) THEN
            WHERE (J)
              x = MAX(x - alpha * g, c)
            END WHERE
          ELSE
            WHERE (J)
              x = MIN(x - alpha * g, c)
            END WHERE
          END IF

          CALL MatrixVectorMultiply(MatA, x, g)
          g = g - b

          ! recompute z and p
          ! Apply diagonal preconditioning with free set masking
          IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
              == 'diagonal' ) THEN
            WHERE (J)
              WHERE (ABS(MatA % Values(MatA % Diag)) > AEPS)
                z = g / MatA % Values(MatA % Diag)
              ELSEWHERE
                z = g
              END WHERE
            ELSEWHERE
              z = 0.0_dp
            END WHERE
          ELSE
            ! No preconditioning, but still mask by free set
            WHERE (J)
              z = g
            ELSEWHERE
              z = 0.0_dp
            END WHERE
          END IF
          p = z

          ! update gf,gc,gr,gp
          gf = MERGE(g, 0.0_dp, J)

          IF (bs == 1) THEN
            gc = 0.0_dp
            WHERE (.NOT. J)
              gc = MIN(g, 0.0_dp)
            END WHERE
            gr = MERGE(MIN(lAl * (x - c), gf), 0.0_dp, J)
          ELSE
            gc = 0.0_dp
            WHERE (.NOT. J)
              gc = MAX(g, 0.0_dp)
            END WHERE
            gr = MERGE(MAX(lAl * (x - c), gf), 0.0_dp, J)
          END IF

          gp = gf + gc
          ne = ne + 1
        END IF
      ELSE
        ! proportioning step
        CALL MatrixVectorMultiply(MatA, gc, Ap)
        pAp = DOT_PRODUCT(gc(1:n), Ap(1:n))
        IF (ABS(pAp) < eps_local) THEN
          CALL Info('my_MPRGP','denominator in proportioning nearly zero, stopping',Level=5)
          EXIT
        END IF
        acg = DOT_PRODUCT(gc(1:n), g(1:n)) / pAp
        x = x - acg * gc
        IF (bs == 1) THEN
          DO i = 1, n
            IF (x(i) < c(i)) x(i) = c(i)
          END DO
        ELSE
          DO i = 1, n
            IF (x(i) > c(i)) x(i) = c(i)
          END DO
        END IF

        J = (bs * x > bs * c)
        g = g - acg * Ap

        ! Apply diagonal preconditioning with free set masking
        IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
            == 'diagonal' ) THEN
          WHERE (J)
            WHERE (ABS(MatA % Values(MatA % Diag)) > AEPS)
              z = g / MatA % Values(MatA % Diag)
            ELSEWHERE
              z = g
            END WHERE
          ELSEWHERE
            z = 0.0_dp
          END WHERE
        ELSE
          ! No preconditioning, but still mask by free set
          WHERE (J)
            z = g
          ELSEWHERE
            z = 0.0_dp
          END WHERE
        END IF
        p = z

        ! update gf,gc,gr,gp
        gf = MERGE(g, 0.0_dp, J)

        IF (bs == 1) THEN
          gc = 0.0_dp
          WHERE (.NOT. J)
            gc = MIN(g, 0.0_dp)
          END WHERE
          gr = MERGE(MIN(lAl * (x - c), gf), 0.0_dp, J)

        ELSE
          gc = 0.0_dp
          WHERE (.NOT. J)
            gc = MAX(g, 0.0_dp)
          END WHERE
          gr = MERGE(MAX(lAl * (x - c), gf), 0.0_dp, J)
        END IF

        gp = gf + gc
        np = np + 1
      END IF
    END DO  ! main loop

    final_norm_gp = ComputeNorm(CurrentModel % Solver, n, gp)
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
      CALL MatrixVectorMultiply(Aptr, v, w)
      normv = ComputeNorm(CurrentModel % Solver, nloc, w)
      IF (normv <= 0.0_dp) EXIT
      v = w / normv
    END DO

    CALL MatrixVectorMultiply(Aptr, v, w)
    out_norm = ComputeNorm(CurrentModel % Solver, nloc, w)

      DEALLOCATE(v, w)
    END SUBROUTINE estimate_matrix_norm

  END SUBROUTINE my_MPRGP
  

!------------------------------------------------------------------------------
SUBROUTINE MyMPRGPSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  ! Initialization function - following ModelPDEevol pattern
  CALL Info('MyMPRGPSolver_init','Initializing MPRGP solver',Level=5)
!------------------------------------------------------------------------------
END SUBROUTINE MyMPRGPSolver_init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE MyMPRGPSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  
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

  LowerLim = ListCheckPresentAnyBodyForce(Model,'Temp Lower Limit') .OR. &
      ListCheckPresentAnyBC(Model,'Temp Lower Limit') 
  UpperLim = ListCheckPresentAnyBodyForce(Model,'Temp Upper Limit') .OR. &
      ListCheckPresentAnyBC(Model,'Temp Upper Limit') 
  IF(LowerLim .AND. UpperLim) THEN
    CALL Warn('MyMPRGPSolver','This code cannot have both upper and lower limit at same time!')
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
    
    CALL my_MPRGP(n_mprgp, x, b, c, 1.0e-8_dp, 100000, 1.0_dp, .FALSE., &
          bound_type, ncg, ne, np, iters, converged_mprgp, final_norm_gp)
    
    CALL Info('MyMPRGPSolver','MPRGP finished: iters='//I2S(iters)// &
              ', ncg='//I2S(ncg)//', ne='//I2S(ne)//', np='//I2S(np), Level=5)

    IF (converged_mprgp .OR. final_norm_gp <= 1.0e-8_dp) THEN
      CALL Info('MyMPRGPSolver','Nonlinear solve: MPRGP converged, exiting outer loop', Level=5)
    END IF

    DEALLOCATE(c)
  END BLOCK

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
      Load(1:n) = GetReal( BodyForce,'Heat Source', Found )
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
      CALL Info('MyMPRGPSolver','Integration points in 1st element: '//I2S(IP % n),Level=8)
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
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight,Lim(n)
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

    IF ( ASSOCIATED(BC) ) THEN
      IF(UpperLim) THEN
        Lim(1:n) = GetReal(BC,'Temp Upper Limit',Found) 
      ELSE IF(LowerLim) THEN
        Lim(1:n) = GetReal(BC,'Temp Lower Limit',Found) 
      END IF
      IF(Found) THEN
        LimVal(Solver % Variable % Perm(Element % NodeIndexes)) = Lim(1:n)
      END IF
    END IF

    
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
END SUBROUTINE MyMPRGPSolver
!------------------------------------------------------------------------------
