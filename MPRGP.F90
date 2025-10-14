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
    
SUBROUTINE my_rpcond(x, r, J)
  REAL(KIND=dp) :: x(:), r(:)
  LOGICAL :: J(:)  ! free-set mask
  TYPE(Matrix_t), POINTER :: A
  LOGICAL :: Found
  INTEGER :: i
  
  IF( ListGetString( CurrentModel % Solver % Values,'Linear System Preconditioning', Found ) &
      == 'diagonal' ) THEN
    A => CurrentModel % Solver % Matrix
    WHERE (J)
      x = r / A % Values(A % Diag)
    ELSEWHERE
      x = 0.0_dp
    END WHERE
  ELSE
    ! No preconditioning, but still mask by free set
    WHERE (J)
      x = r
    ELSEWHERE
      x = 0.0_dp
    END WHERE
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
  
!------------------------------------------------------------------------------
  SUBROUTINE my_MPRGP(n, x, b, c, epsr, maxit, Gamma, adapt, bound, &
                      ncg, ne, np, iters, converged, norm_gp)
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
    REAL(KIND=dp), INTENT(OUT) :: norm_gp

    ! ---------------------------
    ! Local declarations (all here)
    ! ---------------------------
    REAL(KIND=dp), ALLOCATABLE :: g(:), gf(:), gc(:), gr(:), gp(:)
    REAL(KIND=dp), ALLOCATABLE :: z(:), p(:), Ap(:), yy(:), D(:), Agr(:)
    LOGICAL, ALLOCATABLE :: J(:)
    LOGICAL, ALLOCATABLE :: p_mask(:)
    INTEGER :: bs
    REAL(KIND=dp) :: lAl, alpha, a_f
    INTEGER :: i, allocstat
    REAL(KIND=dp) :: rtp, pAp, acg, beta
    REAL(KIND=dp) :: grg, grAgr
    REAL(KIND=dp) :: eps_local
    TYPE(Matrix_t), POINTER :: MatA
    REAL(KIND=dp) :: tol


    ! ---------------------------
    ! Allocate scratch vectors
    ! ---------------------------
    ALLOCATE(g(n), gf(n), gc(n), gr(n), gp(n), z(n), p(n), Ap(n), yy(n), J(n), Agr(n), STAT=allocstat)
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


    ! ---------------------------
    ! Initialization (g = A*x - b)
    ! ---------------------------
    CALL my_matvec(x, g)       ! g = A*x

    g = g - b

    tol = 1.0e-12_dp
    IF (bs == 1) THEN
      J = (x > c + tol) ! J is free set
    ELSE
      J = (x < c - tol)
    END IF
    gf = MERGE(g, 0.0_dp, J)

    IF (bs == 1) THEN
      WHERE (.NOT. J)
        gc = MIN(g, 0.0_dp)
      ELSEWHERE
        gc = 0.0_dp
      END WHERE
    ELSE
      WHERE (.NOT. J)  ! WHERE COMMAND FOR vECTORS
        gc = MAX(g, 0.0_dp)
      ELSEWHERE
        gc = 0.0_dp
      END WHERE
    END IF

    ! estimate matrix norm ||A|| with a small power iteration
    CALL estimate_matrix_norm(n, 10, lAl)
    IF (lAl <= 0.0_dp) lAl = 1.0_dp
    alpha = 1.0_dp / lAl

    ! reduced free gradient gr
    IF (bs == 1) THEN
      gr = MERGE(MIN(lAl * (x - c), gf), 0.0_dp, J) ! this might be wrong,
    ELSE
      gr = MERGE(MAX(lAl * (x - c), gf), 0.0_dp, J)
    END IF

    gp = gf + gc

    ! TODO - write this as a function
    ! preconditioning: z = M^{-1} * g on free set
    CALL my_rpcond(z, g, J)

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
      !WRITE(*,*) "iteration", iters
      !WRITE(*,*) "Norm of gp is:", my_normfun(n, gp)

      IF ( my_dotprodfun(n, gc, gc) <= (Gamma**2) * my_dotprodfun(n, gr, gf) ) THEN
        ! CG-like step
        CALL my_matvec(p, Ap)
        rtp = my_dotprodfun(n, z, g) ! residual * p
        pAp = my_dotprodfun(n, p, Ap)

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
          CALL my_rpcond(z, g, J)

          beta = my_dotprodfun(n, z, Ap) / pAp
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
          p_mask = (bs * p > 0.0_dp) .AND. J ! indexes where p is moving towards the bound
          a_f = MINVAL((x-c) / p,p_mask)
          
          IF (a_f < 0.0_dp) a_f = 0.0_dp
          ! halfstep
          IF (bs == 1) THEN
              x = MAX(x - a_f * p, c)
          ELSE
              x = MIN(x - a_f * p, c)
          END IF

          J = (bs * x > bs * c)
          g = g - a_f * Ap

          ! adaptive alpha
          IF (adapt) THEN
            CALL my_matvec(gr, Agr)
            grg = my_dotprodfun(n, gr, g)
            grAgr = my_dotprodfun(n, gr, Agr)
            IF (ABS(grAgr) < eps_local) THEN
              alpha = 1.0_dp / lAl
            ELSE
              alpha = grg / grAgr
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

          CALL my_matvec(x, g) !calculate Ax, save it in g
          g = g - b

          ! recompute z and p
          CALL my_rpcond(z, g, J)
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
        CALL my_matvec(gc, Ap)
        pAp = my_dotprodfun(n, gc, Ap)
        IF (ABS(pAp) < eps_local) THEN
          CALL Info('my_MPRGP','denominator in proportioning nearly zero, stopping',Level=5)
          EXIT
        END IF
        acg = my_dotprodfun(n, gc, g) / pAp

        x = x - acg * gc

        ! safeguard if we end up outside the bound
        IF (bs == 1) THEN
          x = MAX(x, c)
        ELSE
          x = MIN(x, c)
        END IF

        J = (bs * x > bs * c)
        g = g - acg * Ap

        CALL my_rpcond(z, g, J)
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

    norm_gp = my_normfun(n, gp)
    converged = (norm_gp <= epsr)

    ! cleanup
    IF (ALLOCATED(D)) DEALLOCATE(D)
    DEALLOCATE(g, gf, gc, gr, gp, z, p, Ap, yy, J)

    RETURN
  CONTAINS
    ! power-method norm estimator using A
    SUBROUTINE estimate_matrix_norm(nloc, niter, out_norm)
      INTEGER, INTENT(IN) :: nloc, niter
      REAL(KIND=dp), INTENT(OUT) :: out_norm
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp), ALLOCATABLE :: v(:), w(:)
      INTEGER :: itl, allocstat
      REAL(KIND=dp) :: normv

      ALLOCATE(v(nloc), w(nloc), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('estimate_matrix_norm','alloc fail')
      END IF

      ! Get matrix pointer from CurrentModel like my_matvec does
      A => CurrentModel % Solver % Matrix

      DO itl = 1, nloc
        v(itl) = 1.0_dp / REAL(nloc, dp)
      END DO

      DO itl = 1, niter
        CALL CRS_MatrixVectorMultiply(A, v, w)
        normv = my_normfun(nloc, w)
        IF (normv <= 0.0_dp) EXIT
        v = w / normv
      END DO

      CALL CRS_MatrixVectorMultiply(A, v, w)
      out_norm = my_normfun(nloc, w)

      DEALLOCATE(v, w)
    END SUBROUTINE estimate_matrix_norm

  END SUBROUTINE my_MPRGP
  
END MODULE MyLinearSolver
  

!------------------------------------------------------------------------------
SUBROUTINE MPRGPSolver( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=50) :: filename
  TYPE(Matrix_t), POINTER :: A
  
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
    CALL Warn('MPRGPSolver','This code cannot have both upper and lower limit at same time!')
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

    A => Solver % Matrix

    WRITE(filename, '(A,I0,A)') 'a_mprgp_iter.dat'
    OPEN(1,FILE=TRIM(filename), STATUS='Unknown')
    CALL PrintMatrix(A,.FALSE.,.FALSE.)
    CLOSE(1)

    
    

    ! And finally, solve:
    !--------------------

    ! Use MPRGP for bound-constrained problems
    BLOCK
      INTEGER :: n_mprgp, ncg, ne, np, iters
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp), POINTER :: x(:), b(:)
      REAL(KIND=dp), ALLOCATABLE :: c(:)
      LOGICAL :: converged_mprgp
      REAL(KIND=dp) :: norm_gp
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
      
      ! Set initial guess: u = max(0, c) for lower bounds, u = min(0, c) for upper bounds
      IF(UpperLim) THEN
        ! For upper bounds: u = min(0, c) - start below the upper bound
        x = MIN(0.0_dp, c)
        write(*,*) "Debug: Initial guess x for upper bounds"
      ELSE IF(LowerLim) THEN
        ! For lower bounds: u = max(0, c) - start above the lower bound  
        x = MAX(0.0_dp, c)
        write(*,*) "Debug: Initial guess x for lower bounds"
      ELSE
        ! No bounds, keep default (zeros)
        x = 0.0_dp
        write(*,*) "Debug: Initial guess x for no bounds"
      END IF
      write(*,*) "Debug: Initial guess x =", my_normfun(n_mprgp, x)
      
      CALL my_MPRGP(n_mprgp, x, b, c, 1.0e-8_dp, 100000, 1.0_dp, .FALSE., &
            bound_type, ncg, ne, np, iters, converged_mprgp, norm_gp)
      
      CALL Info('MPRGPSolver','MPRGP finished: iters='//I2S(iters)// &
                ', ncg='//I2S(ncg)//', ne='//I2S(ne)//', np='//I2S(np), Level=5)
      ! Save the x
      OPEN(1, FILE="x_mprgp.dat", STATUS='Unknown')
        DO i=1,SIZE(Solver % Variable % Values)
          WRITE(1,*) Solver % Variable % Values(i)    
        END DO
        CLOSE(1)

      !IF (converged_mprgp .OR. norm_gp <= 1.0e-8_dp) THEN
      !  CALL Info('MPRGPSolver','Nonlinear solve: MPRGP converged, exiting outer loop', Level=5)
      !  EXIT   ! breaks DO iter=1,maxiter
      !END IF

      OPEN(1,FILE="x_mprgp.dat", STATUS='Unknown')
      DO i=1,SIZE(x)
        WRITE(1,*) x(i)    
      END DO
      CLOSE(1)

      DEALLOCATE(c)
    END BLOCK

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
      CALL Info('MPRGPSolver','Integration points in 1st element: '//I2S(IP % n),Level=8)
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
END SUBROUTINE MPRGPSolver
!------------------------------------------------------------------------------