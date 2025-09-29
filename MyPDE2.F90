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
      BLOCK
        INTEGER :: n
        TYPE(Matrix_t), POINTER :: A
        REAL(KIND=dp), POINTER :: x(:), b(:)
        INTEGER :: Rounds, MinIter, OutputInterval, m
        REAL(KIND=dp) :: MinTolerance, MaxTolerance
        LOGICAL :: Converged, Diverged
        REAL(KIND=dp) :: Residual
        
        n = SIZE(Solver % Variable % Values)
        A => Solver % Matrix
        b => Solver % Matrix % Rhs
        x => Solver % Variable % Values

        Rounds = ListGetInteger( Params,'Linear System Max Iterations',Found )
        MinIter = ListGetInteger( Params,'Linear System Min Iterations',Found )
        MinTolerance = ListGetCReal( Params,'Linear System Convergence Tolerance',Found )
        MaxTolerance = 1.0e20
        OutputInterval = ListGetInteger( Params,'Linear System Residual Output', Found )               
        m = ListGetInteger( Params,'Linear System GCR Restart', Found )
        IF(.NOT. Found ) m = MIN(Rounds, 200)
        
        CALL my_GCR( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
            Converged, Diverged, OutputInterval, m, MinIter) 
        
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
