!------------------------------------------------------------------------------
SUBROUTINE DefaultSolver( Model,Solver,dt,TransientSimulation )
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

  LowerLim = ListCheckPresentAnyBodyForce(Model,'Temp Lower Limit')
  UpperLim = ListCheckPresentAnyBodyForce(Model,'Temp Upper Limit')
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

    
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    WRITE(*,*) "Iteration", iter

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
    IF(iter == 1) THEN
    WRITE(filename, '(A,I0,A)') 'a_default_iter', iter, '.dat'
      OPEN(1,FILE=TRIM(filename), STATUS='Unknown')
      CALL PrintMatrix(Solver % Matrix,.FALSE.,.FALSE.)
      CLOSE(1)
    END IF

    IF(iter == 2) THEN
      OPEN(1, FILE="x_default.dat", STATUS='Unknown')
      DO i=1,SIZE(Solver % Variable % Values)
        WRITE(1,*) Solver % Variable % Values(i)    
      END DO
      CLOSE(1)

    END IF

   
    Norm = DefaultSolve()      
    IF( DefaultConverged() ) EXIT    


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
END SUBROUTINE DefaultSolver
!------------------------------------------------------------------------------