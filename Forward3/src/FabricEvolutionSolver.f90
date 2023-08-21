 SUBROUTINE FabricEvolutionSolver( Model,Solver,dt,TransientSimulation )
   !------------------------------------------------------------------------------
   !******************************************************************************
   !  Calculate 2nd order orientation tensor evolution. SemiLagrangian Method
   !
   !  Notes: 
   !  Reference to Martin and Gudmundsson (2012), Effects on nonlinear rheology and anisotropy on the
   !  relationship between age and depth at ice divides,TC, would be appreciated.  
   !
   !  Numerical details:
   !  This solver implement a semilagrangian algoritm with a two-time-level scheme and linear interpolation. 
   !  It is based in Staniforth, Andrew, Jean Côté, 1991: Semi-Lagrangian Integration Schemes for Atmospheric 
   !  Models—A Review. Mon. Wea. Rev., 119, 2206–2223.
   ! 
   !  Changes made by Clara Henry
   !******************************************************************************
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
   Integer :: DIM,N,e,i,j,k,kf,ka,kd,kd2,kk,NMax,NM,stat
   Logical :: GotIt,FirstTime=.TRUE.,IsThere

   INTEGER, POINTER :: a2Perm(:)
 !!!!!  REAL(KIND=dp), POINTER :: a2(:),a20(:)
   Integer :: DOFs

   TYPE(Element_t),POINTER :: Element
   TYPE(Nodes_t)   :: ElementNodes
   INTEGER,  POINTER ::  NodeIndexes(:)
   TYPE(Matrix_t), POINTER :: Systemmatrix

   CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolName
   TYPE(Variable_t), POINTER :: FlowSol
   INTEGER, POINTER :: FlowPerm(:)
   REAL(KIND=dp), POINTER :: Flow(:)
   INTEGER :: FlowDOFs

   INTEGER :: NNMAX
   INTEGER,  ALLOCATABLE :: NoNeigh(:),NeighList(:,:)
   REAL(KIND=dp) :: xa(3),xm(3),xd(3),alpha(3),um(3)
   REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
   REAL(KIND=dp) :: eps=1.0e-6
   REAL(KIND=dp), ALLOCATABLE :: Vector(:,:)
   Integer :: iv,ip,jp,kp,lp

   REAL(KIND=dp) :: AlphaFabric,gama,beta,lambda1,lambda2,EU

   REAL(KIND=dp) :: Rho, Iota

   REAL(KIND=dp) :: fa,fd
   INTEGER :: a4DOFs,etaDOFs,DDOFs,WDOFs,MaxDOFs
   REAL(KIND=dp), ALLOCATABLE :: a2d(:),a2(:), a2a(:),a4(:), a2o(:),&
         Angle(:),FabricGrid(:)&
        ,eta(:,:),etaTemp(:,:),DTay(:),DU(:),D(:),W(:),F(:)
   Integer :: Swap(6)
   CHARACTER(LEN=MAX_NAME_LEN) :: viscosityFile

   !   REAL(KIND=dp) :: at,st,totat,totst,CPUTime,Norm,PrevNorm,RelativeChange
      REAL(KIND=dp) :: at,st,totat,totst,Norm,PrevNorm,RelativeChange

   !Experimental
   LOGICAL, ALLOCATABLE :: Found(:),isBC(:),isa2In(:)
   TYPE(ValueList_t), POINTER :: BC
   REAL(KIND=dp), ALLOCATABLE :: Cond(:)

   INTEGER :: nn,precv
   INTEGER, POINTER :: nlist(:)

   INTEGER :: ierr,gk
   INTEGER :: request(ParEnv % PEs)
   TYPE buffer_t
      INTEGER :: n
      INTEGER, ALLOCATABLE :: gbuff(:)
      REAL(KIND=dp), ALLOCATABLE :: vbuff(:)
   END TYPE buffer_t
   TYPE(buffer_t) :: RequestSend(ParEnv % PEs),RequestRecv(ParEnv % PEs)
   TYPE(buffer_t) :: ReplySend(ParEnv % PEs),ReplyRecv(ParEnv % PEs)

   INTERFACE
      SUBROUTINE IBOF(ai,a4)
        USE Types
        REAL(KIND=dp),intent(in) :: ai(6)
        REAL(KIND=dp),intent(out) :: a4(9)
      END SUBROUTINE IBOF

      Subroutine R2Ro(ai,dim,a2,angle)
        USE Types
        REAL(KIND=dp),intent(in) :: ai(6)
        Integer :: dim
        REAL(KIND=dp),intent(out) :: a2(3), Angle(3)
      End Subroutine R2Ro

      Subroutine OPILGGE_ai_nl(a2,Angle,etaI,eta36)
        USE Types
        REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: a2
        REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: Angle
        REAL(kind=dp), INTENT(in),  DIMENSION(:)   :: etaI
        REAL(kind=dp), INTENT(out), DIMENSION(6,6) :: eta36
      END SUBROUTINE OPILGGE_ai_nl

   END INTERFACE

   SAVE FirstTime,NoNeigh,NeighList,Found,isBC,isa2In

   WRITE(Message,'(a)') 'Start Solver'
   CALL Info('FabricEvolutionSolver', Message, Level=4)

   !------------------------------------------------------------------------------
   ! Get Constants
   !------------------------------------------------------------------------------

   DIM = CoordinateSystemDimension()
   NMAX = Model % MaxElementNodes
   NM = Solver % Mesh % NumberOfNodes

   NNMAX=20!8
   a4DOFs=9
   etaDOFs=36
   DDOFs=6
   WDOFs=3
   MaxDOFs=36


   AlphaFabric = ListGetConstReal(Solver % Values, 'Alpha',GotIt)
   IF(.NOT. GotIt) CALL Fatal( 'FabricEvolutionSolver','Could find value of Alpha, Aborting.' ) 
   Iota = ListGetConstReal(Solver % Values, 'Iota',GotIt)
   IF(.NOT. GotIt) CALL Fatal( 'FabricEvolutionSolver','Could find value of Iota, Aborting.' )
   Rho = ListGetConstReal(Solver % Values, 'Alpha',GotIt) !Swapping to Rho as Fabien
   IF(.NOT. GotIt) CALL Fatal( 'FabricEvolutionSolver','Could find value of Alpha, Aborting.' ) 

   lambda1=2*((gama+2)/(4*gama-1)*beta-1.0)
   lambda2=(1.0-beta)
!   EU = 1.0_dp/(2._dp/5._dp+beta/5.0_dp*(2._dp+3._dp/(4._dp*gama-1._dp)));
   EU = 1.0d0
 !------------------------------------------------------------------------------
  !    Get variables for the solution
  !------------------------------------------------------------------------------

!!!!  a2     => Solver % Variable % Values     ! Nodal values for 
  a2Perm => Solver % Variable % Perm       ! Permutations for 
!!!!  a20 => Solver % Variable % PrevValues(:,1)
  DOFs = Solver % Variable % DOFs

  FlowSolName =  GetString( GetEquation(Solver % Mesh % Elements(1)),'Flow Solution Name', GotIt)
  IF(.NOT.GotIt) THEN        
     CALL WARN('FabricEvolutionSolver','Keyword >Flow Solution Name< not found in section >Equation<')
     CALL WARN('FabricEvolutionSolver','Taking default value >Flow Solution<')
     WRITE(FlowSolName,'(A)') 'Flow Solution'
  END IF
  FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
  IF ( ASSOCIATED( FlowSol ) ) THEN
     FlowPerm     => FlowSol % Perm
     FlowDOFs     =  FlowSol % DOFs
     Flow               => FlowSol % Values
  ELSE
     WRITE(Message,'(A,A,A)') &
          'Convection flag set to >computed<, but no variable >',FlowSolName,'< found'
     CALL FATAL('FabricEvolutionSolver',Message)              
  END IF

  !------------------------------------------------------------------------------
  !    Inicialization
  !-----------------------------------------------------------------------------

  ALLOCATE( ElementNodes % x( NMAX ),      &
       ElementNodes % y( NMAX ),      &
       ElementNodes % z( NMAX ),      &
       Vector(MaxDOFS,NMAX),                  &
       a2d(DOFs),                                 &
       a2(DOFs),                                 &
       a2a(DOFs),                                 &
       a4(a4DOFs),                                 &
       a2o(DOFs),                                 &
       Angle(DOFs),                                 &
       eta(DDOFs,DDOFs),                               &
       etaTemp(DDOFs,DDOFs),                               &
       DTay(DDOFs),                                 &
       DU(DDOFs),                                 &
       D(DDOFs),                                 &
       W(WDOFs),                                 &
       F(DOFs),                                  &
       STAT=stat  )
  IF ( stat /= 0 ) THEN
     CALL Fatal('FabricEvolutionSolver','Memory allocation error, Aborting.')
  END IF

  IF(FirstTime) THEN
     ALLOCATE( NoNeigh(NM),&
          Neighlist(NM,NNMAX),&
          Found(NM),                  &
          Cond(NMAX),                  &
          FabricGrid(4879),                                 &
          STAT=stat  )
     IF ( stat /= 0 ) THEN
        CALL Fatal('FabricEvolutionSolver','Memory allocation error, Aborting.')
     END IF


     NoNeigh=0
     DO e=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(e)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        DO i=1,n
           k=NodeIndexes(i)

           IF(NoNeigh(k)==0) THEN
              NoNeigh(k)=1
              NeighList(k,1)=e
           ELSE
              IsThere=.FALSE.
              DO j=1,NoNeigh(k)
                 IF(NeighList(k,j)==e) THEN
                    IsThere=.TRUE.
                    EXIT
                 END IF
              END DO

              IF(.NOT.IsThere) THEN
                 NoNeigh(k)= NoNeigh(k) + 1
                 NeighList(k, NoNeigh(k))=e                 
              END IF
           END IF
        END DO

     END DO

     ! Flag nodes that are in BC with a2In. 
     !(We will assume that the horizontal gradient of age is null and ignore the horizontal components of velocity.) 
     ALLOCATE(isa2In(NM),                  &
          STAT=stat  )
     IF ( stat /= 0 ) THEN
        CALL Fatal('FabricEvolutionSolver','Memory allocation error, Aborting.')
     END IF
     isa2In=.FALSE.
     DO i=1,GetNofBoundaryElements()
        Element=>GetBoundaryElement(i)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
        BC=>GetBC()
        IF(.NOT.ASSOCIATED(BC)) CYCLE

        IF(ListCheckPresent(BC,'a2In'))  THEN

           DO j=1,Model % NumberOfBCs
              IF ( Element % BoundaryInfo % Constraint == Model % BCs(j) % Tag ) THEN                  
                 !                print *,ListGetLogical(Model % BCs(j) % Values,'a2In', gotIt )
                 isa2In(Nodeindexes)= &
                      ListGetLogical(Model % BCs(j) % Values,'a2In', gotIt )
              END IF
           END DO

        END IF

     END DO




     IF(ParEnv % PEs>1) THEN


        ! Flag nodes that are Dirichlet Boundary conditions

        ALLOCATE(isBC(NM),                  &
             STAT=stat  )
        IF ( stat /= 0 ) THEN
           CALL Fatal('FabricEvolutionSolver','Memory allocation error, Aborting.')
        END IF


        ! Check if the node has Dirichlet BC, in that case we will ignore
        isBC=.FALSE.
        DO i=1,GetNofBoundaryElements()
           Element=>GetBoundaryElement(i)
           n = GetElementNOFNodes()
           NodeIndexes => Element % NodeIndexes
           BC=>GetBC()
           IF(.NOT.ASSOCIATED(BC)) CYCLE

           IF(ListCheckPresent(BC,Solver % Variable % Name))  THEN

              ! Check first if we are using a2 Condition = -1
              IF(ListCheckPresent(BC,Trim(Solver % Variable % Name)//' Condition'))  THEN

                 DO j=1,Model % NumberOfBCs
                    IF ( Element % BoundaryInfo % Constraint == Model % BCs(j) % Tag ) THEN                  
                       isBC(Nodeindexes)= &
                            (ListGetReal(Model % BCs(j) % Values,&
                            Trim(Solver % Variable % Name)//' Condition',n, NodeIndexes, gotIt )>=0.0d0)
                    END IF
                 END DO
              ELSE
                 isBC(Nodeindexes)=.TRUE.
              END IF

           END IF

        END DO
     END IF

     !Read FabricGrid
     viscosityFile = ListGetString(Solver % Values, 'Viscosity File',GotIt)
     IF(.NOT. GotIt) CALL Fatal( 'FabricEvolutionSolver','Could find name of Viscosity File, Aborting.' ) 
     OPEN( 1, File = viscosityFile)
     DO i=1,813
        READ( 1, '(6(e14.8))' ) FabricGrid( 6*(i-1)+1:6*(i-1)+6 )
     END DO
     READ(1 , '(e14.8)' ) FabricGrid(4879)
     CLOSE(1)

  END IF

  SystemMatrix => Solver % Matrix

  totat = 0.0_dp
  totst = 0.0_dp
  at = CPUTime()
  CALL INFO( 'FabricEvolutionSolver', 'start assembly', Level=4 )

  CALL DefaultInitialize()
  CALL DefaultFinishAssembly()

  Norm=0.0d0

  DO k=1,NM
     ka=a2Perm(k)
     kf=FlowPerm(k)

     xa(1) = Solver % Mesh % Nodes % x(k)
     xa(2) = Solver % Mesh % Nodes % y(k)
     xa(3) = Solver % Mesh % Nodes % z(k) 

     IF(DIM==2.AND.FlowDOFs==3) THEN
        alpha(1)=Dt*Flow(FlowDOFs*kf-2)
        alpha(2)=Dt*Flow(FlowDOFs*kf-1)
        alpha(3)=0._dp
     ELSE IF(DIM==3.AND.FlowDOFs==4) THEN
        alpha(1)=Dt*Flow(FlowDOFs*kf-3)
        alpha(2)=Dt*Flow(FlowDOFs*kf-2)
        alpha(3)=Dt*Flow(FlowDOFs*kf-1)
     ELSE
        CALL Fatal('FabricEvolutionSolver','DIM AND FlowDOFS do not combine.  Aborting.')
     END IF

     IF(isa2In(k)) THEN
        DO j=1,DIM-1
           alpha(j)=0.0d0
        END DO
     END IF

     xm(1:3)=xa(1:3)-alpha(1:3)/2._dp


     Found(k)=.FALSE.  
     IsThere=.FALSE.
     DO i=1,NoNeigh(k)
        e=NeighList(k,i)
        Element => Solver % Mesh % Elements(e)

        n = Element % Type % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

        IF ( PointInElement( Element, ElementNodes, xm, LocalCoordinates,NumericEps=eps) )  THEN
           IsThere=.TRUE.
           EXIT
        END IF
        
     END DO

     IF ( IsThere) THEN
        

        DO j=1,DIM
           Vector(j,1:n)=Flow(FlowDOFs*(FlowPerm(NodeIndexes)-1)+j)
        END DO

        IF(isa2In(k)) THEN
           DO j=1,DIM-1
              Vector(j,1:n)=0.0d0
           END DO
        END IF


        !velocity at the middle point
        um=0.0d0
        DO j=1,DIM
           um(j)=InterpolateInElement(Element,Vector(j,1:n),&
                LocalCoordinates(1),LocalCoordinates(2), LocalCoordinates(3))
        END DO
        alpha(1)=Dt*um(1)
        alpha(2)=Dt*um(2)
        alpha(3)=Dt*um(3)           

        ! a2 at the middle point, remember is axx, ayy, azz, axy, axz, ayz
        CALL GetVectorLocalSolution(Vector, 'a2',Element)
        DO j=1,DOFs
           a2(j)=InterpolateInElement(Element,Vector(j,1:n),&
                LocalCoordinates(1),LocalCoordinates(2), LocalCoordinates(3))
        END DO

        ! a4 at the middle point (IBOF closure)  
!!! (a2 rentre dans l'ordre : 11, 22, 33, 12, 23 ,13)
        CALL IBOF((/a2(1),a2(2),a2(3),a2(4),a2(6),a2(5)/),a4)

        !D at the middle point
        CALL GetVectorLocalSolution(Vector, 'VG',Element)
        DO j=1,DDOFs
           DTay(j)=InterpolateInElement(Element,Vector(j,1:n),&
                LocalCoordinates(1),LocalCoordinates(2), LocalCoordinates(3))
        END DO

        !W at the middle point
        DO j=1,WDOFs
           W(j)=InterpolateInElement(Element,Vector(DDOFs+j,1:n),&
                LocalCoordinates(1),LocalCoordinates(2), LocalCoordinates(3))
        END DO

        DU=0.0d0
        IF(Rho>0.0) THEN
           !Getting viscosity at the middle point

           !    fourth order orientation tensor
           ! a4 at the middle point (IBOF closure)  
!!! (a2 rentre dans l'ordre : 11, 22, 33, 12, 23 ,13)
           Swap(1:6)=(/1,2,3,4,6,5/)
           CALL IBOF(a2(Swap(1:6)),a4)

           !     A2 expressed in the orthotropic frame
           call R2Ro(a2(Swap(1:6)),DIM,a2o,Angle)

           !     Get viscosity in Fabien order
           CALL OPILGGE_ai_nl(a2o, Angle, FabricGrid, etaTemp)

           !Swap Fabien index
           Do ip=1,DDOFs
              Do jp=1,DDOFs
                 eta(ip,jp)=etaTemp(Swap(ip),Swap(jp))
              End Do
           EndDo
           ! Also, removing the 2 in the diagonal terms of the compressive part
           Do ip=1,3
              eta(ip,ip)=eta(ip,ip)/2.0;
           END Do

           DU=MATMUL(eta,DTay)
        END IF

        !Alpha Fabric, rho here,  weights between uniform (1) and Taylor (0)
        DO j=1,DDOFs
           D(j)=DTay(j)*(1.0-Rho)+DU(j)*Rho
 !          PRINT *,j,DTay(j),DU(j),D(j)
        END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 
        !  Calculating the forcing in the middle point
        !
        !Algebra thanks to Matlab algebra/da_dt2D.m algebra/da_dt3D.m
        !translation to Elmerish through 
        ! algebra/Forcing2D.f90 and algebra/Forcing3D.f90   
        !

        WRITE( Message, * ) 'Alpha:   : ',Rho
        CALL Info( 'FabricEvolutionSolver', Message, Level=3 )
        WRITE( Message, * ) 'Iota:   : ',Iota
        CALL Info( 'FabricEvolutionSolver', Message, Level=3 )

        ! Fxx
        F(1)=(2*a4(1)*D(1) - 2*a2(1)*D(1) + 4*a4(7)*D(4) - 2*a2(4)*D(4) + &
             4*a4(6)*D(5) - 2*a2(5)*D(5) + 2*a4(3)*D(2) + 4*a4(4)*D(6)) * Iota + &
             2*a2(4)*W(1) + 2*a2(5)*W(2) + (- 2*D(3)*(a4(1) - a2(1) + a4(3))) * Iota

        ! Fyy
        F(2)=(2*a4(3)*D(1) - 2*a2(4)*D(4) + 4*a4(9)*D(4) + 4*a4(5)*D(5) - &
             2*a2(2)*D(2) + 2*a4(2)*D(2) + 4*a4(8)*D(6) - 2*a2(6)*D(6)) * Iota - &
             2*a2(4)*W(1) + 2*a2(6)*W(3) + (- 2*D(3)*(a4(3) - a2(2) + a4(2))) * Iota

           ! Fzz
        F(3)=(2*a2(5)*D(5) - 4*a4(6)*D(5) - 4*a4(4)*D(6) - 4*a4(5)*D(5) - &
             2*a2(1)*D(3) + 2*a4(1)*D(3) + 4*a4(3)*D(3) - 4*a4(8)*D(6) + &
             2*a2(6)*D(6) - 2*a2(2)*D(3) + 2*a4(2)*D(3)) * Iota - 2*a2(5)*W(2) - &
             2*a2(6)*W(3) + (- 2*D(1)*(a4(1) - a2(1) + a4(3)) - 4*D(4)*(a4(7) - &
             a2(4) + a4(9)) - 2*D(2)*(a4(3) - a2(2) + a4(2))) * Iota

           ! Fxy
        F(4)=(2*a4(7)*D(1) - a2(1)*D(4) - a2(4)*D(1) + 4*a4(3)*D(4) + &
             4*a4(4)*D(5) - a2(4)*D(2) - a2(2)*D(4) - a2(5)*D(6) - a2(6)*D(5) + &
             2*a4(9)*D(2) + 4*a4(5)*D(6)) * Iota - a2(1)*W(1) + a2(2)*W(1) + &
             a2(5)*W(3) + a2(6)*W(2) + (- 2*D(3)*(a4(7) - a2(4) + a4(9))) * Iota

        ! Fxz
        F(5)=(2*a4(6)*D(1) + 4*a4(4)*D(4) + 3*a2(1)*D(5) - a2(5)*D(1) - &
             4*a4(1)*D(5) - 4*a4(3)*D(5) - 4*a4(7)*D(6) + 3*a2(4)*D(6) - &
             a2(6)*D(4) - 2*a4(6)*D(3) + 2*a4(5)*D(2) - 4*a4(9)*D(6) + &
             a2(5)*D(3) - a2(3)*D(5) - 2*a4(5)*D(3)) * Iota - a2(1)*W(2) - &
             a2(4)*W(3) + a2(6)*W(1) + a2(3)*W(2)

        ! Fyz
        F(6)=(2*a4(4)*D(1) - 4*a4(7)*D(5) + 3*a2(4)*D(5) - a2(5)*D(4) - &
             4*a4(3)*D(6) + 4*a4(5)*D(4) - 4*a4(9)*D(5) - 2*a4(4)*D(3) + &
             2*a4(8)*D(2) + 3*a2(2)*D(6) - a2(6)*D(2) - 4*a4(2)*D(6) - &
             2*a4(8)*D(3) + a2(6)*D(3) - a2(3)*D(6)) * Iota - a2(4)*W(2) - &
             a2(5)*W(1) - a2(2)*W(3) + a2(3)*W(3)

        xd(1)=xa(1)-alpha(1)
        xd(2)=xa(2)-alpha(2)
        xd(3)=xa(3)-alpha(3)

        IsThere=.FALSE.
        DO i=1,NoNeigh(k)
           e=NeighList(k,i)
           Element => Solver % Mesh % Elements(e)

           n = Element % Type % NumberOfNodes
           NodeIndexes => Element % NodeIndexes

           ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

           IF ( PointInElement( Element, ElementNodes, xd, LocalCoordinates,NumericEps=eps) )  THEN
              IsThere=.TRUE.
              EXIT
           END IF
        END DO

        IF ( IsThere) THEN

           Found(k)=.TRUE.

           ! a2 at the Departure point
           CALL GetVectorLocalSolution(Vector, 'a2',Element)     
           DO j=1,DOFs
              a2d(j)=InterpolateInElement(Element,Vector(j,1:n),&
                   LocalCoordinates(1),LocalCoordinates(2), LocalCoordinates(3))
           END DO

           ! a2 at the arrival point
           DO j=1,DOFs
              a2a(j)=a2d(j)+Dt*F(j)
           END DO

           !Minimal Filtering
!!$           DO j=1,3
!!$              a2a(j) = max(0d0,a2a(j))
!!$              a2a(j) = min(1d0,a2a(j))
!!$           END DO
!!$           DO j=4,DOFs
!!$              a2a(j) = max(-1.0/4.0,a2a(j))
!!$              a2a(j) = min(  1.0/4.0,a2a(j))
!!$           END DO

           !Filtering
           fd=1._dp-27._dp*(a2d(1)*a2d(2)*a2d(3)+2*a2d(4)*a2d(5)*a2d(6)&
                -a2d(1)*a2d(6)*a2d(6)-a2d(2)*a2d(5)*a2d(5)-a2d(3)*a2d(4)*a2d(4))

           fa=1._dp-27._dp*(a2a(1)*a2a(2)*a2a(3)+2*a2a(4)*a2a(5)*a2a(6)&
                -a2a(1)*a2a(6)*a2a(6)-a2a(2)*a2a(5)*a2a(5)-a2a(3)*a2a(4)*a2a(4))

           fa=max(0._dp,fa)
           fa=min(1._dp,fa)

           IF(fa<fd) THEN
              DO j=1,DOFs
                 a2a(j)=a2d(j)
              END DO
           END IF


           !FINALLY we write the matrix for parallel computations
           DO j=1,DOFs
              CALL SetMatrixElement( SystemMatrix,(ka-1)*DOFs+j,(ka-1)*DOFs+j, 1.0d0)
              SystemMatrix % RHS((ka-1)*DOFs+j) =  a2a(j)
           END DO

        END IF

     END IF

  END DO

  CALL DefaultDirichletBCs()

  !------------------------------------------------------------------------------
  !    Solve System  and check for convergence
  !------------------------------------------------------------------------------
  at = CPUTime() - at
  st = CPUTime() 

  PrevNorm = Solver % Variable % Norm

  Norm = DefaultSolve()

  IF ( PrevNorm + Norm /= 0.0_dp ) THEN
     RelativeChange = 2.0_dp * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
  ELSE
     RelativeChange = 0.0_dp
  END IF

  WRITE( Message, * ) 'Result Norm   : ',Norm
  CALL INFO( 'FabricEvolutionSolver', Message, Level=4 )
  WRITE( Message, * ) 'Relative Change : ',RelativeChange
  CALL INFO( 'FabricEvolutionSolver', Message, Level=4 )

  FirstTime=.FALSE.


END SUBROUTINE FabricEvolutionSolver




!!!
!!! Fabian 2006
!!!    
!!! calcul de a4 a partir de a2 par fermeture IBOF (Chung 2002)
!!! a2 rentre dans l'ordre : 11, 22, 33, 12, 23 ,13
!!! a4  sort dans l'ordre : 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
!!!
!------------------------------------------------------------------------------
subroutine IBOF(a2,a4)

  USE Types

  implicit none
  Real(dp),dimension(6),intent(in):: a2  
  Real(dp),dimension(9),intent(out):: a4  
  Real(dp):: a_11,a_22,a_33,a_12,a_13,a_23
  Real(dp):: b_11,b_22,b_12,b_13,b_23
  Real(dp):: aPlusa

  Real(dp),dimension(21) :: vec
  Real(dp),dimension(3,21) :: Mat
  Real(dp),dimension(6) :: beta
  Real(dp) :: Inv2,Inv3
  integer :: i,j      


  a_11=a2(1)
  a_22=a2(2)
  a_33=a2(3)
  a_12=a2(4)
  a_23=a2(5)
  a_13=a2(6)

  !Coefficients 

  Mat(1,1)=0.217774509809788e+02_dp
  Mat(1,2)=-.297570854171128e+03_dp
  Mat(1,3)=0.188686077307885e+04_dp
  Mat(1,4)=-.272941724578513e+03_dp
  Mat(1,5)=0.417148493642195e+03_dp
  Mat(1,6)=0.152038182241196e+04_dp
  Mat(1,7)=-.137643852992708e+04_dp
  Mat(1,8)=-.628895857556395e+03_dp
  Mat(1,9)=-.526081007711996e+04_dp
  Mat(1,10)=-.266096234984017e+03_dp
  Mat(1,11)=-.196278098216953e+04_dp
  Mat(1,12)=-.505266963449819e+03_dp
  Mat(1,13)=-.110483041928547e+03_dp
  Mat(1,14)=0.430488193758786e+04_dp
  Mat(1,15)=-.139197970442470e+02_dp
  Mat(1,16)=-.144351781922013e+04_dp
  Mat(1,17)=-.265701301773249e+03_dp
  Mat(1,18)=-.428821699139210e+02_dp
  Mat(1,19)=-.443236656693991e+01_dp
  Mat(1,20)=0.309742340203200e+04_dp
  Mat(1,21)=0.386473912295113e+00_dp
  Mat(2,1)=-.514850598717222e+00_dp
  Mat(2,2)=0.213316362570669e+02_dp
  Mat(2,3)=-.302865564916568e+03_dp
  Mat(2,4)=-.198569416607029e+02_dp
  Mat(2,5)=-.460306750911640e+02_dp
  Mat(2,6)=0.270825710321281e+01_dp
  Mat(2,7)=0.184510695601404e+03_dp
  Mat(2,8)=0.156537424620061e+03_dp
  Mat(2,9)=0.190613131168980e+04_dp
  Mat(2,10)=0.277006550460850e+03_dp
  Mat(2,11)=-.568117055198608e+02_dp
  Mat(2,12)=0.428921546783467e+03_dp
  Mat(2,13)=0.142494945404341e+03_dp
  Mat(2,14)=-.541945228489881e+04_dp
  Mat(2,15)=0.233351898912768e+02_dp
  Mat(2,16)=0.104183218654671e+04_dp
  Mat(2,17)=0.331489412844667e+03_dp
  Mat(2,18)=0.660002154209991e+02_dp
  Mat(2,19)=0.997500770521877e+01_dp
  Mat(2,20)=0.560508628472486e+04_dp
  Mat(2,21)=0.209909225990756e+01_dp
  Mat(3,1)=0.203814051719994e+02_dp
  Mat(3,2)=-.283958093739548e+03_dp
  Mat(3,3)=0.173908241235198e+04_dp
  Mat(3,4)=-.195566197110461e+03_dp
  Mat(3,5)=-.138012943339611e+03_dp
  Mat(3,6)=0.523629892715050e+03_dp
  Mat(3,7)=0.859266451736379e+03_dp
  Mat(3,8)=-.805606471979730e+02_dp
  Mat(3,9)=-.468711180560599e+04_dp
  Mat(3,10)=0.889580760829066e+01_dp
  Mat(3,11)=-.782994158054881e+02_dp
  Mat(3,12)=-.437214580089117e+02_dp
  Mat(3,13)=0.112996386047623e+01_dp
  Mat(3,14)=0.401746416262936e+04_dp
  Mat(3,15)=0.104927789918320e+01_dp
  Mat(3,16)=-.139340154288711e+03_dp
  Mat(3,17)=-.170995948015951e+02_dp
  Mat(3,18)=0.545784716783902e+00_dp
  Mat(3,19)=0.971126767581517e+00_dp
  Mat(3,20)=0.141909512967882e+04_dp
  Mat(3,21)=0.994142892628410e+00_dp


  ! calcul des invariants
  Inv2=0.5_dp*(1._dp-(a_11*a_11+a_22*a_22+a_33*a_33+ &
       2._dp*(a_12*a_12+a_13*a_13+a_23*a_23)))

  Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
       a_13*(a_12*a_23-a_22*a_13)

  ! polynome complet de degre 5 des 2 invariants.
  vec(1)=1._dp
  vec(2)=Inv2
  vec(3)=vec(2)*vec(2)
  vec(4)=Inv3
  vec(5)=vec(4)*vec(4)
  vec(6)=vec(2)*vec(4)
  vec(7)=vec(3)*vec(4)
  vec(8)=vec(2)*vec(5)
  vec(9)=vec(2)*vec(3)
  vec(10)=vec(5)*vec(4)
  vec(11)=vec(9)*vec(4)
  vec(12)=vec(3)*vec(5)
  vec(13)=vec(2)*vec(10)
  vec(14)=vec(3)*vec(3)
  vec(15)=vec(5)*vec(5)
  vec(16)=vec(14)*vec(4)
  vec(17)=vec(12)*vec(2)
  vec(18)=vec(12)*vec(4)
  vec(19)=vec(2)*vec(15)
  vec(20)=vec(14)*vec(2)
  vec(21)=vec(15)*vec(4)

  ! calcul des beta_bar (cf annexe C Chung)
  ! attention beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
  !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5

  ! calcul des trois beta en fonction du polynome
  beta(:)=0._dp
  Do i=1,3
     Do j=1,21
        beta(i)=beta(i)+Mat(i,j)*vec(j)
     End do
  End do

  ! calcul des 3 autres pour avoir la normalisation
  beta(4)=3._dp*(-1._dp/7._dp+beta(1)*(1._dp/7._dp+4._dp*Inv2/7._dp+8._dp*Inv3/3._dp)/5._dp- &
       beta(2)*(0.2_dp-8._dp*Inv2/15._dp-14._dp*Inv3/15._dp)- &
       beta(3)*(1._dp/35._dp-24._dp*Inv3/105._dp-4._dp*Inv2/35._dp+ &
       16._dp*Inv2*Inv3/15._dp+8._dp*Inv2*Inv2/35._dp))/5._dp

  beta(5)=6._dp*(1._dp-0.2_dp*beta(1)*(1._dp+4._dp*Inv2)+ &
       7._dp*beta(2)*(1._dp/6._dp-Inv2)/5._dp- &
       beta(3)*(-0.2_dp+2._dp*Inv3/3._dp+4._dp*Inv2/5._dp- &
       8._dp*Inv2*Inv2/5._dp))/7._dp

  beta(6)=-4._dp*beta(1)/5._dp-7._dp*beta(2)/5._dp- &
       6._dp*beta(3)*(1._dp-4._dp*Inv2/3._dp)/5._dp

  ! pour avoir les beta_bar
  Do i=1,6
     beta(i)=beta(i)/3._dp
  End do
  beta(2)=beta(2)/2._dp
  beta(5)=beta(5)/2._dp
  beta(6)=beta(6)/2._dp

  !! calcul des 5 b=a.a
  b_11=a_11*a_11+a_12*a_12+a_13*a_13
  b_22=a_22*a_22+a_12*a_12+a_23*a_23
  b_12=a_11*a_12+a_12*a_22+a_13*a_23
  b_13=a_11*a_13+a_12*a_23+a_13*a_33
  b_23=a_12*a_13+a_22*a_23+a_23*a_33

  !Calcul des 9 termes de a4

  a4(1)=3._dp*beta(4)+6._dp*beta(5)*a_11+3._dp*beta(1)*a_11*a_11+&
       6._dp*beta(2)*b_11+6._dp*beta(6)*a_11*b_11+3._dp*beta(3)*b_11*b_11
  a4(2)=3._dp*beta(4)+6._dp*beta(5)*a_22+3._dp*beta(1)*a_22*a_22+&
       6._dp*beta(2)*b_22+6._dp*beta(6)*a_22*b_22+3._dp*beta(3)*b_22*b_22

  a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2._dp*a_12*a_12)+&
       beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4._dp*a_12*b_12)+&
       beta(3)*(b_11*b_22+2._dp*b_12*b_12)


  a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2._dp*a_12*a_13)+beta(2)*b_23+&
       beta(6)*(a_11*b_23+a_23*b_11+2._dp*(a_12*b_13+a_13*b_12))+beta(3)*&
       (b_11*b_23+2._dp*b_12*b_13)
  a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2._dp*a_12*a_23)+beta(2)*b_13+&
       beta(6)*(a_22*b_13+a_13*b_22+2._dp*(a_12*b_23+a_23*b_12))+beta(3)*&
       (b_22*b_13+2._dp*b_12*b_23)


  a4(6)=3._dp*beta(5)*a_13+3._dp*beta(1)*a_11*a_13+3._dp*beta(2)*b_13+&
       3._dp*beta(6)*(a_11*b_13+a_13*b_11)+3._dp*beta(3)*b_11*b_13
  a4(7)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_11*a_12+3._dp*beta(2)*b_12+&
       3._dp*beta(6)*(a_11*b_12+a_12*b_11)+3._dp*beta(3)*b_11*b_12
  a4(8)=3._dp*beta(5)*a_23+3._dp*beta(1)*a_22*a_23+3._dp*beta(2)*b_23+&
       3._dp*beta(6)*(a_22*b_23+a_23*b_22)+3._dp*beta(3)*b_22*b_23
  a4(9)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_22*a_12+3._dp*beta(2)*b_12+&
       3._dp*beta(6)*(a_22*b_12+a_12*b_22)+3._dp*beta(3)*b_22*b_12
End subroutine IBOF
