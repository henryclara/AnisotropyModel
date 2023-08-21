SUBROUTINE AgeSolver( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  !******************************************************************************
  !  Calculate Age. SemiLagrangian Method
  !  Notes:
  !  Reference to Martin and Gudmundsson (2012), Effects on nonlinear rheology and anisotropy on the
  !  relationship between age and depth at ice divides, would be appreciated.
  !
  !  Numerical details:
  !  This solver implement a semilagrangian algoritm with a two-time-level scheme and linear interpolation.
  !  It is based in Staniforth, Andrew, Jean Côté, 1991: Semi-Lagrangian Integration Schemes for Atmospheric
  !  Models—A Review. Mon. Wea. Rev., 119, 2206–2223.
  !
  !  Notes for Reinhard:
  !  I have done 3 shortcuts
  !  -1- As usual, only works under CFL conditions. dt<dx/v. No plans to improve that unless we find a real
  !       case for it.
  !  -2- Inflow BC, it is activated with 'AgeIn=Logical True' in any given BC within the SIF file. It works fine in 2D
  !        BUT in 3D assumes horizontal Advection terms are null (udAgedx=vdAgedy=0) instead of
  !        normal terms (un d Age dn =0). Let's write that in the ToDo list.
  !  -3- In the search of the position of the particle during previous time-step, if the search returns a
  !        position outside the domain it does reposition it to the closest node in the domain. It should be
  !        the closest point in the surface of the domain instead. That is also in the ToDo list.
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
  Integer :: DIM,N,e,i,j,k,kf,km,ka,kd,kd2,kk,NMax,NM,stat
  Logical :: GotIt,FirstTime=.TRUE.,IsThere,ALE

  INTEGER, POINTER :: AgePerm(:)
  REAL(KIND=dp), POINTER :: Age(:),Age0(:)

  TYPE(Element_t),POINTER :: Element
  TYPE(Nodes_t)   :: ElementNodes
  INTEGER,  POINTER ::  NodeIndexes(:)
  TYPE(Matrix_t), POINTER :: Systemmatrix

  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolName
  TYPE(Variable_t), POINTER :: FlowSol,MeshVelSol
  INTEGER, POINTER :: FlowPerm(:),MeshVelPerm(:)
  REAL(KIND=dp), POINTER :: Flow(:),MeshVel(:)
  INTEGER :: FlowDOFs,MeshVelDOFs

  INTEGER :: NNMAX
  INTEGER,  ALLOCATABLE :: NoNeigh(:),NeighList(:,:)


  REAL(KIND=dp) :: xa(3),xm(3),xd(3),alpha(3),um(3)
  REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
  REAL(KIND=dp) :: eps=1.0e-6, eps2=1.0e-6
  INTEGER :: MaxDOFs
  REAL(KIND=dp), ALLOCATABLE :: Vector(:,:)
  REAL(KIND=dp) :: Agea,RcvAgea,TmpAge

  REAL(KIND=dp) :: s1(3),s2(3),PoR(3),d2,d2min
  INTEGER :: emin

  !  REAL(KIND=dp) :: at,st,totat,totst,CPUTime,Norm,PrevNorm,RelativeChange
  REAL(KIND=dp) :: at,st,totat,totst,Norm,PrevNorm,RelativeChange

  !Experimental
  LOGICAL, ALLOCATABLE :: Found(:),isBC(:),isAgeIn(:),isAtBoundary(:)
  LOGICAL, ALLOCATABLE :: SkipThisNodeInThisPartition(:)
  LOGICAL :: Owner=.FALSE., AtInterface=.FALSE.
  TYPE(ValueList_t), POINTER :: BC
  REAL(KIND=dp), ALLOCATABLE :: Cond(:)
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status


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



  SAVE FirstTime,NoNeigh,NeighList,Found,isBC,isAgeIn,isAtBoundary
  SAVE SkipThisNodeInThisPartition
  WRITE(Message,'(a)') 'Start Solver'
  CALL Info('AgeSolver', Message, Level=4)

  !------------------------------------------------------------------------------
  ! Get Constants
  !------------------------------------------------------------------------------

  DIM = CoordinateSystemDimension()
  NMAX = Model % MaxElementNodes
  NM = Solver % Mesh % NumberOfNodes

  IF(DIM==2) THEN
     NNMAX = 10!4
  ELSE
     NNMAX=20!8
  END IF

  MaxDOFs=DIM

  !------------------------------------------------------------------------------
  !    Get variables for the solution
  !------------------------------------------------------------------------------

  Age     => Solver % Variable % Values     ! Nodal values for
  AgePerm => Solver % Variable % Perm       ! Permutations for
  Age0    => Solver % Variable % PrevValues(:,1)

  FlowSolName =  GetString( GetEquation(Solver % Mesh % Elements(1)),'Flow Solution Name', GotIt)
  IF(.NOT.GotIt) THEN
     CALL WARN('AgeSolver','Keyword >Flow Solution Name< not found in section >Equation<')
     CALL WARN('AgeSolver','Taking default value >Flow Solution<')
     WRITE(FlowSolName,'(A)') 'Flow Solution'
  END IF
  FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
  IF ( ASSOCIATED( FlowSol ) ) THEN
     FlowPerm     => FlowSol % Perm
     FlowDOFs     =  FlowSol % DOFs
     Flow         => FlowSol % Values
  ELSE
     WRITE(Message,'(A,A,A)') &
          'Convection flag set to >computed<, but no variable >',FlowSolName,'< found'
     CALL FATAL('AgeSolver',Message)
  END IF

  ! Check if the mesh is moving, if so activate the ALE formulation for the velocity
  MeshVelSol => VariableGet( Solver % Mesh % Variables, "mesh velocity")
  IF ( ASSOCIATED( MeshVelSol ) ) THEN
     ALE = .TRUE.
     CALL INFO( 'AgeSolver', 'ALE formualtion activated', Level=4 )
     MeshVelPerm     => MeshVelSol % Perm
     MeshVelDOFs     =   MeshVelSol % DOFs
     MeshVel         => MeshVelSol % Values
  ELSE
     ALE = .FALSE.
     CALL INFO( 'AgeSolver', 'ALE formualtion NOT activated', Level=4 )
  END IF

  !------------------------------------------------------------------------------
  !    Inicialization
  !-----------------------------------------------------------------------------



  ALLOCATE( ElementNodes % x( NMAX ),      &
       ElementNodes % y( NMAX ),      &
       ElementNodes % z( NMAX ),      &
       Vector(MaxDOFS,NMAX),                  &
       STAT=stat  )
  IF ( stat /= 0 ) THEN
     CALL Fatal('AgeSolver','Memory allocation error, Aborting.')
  END IF

  IF(FirstTime) THEN
     ALLOCATE( NoNeigh(NM),&
          Neighlist(NM,NNMAX),&
          Found(NM),                  &
          Cond(NMAX),                  &
          STAT=stat  )
     IF ( stat /= 0 ) THEN
        CALL Fatal('AgeSolver','Memory allocation error, Aborting.')
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

     ! Flag nodes that are in BC with AgeIn.
     !(We will assume that the horizontal gradient of age is null and ignore the horizontal components of velocity.)
     ALLOCATE(isAgeIn(NM),                  &
          STAT=stat  )
     IF ( stat /= 0 ) THEN
        CALL Fatal('AgeSolver','Memory allocation error, Aborting.')
     END IF
     isAgeIn=.FALSE.
     DO i=1,GetNofBoundaryElements()
        Element=>GetBoundaryElement(i)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
        BC=>GetBC()
        IF(.NOT.ASSOCIATED(BC)) CYCLE

        IF(ListCheckPresent(BC,'AgeIn'))  THEN

           DO j=1,Model % NumberOfBCs
              IF ( Element % BoundaryInfo % Constraint == Model % BCs(j) % Tag ) THEN
                 !                print *,ListGetLogical(Model % BCs(j) % Values,'AgeIn', gotIt )
                 isAgeIn(Nodeindexes)= &
                      ListGetLogical(Model % BCs(j) % Values,'AgeIn', gotIt )
              END IF
           END DO

        END IF

     END DO


     ! Flag nodes that are Dirichlet Boundary conditions

     ALLOCATE(isBC(NM),isAtBoundary(NM),                  &
          STAT=stat  )
     IF ( stat /= 0 ) THEN
        CALL Fatal('AgeSolver','Memory allocation error, Aborting.')
     END IF

     ALLOCATE(SkipThisNodeInThisPartition(NM),                  &
          STAT=stat  )
     IF ( stat /= 0 ) THEN
        CALL Fatal('AgeSolver','Memory allocation error, Aborting.')
     END IF


     ! Check if the node has Dirichlet BC, in that case we will ignore
     isBC=.FALSE.
     isAtBoundary=.FALSE.
     DO i=1,GetNofBoundaryElements()
        Element=>GetBoundaryElement(i)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
        BC=>GetBC()
        IF(.NOT.ASSOCIATED(BC)) CYCLE

        isAtBoundary(Nodeindexes)=.TRUE.

        IF(ListCheckPresent(BC,Solver % Variable % Name))  THEN

           ! Check first if we are using Age Condition = -1
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

  SystemMatrix => Solver % Matrix

  totat = 0.0_dp
  totst = 0.0_dp
  at = CPUTime()
  CALL INFO( 'AgeSolver', 'start assembly', Level=4 )
  CALL DefaultInitialize()

  CALL DefaultFinishAssembly()

  Norm=0.0d0
  SkipThisNodeInThisPartition=.FALSE.

  DO k=1,NM
     Found(k)=.FALSE.

     IF(IsBC(k)) THEN
        Found(k)=.TRUE.
        CYCLE
     END IF


     ka=AgePerm(k)
     kf=FlowPerm(k)
     km=MeshVelPerm(k)

     xa(1) = Solver % Mesh % Nodes % x(k)
     xa(2) = Solver % Mesh % Nodes % y(k)
     xa(3) = Solver % Mesh % Nodes % z(k)

     IF(ALE) THEN
        IF(DIM==2.AND.FlowDOFs==3) THEN
           alpha(1)=Dt*(Flow(FlowDOFs*kf-2)-MeshVel(MeshVelDOFs*km-1))
           alpha(2)=Dt*(Flow(FlowDOFs*kf-1)-MeshVel(MeshVelDOFs*km))
           alpha(3)=0._dp
        ELSE IF(DIM==3.AND.FlowDOFs==4) THEN
           alpha(1)=Dt*Flow(FlowDOFs*kf-3)
           alpha(2)=Dt*Flow(FlowDOFs*kf-2)
           alpha(3)=Dt*Flow(FlowDOFs*kf-1)
        ELSE
           CALL Fatal('AgeSolver','DIM AND FlowDOFS do not combine.  Aborting.')
        END IF
     ELSE
        IF(DIM==2.AND.FlowDOFs==3) THEN
           alpha(1)=Dt*Flow(FlowDOFs*kf-2)
           alpha(2)=Dt*Flow(FlowDOFs*kf-1)
           alpha(3)=0._dp
        ELSE IF(DIM==3.AND.FlowDOFs==4) THEN
           alpha(1)=Dt*Flow(FlowDOFs*kf-3)
           alpha(2)=Dt*Flow(FlowDOFs*kf-2)
           alpha(3)=Dt*Flow(FlowDOFs*kf-1)
        ELSE
           CALL Fatal('AgeSolver','DIM AND FlowDOFS do not combine.  Aborting.')
        END IF
     END IF

     IF(isAgeIn(k)) THEN
        DO j=1,DIM-1
           alpha(j)=0.0d0
        END DO
     END IF

     xm(1:3)=xa(1:3)-alpha(1:3)/2._dp

     IsThere=PointInDomain(Element, ElementNodes, NodeIndexes, n, xm, LocalCoordinates,eps, &
          k,Solver,NoNeigh,NeighList)

     !! Note To Carlos: Putting xm to the closest node is done for all partition boundaries,
     !! (i.e. AtInterface == TRUE) which I don't think is what we want.
     Owner = (Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours(1) == ParEnv % MyPE)
     ! AtInterface = Solver % Mesh % ParallelInfo % INTERFACE(k)   !! [VV]

     IF(.NOT.IsThere) THEN
         IsThere=ClosestPointInBoundary(Element, ElementNodes, NodeIndexes, n, xm, LocalCoordinates,eps2, &
             k,Solver,NoNeigh,NeighList)
     END IF

     DO j=1,DIM
        IF(ALE) THEN
           Vector(j,1:n)=Flow(FlowDOFs*(FlowPerm(NodeIndexes)-1)+j)&
                -MeshVel(MeshVelDOFs*(MeshVelPerm(NodeIndexes)-1)+j)
        ELSE
           Vector(j,1:n)=Flow(FlowDOFs*(FlowPerm(NodeIndexes)-1)+j)
        END IF
     END DO

     IF(isAgeIn(k)) THEN
        DO j=1,DIM-1
           Vector(j,1:n)=0.0d0
        END DO
     END IF


     um=0.0d0
     DO j=1,DIM
        um(j)=InterpolateInElement(Element,Vector(j,1:n),&
             LocalCoordinates(1),LocalCoordinates(2), LocalCoordinates(3))
     END DO

     alpha(1)=Dt*um(1)
     alpha(2)=Dt*um(2)
     alpha(3)=Dt*um(3)

     xd(1)=xa(1)-alpha(1)
     xd(2)=xa(2)-alpha(2)
     xd(3)=xa(3)-alpha(3)

     IsThere=PointInDomain(Element, ElementNodes, NodeIndexes, n, xd, LocalCoordinates,eps, &
          k,Solver,NoNeigh,NeighList)
     IF(.NOT.IsThere) THEN
        IsThere=ClosestPointInBoundary(Element, ElementNodes, NodeIndexes, n, xd, LocalCoordinates,eps2, &
             k,Solver,NoNeigh,NeighList)
     END IF




!!!! Should be the solution in the previous time-step Age0 but let's do it this way for now
     Vector(1,1:n)=Age(AgePerm(Element % NodeIndexes))

     Agea=InterpolateInElement(Element,Vector(1,1:n), LocalCoordinates(1), &
          LocalCoordinates(2), LocalCoordinates(3) )+Dt

    !! Note to Carlos: I noticed that in my specific case different ages were calculated
    !! for nodes belonging to more than one partition. This is a quickfix harmonizing the
    !! ages across paritions. I take the minimum Age because in my case the vertical
    !! streaks were alwasy too old. This helps with the symptom, but is not the cause.
    !!___________________________________________________________________________________________________________________
    IF (1 == 1) THEN
    IF (AtInterface .AND. ParEnv % PEs .gt. 1) THEN
       DO i=1,SIZE(Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
         IF (Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours(i) == ParEnv % MyPE) CYCLE
         !Send Age to neighboring partition(s)
         CALL MPI_BSEND(Agea, 1, MPI_DOUBLE, Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours(i), &
         200,MPI_COMM_WORLD, ierr)

       END DO
       TmpAge=Agea
       DO i=1,SIZE(Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
         IF (Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours(i) == ParEnv % MyPE) CYCLE
         !Receive Age from neighboring partition(s)
         CALL MPI_RECV(RcvAgea, 1, MPI_DOUBLE, Solver %  Mesh % ParallelInfo % NeighbourList(k) % Neighbours(i),  &
         200,MPI_COMM_WORLD,status, ierr)
         IF (RcvAgea .lt. TmpAge) THEN
             TmpAge = RcvAgea
         END IF
       END DO
       Agea = TmpAge
    END IF
    END IF
    !!END________________________________________________________________________________________________________________

     Found(k)=.TRUE.
     CALL SetMatrixElement( SystemMatrix, ka, ka, 1.0d0 )
     SystemMatrix % RHS(ka) = Agea
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
  CALL INFO( 'AgeSolver', Message, Level=4 )
  WRITE( Message, * ) 'Relative Change : ',RelativeChange
  CALL INFO( 'AgeSolver', Message, Level=4 )

  FirstTime=.FALSE.

CONTAINS

   FUNCTION  ClosestPointInBoundary(Element,ElementNodes,NodeIndexes,n, xx, LocalCoordinates,eps&
        ,k,Solver,NoNeigh,NeighList) RESULT(IsThere)
     IMPLICIT NONE
     TYPE(Solver_t) :: Solver
     TYPE(Element_t),POINTER :: Element,ElementOK
     INTEGER,  POINTER ::  NodeIndexes(:)
     Real(KIND=dp) :: xx(:),LocalCoordinates(:),eps,yy(3),yymin(3),dist2,dist2min
     Integer :: k,i,j,e,n,DIM,NumberOfFaces
     Logical :: IsThere
     Integer :: NoNeigh(:),NeighList(:,:)
     TYPE(Nodes_t)   :: ElementNodes

     dist2min=1e20
     ! For now sent to the closest node
     DO i=1,NoNeigh(k)
        e=NeighList(k,i)
        Element => Solver % Mesh % Elements(e)

        n = Element % Type % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        DO j=1,n
           yy=(/Solver % Mesh % Nodes % x( NodeIndexes(j)),&
                Solver % Mesh % Nodes % y( NodeIndexes(j)),&
                Solver % Mesh % Nodes % z( NodeIndexes(j))/)

           dist2=DistanceSquared(xx,yy)
           IF(dist2<dist2min) THEN
              dist2min=dist2
              yymin=yy;
           END IF
        END DO
     END DO

     IsThere=PointInDomain(Element, ElementNodes, NodeIndexes, n, yymin, LocalCoordinates,eps, &
          k,Solver,NoNeigh,NeighList)

     IF (.NOT.IsThere) CALL FATAL('AgeSolver','Cannot find closest node to a point. Stop and cry!')

   END FUNCTION ClosestPointInBoundary


   !Find a point with coordinates xx in an element around node k
  FUNCTION  PointInDomain(Element,ElementNodes,NodeIndexes,n, xx, LocalCoordinates,eps&
       ,k,Solver,NoNeigh,NeighList) RESULT(IsThere)
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    TYPE(Element_t),POINTER :: Element
    INTEGER,  POINTER ::  NodeIndexes(:)
    Real(KIND=dp) :: xx(:),LocalCoordinates(:),eps
    Integer :: k,i,e,n
    Logical :: IsThere
    Integer :: NoNeigh(:),NeighList(:,:)
    TYPE(Nodes_t)   :: ElementNodes

     IsThere=.FALSE.
     DO i=1,NoNeigh(k)
        e=NeighList(k,i)
        Element => Solver % Mesh % Elements(e)

        n = Element % Type % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x( NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y( NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z( NodeIndexes)

        IF ( PointInElement( Element, ElementNodes, xx, LocalCoordinates,NumericEps=eps) )  THEN
           IsThere=.TRUE.
           EXIT
        END IF
     END DO

   END FUNCTION PointInDomain

   !Square of the distance between 2 points P1 and P2
   FUNCTION DistanceSquared(p1,p2) RESULT(d2)
     IMPLICIT NONE
     DOUBLE PRECISION :: p1(:),p2(:),d2

     d2=(p2(1)-p1(1))*(p2(1)-p1(1))+(p2(2)-p1(2))*(p2(2)-p1(2))+(p2(3)-p1(3))*(p2(3)-p1(3))

   END FUNCTION DistanceSquared

END SUBROUTINE AgeSolver
