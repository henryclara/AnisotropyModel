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

  REAL(KIND=dp) :: Rho, Iota

  REAL(KIND=dp) :: fa,fd
  INTEGER :: a4DOFs,etaDOFs,DDOFs,WDOFs,MaxDOFs
  REAL(KIND=dp), ALLOCATABLE :: a2d(:),a2(:), a2a(:),a4(:),a2o(:),Angle(:),FabricGrid(:)&
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


  SAVE FirstTime,NoNeigh,NeighList,Found,isBC,isa2In,FabricGrid

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


  Rho = ListGetConstReal(Solver % Values, 'Alpha',GotIt) !Swapping to Rho as Fabien
  IF(.NOT. GotIt) CALL Fatal( 'FabricEvolutionSolver','Could find value of Alpha, Aborting.' ) 
  Iota = ListGetConstReal(Solver % Values, 'Iota',GotIt)
  IF(.NOT. GotIt) CALL Fatal( 'FabricEvolutionSolver','Could find value of Iota, Aborting.' )
 
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
       a2o(3),                                 &
       Angle(3),                                 &
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

        !    fourth order orientation tensor
           ! a4 at the middle point (IBOF closure)  
!!! (a2 rentre dans l'ordre : 11, 22, 33, 12, 23 ,13)
        Swap(1:6)=(/1,2,3,4,6,5/)
        CALL IBOF(a2(Swap(1:6)),a4)

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

        ! Fxx
        F(1)=2*a4(1)*D(1) - 2*a2(1)*D(1) + 4*a4(7)*D(4) - 2*a2(4)*D(4) &
             +4*a4(6)*D(5) - 2*a2(5)*D(5) + 2*a4(3)*D(2) + 4*a4(4)*D(6) &
             - 2*D(3)*(a4(1) - a2(1) + a4(3))
        F(1)=F(1)*Iota+2*a2(4)*W(1) + 2*a2(5)*W(2)
        
        ! Fyy
        F(2)=2*a4(3)*D(1) - 2*a2(4)*D(4) + 4*a4(9)*D(4) + 4*a4(5)*D(5) &
             -2*a2(2)*D(2) + 2*a4(2)*D(2) + 4*a4(8)*D(6) - 2*a2(6)*D(6) &
             - 2*D(3)*(a4(3) - a2(2) + a4(2))
        F(2)=F(2)*Iota-2*a2(4)*W(1) + 2*a2(6)*W(3)

           ! Fzz
        F(3)=2*a2(5)*D(5) - 4*a4(6)*D(5) - 4*a4(4)*D(6) - 4*a4(5)*D(5) &
             -2*a2(1)*D(3) + 2*a4(1)*D(3) + 4*a4(3)*D(3) - 4*a4(8)*D(6) &
             +2*a2(6)*D(6) - 2*a2(2)*D(3) + 2*a4(2)*D(3)  &
              - 2*D(1)*(a4(1) - a2(1) + a4(3)) - 4*D(4)*(a4(7) &
              -a2(4) + a4(9)) - 2*D(2)*(a4(3) - a2(2) + a4(2)) 
        F(3)=F(3)*Iota- 2*a2(5)*W(2) -2*a2(6)*W(3)
        
           ! Fxy
        F(4)=2*a4(7)*D(1) - a2(1)*D(4) - a2(4)*D(1) + 4*a4(3)*D(4) &
             +4*a4(4)*D(5) - a2(4)*D(2) - a2(2)*D(4) - a2(5)*D(6) - a2(6)*D(5) &
             +2*a4(9)*D(2) + 4*a4(5)*D(6)  &
             - 2*D(3)*(a4(7) - a2(4) + a4(9)) 
        F(4)=F(4)*Iota- a2(1)*W(1) + a2(2)*W(1)+a2(5)*W(3) + a2(6)*W(2)
        
        ! Fxz
        F(5)=2*a4(6)*D(1) + 4*a4(4)*D(4) + 3*a2(1)*D(5) - a2(5)*D(1) &
             -4*a4(1)*D(5) - 4*a4(3)*D(5) - 4*a4(7)*D(6) + 3*a2(4)*D(6) &
             -a2(6)*D(4) - 2*a4(6)*D(3) + 2*a4(5)*D(2) - 4*a4(9)*D(6) &
             +a2(5)*D(3) - a2(3)*D(5) - 2*a4(5)*D(3)
        F(5)=F(5)*Iota- a2(1)*W(2) -a2(4)*W(3) + a2(6)*W(1) + a2(3)*W(2)
        
        ! Fyz
        F(6)=2*a4(4)*D(1) - 4*a4(7)*D(5) + 3*a2(4)*D(5) - a2(5)*D(4) &
             -4*a4(3)*D(6) + 4*a4(5)*D(4) - 4*a4(9)*D(5) - 2*a4(4)*D(3) &
             +2*a4(8)*D(2) + 3*a2(2)*D(6) - a2(6)*D(2) - 4*a4(2)*D(6) &
             -2*a4(8)*D(3) + a2(6)*D(3) - a2(3)*D(6)
        F(6)=F(6)*Iota- a2(4)*W(2)-a2(5)*W(1) - a2(2)*W(3) + a2(3)*W(3)


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


!!$!!! In parallel, look for the lost ones (Ready for Age but not for a2!)
!!$  IF (ParEnv % PES >1) THEN
!!$
!!$     ! Calculate number of nodes we haven't found in this partition and where could they be.
!!$     RequestSend(1:ParEnv % PES)  % n = 0
!!$     DO k=1,NM
!!$        ! WE NEED NODE K if
!!$        IF(Solver % Mesh % ParallelInfo % INTERFACE(k).AND.(.NOT.Found(k)).AND.(.NOT.isBC(k)))  THEN
!!$
!!$           nlist => Solver % Mesh % ParallelInfo % NeighbourList(k) % Neighbours
!!$           nn=SIZE(Solver % Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
!!$           DO i=1,nn
!!$              precv=nlist(i)
!!$              IF(precv==ParEnv % MyPE) CYCLE
!!$              RequestSend(precv+1)  % n = RequestSend(precv+1)  % n + 1
!!$           END DO
!!$        END IF
!!$     END DO
!!$
!!$     !Now get serious and store all the info needed to send a request
!!$     ! Allocate space
!!$     DO i=1,ParEnv % PEs
!!$        ALLOCATE(RequestSend(i) % gbuff( RequestSend(i)  % n))
!!$     END DO
!!$     ! And again but now storing data
!!$     RequestSend(1:ParEnv % PES)  % n = 0
!!$     DO k=1,NM     
!!$        ! WE NEED NODE K if
!!$        IF(Solver % Mesh % ParallelInfo % INTERFACE(k).AND.(.NOT.Found(k)).AND.(.NOT.isBC(k)))  THEN
!!$           nlist => Solver % Mesh % ParallelInfo % NeighbourList(k) % Neighbours
!!$           nn=SIZE(Solver % Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
!!$           DO i=1,nn
!!$              precv=nlist(i)
!!$              IF(precv==ParEnv % MyPE) CYCLE
!!$              RequestSend(precv+1)  % n = RequestSend(precv+1)  % n + 1
!!$              RequestSend(precv+1)  % gbuff(RequestSend(precv+1)  % n)=Solver % Mesh % ParallelInfo % GlobalDOFs(k)
!!$           END DO
!!$        END IF
!!$     END DO

  !Send number of requested nodes to partitions. They are RequestRecv in the other end
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_iRECV(RequestRecv(i) % n, 1, MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,request(i),ierr )       
!!$     END DO
!!$
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_BSEND(RequestSend(i) % n, 1, MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,ierr)
!!$     END DO
!!$
!!$     CALL MPI_WaitAll( ParEnv % PEs,Request, MPI_STATUSES_IGNORE, ierr )

!!$     !Allocate space for requested nodes from partition i-1 to this partition
!!$     DO i=1,ParEnv % PEs
!!$        ALLOCATE(RequestRecv(i) % gbuff( RequestRecv(i)  % n))
!!$     END DO
!!$
!!$     !Send global DOF of the requested nodes across
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_iRECV(RequestRecv(i) % gbuff,RequestRecv(i)  % n, &
!!$             MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,request(i),ierr )       
!!$     END DO
!!$
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_BSEND(RequestSend(i) % gbuff,RequestSend(i)  % n, &
!!$             MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,ierr)
!!$     END DO
!!$
!!$     CALL MPI_WaitAll( ParEnv % PEs,Request, MPI_STATUSES_IGNORE, ierr )
!!$
!!$     ! Now comes the big question. Do we have that info in this partition?     
!!$     DO i=1,ParEnv % PEs
!!$        !(I'm going to be optimistic in the space allocated)
!!$        ALLOCATE(ReplySend(i) % vbuff(RequestRecv(i)  % n),ReplySend(i) % gbuff(RequestRecv(i)  % n))
!!$        ReplySend(i) % n = 0
!!$        DO j=1,RequestRecv(i) % n
!!$           gk = RequestRecv(i) % gbuff(j)
!!$           k=SearchNode( SystemMatrix % ParallelInfo,gk,Order= SystemMatrix % Perm)
!!$           IF(k<0) CYCLE
!!$           IF(Found(k)) THEN
!!$              ReplySend(i) % n = ReplySend(i) % n + 1
!!$              ReplySend(i) % vbuff(ReplySend(i) % n)=SystemMatrix % RHS(a2Perm(k))
!!$              ReplySend(i) % gbuff(ReplySend(i) % n)=gk
!!$           ELSE
!!$              !Warning should be a bit more precisse than this, it could be in a third partition.
!!$!              IF(SIZE(Solver % Mesh % ParallelInfo % NeighbourList(k) % Neighbours)==2) THEN
!!$!                 PRINT *,'Could not find node',gk,k,' For partition' ,i-1
!!$!              END IF
!!$           END IF
!!$        END DO
!!$     END DO
!!$
!!$     !Send number of Replies to partitions. They are ReplyRecv in the other end
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_iRECV(ReplyRecv(i) % n, 1, MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,request(i),ierr )       
!!$     END DO
!!$
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_BSEND(ReplySend(i) % n, 1, MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,ierr)
!!$     END DO
!!$
!!$     CALL MPI_WaitAll( ParEnv % PEs,Request, MPI_STATUSES_IGNORE, ierr )
!!$
!!$     !Send the global DOF of the found nodes
!!$     DO i=1,ParEnv % PEs
!!$        ALLOCATE(ReplyRecv(i) % gbuff(ReplyRecv(i)  % n))
!!$        CALL MPI_iRECV(ReplyRecv(i) % gbuff,ReplyRecv(i)  % n, &
!!$             MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,request(i),ierr )       
!!$     END DO
!!$
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_BSEND(ReplySend(i) % gbuff,ReplySend(i)  % n, &
!!$             MPI_INTEGER,i-1, 910, MPI_COMM_WORLD,ierr)
!!$     END DO
!!$
!!$     CALL MPI_WaitAll( ParEnv % PEs,request, MPI_STATUSES_IGNORE, ierr )
!!$
!!$     !Send the a2 values of the requested nodes
!!$     DO i=1,ParEnv % PEs
!!$        ALLOCATE(ReplyRecv(i) % vbuff(ReplyRecv(i)  % n))
!!$        CALL MPI_iRECV(ReplyRecv(i) % vbuff,ReplyRecv(i)  % n, &
!!$             MPI_DOUBLE_PRECISION,i-1, 910, MPI_COMM_WORLD,request(i),ierr )
!!$     END DO
!!$
!!$     DO i=1,ParEnv % PEs
!!$        CALL MPI_BSEND(ReplySend(i) % vbuff,ReplySend(i)  % n, &
!!$             MPI_DOUBLE_PRECISION,i-1, 910, MPI_COMM_WORLD,ierr)
!!$     END DO
!!$
!!$     CALL MPI_WaitAll( ParEnv % PEs,request, MPI_STATUSES_IGNORE, ierr )
!!$
!!$     !Finally make it happen!
!!$     DO i=1,ParEnv % PEs
!!$        DO j=1,ReplyRecv(i)  % n
!!$           gk=ReplyRecv(i) % gbuff(j)
!!$           k=SearchNode( SystemMatrix % ParallelInfo,gk,Order= SystemMatrix % Perm)
!!$           IF(k<0) CYCLE
!!$           ka=a2Perm(k)
!!$           a2a=ReplyRecv(i) % vbuff(j)
!!$
!!$           CALL SetMatrixElement( SystemMatrix, ka, ka, 1.0d0 ) 
!!$           SystemMatrix % RHS(ka) = a2a
!!$
!!$           !           a2(ka)=a2a
!!$        END DO
!!$     END DO
!!$  END IF



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

