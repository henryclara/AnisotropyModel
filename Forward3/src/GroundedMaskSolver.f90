      SUBROUTINE GroundedMaskSolver( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      LOGICAL :: GL, GRD
      Type(Element_t), POINTER :: Element
      Type(Variable_t), POINTER :: Zb,Bed
      REAL(KIND=dp),POINTER :: Mask(:)
      INTEGER, POINTER :: Perm(:)
      INTEGER :: t,n,i
      REAL(KIND=dp), parameter :: eps=0.5_dp ! originally 0.1, try 1.0


      Mask => Solver % Variable % Values
      Perm => Solver % Variable % Perm
      write(*,*) 'Reading Zb..'
      Zb => VariableGet( Solver % Mesh % Variables, 'Zb' )
      IF ( ASSOCIATED( Zb ) ) THEN
            write(*,*) 'Done Zb.'
      ELSE
            write(*,*) 'Zb NOT ASSOCIATED.'
      END IF

      write(*,*) 'Reading bedrock..'
      Bed => VariableGet( Solver % Mesh % Variables, 'bedrock' )
      IF ( ASSOCIATED( Bed ) ) THEN
            write(*,*) 'Done Bedrock.'
      ELSE
            write(*,*) 'BED NOT ASSOCIATED.'
      END IF

! DEFINE INTIAL MASK
! ice sheet ==  1
! ice shelf == -1
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
         n = GetElementNOFNodes()
         Do i=1,n
            write(*,*) 'Writing Element%NodeIndexes(i)', Element%NodeIndexes(i) 
            if (Zb%Values(Zb%Perm(Element%NodeIndexes(i))).GT.(Bed%Values(Bed%Perm(Element%NodeIndexes(i)))+eps)) then
               Mask(Perm(Element%NodeIndexes(i)))=-1.0_dp
            else
               Mask(Perm(Element%NodeIndexes(i)))=+1.0_dp
            endif
         End do
      End do
     write(*,*) 'Done Initiliazing Mask.'

! This code checks whether a node is ungrounded when the other nodes in
! the element are grounded and then gives the ungrounded node a new
! value so that a basal melt is not applied.

!      DO t=1, Solver % NumberOfActiveElements
!          Element => GetActiveElement(t)
!          IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!          n = GetElementNOFNodes()
!          GRD = .False.
!          IF (ANY(Mask(Perm(Element % NodeIndexes(1:n))).NE.-1.0)) &
!              GRD=.True.
!          IF (ANY(Mask(Perm(Element % NodeIndexes(1:n))).NE.-1.0)) &
!              write (*,*) 'Current element has at least one&
!                  & non-ungrounded node'
!          IF (.NOT.GRD) CYCLE
!          DO i=1, n
!              IF (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0) &
!                  write (*,*) 'This node has groundedmask= -1'
!              IF (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0) &
!                  Mask(Perm(Element%NodeIndexes(i))) = 0.0 ! May need to change this?!
!          END DO
!          GRD = .False.
!      END DO



! GROUNDING LINE NEXT NEIGHBOURS
! on ice sheet side  1/2
! on ice shelf side -1/2
!      DO t=1,Solver % NumberOfActiveElements
!        Element => GetActiveElement(t)
!         write(*,*) 'Initial do loop working'
!         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!         n = GetElementNOFNodes()
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.-1.0_dp)) GL=.True.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.-1.0_dp)) &
!                write (*,*) 'At least one element = -1'
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
!                          Mask(Perm(Element%NodeIndexes(i)))=1.0/2.0
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
!               write (*,*) 'This value should be 1/2: ', Mask(Perm(Element%NodeIndexes(i)))
!         End do
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.1.0/2.0)) GL=.True.
!         write (*,*) GL
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=-1.0/2.0
!         End do
!      End do

!! GROUNDING LINE NEXT NEIGHBOURS SECOND ROUND
!! on ice sheet side  1/3
!! on ice shelf side -1/3
!     DO t=1,Solver % NumberOfActiveElements
!         Element => GetActiveElement(t)
!         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!        n = GetElementNOFNodes()
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.1.0/2.0)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=1.0-1.0/3.0
!         End do
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).GT.-0.9999)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=-1.0+1.0/3.0
!         End do
!      End do

!! GROUNDING LINE NEXT NEIGHBOURS SECOND ROUND
!! on ice sheet side  1/4
!! on ice shelf side -1/4
!      DO t=1,Solver % NumberOfActiveElements
!         Element => GetActiveElement(t)
!         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!         n = GetElementNOFNodes()
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.1.0-1.0/3.0)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=1.0-1.0/4.0
!         End do
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).GT.-0.9999)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=-1.0+1.0/4.0
!         End do
!      End do

      CALL InvalidateVariable( Model % Meshes, Solver%Mesh,Solver%Variable%name )
      End


