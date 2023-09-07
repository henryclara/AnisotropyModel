SUBROUTINE DSolver( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  !******************************************************************************
  !  Calculate D
  !
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

  Integer :: DIM,NMAX,NM,e,n,i,j,k,stat
  Logical :: GotIt

  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(Element_t),POINTER :: Element
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Matrix_t), POINTER :: Systemmatrix

  INTEGER :: DOFs
  INTEGER, POINTER :: DPerm(:)
  REAL(KIND=dp), POINTER :: D(:)

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),FORCE(:)
  Integer(KIND=dp), ALLOCATABLE :: NElements(:)

  REAL(KIND=dp), ALLOCATABLE :: Velocity(:,:)

  Real(KIND=dp) :: Norm,PrevNorm
  !  Real(KIND=dp) :: at,st,totat,totst,CPUTime,RelativeChange
    Real(KIND=dp) :: at,st,totat,totst,RelativeChange

!------------------------------------------------------------------------------
! Say Hello!
!------------------------------------------------------------------------------ 
  WRITE(Message,'(a)') 'Start Assembly'
  CALL Info('DSolver', Message, Level=4)       
  totat = 0.0d0
  totst = 0.0d0
  at = CPUTime()

!------------------------------------------------------------------------------
!    Get constants
!------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  SystemMatrix => Solver % Matrix
  NM = Solver % Mesh % NumberOfNodes

!------------------------------------------------------------------------------
!    Get variables for the solution
!------------------------------------------------------------------------------
  D     => Solver % Variable % Values     ! Nodal values for 
  DPerm => Solver % Variable % Perm       ! Permutations for 
  DOFs = Solver % Variable % DOFs

!------------------------------------------------------------------------------
!    Inicialization
!-----------------------------------------------------------------------------

  NMAX = Model % MaxElementNodes

 ALLOCATE( ElementNodes % x( NMAX ),      &
             ElementNodes % y( NMAX ),    &
             ElementNodes % z( NMAX ),    &
             FORCE(DOFs*NMAX),            &
             STIFF(DOFs*NMAX,DOFs*NMAX ), &
             Velocity(3,NMAX),          &
                                STAT=stat  )
 IF ( stat /= 0 ) THEN
    CALL Fatal('DSolver','Memory allocation error, Aborting.')
 END IF

 CALL DefaultInitialize()

!------------------------------------------------------------------------------
!       Do the assembly
!------------------------------------------------------------------------------

  DO e=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(e)
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes

     CALL GetElementNodes( ElementNodes )
     CALL GetScalarLocalSolution(Velocity(1,1:n), 'Velocity 1')
     CALL GetScalarLocalSolution(Velocity(2,1:n), 'Velocity 2')
     IF(DIM==3) THEN
        CALL GetScalarLocalSolution(Velocity(3,1:n), 'Velocity 3')
     END IF
!------------------------------------------------------------------------------
!       Get element local matrix and vector
!------------------------------------------------------------------------------
     CALL LocalMatrix( STIFF, FORCE,                   &
                       Element,n, ElementNodes, NodeIndexes, &
                       Velocity)
!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------

     CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO

  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
!    Solve System
!------------------------------------------------------------------------------
  at = CPUTime() - at
  st = CPUTime()


  PrevNorm = Solver % Variable % Norm

  Norm = DefaultSolve()

  Solver % Variable % Norm = Norm

  st = CPUTIme()-st
  totat = totat + at
  totst = totst + st

!------------------------------------------------------------------------------
!   Say Good-Bye
!------------------------------------------------------------------------------

  WRITE(Message,'(a,F8.2,F8.2)') 'Assembly: (s)', at, totat
  CALL Info( 'DSolver', Message, Level=4 )
  WRITE(Message,'(a,F8.2,F8.2)') ' Solve:    (s)', st, totst
  CALL Info( 'DSolver', Message, Level=4 )

  ! is the solution converged?
  !---------------------------
  IF ( PrevNorm + Norm /= 0.0d0 ) THEN
     RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
  ELSE
     RelativeChange = 0.0d0
  END IF
  WRITE( Message, * ) 'Result Norm   : ',Norm
  CALL Info( 'D', Message, Level=4 )
  WRITE( Message, * ) 'Relative Change : ',RelativeChange
  CALL Info( 'D', Message, Level=4 )


CONTAINS

  SUBROUTINE LocalMatrix(  STIFF, FORCE,                         &
       Element,n, ElementNodes, NodeIndexes, &
       NVelocity)

    !-----------------------------------------------------------------------
    ! External variables
    !------------------------------------------------------------------------

    REAL(KIND=dp),TARGET :: STIFF(:,:), FORCE(:)
    TYPE(Element_t),POINTER :: Element
    Integer :: n
    TYPE(Nodes_t)   :: ElementNodes
    INTEGER :: NodeIndexes(:)
    Real(KIND=dp) :: NVelocity(:,:)

    !------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),U,V,W,S,SqrtElementMetric
    LOGICAL :: Stat
    INTEGER :: i,j,t,p,q,c
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    REAL(KIND=dp), POINTER :: Sm(:,:),Fm(:)

    REAL(KIND=dp) :: D(DOFs)

    REAL(KIND=dp) :: ss,tt

    INTEGER :: IndexD(3,3),IndexW(3,3)   !Deformations in 3D always as we need fabric in 3D

    !------------------------------------------------------------------------
    ! Inicialization
    !------------------------------------------------------------------------
    FORCE = 0.0d0
    STIFF = 0.0d0

    IntegStuff = GaussPoints( Element )

    c = DOFs
    
    IndexD(1,:)=(/1,4,5/)
    IndexD(2,:)=(/0,2,6/)
    IndexD(3,:)=(/0,0,3/)
    IndexW(1,:)=(/0,7,8/)
    IndexW(2,:)=(/0,0,9/)
    IndexW(3,:)=(/0,0,0/)

    !-------------------------------------------------------------------------
    ! Assembly
    !-------------------------------------------------------------------------

    DO t = 1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       !        Basis function values & derivatives at the integration point:
       !----------------------------------------------------------------------
       stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
            Basis,dBasisdx,Bubbles=.FALSE.)

       !        Correction from metric
       !---------------------------------------------------------------------
       S = S * SqrtElementMetric

       !        Local Values at the integration Point
       !-------------------------------------------------------------
       DO i=1,3
          DO j=i,3
             k=IndexD(i,j)
             IF(dim==2.AND.(i==3.OR.j==3)) THEN
                D(k) =0.0d0
             ELSE
                D(k)= 0.5*(                                &
                     SUM(dBasisdx(1:n,i)*NVelocity(j,1:n))+&
                     SUM(dBasisdx(1:n,j)*NVelocity(i,1:n)))
             END IF
          END DO
       END DO
       DO i=1,3
          DO j=i+1,3
             k=IndexW(i,j)
             IF(dim==2.AND.(i==3.OR.j==3)) THEN
                D(k) =0.0d0
             ELSE
                D(k)= 0.5*(                                &
                     -SUM(dBasisdx(1:n,i)*NVelocity(j,1:n))+&
                     SUM(dBasisdx(1:n,j)*NVelocity(i,1:n)))
             END IF
          END DO
       END DO
       !----------------------------------------------------------------------
       !    Loop over basis functions (of both unknowns and weights)
       !--------------------------------------------------------------
       DO p=1,n
          DO q=1,n

             i = c*(p-1)
             j = c*(q-1)
             Sm => STIFF( i+1:i+c,j+1:j+c )

             !    Stiff Matrix
             !-------------------------------------------------------------
             DO i=1,DOFs
                Sm(i,i) = Sm(i,i) +  s * Basis(q) * Basis(p)
             END DO

             !    Diffusion
             !-------------------------------------------------------------
          END DO

          !    Force Vector
          !-----------------------------------------------------------------
          Fm => FORCE( c*(p-1)+1 : c*(p-1)+c )

          DO i=1,DOFs
             Fm(i) = Fm(i) + s * Basis(p) * D(i)
          END DO

       END DO

    END DO

    !Lumping the stiffness matrix
    ss = 0.d0
    tt = 0.d0
    DO i=1,n*c
       DO j=1,n*c
          ss = ss + STIFF(i,j)
          IF (i /= j) THEN
             STIFF(i,j) = 0.d0
          END IF
       END DO
       tt = tt + STIFF(i,i)
    END DO

    DO i=1,n
       DO j=1,c
          k = c * (i-1) + j
          IF ( tt /= 0.d0 ) THEN
             STIFF(k,k) = STIFF(k,k) * ss / tt
          END IF
       END DO
    END DO
  END SUBROUTINE LocalMatrix

END SUBROUTINE DSolver
