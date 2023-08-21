       FUNCTION BedrockBump( Model, nodenumber, Bedrock) RESULT(BedBump)
          USE types
          USE CoordinateSystems
          USE SolverUtils
          USE ElementDescription
          USE DefUtils
          IMPLICIT NONE
          TYPE(Variable_t), POINTER :: RiseIntoShelfSol
          INTEGER, POINTER :: RiseIntoShelfPerm(:)
          REAL(kind=dp), POINTER :: RiseIntoShelfVal(:)
          TYPE(Model_t) :: Model
          TYPE(Solver_t), TARGET :: Solver
          INTEGER :: nodenumber,  NMAX, i, dim
          REAL(KIND=dp) :: BedBump, Bedrock, x,y,y0,x0,sigmax,sigmay,Ampl,Bump
          REAL(KIND=dp) :: ZbInit, Rise2Shelf, RiseIntoShelf, Incline
          LOGICAL :: FirstTime=.True., UnFoundFatal, GotIt

          SAVE FirstTime
          Ampl = GetConstReal (Model % Constants, 'Bump Amplitude', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No Bump Amplitude specified in sif file.')
          END IF
          sigmax = GetConstReal (Model % Constants, 'Sigmax', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No sigmax specified in sif file.')
          END IF
          sigmay = GetConstReal (Model % Constants, 'Sigmay', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No sigmay specified in sif file.')
          END IF
          x0 = GetConstReal (Model % Constants, 'x0', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No x-center (x0) specified in sif file.')
          END IF
          y0 = GetConstReal (Model % Constants, 'y0', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No y-center (y0) specified in sif file.')
          END IF
          ZbInit = GetConstReal (Model % Constants, 'ZbInit', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No inital constant ZbInit specified in sif file.')
          END IF

          RiseIntoShelf = GetConstReal (Model % Constants,'RiseIntoShelf', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No RiseIntoShelf specified in sif file.')
          END IF

          Incline = GetConstReal (Model % Constants,'Incline', GotIt)
          IF (.NOT.GotIt) THEN
                  CALL FATAL('USF BedrockBump', 'No incline specified in sif file.')
          END IF
          x = Model % Nodes % x(nodenumber) ! get coordinates
          y = Model % Nodes % y(nodenumber) ! get coordinates
          
          !Get the amount the ice rise should inpinge on the ice shelf without catching any error messages if fields don't exist
          !RiseIntoShelfSol => VariableGet( Model % Variables, 'RiseIntoShelf',UnFoundFatal=UnFoundFatal)
          !RiseIntoShelfPerm => RiseIntoShelfSol % Perm
          !RiseIntoShelfVal => RiseIntoShelfSol % Values

          !InclineSol => VariableGet( Model % Variables, 'Incline',UnFoundFatal=UnFoundFatal)
          !IncelinePerm => InlclineSol % Perm
          !InclineVal => InclineSol % Values

          !NORMAL BUMP

          Bump = Ampl * exp(-(((x - x0)**2 + (y - y0)**2)**2)/(2*sigmax**4))
          BedBump = Bedrock+Bump
          !write(*,*) 'BUMP',BedBump

          !if (BedBump .GE. (ZbInit+RiseIntoShelf-20.0)) then
                 ! Rise2Shelf= RiseIntoShelf + (x-x0)* tan(Incline*pi/180) 
                  !write(*,*) 'BUMP',Rise2Shelf
                  !BedBump = Bedrock - (Bedrock-ZbInit) + Rise2Shelf
                  !BedBump = ZbInit
          !end if
       END FUNCTION BedrockBump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FUNCTION ZbAdj( Model, nodenumber, Bedrock) RESULT(ZbBump)
           USE types
           USE CoordinateSystems
           USE SolverUtils
           USE ElementDescription
           USE DefUtils
           IMPLICIT NONE
           TYPE(Variable_t), POINTER :: BedSol
           INTEGER, POINTER :: BedPerm(:)
           REAL(kind=dp), POINTER :: BedVal(:)
           TYPE(Model_t) :: Model
           TYPE(Solver_t), TARGET :: Solver
           INTEGER :: nodenumber,  NMAX, i, dim
           REAL(KIND=dp) :: ZbBump, Bedrock
           LOGICAL :: FirstTime=.True., UnFoundFatal, GotIt

           SAVE FirstTime
           ZbBump = GetConstReal (Model % Constants, 'ZbInit', GotIt)
           IF (.NOT.GotIt) THEN
                  CALL FATAL('USF ZbAdj in BedrockBump.f90', 'No inital constant ZbInit specified in sif file.')
           END IF
           

           if (Bedrock > ZbBump) then
                   ZbBump = Bedrock
           end if
           !ZbBump = ZbInit
        END FUNCTION ZbAdj
        FUNCTION ZsAdj( Model, nodenumber, Bedrock) RESULT(ZsBump)
           USE types
           USE CoordinateSystems
           USE SolverUtils
           USE ElementDescription
           USE DefUtils
           IMPLICIT NONE
           TYPE(Variable_t), POINTER :: BedSol
           INTEGER, POINTER :: BedPerm(:)
           REAL(kind=dp), POINTER :: BedVal(:)
           TYPE(Model_t) :: Model
           TYPE(Solver_t), TARGET :: Solver
           INTEGER :: nodenumber,  NMAX, i, dim
           REAL(KIND=dp) :: ZbBump,Bedrock,ZsInit,ZsBump
           LOGICAL :: FirstTime=.True., UnFoundFatal, GotIt

           SAVE FirstTime
           ZbBump = GetConstReal (Model % Constants, 'ZbInit', GotIt)
           IF (.NOT.GotIt) THEN
                  CALL FATAL('USF ZbAdj in BedrockBump.f90', 'No inital constant ZbInit specified in sif file.')
           END IF
           ZsInit = GetConstReal (Model % Constants, 'ZsInit', GotIt)
           IF (.NOT.GotIt) THEN
                  CALL FATAL('USF ZsAdj in BedrockBump.f90', 'No inital constant ZsInit specified in sif file.')
           END IF
           

           if (Bedrock > ZbBump) then
                   ZsBump = Bedrock + ZsInit - ZbBump
           else
                   ZsBump = ZsInit
           end if
           !ZbBump = ZbInit
        END FUNCTION ZsAdj

        FUNCTION IceFluxAtBack( Model, nodenumber, FluxInit) RESULT(IceFlux)
           USE types
           USE CoordinateSystems
           USE SolverUtils
           USE ElementDescription
           USE DefUtils
           IMPLICIT NONE
           TYPE(Variable_t), POINTER :: DepthSol, HeightSol
           INTEGER, POINTER :: DepthPerm(:), HeightPerm(:)
           REAL(kind=dp), POINTER :: DepthVal(:), HeightVal(:)
           TYPE(Model_t) :: Model
           TYPE(Solver_t), TARGET :: Solver
           INTEGER :: nodenumber,  NMAX, i, dim
           REAL(KIND=dp) :: FluxInit, IceFlux
           LOGICAL :: FirstTime=.True., UnFoundFatal, GotIt

           SAVE FirstTime
          HeightSol => VariableGet( Model % Variables, 'Height',UnFoundFatal=UnFoundFatal)
          HeightPerm => HeightSol % Perm
          HeightVal => HeightSol % Values

          DepthSol => VariableGet( Model % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
          DepthPerm => DepthSol % Perm
          DepthVal => DepthSol % Values
          !IceFlux = FluxInit/(HeightVal(HeightPerm(nodenumber))+DepthVal(DepthPerm(nodenumber)))
          IceFlux = 25.6 * Tanh(HeightVal(HeightPerm(nodenumber))/568.2) + 136.3
          !write(*,*) 'FLUX',IceFlux

        END FUNCTION IceFluxAtBack

        !FUNCTION ZsAdj( Model, nodenumber, ZbBump) RESULT(ZsBump)
           !USE types
           !USE CoordinateSystems
           !USE SolverUtils
           !USE ElementDescription
           !USE DefUtils
           !IMPLICIT NONE
           !TYPE(Variable_t), POINTER :: BedSol
           !INTEGER, POINTER :: BedPerm(:)
           !REAL(kind=dp), POINTER :: BedVal(:)
           !TYPE(Model_t) :: Model
           !TYPE(Solver_t), TARGET :: Solver
           !INTEGER :: nodenumber,  NMAX, i, dim
           !REAL(KIND=dp) :: ZsBump, ZbBump, ZbInit, ZsInit
           !LOGICAL :: FirstTime=.True., UnFoundFatal, GotIt

           !SAVE FirstTime
           !ZbInit = GetConstReal (Model % Constants, 'ZbInit', GotIt)
           !IF (.NOT.GotIt) THEN
                  !CALL FATAL('USF ZsAdj in BedrockBump.f90', 'No inital constant ZbInit specified in sif file.')
           !END IF
           !ZsInit = GetConstReal (Model % Constants, 'ZsInit', GotIt)
           !IF (.NOT.GotIt) THEN
                  !CALL FATAL('USF ZsAdj in BedrockBump.f90', 'No inital constant ZsInit specified in sif file.')
           !END IF
           !!ZsBump=ZsInit+ZbInit
           !ZsBump=ZbBump
           !!if (Zb>ZbInit) then
                   !!ZsBump = Zb + 300
           !!else
                   !!ZsBump=30
           !!end if
           !!ZsBump = Zb + (ZsInit-ZbInit)
           !!ZsBump = Zb + 300
!END FUNCTION ZsAdj
