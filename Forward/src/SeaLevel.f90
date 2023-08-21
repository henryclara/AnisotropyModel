FUNCTION getSeaLevel( Model, nodenumber, t ) RESULT( SeaLevel )
        USE types
        USE CoordinateSystems
        USE SolverUtils
        USE ElementDescription
        USE DefUtils
        IMPLICIT NONE
        !header variables
        TYPE(Model_t) :: Model
        INTEGER :: nodenumber
        REAL(KIND=dp) :: SeaLevel, t
        
        SeaLevel = 0.0

        !IF (t < 4000.0) THEN
        !        SeaLevel = 0.0
        !ELSE IF (t < 8000.0) THEN
        !        SeaLevel = -80.0 + 0.02*t
        !ELSE IF (t < 12000.0) THEN
        !        SeaLevel = 80.0
        !ELSE IF (t < 16000.0) THEN
        !        SeaLevel = 320.0 - 0.02*t
        !ELSE IF (t < 20000.0) THEN
        !        SeaLevel = 0.0
        !ELSE IF (t < 24000.0) THEN
        !        SeaLevel = -400.0 + 0.02*t
        !ELSE IF (t < 28000.0) THEN
        !        SeaLevel = 80.0
        !ELSE IF (t < 32000.0) THEN
        !        SeaLevel = 640.0 - 0.02*t
        !ELSE
        !        SeaLevel = 0.0
        !END IF
 
        ! write(*,*) 'Sea level: ', SeaLevel
END FUNCTION getSeaLevel
