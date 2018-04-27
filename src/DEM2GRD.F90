! COPYRIGHT (C) 2013 MATT BILSKIE / UNIVERSITY OF CENTRAL FLORIDA /
! LOUISIANA STATE UNIVERSITY
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

!	DEVELOPED BY:
!								
!	MATTHEW V. BILSKIE, Ph.D., E.I. (mbilsk3@lsu.edu)
!       LOUISIANA STATE UNIVERSITY
!	    CENTER FOR COASTAL RESILIENCY (lsu.edu/ccr)
!
!   PLEASE EMAIL WITH QUESTIONS/BUGS

!   PORTIONS OF THIS CODE UTILIZE ROUTINES FOUND IN CRYSTAL FULCHERS SURFACE
!   ROUGHNESS CODES - LINES 298 & 405 
!   (MANNINGS_N_FINDER, SURFACE_ROUGHNESS_CALC, AND SURFACE_CANOPY)
!   ALONG WITH USING FORTRAN TYPES THANKS TO ZACH COBELL

!   CITATION (PLEASE CITE THE FOLLOWING PUBLICATION):
!
!   Matthew V. Bilskie, Scott C. Hagen, Topographic accuracy assessment of bare
!   earth lidar-derived unstructured meshes, Advances in Water Resources,
!   Volume 52, February 2013, Pages 165-177

    MODULE MISC
        CONTAINS
        LOGICAL FUNCTION FileExist(fileName)
            CHARACTER(*),INTENT(IN)       :: fileName
            LOGICAL                         :: fExist
            
            INQUIRE(FILE=TRIM(fileName),EXIST=fExist)
            IF(.NOT.fExist)THEN
                WRITE(*,'(A)')'File '//TRIM(fileName)//' not found.'
                FileExist = .FALSE.
                STOP
            ELSE
                FileExist = .TRUE.
            ENDIF
        END FUNCTION

        REAL(8) FUNCTION Distance(lon1,lat1,lon2,lat2,coord)
            IMPLICIT NONE
            REAL(8),INTENT(IN)          :: lon1
            REAL(8),INTENT(IN)          :: lat1
            REAL(8),INTENT(IN)          :: lon2
            REAL(8),INTENT(IN)          :: lat2
            LOGICAL,INTENT(IN)          :: coord

            REAL(8)                     :: a
            REAL(8)                     :: c
            REAL(8),PARAMETER           :: R = 6371.64D0 !km
            REAL(8)                     :: deg2rad
            REAL(8),PARAMETER           :: PI = 3.14159265359D0
            REAL(8)                     :: dlon
            REAL(8)                     :: dlat

            INTRINSIC                   :: SQRT
            INTRINSIC                   :: DSIN
            INTRINSIC                   :: DCOS
            INTRINSIC                   :: DATAN2

            deg2rad = PI / 180.0D0

            IF(coord)THEN !...If geographic coordinates
                ! Haversine Formula (www.movable-type.co.uk_scripts_latlong.html)
                dlon = (lon2 - lon1) * deg2rad
                dlat = (lat2 - lat1) * deg2rad
                a = DSIN(dlat/2.0D0)**2 + DCOS(lat1*deg2rad)*DCOS(lat2*deg2rad) &
                    * DSIN(dlon/2.0D0)**2
                c = 2.D0 * DATAN2(SQRT(a),SQRT(1.D0-a))
                Distance = R * c * 1000 !...Convert km to m
            ELSE !...Cartesian coordinates
                Distance = SQRT((lon2 - lon1)**2 + (lat2 - lat1)**2)
            ENDIF
                
        END FUNCTION
    END MODULE MISC
    
    MODULE ADCGRID
        USE MISC

        INTEGER,ALLOCATABLE             :: Node_InOut(:)
        REAL(8),ALLOCATABLE             :: ZCount(:)
        REAL(8),ALLOCATABLE             :: ZSum(:)

        TYPE BoundaryListing
            INTEGER                     :: NumNodes
            INTEGER                     :: Code
            INTEGER,ALLOCATABLE         :: N1(:)
            INTEGER,ALLOCATABLE         :: N2(:)
            REAL(8),ALLOCATABLE         :: Crest(:)
            REAL(8),ALLOCATABLE         :: Supercritical(:)
            REAL(8),ALLOCATABLE         :: Subcritical(:)
        END TYPE
    
        TYPE MESH
            CHARACTER(200)              :: title
            INTEGER                     :: ne
            INTEGER                     :: np
            REAL(8),ALLOCATABLE         :: nodes(:,:)
            INTEGER,ALLOCATABLE         :: elem(:,:)
            INTEGER                     :: NumOpenBoundaries
            INTEGER                     :: TotNumOpenBoundaryNodes
            INTEGER                     :: NumLandBoundaries
            INTEGER                     :: TotNumLandBoundaryNodes
            TYPE(BoundaryListing),ALLOCATABLE   :: OceanBC(:)
            TYPE(BoundaryListing),ALLOCATABLE   :: LandBC(:)
        END TYPE MESH
        
    CONTAINS
    
    SUBROUTINE ReadMesh(meshName,myMesh)
        IMPLICIT NONE
        CHARACTER(200),INTENT(IN)       :: meshName
        TYPE(MESH),INTENT(OUT)          :: myMesh
        INTEGER                         :: I,J
        INTEGER                         :: tempI
        INTEGER                         :: tempI2
        LOGICAL                         :: exists
        
        exists = FileExist(meshName)
        
        WRITE(*,'(A)')'Reading mesh...'
        OPEN(UNIT=14,FILE=TRIM(meshName),ACTION='READ')
        READ(14,*)myMesh%title
        READ(14,*)myMesh%ne,myMesh%np
        ALLOCATE(myMesh%nodes(1:myMesh%np,1:3))
        ALLOCATE(myMesh%elem(1:myMesh%ne,1:3))
        DO I = 1, myMesh%np
            READ(14,*)tempI,myMesh%nodes(I,1),myMesh%nodes(I,2), & 
                myMesh%nodes(I,3)
            IF(tempI.NE.I)THEN
                WRITE(*,'(A)')'Mesh needs renumbered...'
                WRITE(*,'(A)')''
                STOP
            ENDIF
        ENDDO
        DO I = 1, myMesh%ne
            READ(14,*)tempI,tempI2,myMesh%elem(I,1),myMesh%elem(I,2), &
                myMesh%elem(I,3)
            IF(tempI.NE.I)THEN
                WRITE(*,'(A)')'Mesh needs renumbered...'
                WRITE(*,'(A)')''
                    STOP
            ENDIF
        ENDDO

        !...Read Open BCs
        READ(14,*,END=100) myMesh%NumOpenBoundaries
        READ(14,*,END=100) myMesh%TotNumOpenBoundaryNodes
        ALLOCATE(myMesh%OceanBC(MyMesh%NumOpenBoundaries))
        DO I = 1, myMesh%NumOpenBoundaries
            READ(14,*) myMesh%OceanBC(I)%NumNodes
            ALLOCATE(myMesh%OceanBC(I)%N1(1:myMesh%OceanBC(I)%NumNodes))
            DO J = 1, myMesh%OceanBC(I)%NumNodes
                READ(14,*), myMesh%OceanBC(I)%N1(J)
            ENDDO
        ENDDO

        !...Read Land BCs
        READ(14,*,END=200) myMesh%NumLandBoundaries
        READ(14,*,END=200) myMesh%TotNumLandBoundaryNodes
        ALLOCATE(myMesh%LandBC(myMesh%NumLandBoundaries))
        DO I = 1, myMesh%NumLandBoundaries
            ReAD(14,*) myMesh%LandBC(I)%NumNodes,&
                myMesh%LandBC(I)%Code
            SELECT CASE(myMesh%LandBC(I)%Code)
                CASE(0,1,10,11,12,20,21,22,52)
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*) myMesh%LandBC(I)%N1(J)
                    ENDDO
                CASE(13,23)
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Crest(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Supercritical(1:myMesh%LandBC(I)%NumNodes))
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*) myMesh%LandBC(I)%N1(J), &
                            myMesh%LandBC(I)%Crest(J), &
                            myMesh%LandBC(I)%Supercritical(J)
                    ENDDO
                CASE(24)
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%N2(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Crest(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Supercritical(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Subcritical(1:myMesh%LandBC(I)%NumNodes))
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*)myMesh%LandBC(I)%N1(J),&
                            myMesh%LandBC(I)%N2(J), &
                            myMesh%LandBC(I)%Crest(J), &
                            myMesh%LandBC(I)%Subcritical(J), &
                            myMesh%LandBC(I)%Supercritical(J)
                    ENDDO
                CASE DEFAULT
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    WRITE(*,'(2A,I0)') 'WARNING: Unknown boundary ',&
                        'condition. ADCIRC TYPE = ',&
                        myMesh%LandBC(I)%Code
                    WRITE(*,'(A)') '        READ AS A SINGLE NODE.'
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*) myMesh%LandBC(I)%N1(J)
                    ENDDO
            END SELECT
        ENDDO

        WRITE(*,'(A)')'done!'
        CLOSE(14)

        RETURN
        
100     CONTINUE
        myMesh%NumOpenBoundaries = 0
        myMesh%TotNumOpenBoundaryNodes = 0
200     CONTINUE
        myMesh%NumLandBoundaries = 0
        myMesh%TotNumLandBoundaryNodes = 0
        RETURN

    END SUBROUTINE ReadMesh

    SUBROUTINE WriteMesh(myMesh,fileName,coord)
        IMPLICIT NONE
        TYPE(MESH),INTENT(IN)           :: myMesh
        CHARACTER(200),INTENT(IN)       :: fileName
        LOGICAL,INTENT(IN)              :: coord

        INTEGER                         :: I,J
        
        WRITE(*,'(A)')'Writing new mesh...'
        OPEN(UNIT=14,FILE=TRIM(fileName),ACTION='WRITE')
        WRITE(14,*)'Interpolated from DEM2GRD.F90'
        WRITE(14,'(I9,x,I9)')myMesh%ne,myMesh%np
        DO I = 1, myMesh%np
            IF(.NOT.coord)THEN !...Cartesian
                WRITE(14,'(I9,x,F15.7,5x,F15.7,5x,F15.7)')I,myMesh%nodes(I,1), &
                myMesh%nodes(I,2),myMesh%nodes(I,3)
            ELSE
                WRITE(14,'(I9,x,F14.10,5x,F14.10,5x,F15.7)')I,myMesh%nodes(I,1), &
                myMesh%nodes(I,2),myMesh%nodes(I,3)
            ENDIF
        ENDDO
        DO I = 1, myMesh%ne
            WRITE(14,*)I,'3',myMesh%elem(I,1),myMesh%elem(I,2),myMesh%elem(I,3)  
        ENDDO
        WRITE(14,'(I10)') myMesh%NumOpenBoundaries
        WRITE(14,'(I10)') myMesh%TotNumOpenBoundaryNodes
        DO I = 1, myMesh%NumOpenBoundaries
            WRITE(14,'(I10)') myMesh%OceanBC(I)%NumNodes
            DO J = 1, myMesh%OceanBC(I)%NumNodes
                WRITE(14,*) myMesh%OceanBC(I)%N1(J)
            ENDDO
        ENDDO
        WRITE(14,*) myMesh%NumLandBoundaries
        WRITE(14,*) myMesh%TotNumLandBoundaryNodes
        DO I = 1, myMesh%NumLandBoundaries
            WRITE(14,'(I6,4x,I6,4x,A,I6)') myMesh%LandBC(I)%NumNodes, &
                myMesh%LandBC(I)%Code,"!=seg",I
            SELECT CASE(myMesh%LandBC(I)%Code)
                CASE(0,1,10,11,12,20,21,22,52)
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10)') myMesh%LandBC(I)%N1(J)
                    ENDDO
                CASE(13,23)
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10,2x,F16.3,2x,F16.3)') myMesh%LandBC(I)%N1(J), &
                            myMesh%LandBC(I)%Crest(J),&
                            myMesh%LandBC(I)%Supercritical(J)
                    ENDDO
                CASE(24)
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10,2x,I10,2x,F16.3,2x,F16.3,2x,F16.3,2x,F16.3)') &
                            myMesh%LandBC(I)%N1(J), &
                            myMesh%LandBC(I)%N2(J), &
                            myMesh%LandBC(I)%Crest(J),&
                            myMesh%LandBC(I)%Supercritical(J), &
                            myMesh%LandBC(I)%Subcritical(J)
                    ENDDO
                CASE DEFAULT
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10)') myMesh%LandBC(I)%N1(J)
                    ENDDO
            END SELECT
        ENDDO

        CLOSE(14)
        WRITE(*,'(A)')'done!'

        RETURN

    END SUBROUTINE WriteMesh

    SUBROUTINE ElementSize(myMesh,coord,esize)
        USE MISC
        IMPLICIT NONE
        TYPE(MESH),INTENT(IN)           :: myMesh
        LOGICAL,INTENT(IN)              :: coord
        REAL(8),ALLOCATABLE,INTENT(OUT) :: esize(:)
        
        INTEGER                         :: I
        REAL(8)                         :: DIST(3)

        ALLOCATE(esize(1:myMesh%np))
        WRITE(*,'(A)')'Computing element size...'
        DO I = 1, myMesh%ne
            DIST(1) = Distance(myMesh%nodes(myMesh%elem(I,1),1), &
                myMesh%nodes(myMesh%elem(I,1),2), &
                myMesh%nodes(myMesh%elem(I,2),1), &
                myMesh%nodes(myMesh%elem(I,2),2), coord)
            DIST(2) = Distance(myMesh%nodes(myMesh%elem(I,2),1), &
                myMesh%nodes(myMesh%elem(I,2),2), &
                myMesh%nodes(myMesh%elem(I,3),1), &
                myMesh%nodes(myMesh%elem(I,3),2), coord)
            DIST(3) = Distance(myMesh%nodes(myMesh%elem(I,3),1), &
                myMesh%nodes(myMesh%elem(I,3),2), &
                myMesh%nodes(myMesh%elem(I,1),1), &
                myMesh%nodes(myMesh%elem(I,1),2), coord)
            IF(DIST(1).GT.eSize(MyMesh%elem(I,1)))THEN
                eSize(MyMesh%elem(I,1)) = DIST(1)
            ENDIF
            IF(DIST(2).GT.eSize(MyMesh%elem(I,2)))THEN
                eSize(MyMesh%elem(I,2)) = DIST(2)
            ENDIF
            IF(DIST(3).GT.eSize(MyMesh%elem(I,3)))THEN
                eSize(MyMesh%elem(I,3)) = DIST(3)
            ENDIF

        ENDDO
#ifdef WRITEMESHSIZE
        OPEN(UNIT=14,FILE='ElementSize.grd',ACTION='WRITE')
        WRITE(14,*)myMesh%title
        WRITE(14,'(I9,x,I9)')myMesh%ne,myMesh%np
        DO I = 1, myMesh%np
            WRITE(14,*)I,myMesh%nodes(I,1),myMesh%nodes(I,2),esize(I)
        ENDDO
        DO I = 1, myMesh%ne
            WRITE(14,*)I,'3',myMesh%elem(I,1),myMesh%elem(I,2),myMesh%elem(I,3)  
        ENDDO
        CLOSE(14)
#endif
        WRITE(*,'(A)')'done!'
    END SUBROUTINE ElementSize
    
    END MODULE ADCGrid
    
    MODULE FltRaster
        USE MISC
    
        TYPE RASTER
            INTEGER                     :: NCOLS
            INTEGER                     :: NROWS
            REAL(8)                     :: x_botlt
            REAL(8)                     :: y_botlt
            REAL(8)                     :: x_y_dist
            REAL(8)                     :: nodata_value
            REAL(4),ALLOCATABLE         :: values(:,:)
        END TYPE RASTER
        
    CONTAINS
        
        SUBROUTINE ReadFlt(hdrName,rasterName,myFLT)
            IMPLICIT NONE
            CHARACTER(200),INTENT(IN)   :: hdrName
            CHARACTER(200),INTENT(IN)   :: rasterName
            TYPE(RASTER),INTENT(OUT)    :: myFLT
            LOGICAL                     :: exists
            CHARACTER(10)               :: ndum
            INTEGER                     :: I
            INTEGER                     :: J
            
            exists = FileExist(hdrName)
            exists = FileExist(rasterName)
            
            WRITE(*,'(A)')'Reading raster (*.flt) file...'
            OPEN(UNIT=10,FILE=TRIM(hdrName),ACTION='READ')
            READ(10,*) ndum, myFLT%NCOLS
            READ(10,*) ndum, myFLT%NROWS
            READ(10,*) ndum, myFLT%x_botlt
            READ(10,*) ndum, myFLT%y_botlt
            READ(10,*) ndum, myFLT%x_y_dist
            READ(10,*) ndum, myFLT%nodata_value
            CLOSE(10)

            IF (ALLOCATED(myFLT%values)) DEALLOCATE(myFLT%values)
            ALLOCATE(myFLT%values(1:myFLT%NROWS,1:myFLT%NCOLS))
            
            OPEN(UNIT=10,FILE=TRIM(rasterName),ACTION='READ', &
                FORM='UNFORMATTED',ACCESS='STREAM')
            DO I = 1, myFLT%NROWS
                READ(10) (myFLT%values(I,J),J=1,myFLT%NCOLS)
            ENDDO
            CLOSE(10)
            WRITE(*,'(A)')'done!'
        END SUBROUTINE ReadFLT
    
    END MODULE FltRaster
    
    PROGRAM DEM2GRD
        
        USE ADCGRID
        USE FltRaster

        CHARACTER(200)                  :: CMD
        CHARACTER(100)                  :: inputFile
        CHARACTER(200)                  :: meshFile
        CHARACTER(200),ALLOCATABLE      :: DEMFile(:)
        CHARACTER(200),ALLOCATABLE      :: HDRFile(:)
        CHARACTER(200)                  :: outMeshFile
        CHARACTER(60)                   :: VERSION = 'VERSION 6.0'

        REAL(8)                         :: mult_fac
        REAL(8)                         :: x_inc
        REAL(8)                         :: y_inc
        REAL(8)                         :: x_botrt
        REAL(8)                         :: y_botrt
        REAL(8)                         :: x_toplt
        REAL(8)                         :: y_toplt
        
        INTEGER                         :: I,J
        INTEGER                         :: cs
        INTEGER                         :: NumDEMFiles
        
        LOGICAL                         :: FoundInputFile
        LOGICAL                         :: exists
        LOGICAL                         :: coord

        TYPE(MESH)                      :: myMesh
        TYPE(MESH)                      :: flagMesh
        TYPE(MESH)                      :: zMesh
        TYPE(RASTER)                    :: myFLT

        IF(IARGC().GE.1)THEN
            I = 0
            DO WHILE(I.LT.IARGC())
                I = I + 1
                CALL GETARG(I,CMD)
                SELECT CASE(TRIM(CMD))
                    CASE('-i','-I')
                        I = I + 1
                        CALL GETARG(I,inputFile)
                        FoundInputFile = .TRUE.
                    CASE('-v','-V','-version')
                        WRITE(*,'(A)')VERSION
                        STOP
                    CASE DEFAULT
                        WRITE(*,'(A)')'ERROR: I did not understand argument: '//TRIM(CMD)
                        STOP
                END SELECT
            ENDDO
        ENDIF
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'******************************************************'
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'                 DEM2GRD.F90                          '
        WRITE(*,'(A)')'                 VERSION 6.2                          '
        WRITE(*,'(A)')' PROGRAM TO INTERPOLATE A FLT/GRD RASTER DEM TO       '
        WRITE(*,'(A)')' ADCIRC MESH NODES.                                   '
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'MATTHEW V. BILSKIE, Ph.D.                             '
        WRITE(*,'(A)')'Matt.Bilskie@gmail.com'
        WRITE(*,'(A)')'Copyright M. Bilskie (2013)'
        WRITE(*,'(A)')'Please cite:'
        WRITE(*,'(A)')'Bilsie & Hagen (2013) Adv. Water Res., 52, 165-177,'
        WRITE(*,'(A)')'http://dx.doi.org_10.1016_j.advwatres.2012.09.003'
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'******************************************************'
        WRITE(*,'(A)')''

        IF(.NOT.FoundInputFile)THEN
            WRITE(*,'(A,$)')'Name of input file:'
            READ(*,*)inputFile
            exists = fileExist(TRIM(inputFile))
        ENDIF

        OPEN(UNIT=10,FILE=TRIM(inputFile),ACTION='READ')
        READ(10,*)
        READ(10,*)
        READ(10,*)meshFile
        READ(10,*)cs !... 0 for cartesian, 1 for lat/lon
        READ(10,*)mult_fac
        READ(10,*)outMeshFile
        READ(10,*)NumDEMFiles
        ALLOCATE(DEMFile(1:NumDEMFiles))
        ALLOCATE(HDRFile(1:NumDEMFiles))
        DO I = 1, NumDEMFiles
            READ(10,*)DEMFile(I)
            HDRFile(I) = TRIM(DEMFile(I))//'.hdr'
            DEMFile(I) = TRIM(DEMFile(I))//'.flt'
        ENDDO
        CLOSE(10)

        IF(cs.EQ.0)THEN
            coord = .FALSE. !...cartesian coordinates
        ELSE
            coord = .TRUE.  !...Geographic (lat/lon)
        ENDIF

        CALL ReadMesh(meshFile,myMesh)
        ALLOCATE(Node_InOut(1:myMesh%np))
        ALLOCATE(ZCount(1:myMesh%np))
        ALLOCATE(ZSum(1:myMesh%np))
        ZCount(:) = 0
        ZSum(:) = 0
        DO I = 1, myMesh%np
            IF(myMesh%nodes(I,3).GT.-1000)THEN
                Node_InOut(I) = -1
            ELSE
                Node_InOut(I) = 1
            ENDIF
        ENDDO

        ! Removing the need to have a second mesh of flagged values
        flagMesh = myMesh
        DO I = 1, flagMesh%NP
            IF(flagMesh%nodes(I,3).LT.-1000)THEN
                flagMesh%nodes(I,3) = flagMesh%nodes(I,3) + 1000.0D0
            ENDIF
        ENDDO

        DO I = 1, NumDEMFiles
            CALL ReadFlt(HDRFile(I),DEMFile(I),myFLT)
            !CALL CAA(myMesh,flagMesh,zMesh,myFLT,mult_fac,coord) 
            CALL CAA(myMesh,flagMesh,myFLT,coord) 
            DO J = 1, myMesh%np
                IF(ZCount(J).NE.0) THEN
                    myMesh%nodes(J,3) = mult_fac * (ZSum(J) / ZCount(J))
                ENDIF
            ENDDO
            ZCount(:) = 0
            ZSum(:) = 0
        ENDDO

!       DO I = 1, myMesh%np
!           IF(ZCount(I).NE.0) THEN
!               myMesh%nodes(I,3) = mult_fac * (ZSum(I) / ZCount(I))
!           ENDIF
!       ENDDO
        
        !CALL WriteMesh(zMesh,outMeshFile,coord) 
        CALL WriteMesh(myMesh,outMeshFile,coord) 

        WRITE(*,'(A)')''
        
    END PROGRAM

    !SUBROUTINE CAA(myMesh,flagMesh,newMesh,myRaster,mf,coord)
    SUBROUTINE CAA(myMesh,flagMesh,myRaster,coord)
        USE ADCGRID
        USE FltRaster

        !TYPE(MESH),INTENT(IN)           :: myMesh
        TYPE(MESH)                      :: myMesh
        TYPE(MESH),INTENT(IN)           :: flagMesh
        !TYPE(MESH),INTENT(OUT)          :: newMesh
        TYPE(RASTER),INTENT(IN)         :: myRaster
        !REAL(8),INTENT(IN)              :: mf !...mult. factor
        LOGICAL,INTENT(IN)              :: coord

        INTEGER                         :: I
        INTEGER                         :: J
        INTEGER                         :: K
        INTEGER                         :: L
        INTEGER                         :: std_count
        INTEGER                         :: counter
        INTEGER                         :: no_data_count
        INTEGER                         :: in_pts
        INTEGER                         :: CA

        REAL(8),ALLOCATABLE             :: eSize(:)
        REAL(8)                         :: x
        REAL(8)                         :: y
        REAL(8)                         :: z
        REAL(8)                         :: dx
        REAL(8)                         :: dy
        REAL(8)                         :: col
        REAL(8)                         :: row
        REAL(8)                         :: decimal
        REAL(8)                         :: rastz
        REAL(8)                         :: avgz
        REAL(8),ALLOCATABLE             :: std_vals(:)
        REAL(8)                         :: std_mean
        REAL(8)                         :: std
        REAL(8)                         :: fz !...final z
        REAL(8)                         :: z_check
        REAL(8),PARAMETER               :: PI = 3.14159265359D0
        REAL(8),PARAMETER               :: R = 6371.64D0 !km

        LOGICAL                         :: stdAvg

        INTRINSIC                       :: ABS

        deg2rad = PI / 180.0D0

        CALL ElementSize(myMesh,coord,eSize)

        WRITE(*,'(A)')'Interpolating...'

        !newMesh = myMesh

        x_inc = myRaster%x_y_dist
        y_inc = myRaster%x_y_dist

        x_botrt = myRaster%x_botlt + myRaster%NCOLS * x_inc
        y_botrt = myRaster%y_botlt
        x_toplt = myRaster%x_botlt
        y_toplt = myRaster%y_botlt + myRaster%NROWS * y_inc

        no_data_count = 0
        in_pts = 0
  
        DO I = 1, myMesh%np
            x = myMesh%nodes(I,1)
            y = myMesh%nodes(I,2)
            z = myMesh%nodes(I,3)
            
            IF(z.GT.-1000) GOTO 2000 !...Only interpolate flagged-values
            !IF(Node_InOut(I).EQ.-1) GOTO 2000 !...Only interpolate flagged-values

            IF((x.GT.x_toplt).AND.(x.LT.x_botrt).AND.(y.GT.y_botrt).AND. &
                    (y.LT.y_toplt))THEN

                !Node_InOut(I) = 1

                dx = x_toplt - x
                dy = y_toplt - y

                n_col = ABS(dx / x_inc) + 1
                n_row = ABS(dy / y_inc) + 1

                IF(n_col.EQ.myRaster%NCOLS+1) n_col = myRaster%NCOLS
                IF(n_row.EQ.myRaster%NROWS+1) n_row = myRaster%NROWS

                OPEN(UNIT=17,FILE="debug_output.txt",POSITION="append")
                WRITE(17,*) 'myRaster%NROWS = ', myRaster%NROWS
                WRITE(17,*) 'myRaster%NCOLS = ', myRaster%NCOLS
                WRITE(17,*) 'x_toplt = ', x_toplt
                WRITE(17,*) 'x_botrt = ', x_botrt
                WRITE(17,*) 'y_toplt = ', y_toplt
                WRITE(17,*) 'y_botrt = ', y_botrt
                WRITE(17,*) 'dx = ', dx
                WRITE(17,*) 'dy = ', dy
                WRITE(17,*) 'n_row = ', n_row
                WRITE(17,*) 'n_col = ', n_col
                CLOSE(17)

                rastz = myRaster%values(n_row,n_col)

                counter = 0
                avgz = 0.D0
                IF(rastz.NE.MyRaster%nodata_value)THEN !...Node has value
                    in_pts = in_pts + 1
                    IF(coord)THEN
                        !CA = CEILING((0.25D0 * eSize(I)) / & 
                            !(R * 1000.D0 * PI * myRaster%x_y_dist / 180.D0))
                        CA = INT((0.25D0 * eSize(I)) / & 
                            (R * 1000.D0 * PI * myRaster%x_y_dist / 180.D0))
                    ELSE
                        !CA = CEILING((0.25D0 * eSize(I)) / x_inc)
                        CA = INT((0.25D0 * eSize(I)) / x_inc)
                    ENDIF

                    ! Flagged values are less than 1000. Already added 1000
                    ! after setting flagMesh in main program
                    ! -1001 is not smoothing
                    ! -1005 is a 5x smoothing, for example
                    ! -1105 is a 5x smoothing, but only interpolate "wetted" values
                    stdAvg = .FALSE.
                    IF(flagMesh%nodes(I,3).EQ.-999)THEN
                        !WRITE(*,*)I,CA
                        stdAvg = .TRUE.
                    ELSEIF(flagMesh%nodes(I,3).LT.-100)THEN
                        CA = CA * ABS(flagMesh%nodes(I,3) + 100)
                    ELSEIF(flagMesh%nodes(I,3).LT.0)THEN
                        CA = CA * ABS(flagMesh%nodes(I,3))
                    ENDIF

                    IF(CA.LT.1)THEN
                        !fz = rastz * mf
                        ZCount(I) = ZCount(I) + 1
                        ZSum(I) = ZSum(I) + rastz
                        GOTO 2100
                    ENDIF
                    
                    IF(ALLOCATED(std_vals)) DEALLOCATE(std_vals)
                    ALLOCATE(std_vals(1:(((n_row+CA)-(n_row-CA))+1)*(((n_col+CA)-(ncol-CA))+1)))
                    std_count = 1

                    DO K = n_row - CA, n_row + CA
                        DO L = n_col - CA, n_col + CA
                            !Make sure K is within the DEM rows & L within the columns
                            IF( (K.GE.CA).AND.(K.LE.myRaster%NROWS).AND.(L.GE.CA).AND.(L.LE.myraster%NCOLS) ) THEN
                            !IF((K.GT.CA).AND.(K.LT.myRaster%NROWS))THEN !...If cell window falls outside the DEM
                                z_check = myRaster%values(K,L)
                                IF(z_check.NE.myRaster%nodata_value)THEN
                                    IF(stdAvg)THEN ! Average the two-standard deviation values
                                       ! sum the values and keep the count
                                       !write(*,*)'here',myRaster%values(K,L)
                                       std_vals(std_count) = myRaster%values(K,L)
                                       std_count = std_count + 1
                                    !IF(flagMesh%nodes(I,3).LT.0)THEN
                                    ELSEIF(flagMesh%nodes(I,3).LT.-100)THEN ! only average "wet" values
                                        IF(myRaster%values(K,L).LE.0)THEN
                                            avgz = avgz + myRaster%values(K,L)
                                            counter = counter + 1
                                            ZSum(I) = ZSum(I) + myRaster%values(K,L)
                                            ZCount(I) = ZCount(I) + 1
                                        ENDIF
                                    ELSEIF(flagmesh%nodes(I,3).LT.0)THEN
                                        avgz = avgz + myRaster%values(K,L)
                                        counter = counter + 1
                                        ZSum(I) = ZSum(I) + myRaster%values(K,L)
                                        ZCount(I) = ZCount(I) + 1
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDDO
                    ENDDO

                    IF(stdAvg)THEN
                        avgz = 0
                        counter = 0
                        ! compute the standard deviation
                        std_mean = SUM(std_vals(1:std_count)) / std_count
                        std = SQRT (SUM((std_vals(1:std_count)-std_mean)**2) / std_count)
                        ! Average the values higher than 2x standard deviation
                        DO J = 1, std_count
                            IF(std_vals(J).GE.(2*std))THEN
                                avgz = avgz + std_vals(J)
                                counter = counter + 1
                            ENDIF
                        ENDDO
                        ZSum(I) = avgz
                        ZCount(I) = counter
                    ENDIF

2100                CONTINUE
                    avgz = 0.D0
                    counter = 0
                ENDIF
            ENDIF
2000    CONTINUE
        ENDDO
        WRITE(*,'(A)')'done!'
    END SUBROUTINE CAA
