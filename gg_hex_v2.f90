! ******************************************************************************
!
! Generate a grid in CENTAUR format for generic hexahedral region
!
! Andreas Haselbacher, December 2004
! Modified, December 2004
!
! Manoj Parmar
! Modified, November 2011
!
! ******************************************************************************

  PROGRAM gg_hex
  
  IMPLICIT NONE
  
! ******************************************************************************
! Variables
! ******************************************************************************

  INTEGER, PARAMETER :: SPREAL = SELECTED_REAL_KIND( 6, 37), & 
                        RFREAL = SELECTED_REAL_KIND(15,307), & 
                        CHRLEN = 80
  INTEGER :: bType,choice,errorFlag,fileUnit,fileFormat,i,idum,iFile,il,ih, &
             indx,iq,iRegionGlobal,iw,j,nl,nh,nq,nw,nBounds,nBTris,nCells, &
             nHexs,nPris,nPyrs,nTets,nVars,nVert
  INTEGER :: nBQuads(6)
  INTEGER :: bInfo(3,6)
  INTEGER, ALLOCATABLE :: hex2v(:,:),quad2v(:,:)
  REAL(RFREAL) :: height,length,width,xi,yi,zi,zeta
  REAL(RFREAL), ALLOCATABLE :: xyz(:,:)  
  CHARACTER*(80) :: bName(6),casename,caseTitle,fileName,iFileName
    
! ******************************************************************************
! Start, some data
! ******************************************************************************
    
  casename = 'hex'
    
! ******************************************************************************
! Enter data
! ******************************************************************************
  
!  WRITE(*,*) ' '
!  WRITE(*,*) ' ********************************************** '
!  WRITE(*,*) ' Program to generate grid for hexahedral region ' 
!  WRITE(*,*) '                                                '
!  WRITE(*,*) '               Andreas Haselbacher              '
!  WRITE(*,*) '           Version 2.1, June 29 2005            '
!  WRITE(*,*) '                                                '
!  WRITE(*,*) '                 Manoj Parmar                   '
!  WRITE(*,*) '           Version 3.3, November 16 2011        '
!  WRITE(*,*) ' ********************************************** '
!  WRITE(*,*) ' '
  
  WRITE(*,*) 'Enter length: '
  READ(*,*) length 
 
  WRITE(*,*) 'Enter height: '
  READ(*,*) height
  
  WRITE(*,*) 'Enter width: '
  READ(*,*) width
  
  WRITE(*,*) 'Enter x-initial: '
  READ(*,*) xi 
  
  WRITE(*,*) 'Enter y-initial: '
  READ(*,*) yi 
  
  WRITE(*,*) 'Enter z-initial: '
  READ(*,*) zi
  
  WRITE(*,*) 'Enter number of grid points along length: '
  READ(*,*) nl
 
  WRITE(*,*) 'Enter number of grid points along height: '
  READ(*,*) nh    
  
  WRITE(*,*) 'Enter number of grid points along width: '
  READ(*,*) nw

! TEMPORARY
!  length = Real(0.1)
!  height = Real(5.0)
!  width  = Real(0.1)
  
!  xi = Real(0.0)
!  yi = Real(0.0)
!  zi = Real(0.0)

!  nl = 4
!  nh = 101
!  nw = 2
! END TEMPORARY 
 
! ******************************************************************************
! Grid information
! ******************************************************************************  
  
  WRITE(*,*) ' '
  WRITE(*,*) 'Generating grid...'
  
  nVert   = nl*nw*nh
  nTets   = 0
  nHexs   = (nl-1)*(nw-1)*(nh-1)
  nPris   = 0
  nPyrs   = 0
  nBounds = 6  
  nBTris  = 0
  
  nCells = nHexs  
  
  nBQuads(1) =              (nw-1)*(nh-1)
  nBQuads(2) = nBQuads(1) + (nw-1)*(nh-1) 
  nBQuads(3) = nBQuads(2) + (nl-1)*(nw-1)
  nBQuads(4) = nBQuads(3) + (nl-1)*(nw-1)
  nBQuads(5) = nBQuads(4) + (nl-1)*(nh-1) 
  nBQuads(6) = nBQuads(5) + (nl-1)*(nh-1)          
  
  bType = 300
  bName(1) = 'Solid Wall X'
  bName(2) = 'Solid Wall X'
  bName(3) = 'Solid Wall Y'
  bName(4) = 'Solid Wall Y'
  bName(5) = 'Solid Wall Z'
  bName(6) = 'Solid Wall Z'        
  
  bInfo(1,1:6) = bType
  bInfo(2,1:6) = nbTris
  bInfo(3,1:6) = nBQuads(1:6)
  
! ******************************************************************************
! Generate coordinates
! ******************************************************************************  
  
  ALLOCATE(xyz(3,nVert),STAT=errorFlag)
  IF ( errorFlag /= 0 ) THEN 
    WRITE(*,*) 'ERROR - cannot allocate memory for xyz.'
    STOP
  END IF ! errorFlag
  
  DO ih = 1,nh
    DO iw = 1,nw
      DO il = 1,nl
        indx = il + (iw-1)*nl + (ih-1)*nl*nw
        
! TEMPORARY: Manoj
!        xyz(1,indx) = -0.5D0*length + REAL(il-1)/REAL(nl-1)*length
!        xyz(2,indx) = -0.5D0*height + REAL(ih-1)/REAL(nh-1)*height
!        xyz(3,indx) = -0.5D0*width  + REAL(iw-1)/REAL(nw-1)*width

        xyz(1,indx) = xi + (il-1.0_RFREAL)/REAL(nl-1.0_RFREAL)*length
        xyz(2,indx) = yi + (ih-1.0_RFREAL)/REAL(nh-1.0_RFREAL)*height
        xyz(3,indx) = zi + (iw-1.0_RFREAL)/REAL(nw-1.0_RFREAL)*width

! Stretching in y-direction
!        xyz(1,indx) = xi + REAL(il-1)/REAL(nl-1)*length
!        zeta = REAL(ih-1)/REAL(nh-1)
!        zeta = ERF(2*zeta)*zeta
!        xyz(2,indx) = yi + zeta*height
!        xyz(3,indx) = zi + REAL(iw-1)/REAL(nw-1)*width
! END TEMPORARY        

!        WRITE(*,*) indx,il,iw,ih,xyz(1:3,indx)                 
      END DO ! il
    END DO ! iw
  END DO ! ih
  
! ******************************************************************************
! Connectivity in volume
! ******************************************************************************  
  
  ALLOCATE(hex2v(8,nHexs),STAT=errorFlag)
  IF ( errorFlag /= 0 ) THEN 
    WRITE(*,*) 'ERROR - cannot allocate memory for hex2v.'
    STOP
  END IF ! errorFlag
    
  DO ih = 1,nh-1
    DO iw = 1,nw-1
      DO il = 1,nl-1
        indx = il + (iw-1)*(nl-1) + (ih-1)*(nl-1)*(nw-1)
        
        hex2v(1,indx) = indx + (iw-1) + (ih-1)*(nl*nw - (nl-1)*(nw-1))              
        hex2v(2,indx) = hex2v(1,indx) + nl 
        hex2v(3,indx) = hex2v(1,indx) + nl + 1   
        hex2v(4,indx) = hex2v(1,indx) + 1 
        hex2v(5,indx) = hex2v(1,indx) + nl*nw 
        hex2v(6,indx) = hex2v(1,indx) + nl*nw + nl  
        hex2v(7,indx) = hex2v(1,indx) + nl*nw + nl + 1
        hex2v(8,indx) = hex2v(1,indx) + nl*nw + 1 
        
!        WRITE(*,*) indx,il,iw,ih,hex2v(1:8,indx)                         
      END DO ! il
    END DO ! iw
  END DO ! ih
  
! ******************************************************************************
! Surface connectivity
! ******************************************************************************  
 
  nq = 0
  
  ALLOCATE(quad2v(4,nBQuads(nBounds)),STAT=errorFlag)
  IF ( errorFlag /= 0 ) THEN 
    WRITE(*,*) 'ERROR - cannot allocate memory for quad2v.'
    STOP
  END IF ! errorFlag
  
! x=constant surface at x=xlow  
  
  DO ih = 1,nh-1
    DO iw = 1,nw-1
      nq = nq + 1 
     
      quad2v(1,nq) = (ih-1)*nl*nw + (iw-1)*nl + 1
      quad2v(2,nq) = quad2v(1,nq) + nl
      quad2v(3,nq) = quad2v(2,nq) + nl*nw
      quad2v(4,nq) = quad2v(1,nq) + nl*nw 
    END DO ! iw
  END DO ! ih
  
! x=constant surface at x=xhigh
  
  DO ih = 1,nh-1
    DO iw = 1,nw-1
      nq = nq + 1 
     
      quad2v(1,nq) = ih*nl*nw - (iw-1)*nl
      quad2v(2,nq) = quad2v(1,nq) - nl
      quad2v(3,nq) = quad2v(2,nq) + nl*nw
      quad2v(4,nq) = quad2v(1,nq) + nl*nw 
    END DO ! iw
  END DO ! ih  
    
! y=constant surface at y=ylow
  
  DO iw = 1,nw-1
    DO il = 1,nl-1
      nq = nq + 1 
     
      quad2v(1,nq) = (iw-1)*nl    + il
      quad2v(2,nq) = quad2v(1,nq) + 1
      quad2v(3,nq) = quad2v(2,nq) + nl
      quad2v(4,nq) = quad2v(1,nq) + nl
    END DO ! iw
  END DO ! ih      
  
! y=constant surface at y=yhigh
  
  DO iw = 1,nw-1
    DO il = 1,nl-1
      nq = nq + 1 
     
      quad2v(1,nq) = (nh-1)*nl*nw + (iw-1)*nl + il
      quad2v(2,nq) = quad2v(1,nq) + nl
      quad2v(3,nq) = quad2v(2,nq) + 1
      quad2v(4,nq) = quad2v(1,nq) + 1
    END DO ! iw
  END DO ! ih     
 
! z=constant surface at z=zlow
  
  DO ih = 1,nh-1
    DO il = 1,nl-1
      nq = nq + 1 
     
      quad2v(1,nq) = (ih-1)*nl*nw + nl - (il-1)
      quad2v(2,nq) = quad2v(1,nq) - 1
      quad2v(3,nq) = quad2v(2,nq) + nl*nw
      quad2v(4,nq) = quad2v(1,nq) + nl*nw 
    END DO ! iw
  END DO ! ih    
  
! z=constant surface at z=zhigh
  
  DO ih = 1,nh-1
    DO il = 1,nl-1
      nq = nq + 1 
     
      quad2v(1,nq) = (ih-1)*nl*nw + nl*nw - il
      quad2v(2,nq) = quad2v(1,nq) + 1
      quad2v(3,nq) = quad2v(2,nq) + nl*nw
      quad2v(4,nq) = quad2v(1,nq) + nl*nw 
    END DO ! iw
  END DO ! ih     
  
!  DO iq = 1,nq
!    WRITE(*,*) iq,quad2v(1:4,iq)
!  END DO ! iq
  
  WRITE(*,*) 'Generating grid done.'
 
! ******************************************************************************
! Write grid file for TECPLOT
! ******************************************************************************

  fileName = 'hex.plt'
  fileUnit = 31
  
  WRITE(*,*) ' '
  WRITE(*,*) 'Writing Tecplot grid file...'

  OPEN(fileUnit,FILE=TRIM(fileName),STATUS='UNKNOWN', & 
       FORM='FORMATTED',IOSTAT=errorFlag)
  IF ( errorFlag /= 0 ) THEN 
    WRITE(*,*) 'ERROR - Cannot open file.'
    STOP
  END IF ! errorFlag

  WRITE(fileUnit,*) 'TITLE = "GRID"'
  WRITE(fileUnit,*) 'VARIABLES = "X","Y","Z"'
  WRITE(fileUnit,*) 'ZONE F = FEPOINT, N=',nVert,', E=',nHexs,', ET=BRICK'

  DO j = 1,nVert
    WRITE(fileUnit,'(3(E16.9,2X))') (xyz(i,j),i=1,3)
  END DO ! iVert

  DO j = 1,nHexs
    WRITE(fileUnit,'(8(I10,2X))') (hex2v(i,j),i=1,8)
  END DO ! iCell

  CLOSE(fileUnit,IOSTAT=errorFlag)
  IF ( errorFlag /= 0 ) THEN 
    WRITE(*,*) 'ERROR - Cannot close file.'
    STOP
  END IF ! errorFlag

  WRITE(*,*) 'Writing Tecplot grid file done.'

! ******************************************************************************
! Write file
! ******************************************************************************  
  
  WRITE(*,*) ' '
  WRITE(*,*) 'Enter file format for CENTAUR grid file (Binary - 0, ASCII - 1): '
  READ(*,*) fileFormat
 
  IF (fileFormat == 0) THEN
    WRITE(*,*) ' '
    WRITE(*,*) 'Writing Binary CENTAUR grid file...'
    
    fileUnit = 11
    fileName = 'hex.hyb.bin'
 
    caseTitle = 'hex'
    idum      = 1
  
    OPEN(fileUnit,FILE=fileName,FORM="UNFORMATTED",STATUS="UNKNOWN", & 
      IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN 
      WRITE(*,*) 'ERROR - could not open file.'
      STOP
    END IF ! errorFlag  
    
    WRITE(fileUnit) caseTitle
  
    WRITE(fileUnit) nVert
  
    WRITE(fileUnit) ((xyz(i,j),j=1,nVert),i=1,3)
    WRITE(fileUnit) (idum,j=1,nVert)
  
    WRITE(fileUnit) nTets
  
    WRITE(fileUnit) nHexs
    WRITE(fileUnit) ((hex2v(i,j),j=1,nHexs),i=1,8)
    WRITE(fileUnit) (idum,j=1,nHexs)  
  
    WRITE(fileUnit) nPris
    WRITE(fileUnit) nPyrs
  
    WRITE(fileUnit) nBounds
  
    WRITE(fileUnit) ((bInfo(i,j),j=1,nBounds),i=1,3)
    WRITE(fileUnit) (bName(i),i=1,nBounds)
  
    WRITE(fileUnit) nBTris
    WRITE(fileUnit) nBQuads(nBounds)
    WRITE(fileUnit) ((quad2v(i,j),j=1,nBQuads(nBounds)),i=1,4)
  
  ELSE
    WRITE(*,*) ' '
    WRITE(*,*) 'Writing ASCII CENTAUR grid file...'
    
    fileUnit = 11
    fileName = 'hass.hyb.asc'
 
    caseTitle = 'hex'
    idum      = 1
  
    OPEN(fileUnit,FILE=fileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
      IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN 
      WRITE(*,*) 'ERROR - could not open file.'
      STOP
    END IF ! errorFlag  
    
    WRITE(fileUnit,'(A80)') caseTitle
  
    WRITE(fileUnit,'(I16)') nVert
  
    DO i = 1,3
      WRITE(fileUnit,'(5(E16.9))') (xyz(i,j),j=1,nVert)
    END DO ! i
    WRITE(fileUnit,'(10I16)') (idum,j=1,nVert)
  
    WRITE(fileUnit,'(I16)') nTets
  
    WRITE(fileUnit,'(I16)') nHexs
    DO i=1,8
    WRITE(fileUnit,'(10I16)') (hex2v(i,j),j=1,nHexs)
    END DO
    WRITE(fileUnit,'(10I16)') (idum,j=1,nHexs)  
  
    WRITE(fileUnit,'(I16)') nPris
    WRITE(fileUnit,'(I16)') nPyrs
  
    WRITE(fileUnit,'(I16)') nBounds
  
    DO i = 1,3
      WRITE(fileUnit,'(10I16)') (bInfo(i,j),j=1,nBounds)
    END DO ! i
    WRITE(fileUnit,'(A80)') (bName(i),i=1,nBounds)
  
    WRITE(fileUnit,'(I16)') nBTris
    WRITE(fileUnit,'(I16)') nBQuads(nBounds)
    DO i = 1,4
      WRITE(fileUnit,'(10I16)') (quad2v(i,j),j=1,nBQuads(nBounds))
    END DO ! i
  END IF ! fileFormat

  CLOSE(fileUnit)
  
  WRITE(*,*) 'Writing CENTAUR grid file done.'
  
! ******************************************************************************
! End
! ******************************************************************************
  
  END PROGRAM gg_hex
