!************* THIS HEADER MUST NOT BE REMOVED *******************!
!** Copyright 2013, Lawrence Campbell and Brian McNeil.         **!
!** This program must not be copied, distributed or altered in  **!
!** any way without the prior permission of the above authors.  **!
!*****************************************************************!

MODULE lattice

USE paratype
USE Globals
USE ArrayFunctions
USE ElectronInit
use gtop2
use initConds

implicit none

integer(kind=ip), parameter :: iUnd = 1_ip, &
                               iChic = 2_ip, &
                               iDrift = 3_ip, &
                               iQuad = 4_ip

integer(kind=ip), allocatable :: iElmType(:)
!integer(kind=ip) :: inum_latt_elms

contains

!    ####################################################




  subroutine setupMods(lattFile, taper, sRho)

    implicit none

!     Sets up elements in wiggler lattice. The elements
!     are read in from the file specified.
!
!     Dr Lawrence Campbell
!     University of Strathclyde
!     2015

    character(32_ip), intent(in) :: LattFile 
    real(kind=wp), intent(inout) :: taper
    real(kind=wp), intent(in) :: sRho


    if (lattFile=='') then
      qMod_G = .false.
      if(tProcInfo_G%qRoot) print*, 'There are no dispersive sections'
    else
      qMod_G = .true.
      if(tProcInfo_G%qRoot) print*, 'There are dispersive sections'
    end if


    IF (qMod_G) then

      modNum=numOfMods(lattFile)

      allocate(D(ModNum),zMod(ModNum),delta(modNum))
      allocate(mf(ModNum),delmz(ModNum),tapers(modNum))

!    Latt file name, number of wigg periods converted to z-bar,
!    slippage in chicane in z-bar, 2 dispersive constants, 
!    number of modules

      call readLatt(lattFile,zMod,delta,D,Dfact,ModNum,taper,sRho,sStepSize)
      ModCount = 1

    else 

      modNum = 1
      allocate(iElmType(modNum))
      iElmType(1) = iUnd

    end if

  end subroutine setupMods




!    #####################################################


  SUBROUTINE readLatt(lattFile,zMod,delta,D,Dfact,ModNum,taper,rho,&
                      sStepSize)

  IMPLICIT NONE

! Subroutine to read the information from the lattice file
! and define some of the dispersion parameters.
!
! Modified from Brian McNeil's 1D lattice FEL code 
! to work with this 3D unaveraged model.
!
!                ARGUMENTS
! 
! lattFile      The name of the lattice file (INPUT)
! zMod          Array containing the cumulative 
!               interaction lengths of each module in 
!               z-bar (OUTPUT)
! delta         Array containing lengths of dispersive 
!               sections in slippage lengths (OUTPUT)
! D             Dispersive factor of the chicane, 
!               D=10*Dfact/6*delta (OUTPUT)
! Dfact         Dispersive strength of the chicane (INPUT)
! ModNum        Number of Modules (INPUT)
! rho           FEL parameter

  CHARACTER(32_IP), INTENT(IN) :: lattFile
  REAL(KIND=WP), DIMENSION(:), INTENT(INOUT) :: zMod,delta,D
  REAL(KIND=WP), INTENT(IN) :: Dfact
  REAL(KIND=WP), INTENT(INOUT)  ::  taper
  REAL(KIND=WP), INTENT(IN) :: rho
  INTEGER(KIND=IP), INTENT(IN) :: ModNum
  REAL(KIND=WP), INTENT(OUT)   :: sStepSize

!                LOCAL VARS

  INTEGER(KIND=IP)   :: i,ios,nw,error,ri,NL
  REAL(KIND=WP)      :: c1

  integer(kind=ip) :: nperlam(size(delmz))

  OPEN(1,FILE=lattFile, IOSTAT=ios, ACTION='READ', POSITION ='REWIND')
  IF (ios /= 0_IP) STOP "OPEN(input file) not performed correctly, IOSTAT /= 0"

!   pi = 4.0_WP*ATAN(1.0_WP)
  c1 = 2.0_WP*rho


!     Read whitespace

  NL = 31_IP      !    Number of lines in header
  
  do ri = 1,NL
  
    read (1,*)
  
  end do

!     Read module data from lattice file

  do i=1,ModNum

    read (1,*) nw, delta(i), mf(i), nperlam(i), tapers(i)  !, resFactor(i) ! Wiggler periods, Chicane slippage periods, aw shift, stepsize

    delmz(i) = 4.0_WP * pi * rho / real(nperlam(i),kind=wp)

!     Calculate cumulative interaction length of modules

    if (i==1) then  
      zMod(i) = real(nw,KIND=WP)
      zMod(i) = 2.0_WP*pi*c1*zMod(i)
    else
      zMod(i) = zMod(i-1)+2.0_WP*pi*c1*real(nw,KIND=WP) 
    end if

  end do
	
  close(1, STATUS='KEEP')

!     Convert from # undulator periods to z-bar

  delta = 2.0_WP*pi*c1*delta
  D = Dfact*10.0/6.0*delta ! The Dispersion parameter

  sStepSize =  delmz(1)
  taper = tapers(1)

  END SUBROUTINE readLatt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE disperse(D,delta,i,sStepSize,sZ)

  IMPLICIT NONE

! Subroutine to apply phase changes to electrons beam
! due to passing through a chicane...
! Modified from Brian McNeil's Phase shift routine to 
! work with p2 (electron momentum in z2 frame) as opposed
! to bar{p}_z
!
!              ARGUMENTS
!
! y_e              electron array data
! D                Dispersion parameter for the chicane
!                  (D=10*Dfact/6*delta, Dfact is the
!                  dispersive strength factor of the chicane
! delta            Slippage in resonant wavelengths
	
  REAL(KIND=WP), INTENT(IN) :: D,delta
  INTEGER(KIND=IP), INTENT(IN) :: i
  REAL(KIND=WP), INTENT(OUT) :: sStepSize
  REAL(KIND=WP), INTENT(INOUT) :: sZ

  INTEGER(KIND=IP)  ::  e_tot
 
  REAL(KIND=WP), ALLOCATABLE :: sgamma_j(:),spx0_offset(:),spy0_offset(:), &
                                sx_offset(:),sy_offset(:)

  REAL(KIND=WP)     :: shift_corr, awo, sx_offseto, sx_offsetn, &
                       sy_offseto, sy_offsetn, sZ_new, spx0_offseto, &
                       spy0_offseto, spx0_offsetn, spy0_offsetn, &
                       beta_av, sEta_eff, sl1

  LOGICAL :: qOKL



  e_tot = iGloNumElectrons_G


!     Propagate through chicane

    sElZ2_G = sElZ2_G - 2.0_WP * D *  &
                 (sElGam_G - 1_wp) &
                 + delta



  END SUBROUTINE disperse	



! ##############################################





  subroutine matchOut()



    real(kind=wp), allocatable :: spx0_offset(:),spy0_offset(:), &
                                  sx_offset(:),sy_offset(:)

    real(kind=wp) :: kx, ky


    kx = kx_und_G
    ky = ky_und_G


    allocate(spx0_offset(iNumberElectrons_G), spy0_offset(iNumberElectrons_G))
    allocate(sx_offset(iNumberElectrons_G),sy_offset(iNumberElectrons_G))

!     Get offsets for start of undulator

    sx_offset =    xOffSet(sRho_G, sAw_G, sGammaR_G, sGammaR_G, &
                           sEta_G, sKBeta_G, sFocusfactor_G, spx0_offset, spy0_offset, &
                           fx_G, fy_G, sZ)

    sy_offset =    yOffSet(sRho_G, sAw_G, sGammaR_G, sGammaR_G, &
                           sEta_G, sKBeta_G, sFocusfactor_G, spx0_offset, spy0_offset, &
                           fx_G, fy_G, sZ)


    if (zUndType_G == 'curved') then

! used for curved pole puffin, the 2 order expansion of cosh and sinh
! allows us to simply add a correction term to the intial position
! when calculating initial conditions, this may need change eventually


        spx0_offset = pxOffset(sZ, srho_G, fy_G) & 
            - 0.5_WP * kx**2 * sElX_G**2 &
            -  0.5_WP * kY**2 * sElY_G**2
     
        spy0_offset = -1_wp *  &
                      ( pyOffset(sZ, srho_G, fx_G) &
                      - kx**2 *  sElX_G  * sElY_G)


    else if (zUndType_G == 'planepole') then 

! plane pole initial conditions are calculated as a 2nd order expansion
! and added as a correction term.



        spx0_offset = pxOffset(sZ, srho_G, fy_G) & 
            - 0.5_WP * (sEta_G / (4 * sRho_G**2)) * sElX_G**2 

        spy0_offset = -1_wp * &
                      pyOffset(sZ, srho_G, fx_G) 


    else

! "normal" PUFFIN case with no off-axis undulator
! field variation


        spx0_offset = pxOffset(sZ, srho_G, fy_G) 

        spy0_offset = -1.0_wp * & 
                     pyOffset(sZ, srho_G, fx_G) 


    end if


!     Add on new offset to initialize beam for undulator module


    sElX_G = sElX_G - sx_offset
    sElY_G = sElY_G - sy_offset
    sElPX_G = sElPX_G - spx0_offset
    sElPY_G = sElPY_G - spy0_offset


    deallocate(spx0_offset,spy0_offset,sx_offset,sy_offset)



  end subroutine matchOut




! ###############################################




  subroutine matchIn(sZ)



    real(kind=wp), allocatable :: spx0_offset(:),spy0_offset(:), &
                                  sx_offset(:),sy_offset(:)

    real(kind=wp) :: kx, ky


    kx = kx_und_G
    ky = ky_und_G


    allocate(spx0_offset(iNumberElectrons_G), spy0_offset(iNumberElectrons_G))
    allocate(sx_offset(iNumberElectrons_G),sy_offset(iNumberElectrons_G))

!     Get offsets for start of undulator

    sx_offset =    xOffSet(sRho_G, sAw_G, sGammaR_G, sGammaR_G, &
                           sEta_G, sKBeta_G, sFocusfactor_G, spx0_offset, spy0_offset, &
                           fx_G, fy_G, sZ)

    sy_offset =    yOffSet(sRho_G, sAw_G, sGammaR_G, sGammaR_G, &
                           sEta_G, sKBeta_G, sFocusfactor_G, spx0_offset, spy0_offset, &
                           fx_G, fy_G, sZ)


    if (zUndType_G == 'curved') then

! used for curved pole puffin, the 2 order expansion of cosh and sinh
! allows us to simply add a correction term to the intial position
! when calculating initial conditions, this may need change eventually


        spx0_offset = pxOffset(sZ, srho_G, fy_G) & 
            - 0.5_WP * kx**2 * sElX_G**2 &
            -  0.5_WP * kY**2 * sElY_G**2
     
        spy0_offset = -1_wp *  &
                      ( pyOffset(sZ, srho_G, fx_G) &
                      - kx**2 *  sElX_G  * sElY_G)


    else if (zUndType_G == 'planepole') then 

! plane pole initial conditions are calculated as a 2nd order expansion
! and added as a correction term.



        spx0_offset = pxOffset(sZ, srho_G, fy_G) & 
            - 0.5_WP * (sEta_G / (4 * sRho_G**2)) * sElX_G**2 

        spy0_offset = -1_wp * &
                      pyOffset(sZ, srho_G, fx_G) 


    else

! "normal" PUFFIN case with no off-axis undulator
! field variation


        spx0_offset = pxOffset(sZ, srho_G, fy_G) 

        spy0_offset = -1.0_wp * & 
                     pyOffset(sZ, srho_G, fx_G) 


    end if


!     Add on new offset to initialize beam for undulator module


    sElX_G = sElX_G + sx_offset
    sElY_G = sElY_G + sy_offset
    sElPX_G = sElPX_G + spx0_offset
    sElPY_G = sElPY_G + spy0_offset


    deallocate(spx0_offset,spy0_offset,sx_offset,sy_offset)

  end subroutine matchIn




! #########################################################


  subroutine initUndulator(iM)

    integer(kind=ip), intent(in) :: iM

! Want to update using arrays describing each module...

!     Update undulator parameter:

    n2col0 = mf(i+1)
    n2col = mf(i+1)
    undgrad = tapers(i+1)
    sz0 = sz



!     Update stepsize    

    sStepSize = delmz(i+1) ! Change step size - make sure is still integer
                           ! of 4pirho in input file!!




  end subroutine initUndulator



! #########################################################
	
  FUNCTION lineCount(fname)

!     Function to count the number of lines on a file
!            
!                ARGUMENTS

  INTEGER(KIND=IP)          :: lineCount
  CHARACTER(*), INTENT(IN)  :: fname

!                LOCAL ARGS

  INTEGER :: ios

  OPEN(1,FILE=fname, IOSTAT=ios, ACTION='READ', POSITION ='REWIND')
  IF (ios /= 0) STOP "OPEN(input file) not performed correctly, IOSTAT /= 0"

  lineCount = 0_IP

  DO 
    lineCount=lineCount+1_IP
    READ(1,*,END=10)
  END DO

  10 CLOSE(1, STATUS='KEEP')

  END FUNCTION lineCount

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  FUNCTION numOfMods(fname)

!     Function to count the number of lines on a file
!            
!                ARGUMENTS

  INTEGER(KIND=IP)          :: numOfMods
  CHARACTER(*), INTENT(IN)  :: fname

!                LOCAL ARGS

  INTEGER :: ios

  OPEN(1,FILE=fname, IOSTAT=ios, ACTION='READ', POSITION ='REWIND')
  IF (ios /= 0) STOP "OPEN(input file) not performed correctly, IOSTAT /= 0"

  READ (1,*)
  READ (1,*)
  READ (1,*) numOfMods

  10 CLOSE(1, STATUS='KEEP')

  END FUNCTION numOfMods

END MODULE lattice

