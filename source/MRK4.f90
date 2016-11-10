!************* THIS HEADER MUST NOT BE REMOVED *******************!
!** Copyright 2013, Lawrence Campbell and Brian McNeil.         **!
!** This program must not be copied, distributed or altered in  **!
!** any way without the prior permission of the above authors.  **!
!*****************************************************************!

module RK4int

use ParallelInfoType
use TransformInfoType
!use FFTW_Constants

use Globals
use Derivative
use IO
use ParaField

implicit none

  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: dadz_r0, dadz_i0
  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: dadz_r1, dadz_i1
  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: dadz_r2, dadz_i2

  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: A_localtr0, A_localti0
  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: A_localtr1, A_localti1
  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: A_localtr2, A_localti2
  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: A_localtr3, A_localti3

  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: ac_rfield_in, ac_ifield_in

  REAL(KIND=WP), DIMENSION(:,:),ALLOCATABLE :: sEA_re, sEA_im
  REAL(KIND=WP), DIMENSION(:,:),ALLOCATABLE :: sEdA_re, sEdA_im

  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: dxdx, dydx, dz2dx, dpxdx, dpydx, dpz2dx


  REAL(KIND=WP), DIMENSION(:), ALLOCATABLE :: dxm, dxt, xt    ! *t is 'temp', for use in next rhs call...
  REAL(KIND=WP), DIMENSION(:), ALLOCATABLE :: dym, dyt, yt
  REAL(KIND=WP), DIMENSION(:), ALLOCATABLE :: dpxm, dpxt, pxt
  REAL(KIND=WP), DIMENSION(:), ALLOCATABLE :: dpym, dpyt, pyt
  REAL(KIND=WP), DIMENSION(:), ALLOCATABLE :: dz2m, dz2t, z2t
  REAL(KIND=WP), DIMENSION(:), ALLOCATABLE :: dpz2m, dpz2t, pz2t

!  real(kind=wp), allocatable :: dEdx(:,:), dEm(:,:), dEt(:,:), Et(:,:), Enow(:,:)

contains

subroutine rk4par(sZ,h,qD)

  implicit none
!
! Perform 4th order Runge-Kutta integration, tailored
! to Puffin and its method of parallelization:
! This is NOT a general, all-purpose RK4 routine, it
! is specific to Puffin. Includes MPI_gathers and
! scatters etc between calculation of derivatives for
! use with the parallel field derivative.
!
!                ARGUMENTS
!
! y       INPUT/OUTPUT   Electron values
! SA      INPUT/OUTPUT   Field values
! x       INPUT          Propagation distance zbar
! h       INPUT          Step size in zbar

!  REAL(KIND=WP),  DIMENSION(:), INTENT(INOUT) :: sA, A_local
  REAL(KIND=WP),  INTENT(IN)                  :: sZ
  REAL(KIND=WP),                INTENT(IN)  :: h
  LOGICAL, INTENT(INOUT) :: qD

!               LOCAL ARGS
!
! h6         Step size divided by 6
! hh         Half of the step size
! xh         x position incremented by half a step
! dym        Intermediate derivatives
! dyt        Intermediate derivatives
! yt         Incremental solution
! dAdx       Field derivative
! dydx       Electron derivatives

  INTEGER(KIND=IP) :: iy,idydx,iyout,i,p
  REAL(KIND=WP)    :: h6, hh, szh
  !REAL(KIND=WP), DIMENSION(size(y)) :: dym, dyt, yt




  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: dAdx
  REAL(KIND=WP), DIMENSION(:),ALLOCATABLE :: A_localt
  INTEGER(KIND=IP) :: error, trans, tllen43D

!    Transverse nodes

  trans = (NX_G)*(NY_G)

!    Step sizes

  hh = h * 0.5_WP
  h6 = h / 6.0_WP
  szh = sz + hh

  tllen43D = tllen * ntrndsi_G

!$OMP PARALLEL

  dadz_r0 = 0_wp
  dadz_r1 = 0_wp
  dadz_r2 = 0_wp
  dadz_i0 = 0_wp
  dadz_i1 = 0_wp
  dadz_i2 = 0_wp

  A_localtr0 = 0_wp
  A_localtr1 = 0_wp
  A_localtr2 = 0_wp
  A_localtr3 = 0_wp
  A_localti0 = 0_wp
  A_localti1 = 0_wp
  A_localti2 = 0_wp
  A_localti3 = 0_wp



  xt = 0_wp
  yt = 0_wp
  z2t = 0_wp
  pxt = 0_wp
  pyt = 0_wp
  pz2t = 0_wp

  dxdx = 0.0_wp
  dydx = 0.0_wp
  dz2dx = 0.0_wp
  dpxdx = 0.0_wp
  dpydx = 0.0_wp
  dpz2dx = 0.0_wp

  dxm = 0.0_wp
  dxt = 0.0_wp
  dym = 0.0_wp
  dyt = 0.0_wp
  dpxm = 0.0_wp
  dpxt = 0.0_wp
  dpym = 0.0_wp
  dpyt = 0.0_wp
  dz2m = 0.0_wp
  dz2t = 0.0_wp
  dpz2m = 0.0_wp
  dpz2t = 0.0_wp



  A_localtr0 = ac_rfield_in
  A_localti0 = ac_ifield_in


!if (count(abs(ac_rfield) > 1.0E2) > 0) print*, 'HELP IM RUBBUSH AT START I habve ', &
!              (count(abs(ac_rfield) > 1.0E2) > 0), 'bigger than 100...'

!  allocate(DADx(2*local_rows))
!  allocate(A_localt(2*local_rows))

!    A_local from A_big

  if (qD) then

!    if (tTransInfo_G%qOneD) then
!       A_local(1:local_rows)=sA(fst_row:lst_row)
!       A_local(local_rows+1:2*local_rows)=&
!            sA(fst_row+iNumberNodes_G:lst_row+iNumberNodes_G)
!    ELSE
!       CALL getAlocalFS(sA,A_local)
!    END if
!
    qD = .false.
!
  end if

!    First step
!  iy = size(sElX_G)
!  idydx = size(dxdx)

!    Get derivatives

  call derivs(sZ, A_localtr0, A_localti0, sEA_re, sEA_im, &
              sElX_G, sElY_G, sElZ2_G, sElPX_G, sElPY_G, sElGam_G, &
              dxdx, dydx, dz2dx, dpxdx, dpydx, dpz2dx, &
              dadz_r0, dadz_i0, sEdA_re, sEdA_im)

!call mpi_finalize(error)
!stop

!  allocate(dAm(2*local_rows),dAt(2*local_rows))
  !print*, dpydx

!    Increment local electron and field values

  if (qPArrOK_G) then

!$OMP DO
    do i = 1, iNumberElectrons_G
      xt(i) = sElX_G(i)  +  hh*dxdx(i)
      yt(i) = sElY_G(i)  +  hh*dydx(i)
      z2t(i) = sElZ2_G(i)    +  hh*dz2dx(i)
      pxt(i) = sElPX_G(i)    +  hh*dpxdx(i)
      pyt(i) = sElPY_G(i)    +  hh*dpydx(i)
      pz2t(i) = sElGam_G(i)  +  hh*dpz2dx(i)
    end do
!$OMP END DO

!$OMP DO
    do i = 1, tllen43D
      A_localtr1(i) = A_localtr0(i) + hh * dadz_r0(i)
      A_localti1(i) = A_localti0(i) + hh * dadz_i0(i)
    end do
!$OMP END DO

!    Update large field array with new values
!  call local2globalA(A_localt,sA,recvs,displs,tTransInfo_G%qOneD)

    call upd8a(A_localtr1, A_localti1)

  end if



!    Second step
!    Get derivatives

  if (qPArrOK_G) &
    call derivs(szh, A_localtr1, A_localti1, sEA_re, sEA_im, &
       xt, yt, z2t, pxt, pyt, pz2t, &
       dxt, dyt, dz2t, dpxt, dpyt, dpz2t, &
       dadz_r1, dadz_i1, sEdA_re, sEdA_im)





!    Incrementing with newest derivative value...

  if (qPArrOK_G) then

!$OMP DO
    do i = 1, iNumberElectrons_G
      xt(i) = sElX_G(i)      +  hh*dxt(i)
      yt(i) = sElY_G(i)      +  hh*dyt(i)
      z2t(i) = sElZ2_G(i)    +  hh*dz2t(i)
      pxt(i) = sElPX_G(i)    +  hh*dpxt(i)
      pyt(i) = sElPY_G(i)    +  hh*dpyt(i)
      pz2t(i) = sElGam_G(i)  +  hh*dpz2t(i)
    end do
!$OMP END DO


!$OMP DO
    do i = 1, tllen43D
      A_localtr2(i) = A_localtr0(i) + hh * dadz_r1(i)
      A_localti2(i) = A_localti0(i) + hh * dadz_i1(i)
    end do
!$OMP END DO

!    Update full field array

!  call local2globalA(A_localt,sA,recvs,displs,tTransInfo_G%qOneD)

    call upd8a(A_localtr2, A_localti2)

  end if

!    Third step
!    Get derivatives


  if (qPArrOK_G) &
    call derivs(szh, A_localtr2, A_localti2, sEA_re, sEA_im, &
       xt, yt, z2t, pxt, pyt, pz2t, &
       dxm, dym, dz2m, dpxm, dpym, dpz2m, &
       dadz_r2, dadz_i2, sEdA_re, sEdA_im)

!    Incrementing

  if (qPArrOK_G) then

!$OMP DO
    do i = 1, iNumberElectrons_G
      xt(i) = sElX_G(i)  +  h * dxm(i)
      yt(i) = sElY_G(i)  +  h * dym(i)
      z2t(i) = sElZ2_G(i)    +  h * dz2m(i)
      pxt(i) = sElPX_G(i)    +  h * dpxm(i)
      pyt(i) = sElPY_G(i)    +  h * dpym(i)
      pz2t(i) = sElGam_G(i)  +  h * dpz2m(i)
    end do
!$OMP END DO

!$OMP DO
    do i = 1, tllen43D
      A_localtr3(i) = A_localtr0(i) + h * dadz_r2(i)
      A_localti3(i) = A_localti0(i) + h * dadz_i2(i)
    end do
!$OMP END DO

!  call local2globalA(A_localt, sA, recvs, displs, tTransInfo_G%qOneD)

    call upd8a(A_localtr3, A_localti3)

!$OMP DO
    do i = 1, iNumberElectrons_G
      dxm(i) = dxt(i) + dxm(i)
      dym(i) = dyt(i) + dym(i)
      dz2m(i) = dz2t(i) + dz2m(i)
      dpxm(i) = dpxt(i) + dpxm(i)
      dpym(i) = dpyt(i) + dpym(i)
      dpz2m(i) = dpz2t(i) + dpz2m(i)
    end do
!$OMP END DO


!$OMP DO
    do i = 1, tllen43D
      dadz_r2(i) = dadz_r1(i) + dadz_r2(i)
      dadz_i2(i) = dadz_i1(i) + dadz_i2(i)
      dadz_r1(i) = 0_wp
      dadz_i1(i) = 0_wp
    end do
!$OMP END DO

  end if

!    Fourth step

  szh = sz + h


!    Get derivatives

  if (qPArrOK_G) &
      call derivs(szh, A_localtr3, A_localti3, sEA_re, sEA_im, &
       xt, yt, z2t, pxt, pyt, pz2t, &
       dxt, dyt, dz2t, dpxt, dpyt, dpz2t, &
       dadz_r1, dadz_i1, sEdA_re, sEdA_im)


!    Accumulate increments with proper weights

  if (qPArrOK_G) then

!$OMP DO
    do i = 1, iNumberElectrons_G
      sElX_G(i)    = sElX_G(i) + h6 * ( dxdx(i) + dxt(i) + 2.0_WP * dxm(i) )
      sElY_G(i)    = sElY_G(i) + h6 * ( dydx(i) + dyt(i) + 2.0_WP * dym(i) )
      sElZ2_G(i)   = sElZ2_G(i)  + h6 * ( dz2dx(i)  + dz2t(i)  + 2.0_WP * dz2m(i) )
      sElPX_G(i)   = sElPX_G(i)  + h6 * ( dpxdx(i)  + dpxt(i)  + 2.0_WP * dpxm(i) )
      sElPY_G(i)   = sElPY_G(i)  + h6 * ( dpydx(i)  + dpyt(i)  + 2.0_WP * dpym(i) )
      sElGam_G(i)  = sElGam_G(i) + h6 * ( dpz2dx(i) + dpz2t(i) + 2.0_WP * dpz2m(i))
    end do
!$OMP END DO

!$OMP DO
    do i = 1, tllen43D
      ac_rfield_in(i) = ac_rfield_in(i) + h6 * (dadz_r0(i) + dadz_r1(i) &
                            + 2.0_WP * dadz_r2(i))
      ac_ifield_in(i) = ac_ifield_in(i) + h6 * (dadz_i0(i) + dadz_i1(i) &
                            + 2.0_WP * dadz_i2(i))
    end do
!$OMP END DO

!  if (count(abs(dadz_r0) > 0.0_wp) <= 0) print*, 'HELP IM TOO RUBBUSH'

!  if (count(abs(ac_rfield) > 0.0_wp) <= 0) print*, 'HELP IM RUBBUSH'


!if (count(abs(ac_rfield) > 1.0E2) > 0) print*, 'HELP IM RUBBUSH for I habve ', &
!              count(abs(ac_rfield) > 1.0E2) , 'bigger than 100...', &
!              'and I am', tProcInfo_G%rank



    call upd8a(ac_rfield_in, ac_ifield_in)

  end if

!$OMP END PARALLEL

!  call local2globalA(A_local,sA,recvs,displs,tTransInfo_G%qOneD)

!    Deallocating temp arrays

!  deallocate(dAm,dAt,A_localt)

!  deallocate(DADx)

!  deallocate(dadz_r0, dadz_i0)
!  deallocate(dadz_r1, dadz_i1)
!  deallocate(dadz_r2, dadz_i2)
!
!  deallocate(A_localtr0, A_localti0)
!  deallocate(A_localtr1, A_localti1)
!  deallocate(A_localtr2, A_localti2)
!  deallocate(A_localtr3, A_localti3)
!
!  deallocate(DxDx)
!  deallocate(DyDx)
!  deallocate(DpxDx)
!  deallocate(DpyDx)
!  deallocate(Dz2Dx)
!  deallocate(Dpz2Dx)

!   Set error flag and exit

  GOTO 2000

!   Error Handler - Error log Subroutine in CIO.f90 line 709

1000 CALL Error_log('Error in MathLib:rk4',tErrorLog_G)
  PRINT*,'Error in MathLib:rk4'
2000 CONTINUE

end subroutine rk4par






!> @author
!> Lawrence Campbell
!> University of Strathclyde
!> Glasgow, UK

subroutine allact_rk4_arrs()

  integer(kind=ip) :: tllen43D, iNumElms

  tllen43D = tllen * ntrndsi_G

  iNumElms = tllen * (nspindx-1) * (nspindy-1)


!!        Setup electron structured arrays

! E, dEdx, etc are Electron macroparticles
! phase space coords and associated charge weights.
!
! Think carefully about order of these - z2, gamma, px, py, x, y, chi??
! If want this to be efficient, probably want to end up having big loop 
! through macroparticles???

  allocate(DxDx(iNumberElectrons_G))
  allocate(DyDx(iNumberElectrons_G))
  allocate(DpxDx(iNumberElectrons_G))
  allocate(DpyDx(iNumberElectrons_G))
  allocate(Dz2Dx(iNumberElectrons_G))
  allocate(Dpz2Dx(iNumberElectrons_G))

!  allocate(Enow(7_ip, iNumberElectrons_G))
!
!  allocate(dEdx(7_ip, iNumberElectrons_G))
!
!  allocate(dEm(7_ip, iNumberElectrons_G), &
!           dEt(7_ip, iNumberElectrons_G), &
!           Et(7_ip, iNumberElectrons_G))

!  call popStrArrsE(Enow)


!!           Setup field structured arrays

  allocate(dadz_r0(tllen43D), dadz_i0(tllen43D))
  
  allocate(dadz_r1(tllen43D), dadz_i1(tllen43D))
  allocate(dadz_r2(tllen43D), dadz_i2(tllen43D))

  allocate(A_localtr0(tllen43D), A_localti0(tllen43D))
  allocate(A_localtr1(tllen43D), A_localti1(tllen43D))
  allocate(A_localtr2(tllen43D), A_localti2(tllen43D))
  allocate(A_localtr3(tllen43D), A_localti3(tllen43D))

  allocate(ac_rfield_in(tllen43D), ac_ifield_in(tllen43D))


  allocate(dxm(iNumberElectrons_G), &
    dxt(iNumberElectrons_G), xt(iNumberElectrons_G))
  allocate(dym(iNumberElectrons_G), &
    dyt(iNumberElectrons_G), yt(iNumberElectrons_G))
  allocate(dpxm(iNumberElectrons_G), &
    dpxt(iNumberElectrons_G), pxt(iNumberElectrons_G))
  allocate(dpym(iNumberElectrons_G), &
    dpyt(iNumberElectrons_G), pyt(iNumberElectrons_G))
  allocate(dz2m(iNumberElectrons_G), &
    dz2t(iNumberElectrons_G), z2t(iNumberElectrons_G))
  allocate(dpz2m(iNumberElectrons_G), &
    dpz2t(iNumberElectrons_G), pz2t(iNumberElectrons_G))


    allocate(p_nodes(iNumberElectrons_G))
    call alct_e_srtcts(iNumberElectrons_G)
    if (tTransInfo_G%qOneD) then
      allocate(lis_GR(2,iNumberElectrons_G))
    else
      allocate(lis_GR(8,iNumberElectrons_G))
    end if



  allocate(dadz_w(iNumberElectrons_G))

! Element representation of Areal and Aimag

  allocate(sEA_re(8,iNumElms), sEA_im(8,iNumElms))
  allocate(sEdA_re(8,iNumElms), sEdA_im(8,iNumElms))

  call outer2Inner(ac_rfield_in, ac_ifield_in)

  qInnerXYOK_G = .true.

end subroutine allact_rk4_arrs




!> @author
!> Lawrence Campbell
!> University of Strathclyde
!> Glasgow, UK

subroutine deallact_rk4_arrs()

  call inner2Outer(ac_rfield_in, ac_ifield_in)

  deallocate(ac_rfield_in, ac_ifield_in)

  deallocate(dadz_r0, dadz_i0)
  deallocate(dadz_r1, dadz_i1)
  deallocate(dadz_r2, dadz_i2)

  deallocate(A_localtr0, A_localti0)
  deallocate(A_localtr1, A_localti1)
  deallocate(A_localtr2, A_localti2)
  deallocate(A_localtr3, A_localti3)



!  call popNormArrsE(Enow)

!  allocate(Enow(7_ip, iNumberElectrons_G))
!
!  allocate(dEdx(7_ip, iNumberElectrons_G))
!
!  allocate(dEm(7_ip, iNumberElectrons_G), &
!           dEt(7_ip, iNumberElectrons_G), &
!           Et(7_ip, iNumberElectrons_G))

  deallocate(DxDx)
  deallocate(DyDx)
  deallocate(DpxDx)
  deallocate(DpyDx)
  deallocate(Dz2Dx)
  deallocate(Dpz2Dx)

  deallocate(dxm, &
    dxt, xt)
  deallocate(dym, &
    dyt, yt)
  deallocate(dpxm, &
    dpxt, pxt)
  deallocate(dpym, &
    dpyt, pyt)
  deallocate(dz2m, &
    dz2t, z2t)
  deallocate(dpz2m, &
    dpz2t, pz2t)

  deallocate(lis_GR)
  deallocate(p_nodes)
  call dalct_e_srtcts()
    
  deallocate(dadz_w)


  deallocate(sEA_re, sEA_im)
  deallocate(sEdA_re, sEdA_im)


end subroutine deallact_rk4_arrs


!subroutine RK4_inc



!end subroutine RK4_inc

! Note - the dxdz ad intermediates should all be global,
! if allocating outside of RK4 routine.
!
! They should be passed through to rhs / derivs, and
! be local in there, I think....
!
! All those vars being passed into rhs? GLOBAL.
! Only the arrays should be global.
! Same for the equations module...
!
! Scoop out preamble of rhs, defining temp vars.
! These can all be defined outside this routine to make
! it more readable
!
! label consistently - *_g for main global,
!                      *_rg for rhs global,
!                      *_dg for diffraction global
!
! Check 3D und eqns - are they general? i.e. can kx and ky be anything?
!
! Lj, field4elec (real and imag) and dp2f should be global?
!
! Interface for private var access??
! So have dp2f, field4elec and Lj private arrays in the equations module,
! and interface dp2f and field4elec to rhs.f90 to alter them...
! OR alter them through a subroutine....
! SO Lj is common to both field and e eqns...
! whereas field4elec and dp2f are e only.
! So only making Lj and dp2f private to eqns for now
! In fact, they are just globallay defined ATM until
! we get this working....
! Rename eqns to electron eqns or something...
!
! Make rhs / eqns vars global
! fix dp2f interface in rhs.f90 ... DONE
! only allocate / calc dp2f when it will be used!
!

end module rk4int
