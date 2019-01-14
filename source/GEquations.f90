! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

module Equations


use paratype
use ArrayFunctions
use Globals
use rhs_vars

implicit none



!private real(kind=wp), allocatable :: dp2f(:), sField4ElecReal(:), &
!                                      sField4ElecImag(:) !, &
!                                      !Lj(:)



!real(kind=wp), allocatable :: dp2f(:), sField4ElecReal(:), &
!                              sField4ElecImag(:) !, &
                              !Lj(:)



contains

  subroutine dppdz_r_f(sx, sy, sz2, spr, spi, sgam, &
                       sZ, sdpr, qOK)

  	implicit none


    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)
    real(kind=wp), intent(in) :: sZ
    real(kind=wp), contiguous, intent(out) :: sdpr(:)

    real(kind=wp) :: szt

    logical, intent(inout) :: qOK


    LOGICAL :: qOKL

    qOK = .false.

!$OMP WORKSHARE
    sdpr = sInv2rho * ( n2col * byu  &
                        - sEta_G * sp2 / sKappa_G**2 *    &
                        sField4ElecReal ) &
           + sKappa_G * spi / sgam * (1 + sEta_G * sp2) &
               * n2col * bzu
!$OMP END WORKSHARE

! Set the error flag and exit

    qOK = .true.

    !goto 2000

    !1000 call Error_log('Error in equations:dppdz_r',tErrorLog_G)

    !print*,'Error in equations:dppdz_r'

    !2000 continue

  end subroutine dppdz_r_f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  subroutine dppdz_i_f(sx, sy, sz2, spr, spi, sgam, sZ, &
                       sdpi, qOK)

    implicit none


    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)
    real(kind=wp), intent(in) :: sZ
    real(kind=wp), contiguous, intent(out) :: sdpi(:)

    real(kind=wp) :: szt

    logical, intent(inout) :: qOK

    LOGICAL :: qOKL


    qOK = .false.

!$OMP WORKSHARE
    sdpi = sInv2rho * (  n2col * bxu  &
           - sEta_G * sp2 / sKappa_G**2 * &
                        sField4ElecImag ) &
           - sKappa_G * spr / sgam * (1 + sEta_G * sp2) &
               * n2col * bzu
!$OMP END WORKSHARE

! Set the error flag and exit

    qOK = .true.

    !goto 2000

    !1000 call Error_log('Error in equations:dppdz_i',tErrorLog_G)

    !print*,'Error in equations:dppdz_i'

    !2000 continue


  end subroutine dppdz_i_f







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dgamdz_f(sx, sy, sz2, spr, spi, sgam, &
                      sdgam, qOK)

    implicit none


    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)

    real(kind=wp), contiguous, intent(out) :: sdgam(:)
    real(kind=wp) ran_nums(iNumberElectrons_G)
    real(kind=wp) dsig1
    logical, intent(inout) :: qOK

    LOGICAL :: qOKL

    CALL RANDOM_NUMBER(ran_nums)
    ran_nums=(2.0_WP*ran_nums)-1.0_WP
    qOK = .false.
    dsig1=1.015E-27*(6.283_WP*sAw_G)**2_WP*(1.697_WP*sAw_G+(1_WP/(1_WP+1.88_WP*sAw_G+0.8_WP*sAw_G**2_WP)))
!$OMP WORKSHARE

    sdgam = -sRho_G * ( 1_WP + sEta_G * sp2 ) / sgam * 2_WP *   &
           ( spr * sField4ElecReal + spi * sField4ElecImag ) &
           ! The recoil part - substracts energy from electrons
           - (((2.0_WP/3.0_WP)*2.818E-15*((sgam*sGammaR_G)*(6.283_WP/lam_w_g)*saw_G)**2.0_WP)/sGammaR_G) &
           ! Energy spread increase over the FEL length
           + (ran_nums*DSQRT(dsig1*(sgam*sGammaR_G)**4_WP*(6.283_WP*sAw_G*lam_w_g))*DSQRT(3.0_WP))

           !- ((ran_nums)*(DSQRT((sgam*sGammaR_G)**4.0_WP*18.849_WP*(1.015E-27* &
           !((6.283_WP/lam_w_g)*saw_G)**2.0_WP*((1.697_WP*saw_G) + &
           !1.0_WP/(1.0_WP+1.88_WP*saw_G+0.8_WP*saw_G**2.0_WP)))))/sGammaR_G)


!$OMP END WORKSHARE

    ! print *, sStepSize/2.0_WP
    ! print *,'Routine called !!'
    ! print *,sdgam
    ! print *, -(((2.0_WP/3.0_WP)*2.818E-15*((sgam*sGammaR_G)*(6.283_WP/lam_w_g)*saw_G)**2.0_WP)/sGammaR_G)
    ! print *,'I was called'
    ! Set the error flag and exit
    ! print *,- (((2.0_WP/3.0_WP)*2.818E-15*((sgam*sGammaR_G)*(6.283_WP/lam_w_g)*saw_G)**2.0_WP)/sGammaR_G)
    ! print *,((ran_nums)*(DSQRT((sgam*sGammaR_G)**4_wp*18.849_wp*(1.015E-27*((6.283_WP/lam_w_g)*saw_G)**2_wp*((1.697_wp*saw_G) + &
    ! 1_wp/(1_wp+1.88_wp*saw_G+0.8_wp*saw_G**2_wp)))))/sGammaR_G)
    qOK = .true.

    !goto 2000

    !1000 call Error_log('Error in equations:dp2dz',tErrorLog_G)

    !print*,'Error in equations:dp2dz'

    !2000 continue


  end subroutine dgamdz_f







  subroutine dxdz_f(sx, sy, sz2, spr, spi, sgam, &
                    sdx, qOK)

    implicit none

!   Calculate dx/dz
!
!              Arguments:


    real(kind=wp), contiguous,  intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                              spi(:), sgam(:)
    real(kind=wp), contiguous, intent(out) :: sdx(:)


!    real(kind=wp), intent(in) :: sy(:), Lj(:), nd
!    real(kind=wp), intent(inout) :: sb(:)
    logical, intent(inout) :: qOK

!              Local vars

    logical :: qOKL ! Local error flag


    qOK = .false.


!$OMP WORKSHARE

    sdx = 2 * sRho_G * sKappa_G / sqrt(sEta_G) * &
          (1 + sEta_G * sp2) / sgam *  &
          spr

!$OMP END WORKSHARE

!    sdx = spr * Lj / nd


    qOK = .true.

    !goto 2000

    !1000 call Error_log('Error in equations:dxdz',tErrorLog_G)

    !2000 continue


  end subroutine dxdz_f





  subroutine dydz_f(sx, sy, sz2, spr, spi, sgam, &
                    sdy, qOK)

    implicit none

!   Calculate dy/dz
!
!              Arguments:


    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)

    real(kind=wp), contiguous, intent(out) :: sdy(:)

    logical, intent(inout) :: qOK

!              Local vars

    logical :: qOKL ! Local error flag


    qOK = .false.


!$OMP WORKSHARE

    sdy = - 2 * sRho_G * sKappa_G / sqrt(sEta_G) * &
          (1 + sEta_G * sp2) / sgam *  &
          spi

!$OMP END WORKSHARE

!    sdy = - spi * Lj / nd


    qOK = .true.

    !goto 2000

    !1000 call Error_log('Error in equations:dydz',tErrorLog_G)

    !2000 continue


  end subroutine dydz_f




  subroutine dz2dz_f(sx, sy, sz2, spr, spi, sgam, &
                     sdz2, qOK)

    implicit none

!   Calculate dz2/dz
!
!              Arguments:


    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)

    real(kind=wp), contiguous, intent(out) :: sdz2(:)


!    real(kind=wp), intent(in) :: sy(:), Lj(:), nd
!    real(kind=wp), intent(inout) :: sb(:)
    logical, intent(inout) :: qOK

!              Local vars

    logical :: qOKL ! Local error flag


    qOK = .false.

!$OMP WORKSHARE

    sdz2 = sp2

!$OMP END WORKSHARE

    qOK = .true.

    !goto 2000

    !1000 call Error_log('Error in equations:dz2dz',tErrorLog_G)

    !2000 continue

  end subroutine dz2dz_f







  subroutine alct_e_srtcts(ar_sz)

    implicit none

! Allocate the arrays used in the calculation of
! the electron eqns

    integer(kind=ip), intent(in) :: ar_sz

    allocate(sp2(ar_sz), sField4ElecReal(ar_sz), &
             sField4ElecImag(ar_sz))! , Lj(ar_sz))

    allocate(bxu(ar_sz), byu(ar_sz), bzu(ar_sz))

  end subroutine alct_e_srtcts



  subroutine dalct_e_srtcts()

    implicit none

! Allocate the arrays used in the calculation of
! the electron eqns

    deallocate(sp2, sField4ElecReal, &
             sField4ElecImag)! , Lj(ar_sz))

    deallocate(bxu, byu, bzu)

  end subroutine dalct_e_srtcts



  subroutine adjUndPlace(szl)

    real(kind=wp) :: szl

      if (qUndEnds_G) then

        if (szl < 0) then

          print*, 'undulator section not recognised, sz < 0!!'
          stop

        else if (sZl <= sZFS) then

          iUndPlace_G = iUndStart_G

        else if (sZl >= sZFE) then

          iUndPlace_G = iUndEnd_G

        else if ((sZl > sZFS) .and. (sZl < sZFE)) then

          iUndPlace_G = iUndMain_G

        else

          print*, 'undulator section not recognised, sz > sZFE!!'
          stop

        end if

      else

        iUndPlace_G = iUndMain_G

      end if

  end subroutine adjUndPlace






end module equations
