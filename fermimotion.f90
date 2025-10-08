!------------------------------------------------------------------------------
! Initialize fermi momentum
!------------------------------------------------------------------------------
subroutine InitNucl
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    ! Include all the common blocks
    ! include 'common.f90'

    ! Fill rFM, iZ and iA in function of the Target
    select case (iTg)
    case (0)
        rFM = 0.0
        iZ = 1
        iA = 1
        ! I(RD) extrapolate rFM values between 6Li and 2H
    case (1)
        rFM = 0.07
        iZ = 1
        iA = 2
    case (2)
        rFM = 0.1
        iZ = 1
        iA = 3
    case (3)
        rFM = 0.1
        iZ = 2
        iA = 3
    case (4)
        rFM = 0.120
        iZ = 2
        iA = 4
    case (5)
        rFM = 0.17
        iZ = 3
        iA = 6
    case (6)
        rFM = 0.17
        iZ = 3
        iA = 7
    case (7)
        rFM = 0.221
        iZ = 6
        iA = 12
    case (8)
        rFM = 0.235
        iZ = 13
        iA = 27
    case (9)
        rFM = 0.260
        iZ = 26
        iA = 56
    case (10)
        rFM = 0.260
        iZ = 50
        iA = 120
    case (11)
        rFM = 0.265
        iZ = 82
        iA = 208
    case (12)
        rFM = 0.228
        iZ = 10
        iA = 20
    case (13)
        rFM = 0.260
        iZ = 36
        iA = 84
    case (14)
        rFM = 0.260
        iZ = 54
        iA = 132
    case default
        rFM = 0.0
        iZ = 1
        iA = 1
    end select

end subroutine InitNucl

subroutine InitFM
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    ! Include all the common blocks
    ! include 'common.f90'
    integer :: irho ! dummy variable

    ! Produce the table for FM generation (CS)
    if (iFM == 2 .or. iFM == 3) then
        irho = iFM - 1
        call GenFMtable(irho)
        ! Produce the table for FM generation (RW)
    else if (iFM == 5) then
        call GenRWtable()
    end if

end subroutine InitFM

! Read table from RW
subroutine GenRWtable()
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    ! include 'common.f90'

    integer :: ERROR
    real :: a, b, c, d, e
    real :: sum_n, sum_p
    integer :: i, j
    character(len=100) :: str

    a = 0.0
    b = 0.0
    c = 0.0
    d = 0.0
    e = 0.0
    write (*, *) 'Enter GenRWtable'

    ERROR = 0
    select case (iA)
    case (2)
        if (iZ == 1) then
            open (unit=8, file='../datafiles/fmrw/h2.momentum', status='old')
        else
            ERROR = 1
        end if
    case (3)
        if (iZ == 1) then
            open (unit=8, file='../datafiles/fmrw/h3.momentum', status='old')
        else if (iZ == 2) then
            open (unit=8, file='../datafiles/fmrw/he3.momentum', status='old')
        else
            ERROR = 1
        end if
    case (4)
        if (iZ == 2) then
            open (unit=8, file='../datafiles/fmrw/he4.momentum', status='old')
        else
            ERROR = 1
        end if
    case (6)
        if (iZ == 3) then
            open (unit=8, file='../datafiles/fmrw/lad.momentum', status='old')
        else
            ERROR = 1
        end if
    case (7)
        if (iZ == 3) then
            open (unit=8, file='../datafiles/fmrw/lat.momentum', status='old')
        else
            ERROR = 1
        end if
    case default
        ERROR = 1
    end select

    if (ERROR == 1) then
        write (*, *) 'ERROR'
        close (8)
        stop
    end if

    ! Skip header lines until we find the marker
    do while (str(2:2) /= '*' .or. str(3:3) /= '*')
        read (8, *) str
        write (*, *) str
        write (*, *) str(2:2)
        write (*, *) str(3:3)
    end do

    i = 0
    sum_n = 0.0
    sum_p = 0.0

    if (iZ == 1 .and. iA == 2) then
        do while (a < FMlimit/0.1973269602)
            i = i + 1
            read (8, *) a, b, c, d
            write (*, *) a, b, c, d
            FM_n(i) = d
            FM_p(i) = d
            sum_n = sum_n + d
            sum_p = sum_p + d
        end do
    else if ((iZ == 1 .and. iA == 3) .or. (iZ == 2 .and. iA == 3) &
             .or. (iZ == 3 .and. iA == 7)) then
        do while (a < FMlimit/0.1973269602)
            i = i + 1
            read (8, *) a, b, c, d, e
            write (*, *) a, b, c, d, e
            FM_n(i) = d
            FM_p(i) = b
            sum_n = sum_n + d
            sum_p = sum_p + b
        end do
    else if (iZ == 3 .and. iA == 6) then
        do while (a < FMlimit/0.1973269602)
            i = i + 1
            read (8, *) a, b, c
            write (*, *) a, b, c
            FM_n(i) = b
            FM_p(i) = b
            sum_n = sum_n + b
            sum_p = sum_p + b
        end do
    else if (iZ == 2 .and. iA == 4) then
        do while (a < FMlimit/0.1973269602)
            i = i + 1
            read (8, *) a, b, c, d, e
            write (*, *) a, b, c, d, e
            FM_n(i) = b*(a + 0.05)**2
            FM_p(i) = b*(a + 0.05)**2
            sum_n = sum_n + FM_n(i)
            sum_p = sum_p + FM_p(i)
            FM_i(i) = d/b
        end do
    else
        write (*, *) 'ERROR'
        close (8)
        stop
    end if

    step_size_FM = FMlimit/real(i)
    write (*, *) 'step ', step_size_FM, sum_n

    do j = 1, i
        FM_n(j) = FM_n(j)/sum_n
        FM_p(j) = FM_p(j)/sum_p
    end do

    do j = 2, i
        FM_n(j) = FM_n(j - 1) + FM_n(j)
        FM_p(j) = FM_p(j - 1) + FM_p(j)
        write (*, *) j, FM_n(j), FM_p(j), FM_i(j)
    end do

    write (*, *) 'Exit GenRWtable'
    close (8)

end subroutine GenRWtable

!------------------------------------------------------------------------------
! Initialize kinematics values with fermi motion
!------------------------------------------------------------------------------
subroutine FMParam
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    real, parameter :: a = 2.0
    real :: R

    ! Important variables for simulation
    real :: Ps, Thr, C, Rd ! variables for Kf computation
    integer :: i

    ! Include all the common blocks
    ! include 'common.f90'

    FMintact = 1

    ! Generate Theta and Phi of the particle's fermi momentum
    ThFM = acos(2.0*ranf(0) - 1.0)
    PhiFM = 2.0*pi*ranf(0)

    ! Selector for the FM
    if (iFM == 4) then
        ! Generation of FM in a Fermi sphere
        Kf = rFM*(ranf(0))**(1.0/3.0)

    else if (iFM == 1) then
        ! Thresholds
        Thr = 1.0 - 6.0*(rFM*a/pi)**2
        Ps = 5.0
        ! Constant calculation
        C = 4.0/3.0*pi*rFM**3

        ! Remove events with Pf > 4 GeV/c or negative values
        do while (Ps > FMlimit .or. Ps < 0.0)
            ! Generation of random number
            Kf = ranf(0)
            ! Apply the threshold and produce the tail
            if (Kf <= Thr) then
                Kf = (3.0*C*Kf/4.0/pi/Thr)**(1.0/3.0)
            else
                R = 1.0/(1.0 - rFM/4.0)
                Kf = -rFM/((Kf - Thr)*pi*C/8.0/rFM**5/a**2 - 1.0)
            end if
            Ps = Kf
        end do

        ! Fermi Momentum from Accardi routines
    else if (iFM == 2 .or. iFM == 3) then
        Rd = ranf(0)
        i = 1
        do while (FM_table(i) < Rd)
            i = i + 1
        end do

        Kf = real(i - 1)*step_size_FM + ranf(0)*step_size_FM

        ! Fermi Momentum from R. Wiringa
    else if (iFM == 5) then
        Rd = ranf(0)
        i = 1
        if (nucleon == 2112) then
            do while (FM_n(i) < Rd)
                i = i + 1
            end do
        else if (nucleon == 2212) then
            do while (FM_p(i) < Rd)
                i = i + 1
            end do
        end if

        Kf = real(i - 1)*step_size_FM + ranf(0)*step_size_FM
        if (iTg /= 1) FMintact = FM_i(i)

        if (Kf > FMlimit) then
            write (*, *) 'warning ', i, step_size_FM
            stop
        end if

    else
        write (*, *) 'this iFM is not implemented'
        stop
    end if

end subroutine FMParam

subroutine GenFMtable(irho)
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    integer, intent(in) :: irho
    integer :: i
    integer :: itz, ita
    real(kind=8) :: rhofermi, mom, proba, ptot, step

    ! include 'common.f90'

    itz = iZ
    ita = iA

    if (iTg == 12 .and. iFM == 3) then
        itz = 8
        ita = 16
    else if (iFM == 3 .and. (iTg == 13 .or. iTg == 14)) then
        itz = 26
        ita = 56
    end if

    FMnb = 1000
    step_size_FM = FMlimit/FMnb
    mom = 0.0d0
    proba = 0.0d0
    ptot = 0.0d0

    do i = 1, FMnb
        proba = 3.1415926d0*mom*mom*rhofermi(itz, ita, mom, irho)*step_size_FM
        ptot = ptot + proba
        FM_table(i) = ptot
        mom = mom + step_size_FM
    end do

    do i = 1, FMnb
        FM_table(i) = FM_table(i)/ptot
    end do

end subroutine GenFMtable

!***********************************************************************
!  Subroutines for nucleon Fermi motion                                *
!  programmer: Alberto Accardi                                         *
!  created: 12 Jun 2005                                                *
!                                                                      *
!  SUMMARY:                                                            *
!                                                                      *
!   I. Single nucleon Fermi distributions - main front-end             *
!  II. SVG distributions (Sandel-Vary-Garpman)                         *
! III. CS distributions (Ciofi-Simula)                                 *
!                                                                      *
!  HISTORY:                                                            *
!                                                                      *
!  * fermimotion   : the _first_ one                                   *
!    (12 Jun 05)    --- WORKING VERSION ---                            *
!  * fermimotion2  : Extracts soft CS (normalized and unnormalized     *
!    (31 Aug 06)     and the unnormalized soft SVG                     *
!                                                                      *
!  NEEDS:                                                              *
!   - fermimotion2.SVG.tbl                                             *
!   - fermimotion2.CS.tbl                                              *
!                                                                      *
!  TO DO LIST                                                          *
!  ----------                                                          *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
!   I. Single nucleon Fermi momentum distributions                     *
!                                                                      *
!***********************************************************************

!***********************************************************************
!     Nucleon Fermi momentum distribution [GeV^-3]
function rhofermi(iZ, iA, k, irho)
    !     programmer: Alberto Accardi
    !     date: 25 May 2005
    !
    !  A. COMMENTARY
    !
    !     Returns the nucleon Fermi momentum distribution in a
    !     nucleus (iZ,iA) normalized to 1.
    !
    !       rhofermi [GeV^-3] = (dp) Fermi momentum distribution
    !                                normalized to 1
    !       iZ,iA             = (i)  Atomic and mass numbers
    !       k [GeV]           = (dp) nucleon momentum
    !       irho              = (i)  distribution model
    !                                1=SVG [1]  2=CS [3]
    !
    !     SVG model
    !     ---------
    !
    !     Short range correations (hard part) fitted to computations in [2].
    !     Long range correlations (soft part) taken from [1], and
    !     normalized so that Int(rho)=1.
    !
    !       rho = N*rho_soft + rho_hard
    !
    !     Where rho_soft = rhosoftSVG(p) and rho_hard = rhohardSVG(p)
    !     have [fm^3] dimension and argument p with [fm^-1] dimension
    !
    !     It allows any A>=7. However:
    !     1) The soft part is parametrized only for
    !        a subset of nuclei (7Li, 12C, 160, 20Ne, 27Al, 40Ar, 40Ca, 58Ni,
    !        63Cu, 90Zr, 120Sn, 208Pb). For other nuclei the one in the above
    !        list which is closest in A is used (a warning is printed on screen).
    !     2) The hard tail (short range correlations) is taken from the
    !        parametrization in Ref. [2]. The slope is fixed to 0.22 fm^2.
    !        The height presented in [2] for a few nuclei (2H 3He 12C 16O
    !        40Ca 56Fe an 208Pb - see coments below for 8Be) has been fitted
    !        to a functional dependence on A.
    !        The fit is accurate to the +-3% level, except for 4He, which
    !        is underestimated by ~30%. In this case, the original
    !        value from [1] is used. Caution should be used when considering
    !        light nuclei, in general.
    !
    !     CS model
    !     --------
    !
    !     It is a parametrization taken from [2] of theoretical computation
    !     for a few nuclei (2H 3He 4He 12C 16O 40Ca 56Fe and 208Pb) made
    !     by many authors.
    !     Parameters for 8Be fitted to results in [3] by A.A.
    !     It allows use of the above listed nuclei only.
    !
    !     References:
    !     [1] M.Sandel, J.P.Vary and S.I.A.Garpman, Phys.Rev.C20(1979)744
    !     [2] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
    !     [3] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
    !
    !  B. DECLARATIONS
    !
    implicit none

    real(kind=8) :: rhofermi, k
    integer, intent(in) :: iZ, iA, irho

    ! *** variables
    real(kind=8) :: kk

    ! *** functions
    real(kind=8) :: rhoSVGsh, rhoCSsh

    ! *** constants & initial values
    ! Conversion factors
    real(kind=8), parameter :: FmGeV = 0.1973269602d0
    real(kind=8), parameter :: Fm3GeV3 = FmGeV*FmGeV*FmGeV
    real(kind=8), parameter :: pi = 3.1415926535898d0
    real(kind=8), parameter :: eightpi = 8.0d0*pi

    !
    !  C. ACTION
    !

    kk = k/FmGeV

    if (irho == 1) then
        rhofermi = rhoSVGsh(iZ, iA, kk, 1)/Fm3GeV3
    else if (irho == 2) then
        rhofermi = rhoCSsh(iZ, iA, kk, 1)/Fm3GeV3
    end if

end function rhofermi

!***********************************************************************
!                                                                      *
!  II. SVG distributions                                               *
!                                                                      *
!***********************************************************************

!***********************************************************************
!     SVG parametrization of Fermi momentum distribution [fm^3]
function rhoSVGsh(iZ, iA, k, idist)
    !     programmer: Alberto Accardi
    !     date: 31 Aug 2006
    !
    !  A. COMMENTARY
    !
    !     Front-end which retutrns "soft" and "hard" SVG distributions.
    !
    !     Short range correations (hard part) fitted to computations in [2].
    !     Long range correlations (soft part) taken from [1], and
    !     normalized so that Int(rho)=1.
    !
    !       rho = N*rho_soft + rho_hard
    !
    !     Where rho_soft = rhosoftSVG(p) and rho_hard = rhohardSVG(p)
    !     have [fm^3] dimension and argument p with [fm^-1] dimension
    !
    !     It allows any A>=7. However:
    !     1) The soft part is parametrized only for
    !        a subset of nuclei (7Li, 12C, 160, 20Ne, 27Al, 40Ar, 40Ca, 58Ni,
    !        63Cu, 90Zr, 120Sn, 208Pb). For other nuclei the one in the above
    !        list which is closest in A is used (a warning is printed on screen).
    !     2) The hard tail (short range correlations) is taken from the
    !        parametrization in Ref. [2]. The slope is fixed to 0.22 fm^2.
    !        The height presented in [2] for a few nuclei (2H 3He 12C 16O
    !        40Ca 56Fe an 208Pb - see coments below for 8Be) has been fitted
    !        to a functional dependence on A.
    !        The fit is accurate to the +-3% level, except for 4He, which
    !        is underestimated by ~30%. In this case, the original
    !        value from [1] is used. Caution should be used when considering
    !        light nuclei, in general.
    !
    !     Normalization: \int(d3k rhoCSsh) = 1 (except for idist=2,4)
    !
    !       rhoSVGsh [fm3] = (dp) Fermi momentum distribution
    !                          normalized to 4*pi
    !       iZ,iA          = (i)  Atomic and mass numbers
    !       k [fm^-1]      = (dp) nucleon momentum
    !       idist          = (i)  1 = soft+hard
    !                             2 = soft (unnormalized)
    !                             3 = soft (normalized to 1)
    !                             4 = hard (unnormalized)
    !                            (note: 1=2+4)
    !
    !     References:
    !     [1] M.Sandel, J.P.Vary and S.I.A.Garpman, Phys.Rev.C20(1979)744
    !
    !  B. DECLARATIONS
    !
    implicit none

    real(kind=8) :: rhoSVGsh, k
    integer, intent(in) :: iZ, iA, idist

    ! *** variables
    real(kind=8) :: norm, rs, rh, height, slope

    ! *** functions
    real(kind=8) :: rhosoftSVG, rhohardSVG

    ! *** constants & initial values
    ! Conversion factors
    real(kind=8), parameter :: FmGeV = 0.1973269602d0
    real(kind=8), parameter :: Fm3GeV3 = FmGeV*FmGeV*FmGeV
    real(kind=8), parameter :: pi = 3.1415926535898d0
    real(kind=8), parameter :: eightpi = 8.0d0*pi

    !
    !  C. ACTION
    !

    if (idist == 1) then
        ! ... soft+hard
        rs = rhosoftSVG(iZ, iA, k)
        rh = rhohardSVG(iZ, iA, k, height, slope, norm)
        rhoSVGsh = (norm*rs + rh)
    else if (idist == 2) then
        ! ... soft (unnormalized)
        rs = rhosoftSVG(iZ, iA, k)
        rh = rhohardSVG(iZ, iA, k, height, slope, norm)
        rhoSVGsh = norm*rs
    else if (idist == 3) then
        ! ... soft (normalized)
        rhoSVGsh = rhosoftSVG(iZ, iA, k)
    else if (idist == 4) then
        ! ... hard (unnormalized)
        rhoSVGsh = rhohardSVG(iZ, iA, k, height, slope, norm)
    end if

end function rhoSVGsh

!***********************************************************************
!     Soft Nucleon Fermi momentum distribution [fm^3]
function rhosoftSVG(iZ, iA, k)
    !     programmer: Alberto Accardi
    !     date: 25 May 2005
    !
    !  A. COMMENTARY
    !
    !     Returns the (absolute value) of the soft part of nucleon
    !     Fermi momentum distribution in a nucleus (iZ,iA), according
    !     to Ref.[1]
    !
    !       rhosoftSVG [fm3] = (dp) Fermi momentum distribution
    !                               normalized to 1
    !       iZ,iA            = (i)  Atomic and mass numbers
    !       k [fm^-1]        = (dp) nucleon momentum
    !
    !     NOTE: absolute value of the soft Fermi distribution is taken
    !     because the parametrization chosen in [1] leads to negative
    !     values at very large momenta, well outside the region of
    !     validity of the parametrization itself. At these values of the
    !     nucleon momentum the distribution is nonetheless very small
    !     and completely subdominant compared to the hard tail.
    !
    !     Reference:
    !     [1] M.Sandel, J.P.Vary and S.I.A.Garpman, Phys.Rev.C20(1979)744
    !
    !  B. DECLARATIONS
    !
    implicit none

    real(kind=8) :: rhosoftSVG, k
    integer, intent(in) :: iZ, iA

    ! *** variables
    integer, parameter :: nmax = 20
    integer :: nin, i, j, n
    save n

    real(kind=8) :: a(6, nmax), alpha(nmax), beta(nmax), nj
    integer :: ZZ(nmax), AA(nmax)
    save a, alpha, beta, ZZ, AA

    character(len=100) :: string

    ! *** functions
    integer :: nextunit
    real(kind=8) :: gaussnd

    ! *** constants & initial values
    ! Constants
    real(kind=8), parameter :: pi = 3.1415926535898d0

    ! initialization variables
    logical :: firsttime
    integer :: Zold, Aold, nnuke
    save Zold, Aold, firsttime, nnuke

    data Zold/0/, Aold/0/, firsttime/.true./

    !
    !  C. ACTION
    !

    if (iA < 7) then
        print *
        print *, 'ERROR (rhosoftSVG): called with A<7: ', iA
        print *, ' -- Try using Ciofi-Simula parametrization instead'
        print *
        stop
    end if

    ! *** INITIALIZATION
    if (firsttime) then
        firsttime = .false.
        nin = nextunit()
        open (unit=nin, file='datafiles/fmacc/fermimotion2.SVG.tbl', &
              status='old')
        ! ... skips headers
        do i = 1, 7
            read (nin, *) string
        end do
        ! ... reads values
        j = 0
10      j = j + 1
        read (nin, *) string, ZZ(j), AA(j), a(1, j), a(2, j), a(3, j), &
            a(4, j), a(5, j), a(6, j), alpha(j), beta(j)
        if (ZZ(j) /= 0) goto 10
        nnuke = j - 1
        close (nin)
    end if

    ! ... if new nucleus, locates nearest nucleus
    if ((iZ /= Zold) .or. (iA /= Aold)) then
        Zold = iZ
        Aold = iA
        call ilocatetab(ZZ, nnuke, iZ, j)
        if (j == 0) then
            n = 1
        else if ((ZZ(j) == iZ) .and. (AA(j) == iA)) then
            ! ... if iZ and iA are found in the table OK
            n = j
        else
            ! ... otherwise chooses the nucleus with closest A
            if ((iA - AA(j)) <= (AA(j + 1) - iA)) then
                n = j
            else
                n = j + 1
            end if
        end if
    end if

    ! *** COMPUTES the soft Fermi momentum distribution
    rhosoftSVG = 0.0d0
    do j = 1, 6
        nj = alpha(n)*beta(n)**j
        rhosoftSVG = rhosoftSVG + a(j, n)*(nj/pi)**1.5d0 &
                     *exp(-(nj*k*k))
    end do

    rhosoftSVG = abs(rhosoftSVG)

end function rhosoftSVG

!***********************************************************************
!     Hard nucleon Fermi momentum distribution [fm^3]
function rhohardSVG(iZ, iA, k, height, slope, norm)
    !     programmer: Alberto Accardi
    !     date: 25 May 2005
    !
    !  A. COMMENTARY
    !
    !     Returns the "hard tail" of nucleon Fermi momentum distribution in a
    !     nucleus (iZ,iA), and its height, slope and normalization.
    !
    !       rhofermi [fm^3]   = (dp) Fermi momentum distribution
    !       iZ,iA             = (i)  Atomic and mass numbers
    !       k [fm^-1]         = (dp) nucleon momentum
    !       height            = (dp) [OUT] height of exponential (hard) tail
    !       slope             = (dp) [OUT] slope of exponential (hard) tail
    !       norm              = (dp) [OUT] normalization = int(d3k*rhohard)
    !
    !     Short range correlations (hard tail) are taken from the
    !     parametrization in Ref. [1]. The slope is fixed to 0.22 fm^2.
    !     The height presented in [1] has been fitted to a functional
    !     dependence on A.
    !     The fit included the following nuclei: 2H 3He 12C 16O
    !     40Ca 56Fe an 208Pb, taken from [1] and 8Be fitted by A.Accardi
    !     to the computations of [2]. The fit is accurate to the +-3% level,
    !     except for 4He, which is underestimated by ~30%. In this case,
    !     the original value from [1] is used. Caution should be used when
    !     considering light nuclei, in general.
    !
    !     Reference:
    !     [1] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
    !     [2] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
    !
    !     TO DO LIST
    !     ----------
    !
    !  B. DECLARATIONS
    !
    implicit none

    real(kind=8) :: rhohardSVG, k, height, slope, norm
    integer, intent(in) :: iZ, iA

    ! *** variables
    real(kind=8) :: AA

    ! *** constants & initial values
    real(kind=8), parameter :: pi = 3.1415926535898d0
    real(kind=8), parameter :: fourpi = 4.0d0*pi

    logical :: firsttime
    save firsttime
    data firsttime/.true./

    !
    !  C. ACTION
    !

    if (firsttime) then
        firsttime = .false.
        if (iA < 8) then
            print *
            print *, '**************************************************'
            print *, 'WARNING (rhohardSVG): called with 2<A<8: A=', iA
            print *, ' -- parametrization of hard tail''s height'
            print *, '    should be used with caution'
            print *, '**************************************************'
            print *
        end if
    end if

    if ((iZ == 2) .and. (iA == 3)) then
        slope = 0.234d0
    else
        slope = 0.220d0
    end if

    if ((iZ == 2) .and. (iA == 4)) then
        height = 0.0244d0/fourpi
    else
        AA = real(iA, kind=8)
        height = 0.00623d0/fourpi &
                 *(1.0d0 + 2.69d0*log(AA/2.0d0)**0.695d0/AA**0.133d0)
    end if

    norm = 1.0d0 - height*(pi/slope)**1.5d0

    rhohardSVG = height*exp(-slope*k*k)

end function rhohardSVG

!***********************************************************************
!                                                                      *
!  III. CS distributions                                               *
!                                                                      *
!***********************************************************************

!***********************************************************************
!     Ciofi-Simula parametrization of Fermi momentum distribution [fm^3]
function rhoCS(iZ, iA, k)
    !     programmer: Alberto Accardi
    !     date: 31 Aug 2006
    !
    !  A. COMMENTARY
    !
    !     Front-end to "rhoCSsh" for compatibility with previous versions
    !     and to simplify the syntax. It returns the full soft+hard CS
    !     parametrization taken from [1] of theoretical computation
    !     for a few nuclei (12C 16O 40Ca 56Fe and 208Pb) made by many authors.
    !     Parameters for 8Be fitted to results in [2] by A.A.
    !     It allows use of the above listed nuclei only.
    !
    !     Normalization: \int(d3k rhoCSsh) = 1
    !
    !       rhoCS [fm3] = (dp) Fermi momentum distribution
    !                          normalized to 4*pi
    !       iZ,iA       = (i)  Atomic and mass numbers
    !       k [fm^-1]   = (dp) nucleon momentum
    !
    !     References:
    !     [1] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
    !     [2] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
    !
    !  B. DECLARATIONS
    !
    implicit none

    real(kind=8) :: rhoCS, rhoCSsh, k
    integer, intent(in) :: iZ, iA

    !
    !  C. ACTION
    !

    rhoCS = rhoCSsh(iZ, iA, k, 1)

end function rhoCS

!***********************************************************************
!     Ciofi-Simula parametrization of Fermi momentum distribution [fm^3]
function rhoCSsh(iZ, iA, k, idist)
    !     programmer: Alberto Accardi
    !     date: 25 May 2005
    !           31 Aug 2006
    !
    !  A. COMMENTARY
    !
    !     It is a parametrization taken from [1] of theoretical computation
    !     for a few nuclei (12C 16O 40Ca 56Fe and 208Pb) made by many authors.
    !     Parameters for 8Be fitted to results in [2] by A.A.
    !     It allows use of the above listed nuclei only.
    !
    !     Normalization: \int(d3k rhoCSsh) = 1 (except for idist=2,4)
    !
    !       rhoCS [fm3] = (dp) Fermi momentum distribution
    !                          normalized to 4*pi
    !       iZ,iA       = (i)  Atomic and mass numbers
    !       k [fm^-1]   = (dp) nucleon momentum
    !       idist       = (i)  1 = soft+hard
    !                          2 = soft (unnormalized)
    !                          3 = soft (normalized to 1)
    !                          4 = hard (unnormalized)
    !                          5 = n0, as defined in [1]
    !                          6 = n1, as defined in [1]
    !                          (note: 1=2+4)
    !
    !     References:
    !     [1] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
    !     [2] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
    !
    !  B. DECLARATIONS
    !
    implicit none

    real(kind=8) :: rhoCSsh, k
    integer, intent(in) :: iZ, iA, idist

    ! *** variables
    integer, parameter :: nmax = 20
    integer :: nin, i, j, n
    save n

    real(kind=8) :: a0(nmax), b0(nmax), c0(nmax), d0(nmax), e0(nmax), &
                    f0(nmax), a1(nmax), b1(nmax), b2(nmax), c1(nmax), d1(nmax), &
                    n0, n1, soft, hard, norm, tmp1, k2, k4, k6, k8
    integer :: ZZ(nmax), AA(nmax)
    save a0, b0, c0, d0, e0, f0, a1, b1, b2, c1, d1, ZZ, AA

    character(len=100) :: string

    ! *** functions
    integer :: nextunit

    ! *** constants & initial values
    ! Constants
    real(kind=8), parameter :: pi = 3.1415926535898d0
    real(kind=8), parameter :: fourpi = 4.0d0*pi
    real(kind=8), parameter :: pith = pi*5.5683279968317d0
    real(kind=8), parameter :: normCS = 1.0d0/(4.0d0*pi)

    ! initialization variables
    logical :: firsttime
    integer :: Zold, Aold, nnuke
    save Zold, Aold, firsttime, nnuke

    data Zold/0/, Aold/0/, firsttime/.true./

    !
    !  C. ACTION
    !

    ! *** INITIALIZATION
    if (firsttime) then
        firsttime = .false.
        nin = nextunit()
        open (unit=nin, file='datafiles/fmacc/fermimotion2.CS.tbl', &
              status='old')
        ! ... skips headers
        do i = 1, 9
            read (nin, *) string
        end do
        ! ... reads n0(k) parameters
        j = 0
10      j = j + 1
        read (nin, *) string, ZZ(j), AA(j), a0(j), b0(j), c0(j), d0(j), &
            e0(j), f0(j)
        if (ZZ(j) /= 0) goto 10
        nnuke = j - 1
        ! ... skips text
        do i = 1, 4
            read (nin, *) string
        end do
        ! ... reads n1(k) parameters
        j = 0
20      j = j + 1
        read (nin, *) string, ZZ(j), AA(j), a1(j), b1(j), b2(j), c1(j), d1(j)
        if (ZZ(j) /= 0) goto 20
        nnuke = j - 1
        close (nin)
    end if

    ! *** if new nucleus, locates it in the tables
    if ((iZ /= Zold) .or. (iA /= Aold)) then
        Zold = iZ
        Aold = iA
        call ilocatetab(ZZ, nnuke, iZ, j)
        ! ... if iZ and iA are found in the table OK
        if ((j > 0) .and. (j <= nnuke)) then
            if ((ZZ(j) == iZ) .and. (AA(j) == iA)) then
                n = j
            else if ((ZZ(j - 1) == iZ) .and. (AA(j - 1) == iA)) then
                n = j - 1
            else if ((ZZ(j + 1) == iZ) .and. (AA(j + 1) == iA)) then
                n = j + 1
            else
                ! ... otherwise stops
                print *
                print *, 'ERROR (rhosoftCS): Z,A outside parameter table', iZ, iA
                print *
                stop
            end if
        else
            ! ... otherwise stops
            print *
            print *, 'ERROR (rhosoftCS): Z,A outside parameter table', iZ, iA
            print *
            stop
        end if
    end if

    ! *** COMPUTES the Fermi momentum distribution
    k2 = k*k

    ! ...  n0 (as defined in [1])
    if (iA <= 4) then
        n0 = a0(n)*exp(-b0(n)*k2)/((1.0d0 + c0(n)*k2)**2) &
             + d0(n)*exp(-e0(n)*k2)/((1.0d0 + f0(n)*k2)**2)
    else
        k4 = k2*k2
        k6 = k4*k2
        k8 = k4*k4
        n0 = a0(n)*exp(-b0(n)*k2) &
             *(1.0d0 + c0(n)*k2 + d0(n)*k4 + e0(n)*k6 + f0(n)*k8)
    end if

    ! ... soft, hard, and n1 (as defined in [1])
    tmp1 = a1(n)*exp(-b1(n)*k2)/((1.0d0 + b2(n)*k2)**2)
    hard = c1(n)*exp(-d1(n)*k2)
    soft = n0 + tmp1
    n1 = tmp1 + hard

    ! *** final result
    if (idist == 1) then
        ! ... soft+hard
        rhoCSsh = normCS*(n0 + n1)
    else if (idist == 2) then
        ! ... soft (unnormalized)
        rhoCSsh = normCS*soft
    else if (idist == 3) then
        ! ... soft (normalized to 1)
        norm = 1.0d0 - (c1(n)/fourpi)*(pith/d1(n))
        rhoCSsh = normCS*(soft/norm)
    else if (idist == 4) then
        ! ... hard (unnormalized)
        rhoCSsh = normCS*hard
    else if (idist == 5) then
        ! ... n0 [1]
        rhoCSsh = normCS*n0
    else if (idist == 6) then
        ! ... n1 [1]
        rhoCSsh = normCS*n1
    end if

end function rhoCSsh

!***********************************************************************
!     Search an ordered integer table
subroutine ilocatetab(xx, n, x, j)
    !     from "Numerical Recipes in F77"
    !
    !  A. COMMENTARY
    !
    !     Given an array xx(1:n) and given a value x returns a value j such
    !     that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
    !     either increasing or decreasing. j=0 or j=n is returned to indicate
    !     that x is out of range.
    !
    !  B. DECLARATIONS
    !
    implicit none

    integer, intent(in) :: n, x
    integer, intent(in) :: xx(n)
    integer, intent(out) :: j

    ! *** variables
    integer :: jl, jm, ju

    !
    !  C. ACTION
    !
    jl = 0
    ju = n + 1
10  if (ju - jl > 1) then
        jm = (ju + jl)/2
        if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
            jl = jm
        else
            ju = jm
        end if
        goto 10
    end if

    if (x == xx(1)) then
        j = 1
    else if (x == xx(n)) then
        j = n - 1
    else
        j = jl
    end if

end subroutine ilocatetab

!***********************************************************************
! INPUT/OUTPUT subroutines                                             *
! v 1.0, May 2005                                                      *
! collected or written by A.Accardi                                    *
!                                                                      *
! Contents:                                                            *
!                                                                      *
!  - NextUnit        Returns unallocated i/o unit                      *
!                                                                      *
!***********************************************************************

!***********************************************************************
!     Finds available i/o unit
integer function NextUnit()
    !     Author: CTEQ collab. (taken from "Cteq6Pdf-2004.f")
    !
    !  A. COMMENTARY
    !
    !     Returns an unallocated FORTRAN i/o unit between 10 and 300.
    !
    !  B. DECLARATIONS
    !
    implicit none

    ! *** variables
    logical :: EX
    integer :: N

    !
    !  C. ACTION
    !
    do N = 10, 300
        inquire (unit=N, opened=EX)
        if (.not. EX) then
            NextUnit = N
            return
        end if
    end do
    stop ' There is no available I/O unit. '

end function NextUnit
