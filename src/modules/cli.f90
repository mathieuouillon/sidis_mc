!==============================================================================
! MODULE: cli_module
! Purpose: Command-line interface parsing and utilities
!==============================================================================
module cli_module
    implicit none

contains

subroutine parse_command_line(num_args, nevent, lund_file, iTg, nkin, e0, iColl, eColl, iFM, user_seed)
    implicit none

    integer, intent(in) :: num_args
    integer, intent(out) :: nevent
    character(len=100), intent(out) :: lund_file
    integer, intent(out) :: iTg, nkin, iColl, iFM
    real(kind=4), intent(out) :: e0, eColl
    integer, intent(inout) :: user_seed

    ! Local variables
    integer :: i_arg, iostat
    character(len=100) :: arg_str, next_arg
    logical :: nevent_set, output_set, target_set, nkin_set, e0_set, iColl_set, eColl_set, iFM_set

    ! Initialize flags
    nevent_set = .false.
    output_set = .false.
    target_set = .false.
    nkin_set = .false.
    e0_set = .false.
    iColl_set = .false.
    eColl_set = .false.
    iFM_set = .false.
    nevent = 0
    lund_file = ''
    iTg = -1
    nkin = 0
    e0 = 0.0
    iColl = 0
    eColl = 0.0
    iFM = 5

    if (num_args == 0) then
        write (*, *) 'Usage: ./monte_carlo_simulation --nevent <number_of_events> --output <lund_output_file> &
           --target <target_type> --nkin <nkin_value> --e0 <electron_energy> --iColl <collider_option> --eColl <collider_energy> &
           --iFM <fermi_motion_option>'
        write (*, *) 'Example: ./monte_carlo_simulation --nevent 10000 --output 120k_D.txt --target 1 --nkin 20000 --e0 10.5 &
         --iColl 1 --eColl 0.5 --iFM 5'
        write (*, *) 'Use --help for detailed information about all options'
        stop 1
    end if

    ! Parse named arguments
    i_arg = 1
    do while (i_arg <= num_args)
        call get_command_argument(i_arg, arg_str)

        select case (trim(arg_str))
        case ('--nevent', '-n')
            if (i_arg == num_args) then
                write (*, *) 'Error: --nevent requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) nevent
            if (iostat /= 0 .or. nevent <= 0) then
                write (*, *) 'Error: Invalid number of events. Must be a positive integer.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            nevent_set = .true.
            i_arg = i_arg + 2

        case ('--output', '-o')
            if (i_arg == num_args) then
                write (*, *) 'Error: --output requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, lund_file)
            if (len_trim(lund_file) == 0) then
                write (*, *) 'Error: Lund file name cannot be empty.'
                stop 1
            end if
            output_set = .true.
            i_arg = i_arg + 2

        case ('--target', '-t')
            if (i_arg == num_args) then
                write (*, *) 'Error: --target requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) iTg
            if (iostat /= 0 .or. iTg < 0 .or. iTg > 15) then
                write (*, *) 'Error: Invalid target type. Must be an integer between 0 and 15.'
                write (*, *) 'Target types: 0->p, 1->2H, 2->3H, 3->3He, 4->4He, 5->6Li,'
                write (*, *) '              6->7Li, 7->C, 8->Al, 9->Fe, 10->Sn, 11->Pb,'
                write (*, *) '              12->Ne, 13->Kr, 14->Xe, 15->Cu'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            target_set = .true.
            i_arg = i_arg + 2

        case ('--nkin')
            if (i_arg == num_args) then
                write (*, *) 'Error: --nkin requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) nkin
            if (iostat /= 0 .or. nkin <= 0) then
                write (*, *) 'Error: Invalid nkin value. Must be a positive integer.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            nkin_set = .true.
            i_arg = i_arg + 2

        case ('--e0')
            if (i_arg == num_args) then
                write (*, *) 'Error: --e0 requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) e0
            if (iostat /= 0 .or. e0 <= 0.0) then
                write (*, *) 'Error: Invalid electron energy. Must be a positive number.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            e0_set = .true.
            i_arg = i_arg + 2

        case ('--iColl')
            if (i_arg == num_args) then
                write (*, *) 'Error: --iColl requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) iColl
            if (iostat /= 0 .or. (iColl /= 0 .and. iColl /= 1)) then
                write (*, *) 'Error: Invalid collider option. Must be 0 (no collider) or 1 (collider).'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            iColl_set = .true.
            i_arg = i_arg + 2

        case ('--eColl')
            if (i_arg == num_args) then
                write (*, *) 'Error: --eColl requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) eColl
            if (iostat /= 0 .or. eColl < 0.0) then
                write (*, *) 'Error: Invalid collider energy. Must be a non-negative number.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            eColl_set = .true.
            i_arg = i_arg + 2

        case ('--iFM')
            if (i_arg == num_args) then
                write (*, *) 'Error: --iFM requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) iFM
            if (iostat /= 0 .or. iFM < 0 .or. iFM > 5) then
                write (*, *) 'Error: Invalid Fermi motion option. Must be an integer between 0 and 5.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            iFM_set = .true.
            i_arg = i_arg + 2

        case ('--seed', '-s')
            if (i_arg == num_args) then
                write (*, *) 'Error: --seed requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) user_seed
            if (iostat /= 0 .or. user_seed < 0) then
                write (*, *) 'Error: Invalid seed value. Must be a non-negative integer.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            i_arg = i_arg + 2

        case ('--help', '-h')
            write (*, *) 'Monte Carlo Simulation Program'
            write (*, *) ''
            write (*, *) 'Usage: ./monte_carlo_simulation [OPTIONS]'
            write (*, *) ''
            write (*, *) 'Required Options:'
            write (*, *) '  --nevent, -n    Number of events to generate (positive integer)'
            write (*, *) '  --output, -o    Output file name for Lund format data'
            write (*, *) '  --target, -t    Target type (integer 0-15)'
            write (*, *) '  --nkin          Number of kinematic iterations (positive integer)'
            write (*, *) '  --e0            Electron energy in GeV (positive number)'
            write (*, *) ''
            write (*, *) 'Optional:'
            write (*, *) '  --help, -h      Show this help message'
            write (*, *) '  --iColl         Collider option: 0 (no collider), 1 (collider) [default: 0]'
            write (*, *) '  --eColl        Collider energy in GeV (non-negative number) [default: 0.0]'
            write (*, *) '  --iFM          Fermi motion option (integer 0-5) [default: 5]'
            write (*, *) '  --seed, -s     Random seed for reproducibility (non-negative integer) [default: time-based]'
            write (*, *) ''
            write (*, *) 'Target Types:'
            write (*, *) '  0  -> p      (proton)'
            write (*, *) '  1  -> 2H     (deuterium)'
            write (*, *) '  2  -> 3H     (tritium)'
            write (*, *) '  3  -> 3He    (helium-3)'
            write (*, *) '  4  -> 4He    (helium-4)'
            write (*, *) '  5  -> 6Li    (lithium-6)'
            write (*, *) '  6  -> 7Li    (lithium-7)'
            write (*, *) '  7  -> C      (carbon)'
            write (*, *) '  8  -> Al     (aluminum)'
            write (*, *) '  9  -> Fe     (iron)'
            write (*, *) '  10 -> Sn     (tin)'
            write (*, *) '  11 -> Pb     (lead)'
            write (*, *) '  12 -> Ne     (neon)'
            write (*, *) '  13 -> Kr     (krypton)'
            write (*, *) '  14 -> Xe     (xenon)'
            write (*, *) '  15 -> Cu     (copper)'
            write (*, *) ''
            write (*, *) 'Fermi Motion Options:'
            write (*, *) '  0 -> No Fermi motion'
            write (*, *) '  1 -> like 4 plus a tail from [2] (deut, C, Al, Fe, Sn, Pb, 4He)'
            write (*, *) '  2 -> Accardi SVG (7Li, C, O, Ne, Al, Ar, Ca, Ni, Cu, Zr, Sn, Pb)'
            write (*, *) '  3 -> Accardi CS (2H, 3He, 4He, C, O, Ca, Fe, Pb)'
            write (*, *) '  4 -> Hard sphere with values from [1]'
            write (*, *) '  5 -> R. Wiringa et al. PRC 89, 024305 (2014)'
            write (*, *) ''
            write (*, *) 'Notes:'
            write (*, *) '  - Fermi motion distributions are limited to 1 GeV nucleons.'
            write (*, *) '  - [1] E. J. Moniz et al. PRL 26, 445 (1971)'
            write (*, *) '  - [2] A. Bodek and J. L. Ritchie PRD 23, 1070 (1981)'
            write (*, *) ''
            write (*, *) 'Examples:'
            write (*, *) '  ./monte_carlo_simulation --nevent 10000 --output 120k_D.txt --target 1 --nkin 20000 --e0 10.5 &
            --iColl 1 --eColl 0.5 --iFM 5'
            write (*, *) '  ./monte_carlo_simulation -n 50000 -o results.txt -t 7 --nkin 15000 --e0 12.0'
            write (*, *) '  ./monte_carlo_simulation --target 0 --nevent 25000 --output proton.txt --nkin 10000 --e0 8.5'
            stop 0

        case default
            write (*, *) 'Error: Unknown argument: ', trim(arg_str)
            write (*, *) 'Use --help for usage information'
            stop 1
        end select
    end do

    ! Check if required arguments were provided
    if (.not. nevent_set) then
        write (*, *) 'Error: --nevent argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. output_set) then
        write (*, *) 'Error: --output argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. target_set) then
        write (*, *) 'Error: --target argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. nkin_set) then
        write (*, *) 'Error: --nkin argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. e0_set) then
        write (*, *) 'Error: --e0 argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if
end subroutine parse_command_line


subroutine get_symbol_name(iTg, symbol_name)
    implicit none
    integer, intent(in) :: iTg
    character(len=10), intent(out) :: symbol_name

    select case (iTg)
    case (0)
        symbol_name = 'p'
    case (1)
        symbol_name = '2H'
    case (2)
        symbol_name = '3H'
    case (3)
        symbol_name = 'He3'
    case (4)
        symbol_name = 'He4'
    case (5)
        symbol_name = 'Li6'
    case (6)
        symbol_name = 'Li7'
    case (7)
        symbol_name = 'C'
    case (8)
        symbol_name = 'Al'
    case (9)
        symbol_name = 'Fe'
    case (10)
        symbol_name = 'Sn'
    case (11)
        symbol_name = 'Pb'
    case (12)
        symbol_name = 'Ne'
    case (13)
        symbol_name = 'Kr'
    case (14)
        symbol_name = 'Xe'
    case (15)
        symbol_name = 'Cu'
    case default
        symbol_name = 'Unknown'
    end select
end subroutine get_symbol_name

end module cli_module
