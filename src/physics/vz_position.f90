! ------------------------------------------------------------------------------
! Set target vertex z for Lund (RG-D; 2024)
! ------------------------------------------------------------------------------
subroutine calculate_vz_position(target_type, vz)
    implicit none

    integer, intent(in) :: target_type    ! Input target type
    real, intent(out) :: vz               ! Output vertex z position

    ! Local variables
    real, parameter :: first_pos = -2.5
    real, parameter :: second_pos = -7.5
    real :: target_pos
    real :: target_length
    real :: random_number

    random_number = rand()                 ! Choose a value in [0,1]
    vz = 0.0                              ! Default value for vz

    select case (target_type)
    case (10)
        ! Case for Sn
        target_pos = -3.5
        target_length = 0.018

    case (15)
        ! Case for Cu
        target_pos = -8.5
        target_length = 0.009

    case (7)
        ! Case for Carbon
        target_length = 0.2
        if (random_number < 0.5) then
            target_pos = -3.5
        else
            target_pos = -8.5
        end if

    case (1)
        ! Case for LD2
        target_pos = -5.0
        target_length = 5.0

    case default
        ! Do nothing for unspecified targets
        return
    end select

    vz = target_pos + target_length*(rand() - 0.5)

end subroutine calculate_vz_position
