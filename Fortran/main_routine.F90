module main_routine_mod
    use mod_array
    implicit none
contains
    subroutine process_array(arr, factor, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), intent(inout) :: arr(nx, ny)
        real(8), intent(in) :: factor

        call scale_array(arr, factor, nx, ny)
    end subroutine process_array
end module main_routine_mod
