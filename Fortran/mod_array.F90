module mod_array
    implicit none
contains
    subroutine scale_array(arr, factor, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), intent(inout) :: arr(nx, ny)
        real(8), intent(in) :: factor
        integer :: i, j

        do i = 1, nx
            do j = 1, ny
                arr(i, j) = arr(i, j) * factor
            end do
        end do
    end subroutine scale_array
end module mod_array
