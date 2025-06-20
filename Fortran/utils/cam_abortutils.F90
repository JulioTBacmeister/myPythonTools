module cam_abortutils

   !use shr_kind_mod, only: SHR_KIND_CL
   !use shr_sys_mod,  only: endrun => shr_sys_abort

   implicit none
   private
   save

   public :: endrun
   !public :: handle_allocate_error

CONTAINS

   !!subroutine handle_allocate_error(retval, subname, fieldname)
   subroutine endrun(retval, subname, fieldname)
      ! if <retval> is not zero, generate an error message and abort
      ! Dummy arguments
      integer,          intent(in) :: retval
      character(len=*), intent(in) :: subname
      character(len=*), intent(in) :: fieldname
      ! Local variable
      character(len=40)   :: errmsg

      if (retval /= 0) then
         write(errmsg, '(4a,i0)') trim(subname), ' error allocating ',        &
              trim(fieldname), ', error = ', retval
         write(*,*) errmsg
         STOP
      end if
   end subroutine endrun

end module cam_abortutils
