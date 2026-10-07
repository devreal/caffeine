! Copyright (c), The Regents of the University of California
! Terms of use are as specified in LICENSE.txt

#include "assert_macros.h"

submodule(prif:prif_private_s) critical_s
  ! DO NOT ADD USE STATEMENTS HERE
  ! All use statements belong in prif_private_s.F90
  implicit none

  ! The critical construct is implemented as a lock on the prif_critical_type
  ! coarray element residing on image 1 of the initial team
  integer(c_int), parameter :: critical_image = 1

contains

  module procedure prif_critical
    integer(c_intptr_t) :: lock_ptr
    integer(c_int) :: rc

    call_assert(coarray_handle_check(critical_coarray))

    call base_pointer(critical_coarray, critical_image, lock_ptr)
    rc = caf_lock_acquire(critical_image, lock_ptr, try_only=.false.)
    if (rc /= CAF_LOCK_ACQUIRED) then
      call report_error(PRIF_STAT_LOCKED, "CRITICAL: the critical construct is already active on the executing image", &
                        stat, errmsg, errmsg_alloc)
      return
    end if

    if (present(stat)) stat = 0
  end procedure

  module procedure prif_end_critical
    integer(c_intptr_t) :: lock_ptr
    integer(c_int) :: rc

    call_assert(coarray_handle_check(critical_coarray))

    call base_pointer(critical_coarray, critical_image, lock_ptr)
    rc = caf_lock_release(critical_image, lock_ptr)
    call_assert_describe(rc == CAF_UNLOCK_OK, "END CRITICAL: the critical construct was not active on the executing image")
  end procedure

end submodule critical_s
