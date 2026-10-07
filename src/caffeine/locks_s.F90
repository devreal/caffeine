! Copyright (c), The Regents of the University of California
! Terms of use are as specified in LICENSE.txt

#include "assert_macros.h"

submodule(prif:prif_private_s) locks_s
  ! DO NOT ADD USE STATEMENTS HERE
  ! All use statements belong in prif_private_s.F90
  implicit none

contains

  module procedure prif_lock
    integer(c_intptr_t) :: remote_base

    call_assert(offset >= 0)

    call base_pointer(coarray_handle, image_num, remote_base)
    call prif_lock_indirect( &
        image_num = image_num, &
        lock_var_ptr = remote_base + offset, &
        acquired_lock = acquired_lock, &
        stat = stat, errmsg = errmsg, errmsg_alloc = errmsg_alloc)
  end procedure

  module procedure prif_lock_indirect
    integer(c_int) :: rc

    call_assert(prif_init_called_previously)
    call_assert_describe(image_num > 0 .and. image_num <= initial_team%num_images, "image_num not within valid range")

    rc = caf_lock_acquire(image_num, lock_var_ptr, try_only=present(acquired_lock))
    select case (rc)
    case (CAF_LOCK_ACQUIRED)
      if (present(acquired_lock)) acquired_lock = .true.
    case (CAF_LOCK_BUSY)
      call_assert(present(acquired_lock))
      acquired_lock = .false.
    case default
      call_assert(rc == CAF_LOCK_HELD_BY_ME)
      if (present(acquired_lock)) acquired_lock = .false.
      call report_error(PRIF_STAT_LOCKED, "LOCK: the lock variable is already locked by the executing image", &
                        stat, errmsg, errmsg_alloc)
      return
    end select

    if (present(stat)) stat = 0
  end procedure

  module procedure prif_unlock
    integer(c_intptr_t) :: remote_base

    call_assert(offset >= 0)

    call base_pointer(coarray_handle, image_num, remote_base)
    call prif_unlock_indirect( &
        image_num = image_num, &
        lock_var_ptr = remote_base + offset, &
        stat = stat, errmsg = errmsg, errmsg_alloc = errmsg_alloc)
  end procedure

  module procedure prif_unlock_indirect
    integer(c_int) :: rc

    call_assert(prif_init_called_previously)
    call_assert_describe(image_num > 0 .and. image_num <= initial_team%num_images, "image_num not within valid range")

    rc = caf_lock_release(image_num, lock_var_ptr)
    select case (rc)
    case (CAF_UNLOCK_OK)
      if (present(stat)) stat = 0
    case (CAF_UNLOCK_NOT_LOCKED)
      call report_error(PRIF_STAT_UNLOCKED, "UNLOCK: the lock variable is not locked", &
                        stat, errmsg, errmsg_alloc)
    case default
      call_assert(rc == CAF_UNLOCK_OTHER_OWNER)
      call report_error(PRIF_STAT_LOCKED_OTHER_IMAGE, "UNLOCK: the lock variable is locked by another image", &
                        stat, errmsg, errmsg_alloc)
    end select
  end procedure

end submodule locks_s
