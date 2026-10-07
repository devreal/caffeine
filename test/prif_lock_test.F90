#include "julienne-assert-macros.h"
#include "test-utils.F90"

module prif_lock_test_m
# include "test-uses-alloc.F90"
    use julienne_m, only: test_description_t, test_diagnosis_t, test_result_t, test_t, string_t, usher &
      ,operator(.also.), operator(.equalsExpected.), operator(//)
    use prif
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none
    private
    public :: prif_lock_test_t

    type, extends(test_t) :: prif_lock_test_t
    contains
      procedure, nopass, non_overridable :: subject
      procedure, nopass, non_overridable :: results
    end type

    integer, parameter :: iterations = 50

contains
    pure function subject()
        character(len=:), allocatable :: subject
        subject = "PRIF Locks and Critical"
    end function

    function results() result(test_results)
        type(test_result_t), allocatable :: test_results(:)
        type(prif_lock_test_t) prif_lock_test

        allocate(test_results, source = prif_lock_test%run([ &
              test_description_t("mutual exclusion with a contended lock", usher(check_lock_contended)) &
            , test_description_t("reporting lock and unlock error conditions", usher(check_lock_errors)) &
            , test_description_t("mutual exclusion with a critical construct", usher(check_critical)) &
        ]))
    end function

    ! Allocate a scalar coarray of the given size, returning its handle and local address
    subroutine allocate_scalar_coarray(size_in_bytes, handle, local_addr)
        integer(c_size_t), intent(in) :: size_in_bytes
        type(prif_coarray_handle), intent(out) :: handle
        integer(c_intptr_t), intent(out) :: local_addr
        type(c_ptr) :: mem

        call prif_allocate_coarray( &
                [1_c_int64_t], [integer(c_int64_t)::], &
                size_in_bytes, &
                null_final_proc, &
                coarray_handle = handle, &
                allocated_memory = mem)
        local_addr = transfer(mem, local_addr)
    end subroutine

    ! Non-atomic read-modify-write of counter[1], to be protected by mutual exclusion
    subroutine increment_counter(counter_handle)
        type(prif_coarray_handle), intent(in) :: counter_handle
        integer(c_int64_t), target :: val

        call prif_get(1, counter_handle, 0_c_size_t, c_loc(val), c_sizeof(val))
        val = val + 1
        call prif_put(1, counter_handle, 0_c_size_t, c_loc(val), c_sizeof(val))
    end subroutine

    function read_counter(counter_handle) result(val)
        type(prif_coarray_handle), intent(in) :: counter_handle
        integer(c_int64_t), target :: val

        call prif_get(1, counter_handle, 0_c_size_t, c_loc(val), c_sizeof(val))
    end function

    function check_lock_contended() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_lock_type), pointer :: lock_var
        integer(c_int64_t), pointer :: counter
        integer(c_intptr_t), pointer :: published_addr
        type(prif_coarray_handle) :: lock_handle, counter_handle, addr_handle
        integer(c_intptr_t) :: lock_addr, counter_addr, published_local, image1_lock_addr
        integer(c_intptr_t), target :: tmp_addr
        type(prif_lock_type) :: dummy_lock
        integer(c_int64_t) :: dummy_int
        integer :: num_imgs, i

        diag = .true.
        call prif_num_images(num_images=num_imgs)

        ! type(lock_type) :: lock_var[*]; integer(c_int64_t) :: counter[*]
        ! integer(c_intptr_t) :: published_addr[*]
        call allocate_scalar_coarray(int(storage_size(dummy_lock)/8, c_size_t), lock_handle, lock_addr)
        call allocate_scalar_coarray(int(storage_size(dummy_int)/8, c_size_t), counter_handle, counter_addr)
        call allocate_scalar_coarray(int(storage_size(tmp_addr)/8, c_size_t), addr_handle, published_local)
        call c_f_pointer(transfer(lock_addr, c_null_ptr), lock_var)
        call c_f_pointer(transfer(counter_addr, c_null_ptr), counter)
        call c_f_pointer(transfer(published_local, c_null_ptr), published_addr)
        lock_var = dummy_lock ! default initialization
        counter = 0
        published_addr = lock_addr ! publish the address of my lock variable for indirect access
        call prif_sync_all()

        ! retrieve the address of the lock variable on image 1, as a client would
        ! when following a pointer component of a coarray
        call prif_get(1, addr_handle, 0_c_size_t, c_loc(tmp_addr), c_sizeof(tmp_addr))
        image1_lock_addr = tmp_addr

        do i = 1, iterations
          ! alternate between the direct and indirect interfaces
          if (mod(i, 2) == 0) then
            call prif_lock(1, lock_handle, 0_c_size_t)
            call increment_counter(counter_handle)
            call prif_unlock(1, lock_handle, 0_c_size_t)
          else
            call prif_lock_indirect(1, image1_lock_addr)
            call increment_counter(counter_handle)
            call prif_unlock_indirect(1, image1_lock_addr)
          end if
        end do

        call prif_sync_all()
        ALSO2(read_counter(counter_handle) .equalsExpected. int(iterations * num_imgs, c_int64_t), "lock-protected counter")
        call prif_sync_all()

        call prif_deallocate_coarrays([lock_handle, counter_handle, addr_handle])
    end function

    function check_lock_errors() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_lock_type), pointer :: lock_var
        type(prif_coarray_handle) :: lock_handle
        integer(c_intptr_t) :: lock_addr
        type(prif_lock_type) :: dummy_lock
        integer :: me, num_imgs, peer
        integer(c_int) :: stat
        logical(c_bool) :: acquired

        diag = .true.
        call prif_num_images(num_images=num_imgs)
        call prif_this_image_no_coarray(this_image=me)

        call allocate_scalar_coarray(int(storage_size(dummy_lock)/8, c_size_t), lock_handle, lock_addr)
        call c_f_pointer(transfer(lock_addr, c_null_ptr), lock_var)
        lock_var = dummy_lock
        call prif_sync_all()

        ! lock my own lock variable, then try to lock it again
        call prif_lock(me, lock_handle, 0_c_size_t, stat=stat)
        ALSO2(stat .equalsExpected. 0_c_int, "initial lock")
        call prif_lock(me, lock_handle, 0_c_size_t, stat=stat)
        ALSO2(stat .equalsExpected. PRIF_STAT_LOCKED, "relocking a lock held by the executing image")
        call prif_lock(me, lock_handle, 0_c_size_t, acquired_lock=acquired, stat=stat)
        ALSO2(stat .equalsExpected. PRIF_STAT_LOCKED, "relocking with acquired_lock")

        call prif_sync_all()
        if (num_imgs > 1) then
          ! the lock variable of my peer is held by my peer
          peer = mod(me, num_imgs) + 1
          call prif_lock(peer, lock_handle, 0_c_size_t, acquired_lock=acquired, stat=stat)
          ALSO2(stat .equalsExpected. 0_c_int, "lock with acquired_lock of a lock held by another image")
          ALSO2(logical(acquired) .equalsExpected. .false., "acquired_lock for a lock held by another image")
          call prif_unlock(peer, lock_handle, 0_c_size_t, stat=stat)
          ALSO2(stat .equalsExpected. PRIF_STAT_LOCKED_OTHER_IMAGE, "unlocking a lock held by another image")
        end if
        call prif_sync_all()

        call prif_unlock(me, lock_handle, 0_c_size_t, stat=stat)
        ALSO2(stat .equalsExpected. 0_c_int, "unlock")
        call prif_unlock(me, lock_handle, 0_c_size_t, stat=stat)
        ALSO2(stat .equalsExpected. PRIF_STAT_UNLOCKED, "unlocking an unlocked lock")

        ! acquire with acquired_lock on an unlocked lock
        call prif_lock(me, lock_handle, 0_c_size_t, acquired_lock=acquired, stat=stat)
        ALSO2(stat .equalsExpected. 0_c_int, "lock with acquired_lock of an unlocked lock")
        ALSO2(logical(acquired) .equalsExpected. .true., "acquired_lock for an unlocked lock")
        call prif_unlock(me, lock_handle, 0_c_size_t)

        call prif_sync_all()
        call prif_deallocate_coarray(lock_handle)
    end function

    function check_critical() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_critical_type), pointer :: critical_var
        integer(c_int64_t), pointer :: counter
        type(prif_coarray_handle) :: critical_handle, counter_handle
        integer(c_intptr_t) :: critical_addr, counter_addr
        type(prif_critical_type) :: dummy_critical
        integer(c_int64_t) :: dummy_int
        integer :: num_imgs, i

        diag = .true.
        call prif_num_images(num_images=num_imgs)

        call allocate_scalar_coarray(int(storage_size(dummy_critical)/8, c_size_t), critical_handle, critical_addr)
        call allocate_scalar_coarray(int(storage_size(dummy_int)/8, c_size_t), counter_handle, counter_addr)
        call c_f_pointer(transfer(critical_addr, c_null_ptr), critical_var)
        call c_f_pointer(transfer(counter_addr, c_null_ptr), counter)
        critical_var = dummy_critical
        counter = 0
        call prif_sync_all()

        do i = 1, iterations
          call prif_critical(critical_handle)
          call increment_counter(counter_handle)
          call prif_end_critical(critical_handle)
        end do

        call prif_sync_all()
        ALSO2(read_counter(counter_handle) .equalsExpected. int(iterations * num_imgs, c_int64_t), "critical-protected counter")
        call prif_sync_all()

        call prif_deallocate_coarrays([critical_handle, counter_handle])
    end function

end module prif_lock_test_m
