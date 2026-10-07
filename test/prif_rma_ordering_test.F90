#include "julienne-assert-macros.h"
#include "test-utils.F90"

module prif_rma_ordering_test_m
  ! Tests that coarray accesses issued by one image take effect in program order,
  ! and that image control statements complete outstanding accesses.
  ! These exercise the deferred remote completion of puts in the runtime.
# include "test-uses-alloc.F90"
    use julienne_m, only: test_description_t, test_diagnosis_t, test_result_t, test_t, string_t, usher &
      ,operator(.also.), operator(.equalsExpected.), operator(//)
    use prif
    use, intrinsic :: iso_c_binding, only: c_ptrdiff_t

    implicit none
    private
    public :: prif_rma_ordering_test_t

    type, extends(test_t) :: prif_rma_ordering_test_t
    contains
      procedure, nopass, non_overridable :: subject
      procedure, nopass, non_overridable :: results
    end type

    integer, parameter :: n = 64 ! elements per test array
    integer(c_size_t), parameter :: elem = 8 ! bytes per integer(c_int64_t)

contains
    pure function subject()
        character(len=:), allocatable :: subject
        subject = "PRIF RMA ordering"
    end function

    function results() result(test_results)
        type(test_result_t), allocatable :: test_results(:)
        type(prif_rma_ordering_test_t) prif_rma_ordering_test

        allocate(test_results, source = prif_rma_ordering_test%run([ &
              test_description_t("ordering puts to the same location", usher(check_write_after_write)) &
            , test_description_t("reading after partially overlapping puts", usher(check_read_after_overlapping_writes)) &
            , test_description_t("ordering an atomic after a put to the same location", usher(check_atomic_after_put)) &
            , test_description_t("reading after a strided put", usher(check_get_after_strided_put)) &
            , test_description_t("completing puts at an event post", usher(check_put_then_event)) &
            , test_description_t("completing a put with notify", usher(check_put_with_notify)) &
            , test_description_t("completing puts to every image at sync all", usher(check_puts_to_all_images)) &
        ]))
    end function

    ! integer(c_int64_t) :: a(count)[*], initialized to zero on every image
    subroutine allocate_array(count, handle, a)
        integer, intent(in) :: count
        type(prif_coarray_handle), intent(out) :: handle
        integer(c_int64_t), pointer, intent(out) :: a(:)
        type(c_ptr) :: mem

        call prif_allocate_coarray( &
                [1_c_int64_t], [integer(c_int64_t)::], &
                elem * count, &
                null_final_proc, &
                coarray_handle = handle, &
                allocated_memory = mem)
        call c_f_pointer(mem, a, [count])
        a = 0
        call prif_sync_all()
    end subroutine

    function peer_image() result(peer)
        integer :: peer, me, num_imgs
        call prif_this_image_no_coarray(this_image=me)
        call prif_num_images(num_images=num_imgs)
        peer = mod(me, num_imgs) + 1
    end function

    function check_write_after_write() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: handle
        integer(c_int64_t), pointer :: a(:)
        integer(c_int64_t), target :: val
        integer, parameter :: repetitions = 100
        integer :: peer, k

        diag = .true.
        peer = peer_image()
        call allocate_array(1, handle, a)

        do k = 1, repetitions
          val = k
          call prif_put(peer, handle, 0_c_size_t, c_loc(val), elem)
        end do
        call prif_get(peer, handle, 0_c_size_t, c_loc(val), elem)
        ALSO2(val .equalsExpected. int(repetitions, c_int64_t), "get after repeated puts")

        call prif_sync_all()
        ALSO2(a(1) .equalsExpected. int(repetitions, c_int64_t), "local read after sync all")

        call prif_sync_all()
        call prif_deallocate_coarray(handle)
    end function

    function check_read_after_overlapping_writes() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: handle
        integer(c_int64_t), pointer :: a(:)
        integer(c_int64_t), target :: ones(n), twos(n), result_(2*n), expected(2*n)
        integer :: peer

        diag = .true.
        peer = peer_image()
        call allocate_array(2*n, handle, a)

        ! a(1:n)[peer] = 1; a(n/2+1:n/2+n)[peer] = 2
        ones = 1
        twos = 2
        call prif_put(peer, handle, 0_c_size_t, c_loc(ones), elem * n)
        call prif_put(peer, handle, elem * (n/2), c_loc(twos), elem * n)
        call prif_get(peer, handle, 0_c_size_t, c_loc(result_), elem * 2 * n)

        expected = 0
        expected(1:n/2) = 1
        expected(n/2+1:n/2+n) = 2
        ALSO2(all(result_ == expected) .equalsExpected. .true., "get after overlapping puts")

        call prif_sync_all()
        ALSO2(all(a == expected) .equalsExpected. .true., "local read after sync all")

        call prif_sync_all()
        call prif_deallocate_coarray(handle)
    end function

    function check_atomic_after_put() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: handle
        integer(c_int64_t), pointer :: a(:)
        integer(c_int64_t), target :: val
        integer(PRIF_ATOMIC_INT_KIND) :: old, current
        integer :: peer

        diag = .true.
        peer = peer_image()
        call allocate_array(1, handle, a)

        val = 7
        call prif_put(peer, handle, 0_c_size_t, c_loc(val), elem)
        call prif_atomic_fetch_add(peer, handle, 0_c_size_t, value=1_PRIF_ATOMIC_INT_KIND, old=old)
        ALSO2(old .equalsExpected. 7_PRIF_ATOMIC_INT_KIND, "atomic fetch_add after put")
        call prif_atomic_ref_int(peer, handle, 0_c_size_t, value=current)
        ALSO2(current .equalsExpected. 8_PRIF_ATOMIC_INT_KIND, "atomic ref after fetch_add")

        ! and a put after an atomic
        val = 42
        call prif_put(peer, handle, 0_c_size_t, c_loc(val), elem)
        call prif_atomic_ref_int(peer, handle, 0_c_size_t, value=current)
        ALSO2(current .equalsExpected. 42_PRIF_ATOMIC_INT_KIND, "atomic ref after put")

        call prif_sync_all()
        call prif_deallocate_coarray(handle)
    end function

    function check_get_after_strided_put() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: handle
        integer(c_int64_t), pointer :: a(:)
        integer(c_int64_t), target :: src(n), result_(2*n), expected(2*n)
        integer :: peer, i

        diag = .true.
        peer = peer_image()
        call allocate_array(2*n, handle, a)

        ! a(1:2*n:2)[peer] = src, then a(:)[peer] is read back
        src = [(int(i, c_int64_t), i = 1, n)]
        call prif_put_strided(peer, handle, 0_c_size_t, &
                              remote_stride = [2 * int(elem, c_ptrdiff_t)], &
                              current_image_buffer = c_loc(src), &
                              current_image_stride = [int(elem, c_ptrdiff_t)], &
                              element_size = elem, &
                              extent = [int(n, c_size_t)])
        call prif_get(peer, handle, 0_c_size_t, c_loc(result_), elem * 2 * n)

        expected = 0
        expected(1:2*n:2) = src
        ALSO2(all(result_ == expected) .equalsExpected. .true., "get after strided put")

        ! overwrite part of the strided region with a contiguous put, read back strided
        src(1:n/2) = -1
        call prif_put(peer, handle, 0_c_size_t, c_loc(src), elem * (n/2))
        expected(1:n/2) = -1
        call prif_get(peer, handle, 0_c_size_t, c_loc(result_), elem * 2 * n)
        ALSO2(all(result_ == expected) .equalsExpected. .true., "get after strided then contiguous put")

        call prif_sync_all()
        ALSO2(all(a == expected) .equalsExpected. .true., "local read after sync all")

        call prif_sync_all()
        call prif_deallocate_coarray(handle)
    end function

    function check_put_then_event() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: data_handle, event_handle
        integer(c_int64_t), pointer :: a(:)
        type(prif_event_type), pointer :: evt
        type(prif_event_type) :: dummy_event
        type(c_ptr) :: event_mem
        integer(c_int64_t), target :: vals(n)
        integer :: me, peer, i, round
        integer, parameter :: rounds = 10

        diag = .true.
        call prif_this_image_no_coarray(this_image=me)
        peer = peer_image()
        call allocate_array(n, data_handle, a)

        ! type(event_type) :: evt[*]
        call prif_allocate_coarray( &
                [1_c_int64_t], [integer(c_int64_t)::], &
                int(storage_size(dummy_event)/8, c_size_t), &
                null_final_proc, &
                coarray_handle = event_handle, &
                allocated_memory = event_mem)
        call c_f_pointer(event_mem, evt)
        evt = dummy_event
        call prif_sync_all()

        do round = 1, rounds
          vals = [(int(1000*me + 10*round + i, c_int64_t), i = 1, n)]
          call prif_put(peer, data_handle, 0_c_size_t, c_loc(vals), elem * n)
          call prif_event_post(peer, event_handle, 0_c_size_t)
          ! the data from my predecessor must be visible once its post is received
          call prif_event_wait(event_mem)
          block
            integer :: pred, num_imgs
            integer(c_int64_t) :: expected(n)
            call prif_num_images(num_images=num_imgs)
            pred = mod(me - 2 + num_imgs, num_imgs) + 1
            expected = [(int(1000*pred + 10*round + i, c_int64_t), i = 1, n)]
            ALSO2(all(a == expected) .equalsExpected. .true., "data visible after event wait")
          end block
          call prif_sync_all() ! the predecessor must not overwrite a(:) before it is checked
        end do

        call prif_deallocate_coarrays([data_handle, event_handle])
    end function

    function check_put_with_notify() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: data_handle, notify_handle
        integer(c_int64_t), pointer :: a(:)
        type(prif_notify_type), pointer :: ntf
        type(prif_notify_type) :: dummy_notify
        type(c_ptr) :: notify_mem
        integer(c_int64_t), target :: vals(n)
        integer :: me, peer, pred, num_imgs, i

        diag = .true.
        call prif_this_image_no_coarray(this_image=me)
        call prif_num_images(num_images=num_imgs)
        peer = peer_image()
        pred = mod(me - 2 + num_imgs, num_imgs) + 1
        call allocate_array(n, data_handle, a)

        ! type(notify_type) :: ntf[*]
        call prif_allocate_coarray( &
                [1_c_int64_t], [integer(c_int64_t)::], &
                int(storage_size(dummy_notify)/8, c_size_t), &
                null_final_proc, &
                coarray_handle = notify_handle, &
                allocated_memory = notify_mem)
        call c_f_pointer(notify_mem, ntf)
        ntf = dummy_notify
        call prif_sync_all()

        vals = [(int(100*me + i, c_int64_t), i = 1, n)]
        call prif_put_with_notify(peer, data_handle, 0_c_size_t, c_loc(vals), elem * n, &
                                  notify_handle, 0_c_size_t)
        call prif_notify_wait(notify_mem)
        vals = [(int(100*pred + i, c_int64_t), i = 1, n)]
        ALSO2(all(a == vals) .equalsExpected. .true., "data visible after notify wait")

        call prif_sync_all()
        call prif_deallocate_coarrays([data_handle, notify_handle])
    end function

    function check_puts_to_all_images() result(diag)
        type(test_diagnosis_t) :: diag
        type(prif_coarray_handle) :: handle
        integer(c_int64_t), pointer :: a(:)
        integer(c_int64_t), target :: val
        integer :: me, num_imgs, img, i

        diag = .true.
        call prif_this_image_no_coarray(this_image=me)
        call prif_num_images(num_images=num_imgs)
        call allocate_array(num_imgs, handle, a)

        ! a(me)[img] = me for every image img: one non-overlapping put per target
        val = me
        do img = 1, num_imgs
          call prif_put(img, handle, elem * (me - 1), c_loc(val), elem)
        end do
        call prif_sync_all()
        ALSO2(all(a == [(int(i, c_int64_t), i = 1, num_imgs)]) .equalsExpected. .true., "local read after sync all")

        call prif_sync_all()
        call prif_deallocate_coarray(handle)
    end function

end module prif_rma_ordering_test_m
