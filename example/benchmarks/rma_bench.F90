program rma_bench
  use, intrinsic :: iso_c_binding
  use prif
  implicit none
  integer, parameter :: n = 20000
  integer :: stat, me, np, peer, i
  type(prif_coarray_handle) :: h
  type(c_ptr) :: mem
  integer(c_int64_t), pointer :: a(:)
  integer(c_int64_t), target :: v
  integer(c_int64_t) :: t0, t1, rate
  procedure(prif_coarray_cleanup_interface), pointer :: nofinal => null()

  call prif_init(stat)
  call prif_this_image_no_coarray(this_image=me)
  call prif_num_images(np)
  peer = mod(me, np) + 1
  call prif_allocate_coarray([1_c_int64_t], [int(np,c_int64_t)], int(8*n, c_size_t), nofinal, h, mem)
  call c_f_pointer(mem, a, [n])
  a = 0
  call prif_sync_all()

  call system_clock(t0, rate)
  do i = 1, n
    v = i
    call prif_put(peer, h, int(8*(i-1), c_size_t), c_loc(v), 8_c_size_t)
  end do
  call prif_sync_all()
  call system_clock(t1)
  if (me == 1) print '(a,f8.3,a)', "put, distinct offsets:  ", 1e6*real(t1-t0)/rate/n, " us/op"
  if (any(a /= [(int(i,c_int64_t), i=1,n)])) error stop "put distinct: wrong data"

  call prif_sync_all()
  call system_clock(t0)
  do i = 1, n
    v = i
    call prif_put(peer, h, 0_c_size_t, c_loc(v), 8_c_size_t)
  end do
  call prif_sync_all()
  call system_clock(t1)
  if (me == 1) print '(a,f8.3,a)', "put, same offset:       ", 1e6*real(t1-t0)/rate/n, " us/op"
  if (a(1) /= n) error stop "put same: wrong data"

  call prif_sync_all()
  call system_clock(t0)
  do i = 1, n
    call prif_get(peer, h, int(8*(i-1), c_size_t), c_loc(v), 8_c_size_t)
  end do
  call system_clock(t1)
  if (me == 1) print '(a,f8.3,a)', "get:                    ", 1e6*real(t1-t0)/rate/n, " us/op"

  call prif_sync_all()
  call system_clock(t0)
  do i = 1, n
    v = i
    call prif_put(peer, h, int(8*(i-1), c_size_t), c_loc(v), 8_c_size_t)
    call prif_get(peer, h, int(8*(i-1), c_size_t), c_loc(v), 8_c_size_t)
    if (v /= i) error stop "put/get: read-after-write violated"
  end do
  call prif_sync_all()
  call system_clock(t1)
  if (me == 1) print '(a,f8.3,a)', "put+get same location:  ", 1e6*real(t1-t0)/rate/n, " us/pair"

  call prif_sync_all()
  call system_clock(t0)
  do i = 1, n
    v = i
    call prif_put(me, h, int(8*(i-1), c_size_t), c_loc(v), 8_c_size_t)
    if (a(i) /= i) error stop "self put not visible locally"
  end do
  call system_clock(t1)
  if (me == 1) print '(a,f8.3,a)', "put to self:            ", 1e6*real(t1-t0)/rate/n, " us/op"

  call prif_stop(.true._c_bool)
end program
