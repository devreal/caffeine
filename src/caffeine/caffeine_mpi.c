// Copyright (c), The Regents of the University of California
// Terms of use are as specified in LICENSE.txt
//
// MPI-based communication runtime for Caffeine, selected with -DCAF_RUNTIME_MPI=1
// (./install.sh --runtime=mpi). The GASNet-EX runtime lives in caffeine.c.
// Implements the caf_* C interface used by the Fortran PRIF layer on top of
// MPI-3/4 one-sided communication (RMA) and MPI collectives.
//
// Design overview:
// - One MPI process per image. The initial team uses a private duplicate of
//   MPI_COMM_WORLD so that PRIF traffic is isolated from user MPI traffic.
// - All remotely accessible memory (coarrays and prif_allocate'd memory) lives in
//   a single per-image "segment" created with MPI_Win_allocate over the initial
//   team. A passive-target epoch (MPI_Win_lock_all) is held for the lifetime of
//   the program. Remote addresses are translated to window displacements using
//   the allgathered segment base addresses.
// - Puts only wait for source completion (MPI_Win_flush_local), as permitted by
//   PRIF. The target and the byte range of each incomplete put is recorded
//   ("dirty target tracking"). A later access by this image that overlaps a
//   dirty range on the same target first completes it with MPI_Win_flush,
//   preserving the serial order of this image's accesses. All dirty targets are
//   completed at segment boundaries (image control statements).
// - Gets wait for local completion; atomics are completed at the target.
// - Accesses to the executing image's own memory use memcpy.
// - Teams are MPI communicators created by MPI_Comm_split.

#if CAF_RUNTIME_MPI

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <fenv.h>
#include <time.h>
#include <stdarg.h>
#include <unistd.h>
#include <pthread.h>
#include <mpi.h>
#include <ISO_Fortran_binding.h>
#include "../dlmalloc/dl_malloc_caf.h"
#include "../dlmalloc/dl_malloc.h"
#include "caffeine-internal.h"
#include "caffeine-version.h"

// Ensure assertion enforcement in this file tracks the Caffeine ASSERTIONS setting
#undef NDEBUG
#if !ASSERTIONS
#define NDEBUG 1
#endif
#include <assert.h>

#ifndef CAF_DEBUG_DEFER_PUTS
#define CAF_DEBUG_DEFER_PUTS 0
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

typedef uint8_t byte;

// A team is a communicator plus cached rank translation to the initial team
typedef struct caf_team {
  MPI_Comm comm;
  MPI_Group group;
  int rank, size;
  int *to_initial; // to_initial[team rank] == initial team rank
} caf_team;
typedef caf_team* caf_team_t;

// Client reduction callback, as defined by prif_operation_wrapper_interface
typedef void (*caf_reduce_fn_t)(void *arg1, void *arg2_and_out, size_t count, void *cdata);

static caf_team initial_team;
static MPI_Comm termination_comm = MPI_COMM_NULL; // used only for orderly shutdown
static MPI_Group world_group;
static int myproc, numprocs;
static bool caf_initialized_mpi = false;

static MPI_Win seg_win = MPI_WIN_NULL;
static byte *seg_base;
static size_t seg_size;
static intptr_t *seg_bases; // seg_bases[r] == segment base address on initial-team rank r

static mspace* non_symmetric_heap;
static pthread_mutex_t non_symmetric_heap_lock = PTHREAD_MUTEX_INITIALIZER;

// Dirty target tracking: puts (and deferred event posts) that are not yet
// complete at their target. For each initial-team rank r with incomplete
// operations, [dirty_lo[r], dirty_hi[r]) bounds the window displacements written,
// and r appears in dirty_list[0..num_dirty-1] at position dirty_pos[r].
static MPI_Aint *dirty_lo, *dirty_hi;
static int *dirty_list, *dirty_pos;
static int num_dirty = 0;

// Above this many dirty targets, completion uses one MPI_Win_flush_all
#ifndef CAF_FLUSH_ALL_THRESHOLD
#define CAF_FLUSH_ALL_THRESHOLD 8
#endif

#if CAF_DEBUG_DEFER_PUTS
// Debugging aid that emulates the weakest completion MPI permits:
// contiguous puts to other images are buffered and only issued when their
// target is completed, newest first, so that any missing ordering or
// completion step in the runtime produces stale or reordered data.
typedef struct deferred_put {
  struct deferred_put *next;
  MPI_Aint disp;
  size_t size;
  byte data[];
} deferred_put;
static deferred_put **deferred_puts; // per-rank list, newest first
#endif

static MPI_Op user_reduce_op = MPI_OP_NULL;

static void caf_fatal(const char *fmt, ...) __attribute__((noreturn, format(printf,1,2)));

static void caf_fatal(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  fprintf(stderr, "*** Caffeine FATAL ERROR (image %d): ", myproc + 1);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
  va_end(args);
  fflush(NULL);
  MPI_Abort(MPI_COMM_WORLD, 1);
  abort(); // not reached
}

#define MPI_SAFE(call) do { \
    int _rc = (call); \
    if (_rc != MPI_SUCCESS) { \
      char _msg[MPI_MAX_ERROR_STRING]; int _len = 0; \
      MPI_Error_string(_rc, _msg, &_len); \
      caf_fatal("%s failed at %s:%d: %.*s", #call, __FILE__, __LINE__, _len, _msg); \
    } \
  } while (0)

// ---------------------------------------------------
// Environment helpers

// Parse a memory size with an optional K/M/G/T suffix; bare numbers are in units of `unit`
static size_t getenv_memsize(const char *key, size_t dflt, size_t unit) {
  const char *val = getenv(key);
  if (!val || !*val) return dflt;
  char *end = NULL;
  double num = strtod(val, &end);
  if (end == val || num < 0) caf_fatal("Invalid value for environment variable %s='%s'", key, val);
  while (*end == ' ') end++;
  switch (*end) {
    case 'k': case 'K': unit = (size_t)1 << 10; break;
    case 'm': case 'M': unit = (size_t)1 << 20; break;
    case 'g': case 'G': unit = (size_t)1 << 30; break;
    case 't': case 'T': unit = (size_t)1 << 40; break;
    case 'b': case 'B': unit = 1; break;
    case '\0': break;
    default: caf_fatal("Invalid value for environment variable %s='%s'", key, val);
  }
  return (size_t)(num * unit);
}

static double getenv_double(const char *key, double dflt) {
  const char *val = getenv(key);
  if (!val || !*val) return dflt;
  char *end = NULL;
  double num = strtod(val, &end);
  if (end == val) caf_fatal("Invalid value for environment variable %s='%s'", key, val);
  return num;
}

// ---------------------------------------------------
// Floating-point exception support

#ifndef IEEE_FE_MASK
#define IEEE_FE_MASK   FE_INEXACT
#endif
static fexcept_t fe_flag_save;
void caf_fe_save(void) {
  fegetexceptflag(&fe_flag_save, IEEE_FE_MASK);
}
void caf_fe_restore(void) {
  fesetexceptflag(&fe_flag_save, IEEE_FE_MASK);
}

// ---------------------------------------------------
// Teams

static void team_init(caf_team *t, MPI_Comm comm) {
  t->comm = comm;
  MPI_SAFE(MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN));
  MPI_SAFE(MPI_Comm_rank(comm, &t->rank));
  MPI_SAFE(MPI_Comm_size(comm, &t->size));
  MPI_SAFE(MPI_Comm_group(comm, &t->group));
  t->to_initial = malloc(sizeof(int) * t->size);
  int *ranks = malloc(sizeof(int) * t->size);
  if (!t->to_initial || !ranks) caf_fatal("out of memory in team_init");
  for (int i = 0; i < t->size; i++) ranks[i] = i;
  MPI_SAFE(MPI_Group_translate_ranks(t->group, t->size, ranks, world_group, t->to_initial));
  free(ranks);
}

int caf_this_image(caf_team_t tm) {
  return tm->rank + 1;
}
int caf_num_images(caf_team_t tm) {
  return tm->size;
}

// Given team and corresponding image_num, return image number in the initial team
int caf_image_to_initial(caf_team_t tm, int image_num) {
  assert(image_num >= 1);
  assert(image_num <= tm->size);
  return tm->to_initial[image_num-1] + 1;
}
// Given image number in the initial team, return image number corresponding to given team
int caf_image_from_initial(caf_team_t tm, int image_num) {
  assert(image_num >= 1);
  assert(image_num <= numprocs);
  int world_rank = image_num - 1, team_rank;
  MPI_SAFE(MPI_Group_translate_ranks(world_group, 1, &world_rank, tm->group, &team_rank));
  // MPI_UNDEFINED indicates the provided image_num in initial team is not part of tm
  assert(team_rank != MPI_UNDEFINED);
  return team_rank + 1;
}

void caf_form_team(caf_team_t current_team, caf_team_t* new_team, int64_t team_number, int new_index)
{
  // MPI color argument is a non-negative int, check for value truncation:
  if (team_number < 0 || team_number > INT_MAX)
    caf_fatal("FORM TEAM: team_number %lld is not supported (must be in [1,%d])",
              (long long)team_number, INT_MAX);
  MPI_Comm newcomm;
  MPI_SAFE(MPI_Comm_split(current_team->comm, (int)team_number, new_index, &newcomm));
  caf_team *t = malloc(sizeof(caf_team));
  if (!t) caf_fatal("out of memory in caf_form_team");
  team_init(t, newcomm);
  *new_team = t;
}

// ---------------------------------------------------
// Program startup and shutdown

static void reduce_trampoline(void *invec, void *inoutvec, int *len, MPI_Datatype *dt);

static bool caf_finalized = false;

// Release MPI resources held by Caffeine and finalize MPI if we initialized it
static void caf_finalize(void) {
  if (caf_finalized) return;
  caf_finalized = true;
  MPI_Win_unlock_all(seg_win);
  if (caf_initialized_mpi) MPI_Finalize();
}

// Invoked when the process exits without passing through prif_stop/prif_error_stop,
// e.g. when a program not compiled for multi-image execution reaches END PROGRAM
static void caf_atexit(void) {
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized) return;
  fflush(NULL);
  caf_finalize();
}

void caf_caffeinate(
  intptr_t* total_heap_size,
  mspace* symmetric_heap,
  intptr_t* symmetric_heap_start,
  intptr_t* symmetric_heap_size,
  caf_team_t* initial_team_out
) {
  // The MPI library may raise floating-point exceptions during initialization
  // (observed: FE_OVERFLOW with Open MPI). These must not be visible to the
  // program, so restore the exception flags on exit.
  fexcept_t fe_flags_on_entry;
  fegetexceptflag(&fe_flags_on_entry, FE_ALL_EXCEPT);

  int flag = 0;
  MPI_SAFE(MPI_Initialized(&flag));
  if (!flag) {
    int provided;
    MPI_SAFE(MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided));
    caf_initialized_mpi = true;
  }

  MPI_Comm world;
  MPI_SAFE(MPI_Comm_dup(MPI_COMM_WORLD, &world));
  MPI_SAFE(MPI_Comm_dup(MPI_COMM_WORLD, &termination_comm));
  MPI_SAFE(MPI_Comm_group(world, &world_group));
  MPI_SAFE(MPI_Comm_rank(world, &myproc));
  MPI_SAFE(MPI_Comm_size(world, &numprocs));
  team_init(&initial_team, world);
  *initial_team_out = &initial_team;

  // Segment sizing
  long pagesize_l = sysconf(_SC_PAGESIZE);
  const size_t pagesize = (size_t)(pagesize_l > 0 ? pagesize_l : 4096);
  #define PAGE_ALIGNUP(sz) (((sz) + pagesize - 1) & ~(pagesize - 1))

  #ifndef CAF_DEFAULT_HEAP_SIZE
  #define CAF_DEFAULT_HEAP_SIZE (128*1024*1024) // 128 MiB
  #endif
  // retrieve user preference, defaulting to the above and units of MiB
  size_t segsz = getenv_memsize("CAF_HEAP_SIZE", CAF_DEFAULT_HEAP_SIZE, 1024*1024);
  segsz = MAX(segsz, 2*pagesize);
  segsz = PAGE_ALIGNUP(segsz);
  assert(segsz % pagesize == 0);

  // All images must agree on the segment size, use the minimum request
  unsigned long long segsz_ull = segsz;
  MPI_SAFE(MPI_Allreduce(MPI_IN_PLACE, &segsz_ull, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, world));
  segsz = (size_t)segsz_ull;

  MPI_Info info;
  MPI_SAFE(MPI_Info_create(&info));
  MPI_SAFE(MPI_Info_set(info, "accumulate_ordering", "none"));
  MPI_SAFE(MPI_Info_set(info, "same_disp_unit", "true"));
  void *base = NULL;
  int rc = MPI_Win_allocate((MPI_Aint)segsz, 1, info, world, &base, &seg_win);
  MPI_SAFE(MPI_Info_free(&info));
  if (rc != MPI_SUCCESS)
    caf_fatal("Failed to allocate a shared heap segment of %zu bytes. "
              "Consider setting the CAF_HEAP_SIZE environment variable to request a smaller heap.", segsz);
  MPI_SAFE(MPI_Win_set_errhandler(seg_win, MPI_ERRORS_RETURN));
  seg_base = base;
  seg_size = segsz;

  // Local load/stores to coarray memory must observe the results of RMA
  {
    int *model = NULL;
    int found = 0;
    MPI_SAFE(MPI_Win_get_attr(seg_win, MPI_WIN_MODEL, &model, &found));
    if (!found || *model != MPI_WIN_UNIFIED)
      caf_fatal("The MPI library does not provide the MPI_WIN_UNIFIED memory model required by Caffeine");
  }

  seg_bases = malloc(sizeof(intptr_t) * numprocs);
  dirty_lo = malloc(sizeof(MPI_Aint) * numprocs);
  dirty_hi = malloc(sizeof(MPI_Aint) * numprocs);
  dirty_list = malloc(sizeof(int) * numprocs);
  dirty_pos = malloc(sizeof(int) * numprocs);
  if (!seg_bases || !dirty_lo || !dirty_hi || !dirty_list || !dirty_pos)
    caf_fatal("out of memory allocating segment tables");
  for (int r = 0; r < numprocs; r++) dirty_lo[r] = dirty_hi[r] = 0;
#if CAF_DEBUG_DEFER_PUTS
  deferred_puts = calloc(numprocs, sizeof(deferred_put*));
  if (!deferred_puts) caf_fatal("out of memory allocating deferred put table");
#endif
  intptr_t my_base = (intptr_t)seg_base;
  MPI_SAFE(MPI_Allgather(&my_base, sizeof(intptr_t), MPI_BYTE, seg_bases, sizeof(intptr_t), MPI_BYTE, world));

  MPI_SAFE(MPI_Win_lock_all(MPI_MODE_NOCHECK, seg_win));

  *symmetric_heap_start = (intptr_t)seg_base;
  *total_heap_size = (intptr_t)seg_size;

  #ifndef CAF_DEFAULT_COMP_FRAC
  #define CAF_DEFAULT_COMP_FRAC 0.1f // 10%
  #endif
  float default_comp_frac = MAX(MIN(0.99f, CAF_DEFAULT_COMP_FRAC), 0.01f);
  double non_symmetric_fraction = getenv_double("CAF_COMP_FRAC", default_comp_frac);
  if (non_symmetric_fraction <= 0 || non_symmetric_fraction >= 1) {
    caf_fatal("If used, environment variable 'CAF_COMP_FRAC' must be a valid floating point value or fraction between 0 and 1.");
  }

  size_t non_symmetric_heap_size = (size_t)(seg_size * non_symmetric_fraction);
  non_symmetric_heap_size = PAGE_ALIGNUP(non_symmetric_heap_size);
  if (non_symmetric_heap_size >= seg_size) non_symmetric_heap_size = seg_size - pagesize;
  *symmetric_heap_size = (intptr_t)(seg_size - non_symmetric_heap_size);
  assert(non_symmetric_heap_size > 0);
  assert(non_symmetric_heap_size % pagesize == 0);
  assert(*symmetric_heap_size > 0);
  assert(*symmetric_heap_size % pagesize == 0);
  intptr_t non_symmetric_heap_start = *symmetric_heap_start + *symmetric_heap_size;

  if (myproc == 0) {
    *symmetric_heap = create_mspace_with_base((void*)*symmetric_heap_start, *symmetric_heap_size, 0);
    assert(*symmetric_heap);
    mspace_set_footprint_limit(*symmetric_heap, *symmetric_heap_size);
  } else {
    *symmetric_heap = NULL;
  }
  non_symmetric_heap = create_mspace_with_base((void*)non_symmetric_heap_start, non_symmetric_heap_size, 0);
  assert(non_symmetric_heap);
  mspace_set_footprint_limit(non_symmetric_heap, non_symmetric_heap_size);
  #undef PAGE_ALIGNUP

  MPI_SAFE(MPI_Op_create(reduce_trampoline, 1, &user_reduce_op));

  atexit(caf_atexit);

  // ensure every image has its window epoch established before any RMA
  MPI_SAFE(MPI_Barrier(world));

  fesetexceptflag(&fe_flags_on_entry, FE_ALL_EXCEPT);
}

void caf_acquire_exit_lock() {
  static pthread_mutex_t exit_lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_lock(&exit_lock);
}

static void caf_segment_release(void);
static void complete_all_targets(void);

// Normal termination: invoked once all images have synchronized in prif_stop
void caf_decaffeinate(int exit_code) {
  fflush(NULL);
  caf_segment_release();
  // Wait for every image (including any failed images) to reach termination,
  // so no image tears down its window while RMA may still target it.
  MPI_SAFE(MPI_Barrier(termination_comm));
  caf_finalize();
  exit(exit_code);
}

// Error termination: terminate all images
void caf_abort(int exit_code) {
  fflush(NULL);
  MPI_Abort(MPI_COMM_WORLD, exit_code);
  exit(exit_code); // not reached
}

void caf_fail_image() {
  fprintf(stderr,"FAIL IMAGE on image %d\n", myproc+1);
  fflush(NULL);
  complete_all_targets();

  if (numprocs > 1) {
    // Failed images are not detectable by other images in this implementation.
    // Cease participating in program execution while still servicing RMA
    // (which progresses inside the blocking barrier) until all other images
    // have reached normal termination.
    MPI_Barrier(termination_comm);
  }
  caf_finalize();

  exit(1);
}

void caf_fatal_error( const CFI_cdesc_t* Fstr )
{
  const char *msg = (char *)Fstr->base_addr;
  int len = (int)Fstr->elem_len;
  caf_fatal("%.*s", len, msg);
}

// ---------------------------------------------------
// Memory management

void* caf_allocate(mspace heap, size_t bytes) {
  void* allocated_space = mspace_memalign(heap, 8, bytes);
  return allocated_space;
}

void* caf_allocate_non_symmetric(size_t bytes) {
  pthread_mutex_lock(&non_symmetric_heap_lock);
  void* allocated_space = caf_allocate(non_symmetric_heap, bytes);
  pthread_mutex_unlock(&non_symmetric_heap_lock);
  return allocated_space;
}

void caf_allocate_remaining(mspace heap, void** allocated_space, size_t* allocated_size)
{
  // The following doesn't necessarily give us all remaining space
  // nor necessarily the largest open space, but in practice is likely
  // to work out that way
  struct mallinfo heap_info = mspace_mallinfo(heap);

  // clang's implementation of nearbyint() raises FE_INEXACT,
  // in direct contradiction to its specified purpose.
  // Workaround this defect by saving and restoring the FE flags
  caf_fe_save();
  *allocated_size = (size_t)nearbyint(heap_info.keepcost * 0.9f);
  caf_fe_restore();

  *allocated_space = mspace_memalign(heap, 8, *allocated_size);
  if (!*allocated_space) // uh-oh, something went wrong..
    caf_fatal("caf_allocate_remaining failed to mspace_memalign(%zu)", *allocated_size);
}

void caf_deallocate(mspace heap, void* mem) {
  mspace_free(heap, mem);
}

void caf_deallocate_non_symmetric(void* mem) {
  pthread_mutex_lock(&non_symmetric_heap_lock);
  caf_deallocate(non_symmetric_heap, mem);
  pthread_mutex_unlock(&non_symmetric_heap_lock);
}

void caf_establish_mspace(mspace* heap, void* heap_start, size_t heap_size)
{
  *heap = create_mspace_with_base(heap_start, heap_size, 0);
  mspace_set_footprint_limit(*heap, heap_size);
}

// take address in a segment and convert to an address on given image
intptr_t caf_convert_base_addr(void* addr, int image)
{
  assert(image >= 1 && image <= numprocs);
  ptrdiff_t offset = (byte*)addr - seg_base;
  assert(offset >= 0 && (size_t)offset <= seg_size);
  return seg_bases[image-1] + offset;
}

// Convert an address on initial-team rank `rank` to a displacement in the segment window
static inline MPI_Aint to_disp(int rank, intptr_t addr, size_t len) {
  assert(rank >= 0 && rank < numprocs);
  MPI_Aint disp = (MPI_Aint)(addr - seg_bases[rank]);
  if (disp < 0 || (size_t)disp + len > seg_size)
    caf_fatal("Remote address %p on image %d (length %zu) is not within the shared heap",
              (void*)addr, rank + 1, len);
  return disp;
}

// _______________________ Dirty target tracking ____________________________

static inline bool target_is_dirty(int rank) {
  return dirty_hi[rank] > dirty_lo[rank];
}

// Record that [lo, hi) on rank was written by an operation not yet complete at the target
static inline void mark_dirty(int rank, MPI_Aint lo, MPI_Aint hi) {
  assert(lo < hi);
  if (target_is_dirty(rank)) {
    dirty_lo[rank] = MIN(dirty_lo[rank], lo);
    dirty_hi[rank] = MAX(dirty_hi[rank], hi);
  } else {
    dirty_lo[rank] = lo;
    dirty_hi[rank] = hi;
    dirty_pos[rank] = num_dirty;
    dirty_list[num_dirty++] = rank;
  }
}

static inline void mark_clean(int rank) {
  if (!target_is_dirty(rank)) return;
  dirty_lo[rank] = dirty_hi[rank] = 0;
  // swap-remove rank from dirty_list
  int pos = dirty_pos[rank], last = dirty_list[--num_dirty];
  dirty_list[pos] = last;
  dirty_pos[last] = pos;
}

// Complete all outstanding operations to rank at the target
static inline void flush_target(int rank) {
#if CAF_DEBUG_DEFER_PUTS
  for (deferred_put *p = deferred_puts[rank], *next; p; p = next) {
    next = p->next;
    MPI_SAFE(MPI_Put(p->data, (int)p->size, MPI_BYTE, rank, p->disp, (int)p->size, MPI_BYTE, seg_win));
    MPI_SAFE(MPI_Win_flush(rank, seg_win));
    free(p);
  }
  deferred_puts[rank] = NULL;
#endif
  MPI_SAFE(MPI_Win_flush(rank, seg_win));
  mark_clean(rank);
}

// Before accessing [lo, hi) on rank: complete any incomplete write that may overlap,
// so this image's accesses take effect in program order
static inline void order_access(int rank, MPI_Aint lo, MPI_Aint hi) {
  if (target_is_dirty(rank) && lo < dirty_hi[rank] && dirty_lo[rank] < hi)
    flush_target(rank);
}

// Complete all outstanding operations at all targets
static void complete_all_targets(void) {
  if (num_dirty == 0) return;
  if (num_dirty > CAF_FLUSH_ALL_THRESHOLD && !CAF_DEBUG_DEFER_PUTS) {
    MPI_SAFE(MPI_Win_flush_all(seg_win));
    while (num_dirty > 0) mark_clean(dirty_list[num_dirty-1]);
  } else {
    while (num_dirty > 0) flush_target(dirty_list[num_dirty-1]);
  }
}

// _______________________ Contiguous RMA ____________________________

// largest single transfer issued in one MPI call
#define CAF_MAX_CHUNK ((size_t)1 << 30)

void caf_put(int image, intptr_t dest, void* src, size_t size)
{
  if (size == 0) return;
  const int rank = image - 1;
  MPI_Aint disp = to_disp(rank, dest, size);
  order_access(rank, disp, disp + (MPI_Aint)size);
  if (rank == myproc) { // local load/stores must observe the result immediately
    memcpy((void*)dest, src, size);
    return;
  }
#if CAF_DEBUG_DEFER_PUTS
  if (size <= CAF_MAX_CHUNK) {
    deferred_put *p = malloc(sizeof(deferred_put) + size);
    if (!p) caf_fatal("out of memory deferring a put");
    p->disp = disp;
    p->size = size;
    memcpy(p->data, src, size);
    p->next = deferred_puts[rank];
    deferred_puts[rank] = p;
    mark_dirty(rank, disp, disp + (MPI_Aint)size);
    return;
  }
#endif
  for (size_t off = 0; off < size; off += CAF_MAX_CHUNK) {
    int n = (int)MIN(CAF_MAX_CHUNK, size - off);
    MPI_SAFE(MPI_Put((byte*)src + off, n, MPI_BYTE, rank, disp + (MPI_Aint)off, n, MPI_BYTE, seg_win));
  }
  // PRIF only requires source completion; remote completion is deferred
  MPI_SAFE(MPI_Win_flush_local(rank, seg_win));
  mark_dirty(rank, disp, disp + (MPI_Aint)size);
}

void caf_get(int image, void* dest, intptr_t src, size_t size)
{
  if (size == 0) return;
  const int rank = image - 1;
  MPI_Aint disp = to_disp(rank, src, size);
  order_access(rank, disp, disp + (MPI_Aint)size);
  if (rank == myproc) {
    memcpy(dest, (void*)src, size);
    return;
  }
  for (size_t off = 0; off < size; off += CAF_MAX_CHUNK) {
    int n = (int)MIN(CAF_MAX_CHUNK, size - off);
    MPI_SAFE(MPI_Get((byte*)dest + off, n, MPI_BYTE, rank, disp + (MPI_Aint)off, n, MPI_BYTE, seg_win));
  }
  // local completion of a get means the data has arrived in dest
  MPI_SAFE(MPI_Win_flush_local(rank, seg_win));
}

// _______________________ Strided RMA ____________________________

// Build a datatype describing a dims-dimensional strided region of
// element_size-byte contiguous blocks. Dimension 0 varies fastest.
// Returns false if the region is empty.
// *lo and *hi receive the lowest and one-past-highest byte offsets touched.
static bool make_strided_type(int dims, const ptrdiff_t *stride, const size_t *extent,
                              size_t element_size, MPI_Datatype *type,
                              ptrdiff_t *lo, ptrdiff_t *hi) {
  if (element_size > INT_MAX) caf_fatal("strided element_size %zu too large", element_size);
  *lo = 0; *hi = (ptrdiff_t)element_size;
  for (int d = 0; d < dims; d++) {
    if (extent[d] == 0) return false;
    if (extent[d] > INT_MAX) caf_fatal("strided extent %zu too large", extent[d]);
    ptrdiff_t span = stride[d] * (ptrdiff_t)(extent[d] - 1);
    if (span < 0) *lo += span; else *hi += span;
  }
  MPI_Datatype t, nt;
  MPI_SAFE(MPI_Type_contiguous((int)element_size, MPI_BYTE, &t));
  for (int d = 0; d < dims; d++) {
    MPI_SAFE(MPI_Type_create_hvector((int)extent[d], 1, (MPI_Aint)stride[d], t, &nt));
    MPI_SAFE(MPI_Type_free(&t));
    t = nt;
  }
  MPI_SAFE(MPI_Type_commit(&t));
  *type = t;
  return true;
}

static void strided_rma(bool is_put, int dims, int image_num,
                        intptr_t remote_ptr, const ptrdiff_t* remote_stride,
                        void *current_image_buffer, const ptrdiff_t *current_image_stride,
                        size_t element_size, const size_t *extent) {
  MPI_Datatype rtype, ltype;
  ptrdiff_t rlo, rhi, llo, lhi;
  if (!make_strided_type(dims, remote_stride, extent, element_size, &rtype, &rlo, &rhi))
    return; // empty transfer
  make_strided_type(dims, current_image_stride, extent, element_size, &ltype, &llo, &lhi);

  const int rank = image_num - 1;
  // validate the full extent of the remote region, then use the region start as displacement
  (void)to_disp(rank, remote_ptr + rlo, (size_t)(rhi - rlo));
  MPI_Aint disp = (MPI_Aint)(remote_ptr - seg_bases[rank]);
  order_access(rank, disp + rlo, disp + rhi);
  if (is_put) {
    MPI_SAFE(MPI_Put(current_image_buffer, 1, ltype, rank, disp, 1, rtype, seg_win));
    if (rank == myproc) { // local load/stores must observe the result immediately
      flush_target(rank);
    } else {
      MPI_SAFE(MPI_Win_flush_local(rank, seg_win));
      mark_dirty(rank, disp + rlo, disp + rhi);
    }
  } else {
    MPI_SAFE(MPI_Get(current_image_buffer, 1, ltype, rank, disp, 1, rtype, seg_win));
    MPI_SAFE(MPI_Win_flush_local(rank, seg_win));
  }
  MPI_SAFE(MPI_Type_free(&rtype));
  MPI_SAFE(MPI_Type_free(&ltype));
}

void caf_put_strided(int dims, int image_num,
                     intptr_t remote_ptr, void* remote_stride,
                     void *current_image_buffer, void * current_image_stride,
                     size_t element_size, void *extent) {
  strided_rma(true, dims, image_num, remote_ptr, remote_stride,
              current_image_buffer, current_image_stride, element_size, extent);
}

void caf_get_strided(int dims, int image_num,
                     intptr_t remote_ptr, void* remote_stride,
                     void *current_image_buffer, void * current_image_stride,
                     size_t element_size, void *extent) {
  strided_rma(false, dims, image_num, remote_ptr, remote_stride,
              current_image_buffer, current_image_stride, element_size, extent);
}

//-------------------------------------------------------------------

// caf_segment_release() is invoked whenever this image is ending a
// segment, to flush any pending actions that are specified to be
// ordered before a subsequent segment.
static void caf_segment_release(void) {
  complete_all_targets();
  MPI_SAFE(MPI_Win_sync(seg_win));
}

void caf_sync_memory() {
  caf_segment_release();
}

void caf_sync_team( caf_team_t team ) {
  caf_segment_release();
  MPI_SAFE(MPI_Barrier(team->comm));
  MPI_SAFE(MPI_Win_sync(seg_win));
}

// _______________________ Events ____________________________

static inline int64_t atomic_fetch_op_local(void *addr, int64_t operand, MPI_Op op) {
  int64_t result = 0;
  MPI_Aint disp = to_disp(myproc, (intptr_t)addr, sizeof(int64_t));
  order_access(myproc, disp, disp + (MPI_Aint)sizeof(int64_t));
  MPI_SAFE(MPI_Fetch_and_op(&operand, &result, MPI_INT64_T, myproc, disp, op, seg_win));
  flush_target(myproc);
  return result;
}

void caf_event_post(int image, intptr_t event_var_ptr, int segment_boundary, int release_fence) {
  assert(event_var_ptr);
  const int rank = image - 1;

  if (segment_boundary) {
    caf_segment_release();
  } else if (release_fence && target_is_dirty(rank)) {
    // MPI does not order a put and an accumulate to different locations, so
    // earlier puts to this target (e.g. the data of a put with NOTIFY=) must
    // be complete before the post can become visible
    flush_target(rank);
  }

  int64_t one = 1;
  MPI_Aint disp = to_disp(rank, event_var_ptr, sizeof(int64_t));
  order_access(rank, disp, disp + (MPI_Aint)sizeof(int64_t));
  MPI_SAFE(MPI_Accumulate(&one, 1, MPI_INT64_T, rank, disp, 1, MPI_INT64_T, MPI_SUM, seg_win));

  if (segment_boundary || release_fence) {
    flush_target(rank);
  } else {
    // Defer completion; this will later be synchronized in caf_segment_release()
    // or before any subsequent wait
    mark_dirty(rank, disp, disp + (MPI_Aint)sizeof(int64_t));
  }
}

void caf_event_query(void *event_var_ptr, int64_t *count) {
  assert(event_var_ptr);
  assert(count);
  *count = atomic_fetch_op_local(event_var_ptr, 0, MPI_NO_OP);
}

void caf_event_wait(void *event_var_ptr, int64_t threshold,
                    int segment_boundary, int acquire_fence, int maybe_concurrent) {
  static pthread_mutex_t notify_wait_lock = PTHREAD_MUTEX_INITIALIZER;
  assert(event_var_ptr);
  assert(threshold >= 1);

  // complete any of our own deferred posts, to ensure the waited-for peers can progress
  if (segment_boundary || num_dirty > 0) caf_segment_release();

  int64_t cnt = 0;
  while (1) {
    while (caf_event_query(event_var_ptr, &cnt), cnt < threshold) {
      // caf_event_query enters the MPI library, which provides progress
    }
    if (maybe_concurrent) pthread_mutex_lock(&notify_wait_lock);
    caf_event_query(event_var_ptr, &cnt);
    if (cnt >= threshold) {
      cnt = atomic_fetch_op_local(event_var_ptr, -threshold, MPI_SUM);
      assert(cnt >= threshold);
      if (maybe_concurrent) pthread_mutex_unlock(&notify_wait_lock);
      break;
    }
    if (maybe_concurrent) pthread_mutex_unlock(&notify_wait_lock);
  }

  if (segment_boundary || acquire_fence) {
    MPI_SAFE(MPI_Win_sync(seg_win));
  }
}

// _______________________ Atomics ____________________________

void caf_atomic_int(int opcode, int image, void* addr, int64_t *result, int64_t op1, int64_t op2) {
  assert(addr);
  const int rank = image - 1;
  MPI_Aint disp = to_disp(rank, (intptr_t)addr, sizeof(int64_t));
  order_access(rank, disp, disp + (MPI_Aint)sizeof(int64_t));

  switch (opcode) {
    case CAF_OP_GET:
      MPI_SAFE(MPI_Fetch_and_op(NULL, result, MPI_INT64_T, rank, disp, MPI_NO_OP, seg_win));
      break;
    case CAF_OP_SET:
      MPI_SAFE(MPI_Accumulate(&op1, 1, MPI_INT64_T, rank, disp, 1, MPI_INT64_T, MPI_REPLACE, seg_win));
      break;
    case CAF_OP_ADD:
      MPI_SAFE(MPI_Accumulate(&op1, 1, MPI_INT64_T, rank, disp, 1, MPI_INT64_T, MPI_SUM, seg_win));
      break;
    case CAF_OP_AND:
      MPI_SAFE(MPI_Accumulate(&op1, 1, MPI_INT64_T, rank, disp, 1, MPI_INT64_T, MPI_BAND, seg_win));
      break;
    case CAF_OP_OR:
      MPI_SAFE(MPI_Accumulate(&op1, 1, MPI_INT64_T, rank, disp, 1, MPI_INT64_T, MPI_BOR, seg_win));
      break;
    case CAF_OP_XOR:
      MPI_SAFE(MPI_Accumulate(&op1, 1, MPI_INT64_T, rank, disp, 1, MPI_INT64_T, MPI_BXOR, seg_win));
      break;
    case CAF_OP_FADD:
      MPI_SAFE(MPI_Fetch_and_op(&op1, result, MPI_INT64_T, rank, disp, MPI_SUM, seg_win));
      break;
    case CAF_OP_FAND:
      MPI_SAFE(MPI_Fetch_and_op(&op1, result, MPI_INT64_T, rank, disp, MPI_BAND, seg_win));
      break;
    case CAF_OP_FOR:
      MPI_SAFE(MPI_Fetch_and_op(&op1, result, MPI_INT64_T, rank, disp, MPI_BOR, seg_win));
      break;
    case CAF_OP_FXOR:
      MPI_SAFE(MPI_Fetch_and_op(&op1, result, MPI_INT64_T, rank, disp, MPI_BXOR, seg_win));
      break;
    case CAF_OP_FCAS: // op1 = compare, op2 = new value
      MPI_SAFE(MPI_Compare_and_swap(&op2, &op1, result, MPI_INT64_T, rank, disp, seg_win));
      break;
    default:
      caf_fatal("caf_atomic_int: invalid opcode %d", opcode);
  }
  // PRIF atomics are fully blocking
  flush_target(rank);
}

void caf_atomic_logical(int opcode, int image, void* addr, int64_t *result, int64_t op1, int64_t op2) {
  caf_atomic_int(opcode, image, addr, result, op1, op2);
}

//-------------------------------------------------------------------
// gfortran 13.2 .. 15.2 : c_funloc is non-compliant:
//   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=124652
// it erroneously generates a non-callable pointer to a pointer to the subroutine
// This helper is used to undo that incorrect extra level of indirection
typedef void (*funloc_t)(void);
funloc_t caf_c_funloc_deref(funloc_t funloc) {
  return *(funloc_t *)funloc;
}

// Type-erased collective subroutines
//-------------------------------------------------------------------

// Context for the reduction currently in progress.
// PRIF collectives are blocking, so a single active context suffices.
static struct {
  caf_reduce_fn_t fn;
  void *cdata;
} reduce_ctx;

static void reduce_trampoline(void *invec, void *inoutvec, int *len, MPI_Datatype *dt) {
  (void)dt;
  assert(reduce_ctx.fn);
  reduce_ctx.fn(invec, inoutvec, (size_t)*len, reduce_ctx.cdata);
}

// Reduce num_elements values of type dt (each elem_bytes bytes) in place at a_ptr
static void reduce_in_place(void *a_ptr, int result_image, size_t num_elements, size_t elem_bytes,
                            MPI_Datatype dt, MPI_Op op, caf_team_t team) {
  const size_t chunk = (size_t)INT_MAX;
  const int root = result_image - 1;
  for (size_t off = 0; off < num_elements; off += chunk) {
    int n = (int)MIN(chunk, num_elements - off);
    byte *p = (byte*)a_ptr + off * elem_bytes;
    if (result_image) {
      if (team->rank == root) {
        MPI_SAFE(MPI_Reduce(MPI_IN_PLACE, p, n, dt, op, root, team->comm));
      } else {
        MPI_SAFE(MPI_Reduce(p, NULL, n, dt, op, root, team->comm));
      }
    } else {
      MPI_SAFE(MPI_Allreduce(MPI_IN_PLACE, p, n, dt, op, team->comm));
    }
  }
}

void caf_co_reduce_cptr( void *a_ptr, int result_image, size_t num_elements, size_t element_size,
                         caf_reduce_fn_t op_wrapper, void* client_data, caf_team_t team) {
  assert(result_image >= 0);
  assert(num_elements > 0);
  assert(op_wrapper);
  if (element_size > INT_MAX) caf_fatal("co_reduce: element size %zu too large", element_size);

  MPI_Datatype dt;
  MPI_SAFE(MPI_Type_contiguous((int)element_size, MPI_BYTE, &dt));
  MPI_SAFE(MPI_Type_commit(&dt));

  reduce_ctx.fn = op_wrapper;
  reduce_ctx.cdata = client_data;
  reduce_in_place(a_ptr, result_image, num_elements, element_size, dt, user_reduce_op, team);
  reduce_ctx.fn = NULL;
  reduce_ctx.cdata = NULL;

  MPI_SAFE(MPI_Type_free(&dt));
}

void caf_co_reduce( CFI_cdesc_t* a_desc, int result_image, size_t num_elements,
                    caf_reduce_fn_t op_wrapper, void* client_data, caf_team_t team) {
  assert(a_desc);
  char* a_ptr = (char*) a_desc->base_addr;
  size_t element_size = a_desc->elem_len;
  caf_co_reduce_cptr(a_ptr, result_image, num_elements, element_size,
                     op_wrapper, client_data, team);
}

void caf_co_broadcast_cptr(void *a_ptr, int source_image, size_t nbytes, caf_team_t team) {
  assert(a_ptr);
  assert(source_image >= 1);
  assert(nbytes > 0);

  for (size_t off = 0; off < nbytes; off += CAF_MAX_CHUNK) {
    int n = (int)MIN(CAF_MAX_CHUNK, nbytes - off);
    MPI_SAFE(MPI_Bcast((byte*)a_ptr + off, n, MPI_BYTE, source_image - 1, team->comm));
  }
}

void caf_co_broadcast(CFI_cdesc_t * a_desc, int source_image, int num_elements, caf_team_t team) {
  assert(a_desc);
  char* a_ptr = (char*) a_desc->base_addr;
  size_t element_size = a_desc->elem_len;
  size_t nbytes = (size_t)num_elements * element_size;
  caf_co_broadcast_cptr(a_ptr, source_image, nbytes, team);
}

//-------------------------------------------------------------------
// Typed computational collective subroutines
//-------------------------------------------------------------------

static MPI_Datatype int_type_for_size(size_t sz, CFI_type_t cfi_type) {
  switch (sz) {
    case 1: return MPI_INT8_T;
    case 2: return MPI_INT16_T;
    case 4: return MPI_INT32_T;
    case 8: return MPI_INT64_T;
  }
  caf_fatal("Unsupported integer type: %d (size %zu)", (int)cfi_type, sz);
}

// Convert CFI_type_t to the corresponding MPI datatype.
// Sets *is_complex for complex types.
static MPI_Datatype CFI_to_MPI_DT(CFI_type_t cfi_type, size_t elem_len, bool *is_complex) {
  *is_complex = false;

  // real and complex cases
  // (an if-chain because long double may be a duplicate of double on some platforms)
  if (cfi_type == CFI_type_float)  return MPI_FLOAT;
  if (cfi_type == CFI_type_double) return MPI_DOUBLE;
  if (cfi_type == CFI_type_long_double) return MPI_LONG_DOUBLE;
  *is_complex = true;
  if (cfi_type == CFI_type_float_Complex)       return MPI_C_FLOAT_COMPLEX;
  if (cfi_type == CFI_type_double_Complex)      return MPI_C_DOUBLE_COMPLEX;
  if (cfi_type == CFI_type_long_double_Complex) return MPI_C_LONG_DOUBLE_COMPLEX;
  *is_complex = false;

  // integer types
  // these must be handled outside the switch because there are duplicates
  // for the same reason, start with the most likely candidates
  if (cfi_type == CFI_type_int64_t ||
      cfi_type == CFI_type_int32_t ||
      cfi_type == CFI_type_int16_t ||
      cfi_type == CFI_type_int8_t ||
      cfi_type == CFI_type_Bool ||
      cfi_type == CFI_type_char ||
      cfi_type == CFI_type_signed_char ||
      cfi_type == CFI_type_short ||
      cfi_type == CFI_type_int ||
      cfi_type == CFI_type_long ||
      cfi_type == CFI_type_long_long ||
      cfi_type == CFI_type_size_t ||
      cfi_type == CFI_type_int_least8_t ||
      cfi_type == CFI_type_int_least16_t ||
      cfi_type == CFI_type_int_least32_t ||
      cfi_type == CFI_type_int_least64_t ||
      cfi_type == CFI_type_int_fast8_t ||
      cfi_type == CFI_type_int_fast16_t ||
      cfi_type == CFI_type_int_fast32_t ||
      cfi_type == CFI_type_int_fast64_t ||
      cfi_type == CFI_type_intmax_t ||
      cfi_type == CFI_type_intptr_t ||
      cfi_type == CFI_type_ptrdiff_t) {
    return int_type_for_size(elem_len, cfi_type);
  }

  caf_fatal("Unrecognized type: %d", (int)cfi_type);
}

static void caf_co_common(CFI_cdesc_t* a_desc, int result_image, size_t num_elements, caf_team_t team, MPI_Op op) {
  bool is_complex;
  MPI_Datatype dt = CFI_to_MPI_DT(a_desc->type, a_desc->elem_len, &is_complex);
  if (is_complex && op != MPI_SUM)
    caf_fatal("This operation does not support complex types");
  reduce_in_place(a_desc->base_addr, result_image, num_elements, a_desc->elem_len, dt, op, team);
}

void caf_co_max(CFI_cdesc_t* a_desc, int result_image, size_t num_elements, caf_team_t team) {
  caf_co_common(a_desc, result_image, num_elements, team, MPI_MAX);
}

void caf_co_min(CFI_cdesc_t* a_desc, int result_image, size_t num_elements, caf_team_t team) {
  caf_co_common(a_desc, result_image, num_elements, team, MPI_MIN);
}

void caf_co_sum(CFI_cdesc_t* a_desc, int result_image, size_t num_elements, caf_team_t team) {
  caf_co_common(a_desc, result_image, num_elements, team, MPI_SUM);
}

#endif // CAF_RUNTIME_MPI
