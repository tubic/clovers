# cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3
# distutils: extra_compile_args = -O3 -mavx -mfma -fopenmp
# distutils: extra_link_args = -fopenmp
"""Cython implementation of the Z-curve encoder.

Encodes a list of DNA sequences into 189-dimensional Z-curve feature
vectors (mono 9D + di 36D + tri 144D), with OpenMP parallelism.

Public API (identical to the previous C extension):
    encode(records: List[str], n_jobs: int = 0) -> np.ndarray (n, 189) float32
"""
import numpy as np
cimport numpy as cnp
from cython.parallel import prange, parallel
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

cnp.import_array()

cdef extern from *:
    int omp_get_max_threads() nogil
    void omp_set_num_threads(int) nogil

DEF DIM_A = 189
DEF PHASE = 3
DEF N_BASE = 128   # ASCII table size
DEF ONE_HOT_SIZE = 512   # N_BASE * 4
DEF Z_COORD_SIZE = 384   # N_BASE * 3

cdef float V = 1.0 / 2
cdef float W = 1.0 / 3
cdef float Q = 1.0 / 4

# ---------------------------------------------------------------------------
# Lookup tables (identical to the original Zcurve.cpp)
# C arrays are declared uninitialized and filled from the Python lists at
# import time (Cython does not allow init lists on module-level cdef arrays)
# ---------------------------------------------------------------------------
cdef float ONE_HOT[ONE_HOT_SIZE]
cdef float Z_COORD[Z_COORD_SIZE]

_ONE_HOT_PY = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 1.0/3, 1.0/3, 1.0/3, 0, 0, 1, 0, 1.0/3, 1.0/3, 0, 1.0/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    1.0/3, 0, 1.0/3, 1.0/3, 1.0/4, 1.0/4, 1.0/4, 1.0/4, 0, 0, 0, 0, 0, 1.0/2, 0, 1.0/2, 0, 0, 0, 0, 1.0/2, 0, 1.0/2, 0, 1.0/4, 1.0/4, 1.0/4, 1.0/4, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1.0/2, 1.0/2, 0, 0, 0, 1.0/2, 1.0/2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1.0/3, 1.0/3, 1.0/3, 0, 1.0/2, 0, 0, 1.0/2,
    0, 0, 0, 0, 0, 0, 1.0/2, 1.0/2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 1.0/3, 1.0/3, 1.0/3, 0, 0, 1, 0, 1.0/3, 1.0/3, 0, 1.0/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    1.0/3, 0, 1.0/3, 1.0/3, 1.0/4, 1.0/4, 1.0/4, 1.0/4, 0, 0, 0, 0, 0, 1.0/2, 0, 1.0/2, 0, 0, 0, 0, 1.0/2, 0, 1.0/2, 0, 1.0/4, 1.0/4, 1.0/4, 1.0/4, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1.0/2, 1.0/2, 0, 0, 0, 1.0/2, 1.0/2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1.0/3, 1.0/3, 1.0/3, 0, 1.0/2, 0, 0, 1.0/2,
    0, 0, 0, 0, 0, 0, 1.0/2, 1.0/2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]

# Z_COORD: ASCII char -> Z-curve coordinates (degenerate bases by weighted sum)
_Z_COORD_PY = [
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +1, +1, +1, -1.0/3, -1.0/3, -1.0/3, -1, +1, -1, +1.0/3, -1.0/3, +1.0/3, +0, +0, +0, +0, +0, +0, +1, -1, -1,
    -1.0/3, +1.0/3, +1.0/3, +0, +0, +0, +0, +0, +0, +0, -1, +0, +0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +1, +0, +0, +0, +0, -1, -1, -1, +1, -1, -1, +1, +1.0/3, +1.0/3, -1.0/3, +0, +0, +1,
    +0, +0, +0, +1, +0, +0, +1, +1, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +1, +1, +1, -1.0/3, -1.0/3, -1.0/3, -1, +1, -1, +1.0/3, -1.0/3, +1.0/3, +0, +0, +0, +0, +0, +0, +1, -1, -1,
    -1.0/3, +1.0/3, +1.0/3, +0, +0, +0, +0, +0, +0, +0, -1, +0, +0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +1, +0, +0, +0, +0, -1, -1, -1, +1, -1, -1, +1, +1.0/3, +1.0/3, -1.0/3, +0, +0, +1,
    +0, +0, +0, +1, +0, +0, +1, +1, +1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]

cdef void _init_tables() noexcept:
    cdef int i
    for i in range(ONE_HOT_SIZE):
        ONE_HOT[i] = _ONE_HOT_PY[i]
    for i in range(Z_COORD_SIZE):
        Z_COORD[i] = _Z_COORD_PY[i]

_init_tables()

# ---------------------------------------------------------------------------
# Core transforms
# ---------------------------------------------------------------------------
cdef inline void mono_trans(const char* seq, int length, float* params) noexcept nogil:
    cdef float counts[PHASE][3]
    cdef int i, p, b
    for p in range(PHASE):
        for b in range(3):
            counts[p][b] = 0.0
    for i in range(length):
        p = i % PHASE
        counts[p][0] += Z_COORD[seq[i] * 3 + 0]
        counts[p][1] += Z_COORD[seq[i] * 3 + 1]
        counts[p][2] += Z_COORD[seq[i] * 3 + 2]
    for p in range(PHASE):
        params[p * 3 + 0] = counts[p][0] / length * PHASE
        params[p * 3 + 1] = counts[p][1] / length * PHASE
        params[p * 3 + 2] = counts[p][2] / length * PHASE

cdef inline void di_trans(const char* seq, int length, float* params) noexcept nogil:
    cdef float counts[PHASE][4][3]
    cdef int i, p, b
    for p in range(PHASE):
        for b in range(4):
            counts[p][b][0] = 0.0
            counts[p][b][1] = 0.0
            counts[p][b][2] = 0.0
    for i in range(length - 1):
        p = i % PHASE
        for b in range(4):
            counts[p][b][0] += ONE_HOT[seq[i] * 4 + b] * Z_COORD[seq[i + 1] * 3 + 0]
            counts[p][b][1] += ONE_HOT[seq[i] * 4 + b] * Z_COORD[seq[i + 1] * 3 + 1]
            counts[p][b][2] += ONE_HOT[seq[i] * 4 + b] * Z_COORD[seq[i + 1] * 3 + 2]
    for p in range(PHASE):
        for b in range(4):
            params[(p * 4 + b) * 3 + 0] = counts[p][b][0] / length * PHASE
            params[(p * 4 + b) * 3 + 1] = counts[p][b][1] / length * PHASE
            params[(p * 4 + b) * 3 + 2] = counts[p][b][2] / length * PHASE

cdef inline void tri_trans(const char* seq, int length, float* params) noexcept nogil:
    cdef float counts[PHASE][4][4][3]
    cdef int i, p, s, b
    for p in range(PHASE):
        for s in range(4):
            for b in range(4):
                counts[p][s][b][0] = 0.0
                counts[p][s][b][1] = 0.0
                counts[p][s][b][2] = 0.0
    for i in range(length - 2):
        p = i % PHASE
        for s in range(4):
            for b in range(4):
                counts[p][s][b][0] += ONE_HOT[seq[i] * 4 + s] * ONE_HOT[seq[i + 1] * 4 + b] * Z_COORD[seq[i + 2] * 3 + 0]
                counts[p][s][b][1] += ONE_HOT[seq[i] * 4 + s] * ONE_HOT[seq[i + 1] * 4 + b] * Z_COORD[seq[i + 2] * 3 + 1]
                counts[p][s][b][2] += ONE_HOT[seq[i] * 4 + s] * ONE_HOT[seq[i + 1] * 4 + b] * Z_COORD[seq[i + 2] * 3 + 2]
    for p in range(PHASE):
        for s in range(4):
            for b in range(4):
                params[((p * 4 + s) * 4 + b) * 3 + 0] = counts[p][s][b][0] / length * PHASE
                params[((p * 4 + s) * 4 + b) * 3 + 1] = counts[p][s][b][1] / length * PHASE
                params[((p * 4 + s) * 4 + b) * 3 + 2] = counts[p][s][b][2] / length * PHASE

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def encode(records, int n_jobs=0):
    """Encode DNA sequences into 189D Z-curve feature vectors.

    Parameters
    ----------
    records : List[str]
        DNA sequences (ASCII).
    n_jobs : int, optional
        Number of OpenMP threads; <= 0 means all available cores.

    Returns
    -------
    np.ndarray of shape (n, 189), dtype float32
    """
    cdef int n, i
    cdef bytes b
    cdef list blobs
    cdef char** seqs
    cdef int* lens
    cdef cnp.ndarray[cnp.float32_t, ndim=2, mode='c'] out
    cdef float* params

    if isinstance(records, (str, bytes)):
        raise TypeError("records must be List[str]")
    try:
        n = len(records)
    except TypeError:
        raise TypeError("records must be List[str]")
    if n == 0:
        return np.empty((0, DIM_A), dtype=np.float32)

    blobs = []
    for r in records:
        if not isinstance(r, str):
            raise TypeError("records must be List[str]")
        blobs.append(r.encode('ascii'))

    seqs = <char**>malloc(n * sizeof(char*))
    lens = <int*>malloc(n * sizeof(int))
    if seqs == NULL or lens == NULL:
        free(seqs); free(lens)
        raise MemoryError("Failed to allocate memory")

    for i in range(n):
        seqs[i] = blobs[i]
        lens[i] = len(blobs[i])

    out = np.empty((n, DIM_A), dtype=np.float32)
    params = &out[0, 0]

    if n_jobs <= 0:
        n_jobs = omp_get_max_threads()
    omp_set_num_threads(n_jobs)

    with nogil, parallel():
        for i in prange(n):
            mono_trans(seqs[i], lens[i], params + i * DIM_A)
            di_trans(seqs[i], lens[i], params + i * DIM_A + 9)
            tri_trans(seqs[i], lens[i], params + i * DIM_A + 45)

    free(seqs)
    free(lens)
    return out
