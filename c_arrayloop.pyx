cimport numpy as np
def array2d_3d_loop (np.ndarray[np.double_t, ndim=2] A,
                     np.ndarray[np.double_t, ndim=3] B,
                     NA,
                     NB):
    # Defining variables
    cdef int c = 0;
    cdef int i, j, k;
    for i from 0 <= i < NA[0]:
        for j from 0 <= j < NA[1]:
            A[i, j] = c
            c += 1
    for i from 0 <= i < NB[0]:
        for j from 0 <= j < NB[1]:
            for k from 0 <= k < NB[2]:
                B[i, j, k] = c
                c += 1
