cimport numpy as np

def hik_kernel(np.ndarray[np.double_t, ndim=2] X,
               np.ndarray[np.double_t, ndim=2] Y,
               n_rowsX, n_rowsY, n_features,
               np.ndarray[np.double_t, ndim=2] K):
    """
        Hik kernel:
            K(x, y) = Sum (min(x, y))
    """
    # Defining variables
    cdef double ans
    cdef Py_ssize_t idx_rowX, idx_rowY, idx_feat
    # Calculating the HIK Kernel
    for idx_rowX from 0 <= idx_rowX < n_rowsX:
        for idx_rowY from 0 <= idx_rowY < n_rowsY:
            ans = 0
            for idx_feat from  0 <= idx_feat < n_features:
                ans += min(X[idx_rowX, idx_feat], Y[idx_rowY, idx_feat])
            K[idx_rowX, idx_rowY] = ans

def chi_square_kernel_base(np.ndarray[np.double_t, ndim=2] X,
                           np.ndarray[np.double_t, ndim=2] Y,
                           n_rowsX, n_rowsY, n_features,
                           np.ndarray[np.double_t, ndim=2] B):
    """
        Chi2 kernel base:
            B(x, y) = Sum ((x - y)^2 / (x + y))
        To get the chi2 kernel:
            K = exp(-gamma * B)
    """
    # Defining variables
    cdef double eps = 1E-8
    cdef double ans, tmp
    cdef Py_ssize_t idx_rowX, idx_rowY, idx_feat
    # Calculating the HIK Kernel
    for idx_rowX from 0 <= idx_rowX < n_rowsX:
        for idx_rowY from 0 <= idx_rowY < n_rowsY:
            ans = 0
            for idx_feat from  0 <= idx_feat < n_features:
                tmp = X[idx_rowX, idx_feat] - Y[idx_rowY, idx_feat]
                ans += (tmp * tmp) / (eps + X[idx_rowX, idx_feat] + Y[idx_rowY, idx_feat])
            B[idx_rowX, idx_rowY] = ans
