cimport numpy as np
from libc.math cimport floor
def calculate_histograms (np.ndarray[np.double_t, ndim=2] mag,
                          np.ndarray[np.double_t, ndim=2] ang,
                          cx, cy, sx, sy, n_cellsx, n_cellsy, n_orientations,
                          np.ndarray[np.double_t, ndim=3] orientation_histogram):
    # Defining variables
    cdef double cell_ang, cell_mag
    cdef double dist_cur_center, dist_prv_center, dist_nxt_center
    cdef double bin_width = 180. / n_orientations
    cdef int cur_bin, prv_bin, nxt_bin
    cdef Py_ssize_t idx_X, idx_Y, x, y, to_x, to_y
    #cdef int idx_X, idx_Y, x, y, to_x, to_y
    # Calculating the histograms for each cell
    for idx_X from 0 <= idx_X < n_cellsx:
        for idx_Y from 0 <= idx_Y < n_cellsy:
            to_x = (idx_X + 1) * cx
            to_y = (idx_Y + 1) * cy
            for x from (idx_X * cx) <= x < to_x:
                for y from (idx_Y * cy) <= y < to_y:
                    cell_ang = ang[y, x]
                    cell_mag = mag[y, x]
                    cur_bin = int(floor(cell_ang / bin_width))
                    # Check when ang is 180
                    if cur_bin == n_orientations:
                        cur_bin -= 1
                    dist_cur_center = cell_ang - (cur_bin * bin_width + bin_width * 0.5)
                    if 0 < dist_cur_center:
                        # share with next bin
                        nxt_bin = cur_bin + 1
                        dist_nxt_center = (nxt_bin * bin_width + bin_width * 0.5) - cell_ang
                        if nxt_bin == n_orientations:
                            nxt_bin = 0
                        orientation_histogram[idx_Y, idx_X, cur_bin] += (dist_nxt_center * 1.0 / bin_width) * cell_mag
                        orientation_histogram[idx_Y, idx_X, nxt_bin] += (dist_cur_center * 1.0 / bin_width) * cell_mag
                    else:
                        # share with prev bin
                        prv_bin = cur_bin - 1
                        dist_prv_center = cell_ang - (prv_bin * bin_width + bin_width * 0.5)
                        if prv_bin == -1:
                            prv_bin = n_orientations - 1
                        orientation_histogram[idx_Y, idx_X, cur_bin] += (dist_prv_center * 1.0 / bin_width) * cell_mag
                        orientation_histogram[idx_Y, idx_X, prv_bin] += (-dist_cur_center * 1.0 / bin_width) * cell_mag
