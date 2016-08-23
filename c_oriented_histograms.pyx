cimport numpy as np
from libc.math cimport floor
def calculate_histograms (np.ndarray[np.double_t, ndim=2] mag,
                          np.ndarray[np.double_t, ndim=2] ang,
                          cx, cy, sx, sy, n_cellsx, n_cellsy, n_orientations, ovx, ovy,
                          np.ndarray[np.double_t, ndim=3] orientation_histogram):
    # Defining variables
    cdef double cell_ang, cell_mag
    cdef double dist_cur_center, dist_prv_center, dist_nxt_center
    cdef double bin_width = 180. / n_orientations
    cdef int cur_bin, prv_bin, nxt_bin
    cdef Py_ssize_t idx_X, idx_Y, x, y, from_x, from_y
    #cdef int idx_X, idx_Y, x, y, to_x, to_y
    # Calculating the histograms for each cell
    cdef int gap_cx = cx - ovx;
    cdef int gap_cy = cy - ovy;
    for idx_X from 0 <= idx_X < n_cellsx:
        for idx_Y from 0 <= idx_Y < n_cellsy:
            from_x = idx_X * gap_cx 
            from_y = idx_Y * gap_cy
            for x from from_x <= x < (from_x + cx):
                for y from from_y <= y < (from_y + cy):
                    cell_ang = ang[y, x]
                    cell_mag = mag[y, x]
                    cur_bin = int(floor(cell_ang / bin_width))
                    # Check when ang is 180
                    if cur_bin == n_orientations:
                        cur_bin -= 1
                    #assert(0 <= cur_bin and cur_bin <= n_orientations)
                    dist_cur_center = cell_ang - (cur_bin * bin_width + bin_width * 0.5)
                    if 0 < dist_cur_center:
                        # share with next bin
                        nxt_bin = cur_bin + 1
                        dist_nxt_center = (nxt_bin * bin_width + bin_width * 0.5) - cell_ang
                        if nxt_bin == n_orientations:
                            nxt_bin = 0
                        #continue
                        orientation_histogram[idx_Y, idx_X, cur_bin] += (dist_nxt_center * 1.0 / bin_width) * cell_mag
                        orientation_histogram[idx_Y, idx_X, nxt_bin] += (dist_cur_center * 1.0 / bin_width) * cell_mag
                    else:
                        # share with prev bin
                        prv_bin = cur_bin - 1
                        dist_prv_center = cell_ang - (prv_bin * bin_width + bin_width * 0.5)
                        #continue
                        if prv_bin == -1:
                            prv_bin = n_orientations - 1
                        assert(0 <= cur_bin and cur_bin < n_orientations)
                        assert(0 <= prv_bin and prv_bin < n_orientations)
                        orientation_histogram[idx_Y, idx_X, cur_bin] += (dist_prv_center * 1.0 / bin_width) * cell_mag
                        orientation_histogram[idx_Y, idx_X, prv_bin] += (-dist_cur_center * 1.0 / bin_width) * cell_mag
