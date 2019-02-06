# cython: boundscheck=False
# cython: language_level=2

from libc.math cimport sqrt
cimport numpy as np
import numpy as np


def LSE(double[:] x, double[:] y, int n, float bad):
    """
    This is a Least Square Estimate subroutine to fit a linear
    equation for (xi, yi) (i = 1, ..., n), so that yi = a * xi + b
    INPUTs: x(i), y(i), n, (i = 1, ..., n ).
    OUTPUTs: a (slope), b (intercept).
    Li Liu - Sep. 23, 92
    """
    cdef:
        float a
        float b
        float xsum = 0.0
        float ysum = 0.0
        float xxsum = 0.0
        float xysum = 0.0
        float det
        float total
        Py_ssize_t i

    total = float(n)
    for i in range(n):
        if x[i] == bad or y[i] == bad:
            total -= 1.0
        else:
            xsum += x[i]
            ysum += y[i]
            xxsum += x[i] * x[i]
            xysum += x[i] * y[i]

    det = total * xxsum - xsum**2
    a = (total * xysum - xsum * ysum) / det
    b = (ysum * xxsum - xsum * xysum) / det
    return a, b


def calc_kdp_ray_fir(
        int ngates, float[:] dp, float[:] dz, float[:] rng,
        float[:] thsd, int nf, float bad, int fir_order, float fir_gain,
        double[:] fir_coeff, int std_gate):
    """
    Calculate KDP along a 1-D ray of radar data, using an FIR filter to first
    smooth the differential phase. Standard deviation of phase is also
    calculated.

    Arguments
    ---------
    dp = 1D ray of differential phase
    dz = 1D ray of reflectivity
    rng = 1D ray of range
    thsd = Scalar or 1D ray of diff phase stddev thresholds
    nf = Number of times to filter the data
    bad = Bad/missing data value
    fir_order = FIR filter order
    fir_gain = FIR filter gain
    fir_coeff = FIR filter coefficients
    std_gate = Number of gates to use for diff phase stddev calc

    Returns
    -------
    kd_lin = Specific differential phase (deg/km, 1D array)
    dp_lin = Filtered differential phase (deg, 1D array)
    sd_lin = Standard deviation of diff. phase (deg, 1D array)
    """
    # Define needed variables
    cdef:
        np.ndarray[np.float64_t, ndim=1] xx = bad + np.zeros(500, dtype=np.float64)
        np.ndarray[np.float64_t, ndim=1] yy = bad + np.zeros(500, dtype=np.float64)
        np.ndarray[np.float32_t, ndim=1] y = np.zeros(ngates, dtype=np.float32)
        np.ndarray[np.float32_t, ndim=1] z = np.zeros(ngates, dtype=np.float32)
        np.ndarray[np.float32_t, ndim=1] kd_lin = np.zeros(ngates, dtype=np.float32)
        np.ndarray[np.float32_t, ndim=1] dp_lin = np.zeros(ngates, dtype=np.float32)
        np.ndarray[np.float32_t, ndim=1] sd_lin = np.zeros(ngates, dtype=np.float32)
        int half_std_win, nadp, half_fir_win
        int index1, index2, N, half_nadp
        double X, A, V, W, aa, bb
        Py_ssize_t i, j, mloop

    # Half window size for calculating stdev phase
    half_std_win = (std_gate - 1) / 2
    # Half window size for FIR filtering
    half_fir_win = fir_order / 2

    #####################################################################
    # Calculate standard deviation of phidp
    for i in range(ngates):
        kd_lin[i] = bad
        sd_lin[i] = 100.0
        y[i] = bad
        z[i] = dp[i]
        index1 = i - half_std_win
        index2 = i + half_std_win
        if index1 >= 0 and index2 <= ngates - 1:
            N = 0
            A = 0.0
            V = 0.0
            for j in range(index1, index2 + 1, 1):
                if dp[j] != bad:
                    # Standard deviation algorithm
                    X = dp[j]
                    if N <= 0:
                        W = X
                    N += 1
                    D = X - W
                    V = (N - 1) * (D - A)**2 / N + V
                    A = (D - A) / N + A
            if N > half_std_win:
                sd_lin[i] = sqrt(V / N)

    # ------------- MAIN LOOP of Phidp Adaptive Filtering ------------------
    # FIR FILTER SECTION
    for mloop in range(nf):
        for i in range(half_fir_win, ngates - half_fir_win, 1):
            if sd_lin[i] <= thsd[i] and z[i] != bad:
                index1 = i - half_fir_win
                index2 = i + half_fir_win
                N = 0
                for j in range(index1, index2 + 1, 1):
                    if sd_lin[j] <= thsd[j] and z[j] != bad:
                        yy[N] = z[j]
                        xx[N] = rng[j]
                        N += 1

                # Now fill in gaps if they aren't too big
                if float(N) > (0.8 * float(fir_order)):
                    if N < fir_order + 1:
                        aa, bb = LSE(xx, yy, N, bad)
                        for j in range(index1, index2 + 1, 1):
                            if z[j] == bad:
                                z[j] = aa * rng[j] + bb

                    # Now do the FIR filtering
                    A = 0.0
                    for j in range(index1, index2 + 1, 1):
                        A += fir_coeff[j - index1 + 1] * z[j]
                    y[i] = A * fir_gain
        for i in range(ngates):
            z[i] = y[i]  # Enables re-filtering of processed phase
    for i in range(ngates):
        dp_lin[i] = z[i]

    # *****************END LOOP for Phidp Adaptive Filtering******************

    # CALCULATE KDP
    # Default value for nadp is half_fir_win, but varies based on Zh
    for i in range(ngates):
        if dz[i] != bad:
            if dz[i] >= 45.0:
                nadp = half_fir_win
            if dz[i] >= 35.0 and dz[i] < 45.0:
                nadp = 2 * half_fir_win
            if dz[i] < 35.0:
                nadp = 3 * half_fir_win
            half_nadp = nadp / 2
            index1 = i - half_nadp
            index2 = i + half_nadp
            N = 0
            for j in range(index1, index2, 1):
                if index1 >= 0 and index2 <= ngates - 1:
                    if dp_lin[j] != bad:
                        yy[N] = dp_lin[j]
                        xx[N] = rng[j]
                        N += 1
            if float(N) > (0.8 * float(nadp)):
                aa, bb = LSE(xx, yy, N, bad)
                kd_lin[i] = 0.5 * aa

    # *******************END KDP CALCULATION****************************
    return kd_lin, dp_lin, sd_lin


def hid_beta_f(int ngates, float[:] x_arr, float a, float b, float m):
    """
    Beta function calculator
    ngates = Number of gates
    x_arr = Array to process
    a, b , m = Beta function parameters
    """

    cdef:
        np.ndarray[np.float32_t, ndim=1] beta = np.zeros(ngates, dtype=np.float32)
        Py_ssize_t i

    for i in range(ngates):
        beta[i] = 1.0 / (1.0 + (((x_arr[i] - m) / a)**2.0)**b)

    return beta
