# Purpose: Cython wrapper for the CCF module
# Author: Attila Bódi
# Version: 2022MARC06

import cython
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern double ccfpix_c(
        double* M_L,
        double* M_H,
        double* X,
        double* THAR,
        const double DELTA,
        const int N,
        const int M)

cdef extern void ccfcos_c(
        double* M_L,
        double* M_H,
        double* WAV,
        double* SPEC,
        double* WEIGHT,
        double* SN,
        const double V_R,
        const int N,
        const int M,
        double* CCFCOS,
        double* SNW)

@cython.boundscheck(False)
@cython.wraparound(False)
def ccfpix( np.ndarray[np.double_t, ndim=1, mode="c"] m_l not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] m_h not None,
            np.ndarray xin not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] thar not None,
            double delta
    ):
    """
    This function calculates the cross-correlation function.

    Parameters
    ----------
    m_l: array
    m_h: array
    x: array
    thar: array
    delta: float

    Returns
    -------
    ccfpix : float
    """

    cdef np.ndarray[np.double_t, ndim=1, mode="c"] x = xin.astype(np.float64)

    del xin

    cdef int n = m_l.size
    cdef int m = x.size
    cdef double ccfpix

    # --- Call C funtion ---
    ccfpix = ccfpix_c(
              &m_l[0],
              &m_h[0],
              &x[0],
              &thar[0],
              delta,
              n,
              m)

    return ccfpix


@cython.boundscheck(False)
@cython.wraparound(False)
def ccfcos( np.ndarray[np.double_t, ndim=1, mode="c"] m_l not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] m_h not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] WAV not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] SPEC not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] WEIGHT not None,
            np.ndarray[np.double_t, ndim=1, mode="c"] SN not None,
            double V_R
    ):
    """
    This function calculates the cross-correlation function
    of a spectrum (SPE, at wavelengths WAV) with a mask.

    Parameters
    ----------
    m_l: array
    m_h: array
    WAV: array
    SPEC: array
    WEIGHT: array
    SN: array
    V_R: float

    Returns
    -------
    ccfcos : float
    snw : float
    """

    cdef int n = m_l.size
    cdef int m = WAV.size
    cdef double ccfcos
    cdef double snw

    # --- Call C funtion ---
    ccfcos_c(
            &m_l[0],
            &m_h[0],
            &WAV[0],
            &SPEC[0],
            &WEIGHT[0],
            &SN[0],
            V_R,
            n,
            m,
            &ccfcos,
            &snw)

    return ccfcos,snw
