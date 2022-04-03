import math
import numpy as np


def cwt(data, wavelet_name, sampling_frequency=1.):

    # Currently only supported for Morlet wavelets
    if wavelet_name == 'morl':
        n_orig = data.size
        nv = 10
        ds = 1 / nv
        fs = sampling_frequency
        dt = 1 / fs
        
        # Pad data symmetrically
        padvalue = n_orig // 2
        print(n_orig, " ds : ", ds, " dt : ", padvalue)
        
       
        x = np.concatenate((np.flipud(data[0:padvalue]), data, np.flipud(data[-padvalue:])))
        n = x.size
        # print(x)
        # Define scales
        _, _, wavscales = getDefaultScales(wavelet_name, n_orig, ds)
        num_scales = wavscales.size
        # Frequency vector sampling the Fourier transform of the wavelet
        omega = np.arange(1, math.floor(n / 2) + 1, dtype=np.float64)
        omega *= (2 * np.pi) / n
        
        omega = np.concatenate((np.array([0]), omega, -omega[np.arange(math.floor((n - 1) / 2), 0, -1, dtype=int) - 1]))
        print("np.array([0]) : ",np.concatenate((np.array([0]), omega)))
        # print("Omega : ",omega)
        # Compute FFT of the (padded) time series
        f = np.fft.fft(x)

        # Loop through all the scales and compute wavelet Fourier transform
        psift, freq = waveft(wavelet_name, omega, wavscales)

        # Inverse transform to obtain the wavelet coefficients.
        cwtcfs = np.fft.ifft(np.kron(np.ones([num_scales, 1]), f) * psift)
        cfs = cwtcfs[:, padvalue:padvalue + n_orig]
        freq = freq * fs

        return cfs, freq
    else:
        raise Exception


def getDefaultScales(wavelet, n, ds):
    """
    getDefaultScales(wavelet, n, ds)

    Calculate default scales given a wavelet and a signal length.

    Parameters
    ----------
    wavelet : string
        Name of wavelet
    n : int
        Number of samples in a given signal
    ds : float
        Scale resolution (inverse of number of voices in octave)

    Returns
    -------
    s0 : int
        Smallest useful scale
    ds : float
        Scale resolution (inverse of number of voices in octave). Here for legacy reasons; implementing more wavelets
        will need this output.
    scales : array_like
        Array containing default scales.
    """
    wname = wavelet
    nv = 1 / ds

    if wname == 'morl':

        # Smallest useful scale (default 2 for Morlet)
        s0 = 2

        # Determine longest useful scale for wavelet
        max_scale = n // (np.sqrt(2) * s0)
        
        if max_scale <= 1:
            max_scale = n // 2
        max_scale = np.floor(nv * np.log2(max_scale)) 
        a0 = 2 ** ds
        print("max a0 = ",a0)
        scales = s0 * a0 ** np.arange(0, max_scale + 1)
        # print((scales))
    else:
        raise Exception

    return s0, ds, scales


def waveft(wavelet, omega, scales):
    """
    waveft(wavelet, omega, scales)

    Computes the Fourier transform of a wavelet at certain scales.

    Parameters
    ----------
    wavelet : string
        Name of wavelet
    omega : array_like
        Array containing frequency values in Hz at which the transform is evaluated.
    scales : array_like
        Vector containing the scales used for the wavelet analysis.

    Returns
    -------
    wft : array_like
        (num_scales x num_freq) Array containing the wavelet Fourier transform
    freq : array_like
        Array containing frequency values
    """
    wname = wavelet
    num_freq = omega.size
    num_scales = scales.size
    wft = np.zeros([num_scales, num_freq])

    if wname == 'morl':
        gC = 6
        mul = 2
        for jj, scale in enumerate(scales):
            expnt = -(scale * omega - gC) ** 2 / 2 * (omega > 0)
            wft[jj, ] = mul * np.exp(expnt) * (omega > 0)

        fourier_factor = gC / (2 * np.pi)
        frequencies = fourier_factor / scales

    else:
        raise Exception

    return wft, frequencies
