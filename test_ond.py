import numpy as np

#%matplotlib widget
from matplotlib import pyplot as plt
#from qdm.auto import trange, tqdm
import scipy.signal
from scipy import signal 
import scipy
from matplotlib import *
import pywt
from ssqueezepy import cwt
from ssqueezepy.visuals import plot, imshow
import pywt 

atest = pywt.scale2frequency('db2', 8, precision=8)


print('World')
print(pywt.scale2frequency('gaus1', [273]) / 40e6)

signal_bw = 160  # MHz
fft_len = 2**10

# Load from cfloat fie
sig0 = np.fromfile(r"C:\Users\HP\Desktop\wavelet\pasnoisy.iq", dtype=np.double)
za=[]
zb=[]
zz = []
for i in range(0,len(sig0)):
    if i  % 2 != 0:
        zb.append(sig0[i])
    else:
        za.append(sig0[i])
        
zza = np.array(za)
zzb = np.array(zb)
zzz = zza + 1j*zzb
sig0 = zzz
plt.specgram(sig0,Fs=40e6)
plt.show()
#sig0 = sig0.astype(np.complex64) # Convert to 64
#sig0.tofile('bpsk_in_noise.txt') # Save to file



# sig0 = np.fromfile(r"C:\Users\HP\Desktop\wavelet\bpsk_in_noise.iq", dtype=np.complex64)
Zxx = (np.fft.fft(sig0))
fft0 = np.square(np.abs(Zxx))
# fft_mean = 10*np.log10(np.mean(fft0, axis=0)/(fft_len))  # Should be a time average

xf = np.fft.fftfreq(len(sig0))
# plt.plot(xf,fft0.real**2+fft0.imag**2)
# plt.show()  

# =============================================================================
# From ssqueezepy
# Wx, scales = cwt(sig0, 'morlet')
# =============================================================================

#from scipy
witdhs = np.arange(2, 10000)
scales = [2**4, 2**6, 2**8]  # Wavelet scales

# cwtmatr = signal.cwt(sig0, signal.morlet, widths)
Zxx = (np.fft.fft(sig0))
freqs = np.square(np.abs(Zxx))
Wx, scaless = cwt(freqs/np.max(freqs), 'morlet')
 

[cfs, frequencies] = pywt.cwt(freqs/np.max(freqs), scales, 'gaus1')

imshow(cfs*np.abs(cfs), yticks=scales, abs=1,
       title="abs(CWT) | Morlet wavelet",
       ylabel="scales", xlabel="samples")

imshow(Wx, yticks=scaless, abs=1,
       title="abs(CWT) | Morlet wavelet",
       ylabel="scales", xlabel="samples")