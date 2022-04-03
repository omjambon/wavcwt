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
import cwt

Zxx = (np.fft.fft(sig0))
freqs = np.square(np.abs(Zxx))
cfs, f = cwt.cwt(freqs/np.max(freqs), 'morl', sampling_frequency=40e6)
scales = [2**4, 2**6, 2**8]  # Wavelet scales
imshow(cfs*np.abs(cfs), yticks=scales, abs=1,
       title="abs(CWT) | Morlet wavelet",
       ylabel="scales", xlabel="samples")