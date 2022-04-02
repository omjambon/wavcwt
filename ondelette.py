import numpy as np
import pywt

#%matplotlib widget
from matplotlib import pyplot as plt
#from qdm.auto import trange, tqdm
import scipy.signal
import scipy
from ssqueezepy.visuals import plot, imshow
from scipy import signal

#from scipy.ndimage.filters import uniform_filter1d
#from multiprocessing import Pool
#import cv2
#import imutils
print('World')


signal_bw = 160  # MHz
fft_len = 2**10


# Load from cfloat file
sig0 = np.fromfile(r"C:\Users\HP\Desktop\wavelet\bpsk_in_noise.iq", dtype=np.complex64)
#sig0 = sig0.astype(np.complex64) # Convert to 64
#sig0.tofile('bpsk_in_noise.txt') # Save to file
# Load from SigMF file
# signal = sigmffile.fromfile(os.path.join(folder, sigmf_files[0]))
# file_sample_count = signal.sample_count
# sig0 = signal.read_samples(0, -1)
# annotations = signal.get_annotations()


#plt.figure(figsize=[12, 5])
Zxx = (np.fft.fft(sig0))
fft0 = np.square(np.abs(Zxx))
fft_mean = 10*np.log10(np.mean(fft0, axis=0)/(fft_len))  # Should be a time average

xf = np.fft.fftfreq(len(sig0))
plt.plot(xf,fft0.real**2+fft0.imag**2)
plt.show()  

scales = [2**4, 2**6, 2**8]  # Wavelet scales
scaling = np.array([[2**1], [2**1], [2**0]])  # Weighting of scales
side_cut = 1500  # Limit to peak detection (reduce side artefacts)

freqs = fft0

[cfs, frequencies] = pywt.cwt(freqs/np.max(freqs), scales, 'gaus1')

cfs *= scaling 
sum_cwt = np.sum(cfs*np.abs(cfs), axis=0)  # We want to square the signal to better distinguish peaks from noise, but keeping sign information to get direction

# Peak detection (Peaks and trophs)
peaks, _ = scipy.signal.find_peaks(sum_cwt[side_cut:-side_cut], prominence=0.5)
neg_peaks, _ = scipy.signal.find_peaks(-sum_cwt[side_cut:-side_cut], prominence=0.5)
peaks += side_cut  # Translate back to original place
neg_peaks += side_cut
occupied = np.zeros_like(freqs, dtype=bool)
id_l1 = 0
#id_l2 = 0

## Use the in and out edges (sign information) to fill in occupancy
#while id_l1 < len(peaks) and id_l2 < len(neg_peaks):
#    if peaks[id_l1] < neg_peaks[id_l2]:
#        occupied[peaks[id_l1]:neg_peaks[id_l2]] = False
#        id_l1 += 1
#    elif peaks[id_l1] == neg_peaks[id_l2]:
#        print("Shouldn't happen!!!")
#    else:
#        occupied[neg_peaks[id_l2]:peaks[id_l1]] = True
#        id_l2 += 1
#    if id_l1 == len(peaks):
#        occupied[neg_peaks[id_l2]:] = False
#    if id_l2 == len(neg_peaks):
#        occupied[peaks[id_l1]:] = False
widths = np.arange(1, 4)
cwtmatr = signal.cwt(freqs/np.max(freqs), signal.morlet2, widths)
imshow(cwtmatr *np.abs(cwtmatr ), yticks=scales, abs=1,
       title="abs(cwtmatr ) | Morlet wavelet",
       ylabel="scales", xlabel="samples")

imshow(cfs*np.abs(cfs), yticks=scales, abs=1,
       title="abs(CWT) | Morlet wavelet",
       ylabel="scales", xlabel="samples")




import scipy.io
mdic = {"sig0": sig0, "label": "experiment"}
scipy.io.savemat('sig0.mat',  mdict={'sig0': sig0})

