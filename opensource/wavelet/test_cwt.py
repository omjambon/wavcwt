import numpy as np 
import matplotlib.pyplot as plt
import cwt

from scipy import signal 

fs = 4e7
t = np.linspace(0, 1, fs+1, endpoint=True)
x = np.cos(2*np.pi*32*t) * np.logical_and(t >= 0.1, t < 0.3) + np.sin(2*np.pi*64*t) * (t > 0.7)
wgnNoise = 0.05 * np.random.standard_normal(t.shape)
#x += wgnNoise



sig0 = np.fromfile("/home/rpd316/Documents/Ababil3/EES300821/filtresolodaltons.iq", dtype=np.complex64)

c, f = cwt.cwt(sig0, 'morl', sampling_frequency=fs)

fig, ax = plt.subplots()
ax.imshow(np.absolute(c), aspect='auto')



#==============================================================================
# 
#==============================================================================

cwtmatr = signal.cwt(sig0,signal.morlet,np.arange(1,365))
fig1, ax1 = plt.subplots()
ax1.imshow(np.absolute(cwtmatr), aspect='auto')