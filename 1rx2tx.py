import adi
import time
import matplotlib.pyplot as plt
import numpy as np
from arrayfactor import af as af


sfasamento_sistema = -6          #deg


# Create radio
sdr = adi.ad9361(uri='ip:192.168.2.1')
samp_rate = 30e6    
num_samps = 2**12     
rx_lo = 2.4e9
rx_mode = "manual"
rx_gain0 = 10
rx_gain1 = -88
tx_lo = rx_lo
tx_gain = -3
c = 299792458  


sdr.rx_enabled_channels = [0, 1]
sdr.sample_rate = int(samp_rate)
sdr.rx_lo = int(rx_lo)
sdr.gain_control_mode = rx_mode
sdr.rx_hardwaregain_chan0 = int(rx_gain0)
sdr.rx_hardwaregain_chan1 = int(rx_gain1)
sdr.rx_buffer_size = int(num_samps)


sdr.tx_rf_bandwidth = int(samp_rate)
sdr.tx_lo = int(tx_lo)
sdr.tx_cyclic_buffer = True
sdr.tx_hardwaregain_chan0 = int(tx_gain)
sdr.tx_hardwaregain_chan1 = int(tx_gain)
sdr.tx_buffer_size = int(num_samps)


d_wavelength = 5              
wavelength = c/rx_lo            
d = d_wavelength*wavelength   
time_max = d / c  


def dbfs(raw_data):
    NumSamples = len(raw_data)
    win = np.hamming(NumSamples)
    y = raw_data * win
    s_fft = np.fft.fft(y) / np.sum(win)
    s_shift = np.fft.fftshift(s_fft)
    s_dbfs = 20*np.log10(np.abs(s_shift)/(2**11))  
    return s_dbfs

def dbm(raw_data):
    NumSamples = len(raw_data)
    #win = np.hamming(NumSamples)
    #y = raw_data * win
    s_fft = np.fft.fft(raw_data) #/ np.sum(win)
    s_shift = np.fft.fftshift(s_fft)
    s_db = 10*np.log10(np.abs(s_shift)/(2**11)) + 30 
    return s_db

# Example read properties
print("RX LO %s" % (sdr.rx_lo))

# Program the Tx with some data
fs = int(sdr.sample_rate)
fc0 = int(200e3)
N = 2**16
ts = 1 / float(fs)
t = np.arange(0, N * ts, ts)
i0 = np.cos(2 * np.pi * t * fc0) * 2 ** 14
q0 = np.sin(2 * np.pi * t * fc0) * 2 ** 14
iq0 = i0 + 1j * q0
sp = np.fft.fft(iq0)

dataFin = []
dataFinNotFiltered = []
result = af.af_asym_phasescannig(d_wavelength,0,0,tx_lo,tx_lo,sfasamento_sistema,2,0,0,0,'E')
xAxes = []
teorico = []

for i in range(-450,450):
    xAxes.append(i/5)
    teorico.append(result[1][900 + i])

for angle in range(-450,450):
    tf = (sp/2) * np.exp(1j * 2*np.pi*rx_lo*(np.sin(np.deg2rad(angle/5))*time_max))
    iq1 = np.fft.ifft(tf)*2
    #plt.clf()
    #plt.plot(iq0)
    #plt.plot(iq1)
    #plt.xlim(0,800)
    #plt.pause(0.1)
    #time.sleep(0.1)
    sdr.tx_destroy_buffer() 
    sdr.tx([iq0, iq1])   # Send Tx data.

    mean = 0
    for times in range(30):
        data = sdr.rx() 
        Rx_0=data[0]
        mean += np.max(dbm(Rx_0))/30
        if times == 0:
            dataFinNotFiltered.append(np.max(dbm(Rx_0)))

    dataFin.append(mean)

norm =  [(float(j)-np.min(dataFin))/np.max(dataFin-np.min(dataFin)) for j in dataFin]
norm1 = [(float(j)-np.min(dataFinNotFiltered))/np.max(dataFinNotFiltered-np.min(dataFinNotFiltered)) for j in dataFinNotFiltered]


plt.plot(xAxes,norm1)
plt.plot(xAxes,teorico)
plt.plot(xAxes,norm)
plt.legend(['Segnale ricevuto istantaneo', 'segnale teorico', 'Segnale mediato'])
plt.show()




