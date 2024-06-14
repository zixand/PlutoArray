import adi
import time
import matplotlib.pyplot as plt
import numpy as np
from arrayfactor import af as af


sfasamento_sistema = -3          #deg


# Create radio
sdr = adi.ad9361(uri='ip:192.168.2.1')
samp_rate = 30e6    
num_samps = 2**12     
rx_lo = 2.4e9
rx_mode = "manual"
rx_gain0 = 20
rx_gain1 = 20
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
sdr.tx_hardwaregain_chan1 = int(-88)
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
sdr.tx([iq0, iq0])   # Send Tx data.


result = af.af_asym_phasescannig(d_wavelength,0,0,tx_lo,tx_lo,sfasamento_sistema,2,0,0,0,'E')




# Collect data
for r in range(20):    # grab several buffers to give the AGC time to react (if AGC is set to "slow_attack" instead of "manual")
    data = sdr.rx()


for i in range(1000):


    dataFin = []
    data = sdr.rx() 
    Rx_0=data[0]
    Rx_1=data[1]


    NumSamples = len(Rx_1)
    y = Rx_1 
    sp = np.fft.fft(y)


    for angle in range(-450,450):
        tf = (sp/2) * np.exp(1j * 2*np.pi*rx_lo*(np.sin(np.deg2rad(angle/5))*time_max))
        tf = np.fft.ifft(tf)*4
        Rx_total = tf + Rx_0
        if i == 0:
            plt.clf()
            plt.plot(Rx_total)
            plt.plot(Rx_0)
            plt.plot(Rx_1)
            plt.plot(tf)
            plt.xlim(0,800)
            plt.ylim(-2000,2000)
            plt.legend(['Somma Segnali', 'Antenna 1', 'Antenna 2', 'Antenna 2 (Traslata)'])
            plt.pause(0.05)
            time.sleep(0.05)
        Rx_total = dbm(Rx_total)
        dataFin.append(np.max(Rx_total))


    xAxes = []
    teorico = []

    for i in range(-450,450):
        xAxes.append(i/5)
        teorico.append(result[1][900 + i])

    dataFinNorm = [float(z) - np.min(dataFin) for z in dataFin]
    norm =  [float(j)/max(dataFinNorm) for j in dataFinNorm]


    plt.clf()
    plt.plot(xAxes, norm)
    #plt.plot(xAxes, dataFin)
    plt.plot(xAxes, teorico)
    plt.legend(['Sperimentale', 'Teorico'])
    plt.grid(True)
    plt.xlim(95,-95)
    plt.pause(0.05)

plt.show()




sp = sp[1:-1]
sp = np.fft.fftshift(sp)
s_mag = np.abs(sp) /2   #
s_dbfs = 20*np.log10(s_mag/(2**12))     


xf = np.fft.fftfreq(NumSamples, ts)
xf = np.fft.fftshift(xf[1:-1])/1e6

plt.plot(xf, s_dbfs)
plt.xlabel("frequency [MHz]")
plt.ylabel("dBfs")
plt.draw()
plt.show()








