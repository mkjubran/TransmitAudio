
import matplotlib.pyplot as plt 
import numpy as np 
import random 
import scipy 
import wave 
import struct 
import pylab
import pdb
from scipy.io import wavfile

################################################
filenameWave='./summer.wav'
#filenameWave='./counting.wav'
#filenameWave='./athan1.wav'
#filenameWave='./sunnyDay.wav'
filenameWavefiltered='./Filtered.wav'
filenameWavewithoutfilter='./Withoutfilter.wav'
PowerBW=0.85
write=0
read=0
################################################

rate, data = wavfile.read(filenameWave)
#pdb.set_trace()
if len(data.shape) > 1:
   data=data[:,0]

filtereddata = np.fft.rfft(data, axis=0)

filteredwrite = np.fft.irfft(filtereddata, axis=0)
##############################################
## Generate Signal and add noise to it
Fs=rate;
Ts = 1.0/Fs; # sampling interval
t = np.arange(0,1,Ts) # time vector
t = np.arange(0,len(data)*Ts,Ts) # time vector

amplitude = 30
noise = amplitude * np.random.normal(scale=0.1, size=len(t))
y=[float(x) for x in data]
## Write values to a file
#Open new data file
if write!=0:
   f = open("Signal_in_text.txt", "w")
   for i in range(len(y)):
       f.write( str(y[i]) + " " + str(float(t[i])) + "\n"  )
   f.close()


## Read values from a file
if read !=0:
   with open('Signal_in_text.txt') as f:
      w=f.read()
   y=[];
   t=[];
   for x in w.split('\n'):
      if x != '':
         y.append(float(x.split()[0]))
         t.append(float(x.split()[1]))

n = len(y) # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T # two sides frequency range
fcen=frq[len(frq)/2]
frq_DS=frq-fcen
frq_SS = frq[range(n/2)] # one side frequency range

Y = np.fft.fft(y) # fft computing and normalization
yinv= np.fft.ifft(Y).real # ifft computing and normalization
Y_DS=np.roll(Y,n/2)
Y_SS = Y[range(n/2)]

fig, ax = plt.subplots(7, 1)
ax[0].plot(t,y)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[1].plot(frq_SS,abs(Y_SS),'r') # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[2].plot(frq_DS[n/4:3*n/4],abs(Y_DS[n/4:3*n/4]),'r') # plotting the spectrum
ax[2].set_xlabel('Freq (Hz)')
ax[2].set_ylabel('|Y(freq)|')
ax[3].plot(t,yinv,'g') # plotting the spectrum
ax[3].set_xlabel('Time')
ax[3].set_ylabel('Amplitude')

#plt.show()
y=np.array(y)
y_int=y.astype(np.int16)

yinv=np.array(yinv)
yinv_int=yinv.astype(np.int16)
wavfile.write(filenameWavewithoutfilter, rate, y_int)

Energy_time=np.sum(data.astype(float)**2)
Power_time=1.0/(2*(data.size)+1)*np.sum(data.astype(float)**2)/rate
Energy_freq_DS=np.real(sum(Y_DS*np.conj(Y_DS)))/n
Power_freq_DS=np.real(sum(Y_DS*np.conj(Y_DS)))/(2*(n**2)/(Ts))
print('Duration={} sec'.format(len(y)*Ts))
print('Energy of wave (time domain)={} KJouls'.format(Energy_time/1000))
print('Power of wave (time domain)={} Watts'.format(Power_time))
print('Energy of wave (DS frequency domain)={} KJouls'.format(Energy_freq_DS/1000))
print('Power of wave (DS frequency domain)={} Watts'.format(Power_freq_DS))

Mask_DS=np.ones(len(frq_DS))
Yf_DS=np.copy(Y_DS)
Bmax=frq_DS[len(frq_DS)-1]
Bmin=0
B=Bmax
Bold=0
Energy_time_filtered=Energy_time
while (abs(Bold - B) > 0) and B > 1:
   print('\n============ BW={} Hz'.format(B))
   del yinv_int
   Mask_DS=np.ones(len(frq_DS))
   Yf_DS=np.copy(Y_DS)
   for cnt in range(len(frq_DS)):
      if ~(((frq_DS[cnt])>-1*B) and ((frq_DS[cnt])<B)):
          Mask_DS[cnt]=0;
          Yf_DS[cnt]=Y_DS[cnt]*0;

   Yf=np.roll(Yf_DS,n/2)
   yinv= np.fft.ifft(Yf).real # ifft computing and normalization
   yinv=np.array(yinv)
   yinv_int=yinv.astype(np.int16)

   Energy_time_filtered=np.sum(yinv.astype(float)**2)
   Power_time_filtered=1.0/(2*(yinv.size)+1)*np.sum(yinv.astype(float)**2)/rate
   Energy_freq_DS_filtered=np.real(sum(Yf_DS*np.conj(Yf_DS)))/n
   Power_freq_DS_filtered=np.real(sum(Yf_DS*np.conj(Yf_DS)))/(2*(n**2)/(Ts))
   print('Energy of the filtered wave (time domain)={} KJouls'.format(Energy_time_filtered/1000))
   print('Power of the filtered wave (time domain)={} Watts'.format(Power_time_filtered))
   print('Energy of the filtered wave (DS frequency domain)={} KJouls'.format(Energy_freq_DS_filtered/1000))
   print('Power of the filtered wave (DS frequency domain)={} Watts'.format(Power_freq_DS_filtered))
   print('Energy Ratio (filtered / original) = {}%'.format(((Energy_time_filtered/Energy_time))*100))

   ax[4].clear()
   ax[4].plot(frq_DS[n/4:3*n/4],abs(Mask_DS[n/4:3*n/4]),'r') # plotting the spectrum
   ax[4].set_xlabel('Freq (Hz)')
   ax[4].set_ylabel('|Y(freq)|')
   ax[5].clear()
   ax[5].plot(frq_DS[n/4:3*n/4],abs(Yf_DS[n/4:3*n/4]),'r') # plotting the spectrum
   ax[5].set_xlabel('Freq (Hz)')
   ax[5].set_ylabel('|Y(freq)|')
   ax[6].clear()
   ax[6].plot(t,yinv,'g') # plotting the spectrum
   ax[6].set_xlabel('Time')
   ax[6].set_ylabel('Amplitude')
   plt.pause(0.3)

   if (Energy_time_filtered > PowerBW*Energy_time):
      Bmax=B;
   else:
      Bmin=B;
   Bold=B
   B=Bmin+(Bmax-Bmin)/2

wavfile.write(filenameWavefiltered, rate, yinv_int)

ax[4].plot(frq_DS[n/4:3*n/4],abs(Mask_DS[n/4:3*n/4]),'r') # plotting the spectrum
ax[4].set_xlabel('Freq (Hz)')
ax[4].set_ylabel('|Y(freq)|')
ax[5].plot(frq_DS[n/4:3*n/4],abs(Yf_DS[n/4:3*n/4]),'r') # plotting the spectrum
ax[5].set_xlabel('Freq (Hz)')
ax[5].set_ylabel('|Y(freq)|')
ax[6].plot(t,yinv,'g') # plotting the spectrum
ax[6].set_xlabel('Time')
ax[6].set_ylabel('Amplitude')

Energy_time=np.sum(data.astype(float)**2)
Power_time=1.0/(2*(data.size)+1)*np.sum(data.astype(float)**2)/rate
Energy_freq_DS=np.real(sum(Y_DS*np.conj(Y_DS)))/n
Power_freq_DS=np.real(sum(Y_DS*np.conj(Y_DS)))/(2*(n**2)/(Ts))
Energy_time_filtered=np.sum(yinv.astype(float)**2)
Power_time_filtered=1.0/(2*(yinv.size)+1)*np.sum(yinv.astype(float)**2)/rate
Energy_freq_DS_filtered=np.real(sum(Yf_DS*np.conj(Yf_DS)))/n
Power_freq_DS_filtered=np.real(sum(Yf_DS*np.conj(Yf_DS)))/(2*(n**2)/(Ts))

print('\n\n============ Summary ======'.format(B))
print('Duration = {} sec'.format(len(y)*Ts))
print('Energy of wave (time domain) = {} KJouls'.format(Energy_time/1000))
print('Power of wave (time domain) = {} Watts'.format(Power_time))
print('Energy of wave (DS frequency domain) = {} KJouls'.format(Energy_freq_DS/1000))
print('Power of wave (DS frequency domain) = {} Watts'.format(Power_freq_DS))
print('Bandiwdth of filtered wave = {} KHz'.format(float(B)/1000))
print('Energy of the filtered wave (time domain) = {} KJouls'.format(Energy_time_filtered/1000))
print('Power of the filtered wave (time domain) = {} Watts'.format(Power_time_filtered))
print('Energy of the filtered wave (DS frequency domain) = {} KJouls'.format(Energy_freq_DS_filtered/1000))
print('Power of the filtered wave (DS frequency domain) = {} Watts'.format(Power_freq_DS_filtered))
print('Energy Ratio (filtered / original) = {}%'.format(((Energy_time_filtered/Energy_time))*100))

#plt.show()
