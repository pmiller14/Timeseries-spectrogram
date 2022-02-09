#plot spectrogram of acoustic waveforms
#import relevant packages
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import scipy 
from scipy import signal 
from collections import OrderedDict
#This is different for python 3
from scipy.signal import kaiserord, lfilter, firwin, freqz

#Create frequency content class
class fc():
    def __init__(self):
        self.name = 'makeplot'
        
    #function to load variables
    def loadvar(self,max_ind,stress,waveform_ar,time,freqcont,freqvec,fs): #stress,    
        self.Waveform_array = waveform_ar  
        self.max_index = max_ind   
        self.stress = stress
        self.time = time
        self.frequency_content = freqcont
        self.freq_vector= freqvec
        self.fs = fs								
    
    #plot figure    
    def plotfigure(self,n):
    #Set local variable for waveform number
    	self.n=n
        fig = plt.figure(figsize=(13,10))
        #create gridspec instance
        gs = gridspec.GridSpec(2, 2,width_ratios=[1,3],height_ratios=[1,2]) 
        
        ax = fig.add_subplot(gs[1])
        ax.plot(self.time[0:self.max_index],self.Waveform_array[0:self.max_index,self.n],label = '%1.1f MPa'%float(self.stress[self.n]), linewidth = 2.5,c='0.05')
        ax.set_xlim([0,120])
        plt.legend(loc='upper right', fontsize = 14)
        ax.set_ylabel('Amplitude [$mV$]',fontsize = 20,labelpad=-8)
        plt.yticks(np.linspace(np.min(self.Waveform_array[:self.max_index,self.n]), np.max(self.Waveform_array[:self.max_index,self.n]), 5))
        ax.tick_params(axis='y', which='major', labelsize=18)
        ax.tick_params(axis='x', which='major', labelsize=18)

        ax = fig.add_subplot(gs[2])
        ax.plot(self.frequency_content[:self.freq_vector.size,self.n]/np.max(self.frequency_content[:self.freq_vector.size,self.n]),self.freq_vector/1e3, label = '50 MPa' ,linewidth=1.5)
        ax.invert_xaxis()
        ax.set_ylim([0,1000])
        ax.set_xlabel('Normalized Amplitude',fontsize = 20)
        ax.tick_params(axis='x', which='major', labelsize=18)
        plt.xticks(np.linspace(0, 1.0, 3))
        ax.axes.get_yaxis().set_visible(False)

        ax = fig.add_subplot(gs[3])
        self.f, self.t, self.Sxx = signal.spectrogram(self.Waveform_array[:self.max_index,self.n], self.fs,window='hamming', nfft=1500, noverlap=255)
        im=ax.pcolormesh(self.t*1e6, self.f/1e3, self.Sxx)
        fig.colorbar(im)
        ax.set_ylabel('frequency [$kHz$] ',fontsize = 20,labelpad=-12)
        ax.set_xlabel('Time [$\mu s$]',fontsize = 20)
        ax.set_ylim([0,1e3])
        ax.set_xlim([0,120])
        ax.tick_params(axis='y', which='major', labelsize=18)
        ax.tick_params(axis='x', which='major', labelsize=18) 

        return self.f, self.t, self.Sxx
    
    def findamp(self,f,t,Sxx,tmin,tmax):
    	self.tmin=tmin/1e6
    	self.tmax=tmax/1e6
    	self.f,self.t,self.Sxx = f,t,Sxx
    	tminidx = (np.abs(self.t - self.tmin)).argmin()
    	tmaxidx = (np.abs(self.t - self.tmax)).argmin()
    	self.tseg=self.t[tminidx:tmaxidx]
    	x,y = np.where(self.Sxx[:,tminidx:tmaxidx]==np.max(self.Sxx[:,tminidx:tmaxidx]))
    	return self.f[x]/1e3,self.tseg[y]*1e6,np.max(self.Sxx[:,tminidx:tmaxidx])
		
    

