#	Script name: 	fourtran.py
#	Written by:		Troels Suhr Skovgaard
#	Date: 			25.05.2012
#
#	Takes equividistant time series and Fourier transforms using either fast-fourier
#	transform, continous fourier transform or harmonic inversion. Outputs transformed data
#	and estimated Q factors.
#

import argparse
import numpy as np
from scipy import optimize
#from scipy import *
from math import *
import matplotlib.pyplot as plt

def FT(t,x):
	# Continous Fourier Transform module
	
	print "-"*40 + "\n"
	print "Fourier transforming data"
	print "Using continouos transform\n"
	
	dt   = t[1]-t[0]
	Nf   = 250
	freq = np.linspace(0.7,1.3,Nf)
	fx   = np.zeros(np.size(freq),dtype=complex)
	
	print "Progress:"

	for i,f in enumerate(freq):
		if i%(Nf/10) is 0:
			print str(int(100*i/Nf)) + "%"
		fx[i] = dt*sum(x*np.exp(2*pi*1j*f*t))
	
	print "100%"
	print "Done\n"
	print "-"*40 + "\n"
	
	plt.figure()
	plt.plot(freq,abs(fx)**2)
	plt.xlabel('frequency')
	plt.ylabel('amplitude')
	plt.title('Continouos Fourier Transform')
	plt.show()
	
	return freq,fx

def FFT(t,x):
	# Fast Fourier Transform module
	
	print "-"*40 + "\n"
	print "Fourier transforming data"
	print "Using fast fourier transform\n"
	
	fx   = np.fft.fft(x)
	freq = np.fft.fftfreq(t.shape[-1], d=t[1]-t[0])
	
	N    = len(fx)
	fx   = fx[:int(N/2)]
	freq = freq[:int(N/2)]
	
	print "Done\n"
	print "-"*40 + "\n"
	
	plt.figure()
	plt.plot(freq,abs(fx)**2)
	plt.xlabel('frequency')
	plt.ylabel('amplitude')
	plt.title('Fast Fourier Transform')
	plt.show()
	
	return freq,fx

def HI():
	# Harmonic Inversion module
	
	print "-"*40 + "\n"
	print "Fourier transforming data"
	print "Using harmonic inversion method\n"
	
	print "Not yet implemented\n"
	
	print "-"*40 + "\n"	

def findQ(freq, fx, peaks=1):
	# Fits lorentzian to input data using Scipy optimize function
	
	print "-"*40 + "\n"
	print "Fitting lorentzian to spectrum\n"
	
	fitfunc = lambda p, x: p[0]*(p[1]/(4*pi*p[2]))**2/((x-p[1])**2 + (p[1]/(4*pi*p[2]))**2)
	errfunc = lambda p, x, y: fitfunc(p,x) - y
	
	P = []
	for i,p in enumerate(peaks):
		p0 = [1e6, p, 1000]
		p, success = optimize.leastsq(errfunc, p0[:], args=(freq, fx))
		P.append(p)
		
	if len(P) > 1:
		fmin  = min(P[i][1] for i in range(len(P)))
		fmax  = max(P[i][1] for i in range(len(P)))
		fmean = (fmin + fmax)/2
		df    = (fmax - fmin)/2
		
	if len(P) is 1 or df < 1e-2:
		fmean = P[0][1]
		df    = 0.2*fmean
	
	plt.figure()
	plt.plot(freq, fx, 'b')
	plt.xlim([fmean - 1.1*df,fmean + 1.1*df])
	plt.xlabel('freq')
	plt.ylabel('amplitude')
	plt.title('Fitted data')
	plt.hold('on')
	for i in range(len(P)):
		plt.plot(freq,fitfunc(P[i],freq),'rx')
	plt.show()
	
	print "Found resonances:"
	for i,v in enumerate(P):
		print "Peak No " + str(i)
		print "Amplitude\t = " + str(v[0])
		print "Resonance\t = " + str(v[1])
		print "Q-factor\t = "  + str(abs(v[2])) + "\n"
	
	print "-"*40 + "\n"
	
	return p
	
def findpeaks(freq,fx):
	# Identifies peaks of a given data by inspecting the 1st derivative
	
	print "-"*40 + "\n"
	print "Identifying spectrum peaks\n"
	
	dfx = np.diff(fx)
	
	npeaks = []
	for i in range(1,len(dfx)-1):
		if dfx[i] > dfx[i-1] and dfx[i] > dfx[i+1]:
			npeaks.append(i+2)
	
	fpeaks = freq[npeaks]
	
	plt.figure()
	plt.plot(freq,fx,'b',fpeaks,fx[npeaks],'ro')
	plt.xlabel('frequency')
	plt.ylabel('amplitude')
	plt.title('Peaks')
	plt.show()
	
	print "-"*40 + "\n"
	print "Found " + str(len(fpeaks)) + " peaks:"
	print "Peak\tfrequency"
	for i,p in enumerate(fpeaks):
		print str(i) + "\t" + str(p)
	
	print "\n"	
	print "-"*40 + "\n"
	
	return fpeaks

def main():
	
	plt.close('all')
	
	t = np.linspace(0,1e3,2e4)
	#x = np.sin(2*pi*t)*np.exp(-0.01*t)
	x = np.sin(2*pi*t)*np.exp(-0.01*t) + np.sin(1.01*2*pi*t)*np.exp(-0.01*t)
	
	plt.figure()
	plt.plot(t,x)
	plt.xlabel('time')
	plt.ylabel('amplitude')
	plt.title('Time series')
	plt.show()
	
	#freq,fx = FT(t,x)
	
	freq2,fx2 = FFT(t,x)
	
	fpeaks = findpeaks(freq2,abs(fx2)**2)
	
	#A,res,Q = findQ(freq,abs(fx)**2)
	
	A2,res2,Q2 = findQ(freq2,abs(fx2)**2,peaks=fpeaks)

if __name__ == '__main__':
	main()