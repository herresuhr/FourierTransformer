#	Script name: 	fourtran.py
#	Written by:		Troels Suhr Skovgaard
#	Date: 			25.05.2012
#
#	Takes equividistant time series and Fourier transforms using either fast-fourier
#	transform, continous fourier transform or harmonic inversion. Outputs transformed data
#	and estimated Q factors.
#
#   Added routine for fitting multiple peaks to a truncated harmonic oscillation. 
#   Generally finds peaks nicely, but fitting routine occasionally greatly overestimates 
#   Q if the fitted FWHM is below the spectral resolution.

import argparse
import numpy as np
from scipy import optimize
#from scipy import *
from math import *
import matplotlib.pyplot as plt
import movingaverage as ma

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

def findTQ(freq, fx, peaks=1, t=-1):
	# Fits truncated lorentzian to input data using Scipy optimize function
	
	print "-"*40 + "\n"
	print "Fitting lorentzian to spectrum\n"
	
	if t is -1:
		t = 1/(freq[1]-freq[0])
	
	fitfunc = lambda p, x: p[0]/4.0*(1 - 2*np.exp(-abs(p[2])*t)*np.cos((x-p[1])*t) + np.exp(-2*abs(p[2])*t))/((x-p[1])**2 + p[2]**2)
	errfunc = lambda p, x, y: (fitfunc(p,x) - y)
	
	if type(peaks) is not np.ndarray:
		peaks = np.array([peaks])
	
	P = {}
	for i in range(len(peaks)):
		p0 = [1e3, peaks[i], 1e-3]

		p, success = optimize.leastsq(errfunc, p0[:], args=(freq, fx))

		P[i] = (p)
		print t*p[2]
		
	fmean = P[0][1]
	df = 0.2*fmean
				
	if len(P) > 1:
		fmin  = min(P[i][1] for i in P.keys())
		fmax  = max(P[i][1] for i in P.keys())
		fmean = (fmin + fmax)/2
		df    = (fmax - fmin)/2
		if df < 1e-2:
			df = 0.2*fmean
		
	ff = np.linspace(fmean - 2*df,fmean + 2*df,1e4)
	
	plt.figure()
	plt.plot(freq, fx, 'bx')
	plt.xlim([fmean - 2*df,fmean + 2*df])
	plt.xlabel('freq')
	plt.ylabel('amplitude')
	plt.title('Fitted data')
	plt.hold('on')
	for i in P.keys():
		plt.plot(ff,fitfunc(P[i],ff),'r')
	plt.show()
	
	print "Found resonances:\n"
	for key in P.keys():
		print "Peak No " + str(key)
		print "Amplitude\t = " + str(P[key][0])
		print "Resonance\t = " + str(P[key][1])
		print "Q-factor\t = "  + str(abs(P[key][1]/(2*pi*P[key][2]))) + "\n"
	
	print "-"*40 + "\n"
	
	return P

def findQ(freq, fx, peaks=1):
	# Fits lorentzian to input data using Scipy optimize function
	
	print "-"*40 + "\n"
	print "Fitting lorentzian to spectrum\n"
	
	fitfunc = lambda p, x: p[0]**2*(p[1]/(2*p[2]))**2/((x-p[1])**2 + (p[1]/(2*p[2]))**2)
	errfunc = lambda p, x, y: (fitfunc(p,x) - y)*fitfunc(p,x)
	
	if type(peaks) is not np.ndarray:
		peaks = np.array([peaks])
	
	P = {}
	for i in range(len(peaks)):
		p0 = [1e3, peaks[i], 100]

		p, success = optimize.leastsq(errfunc, p0[:], args=(freq, fx))
		P[i] = (p)
		
	for key in P.keys():
		if abs(P[key][0]) < 1e-5:
			print P[key][0]
			P.pop(key)
				
	if len(P) is 1:
		fmean = P[0][1]
		df    = 0.2*fmean

	if len(P) > 1:
		fmin  = min(P[i][1] for i in range(len(P)))
		fmax  = max(P[i][1] for i in range(len(P)))
		fmean = (fmin + fmax)/2
		df    = (fmax - fmin)/2
		if df < 1e-2:
			df = 0.2*fmean
		
	

	ff = np.linspace(0,17,1e4)
	
	plt.figure()
	plt.plot(freq, fx, 'bx')
	plt.xlim([fmean - 1.1*df,fmean + 1.1*df])
	plt.xlabel('freq')
	plt.ylabel('amplitude')
	plt.title('Fitted data')
	plt.hold('on')
	for i in range(len(P)):
		plt.plot(ff,fitfunc(P[i],ff),'rx')
	plt.show()
	
	print "Found resonances:\n"
	for key in P.keys():
		print "Peak No " + str(key)
		print "Amplitude\t = " + str(P[key][0])
		print "Resonance\t = " + str(P[key][1])
		print "Q-factor\t = "  + str(abs(P[key][2])) + "\n"
	
	print "-"*40 + "\n"
	
	return P

def findMTQ(freq, fx, peaks=1, t=-1):
	# Fits truncated double lorentzian to input data using Scipy optimize function
	
	print "-"*40 + "\n"
	print "Fitting multiple truncated lorentzians to spectrum\n"

	if t is -1:
		t = 1/(freq[1]-freq[0])
	
	def fitfunc(p,x):
		N = len(p)/3
		
		# Pull out amplitudes, resonance frequencies and decay constants
		# Absolute value due to fitting crossing zero
		A = abs(p[0::3])
		w = 2*pi*abs(p[1::3])
		a = 2*pi*abs(p[2::3])
		x = 2*pi*x
	
		fit = 0
		for n in range(N):
			for m in range(N):
				NR = (a[n]*a[m]+(x-w[n])*(x-w[m]))*(1 + np.exp(-(a[n]+a[m])*t)*np.cos((w[m]-w[n])*t) - np.exp(-a[n]*t)*np.cos((x-w[n])*t) - np.exp(-a[m]*t)*np.cos((x-w[m])*t))
				NI = (a[n]*(x-w[m]) - a[m]*(x-w[n])) * (np.exp(-(a[n]+a[m])*t)*np.sin((w[m]-w[n])*t) - np.exp(-a[n]*t)*np.sin((x-w[n])*t) + np.exp(-a[m]*t)*np.sin((x-w[m])*t))
				DN = (a[n]*a[m] + (x-w[n])*(x-w[m]))**2 + (a[n]*(x-w[m])-a[m]*(x-w[n]))**2
				fit += A[n]*A[m]/2.0*a[n]*a[m]*(NR+NI)/DN
				
		return fit
		
	errfunc = lambda p, x, y: (fitfunc(p,x) - y)
	
	if type(peaks) is not np.ndarray:
		peaks = np.array([peaks])
	
	p0 = 3*len(peaks)*[0]
	p0[0::3] = len(peaks)*[1e3]
	for n in range(len(peaks)): 
		p0[1+3*n] = peaks[n]
	p0[2::3] = len(peaks)*[1e-3]
	
	P, success = optimize.leastsq(errfunc, p0[:], args=(freq, fx))

	#P_max = max(P[::3])
	#P[::3] = P[::3]/P_max
	#fx = fx/P_max
		
	fmean = (peaks[-1] + peaks[0])/2
	df    = (peaks[-1] - peaks[0])/2
		
	ff = np.linspace(fmean - 2*df,fmean + 2*df,1e4)
	
	plt.figure()
	plt.plot(freq, fx, 'bx',ff,fitfunc(P,ff),'r')
	plt.xlim([fmean - 2*df,fmean + 2*df])
	plt.xlabel('freq [c/a]')
	plt.ylabel('amplitude')
	plt.title('Fitted data')
	plt.show()
	
	print "Found resonances:\n"
	for n in range(len(peaks)):
		print "Peak No " + str(n)
		print "Amplitude\t = " + str(round(P[3*n],4))
		print "Resonance\t = " + str(round(P[3*n+1],4))
		print "Q-factor\t = "  + str(round(abs(P[3*n+1]/P[3*n+2]),1)) + "\n"
	
	print "-"*40 + "\n"
	
	return P	
	
def findpeaks(freq,fx,ran=-1):
	# Identifies peaks of a given data by inspecting the 1st derivative
	
	print "-"*40 + "\n"
	print "Identifying spectrum peaks\n"
	
	if ran is -1:
		ran = [0,freq[-1]]
	
	dfx = np.diff(fx)
	
	npeaks = []
	for i in range(1,len(dfx)):
		if freq[i] >= ran[0] and freq[i] <= ran[1] and dfx[i-1] > 0 and dfx[i] < 0:
			npeaks.append(i)
	
	A_max = max(fx[npeaks])
	tmp = npeaks
	npeaks = []
	for n in tmp:
		if fx[n]/A_max > 1e-2:
			npeaks.append(n)
	
	peaks = freq[npeaks]

#	plt.figure()
#	plt.plot(freq,fx,'b',fpeaks,fx[npeaks],'ro')
#	plt.xlabel('frequency')
#	plt.ylabel('amplitude')
#	plt.title('Peaks')
#	plt.show()
	
	print "Found " + str(len(peaks)) + " peaks:\n"
	print "Peak\tfrequency"
	for i,p in enumerate(peaks):
		print str(i) + "\t" + str(round(p,4))
	
	print "\n"	
	print "-"*40 + "\n"
	
	return peaks

def avdata(data,win):

	avd = ma.smoothListGaussian(data,win)
	
	plt.figure()
	plt.plot(data,'r')
	plt.hold()
	plt.plot(avd,'b')
	plt.show()

def main():
	
	plt.close('all')
		
	# Constructed test data
#	t = np.linspace(0,1e3,1e4)
#	x = np.sin(2*pi*t)*np.exp(-0.01*t)
#	x = np.sin(2*pi*0.23*t)*np.exp(-2*pi*0.23/100*t) + np.sin(2*pi*0.26*t)*np.exp(-2*pi*0.26/100*t)
	
	
	
#	fid = open('data.txt')
#	fid = open('L3mod_f02516_FFcen_Ey_re.txt')
	fid = open('BHS_r030_M0_re.txt')
	fid = open('BHS_r034_longrun_M0_re.txt')
		
	x = []
	for line in fid.readlines():
		x.append(float(line))

	fid.close()

	x = x[1500:]
#	print len(x)
	t = np.linspace(0,len(x)/32.,len(x))
	
	plt.figure()
	plt.plot(t,x)
	plt.xlabel('time')
	plt.ylabel('amplitude')
	plt.title('Time series')
	plt.show()
	
	freq,fx = FFT(t,x)
	
	peaks = findpeaks(freq,abs(fx)**2,ran=[0.19,0.31])
	
	P = findMTQ(freq,abs(fx)**2,peaks=peaks,t=2*t[-1])
	
	newpeaks = np.array([])
	for n in range(2,len(P),3):
		if P[n] > (freq[1] - freq[0])/10.0:
			newpeaks = np.append(newpeaks,P[n-1])


	if len(peaks) > len(newpeaks):
		Q = findMTQ(freq,abs(fx)**2,peaks=newpeaks,t=2*t[-1])
	
			
if __name__ == '__main__':
	main()
