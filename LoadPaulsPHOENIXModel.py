# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 13:59:34 2018

@author: andrew
"""

import numpy as np
import matplotlib.pyplot as plt 
from astropy.modeling.models import Lorentz1D
from scipy.interpolate import interp1d

plt.close('all')

c_ms = 299792458

sigma_cgs = 5.6704e-5  # cgs units (erg cm^(−2) s^(−1) K^(-4).

Teff_K = 3820

h_erg_s = 6.626070040e-27 # in erg s (1 erg  = 1e-7 J)


#LambdaPeak_m = 0.7626315789473685e-6 
#
#FreqPeak_Hz = c_ms/LambdaPeak_m
#print 'Freq Preak Hz: %.3e'%(FreqPeak_Hz)
#
#FreqPeak_Hz2 = 5.879e10*3820
#
#l2 = c_ms/FreqPeak_Hz2
#
#
#
d = np.loadtxt('phoenix_teff_3820_k.dat', skiprows=1)
#

UnscaledF =  d[:,1]

w_cm = d[:,0]
w_m = d[:,0]*1e-2

nu_Hz = c_ms/w_m

nu_104nm = c_ms/104e-9

nu_uband1 = c_ms/(365e-9 - 2*66e-9)
nu_uband2 = c_ms/(365e-9 + 2*66e-9)



ionizenu_index = np.argmin(np.abs(nu_Hz - nu_104nm))  #all values from here and higher (to the right) can ioinize 

uband_index1 = np.argmin(np.abs(nu_Hz - nu_uband1))  #all values from here and higher (to the right) can ioinize 
uband_index2 = np.argmin(np.abs(nu_Hz - nu_uband2))  #all values from here and higher (to the right) can ioinize 



nu_Hz = np.flipud(nu_Hz)
UnscaledF = np.flipud(UnscaledF)

#plt.plot(nu_Hz,UnscaledF)

IntegralOverNu = np.trapz(UnscaledF,nu_Hz) 

#print '%.3e'%(IntegralOverNu)

ScaledF = sigma_cgs*(Teff_K**4)*UnscaledF/IntegralOverNu

ScaledFIonize = np.trapz(ScaledF[ionizenu_index:],nu_Hz[ionizenu_index:])

ScaledFuband = np.trapz(ScaledF[uband_index1:uband_index2],nu_Hz[uband_index1:uband_index2])



#plt.figure()
#plt.plot(nu_Hz,ScaledF)

#w0_m = 3933.66e-10
#nu0_Hz = c_ms/w0_m
#A = 1.47e8 # FWHM in angular frequency units  
#OscillatorStrength = 0.682

#w0_m = 3968.47e-10
#nu0_Hz = c_ms/w0_m
#A = 1.4e8 # FWHM in angular frequency units  
#OscillatorStrength = 0.33

#######################

### Just assume this transition does not occur 

#w0_m = 7291.47e-10
#nu0_Hz = c_ms/w0_m
#A = 1.3 # FWHM in angular frequency units  
#OscillatorStrength = ????
#
#w0_m = 7323.89-10
#nu0_Hz = c_ms/w0_m
#A = 1.3 # FWHM in angular frequency units  
#OscillatorStrength = ????

###################################


#w0_m = 8498.02e-10
#nu0_Hz = c_ms/w0_m
#A = 1.11e6 # FWHM in angular frequency units  
#OscillatorStrength = 1.2e-2


#w0_m = 8542.09e-10
#nu0_Hz = c_ms/w0_m
#A = 9.9e6 # FWHM in angular frequency units  
#OscillatorStrength = 7.2e-2

w0_m = 8662.14e-10
nu0_Hz = c_ms/w0_m
A = 1.06e7 # FWHM in angular frequency units  
OscillatorStrength = 5.97e-2

#############

FWHM_freq_Hz = A/(2*np.pi)

## Consider the line profile from -100 ... 100 FWHM/2

LineProfileFreq_HalfExtent_Hz = (FWHM_freq_Hz/2)*100

LineProfileFreq_Hz = np.linspace(nu0_Hz-LineProfileFreq_HalfExtent_Hz,nu0_Hz+LineProfileFreq_HalfExtent_Hz,1e4)

RstarOverD = 1/3.3

######## Check Pauls calculation for Na to find the value of that factor 

fNa = 0.641
ANa = 6.16e7

Factor = 1.1e-9*ANa/fNa

#Nas2 = Factor*fNa/ANa  ### By inferring the factor From Paul's calculation 

PeakAbsorptionCrossSection_cm2 = Factor*OscillatorStrength/A

#PeakAbsorptionCrossSection_m2 = PeakAbsorptionCrossSection_cm2*1e-4

#################################

LineProfileFunction = Lorentz1D(amplitude=PeakAbsorptionCrossSection_cm2,x_0 = nu0_Hz, fwhm = FWHM_freq_Hz)

LineProfile = LineProfileFunction(LineProfileFreq_Hz)

ModelSpecRefIndex = np.argmin(np.abs(nu0_Hz - nu_Hz))

RelevantModelScaledF = ScaledF[ModelSpecRefIndex-1 : ModelSpecRefIndex+2] ## +2 because up to but not including 
RelevantModelFreq_Hz = nu_Hz[ModelSpecRefIndex-1 : ModelSpecRefIndex+2]

RelevantModelScaledFInterpFunc = interp1d(RelevantModelFreq_Hz, RelevantModelScaledF)

InterpolatedRelevantModelScaledF = RelevantModelScaledFInterpFunc(LineProfileFreq_Hz)

##########################

#plt.figure()
plt.figure()
plt.plot(1e6*c_ms/RelevantModelFreq_Hz,RelevantModelScaledF*(RstarOverD**2)*(1/(h_erg_s*RelevantModelFreq_Hz)),'b')
plt.plot(1e6*c_ms/LineProfileFreq_Hz,InterpolatedRelevantModelScaledF*(RstarOverD**2)*(1/(h_erg_s*LineProfileFreq_Hz)),'ro')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

plt.xlabel('wavelength (um)')
plt.ylabel(r'photons/s/cm$^2$/Hz')
plt.title('photons per unit frequency at location\n of K2-22 b from Teff = 3830 K model spectrum\n at wavelength of absorption line')


plt.figure()
plt.plot(nu0_Hz-LineProfileFreq_Hz,LineProfile,'.-')
plt.xlabel('delta frequency (Hz)')
plt.ylabel(r'absorption cross-section (cm$^2$)')
plt.title('Absorption cross-section')

################################


Eph = h_erg_s*LineProfileFreq_Hz

ProductFunction = InterpolatedRelevantModelScaledF*(RstarOverD**2)*(1/Eph)*LineProfile

IntegralOfProductFunction = np.trapz(ProductFunction,LineProfileFreq_Hz)

plt.figure()
plt.plot(nu0_Hz-LineProfileFreq_Hz,ProductFunction)
plt.ylabel('photons/s/Hz')
plt.xlabel('delta frequency (Hz)')
plt.title('Function is integrated over frequency\n to get photons/s')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

print('number of potentially absorbed photons per second')
print('%.3e'%(IntegralOfProductFunction))

NumberOfPhotonsPer_s_per_cm2_perHz = InterpolatedRelevantModelScaledF*(RstarOverD**2)*(1/Eph)


Eph_all = h_erg_s*nu_Hz

GammaAllNu = ScaledF*(RstarOverD**2)*(1/Eph_all)

Allplotlambda = 1e6*c_ms/nu_Hz

plt.figure()
plt.semilogx(Allplotlambda,GammaAllNu)
plt.xlabel('wavelength (um)')
plt.ylabel(r'photons/s/cm$^2$/Hz')
plt.title('photons per unit frequency at location\n of K2-22 b from Teff = 3830 K model spectrum')
plt.xlim((1e-1,1e2))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot([0.299,0.299],[0,4e6],'r--')
plt.plot([0.431,0.431],[0,4e6],'r--')

plt.plot([0.104,0.104],[0,4e6],'k--')

plt.figure()
plt.semilogx(nu_Hz,GammaAllNu)
plt.xlabel('wavelength (um)')
plt.ylabel(r'photons/s/cm$^2$/Hz')
plt.title('photons per unit frequency at location\n of K2-22 b from Teff = 3830 K model spectrum')
#plt.xlim((1e-1,1e2))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot([nu_uband1,nu_uband1],[0,4e6],'r--')
plt.plot([nu_uband2,nu_uband2],[0,4e6],'r--')

plt.plot([nu_104nm,nu_104nm],[0,4e6],'k--')



#d_wA = d[:,0]*1e8
#d_wum = d_wm*1e6

#Fl = d[:,1]*c_ms/d_wm**2

#Fl *= d_wm

#

#Teff = 3800
#sigma = 5.670367e-8 #W⋅m**(−2)⋅K**(−4)

#ScaledFl = sigma*(Teff**4)*Fl/IntegralOverWavelength

#
#d_wA = d[:,0]*1e8
#
#d_wum = d_wm*1e6
#
#
#
#f_Hz = c_ms/d_wm
#
#plt.plot(f_Hz,d[:,1])
#plt.plot([FreqPeak_Hz2,FreqPeak_Hz2],[0,np.max(d[:,1])],'r')
#
#plt.figure()
#plt.plot(d_wA,ScaledFl)

