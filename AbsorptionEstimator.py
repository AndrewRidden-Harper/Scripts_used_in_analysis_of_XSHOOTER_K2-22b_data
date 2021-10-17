# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 23:51:54 2018

@author: andrew
"""

import numpy as np
import matplotlib.pyplot as plt 
#from PyAstronomy import pyasl
from astropy.modeling.models import Lorentz1D
from scipy import interpolate
import broadChangedByAndrew 
import time

def EstimateNumberOfAtoms(DustMassPerOrbit_g,WeightFraction,AtomicWeight):
    print('total mass lost per orbit: %.3e g'%(DustMassPerOrbit_g))    
    
    MassOfElement = DustMassPerOrbit_g*WeightFraction/100.0
    print('Mass of element: %.3e g'%(MassOfElement))

    
    NumberOfMoles = MassOfElement/AtomicWeight
    print('number of moles of element: %.3e'%(NumberOfMoles))

        
    NumberOfAtoms = NumberOfMoles*NA
    
    print('number of atoms of element: %.3e'%(NumberOfAtoms))
    print()
    
    return NumberOfAtoms
    
    
def EstimateColumnDensity(NumberOfAtoms,l_tail_m,rp_m):
    
    TailArea_cm2 = l_tail_m*100.0*rp_m*100.0 ## in cm**2
    print('tail area: %.3e (cm**2)'%(TailArea_cm2))
    print('Tail area (star area): %.3e'%(TailArea_cm2/(np.pi*(Rstar*100)**2)))
    
    NatomsPerCm2 = NumberOfAtoms/TailArea_cm2
    print('Column density: %.3e (number of atoms per cm^2)'%(NatomsPerCm2))
    print()
    
    
    return NatomsPerCm2
    
def CalculateAbsorptionInCore(f,wavelength,NatomsPerCm2):
    
    '''
    Where f is the oscillator strength from NIST 
    '''
    
    tau = 2.654e-15*f*wavelength*NatomsPerCm2*FractionInCentralVelocityBin
    print('optical depth:', tau)
 
    T = np.exp(-tau)    
    print('transmission', T)
    
    A = 1 - T
    print('absorption', A)
    
    return tau
    
def EstimateAbsorption(DustMassPerOrbit_g,WeightFraction,AtomicWeight,l_tail_m,rp_m,f,wavelength):
    
    NumberOfAtoms = EstimateNumberOfAtoms(DustMassPerOrbit_g,WeightFraction,AtomicWeight)
    
    NatomsPerCm2 = EstimateColumnDensity(NumberOfAtoms,l_tail_m,rp_m)
    
    tau = CalculateAbsorptionInCore(f,wavelength,NatomsPerCm2)
    
    return tau
    
    

def EstimateOpticalDepth(FractionOfMassAsGas,timescale_s,TotalMasslossRate,TotalAbsorbingArea_cm2,AtomicWeight,OscillatorStrength,wavelength,FWHM_kms,NormalisationFactor):
    
    print('total mass lost per orbit: %.10e g'%(TotalMasslossRate*P_sec))    
    
    MassOfElement = TotalMasslossRate*FractionOfMassAsGas*timescale_s
    print('Mass of element: %.10e g'%(MassOfElement))
    
    NumberOfMoles = MassOfElement/AtomicWeight
    print('number of moles of element: %.10e'%(NumberOfMoles))
        
    NumberOfAtoms = NumberOfMoles*NA    
    print('number of atoms of element: %.10e'%(NumberOfAtoms))
    print()
    print('Total absorbing area: %.10e'%(TotalAbsorbingArea_cm2))
   
    NatomsPerCm2 = NumberOfAtoms/TotalAbsorbingArea_cm2
    print('Column density: %.10e (number of atoms per cm^2)'%(NatomsPerCm2))
    print()
    
    NPointsInProfile = 3e5

    VelocityVect_kms = np.linspace(-1,1,NPointsInProfile)
    
    WavelengthVect = wavelength*((VelocityVect_kms/c_kms) + 1.0)
    #WavelengthVect = wavelength*(VelocityVect_kms/c_kms) + wavelength
    #RegularWavelengthVect = np.linspace(np.min(WavelengthVect),np.max(WavelengthVect),NPointsInProfile)
    
#    FWHM_ms = 5.769206025
#    
#    NormalisationFactor = 1.0/9.0456058096719438  ### Found by indepently integrating this curve with amplitude = 1 
    
    LineProfileFunction = Lorentz1D(amplitude=NormalisationFactor,x_0 = 0, fwhm = FWHM_kms)
    
    LineProfile = LineProfileFunction(VelocityVect_kms)        
    
    tau = 2.654e-15*OscillatorStrength*wavelength*NatomsPerCm2*LineProfile
    

    
#    tauFunctionForInterp = interpolate.interp1d(WavelengthVect,tau,bounds_error=False)  
#        
#    tauInterpolatedToRegularWavelengthScale = tauFunctionForInterp(RegularWavelengthVect)
    
    
    ConvolvedTau, fwhm = broadChangedByAndrew.instrBroadGaussFast(WavelengthVect, tau, 11400,
          edgeHandling="firstlast", fullout=True)   
    
    TransmittedFlux = np.exp(-ConvolvedTau)
    
    FluxAbsorption = (1 - TransmittedFlux)*100  ## as percent 
    

    
    #return VelocityVect_ms,ConvolvedTau##tau #FluxAbsorptioProfileNormFactor_vectn,ConvolvedFlux
    return np.max(FluxAbsorption)
    #return VelocityVect_kms, tau
    

def AbsorptionFuncMass(MassOfElement,TotalAbsorbingArea_cm2,AtomicWeight,OscillatorStrength,wavelength,FWHM_kms,NormalisationFactor):
    
    print('Mass of element: %.10e g'%(MassOfElement))
    
    NumberOfMoles = MassOfElement/AtomicWeight
    print('number of moles of element: %.10e'%(NumberOfMoles))
        
    NumberOfAtoms = NumberOfMoles*NA    
    print('number of atoms of element: %.10e'%(NumberOfAtoms))
    print()
    print('Total absorbing area: %.10e'%(TotalAbsorbingArea_cm2))
   
    NatomsPerCm2 = NumberOfAtoms/TotalAbsorbingArea_cm2
    print('Column density: %.10e (number of atoms per cm^2)'%(NatomsPerCm2))
    print()
    
    NPointsInProfile = 3e5

    VelocityVect_kms = np.linspace(-1,1,NPointsInProfile)
    
    WavelengthVect = wavelength*((VelocityVect_kms/c_kms) + 1.0)
    
    LineProfileFunction = Lorentz1D(amplitude=NormalisationFactor,x_0 = 0, fwhm = FWHM_kms)
    
    LineProfile = LineProfileFunction(VelocityVect_kms)        
    
    tau = 2.654e-15*OscillatorStrength*wavelength*NatomsPerCm2*LineProfile
    

    
    ###  we convolved the optical depth profile with the instrument resolution 
    ### to broaden the line. This maximised the broadening while minimizing 
    ### optical thickness effects.      
    
    ConvolvedTau, fwhm = broadChangedByAndrew.instrBroadGaussFast(WavelengthVect, tau, 11400,
          edgeHandling="firstlast", fullout=True)   
    
    TransmittedFlux = np.exp(-ConvolvedTau)
    
    FluxAbsorption = (1 - TransmittedFlux)*100  ## as percent 
    

    
    return np.max(FluxAbsorption)    
    

    
def EstimateOpticalDepthConvoleIntensity(FractionOfMassAsGas,timescale_s,TotalMasslossRate,TotalAbsorbingArea_cm2,AtomicWeight,OscillatorStrength,wavelength,FWHM_kms,NormalisationFactor):
    
    print('total mass lost per orbit: %.10e g'%(TotalMasslossRate*P_sec))    
    
    MassOfElement = TotalMasslossRate*FractionOfMassAsGas*timescale_s
    print('Mass of element: %.10e g'%(MassOfElement))
    
    NumberOfMoles = MassOfElement/AtomicWeight
    print('number of moles of element: %.10e'%(NumberOfMoles))
        
    NumberOfAtoms = NumberOfMoles*NA    
    print('number of atoms of element: %.10e'%(NumberOfAtoms))
    print()
    print('Total absorbing area: %.10e'%(TotalAbsorbingArea_cm2))
   
    NatomsPerCm2 = NumberOfAtoms/TotalAbsorbingArea_cm2
    print('Column density: %.10e (number of atoms per cm^2)'%(NatomsPerCm2))
    print()
    
    NPointsInProfile = 3e5

    VelocityVect_kms = np.linspace(-1,1,NPointsInProfile)
    
    WavelengthVect = wavelength*((VelocityVect_kms/c_kms) + 1.0)
    
    LineProfileFunction = Lorentz1D(amplitude=NormalisationFactor,x_0 = 0, fwhm = FWHM_kms)
    
    LineProfile = LineProfileFunction(VelocityVect_kms)        
    
    tau = 2.654e-15*OscillatorStrength*wavelength*NatomsPerCm2*LineProfile
    
    
    TransmittedFlux = np.exp(-tau)
    
    FluxAbsorption = (1 - TransmittedFlux)*100  ## as percent 
    

    ConvolvedAbsorption, fwhm = broadChangedByAndrew.instrBroadGaussFast(WavelengthVect, FluxAbsorption, 11400,
          edgeHandling="firstlast", fullout=True)   
    

    return np.max(ConvolvedAbsorption)

    
def EstimateOpticalDepthTesting(FractionOfMassAsGas,timescale_s,TotalMasslossRate,TotalAbsorbingArea_cm2,AtomicWeight,OscillatorStrength,wavelength,FWHM_kms,NormalisationFactor):
    
    print('total mass lost per orbit: %.10e g'%(TotalMasslossRate*P_sec))    
    
    MassOfElement = TotalMasslossRate*FractionOfMassAsGas*timescale_s
    print('Mass of element: %.10e g'%(MassOfElement))
    
    NumberOfMoles = MassOfElement/AtomicWeight
    print('number of moles of element: %.10e'%(NumberOfMoles))
        
    NumberOfAtoms = NumberOfMoles*NA    
    print('number of atoms of element: %.10e'%(NumberOfAtoms))
    print()
    print('Total absorbing area: %.10e'%(TotalAbsorbingArea_cm2))
   
    NatomsPerCm2 = NumberOfAtoms/TotalAbsorbingArea_cm2
    print('Column density: %.10e (number of atoms per cm^2)'%(NatomsPerCm2))
    print()
    
    NPointsInProfile = 3e5

    #VelocityVect_kms = np.linspace(-1,1,NPointsInProfile)
    VelocityVect_kms = np.linspace(-50,50,NPointsInProfile)

    
    WavelengthVect = wavelength*((VelocityVect_kms/c_kms) + 1.0)

    
    LineProfileFunction = Lorentz1D(amplitude=NormalisationFactor,x_0 = 0, fwhm = FWHM_kms)
    
    LineProfile = LineProfileFunction(VelocityVect_kms)        
    
    tau = 2.654e-15*OscillatorStrength*wavelength*NatomsPerCm2*LineProfile
    
    return WavelengthVect,tau
    
def ConvertAbsorptionToTau(AbsorptionPercent):
    
    Afraction = AbsorptionPercent/100
    
    T = 1.0 - Afraction
    
    tau = -np.log(T)
    
    return tau
    
def CalculateColumnDensity(AbsorptionPercent,f,w,phi0):
    
    const = 2.654e-15
    
    tau = ConvertAbsorptionToTau(AbsorptionPercent)
    
    N = tau/(const*f*w*phi0)
    
    return N 
    
    
    
    


c_kms = 299792.458 

NA = 6.0221409e+23  #Avogadro's number 

Rearth = 6371e3 # in m
Rmars = 3390e3 # in m
Rmercury = 2440e3 # in m
Rmoon = 1737e3 # in m 
Rsun = 695508e3 # in m

Rstar = 0.57*Rsun

Rstar_cm = Rstar*100

Astar_cm2 = np.pi*Rstar_cm**2

Atail_cm2 = 0.013*Astar_cm2


Lstar = 0.063  ## in solar luminosity 
Orbital_distance = 0.0088 ## in au 


wavelengths = np.array([5889.95, 5895.92, 8498.02, 8542.09, 8662.14])
oscillatorStrengths = np.array([6.41e-01, 3.20e-01, 1.20e-02, 7.2e-02, 5.97e-02])    # f 
#oscillatorStrengths = 10**np.array([0.108, -0.194, -1.318, -0.36, -0.622])    # f 



NaAtomicWeight_vect = np.array([22.90, 22.90, 40.078, 40.078, 40.078])
OLDtimescale_vect = np.array([20.0,20.0,1e3,1e3,1e3])
timescale_vect = np.array([7036.45491004,7036.45491004,1.25651e+05,1.25651e+05,1.25651e+05])
OldToNewTimeScalingFactorVect =  OLDtimescale_vect/timescale_vect


FWHM_kms_vect = np.array([0.005774474287514911,
                          0.005761559946136618,
                          0.00014363310437195365,
                          0.00134592068935749,
                          0.0014613397426792723])


## More accurate normalisation factors 
ProfileNormFactor_vect = np.array([112.02986546200165,
                                       112.2809766090983,
                                       4309.0884333819349,
                                       480.64762111124384, 
                                       442.68526931184999])
NaAbundance = 2.5e-2
CaAbundance = 2.96e-2

### For Na 
i = 3
#mass = 0.9991673605328892*0.9730238857903485*2.3484449720e+09*1.2 ##g
#mass = 1.7382671480144405e+09*1.0006671114076051
#TerrestrailFrac = NaAbundance 

### For Ca 
#i = 3
#mass = 0.9995456610631531*0.9972801450589303*8.8320590504e+09 ##g
#mass = 1.053e9
mass = 5.587e+09*0.9992862241256245
TerrestrailFrac = CaAbundance 

A = AbsorptionFuncMass(mass,Astar_cm2,NaAtomicWeight_vect[i],oscillatorStrengths[i],wavelengths[i],FWHM_kms_vect[i],ProfileNormFactor_vect[i])   
print('Mass %.5e gives %.3f %% absorption in line %.2f'%(mass,A,wavelengths[i]))

ImpliedTotalM = mass/TerrestrailFrac
print('assuming terrestrial abundances, gives a total mass of %.5e g'%(ImpliedTotalM))



Mdot = 2e11 # g/s 

P_hr = 9.145872

P_sec = P_hr*3600.0

MassLostPerOrbit_g = Mdot*P_sec

MassLostPerOrbit_kg = MassLostPerOrbit_g/1e3



