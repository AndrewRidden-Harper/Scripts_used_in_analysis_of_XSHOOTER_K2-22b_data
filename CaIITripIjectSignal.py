# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:30:26 2015

@author: arh

Think need an interval of 800 pixels either side 

"""

from CaII_inject_model_spectrum import call_model_injector
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt 
import scipy
from DopplerShiftToPlanetFrame import shift_to_planet_frame, findNANedges, whole_pixel_shift_to_planet_frame
from K2_22bJDToPhaseFuncCustomTransLims import JDtophase_radial_vel
import os
from pcasub import pcasub 
from CaII_VALD_MakeSpec import makespec
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from astropy.stats import sigma_clip, mad_std
import pandas as pd
from astropy.modeling import models, fitting

plt.close('all')

vsys = 27.58  # From SIMBAD 

c = 299792.458

def binspec(spec,pointsinbin=5):
    
    '''
    A way to bin every N points in a spectrum based on 
    doing a clever reshape then average along a row
    
    It returns a spectrum with representative bin sizes in a simple way.    
    
    '''   
        
    numcols = len(spec)
    
    remainderPoints = numcols%pointsinbin 
    
    extraRequiredPoints = pointsinbin-remainderPoints
    
    NANfiller_array = np.empty((numcols+extraRequiredPoints))
    NANfiller_array[:] = np.NAN
    
    NANfiller_array[:numcols] = spec
    
    binnedspec = np.nanmean(NANfiller_array.reshape(-1, pointsinbin), axis=1)
    
    full_length_binned_spec = np.empty_like(spec)
    
    for i in range(len(binnedspec)):
        
        full_length_binned_spec[i*pointsinbin:(i+1)*pointsinbin] = binnedspec[i]
        
    return full_length_binned_spec
    
def BinnnerCentrePoint(VectToBin,PointsInBin,CentrePoint):
    
    '''
    centre point should be an index (starting at 0)
    
    This is a function that depends on binspec which adds the extra 
    functionality to allow the point that the bins should be centred around.
    It places a central bin around this central point, then uses binspec 
    on the left and right parts of the vector around the centred point    
    '''
    
    binned_vect = np.empty_like(VectToBin)
    binned_vect[:] = np.NAN
    
    halfbin = np.floor(PointsInBin/2.0)
    
    if PointsInBin % 2.0 != 0: #is an odd number 
    
        print('odd number of points in bin')
    
        binned_vect[CentrePoint-halfbin:CentrePoint+halfbin+1] =  np.mean(VectToBin[CentrePoint-halfbin:CentrePoint+halfbin+1])    
        
        leftpart = VectToBin[:CentrePoint-halfbin]
        rightpart = VectToBin[CentrePoint+halfbin+1:]
        
        flippedleftpart = leftpart[::-1]    
        flippedbinned_leftpart = binspec(flippedleftpart,PointsInBin)    
        binned_leftpart = flippedbinned_leftpart[::-1]
        
        binned_vect[:CentrePoint-halfbin] = binned_leftpart    
           
        ## put the binned right part in the binned vector    
        binned_vect[CentrePoint+halfbin+1:] = binspec(rightpart,PointsInBin)
    
    if PointsInBin % 2.0 == 0: #is an even number 
    
        print('even number of points in bin')
    
    
        binned_vect[CentrePoint-halfbin:CentrePoint+halfbin] =  np.mean(VectToBin[CentrePoint-halfbin:CentrePoint+halfbin])
            
        leftpart = VectToBin[:CentrePoint-halfbin]
        rightpart = VectToBin[CentrePoint+halfbin:]
    
        flippedleftpart = leftpart[::-1]
        
        flippedbinned_leftpart = binspec(flippedleftpart,PointsInBin)
        
        binned_leftpart = flippedbinned_leftpart[::-1]
        
        binned_vect[:CentrePoint-halfbin] = binned_leftpart    
        
        binned_vect[CentrePoint+halfbin:] = binspec(rightpart,PointsInBin)
        
    return binned_vect
    

def makeHistogram(data,numbins=50):

    # the histogram of the data with histtype='step'
    n, bins, patches = P.hist(data, numbins, normed=1, histtype='bar') #histtype='stepfilled'
    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)    
    P.xlabel('noise value bin')
    P.ylabel('occurance count')
    P.savefig('../../papers/55Cnce_Na_paper/plots/combined_no_injected_noise_hist.png')
    
    
    # add a line showing the expected distribution
    #y = P.normpdf( bins, mu, sigma)
    #l = P.plot(bins, y, 'k--', linewidth=1.5)
    
    


def calculate_SNR(spectrum,minindex=402,signal_extent=22,noise_extent=230,toplot=False):
#This extent used to be 10 based on the narrow injected signal but the binned 
#combined signal indicates a wider signal 
    
    recoveredsignal = -np.min(spectrum[minindex-signal_extent:minindex+signal_extent])
    signal = InjectionStrength#*(recoveredsignal/0.048271636788569633)   ## option to scale relative to the 0 pca minimum
    
    leftnoise = np.std(spectrum[minindex-noise_extent:minindex-signal_extent])
    
    rightnoise = np.std(spectrum[minindex+signal_extent:minindex+noise_extent])
    
    avgnoise = (leftnoise+rightnoise)/2.0
    
    if toplot == True:
        
        xvect = np.arange(len(spectrum))
        plt.figure()
        plt.plot(xvect[minindex-noise_extent:minindex+noise_extent],spectrum[minindex-noise_extent:minindex+noise_extent],'b')
        plt.plot(xvect[minindex-signal_extent:minindex+signal_extent],spectrum[minindex-signal_extent:minindex+signal_extent],'r')
        
                
    #return signal/avgnoise 
    
    return [signal/avgnoise,signal,avgnoise]
    
def calculate_SNR_BetterPlot(spectrum,minindex=402,signal_extent=22,noise_extent=230,toplot=False):
#This extent used to be 10 based on the narrow injected signal but the binned 
#combined signal indicates a wider signal 
    
    recoveredsignal = -np.min(spectrum[minindex-signal_extent:minindex+signal_extent])    
    signal = InjectionStrength#*(recoveredsignal/0.048271636788569633)  ## option to scale relative to the 0 pca minimum 
    
    leftnoise = np.std(spectrum[minindex-noise_extent:minindex-signal_extent])
    
    rightnoise = np.std(spectrum[minindex+signal_extent:minindex+noise_extent])
    
    avgnoise = (leftnoise+rightnoise)/2.0
    
    
    if toplot == True:

        xvect = plot_vx
        plt.figure()
        plt.title('min of red = signal, std dev of blue = noise\n signal = %f, noise = %f, S/n = %f'%(signal,avgnoise,signal/avgnoise))
        plt.plot(xvect[minindex-noise_extent:minindex+noise_extent],spectrum[minindex-noise_extent:minindex+noise_extent],'b')
        plt.plot(xvect[minindex-signal_extent:minindex+signal_extent],spectrum[minindex-signal_extent:minindex+signal_extent],'r') 
        plt.xlabel('radial velocity (km/s)')           
        plt.ylabel('absorption relative to stellar spectrum')
        plt.ylim(PlotYLimsCombined)
        
        plt.figure()
        plt.title('part shown in plot')
        plt.plot(xvect[minindex-noise_extent:minindex+noise_extent],spectrum[minindex-noise_extent:minindex+noise_extent],'b')
        plt.plot(xvect[minindex-signal_extent:minindex+signal_extent],spectrum[minindex-signal_extent:minindex+signal_extent],'r')   
        plt.xlim(PlotXLims)
        plt.ylim(PlotYLimsCombined)
        plt.xlabel('radial velocity (km/s)')           
        plt.ylabel('absorption relative to stellar spectrum')
                
    #return signal/avgnoise 
    
    return [signal/avgnoise,signal,avgnoise]
    
    
def MakeVariableKpMatrix(name,fitsfile,bjd,componets_to_remove,width,vorbital_range,points_in_bin):
    
    ResultsArray = np.empty(((len(vorbital_range)),800))
    ResultsArray[:] = np.NAN
    
    BinnnedResultsArray = np.empty(((len(vorbital_range)),800))
    BinnnedResultsArray[:] = np.NAN
    
    for i in range(len(vorbital_range)):
        
        Kp = vorbital_range[i]
    
        ResultsArray[i,:] = process(name,fitsfile,bjd,componets_to_remove,width,vorbital=Kp)
        BinnnedResultsArray[i,:] = binspec(ResultsArray[i,:],points_in_bin)
        
    pyfits.writeto('VariableKp.fits',ResultsArray,clobber=True)
    pyfits.writeto('BinnedVariableKp.fits',BinnnedResultsArray,clobber=True)
    
    return ResultsArray, BinnnedResultsArray
        
        
def process(name,data,wavelength,bjd,componets_to_remove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=5E-5,ModelRadialVelocity = 0):


## New line positions
#### These are from visual inspection but form the wavelenght scale it should be this +1 \
#### By looking at the 2D spectra, it almost needs to be +2 
    if name == 'night1':       
        CaTrip1PixCenter = 1000 
        CaTrip2PixCenter = 1221 
        CaTrip3PixCenter = 1821
        
        I1cols = [993,1213,1813]        
        I2cols = [998,1218,1819]        
        I3cols = [1003,1223,1824]        
        I4cols = [1010,1230,1831]    
        
    if name == 'night2':       
        CaTrip1PixCenter = 1001
        CaTrip2PixCenter = 1222 
        CaTrip3PixCenter = 1822
        
        I1cols = [994,1214,1814]        
        I2cols = [999,1219,1819]        
        I3cols = [1004,1224,1825]        
        I4cols = [1010,1231,1831]       
        
    if name == 'night3':       
        CaTrip1PixCenter = 999 
        CaTrip2PixCenter = 1220 
        CaTrip3PixCenter = 1820
        
        I1cols = [992,1212,1812]        
        I2cols = [997,1217,1818]        
        I3cols = [1002,1223,1823]        
        I4cols = [1007,1228,1828]  

    if name == 'night4':       
        CaTrip1PixCenter = 1000 
        CaTrip2PixCenter = 1220 
        CaTrip3PixCenter = 1821  
        
        I1cols = [992,1212,1813]        
        I2cols = [997,1217,1818]        
        I3cols = [1004,1224,1825]        
        I4cols = [1009,1229,1830]  
    
    
    LinePartHalfExtent = 900
    
    numrows,numcols = np.shape(data)        
    
    if not os.path.exists('processing/%s/pca%d'%(name,pcatoremove)):
        os.makedirs('processing/%s/pca%d'%(name,pcatoremove))
    
    if ToInject == True:
        toprocess = 'injected'
    
    if ToInject == False:
        toprocess = 'no_injected'
    
    species = 'ToyCaIITriplet'
    
    flux = injectStrength
   
       
    #Create the model template if it does not exist already 
    #if not os.path.exists('model_spectra/%s/ToyCaIITriplet_F%.2E_s%f' % (species,flux,width)):
    makespec(species='%s'%(species),flux=flux,sigma=width)    
    
    d = data   
    w = wavelength
    
   
    if toprocess == 'injected':
        #d = call_model_injector(d,w,jd,species,flux,width,vshift=vsys+8)  # This vshift includes the system velocity and and shift to seperate 
        d, injected_spec = call_model_injector(d,w,jd,species,flux,width,SpecialTransLims,vshift=ModelRadialVelocity)        
        #d = injected_spec
        
        writedir = 'CaIItriplet/%s/pca%d/Injected_F%.2E' % (name,pcatoremove,flux)     
        
        if not os.path.exists(writedir):
            os.makedirs(writedir)
            
    if toprocess == 'no_injected':
        
        writedir = 'CaIItriplet/%s/pca%d' % (name,pcatoremove)
        
        if not os.path.exists(writedir):
            os.makedirs(writedir)
            
    if ToInject == True:   
        pyfits.writeto('%s/InjectedSpectrum.fits'%(writedir),injected_spec,clobber=True)    
        
    ##### Divide by the 5-sigma clipped mean spectrum 
    d2 = np.copy(d)

    for i in range(numcols):
        
        col = d2[:,i]
        
        mask = np.zeros(len(col))
        mask[SpecialTransLims[0]:SpecialTransLims[1]+1] = 1  ### 1 is a masked value 
        
        
        MaskedArray = np.ma.masked_array(col,mask)        
                
        filtered_data = sigma_clip(MaskedArray, sigma=5, iters=1, stdfunc=mad_std)
        #filtered_data = sigma_clip(col, sigma=5, iters=1, stdfunc=mad_std)
        
        
        #d2[:,i] = d[:,i]/np.mean(filtered_data)
        #d2[:,i] = d[:,i]/np.mean(MaskedArray)
        d2[:,i] = d[:,i]/np.mean(col)       
    
    d2[np.isnan(d2)] = 1.0
    d2[np.isinf(d2)] = 1.0
    
    pyfits.writeto('%s/d2.fits'%(writedir),d2,clobber=True)
    print('Writedir: %s'%(writedir))
    
    d1b = np.copy(d2)    
   
    d2_planetframe = shift_to_planet_frame(d2,w,jd,vorbital='calculate')
    pyfits.writeto('%s/d2PlanetFrame.fits'%(writedir),d2_planetframe,clobber=True)    
    
    d2b = d1b    
     
    pyfits.writeto('%s/d2b.fits'%(writedir),d2b,clobber=True)
    
    d3 = np.copy(d2b)    
   
    pyfits.writeto('%s/d3_prepcasub.fits'%(writedir),d3,clobber=True)
    
    
    #### Optionally use PCA         
    if componets_to_remove>0:
       
        d3 = pcasub(d2b,componets_to_remove)
        
    else: 

        d3 = d3 - np.mean(d3)
    
    pyfits.writeto('%s/d3_pcasub.fits'%(writedir),d3,clobber=True)
        
    d3preSNRweight = np.copy(d3)
    
    def GiveAverageSNROfACol(ColNumber):
        
        return 1/np.mean([np.std(d3[:,ColNumber-1]),np.std(d3[:,ColNumber]),np.std(d3[:,ColNumber+1])])        
        
    
    
    ### Variables are called noise but they are actually S/N    
    I1anoise = GiveAverageSNROfACol(I1cols[0])    
    I1bnoise = GiveAverageSNROfACol(I1cols[1])
    I1cnoise = GiveAverageSNROfACol(I1cols[2])
    
    I2anoise = GiveAverageSNROfACol(I2cols[0])
    I2bnoise = GiveAverageSNROfACol(I2cols[1])
    I2cnoise = GiveAverageSNROfACol(I2cols[2])
    
    I3anoise = GiveAverageSNROfACol(I3cols[0])
    I3bnoise = GiveAverageSNROfACol(I3cols[1])
    I3cnoise = GiveAverageSNROfACol(I3cols[2])
    
    I4anoise = GiveAverageSNROfACol(I4cols[0])
    I4bnoise = GiveAverageSNROfACol(I4cols[1])
    I4cnoise = GiveAverageSNROfACol(I4cols[2])
    
    I1meannoise = np.mean([I1anoise,I1bnoise,I1cnoise])
    I2meannoise = np.mean([I2anoise,I2bnoise,I2cnoise])
    I3meannoise = np.mean([I3anoise,I3bnoise,I3cnoise])
    I4meannoise = np.mean([I4anoise,I4bnoise,I4cnoise])
    
    pyfits.writeto('%s/pre_col_std_scaled_d3.fits'%(writedir),d3,clobber=True)
    
    residual_col_std = np.std(d3,axis=0)
    
    np.savetxt('%s/CaSNR.txt'%(writedir),1/residual_col_std)
    
    scaled_resid_col_std = residual_col_std*(1.0/np.median(residual_col_std))
    
    residual_matrix_mean = np.mean(d3)
    

    
    ########  DON'T DO: Weight columns according to noise in column 
#    for i in range(numcols):
#        
#        if scaled_resid_col_std[i] != 0.0:
#        
#            d3[:,i] = d3[:,i]/scaled_resid_col_std[i]
        
    pyfits.writeto('%s/d3_pcasub_colweighted.fits'%(writedir),d3,clobber=True)
    
    #### Combine the Na D lines     
    combined_interval = 500
    
    pyfits.writeto('%s/CaTripLine1Part.fits'%(writedir),d3[:,CaTrip1PixCenter-LinePartHalfExtent:CaTrip1PixCenter+LinePartHalfExtent],clobber=True)
    pyfits.writeto('%s/CaTripLine2Part.fits'%(writedir),d3[:,CaTrip2PixCenter-LinePartHalfExtent:CaTrip2PixCenter+LinePartHalfExtent],clobber=True)
    pyfits.writeto('%s/CaTripLine3Part.fits'%(writedir),d3[:,CaTrip3PixCenter-LinePartHalfExtent:CaTrip3PixCenter+LinePartHalfExtent],clobber=True)
    
    
    CaTrip1 = d3[:,CaTrip1PixCenter-LinePartHalfExtent:CaTrip1PixCenter+LinePartHalfExtent]
    CaTrip2 = d3[:,CaTrip2PixCenter-LinePartHalfExtent:CaTrip2PixCenter+LinePartHalfExtent]
    CaTrip3 = d3[:,CaTrip3PixCenter-LinePartHalfExtent:CaTrip3PixCenter+LinePartHalfExtent]

     
    
    #### Weights proportional to oscillator strengths (like Paul Molliere did)
    Weight1 = 0.16666667**2
    Weight2 = 1.0**2
    Weight3 = 0.82916667**2    
    
#    Weight1 = 1.0
#    Weight2 = 1.0
#    Weight3 = 1.0
    
    
    combined = (Weight1*CaTrip1 + Weight2*CaTrip2 + Weight3*CaTrip3)/(Weight1 + Weight2 + Weight3)
    
    
    pyfits.writeto('%s/CaIItrip_combined.fits'%(writedir),combined,clobber=True)
    
    wcombined = w[CaTrip2PixCenter-LinePartHalfExtent:CaTrip2PixCenter+LinePartHalfExtent]
    
    np.save('%s/%s_combined.npy'%(writedir,name),wcombined)
    
    d4 = shift_to_planet_frame(combined,wcombined,jd,vorbital='calculate')
    
    combined_numcols = np.shape(combined)[1]
    
    combined_Na_centre = (combined_numcols/2.0)-1 # -1 to get a python index 
    
    pyfits.writeto('%s/not_cropped_edge_nans.fits' % (writedir),d4,clobber=True)
    
    d4 = d4[:,combined_Na_centre-combined_interval:combined_Na_centre+combined_interval]
    
    pyfits.writeto('%s/d4.fits'%(writedir),d4,clobber=True)
    
    
    phase,radv,translims,numorbits = JDtophase_radial_vel(jd,SpecialTransLims)
    
    
    ### Weights for the weighted average planet transmission spectrum by previously reading off the columns 
    starttrans = translims[0]
    endtrans = translims[1]
    
   
    SumOfWeights = I1meannoise + I2meannoise + I3meannoise + I4meannoise
    
    InjectedSpectrumWeights = np.ones_like(d4[starttrans:endtrans+1,:])
    InjectedSpectrumWeights[0,:] *= I1meannoise
    InjectedSpectrumWeights[1,:] *= I2meannoise
    InjectedSpectrumWeights[2,:] *= I3meannoise
    InjectedSpectrumWeights[3,:] *= I4meannoise    
    
    np.savetxt('SummingWeights.txt',InjectedSpectrumWeights)   
    
    #planetspec = np.sum(d4[starttrans:endtrans+1,:]*InjectedSpectrumWeights,axis=0)/SumOfWeights  ## To include the last spectrum in SpecialTransitLimits 
    
    ####planetspec = np.sum(d4[starttrans:endtrans+1,:],axis=0)  ## old To include the last spectrum in SpecialTransitLimits 
    
    ### Unweighted alternative 
    planetspec = np.sum(d4[starttrans:endtrans+1,:],axis=0)/4  ## To include the last spectrum in SpecialTransitLimits    
    
    sinday = 86164.09054 #exact number of seconds in a sidereal day from WolframAlpha
    
    cadence = (jd[1:]-jd[0:-1])*sinday
    
    print('median cadence %f' % (np.median(cadence)))

    
    mbjd = jd-2400000.5
    
    print('start bjd %.5f'%mbjd[0])
    print('end bjd %.5f'%mbjd[-1])
    print('median cadence %.1f:'%np.median(cadence))
    
    numspecbeforetrans = translims[0]
    numspecduringtrans = translims[1]-translims[0] + 1 
    numspecaftertrans = (numrows)-translims[1]
    
    
    print('Nr. spec before transit (index)',numspecbeforetrans)
    print('Nr. spec during transit',numspecduringtrans)
    print('Nr. spec after transit',numspecaftertrans)
    print('start phase %.3f'% phase[0])
    print('end phase %.3f' % phase[-1])
    print('end trans lim', translims[1])
    
    print('total number of specta:', numrows)
    print('sum of the pre,during, post:',(numspecbeforetrans+numspecduringtrans+numspecaftertrans))
    
    print('wstep',wcombined[2]-wcombined[1])
    
  
    InjectedPartsSTD = np.array([np.mean(residual_col_std[987:1015]),np.mean(residual_col_std[1207:1235]),np.mean(residual_col_std[1808:1835])])


    return planetspec,phase,radv,translims,numorbits,[I1anoise,I1bnoise,I1cnoise,I2anoise, I2bnoise, I2cnoise, I3anoise, I3bnoise, I3cnoise, I4anoise, I4bnoise, I4cnoise], InjectedPartsSTD
    
    


#################################


notbinned_snrlist = []
binned_snrlist = []



No1InjectedSNR = []
No1BinnedInjectedSNR = []

pcavect = [0]  



#width = 0.01
#width = 0.34#*5
width = 0.34*1.2878236525174087
#width = 0.25
#vorbital = 240



ToInjectTrueFals = True  
#InjectionStrength = 4.0e-2
#InjectionStrength = 3e-2
#InjectionStrength = 3e-2
#InjectionStrength = 1.5e-2

#InjectionStrength = 1.7e-2  ### 5 sigma with wider signal 
#InjectionStrength = 1.1e-2
#InjectionStrength = 5.5e-2
#InjectionStrength = 1.1e-2

#InjectionStrength = 1.4E-2
#InjectionStrength = 2.15e-2 
#InjectionStrength = 4.4e-2 

#InjectionStrength = 0.56e-2 ### # !! 
#InjectionStrength = 1.1e-2

#InjectionStrength = 1e-2
#InjectionStrength = 1.5e-2

#InjectionStrength = 0.7e-2
#InjectionStrength = 0.8e-2

#InjectionStrength = 2e-2
#InjectionStrength = 1.6e-2
#InjectionStrength = 1.5e-2

#InjectionStrength = 0.488e-2

#InjectionStrength = 0.55627079774341615e-2
#InjectionStrength = 10e-2

#InjectionStrength = 1e-2
#### What it should be based on the combined noise 
InjectionStrength = 0.0067   
EffectiveInjectionStrength = 0.0067# 6e-3

#InjectionStrength = 1e-2
#EffectiveInjectionStrength = InjectionStrength

#InjectionStrength = 6e-3 
#EffectiveInjectionStrength = 6e-3


#InjectionStrength = 0.6e-2
#InjectionStrength = 6.6e-3

#InjectionStrength = 99e-2
#InjectionStrength = 10e-2

#InjectionStrength = 1e-2 
#EffectiveInjectionStrength = 1e-2

#PlotYLims = (-0.07,0.04)
PlotYLims = (-0.012,0.012)
#PlotYLims = (-0.15,0.15)
PlotYLimsCombined = PlotYLims#(-0.05,0.02)
#PlotXLims = (-500,500)
#PlotXLims = (-1300,1300)
PlotXLims = (-500,500)

EmpiricalWavelenghtShift = 0.0#-0.22
EmpiricalWavelenghtShift = 0.0#-0.5 ## required to get the planet frame transformed to have zero radial velocity 



#Weights = np.array([66.033222,70.443019,68.937256,66.469294])
#Weights = np.array([1.0,1.0,1.0,1.0])



#Weights = np.array([96.56360597, 107.67688262, 109.68346268,106.68958225])


#Weights = np.array([94.104425,108.875060,99.896728,98.926709])  ## From 1/std of a part without lines in the AlignedNoCosmics residuals 

#DataSetWeight = 1.0/(np.array([0.0063426417309907614,0.0031646285957089198,0.0037694087799296839,0.0031563813481578621]))
DataSetWeight = np.array([1.0, 1.0, 1.0, 1.0])

#DataSetWeight = 1.0/(np.array([0.0063426417309907614,0.0031646285957089198,0.0037694087799296839,0.0031563813481578621])**2)



BinCentralPixel = 501
BinSizePixels = 4

#SNRCalcParams  = (500,10,90)
SNRCalcParams  = (500,10,250)

#
for pcatoremove in pcavect:       
    

    name = 'night1'
    DopplerShiftForModel = -7.68873362
    SpecialTransLims = (5,8) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  
    
    fitsfile = pyfits.open('processing/night1/AlignedNoCosmics/night1MedNormNoCosmics.fits')

    
    data = fitsfile[0].data    

    full_mjd_vect = np.loadtxt('IgnasReduced/night1/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night1/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/night1/Night1WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n1spec,n1phase,n1radv,n1translims,n1numorbits,n1linenoise,n1InjectedPartsSTD = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    
    
    n1specNoInject,n1phaseNoInject,n1radvNoInject,n1translimsNoInject,n1numorbitsNoInject,n1linenoiseNoInject,n1InjectedPartsSTDNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    
#
#
#
    name = 'night2'
    DopplerShiftForModel = -1.48105377
    SpecialTransLims = (5,8) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  
    
    fitsfile = pyfits.open('processing/night2/AlignedNoCosmics/night2MedNormNoCosmics.fits')

    data = fitsfile[0].data

    full_mjd_vect = np.loadtxt('IgnasReduced/night2/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night2/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/night2/Night2WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n2spec,n2phase,n2radv,n2translims,n2numorbits,n2linenoise,n2InjectedPartsSTD = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    
    n2spec = np.roll(n2spec,1)    
    
    n2specNoInject,n2phaseNoInject,n2radvNoInject,n2translimsNoInject,n2numorbitsNoInject,n2linenoiseNoInject,n2InjectedPartsSTDNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    
    n2specNoInject = np.roll(n2specNoInject,1)    
#    
#
    name = 'night3'
    DopplerShiftForModel = -12.59669067 
    SpecialTransLims = (4,7) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  
    
    fitsfile = pyfits.open('processing/night3/AlignedNoCosmics/night3MedNormNoCosmics.fits')

    data = fitsfile[0].data
    
    full_mjd_vect = np.loadtxt('IgnasReduced/night3/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night3/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/night3/Night3WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n3spec,n3phase,n3radv,n3translims,n3numorbits,n3linenoise,n3InjectedPartsSTD = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    
    n3spec = np.roll(n3spec,-1)   
    
    n3specNoInject,n3phaseNoInject,n3radvNoInject,n3translimsNoInject,n3numorbitsNoInject,n3linenoiseNoInject,n3InjectedPartsSTDNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    
    n3specNoInject = np.roll(n3specNoInject,-1)   
    
    
    
    name = 'night4'

    DopplerShiftForModel = -8.80273704
    SpecialTransLims = (6,9) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  
    
    fitsfile = pyfits.open('processing/night4/AlignedNoCosmics/night4MedNormNoCosmics.fits')

    data = fitsfile[0].data
    
    full_mjd_vect = np.loadtxt('IgnasReduced/night4/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night4/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/night4/Night4WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n4spec,n4phase,n4radv,n4translims,n4numorbits,n4linenoise,n4InjectedPartsSTD = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    

    n4specNoInject,n4phaseNoInject,n4radvNoInject,n4translimsNoInject,n4numorbitsNoInject,n4linenoiseNoInject,n4InjectedPartsSTDNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    

    deltaw = wavelength[1] - wavelength[0]
    
    ##################
    ### Divide by rolling average 
    
    RollingMeanWindow = 20 
    
    n1NIRM = pd.rolling_mean(pd.DataFrame(n1specNoInject),window=RollingMeanWindow).values.flatten()
    n2NIRM = pd.rolling_mean(pd.DataFrame(n2specNoInject),window=RollingMeanWindow).values.flatten()
    n3NIRM = pd.rolling_mean(pd.DataFrame(n3specNoInject),window=RollingMeanWindow).values.flatten()
    n4NIRM = pd.rolling_mean(pd.DataFrame(n4specNoInject),window=RollingMeanWindow).values.flatten()

#    n1specNoInject = ((n1specNoInject + 1)/(n1NIRM + 1)) - 1 
#    n1spec = ((n1spec + 1)/(n1NIRM + 1)) - 1 
#    
#    n2specNoInject = ((n2specNoInject + 1)/(n2NIRM + 1)) - 1 
#    n2spec = ((n2spec + 1)/(n2NIRM + 1)) - 1 
#    
#    n3specNoInject = ((n3specNoInject + 1)/(n3NIRM + 1)) - 1 
#    n3spec = ((n3spec + 1)/(n3NIRM + 1)) - 1 
#    
#    n4specNoInject = ((n4specNoInject + 1)/(n4NIRM + 1)) - 1 
#    n4spec = ((n4spec + 1)/(n4NIRM + 1)) - 1 


    
    #### Now the binning 
    
    n1specBinned = BinnnerCentrePoint(n1spec,BinSizePixels,BinCentralPixel)    

    n1specBinnedNoInject = BinnnerCentrePoint(n1specNoInject,BinSizePixels,BinCentralPixel)    

    n2specBinned = BinnnerCentrePoint(n2spec,BinSizePixels,BinCentralPixel) 

    n2specBinnedNoInject = BinnnerCentrePoint(n2specNoInject,BinSizePixels,BinCentralPixel) 

    n2specBinnedNoInject = BinnnerCentrePoint(n2specNoInject,BinSizePixels,BinCentralPixel)

    n3specBinned = BinnnerCentrePoint(n3spec,BinSizePixels,BinCentralPixel)  

    n3specBinnedNoInject = BinnnerCentrePoint(n3specNoInject,BinSizePixels,BinCentralPixel) 

    n4specBinned = BinnnerCentrePoint(n4spec,BinSizePixels,BinCentralPixel) 

    n4specBinnedNoInject = BinnnerCentrePoint(n4specNoInject,BinSizePixels,BinCentralPixel)   
    
    ## Combine the nights  
  
  
    combined = (DataSetWeight[0]*n1spec+DataSetWeight[1]*n2spec+DataSetWeight[2]*n3spec+DataSetWeight[3]*n4spec)/np.sum(DataSetWeight)
    combined_binned = BinnnerCentrePoint(combined,BinSizePixels,BinCentralPixel)   
    
    combinedWithout1 = (DataSetWeight[1]*n2spec+DataSetWeight[2]*n3spec+DataSetWeight[3]*n4spec)/np.sum(DataSetWeight[1:])
    combinedWithout1_binned = BinnnerCentrePoint(combinedWithout1,BinSizePixels,BinCentralPixel)   
    
    NoInjectedCombined = (DataSetWeight[0]*n1specNoInject + DataSetWeight[1]*n2specNoInject + DataSetWeight[2]*n3specNoInject + DataSetWeight[3]*n4specNoInject)/np.sum(DataSetWeight)
    NoInjectedCombined_binned = BinnnerCentrePoint(NoInjectedCombined,BinSizePixels,BinCentralPixel) 
    
    NoInjectedCombinedWithout1 = (DataSetWeight[1]*n2specNoInject + DataSetWeight[2]*n3specNoInject + DataSetWeight[3]*n4specNoInject)/np.sum(DataSetWeight[1:])
    NoInjectedCombinedWithout1_binned = BinnnerCentrePoint(NoInjectedCombinedWithout1,BinSizePixels,BinCentralPixel) 
    
    
### No injected combined 
    
    combinedNoInject = (DataSetWeight[0]*n1specNoInject+DataSetWeight[1]*n2specNoInject+DataSetWeight[2]*n3specNoInject+DataSetWeight[3]*n4specNoInject)/np.sum(DataSetWeight)
    combined_binnedNoInject = BinnnerCentrePoint(combinedNoInject,BinSizePixels,BinCentralPixel)   
    
    combinedWithout1NoInject = (DataSetWeight[1]*n2specNoInject+DataSetWeight[2]*n3specNoInject+DataSetWeight[3]*n4specNoInject)/np.sum(DataSetWeight[1:])
    combinedWithout1_binnedNoInject = BinnnerCentrePoint(combinedWithout1NoInject,BinSizePixels,BinCentralPixel)   

    
    ### subtracted no injected 
    
    n1specNoInjectSub = n1spec - n1specNoInject     
    n1specBinnedNoInjectSub = BinnnerCentrePoint(n1specNoInjectSub,BinSizePixels,BinCentralPixel)   
    
    n2specNoInjectSub = n2spec - n2specNoInject    
    n2specBinnedNoInjectSub = BinnnerCentrePoint(n2specNoInjectSub,BinSizePixels,BinCentralPixel)   
    
    n3specNoInjectSub = n3spec - n3specNoInject    
    n3specBinnedNoInjectSub = BinnnerCentrePoint(n3specNoInjectSub,BinSizePixels,BinCentralPixel)   
       
    n4specNoInjectSub = n4spec - n4specNoInject     
    n4specBinnedNoInjectSub = BinnnerCentrePoint(n4specNoInjectSub,BinSizePixels,BinCentralPixel)   
    
    CombinedNoInjectSubbed = (DataSetWeight[0]*n1specNoInjectSub + DataSetWeight[1]*n2specNoInjectSub + DataSetWeight[2]*n3specNoInjectSub + DataSetWeight[3]*n4specNoInjectSub)/np.sum(DataSetWeight)
    CombinedNoInjectSubbedWithout1 = (DataSetWeight[1]*n2specNoInjectSub + DataSetWeight[2]*n3specNoInjectSub + DataSetWeight[3]*n4specNoInjectSub)/np.sum(DataSetWeight[1:])
    
    CombinedNoInjectSubbedBinned = BinnnerCentrePoint(CombinedNoInjectSubbed,BinSizePixels,BinCentralPixel)   
    CombinedNoInjectSubbedWithout1Binned = BinnnerCentrePoint(CombinedNoInjectSubbedWithout1,BinSizePixels,BinCentralPixel)   

      
####

    avg_CaII_wavelength = (8.49801800e+03+8.54208900e+03+8.66214000e+03)/3.0  
   
    #plot_wx = (np.arange(len(planetspec2012_01_27))-399.5)*0.01
    #plot_wx = (np.arange(len(Feb20))-len(Feb20)/2.0)*deltaw
   
    plot_wx = (np.arange(len(n4spec))-len(n4spec)/2.0)*deltaw
    #plot_wx = ((np.arange(len(n4spec))-len(n4spec)/2.0)-2)*deltaw
    plot_vx = (plot_wx/avg_CaII_wavelength)*c

    snrlist = []
    signallist = []
    noiselist = []
    
    origsnrlist = []
    


    
    ## Divide by about 6 to get A 
    
    plt.figure()
    plt.title('n1')
    plt.plot(plot_vx,n1spec,'k-')
    plt.plot(plot_vx,n1specBinned,'k:')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)    
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n1 no injected signal')
    plt.plot(plot_vx,n1specNoInject,'k:')
    plt.plot(plot_vx,n1specBinnedNoInject,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)    
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n2')
    plt.plot(plot_vx,n2spec,'k:')
    plt.plot(plot_vx,n2specBinned,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n2 no injected signal')
    plt.plot(plot_vx,n2specNoInject,'k:')
    plt.plot(plot_vx,n2specBinnedNoInject,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n3')
    plt.plot(plot_vx,n3spec,'k:')
    plt.plot(plot_vx,n3specBinned,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n3 no injected signal')
    plt.plot(plot_vx,n3specNoInject,'k:')
    plt.plot(plot_vx,n3specBinnedNoInject,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n4')
    plt.plot(plot_vx,n4spec,'k:')
    plt.plot(plot_vx,n4specBinned,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n4 no injected signal')
    plt.plot(plot_vx,n4specNoInject,'k:')
    plt.plot(plot_vx,n4specBinnedNoInject,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
###########   No injected subtracted 
    
    plt.figure()
    plt.title('n1 no injected subtracted')
    plt.plot(plot_vx,n1specNoInjectSub,'k:')
    plt.plot(plot_vx,n1specBinnedNoInjectSub,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)    
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n2 no injected subtracted')
    plt.plot(plot_vx,n2specNoInjectSub,'k:')
    plt.plot(plot_vx,n2specBinnedNoInjectSub,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n3 no injected subtracted')
    plt.plot(plot_vx,n3specNoInjectSub,'k:')
    plt.plot(plot_vx,n3specBinnedNoInjectSub,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('n4 no injected subtracted')
    plt.plot(plot_vx,n4specNoInjectSub,'k:')
    plt.plot(plot_vx,n4specBinnedNoInjectSub,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    



    plt.figure()
    plt.title('combined')
    plt.plot(plot_vx,combined,'k:')
    plt.plot(plot_vx,combined_binned,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    
    plt.figure()
    plt.title('combined without night 1')
    plt.plot(plot_vx,combinedWithout1,'k:')
    plt.plot(plot_vx,combinedWithout1_binned,'k-')
    plt.xlim(PlotXLims)
    plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
   
    n1SNR = calculate_SNR(n1spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n2SNR = calculate_SNR(n2spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n3SNR = calculate_SNR(n3spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n4SNR = calculate_SNR(n4spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)

     
    CombinedSNR = calculate_SNR(combined,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    
    CombinedNo1SNR = calculate_SNR_BetterPlot(combinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    CombinedNoNoInjectedSNR = calculate_SNR_BetterPlot(NoInjectedCombinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    
    ## The binned SNRS 
    
    n1SNRBinned = calculate_SNR(n1specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n2SNRBinned = calculate_SNR(n2specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n3SNRBinned = calculate_SNR(n3specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n4SNRBinned = calculate_SNR(n4specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)

     
    CombinedSNRBinned = calculate_SNR(combined_binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    
    CombinedNo1SNRBinned = calculate_SNR_BetterPlot(combinedWithout1_binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    
    
    ##### Subtracting the spectrum without an injected signal     
    CombinedNo1NoInjectSubbedSNR = calculate_SNR_BetterPlot(CombinedNoInjectSubbedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    
    
    ActualNoInjectedSubbedSNR = CombinedNo1NoInjectSubbedSNR[1]/CombinedNo1SNR[2]
    
    CombinedNoInjectSubbedWithout1BinnedSNR =  calculate_SNR(CombinedNoInjectSubbedWithout1Binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
 
    #ActualBinnedNoInjectedSubbedSNR = CombinedNoInjectSubbedWithout1BinnedSNR[1]/(2*CombinedSNRBinned[2])  ## *2 noise since the signal is subtracted 
    ActualBinnedNoInjectedSubbedSNR = CombinedNoInjectSubbedWithout1BinnedSNR[1]/(CombinedSNRBinned[2])  ## *2 noise since the signal is subtracted 
 
 
    ### I think this combined SNR is the best measure of noise 
    print(CombinedNo1SNR)
    
    No1InjectedSNR.append(CombinedNo1SNR[0])
    No1BinnedInjectedSNR.append(CombinedNo1SNRBinned[0])
    
    
STDOfInjectionPartBeforeDownWeightingNoisyColumns = np.mean([np.mean(n1InjectedPartsSTDNoInject),np.mean(n2InjectedPartsSTDNoInject),np.mean(n3InjectedPartsSTDNoInject),np.mean(n4InjectedPartsSTDNoInject)])/np.sqrt(2)  ## Since the two lines are combined 


##### Plot binned and unbinned columns with guiding lines  

dashparams = (4,1.5)

fig, axes = plt.subplots(nrows=4, ncols=2,figsize=(7.27,10.69))



axes[0][0].plot(plot_vx,n2spec,'k--',dashes=(dashparams[0], dashparams[1]))
axes[0][0].plot(plot_vx,n2specNoInject,'k')
axes[0][0].plot([PlotXLims[0],PlotXLims[1]],[-EffectiveInjectionStrength,-EffectiveInjectionStrength],'b')
axes[0][0].plot([PlotXLims[0],PlotXLims[1]],[-np.std(n2specNoInject[500-5:500+5]),-np.std(n2specNoInject[500-5:500+5])],'r')
axes[0][0].plot([PlotXLims[0],PlotXLims[1]],2*np.array([-np.std(n2specNoInject[500-5:500+5]),-np.std(n2specNoInject[500-5:500+5])]),'r')
axes[0][0].plot([PlotXLims[0],PlotXLims[1]],3*np.array([-np.std(n2specNoInject[500-5:500+5]),-np.std(n2specNoInject[500-5:500+5])]),'r')
#axes[0][0].plot([PlotXLims[0],PlotXLims[1]],[np.std(n2specNoInject[500-5:500+5]),np.std(n2specNoInject[500-5:500+5])],'r')
axes[0][0].set_xlim(PlotXLims)
axes[0][0].set_ylim(PlotYLims)


axes[0][1].plot(plot_vx,n2specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[0][1].plot(plot_vx,n2specBinnedNoInject,'k')
axes[0][1].set_xlim(PlotXLims)
axes[0][1].set_ylim(PlotYLims)

axes[1][0].plot(plot_vx,n3spec,'k--',dashes=(dashparams[0], dashparams[1]))
axes[1][0].plot(plot_vx,n3specNoInject,'k')
axes[1][0].plot([PlotXLims[0],PlotXLims[1]],[-EffectiveInjectionStrength,-EffectiveInjectionStrength],'b')
axes[1][0].plot([PlotXLims[0],PlotXLims[1]],[-np.std(n3specNoInject[500-5:500+5]),-np.std(n3specNoInject[500-5:500+5])],'r')
axes[1][0].plot([PlotXLims[0],PlotXLims[1]],2*np.array([-np.std(n3specNoInject[500-5:500+5]),-np.std(n3specNoInject[500-5:500+5])]),'r')
axes[1][0].plot([PlotXLims[0],PlotXLims[1]],3*np.array([-np.std(n3specNoInject[500-5:500+5]),-np.std(n3specNoInject[500-5:500+5])]),'r')

#axes[1][0].plot([PlotXLims[0],PlotXLims[1]],[np.std(n3specNoInject[500-5:500+5]),np.std(n3specNoInject[500-5:500+5])],'r')
axes[1][0].set_xlim(PlotXLims)
axes[1][0].set_ylim(PlotYLims)

axes[1][1].plot(plot_vx,n3specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[1][1].plot(plot_vx,n3specBinnedNoInject,'k')
axes[1][1].set_xlim(PlotXLims)
axes[1][1].set_ylim(PlotYLims)

axes[2][0].plot(plot_vx,n4spec,'k--',dashes=(dashparams[0], dashparams[1]))
axes[2][0].plot(plot_vx,n4specNoInject,'k')
axes[2][0].plot([PlotXLims[0],PlotXLims[1]],[-EffectiveInjectionStrength,-EffectiveInjectionStrength],'b')
axes[2][0].plot([PlotXLims[0],PlotXLims[1]],[-np.std(n4specNoInject[500-5:500+5]),-np.std(n4specNoInject[500-5:500+5])],'r')
axes[2][0].plot([PlotXLims[0],PlotXLims[1]],2*np.array([-np.std(n4specNoInject[500-5:500+5]),-np.std(n4specNoInject[500-5:500+5])]),'r')
axes[2][0].plot([PlotXLims[0],PlotXLims[1]],3*np.array([-np.std(n4specNoInject[500-5:500+5]),-np.std(n4specNoInject[500-5:500+5])]),'r')
#axes[2][0].plot([PlotXLims[0],PlotXLims[1]],[np.std(n4specNoInject[500-5:500+5]),np.std(n4specNoInject[500-5:500+5])],'r')
axes[2][0].set_xlim(PlotXLims)
axes[2][0].set_ylim(PlotYLims)

axes[2][1].plot(plot_vx,n4specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[2][1].plot(plot_vx,n4specBinnedNoInject,'k')
axes[2][1].set_xlim(PlotXLims)
axes[2][1].set_ylim(PlotYLims)

axes[3][0].plot(plot_vx,combinedWithout1,'k--',dashes=(dashparams[0], dashparams[1]))
axes[3][0].plot(plot_vx,NoInjectedCombinedWithout1,'k')
axes[3][0].plot([PlotXLims[0],PlotXLims[1]],[-EffectiveInjectionStrength,-EffectiveInjectionStrength],'b')
axes[3][0].plot([PlotXLims[0],PlotXLims[1]],[-np.std(NoInjectedCombinedWithout1[500-5:500+5]),-np.std(NoInjectedCombinedWithout1[500-5:500+5])],'r')
axes[3][0].plot([PlotXLims[0],PlotXLims[1]],2*np.array([-np.std(NoInjectedCombinedWithout1[500-5:500+5]),-np.std(NoInjectedCombinedWithout1[500-5:500+5])]),'r')
axes[3][0].plot([PlotXLims[0],PlotXLims[1]],3*np.array([-np.std(NoInjectedCombinedWithout1[500-5:500+5]),-np.std(NoInjectedCombinedWithout1[500-5:500+5])]),'r')
#axes[3][0].plot([PlotXLims[0],PlotXLims[1]],[np.std(NoInjectedCombinedWithout1[500-5:500+5]),np.std(NoInjectedCombinedWithout1[500-5:500+5])],'r')
axes[3][0].set_xlim(PlotXLims)
axes[3][0].set_ylim(PlotYLimsCombined)


axes[3][1].plot(plot_vx,combinedWithout1_binned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[3][1].plot(plot_vx,NoInjectedCombinedWithout1_binned,'k')
axes[3][1].set_xlim(PlotXLims)
axes[3][1].set_ylim(PlotYLimsCombined)

for i in range(3):
    plt.setp(axes[i][0].get_xticklabels(), visible=False)
    plt.setp(axes[i][1].get_xticklabels(), visible=False)
    plt.setp(axes[i][1].get_yticklabels(), visible=False)
plt.setp(axes[3][1].get_yticklabels(), visible=False)

plt.text(0.5, 1.2, 'Unbinned',
         horizontalalignment='center',
         fontsize=12,
         transform = axes[0][0].transAxes)
         
plt.text(0.5, 1.2, 'Binned',
         horizontalalignment='center',
         fontsize=12,
         transform = axes[0][1].transAxes)
         
plt.text(1, 1.08, 'Night 2',
         horizontalalignment='center',
         fontsize=12,
         transform = axes[0][0].transAxes)
         
plt.text(1, 1.08, 'Night 3',
         horizontalalignment='center',
         fontsize=12,
         transform = axes[1][0].transAxes)
         
plt.text(1, 1.08, 'Night 4',
         horizontalalignment='center',
         fontsize=12,
         transform = axes[2][0].transAxes)
         
        
plt.text(1, 1.08, 'Weigted mean',
         horizontalalignment='center',
         fontsize=12,
         transform = axes[3][0].transAxes)
         
         
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.ylabel('absorption relative to stellar spectrum',labelpad=20)
plt.xlabel(r'radial velocity (km s$^{-1}$)')
         
fig.subplots_adjust(hspace=0.3, wspace=0.0)


#########################

#### Checking the width of the signal 

plt.figure()
plt.plot(n1spec,label='1'),plt.plot(n2spec,label='2'),plt.plot(n3spec,label='3'),plt.plot(n4spec,label='4'),plt.legend()

#### Fit the damn gaussians 

def FitTheDamnGaussians(spec):
    
      
    minindex = np.argmin(spec)
    
    fitlinewidth = 8
    
    linetofit = spec[minindex - fitlinewidth:minindex + fitlinewidth+1]
    
    x = np.arange(len(linetofit))
    
    g_init = models.Gaussian1D(amplitude=0.25, mean=8, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, linetofit)
    
    xx = np.linspace(x[0],x[-1],1e4)
    
    plt.figure()
    plt.plot(x,linetofit,'b')
    plt.plot(xx,g(xx),'r')
    
    std = g.stddev.value
    
    plt.title('standard deviation: %f px' % (std))
    
    return fit_g
    
    
    
g1 = FitTheDamnGaussians(n1spec)
    
#spec = n1spec
#    
#minindex = np.argmin(spec)
#
#fitlinewidth = 8
#
#linetofit = spec[minindex - fitlinewidth:minindex + fitlinewidth+1]
#
#x = np.arange(len(linetofit))
#
#g_init = models.Gaussian1D(amplitude=0.25, mean=8, stddev=1.)
#fit_g = fitting.LevMarLSQFitter()
#g = fit_g(g_init, x, linetofit)
#
#xx = np.linspace(x[0],x[-1],1e4)
#
#plt.figure()
#plt.plot(x,linetofit,'b')
#plt.plot(xx,g(xx),'r')
#
#std = g.stddev.value
    

    

    
    
    
    
    
    

    
    
    
    
    




#########################################

###### Plot binned and unbinned columns 
#
#dashparams = (4,1.5)
#
#fig, axes = plt.subplots(nrows=4, ncols=2,figsize=(7.27,10.69))
#
#
#
#axes[0][0].plot(plot_vx,n2spec,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[0][0].plot(plot_vx,n2specNoInject,'k')
#axes[0][0].set_xlim(PlotXLims)
#axes[0][0].set_ylim(PlotYLims)
#
#
#axes[0][1].plot(plot_vx,n2specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[0][1].plot(plot_vx,n2specBinnedNoInject,'k')
#axes[0][1].set_xlim(PlotXLims)
#axes[0][1].set_ylim(PlotYLims)
#
#axes[1][0].plot(plot_vx,n3spec,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[1][0].plot(plot_vx,n3specNoInject,'k')
#axes[1][0].set_xlim(PlotXLims)
#axes[1][0].set_ylim(PlotYLims)
#
#axes[1][1].plot(plot_vx,n3specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[1][1].plot(plot_vx,n3specBinnedNoInject,'k')
#axes[1][1].set_xlim(PlotXLims)
#axes[1][1].set_ylim(PlotYLims)
#
#axes[2][0].plot(plot_vx,n4spec,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[2][0].plot(plot_vx,n4specNoInject,'k')
#axes[2][0].set_xlim(PlotXLims)
#axes[2][0].set_ylim(PlotYLims)
#
#axes[2][1].plot(plot_vx,n4specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[2][1].plot(plot_vx,n4specBinnedNoInject,'k')
#axes[2][1].set_xlim(PlotXLims)
#axes[2][1].set_ylim(PlotYLims)
#
#axes[3][0].plot(plot_vx,combinedWithout1,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[3][0].plot(plot_vx,NoInjectedCombinedWithout1,'k')
#axes[3][0].set_xlim(PlotXLims)
#axes[3][0].set_ylim(PlotYLimsCombined)
#
#
#axes[3][1].plot(plot_vx,combinedWithout1_binned,'k--',dashes=(dashparams[0], dashparams[1]))
#axes[3][1].plot(plot_vx,NoInjectedCombinedWithout1_binned,'k')
#axes[3][1].set_xlim(PlotXLims)
#axes[3][1].set_ylim(PlotYLimsCombined)
#
#for i in range(3):
#    plt.setp(axes[i][0].get_xticklabels(), visible=False)
#    plt.setp(axes[i][1].get_xticklabels(), visible=False)
#    plt.setp(axes[i][1].get_yticklabels(), visible=False)
#plt.setp(axes[3][1].get_yticklabels(), visible=False)
#
#plt.text(0.5, 1.2, 'Unbinned',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][0].transAxes)
#         
#plt.text(0.5, 1.2, 'Binned',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][1].transAxes)
#         
#plt.text(1, 1.08, 'Night 2',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 3',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[1][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 4',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[2][0].transAxes)
#         
#        
#plt.text(1, 1.08, 'Weigted mean',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[3][0].transAxes)
#         
#         
#fig.add_subplot(111, frameon=False)
## hide tick and tick label of the big axes
#plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#plt.grid(False)
#plt.ylabel('absorption relative to stellar spectrum',labelpad=20)
#plt.xlabel(r'radial velocity (km s$^{-1}$)')
#         
#fig.subplots_adjust(hspace=0.3, wspace=0.0)




################### 
#### Without the combined night 1   Old plot with binned and not binned column 

#fig, axes = plt.subplots(nrows=5, ncols=2,figsize=(7.27,10.69))
#
#axes[0][0].plot(plot_vx,n1specNoInject,'k:')
#axes[0][0].plot(plot_vx,n1spec,'k')
#axes[0][0].set_xlim(PlotXLims)
#axes[0][0].set_ylim(PlotYLims)
#
#
#axes[0][1].plot(plot_vx,n1specBinnedNoInject,'k:')
#axes[0][1].plot(plot_vx,n1specBinned,'k')
#axes[0][1].set_xlim(PlotXLims)
#axes[0][1].set_ylim(PlotYLims)
#
#axes[1][0].plot(plot_vx,n2specNoInject,'k:')
#axes[1][0].plot(plot_vx,n2spec,'k')
#axes[1][0].set_xlim(PlotXLims)
#axes[1][0].set_ylim(PlotYLims)
#
#
#axes[1][1].plot(plot_vx,n2specBinnedNoInject,'k:')
#axes[1][1].plot(plot_vx,n2specBinned,'k')
#axes[1][1].set_xlim(PlotXLims)
#axes[1][1].set_ylim(PlotYLims)
#
#axes[2][0].plot(plot_vx,n3specNoInject,'k:')
#axes[2][0].plot(plot_vx,n3spec,'k')
#axes[2][0].set_xlim(PlotXLims)
#axes[2][0].set_ylim(PlotYLims)
#
#
#axes[2][1].plot(plot_vx,n3specBinnedNoInject,'k:')
#axes[2][1].plot(plot_vx,n3specBinned,'k')
#axes[2][1].set_xlim(PlotXLims)
#axes[2][1].set_ylim(PlotYLims)
#
#axes[3][0].plot(plot_vx,n4specNoInject,'k:')
#axes[3][0].plot(plot_vx,n4spec,'k')
#axes[3][0].set_xlim(PlotXLims)
#axes[3][0].set_ylim(PlotYLims)
#
#axes[3][1].plot(plot_vx,n4specBinnedNoInject,'k:')
#axes[3][1].plot(plot_vx,n4specBinned,'k')
#axes[3][1].set_xlim(PlotXLims)
#axes[3][1].set_ylim(PlotYLims)
#
#axes[4][0].plot(plot_vx,NoInjectedCombinedWithout1,'k:')
#axes[4][0].plot(plot_vx,combinedWithout1,'k')
#axes[4][0].set_xlim(PlotXLims)
#axes[4][0].set_ylim(PlotYLims)
#
#
#axes[4][1].plot(plot_vx,NoInjectedCombinedWithout1_binned,'k:')
#axes[4][1].plot(plot_vx,combinedWithout1_binned,'k')
#axes[4][1].set_xlim(PlotXLims)
#axes[4][1].set_ylim(PlotYLims)
#
#for i in range(4):
#    plt.setp(axes[i][0].get_xticklabels(), visible=False)
#    plt.setp(axes[i][1].get_xticklabels(), visible=False)
#    plt.setp(axes[i][1].get_yticklabels(), visible=False)
#plt.setp(axes[4][1].get_yticklabels(), visible=False)
#
#plt.text(0.5, 1.2, 'Unbinned',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][0].transAxes)
#         
#plt.text(0.5, 1.2, 'Binned',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][1].transAxes)
#
#
#plt.text(1, 1.08, 'Night 1',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 2',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[1][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 3',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[2][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 4',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[3][0].transAxes)
#         
#        
#plt.text(1, 1.08, 'Weigted mean of all nights excluding Night 1',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[4][0].transAxes)
#         
#plt.text(-0.3, 0.55, 'transmission spectrum (arbitrary units)',
#         verticalalignment='center',
#         fontsize=12,
#         rotation=90,
#         transform = axes[2][0].transAxes)
#         
#plt.text(1,-0.4 , r'radial velocity relative to combined line (km s$^{-1}$)',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[4][0].transAxes)
#         
#fig.subplots_adjust(hspace=0.3, wspace=0.0)

##############################################################




###############################################################

#### Including a combined plot with night 1 

#fig, axes = plt.subplots(nrows=6, ncols=2,figsize=(7.27,10.69))
#
#axes[0][0].plot(plot_vx,n1specNoInject,'k')
#axes[0][0].plot(plot_vx,n1spec,'k:')
#axes[0][0].set_xlim(PlotXLims)
#axes[0][0].set_ylim(PlotYLims)
#
#
#axes[0][1].plot(plot_vx,n1specBinnedNoInject,'k')
#axes[0][1].plot(plot_vx,n1specBinned,'k:')
#axes[0][1].set_xlim(PlotXLims)
#axes[0][1].set_ylim(PlotYLims)
#
#axes[1][0].plot(plot_vx,n2specNoInject,'k')
#axes[1][0].plot(plot_vx,n2spec,'k:')
#axes[1][0].set_xlim(PlotXLims)
#axes[1][0].set_ylim(PlotYLims)
#
#
#axes[1][1].plot(plot_vx,n2specBinnedNoInject,'k')
#axes[1][1].plot(plot_vx,n2specBinned,'k:')
#axes[1][1].set_xlim(PlotXLims)
#axes[1][1].set_ylim(PlotYLims)
#
#axes[2][0].plot(plot_vx,n3specNoInject,'k')
#axes[2][0].plot(plot_vx,n3spec,'k:')
#axes[2][0].set_xlim(PlotXLims)
#axes[2][0].set_ylim(PlotYLims)
#
#
#axes[2][1].plot(plot_vx,n3specBinnedNoInject,'k')
#axes[2][1].plot(plot_vx,n3specBinned,'k:')
#axes[2][1].set_xlim(PlotXLims)
#axes[2][1].set_ylim(PlotYLims)
#
#axes[3][0].plot(plot_vx,n4specNoInject,'k')
#axes[3][0].plot(plot_vx,n4spec,'k:')
#axes[3][0].set_xlim(PlotXLims)
#axes[3][0].set_ylim(PlotYLims)
#
#axes[3][1].plot(plot_vx,n4specBinnedNoInject,'k')
#axes[3][1].plot(plot_vx,n4specBinned,'k:')
#axes[3][1].set_xlim(PlotXLims)
#axes[3][1].set_ylim(PlotYLims)
#
#
#axes[4][0].plot(plot_vx,combinedNoInject,'k')
#axes[4][0].plot(plot_vx,combined,'k:')
#axes[4][0].set_xlim(PlotXLims)
#axes[4][0].set_ylim(PlotYLims)
#
#
#axes[4][1].plot(plot_vx,combined_binnedNoInject,'k')
#axes[4][1].plot(plot_vx,combined_binned,'k:')
#axes[4][1].set_xlim(PlotXLims)
#axes[4][1].set_ylim(PlotYLims)
#
#
#axes[5][0].plot(plot_vx,NoInjectedCombinedWithout1,'k')
#axes[5][0].plot(plot_vx,combinedWithout1,'k:')
#axes[5][0].set_xlim(PlotXLims)
#axes[5][0].set_ylim(PlotYLims)
#
#
#axes[5][1].plot(plot_vx,NoInjectedCombinedWithout1_binned,'k')
#axes[5][1].plot(plot_vx,combinedWithout1_binned,'k:')
#axes[5][1].set_xlim(PlotXLims)
#axes[5][1].set_ylim(PlotYLims)
#
#for i in range(5):
#    plt.setp(axes[i][0].get_xticklabels(), visible=False)
#    plt.setp(axes[i][1].get_xticklabels(), visible=False)
#    plt.setp(axes[i][1].get_yticklabels(), visible=False)
#plt.setp(axes[5][1].get_yticklabels(), visible=False)
#
#plt.text(0.5, 1.2, 'Unbinned',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][0].transAxes)
#         
#plt.text(0.5, 1.2, 'Binned',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][1].transAxes)
#
#plt.text(1, 1.08, 'Night 1',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[0][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 2',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[1][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 3',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[2][0].transAxes)
#         
#plt.text(1, 1.08, 'Night 4',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[3][0].transAxes)
#         
#plt.text(1, 1.08, 'Weigted mean of all nights',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[4][0].transAxes)
#         
#plt.text(1, 1.08, 'Weigted mean of all nights excluding Night 1',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[5][0].transAxes)
#         
#plt.text(-0.3, 1.0, 'transmission spectrum (arbitrary units)',
#         verticalalignment='center',
#         fontsize=12,
#         rotation=90,
#         transform = axes[3][0].transAxes)
#         
#plt.text(1,-0.4 , 'radial velocity relative to combined line (km/s)',
#         horizontalalignment='center',
#         fontsize=12,
#         transform = axes[5][0].transAxes)
#         
#fig.subplots_adjust(hspace=0.3, wspace=0.0)



    
