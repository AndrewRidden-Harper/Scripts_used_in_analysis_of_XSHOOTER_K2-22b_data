# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:30:26 2015

@author: arh

Think need an interval of 800 pixels either side 

"""

from NaD_inject_model_spectrum import call_model_injector
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt 
import scipy
from DopplerShiftToPlanetFrame import shift_to_planet_frame, findNANedges, whole_pixel_shift_to_planet_frame
from K2_22bJDToPhaseFuncCustomTransLims import JDtophase_radial_vel
import os
from pcasub import pcasub 
from NaD_VALD_MakeSpec import makespec
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from astropy.stats import sigma_clip, mad_std
import pandas as pd 

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
    signal = InjectionStrength#*(recoveredsignal/0.073134605451126752)  ## option to scale relative to the 0 pca minimum
    
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
    signal = InjectionStrength#*(recoveredsignal/0.073134605451126752)  ## option to scale relative to the 0 pca minimum 
    
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


    if name == 'night1':
        CaTrip1PixCenter = 396
        CaTrip2PixCenter = 425
        
        L1IndicesForIteration = np.array([375, 377, 380, 382, 387, 390, 394, 397, 402, 405, 408, 411, 414])
        L2IndicesForIteration = np.array([405, 407, 409, 412, 417, 420, 423, 427, 432, 435, 438, 441, 444])
              
        
    if name == 'night2':
        CaTrip1PixCenter = 396
        CaTrip2PixCenter = 426
        
        L1IndicesForIteration = np.array([376, 378, 380, 383, 387, 391, 394, 398, 402, 405, 409, 411, 414])
        L2IndicesForIteration = np.array([405, 408, 410, 413, 417, 421, 424, 428, 432, 435, 438, 441, 444])

        
    if name == 'night3':
        CaTrip1PixCenter = 395
        CaTrip2PixCenter = 425
        
        L1IndicesForIteration = np.array([377, 379, 382, 385, 389, 393, 397, 400, 405, 408, 411, 413, 416])
        L2IndicesForIteration = np.array([406, 409, 412, 415, 419, 423, 426, 430, 435, 438, 441, 443, 446])


        
    if name == 'night4':
        CaTrip1PixCenter = 395
        CaTrip2PixCenter = 425
        
        L1IndicesForIteration = np.array([373, 375, 377, 379, 383, 386, 389, 393, 398, 401, 405, 407, 411])
        L2IndicesForIteration = np.array([403, 404, 406, 409, 413, 416, 419, 423, 428, 432, 434, 438, 441])
        
    LinePartHalfExtent = 380
    
    numrows,numcols = np.shape(data)    
    #Qpcatoremove = 0    
    
    
    if not os.path.exists('processing/%s/pca%d'%(name,pcatoremove)):
        os.makedirs('processing/%s/pca%d'%(name,pcatoremove))
    
    if ToInject == True:
        toprocess = 'injected'
    
    if ToInject == False:
        toprocess = 'no_injected'
    
    species = 'VALD_Na'
    
    flux = injectStrength
       
    #Create the model template if it does not exist already 
    if not os.path.exists('model_spectra/%s/NaD_F%.2E_s%f' % (species,flux,width)):
        makespec(species='%s'%(species),flux=flux,sigma=width)     
    
    d = data    
    w = wavelength
    
    if toprocess == 'injected':
        #d = call_model_injector(d,w,jd,species,flux,width,vshift=vsys+8)  # This vshift includes the system velocity and and shift to seperate 
        d, injected_spec = call_model_injector(d,w,jd,species,flux,width,SpecialTransLims,vshift=ModelRadialVelocity)
        #d = injected_spec
        
        
        writedir = 'NaDLines/%s/pca%d/Injected_F%.2E' % (name,pcatoremove,flux)     
        
        if not os.path.exists(writedir):
            os.makedirs(writedir)
            
    if toprocess == 'no_injected':
        
        writedir = 'NaDLines/%s/pca%d' % (name,pcatoremove)
        
        if not os.path.exists(writedir):
            os.makedirs(writedir)
            
    if ToInject == True:   
        pyfits.writeto('%s/InjectedSpectrum.fits'%(writedir),injected_spec,clobber=True)
    
    
    d2 = np.copy(d)
    
    for i in range(numcols):
        
        col = d2[:,i]
        
        mask = np.zeros(len(col))
        mask[SpecialTransLims[0]:SpecialTransLims[1]+1] = 1  ### 1 is a masked value 
                
        MaskedArray = np.ma.masked_array(col,mask)        
                
        filtered_data = sigma_clip(MaskedArray, sigma=5, iters=1, stdfunc=mad_std)
        
        #filtered_data = sigma_clip(d2[:,i], sigma=5, iters=1, stdfunc=mad_std)
        d2[:,i] = d[:,i]/np.mean(filtered_data)   
        
        #d2[:,i] = d[:,i]/np.mean(col)
    
    
    d2[np.isnan(d2)] = 1.0
    d2[np.isinf(d2)] = 1.0
    
    pyfits.writeto('%s/d2.fits'%(writedir),d2,clobber=True)
    print('Writedir: %s'%(writedir))
    
    d1b = np.copy(d2)
    
    d2_planetframe = shift_to_planet_frame(d2,w,jd,vorbital='calculate')
    pyfits.writeto('%s/d2PlanetFrame.fits'%(writedir),d2_planetframe,clobber=True)
    
    
    d2b = d1b#/np.median(d1b,axis=0)
    
    residual_std = np.std(d2b,axis=0)
    SNR = 1.0/residual_std
    
       
    pyfits.writeto('%s/d2b.fits'%(writedir),d2b,clobber=True)
    
       
    d3 = np.copy(d2b)
    
   
    #### Use PCA to remove the telluric lines and sublte trends (either do this or fit linear to each column)    
    pyfits.writeto('%s/d3_prepcasub.fits'%(writedir),d3,clobber=True)    
        
    if componets_to_remove>0:
      
        d3 = pcasub(d2b,componets_to_remove)
        
    else: 

        d3 = d3 - np.mean(d3)
    
    pyfits.writeto('%s/d3_pcasub.fits'%(writedir),d3,clobber=True)
        
    pyfits.writeto('%s/pre_col_std_scaled_d3.fits'%(writedir),d3,clobber=True)
    
    residual_col_std = np.std(d3,axis=0)
    scaled_resid_col_std = residual_col_std*(1.0/np.median(residual_col_std))
    
    d3_precolweighting = np.copy(d3)  
    


    pyfits.writeto('%s/d3_pcasub_colweighted.fits'%(writedir),d3,clobber=True)
    
    #### Combine the Na D lines  

    
    pyfits.writeto('%s/CaTripLine1Part.fits'%(writedir),d3[:,CaTrip1PixCenter-LinePartHalfExtent:CaTrip1PixCenter+LinePartHalfExtent],clobber=True)
    pyfits.writeto('%s/CaTripLine2Part.fits'%(writedir),d3[:,CaTrip2PixCenter-LinePartHalfExtent:CaTrip2PixCenter+LinePartHalfExtent],clobber=True)
    
    
    CaTrip1 = d3[:,CaTrip1PixCenter-LinePartHalfExtent:CaTrip1PixCenter+LinePartHalfExtent]
    CaTrip2 = d3[:,CaTrip2PixCenter-LinePartHalfExtent:CaTrip2PixCenter+LinePartHalfExtent]
    
    ###########################################################
    ############### New noise scaling 
    
    nrows,ncols = np.shape(CaTrip1)
    
    ScaledCaTrip1 = np.copy(CaTrip1)
    
    L1sigma = np.std(CaTrip1,axis = 0)
    L2sigma = np.std(CaTrip2,axis = 0)
    

    L1SNR = (1/L1sigma)#**2
    L2SNR = (1/L2sigma)#**2
    
    RollingMeanWindow = 25 
    
    L1SNR = pd.rolling_mean(pd.DataFrame(L1SNR),window=RollingMeanWindow).values.flatten()
    L2SNR = pd.rolling_mean(pd.DataFrame(L2SNR),window=RollingMeanWindow).values.flatten()
    
    L1SNR[0:RollingMeanWindow-1] = L1SNR[RollingMeanWindow-1]
    L2SNR[0:RollingMeanWindow-1] = L2SNR[RollingMeanWindow-1]
    
    DeltaPix = np.diff(L1SNR)
    
#    plt.figure()
#    plt.title('difference in SNR between neighbouring pixels')
#    plt.plot(DeltaPix)
#    plt.ylabel('difference in SNR')
#    
#    
#    
#    plt.figure()
#    plt.plot(L1SNR)
#    plt.title('signal to noise of line 1')    
        
        
#    L1SNRSquared = 1/np.std(CaTrip1,axis = 0)
#    L2SNRSquared = 1/np.std(CaTrip2,axis = 0)
    
    plt.plot()
    
    #VectOfInjectedRows = np.arange(SpecialTransLims[0],SpecialTransLims[1]+1) ### The special transit limits includes the last one so adding one to range since up to but not including 
    VectOfInjectedRows = np.arange(13) ### The special transit limits includes the last one so adding one to range since up to but not including 

    
    
#    plt.figure()
#    #for i in range(nrows):
#    for i in range(6,10,1):
#        plt.plot(CaTrip1[i,:],label=i)
#    plt.legend()     
#    plt.title('before scaling')
    
    sflist = []
    refpixellist = []
    
    
    
    ScalingFunctions = np.zeros((13,ncols))
    
    
    
    for InjectedRowCounter in range(13):  ## Always 4 in transit spectra         
        
        InjectedRow = VectOfInjectedRows[InjectedRowCounter]  
        
        MiddleIndex = 379
        
        ### For L1 
        
        L1OffsetOfPlanetSignal = L1IndicesForIteration[InjectedRowCounter] - CaTrip1PixCenter          
        
        L1PixelToScaleToOne = MiddleIndex + L1OffsetOfPlanetSignal
        refpixellist.append(L1PixelToScaleToOne)
        
        L1SNRScalingFactor = 1/L1SNR[L1PixelToScaleToOne]
        
        sflist.append(L1SNRScalingFactor)
        
        ScaledL1SNR = np.copy(L1SNR)
        
        ScaledL1SNR *= L1SNRScalingFactor    
        
        ScalingFunctions[InjectedRowCounter,:] = ScaledL1SNR
        
#        plt.figure()
#        xx = np.arange(ncols)
#        plt.plot(ScaledL1SNR,label=InjectedRow)
#        plt.plot([xx[0],xx[-1]],[1,1],'k--')
#        plt.plot([L1PixelToScaleToOne,L1PixelToScaleToOne],[np.min(ScaledL1SNR),np.max(ScaledL1SNR)],'k--')
#        plt.title('injected row %d'%(InjectedRow))
        
        
        
    #    np.save('scaleddump.npy',ScaledL1SNR)
    #    np.save('scaleddumppix.npy',L1PixelToScaleToOne)
        
        ScaledCaTrip1[InjectedRow,:] *= ScaledL1SNR
        
        CaTrip1 = ScaledCaTrip1
    
#    plt.figure()
#    for i in range(4):
#        plt.plot(ScaledCaTrip1[i,:])    
        
        ### For L2 
        
        L2OffsetOfPlanetSignal = L2IndicesForIteration[InjectedRowCounter] - CaTrip2PixCenter          
          
        L2PixelToScaleToOne = MiddleIndex + L2OffsetOfPlanetSignal
        
        L2SNRScalingFactor = 1/L2SNR[L2PixelToScaleToOne]
        
        ScaledL2SNR = np.copy(L2SNR)
        
        ScaledL2SNR *= L2SNRScalingFactor       
        
        CaTrip2[InjectedRow,:] *= ScaledL2SNR      
   
    
    
    ########### End new noise scaling 
    #######################     
    
    
    Weight1 = 2.0**2
    Weight2 = 1.0**2
     
#    Weight1 = 1.0
#    Weight2 = 1.0    
    
    combined = (Weight1*CaTrip1+Weight2*CaTrip2)/(Weight1+Weight2)    
    
    pyfits.writeto('%s/CaIItrip_combined.fits'%(writedir),combined,clobber=True)
    
    wcombined = w[CaTrip2PixCenter-LinePartHalfExtent:CaTrip2PixCenter+LinePartHalfExtent]
    
    np.save('%s/%s_combined.npy'%(writedir,name),wcombined)
    
    d4 = shift_to_planet_frame(combined,wcombined,jd,vorbital='calculate')
    
    combined_interval = 360
    
    combined_numcols = np.shape(combined)[1]
    
    combined_Na_centre = (combined_numcols/2.0)-1 # -1 to get a python index 
    
    pyfits.writeto('%s/not_cropped_edge_nans.fits' % (writedir),d4,clobber=True)
    
    d4 = d4[:,combined_Na_centre-combined_interval:combined_Na_centre+combined_interval]
    
    pyfits.writeto('%s/d4.fits'%(writedir),d4,clobber=True)
    
    
    phase,radv,translims,numorbits = JDtophase_radial_vel(jd,SpecialTransLims)
    
    starttrans = translims[0]
    endtrans = translims[1]
    
    planetspec = np.sum(d4[starttrans:endtrans+1,:],axis=0)/4.0  ## To include the last spectrum in SpecialTransitLimits     

     ##### Estimate SNR in the matrix by taking the stddev of a region without residuals from tellurics and stellar 
    #noisesample = np.std(d3_precolweighting[:,12925:12947],axis=0) #old interval in heliocentric frame spectra
    noisesample = np.std(d3_precolweighting[:,47900:47990],axis=0) 
    stdcol = np.median(noisesample)
    snrperspecperNaline = 1.0/stdcol  # Note these limits have been adjusted to be for Ca II but name kept same.  Na limits were d3[:,5405:5525]
    SNRaverage_over_transit_and_Na_lines = (1.0/stdcol)*np.sqrt(2.0)*np.sqrt(endtrans-starttrans)
    
    snrarray = np.array([stdcol,snrperspecperNaline,endtrans-starttrans,SNRaverage_over_transit_and_Na_lines])
    np.savetxt('%s/%s_SNR.txt'%(writedir,name),snrarray)
    print('snr',snrarray)
    
    
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
    
    InjectedPartsSTD = np.array([np.mean(residual_col_std[385:405]),np.mean(residual_col_std[415:435])])

    return planetspec,phase,radv,translims,numorbits
    
    


#################################


notbinned_snrlist = []
binned_snrlist = []


pcavect = [0]  

#width = 0.01
#width = 0.1*2
width = 0.25*1.191305315878679

#vorbital = 240



ToInjectTrueFals = True  
#InjectionStrength = 4.0e-2
#InjectionStrength = 3e-2
#InjectionStrength = 3e-2

#InjectionStrength = 1.5e-2

#InjectionStrength = 50e-2
#InjectionStrength = 18e-2
#InjectionStrength = 8e-2  ## 5 sigma limit with wide signal 

#InjectionStrength = 3.6e-2 
#InjectionStrength = 36e-2 
#InjectionStrength = 18e-2 

#InjectionStrength = 5.4e-2  ## 5.3e-2 gives an effective injection strenght of  5.1805802332657747e-2 needed for 3 sig 
#InjectionStrength = 10.8e-2  ## 10.8e-2 gives an effective injection strenght of  10.36e-2 needed for 6 sig 
#InjectionStrength = 5.2e-2

#InjectionStrength = 6e-2
#InjectionStrength = 10e-2
#InjectionStrength = 3e-2
#InjectionStrength = 2e-2
#InjectionStrength = 1.5e-2
#InjectionStrength = 12e-2

#InjectionStrength = 12e-2
#InjectionStrength = 30e-2

#InjectionStrength = 10e-2 


#InjectionStrength = 5e-2
#InjectionStrength = 6.7e-2
#InjectionStrength = 9.3e-2
#InjectionStrength = 9.7e-2
#InjectionStrength = 8.4e-2
#InjectionStrength = 5.9e-2
#InjectionStrength = 9e-2

InjectionStrength = 50e-2



EmpiricalWavelenghtShift = 0.0#-0.22
EmpiricalWavelenghtShift = 0.0#-0.5 ## required to get the planet frame transformed to have zero radial velocity 


############# end not as function  
##############################################################

#### This is all done further down by taking the standard deviation of the central 50 pixels now 
#Weights = np.array([66.033222,70.443019,68.937256,66.469294])

#Weights = np.array([29.60602097, 36.67902957, 34.73562511, 27.0476788])**2

#Weights = np.array([1.0,1.0,1.0,1.0])

BinCentralPixel = 501
#BinSizePixels = 3
BinSizePixels = 4

TransitShift = 0

SNRCalcParams  = (360,50,250)


for pcatoremove in pcavect:       
    

    name = 'night1'
    DopplerShiftForModel = -7.68873362
    SpecialTransLims = (5+TransitShift,8+TransitShift) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  

    
    fitsfile = pyfits.open('processing/Na/night1/AlignedNoCosmics/night1MedNormNoCosmics.fits')

    
    data = fitsfile[0].data
    
    full_mjd_vect = np.loadtxt('IgnasReduced/night1/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night1/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/Na/night1/Night1WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n1spec,n1phase,n1radv,n1translims,n1numorbits = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    
    
    n1specNoInject,n1phaseNoInject,n1radvNoInject,n1translimsNoInject,n1numorbitsNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    

    name = 'night2'
    DopplerShiftForModel = - 1.48105377
    SpecialTransLims = (5+TransitShift,8+TransitShift) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  


    fitsfile = pyfits.open('processing/Na/night2/AlignedNoCosmics/night2MedNormNoCosmics.fits')

    data = fitsfile[0].data
    
    full_mjd_vect = np.loadtxt('IgnasReduced/night2/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night2/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/Na/night2/Night2WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n2spec,n2phase,n2radv,n2translims,n2numorbits = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    
    
    n2specNoInject,n2phaseNoInject,n2radvNoInject,n2translimsNoInject,n2numorbitsNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    
#    
#    
#    
#
    name = 'night3'
    DopplerShiftForModel = -12.59669067 
    SpecialTransLims = (4+TransitShift,7+TransitShift) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  

    fitsfile = pyfits.open('processing/Na/night3/AlignedNoCosmics/night3MedNormNoCosmics.fits')

    data = fitsfile[0].data
    
    full_mjd_vect = np.loadtxt('IgnasReduced/night3/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night3/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/Na/night3/Night3WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n3spec,n3phase,n3radv,n3translims,n3numorbits = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    
    
    n3specNoInject,n3phaseNoInject,n3radvNoInject,n3translimsNoInject,n3numorbitsNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    
    
    
    
    name = 'night4'
    
    DopplerShiftForModel = -8.80273704
    SpecialTransLims = (6+TransitShift,9+TransitShift) ## including 8 so a +1 was added when summing over the in-transit spectra.  The model injector already injects up to and including  

    
    fitsfile = pyfits.open('processing/Na/night4/AlignedNoCosmics/night4MedNormNoCosmics.fits')


    data = fitsfile[0].data
    
    full_mjd_vect = np.loadtxt('IgnasReduced/night4/MJDsFromArchive.txt') 
    exps = np.loadtxt('IgnasReduced/night4/EXPsFromArchive.txt')

    mjd_exps = full_mjd_vect + 0.5*exps/(24*3600.)
    
    mjd = np.mean(mjd_exps.reshape(-1, 2), axis=1)  ## Averaging every two values 
  
    jd = mjd + 2400000.5      
    
    wavelength = np.loadtxt('processing/Na/night4/Night4WavelengthNoPadding.txt') - EmpiricalWavelenghtShift
    
    n4spec,n4phase,n4radv,n4translims,n4numorbits = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=ToInjectTrueFals,injectStrength=InjectionStrength,ModelRadialVelocity=DopplerShiftForModel)    

    n4specNoInject,n4phaseNoInject,n4radvNoInject,n4translimsNoInject,n4numorbitsNoInject = process(name,data,wavelength,jd,pcatoremove,width,SpecialTransLims,vorbital='calculate',ToInject=False,injectStrength=InjectionStrength)    

    deltaw = wavelength[1] - wavelength[0]

    ##################
    ### Divide by rolling average 
    
    RollingMeanWindow = 20 
    
    n1NIRM = pd.rolling_mean(pd.DataFrame(n1specNoInject),window=RollingMeanWindow).values.flatten()
    n2NIRM = pd.rolling_mean(pd.DataFrame(n2specNoInject),window=RollingMeanWindow).values.flatten()
    n3NIRM = pd.rolling_mean(pd.DataFrame(n3specNoInject),window=RollingMeanWindow).values.flatten()
    n4NIRM = pd.rolling_mean(pd.DataFrame(n4specNoInject),window=RollingMeanWindow).values.flatten()

    n1specNoInject = ((n1specNoInject + 1)/(n1NIRM + 1)) - 1 
    n1spec = ((n1spec + 1)/(n1NIRM + 1)) - 1 
    
    n2specNoInject = ((n2specNoInject + 1)/(n2NIRM + 1)) - 1 
    n2spec = ((n2spec + 1)/(n2NIRM + 1)) - 1 
    
    n3specNoInject = ((n3specNoInject + 1)/(n3NIRM + 1)) - 1 
    n3spec = ((n3spec + 1)/(n3NIRM + 1)) - 1 
    
    n4specNoInject = ((n4specNoInject + 1)/(n4NIRM + 1)) - 1 
    n4spec = ((n4spec + 1)/(n4NIRM + 1)) - 1 


    
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
    
    DataSetWeight = 1.0/(np.array([np.std(n1specBinnedNoInject[500-25:500+25]),np.std(n2specBinnedNoInject[500-25:500+25]),np.std(n3specBinnedNoInject[500-25:500+25]),np.std(n4specBinnedNoInject[500-25:500+25])])**2)



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

    avg_CaII_wavelength = np.mean(np.array([5889.95, 5895.92]))
   

    plot_wx = (np.arange(len(n4spec))-len(n4spec)/2.0)*deltaw
    #plot_wx = ((np.arange(len(n4spec))-len(n4spec)/2.0)-2)*deltaw
    plot_vx = (plot_wx/avg_CaII_wavelength)*c

    snrlist = []
    signallist = []
    noiselist = []
    
    origsnrlist = []
    
    #PlotYLims = (-0.15,0.15)
    #PlotYLims = (-0.03,0.03)
    #PlotYLims = (-0.040,0.005)
    #PlotYLims = (-0.12,0.08)    
    PlotYLims = (-0.14,0.12)   
    #PlotYLims = (-1,1)
    #PlotYLims = (-0.05,0.05)  
    #PlotYLims = (-0.03,0.03)  
    
    PlotYLimsCombined = PlotYLims#(-0.08,0.08)
    PlotXLims = (-500,500)
    #PlotXLims = (-1000,1000)
    
    ## Divide by about 6 to get A 
    
    plt.figure()
    plt.title('n1')
    plt.plot(plot_vx,n1spec,'k:')
    plt.plot(plot_vx,n1specBinned,'k-')
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
    #plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
    
    plt.figure()
    plt.title('No injected combined without night 1')
    plt.plot(plot_vx,NoInjectedCombinedWithout1,'k:')
    plt.plot(plot_vx,NoInjectedCombinedWithout1_binned,'k-')
    plt.xlim(PlotXLims)
    #plt.ylim(PlotYLims)
    plt.xlabel('systematic velocity (km/s)')
    plt.ylabel('absorption (arbitrary units)')
   
    n1SNR = calculate_SNR(n1spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n2SNR = calculate_SNR(n2spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n3SNR = calculate_SNR(n3spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n4SNR = calculate_SNR(n4spec,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)

     
    CombinedSNR = calculate_SNR(combined,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    
    #CombinedNo1SNR = calculate_SNR(combinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    CombinedNo1SNR = calculate_SNR_BetterPlot(combinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    CombinedNoInjecteNo1SNR = calculate_SNR_BetterPlot(NoInjectedCombinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    
    CombinedNoInjecteNo1SNR_binned = calculate_SNR_BetterPlot(NoInjectedCombinedWithout1_binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)

    
    ###############################################
    
    spectrum = combinedWithout1
    minindex = SNRCalcParams[0]
    signal_extent = SNRCalcParams[1]
    noise_extent = SNRCalcParams[2]
    toplot = True 
    
    
    signal = -np.min(spectrum[minindex-signal_extent:minindex+signal_extent])
    
    leftnoise = np.std(spectrum[minindex-noise_extent:minindex-signal_extent])
    
    rightnoise = np.std(spectrum[minindex+signal_extent:minindex+noise_extent])
    
    avgnoise = (leftnoise+rightnoise)/2.0
    
    
    if toplot == True:
        
        xvect = plot_vx
        plt.figure()
        plt.plot(xvect[minindex-noise_extent:minindex+noise_extent],spectrum[minindex-noise_extent:minindex+noise_extent],'b')
        plt.plot(xvect[minindex-signal_extent:minindex+signal_extent],spectrum[minindex-signal_extent:minindex+signal_extent],'r') 
        plt.xlabel('radial velocity (km/s)')           
        plt.ylabel('absorption relative to stellar spectrum')
        
        plt.figure()
        plt.plot(xvect[minindex-noise_extent:minindex+noise_extent],spectrum[minindex-noise_extent:minindex+noise_extent],'b')
        plt.plot(xvect[minindex-signal_extent:minindex+signal_extent],spectrum[minindex-signal_extent:minindex+signal_extent],'r')   
        plt.xlim(PlotXLims)
        plt.xlabel('radial velocity (km/s)')           
        plt.ylabel('absorption relative to stellar spectrum')
        
    plt.figure()
    plt.plot(plot_vx,combinedWithout1)
    plt.xlim(PlotXLims)
    plt.xlabel('radial velocity (km/s)')           
    plt.ylabel('absorption relative to stellar spectrum')
    plt.title('injected signal') 
    
    plt.figure()
    plt.plot(plot_vx,NoInjectedCombinedWithout1)
    plt.xlim(PlotXLims)
    plt.xlabel('radial velocity (km/s)')           
    plt.ylabel('absorption relative to stellar spectrum')
    plt.title('no injected signal')
    
    plt.figure()
    plt.plot(plot_vx,combinedWithout1,'r')
    plt.plot(plot_vx,NoInjectedCombinedWithout1)
    plt.xlim(PlotXLims)
    plt.xlabel('radial velocity (km/s)')           
    plt.ylabel('absorption relative to stellar spectrum')
       
    
    plt.plot(plot_vx,NoInjectedCombinedWithout1,'k:')
    plt.plot(plot_vx,combinedWithout1,'k-')
    plt.xlim(PlotXLims)
        
        
    
    ######################################################
    
    
    ## The binned SNRS 
    
    n1SNRBinned = calculate_SNR(n1specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n2SNRBinned = calculate_SNR(n2specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n3SNRBinned = calculate_SNR(n3specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    n4SNRBinned = calculate_SNR(n4specBinned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)

     
    CombinedSNRBinned = calculate_SNR(combined_binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
    
    CombinedNo1SNR = calculate_SNR_BetterPlot(combinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    CombinedNo1NoInjectedSNR = calculate_SNR_BetterPlot(NoInjectedCombinedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    
    CombinedNo1SNR_Binned = calculate_SNR_BetterPlot(combinedWithout1_binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    CombinedNo1NoInjectedSNR_Binned = calculate_SNR_BetterPlot(NoInjectedCombinedWithout1_binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)           

    
    ##### Subtracting the spectrum without an injected signal
    CombinedNo1NoInjectSubbedSNR = calculate_SNR_BetterPlot(CombinedNoInjectSubbedWithout1,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],True)
    
    ##### Signal to noise after subtracting 
    ActualNoInjectedSubbedSNR = CombinedNo1NoInjectSubbedSNR[1]/CombinedNo1SNR[2]
    
    CombinedNoInjectSubbedWithout1BinnedSNR =  calculate_SNR(CombinedNoInjectSubbedWithout1Binned,SNRCalcParams[0],SNRCalcParams[1],SNRCalcParams[2],False)
 
    #ActualBinnedNoInjectedSubbedSNR = CombinedNoInjectSubbedWithout1BinnedSNR[1]/(2*CombinedSNRBinned[2])  ## *2 noise since the signal is subtracted 
    ActualBinnedNoInjectedSubbedSNR = CombinedNoInjectSubbedWithout1BinnedSNR[1]/(CombinedSNRBinned[2])  ## *2 noise since the signal is subtracted 

 
    ### I think this combined SNR is the best measure of noise 
    print(CombinedNo1SNR)
    
yI = combinedWithout1 + 1 
yNI = NoInjectedCombinedWithout1 + 1

AbsI = 1 - np.min(yI)
AbsNoI = 1 - np.min(yNI)

RecoveredInjectionStrength = 1 - (np.min(yI)/np.min(yNI))
#RecoveredInjectionStrength = 1 - np.min(yI)/1.0

RecoveryEfficiency = RecoveredInjectionStrength/InjectionStrength

CombinedInjectedPartSTD = np.std(NoInjectedCombinedWithout1[360-4:360+4])

StdOfInjectedPartWithoutWeighting = 0.017176140762723384  ## does use the 2:1 weighting of the Na lines though.  Obtained from np.std(NoInjectedCombinedWithout1[360-4:360+4]) of this code without down weighting noisy columns and weighting the planet spec average 
  
#SNRFromNoiseOfInjectedPart = InjectionStrength/STDOfInjectionPartBeforeDownWeightingNoisyColumns

SNRFromNoiseOfInjectedPart =  5.1805802332657747e-2/StdOfInjectedPartWithoutWeighting


print(SNRFromNoiseOfInjectedPart)    

#######################################
##### Plot binned and unbinned columns 

dashparams = (4,1.5)

fig, axes = plt.subplots(nrows=4, ncols=2,figsize=(7.27,10.69))



axes[0][0].plot(plot_vx,n2spec,'k--',dashes=(dashparams[0], dashparams[1]))
axes[0][0].plot(plot_vx,n2specNoInject,'k')
axes[0][0].set_xlim(PlotXLims)
axes[0][0].set_ylim(PlotYLims)


axes[0][1].plot(plot_vx,n2specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[0][1].plot(plot_vx,n2specBinnedNoInject,'k')
axes[0][1].set_xlim(PlotXLims)
axes[0][1].set_ylim(PlotYLims)

axes[1][0].plot(plot_vx,n3spec,'k--',dashes=(dashparams[0], dashparams[1]))
axes[1][0].plot(plot_vx,n3specNoInject,'k')
axes[1][0].set_xlim(PlotXLims)
axes[1][0].set_ylim(PlotYLims)

axes[1][1].plot(plot_vx,n3specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[1][1].plot(plot_vx,n3specBinnedNoInject,'k')
axes[1][1].set_xlim(PlotXLims)
axes[1][1].set_ylim(PlotYLims)

axes[2][0].plot(plot_vx,n4spec,'k--',dashes=(dashparams[0], dashparams[1]))
axes[2][0].plot(plot_vx,n4specNoInject,'k')
axes[2][0].set_xlim(PlotXLims)
axes[2][0].set_ylim(PlotYLims)

axes[2][1].plot(plot_vx,n4specBinned,'k--',dashes=(dashparams[0], dashparams[1]))
axes[2][1].plot(plot_vx,n4specBinnedNoInject,'k')
axes[2][1].set_xlim(PlotXLims)
axes[2][1].set_ylim(PlotYLims)

axes[3][0].plot(plot_vx,combinedWithout1,'k--',dashes=(dashparams[0], dashparams[1]))
axes[3][0].plot(plot_vx,NoInjectedCombinedWithout1,'k')
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

    

    


################### 
#### Without the combined night 1 

#fig, axes = plt.subplots(nrows=5, ncols=2,figsize=(7.27,10.69))
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
#axes[4][0].plot(plot_vx,NoInjectedCombinedWithout1,'k')
#axes[4][0].plot(plot_vx,combinedWithout1,'k:')
#axes[4][0].set_xlim(PlotXLims)
#axes[4][0].set_ylim(PlotYLims)
#
#
#axes[4][1].plot(plot_vx,NoInjectedCombinedWithout1_binned,'k')
#axes[4][1].plot(plot_vx,combinedWithout1_binned,'k:')
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

#################################################

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
##axes[0][0].set_title('Unbinned \n Night 1')
##axes[0][1].set_title('Binned \n Night 1')
##
##for i in range(1,4):
##
##    axes[i][0].set_title('Night %d'%(i+1))
##    axes[i][1].set_title('Night %d'%(i+1))
##
##axes[4][0].set_title('Weighted mean of all nights')  
###axes[4][1].set_title('mean of all nights')
##
##axes[5][0].set_title('Weighted mean excluding Night 1')  
###axes[5][1].set_title('mean excluding Night 1')
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
    
