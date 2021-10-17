# -*- coding: utf-8 -*-
def makespec(species,flux=0.01,sigma='instr_res'):

    """
    Created on Wed Feb 25 18:28:20 2015
    
    @author: arh
    
    This script loads a .npy file of the lines for Fe that were obtained from the VALD 
    spectral line database, convolves the delta function lines with a Gaussian of width 
    equal to the resolution of the spectrograph as obtained by measuring the width 
    of a telluric line in order u2 by fitting a Gaussian and using that value of sigma.
    It then saves the resulting spectrum to be used as a model spectrum in the 
    correlation.
    
    The .npy file that this script takes as input was generated using a 
    different file which did the convoluted reading the VALD text file.
    The VALD file reading script is VALD_spec_maker.py   
    
    Need about 150 km/s to shift both ways    
    
    takes as default parameters:
    
    cdelt = 0.0168819795611807    # from u2 in A/pixel
    pixsigma = 1.8378453090531477  #pixel (found by measuing the width of a telluric line in the u2 order)
    
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from astropy.convolution import Gaussian1DKernel, convolve
    import os
    
    
    path = 'model_spectra/%s/' % (species) #the path to the Fe model 
    
    fname = path+'ToyCaIITriplet.npy' 
    
    spec_matrix = np.load(fname)
    
    warray = spec_matrix[1,:]
    relarray = flux*spec_matrix[0,:]/max(spec_matrix[0,:]) # the relative intensity    
        
    ## The default cdelt and pixsigma for calculating the resolution width 
        
    cdelt = 0.0168819795611807
    pixsigma = 1.8378453090531477
    
    if sigma == 'instr_res':
        wsigma = cdelt*pixsigma #width in A 
    else:
        wsigma = sigma 
        
    zeroedgepadding = 5.0 # in A, the amount to set as zeros on either side of the wavelengths
    
    
    # wavelength lims are 5830 - 5911.01
    

    leftwlim = 8.49801800e+03 - 1000
    rightwlim = 8.66214000e+03 + 1000
        
            
    Fe_wave = warray
    Fe_rel = relarray
    
    Ninterp = ((rightwlim+zeroedgepadding)-(leftwlim-zeroedgepadding))/(wsigma/32.0)
    
    ww = np.linspace(leftwlim-zeroedgepadding,rightwlim+zeroedgepadding,num=Ninterp)
    yy =np.zeros_like(ww)
        
        #######################
        
    print('making the tuple list')
    
    tuplist = []
    for i in range(len(ww)):
        
        tuplist.append((ww[i],yy[i]))
        
    for i in range(len(Fe_wave)):
        
        tuplist.append((Fe_wave[i],Fe_rel[i]))
    
    print('sorting the tuple list')
    
    tuplist.sort()
    
    fullwarray = np.empty((len(tuplist)))
    fullyarray = np.empty((len(tuplist)))
    
    #taking the tuples out of the sorted tuplist and putting them into vectors
    
    print('making the vectors')
    
    for i in range(len(tuplist)):
        
        fullwarray[i]=tuplist[i][0]
        fullyarray[i]=tuplist[i][1]  
        
    localfluxmax = max(fullyarray) # rescale the flux after convolution with the gaussian 

    
    templatecdelt = np.median(fullwarray[1:]-fullwarray[0:-1])
    
    tempatepixsigma = wsigma/templatecdelt 
    print('wsigma',wsigma) 
    
    tempatepixsigma = np.round(tempatepixsigma)
    print('template pix sigma', tempatepixsigma)
    
#    if tempatepixsigma % 2.0 != 0.0:
#        tempatepixsigma = tempatepixsigma + 1.0        
        
   
    g = Gaussian1DKernel(stddev=tempatepixsigma)
   
    # Convolve data
    z = convolve(fullyarray, g, boundary='extend')  
    
    rescaling_factor = max(fullyarray)/localfluxmax
    
    z=z*rescaling_factor
    
    results_array = np.empty((2,len(z)))
    
    results_array[0,:]=z*(flux/np.max(z))
    results_array[1,:]=fullwarray
    print(max(z))

    savedir = '%s/ToyCaIITriplet_F%.2E_s%f/' % (path,flux,wsigma)    
#    
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    
    np.save('%s/ToyCaIITriplet.npy' % (savedir),results_array)    
              ###/ (path ends in a /)
   
    plt.figure()
    #plt.plot(fullwarray,fullyarray,'b')
    #plt.plot(fullwarray,z,'r')
    plt.plot(results_array[1,:],results_array[0,:],'.')
    #plt.plot(fullwarray[54919:54919+257],g.array/max(g.array),'k')

    #plt.figure()
    #plt.plot(g.array/max(g.array))

    return results_array  

#makespec(species='VALD_Na',flux=0.010)
#q = makespec(species='VALD_CaII',flux=0.01,sigma=0.09261)     
