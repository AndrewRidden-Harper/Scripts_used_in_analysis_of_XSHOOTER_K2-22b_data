# -*- coding: utf-8 -*-
def pcasub(data,components_to_remove=4):


    """
    Created on Wed Feb  4 20:52:49 2015
    
    @author: andrew
    
   
    """
    import numpy as np
    import astropy.io.fits as pyfits
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import time
    
       
    starttime = time.clock()
    
    plt.close('all')
    
    hdu=pyfits.PrimaryHDU(np.array([]))
    
    hdulist_removed=pyfits.HDUList([hdu])    
    hdulist_data = pyfits.HDUList([hdu])  
    

    ns = components_to_remove
    
   
    if ns>0:
     

        d0 = np.copy(data)

        
        nx = np.shape(d0)[1]
        ny = np.shape(d0)[0]
        
        t=np.array((nx))
        
        
        d0=d0-np.mean(d0,axis=0)
    
        removed = np.empty((ny,nx,ns))    
        removed[:]=np.NAN
        
        d0subbed = np.empty((ny,nx,ns))    
        d0subbed[:]=np.NAN
        
        for j in range(1,ns+1):
            
            p=np.random.rand(nx)
            
            for i in range(1,31):
                
                
                a = d0*p
                b=np.sum(a,axis=1)
                
                c=np.empty((nx,ny))
                for z in range(nx):
                    c[z,:]=b
                    
                d = c.T*d0
                
                t=np.sum(d,axis=0)
                
                p=t/np.sqrt(np.sum(t**2))
                
            d0subbed[:,:,j-1] = d0 - removed[:,:,j-1]
            
            for i in range(ny):
                tosubtract = np.sum(d0[i,:]*p)*p
                removed[i,:,j-1] =  tosubtract
                
                
                d0[i,:]=d0[i,:] - tosubtract
                
            newhdu_removed=pyfits.ImageHDU(removed[:,:,j-1])        
            hdulist_removed.append(newhdu_removed)
        
            newhdu_data = pyfits.ImageHDU(d0subbed[:,:,j-1])
            hdulist_data.append(newhdu_data)  
            
        ## To save the componets that are removed 
        #hdulist_removed.writeto('UVES_4pca2.fits',clobber=True)
            

                
    endtime = time.clock()
    
    
 
#        
    return d0 #so it's still cropped 

       

            
            