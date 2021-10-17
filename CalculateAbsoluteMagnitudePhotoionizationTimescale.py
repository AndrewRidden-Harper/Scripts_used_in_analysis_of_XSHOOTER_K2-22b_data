# -*- coding: utf-8 -*-
"""
Created on Sun Jul  8 21:14:51 2018

@author: andrew
"""

import numpy as np

## K2-22 V bad 
#m = 16.44  # K2-22 apparent visual 
#Msun = 4.81  # Sun absolute visual ## taken from Wilmer 2018  Johnson_V


## U band 
## K2-22 U band 
m = 19.07  # K2-22 apparent u  
Msun = 5.61  # Sun absolute u ## taken from Wilmer 2018  Johnson_U

### For CoRoT-7b
#m = 11.668
d = 225  ### For K2-22 b
#d = 520  ## For CoRoT-7b 

M = m - 5*(np.log10(d)-1)
#
dM = Msun - M
#dM = -3.5 ## Ignas estiamte 

#
#print dM
#
#

#
F2overF1 = 10**(0.4*dM)
#F2overF1 = 1


d1 = 1       ### for 1 au to compare to the timescales of comets at 1 au in the solar system
d2 = 0.0088  ## For K2-22 b
#d2 = 0.01  ## Ignas estimate for K2-22 b 

#d2 = 0.0172  ## for CoRoT-7


dfluxratio = (d1/d2)**2


OverallFluxIncreaseFactor = dfluxratio*F2overF1

print('Overall UV Flux Increase Factor for K2-22 b (relative to Solar system at 1 au)')
print(OverallFluxIncreaseFactor)
print() 

#NaT1au = 1.69e5  #s 
NaTFulle = 1.9e5

CoRoT7NaTime = 56
CoRoT7CaITime = 1e3

K2NaTime = NaTFulle/OverallFluxIncreaseFactor


print('Scaled Na photoionisation time')
print(K2NaTime)

RatioK2toCorot = K2NaTime/CoRoT7NaTime

K222CaITime = RatioK2toCorot*CoRoT7CaITime

print('The Na ionisation timescale of K2-22 is %f times longer than for CoRoT-7 (due to its low UV flux)'%(RatioK2toCorot))
print() 
print('Assuming this same scaling holds for Ca+, the timescale of for K2-22 is %.5e'%(K222CaITime)) 
