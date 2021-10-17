# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 20:46:26 2018

@author: andrew
"""

import numpy as np

D1 = 0.387098 ## AU
L1 = 1        ## luminosity in solar units

D2 = 0.0088 ## AU 
L2 = 0.063  ## luminosity in solar units 


ratio = (L1/L2)*(D2/D1)**2  #P1/P2 


## So K2-22 b has 1/ratio times more energy density 

1/ratio

### It would increase the velocity proportionally to sqrt(1/ratio)

np.sqrt(1/ratio)

### Killen 2002 gives an average Na velocity on Mercury of 7 km/s so that becomes (77 km/s):

np.sqrt(1/ratio)*7