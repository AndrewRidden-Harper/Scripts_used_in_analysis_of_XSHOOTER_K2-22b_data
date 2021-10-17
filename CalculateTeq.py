# -*- coding: utf-8 -*-
"""
Created on Thu May 17 14:11:41 2018

@author: andrew
"""

import numpy as np

R_sun_m = 6.957e8 # meters
M_sun_kg = 1.988435e30# kg (kilograms

AU_in_m = 1.496e11 # meters


def calculate_Teq(Tstar,A,Rstar,D):
    
    Teq = Tstar*((1-A)**(1.0/4.0))*(((Rstar/(2*D)))**(1.0/2.0))
    
    return Teq 
    
Tstar = 3830.0

Rstar = 0.57*R_sun_m 
Mstar = 0.60*M_sun_kg
D = 0.0088*AU_in_m  ## AU
A = 0.0 


Teq = calculate_Teq(Tstar,A,Rstar,D)

print('Teq', Teq)