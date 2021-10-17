# -*- coding: utf-8 -*-
"""
Created on Sun May  6 16:19:01 2018

@author: andrew
"""

#def JDtophase_radial_vel(jd_vect,vorbital='calculate'):
#    
#    
#    
#    '''
#    Returns the planet radial velocities and transit limits based on the 
#    published data of period, mid-transit time and observation times.
#    The default orbital velocity is the one calculated from orbital period 
#    and the assumed circular orbital semi-major axis based on a calculation 
#    of semi-major axis based on Kepler's law 
#    '''
#
#    #Get the times of observations (in MJD)
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt

from  PyAstronomy import pyasl 

from astropy.time import Time

from astropy import time, coordinates as coord, units as u

import pytz

#RA = 11 h 17 min 55.856  
## 1 hour = 15 degrees 

RA = 15.0*11.0 + (17.0/60.0)*15.0 + (55.856/3600.)*15.0 



dec = 2.0 + (37.0/60.0) + (06.79/3600.)

K2_22_coords = coord.SkyCoord("11:17:55.856", "+02:37:06.79",unit=(u.hourangle, u.deg), frame='icrs')


vlt_coords = coord.EarthLocation.of_site('cerro paranal')



#mjd_vect = np.loadtxt('18March2017mjd.txt')
#mjd_vect = np.loadtxt('3April2017mjd.txt')

#SciMjd = np.array([57831.083729,
#57831.088236,
#57831.092681,
#57831.097208,
#57831.101646,
#57831.106166,
#57831.110604,
#57831.115115])





#mjd_vect = SciMjd

#mjd_vect = np.loadtxt('ScienceObMJDs/19March2017VIS.txt')
#mjd_vect = np.loadtxt('ScienceObMJDs/3April2017NIR.txt')
#mjd_vect = np.loadtxt('ScienceObMJDs/3April2017VIS.txt')
#mjd_vect = np.loadtxt('ScienceObMJDs/3April2017UVB.txt')

#mjd_vect = np.loadtxt('ScienceObMJDs/10March2018VIS.txt')
#mjd_vect = np.loadtxt('ScienceObMJDs/18March2018VIS.txt')

#mjd_vect = np.loadtxt('ScienceObMJDs/18March2017_2.txt')
#mjd_vect_full = np.loadtxt('IgnasReduced/night1/MJDsFromArchive.txt')

nightnum = 1

mjd_vect_full = np.loadtxt('IgnasReduced/night%d/MJDsFromArchive.txt'%(nightnum))
exps = np.loadtxt('IgnasReduced/night%d/EXPsFromArchive.txt'%(nightnum))




#np.savetxt('VradsAtStartofObsN4.txt',v_rad)
#np.savetxt('VradsAfterExposureN4.txt',v_rad)

mjd_exps = mjd_vect_full #+ exps/(24*3600.)


mjd_vect = mjd_vect_full#np.mean(mjd_exps.reshape(-1, 2), axis=1)

#mjd_vect = mjd_vect_full
#mjd_vect = mjd_exps


vorbital = 'calculate'

#mjd_vect = np.array([ 2460240.7275305,  2460240.8228   ,  2460240.7275305])


####
# In this script, phase 0=1 is defined to be the mid-transit point 
####    

jd_vect = mjd_vect + 2400000.5  
jd_vect_full = mjd_vect_full + 2400000.5  
deltajd_s = np.mean(jd_vect_full[1:] - jd_vect_full[0:-1])*(24.0*3600.0)


jdvect_as_times = Time(jd_vect_full,format='jd',scale='utc',location=vlt_coords) ## time at barycentre 
llt_bary = jdvect_as_times.light_travel_time(K2_22_coords)

time_at_observatory = jdvect_as_times - llt_bary 


#
#local_time = pytz.utc.localize(midtransittime.datetime).astimezone(pytz.timezone('Chile/Continental'))
#utc_time_observatory = pytz.utc.localize(jdvect_as_times.datetime).astimezone(pytz.timezone('UTC'))


Rsun = 695508 # in km 
                  
sinday = 86164.09054 #exact number of seconds in a sidereal day from WolframAlpha
G=6.67384E-11 #in m^3kg^-1 s^-2 from wiki
Msun = 1.988435E30 # kg  From WolframAlpha

Rstar = 0.57*Rsun

M55Cnc = 0.6*Msun  ## Actually the mass of K2-22 
   
Perror = 0.000001 


T0=2456811.1208 + Perror
P = 0.381078   

#TestNumOrbs = 9000.0 - 0.4
#jd_vect = np.array([T0+(TestNumOrbs*P-0.25*P),T0+(TestNumOrbs*P),T0+(TestNumOrbs*P+0.25*P)]) mjd_vect = jd_vect - 2400000.5  
   
Psec = P*sinday   

inclination=np.radians(79.52)

#calculating a with Keplers third law P^2=a^3
#a_AU=(P/365.256360417)**(2.0/3.0) #dividing P by number of days in sidereal year

#a_km=a_AU*149597870.700 #big number is 1 AU in km

a_metres = ((G*M55Cnc/(4*(np.pi**2)))*(Psec**2))**(1.0/3.0)
a_km = a_metres/1000.0

#assuming circular orbit
#v_orb=2*np.pi*a_km/(P*86400.0) #converting period of days to seconds
v_orb=2*np.pi*a_km/(Psec)

if vorbital != 'calculate':
    v_orb = vorbital 

#b=a*np.cos(inclination)

Ttot=48.0/(60.0*24.0)  # 48 minutes 
#b_trans = (a*cos(i)/R_star)*((1-e^2)/(1+e^2))
#so with e=0 (circular orbit),
#This value for b does not agree with G12 values

#Times at mid transit should give zero remainder
remainder=np.mod(jd_vect-T0,P)

#index of minimum value of remainder
#(ie, the index in this data that is closest to an exact mid-transit point,
# t=N*P+T0 
minindex = np.argmin(remainder)
print(minindex)
#The number of orbits between T0 and this data's time closest to mid-transit
#(found by taking the floor of the fractional number of orbits until the actual
#time closest to the mid-transit time) so that it is an integer number of orbits
#to reference to

num_orbits=np.floor((jd_vect[minindex]-T0)/P) 
#num_orbits = 2678

transit_start=T0+(num_orbits*P)-(Ttot/2) 
transit_end=T0+(num_orbits*P)+(Ttot/2)

midtransittime = T0+(num_orbits*P)


StartMidEnd = np.array([transit_start,midtransittime,transit_end])


jdvect_as_timesStartMidEnd = Time(StartMidEnd,format='jd',scale='utc',location=vlt_coords) ## time at barycentre 
llt_baryStartMidEnd = jdvect_as_timesStartMidEnd.light_travel_time(K2_22_coords)

time_at_observatoryStartMidEnd = jdvect_as_timesStartMidEnd - llt_baryStartMidEnd 


#np.savetxt('n%d_TransitTiming.txt'%(nightnum),np.array([transit_start,midtransittime,transit_end])-2400000.5)

print('midtransit time', midtransittime)

StartClosestIndex = np.argmin(np.abs(jd_vect - transit_start))
EndClosestIndex = np.argmin(np.abs(jd_vect - transit_end))

if jd_vect[StartClosestIndex] > transit_start:
    #### This exposure starts just after the transit starts 

    print() 
    print('at start: exposures starts just after the transit starts')
    
    #if (jd_vect[StartClosestIndex-1] + exps/(24*3600))  > transit_start:
        
    StartTransitLimit = StartClosestIndex 
    
if jd_vect[StartClosestIndex] < transit_start:
    #### this exposure starts just before the transit starts 
    print('at start: the exposure starts just before the transit starts')

    if (transit_start - jd_vect[StartClosestIndex]) < (exps[StartClosestIndex]/(24*3600)):

        StartTransitLimit = StartClosestIndex
        print('but it is still exposing as the transit starts')
    
    else:
        ## The transit starts during the dead time between the end of expsoure and the start of the next one 
        StartTransitLimit = StartClosestIndex + 1 
        print('and is during the dead time after the exposure has finished')
        
print() 
    
################    
    
if jd_vect[EndClosestIndex] > transit_end:
    ## this exposure starts just after the transit ends 
    print('at end: the exposure starts just after the transit ends')

    EndTransitLimit = EndClosestIndex -1 

if jd_vect[EndClosestIndex] < transit_end:
    #### this exposure starts just before the transit ends 

    print('at end: the exposure starts just before the transit ends')

    EndTransitLimit = EndClosestIndex
    
print() 










#    if transit_start<T0+(1261*P):
#        print 'Closest observation to mid-transit occurs after the exact transit mid point'
#    if transit_start>T0+(1261*P):
#        'Closest observation to mid-transit occurs before the exact transit mid point so please check if it works properly'
#        raise Exception

#    phase=np.empty(len(jd_vect))
#    phase[:]=np.NAN
#    
#    for i in range(len(jd_vect)):
#        if i>=minindex:
#            phase[i]=(jd_vect[i]-(P*num_orbits+T0))/P
#    
#        if i<minindex:  
#            phase[i]=(jd_vect[i]-(P*(num_orbits-1)+T0))/P
#            #calculate the phase from the orbit which finishes during these observations 
        
##Can also just define phase as:
phase=((jd_vect - (T0+num_orbits*P))/P)%1.0

transitlims=np.where(((transit_start<=jd_vect)&(transit_end>=jd_vect)))       
        
print(phase)

avgphase = np.mean(phase.reshape(-1, 2), axis=1)

np.savetxt('AvgphaseN%d.txt'%(nightnum),avgphase)


phase_angle = 2*np.pi*phase #convert the phase to an angle in radians

#This radial velocity may very well be the negative of what is calculated here,
#depedning on the direction that the planet orbits (clockwise or anticlockwise)
v_rad=(-1.0)*v_orb*np.sin(inclination)*np.sin(phase_angle)
EverySecond = v_rad[::2]

np.mean(np.diff(EverySecond))




results_array=np.empty((len(jd_vect),3))
results_array[:,0]=mjd_vect
results_array[:,1]=phase
results_array[:,2]=v_rad
##Transit limits are array indices

print(v_orb)
    #return phase,v_rad,transitlims,num_orbits 

print('phase angle (rad)', phase_angle)
print('radial v',v_rad)
print() 
print() 
    
print() 
print() 
print(transitlims)
print('alternative start lim:', StartTransitLimit)
print('alternative end lim:', EndTransitLimit)

print() 
print('spectra in transit', EndTransitLimit - StartTransitLimit +1)

baryvel_list = []

for i in range(len(jd_vect)):
    heliovel,baryvel = pyasl.baryCorr(jd_vect[i],RA,dec)    

    baryvel_list.append(baryvel)

print() 
print('barycentric velocity = %f' % (np.mean(baryvel_list)))

    
    