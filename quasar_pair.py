import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import chisquare
from astropy.io.fits import getdata
from astropy.stats import sigma_clip
from astropy.io import fits
from scipy.signal import argrelextrema
from scipy import ndimage
from PyAstronomy import pyasl
from astropy import constants as const
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u

#cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
data, header = getdata("/home/priyankajalan/Suv/dr14q_spec_prop.fits",1,header=True) #sim_image 

idx = np.where(data['REDSHIFT'] <= 0.5)[0]
ra = data['RA'][idx]
dec = data['DEC'][idx]
z = data['REDSHIFT'][idx]


f = open('pairs.dat','w')
pair = ([-1,-1])
for i in range(len(idx)):
    print('-------- '+str(i) + ' source out of '+str(len(idx)))
    ang = pyasl.getAngDist(ra[i],dec[i],ra,dec)
    idz = np.where((ang >0) & (ang <= 10/60.))[0] # 10 arcmin
    for c in range(len(idz)):
        j = idz[c]
        delz = abs(z[i]-z[j])
        H_z = cosmo.H(z[i]) 
        r_parr  = const.c.to('km/s')*delz/((1+z[i])*H_z)
        d_A     = cosmo.angular_diameter_distance(z[i])
        r_perp  = d_A*(ang[j]*np.pi/180.)
        r       = np.sqrt((r_parr**2+r_perp**2))
        if (r < 10*u.Mpc): # 10 Mpc #(ang <= 50/3600.): # & (const.c*abs((z[i]-z[j])/(1.+z[i])) <= 500000):
            
            f.write(str(data['SDSS_NAME'][idx][i]) + '\t' + str(round(ra[i],5)) + '\t' + str(round(dec[i],5)) + '\t' + str(round(z[i],5)) +'\t' + str(data['SDSS_NAME'][idx][j]) + '\t' + str(round(ra[j],5)) + '\t' + str(round(dec[j],5)) + '\t' + str(round(z[j],5)) + '\t' + str(round(ang[j]*3600,5))+ '\t' + str(round(r.value,5))+ '\n')
            '''
            for k in range(len(pair)):
            if [j,i] != pair[k]: # removing duplicate pairs
                pair.append([i,j])
            '''
f.close()
