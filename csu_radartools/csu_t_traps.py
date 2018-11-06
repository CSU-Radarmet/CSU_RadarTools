"""
csu_t_traps.py

Brenda Dolan, CSU, October 2018
bdolan@atmos.colostate.edu

Creating Trapazoidal functions for the temperature fuzzy logic.
In this sense we can penalize certain temperatures for certain HID types, and 
therefore reduce some of the weighting on temperature while also reducing the harsh
boundaries between some of the species that are directly related to temperature.
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import numpy as np


def get_hid_traps(c,x_arr):
    if c == 0 :
        #Drizzle
        x1 = -3
        x2 = -1
        beta = t_trap_dz(x_arr,x1,x2)
    elif c==1:
        #Rain
        x1 = -10
        x2 = -2
        beta = t_trap_rain(x_arr,x1,x2)
    elif c==2:
        #Ice Crystals
        x1 = 0
        x2 = 5
        beta = t_trap_ic(x_arr,x1,x2)
    elif c ==3:
        #Aggregates
        x1 = -35
        x2 = -25
        x3 = 0
        x4 = 5
        beta = t_trap_agg(x_arr,x1,x2,x3,x4)
    elif c==4:
        #Wet Snow
        x1 = -3
        x2 = 3
        beta = t_trap_ws(x_arr,x1,x2)
    elif c==5:
        #Vertical ICe. Should exactly match Ice Crystals
        x1 = 0
        x2 = 5
        beta = t_trap_ic(x_arr,x1,x2)
    elif c ==6:
        #Low Density Graupel
        x1 = -40
        x2 = -30
        x3 = 0
        x4 = 5
        beta = t_trap_ldg(x_arr,x1,x2,x3,x4)
    elif c == 7:
        #High Density Graupel
        x1 = -30
        x2 = -20
        x3 = 5
        x4 = 30
        beta = t_trap_hdg(x_arr,x1,x2,x3,x4)
    elif c == 8:
        #Hail
        x1 = -30
        x2 = -20
        x3 = 5
        x4 = 15
        beta = t_trap_ha(x_arr,x1,x2,x3,x4)
    elif c == 9:
        #Big Drops
        x1 = -2
        x2 = 1
        beta = t_trap_bd(x_arr,x1,x2)
    else:
        print('HID type is out of range. Returning')
        return
    return beta

def t_trap_agg(x_arr,x1,x2,x3,x4):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = 0.5
    wh2 = np.where(np.logical_and(x_arr>=x2,x_arr<=x3))
    beta[wh2] = 1
    wh3 = np.where(x_arr >= x4)
    beta[wh3] = -1
    
    slope_x1_x2 = (.5)/np.abs(np.abs(x1)-np.abs(x2))
    slope_x3_x4 = (1.)/np.abs(np.abs(x4)-np.abs(x3))
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = slope_x1_x2*(np.abs(x1)+x_arr)+0.5
    beta[wh4] = test1[wh4]
    wh5 = np.where(np.logical_and(x_arr >=x3,x_arr<=x4))
    test2 = slope_x3_x4*(np.abs(x4)-x_arr)
    beta[wh5] =test2[wh5]
    return beta

def t_trap_bd(x_arr,x1,x2):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = -1
    wh2 = np.where(x_arr >=x2)
    beta[wh2] = 1
    
    slope_x1_x2 = (2.)/np.abs(np.abs(x1)+np.abs(x2))
#    print(slope_x1_x2)
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = np.abs(slope_x1_x2*(x_arr-x1)/(x2-x1))
    beta[wh4] = test1[wh4]
    return beta

def t_trap_dz(x_arr,x1,x2):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = -1
    wh2 = np.where(x_arr >=x2)
    beta[wh2] = 1
    
    slope_x1_x2 = (2.)/np.abs(np.abs(x1)+np.abs(x2))
    print(slope_x1_x2)
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = np.abs(slope_x1_x2*(x_arr-x1)/(x2-x1))
    beta[wh4] = test1[wh4]
    return beta
    
def t_trap_ha(x_arr,x1,x2,x3,x4):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = 0.2
    wh2 = np.where(np.logical_and(x_arr>=x2,x_arr<=x3))
    beta[wh2] = 0.8
    wh3 = np.where(x_arr >= x4)
    beta[wh3] = 0.2
    
    slope_x1_x2 = (0.6)/np.abs(np.abs(x1)-np.abs(x2))
    slope_x3_x4 = (0.6)/np.abs(np.abs(x4)-np.abs(x3))
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = slope_x1_x2*(np.abs(x1)+x_arr)+0.2
    beta[wh4] = test1[wh4]
    wh5 = np.where(np.logical_and(x_arr >=x3,x_arr<=x4))
    test2 = slope_x3_x4*(np.abs(x4)-x_arr)+0.2
    beta[wh5] =test2[wh5]
    return beta    
    
def t_trap_hdg(x_arr,x1,x2,x3,x4):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = 0.3
    wh2 = np.where(np.logical_and(x_arr>=x2,x_arr<=x3))
    beta[wh2] = 1
    wh3 = np.where(x_arr >= x4)
    beta[wh3] = -1
    
    slope_x1_x2 = (0.7)/np.abs(np.abs(x1)-np.abs(x2))
    slope_x3_x4 = (2)/np.abs(np.abs(x4)-np.abs(x3))
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = slope_x1_x2*(np.abs(x1)+x_arr)+0.3
    beta[wh4] = test1[wh4]
    wh5 = np.where(np.logical_and(x_arr >=x3,x_arr<=x4))
    test2 = slope_x3_x4*(np.abs(x4)-x_arr)-1
    beta[wh5] =test2[wh5]
    return beta 
    
    
def t_trap_ic(x_arr,x1,x2):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = 1
    wh3 = np.where(x_arr >= x2)
    beta[wh3] = -1
    
    slope_x1_x2 = (1.)/np.abs(x2+x1)
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = 1-np.abs(slope_x1_x2*(x_arr-x1))
    beta[wh4] = test1[wh4]

    return beta           
    
def t_trap_ldg(x_arr,x1,x2,x3,x4):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = 0.5
    wh2 = np.where(np.logical_and(x_arr>=x2,x_arr<=x3))
    beta[wh2] = 1
    wh3 = np.where(x_arr >= x4)
    beta[wh3] = -1
    
    slope_x1_x2 = (0.5)/np.abs(np.abs(x1)-np.abs(x2))
    slope_x3_x4 = (1)/np.abs(np.abs(x4)-np.abs(x3))
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = slope_x1_x2*(np.abs(x1)+x_arr)+0.5
    beta[wh4] = test1[wh4]
    wh5 = np.where(np.logical_and(x_arr >=x3,x_arr<=x4))
    test2 = slope_x3_x4*(np.abs(x4)-x_arr)
    beta[wh5] =test2[wh5]
    return beta 
   
def t_trap_rain(x_arr,x1,x2):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = -1
    wh3 = np.where(x_arr >= x2)
    beta[wh3] = 1
    
    slope_x1_x2 = (1.)/np.abs(x2+x1)
    wh4 = np.where(np.logical_and(x_arr >= x1,x_arr<=x2))
    test1 = np.abs(slope_x1_x2*(x_arr-x1))
    beta[wh4] = test1[wh4]

    return beta           

def t_trap_ws(x_arr,x1,x2):
    beta = np.zeros_like(x_arr)
    wh1 = np.where(x_arr <= x1)
    beta[wh1] = -1
    wh3 = np.where(x_arr >= x2)
    beta[wh3] = -1
    wh2 = np.where(np.logical_and(x_arr>=x1,x_arr<=x2))
    beta[wh2] = 1
    return beta 
