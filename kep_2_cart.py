#!/usr/bin/python3

import numpy as np
import sys


mu=0.9
t=0
T=0
EA=0
nu=0
a=1
omega_LAN = 0
omega_AP = 0
i=0

def kep_2_cart(e):
   
   
    n = np.sqrt(mu/(a**3))
    M = n*(t - T)
    #2
    MA = EA - e*np.sin(EA)
    #3
    #
    nu = 2*np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(EA/2))
    #4
    r = a*(1 - e*np.cos(EA))
    #5
    
    h = np.sqrt(mu*a * (1 - e**2))
    #6
    Om = omega_LAN
    w =  omega_AP
    
    X = r*(np.cos(Om)*np.cos(w+nu) - np.sin(Om)*np.sin(w+nu)*np.cos(i))
    Y = r*(np.sin(Om)*np.cos(w+nu) + np.cos(Om)*np.sin(w+nu)*np.cos(i))
    
    #7
    p = a*(1-e**2)
    
    V_X = (X*h*e/(r*p))*np.sin(nu) - (h/r)*(np.cos(Om)*np.sin(w+nu) + np.sin(Om)*np.cos(w+nu)*np.cos(i))
    V_Y = (Y*h*e/(r*p))*np.sin(nu) - (h/r)*(np.sin(Om)*np.sin(w+nu) - np.cos(Om)*np.cos(w+nu)*np.cos(i))
    return [X,Y],[V_X,V_Y]


if len(sys.argv)!=2:
    print("Wrong number of arguments! Exiting")
    exit()



e = float(sys.argv[1])
[x,y], [vx,vy] = kep_2_cart(e)
print(e,x,vy)
