# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

fig= plt.figure()

#read in necessary information
lcnd=np.genfromtxt("loading_condition.txt",unpack=True)

#data calculated using fortran
fstrain,fstress=np.genfromtxt("af_stress_strain.txt", unpack=True)

#initialise values
j=1
k=0
epsp=0
xsum=0
xmod=69000
sigy=200
depsp=0
max_it=10000000
toler = 1e-8
epsp_new=0
sig=0



a=np.array([7500,1000,30])
c=np.array([310,300,4])

#step size
steps=100

#initialise arrays
estrain=np.zeros(len(lcnd))
nu=np.zeros(len(lcnd))
dsig=np.zeros(len(lcnd))
inc=[0 for i in range(10000)]
#dxback=np.zeros(3)
#xback=np.zeros(3)
xbackold=np.zeros(3)
lam=0
enew=np.zeros(10000)
snew=np.zeros(10000)

etotal=np.zeros(10000)
#sig=np.zeros(10000)
el=np.zeros(10000)



estrain[0] = lcnd[1]

#xback=np.empty([3,10000])
#dxback=np.empty([3,10000])

xbackprev=0
plas=0

#inc[0]=0

"nonlinear backstress calculation"

def back(a,c,xbackprev,plas,epsp):
    xback=a*epsp -c*epsp*xbackprev
    
    return xback
    
"derivative of backstress w.r.t plastic strain "   
def dback(a,c,xbackprev,plas):

        dxback=a - c*xbackprev
        return dxback

#newton-raphson method
def newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback):


    for n in range(0,max_it):
        
        sig=xmod*(el-epsp)
      
       
        
        xsum=np.sum(xback)     
    #check to see if the point remains in the yield surface
        func=sig-xsum-sigy
        if(abs(func)<toler):
            return epsp
        else:

            dxback=dback(a,c,xbackprev,plas)
            dxsum=np.sum(dxback)
            dfunc = nu*(-xmod - dxsum)
            depsp=-func/dfunc
            epsp=epsp+depsp
            xback=back(a,c,xbackprev,plas,epsp)
            
    return epsp
        

#calculate the increment in loading
#for i in range(1,len(lcnd)):
for i in range(1,1):
        
    if i+1 < len(lcnd):
        estrain[j]=lcnd[i+1]-lcnd[i]
        j +=1

#removing end zeros
estrain=np.trim_zeros(estrain,'b')

#create a loop to loop through the increments of strain
#for i in range(0,len(estrain)):
for i in range(0,len(estrain)):
    
    "calculate the current increment in stress from increment of strain"
    dsig=xmod*estrain[i]

    "calculate the backstress"
    xback = back(a,c,xbackprev,plas,epsp)

    "loading direction provided by sign of the strain increment"
    nu=np.sign(estrain[i])

#NEED TO CALCULATE THE BACKSTRESS HERE

    " now we know need to check with the current increment in stress is greater"
    "than the yield taking into account the shift as advised by the backstress"
    "eg plasticity will occur when sigma_old+dsig > sig_yield + backstress"

    "total backstress"
    xsum=np.sum(xback)
    
    "check to see if there is plasticity"
    lam = (nu*sigy + xsum - sig)/dsig

    #if lam> 1 then the increment is elastic, therefore the total stress can be
    #calculated directly from this increment

    if (lam > 1.0) or (np.abs(lam-1.0) < 0.005):
        etotal=etotal + estrain[i]
        sig=sig + xmod* estrain[i]
        continue
      
    "caclulate the stress and strain at yield point"

    etotal[i] = lcnd[i] + lam*estrain[i]
    sig = sig + lam*xmod*estrain[i]    
    
    #Plasticity has occured therefore the increment has to be separated in to parts
    #and calculated.
    
    #calculate the increment of strain by separating the loading increment into 
    #smaller increments
    
    de=(estrain[i]-etotal[i])/steps
    
    for k in range(0,(steps)):
        #develop increments of strain starting from previous total strain
        if (k==steps):
            el=lcnd[i+1]
        else:
            el=etotal[i]+de*(k)
            
            inc[i] = inc[i]+ 1
            
            #xback=back(a,c,xbackprev,plas,epsp)
        
            #use Netwon-Raphson to determine the increment of plastic strain
            epsp=newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback)
        
            #Use this plastic strain to calculate the stress. Note: since this is a 
            #strain-controlled analysis, the increment in strain is the increment of
            #used in the Newton-Raphson input.
            
            xbackprev=back(a,c,xbackprev,plas,epsp)
            enew[k]=el
            snew[k]=xmod*(el-epsp)
            
            
            
             
#remove end zeros from array of stress and strain
enew=np.trim_zeros(enew,'b')
snew=np.trim_zeros(snew,'b')

        
plt.plot(enew,snew)
plt.plot(fstrain,fstress,'ro')
        