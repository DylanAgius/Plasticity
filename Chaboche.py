# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy import optimize

class oneD_plasticity:
   
    def __init__(self):
        pass

    def AF_Model():
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
        max_it=100000000000000
        toler = 1e-10
        epsp_new=0
        sig=0
        
        
        depsp=0
        a=np.array([7500,1000,30])
        c=np.array([310,300,4])
        Qs=100.0
        bs=6.8
        sigy0=200.0
        
        #step size
        steps=10
        
        #initialise arrays
        estrain=np.zeros(len(lcnd))
        nu=np.zeros(len(lcnd))
        dsig=np.zeros(len(lcnd))
        inc=[0 for i in range(10000)]
        #dxback=np.zeros(3)
        #xback=np.zeros(3)
        xbackold=np.zeros(3)
        lam=0
        #enew=np.zeros(10000)
        #snew=np.zeros(10000)
        
        enew=[]
        snew=[]
        
        etotal=np.zeros(10000)
        #sig=np.zeros(10000)
        el=np.zeros(10000)
        enewtotal=[]
        snewtotal=[]
        epspprev=0
        plasprev=0
        sigend=0
        
        
        estrain[0] = lcnd[1]
        
        #xback=np.empty([3,10000])
        #dxback=np.empty([3,10000])
        
        xbackprev=np.zeros(3)
        plas=0
        xbackprevcy=0
        rsprev=0
        plasprev=0
        
        #inc[0]=0
        
        "nonlinear backstress calculation"
        
        def back(a,c,xbackprev,plas,epsp,plasprev,epspprev,xbackprevcy,nu):
            xback=(a*(epsp-epspprev) -c*nu*(xbackprev)*(epsp-epspprev))+xbackprev
            #xback=a*(epsp) -c*(xbackprev)*(epsp)
            return xback
            
        "derivative of backstress w.r.t plastic strain "   
        def dback(a,c,xbackprev,plas,nu,xbackprevcy):
        
            dxback=a - c*(xbackprev)*nu
            return dxback
        
        "nonlinear isotropic hardening"
        def iso(Qs,bs,plas,plasprev,rsprev):
            rs=bs*(Qs-rsprev)*(plas-plasprev)+rsprev
            return rs
        
        def diso(Qs,bs,rsprev):
            drs=bs*(Qs-rsprev)
            return drs
        
        #newton-raphson method
        def newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,xbackprevcy,rs,rsprev):
        
        
            for n in range(0,max_it):
                
                sig=xmod*(el-epsp)
                "von mises stress invariant"
              
              
                strvm=abs(sig-np.sum(xback))
               #need to make this state variable
                sigy=200+rs
            #check to see if the point remains in the yield surface
                func=strvm-sigy
           
                #func=strvm-sigy
                if(abs(func)<toler):
                    return epsp,plas,xback,rs
                else:
                    
                    dxback=dback(a,c,xbackprev,plas,nu,xbackprevcy)
                    dxsum=np.sum(dxback)
                    drs=diso(Qs,bs,rsprev)
                    dfunc = nu*(-xmod - dxsum)-drs
                    depsp=-func/dfunc
                    
                    #epsp=(epsp)+flow*depsp
                    epsp += depsp
                    plas += nu*depsp
                    
                    #epsp=epsp+depsp
                    xback=back(a,c,xbackprev,plas,epsp,plasprev,epspprev,xbackprevcy,nu)
                    rs= iso(Qs,bs,plas,plasprev,rsprev)
                    
            return epsp,plas,xback,rs
                
        estrain=np.diff(lcnd)
        
        
        "create a loop to loop through the increments of strain"
        for i in range(0,len(estrain)):
            
            "starting strain"
            estart=lcnd[i+1]
            
        
            "calculate the current increment in stress from increment of strain"
            dsig=xmod*estrain[i]
        
            "calculate the backstress"
        
            xback=xbackprev
            sig=sigend
            rs=rsprev
        
            "loading direction provided by sign of the strain increment"
            nu=np.sign(estrain[i])
        
        #NEED TO CALCULATE THE BACKSTRESS HERE
        
            "now we know need to check with the current increment in stress is greater"
            "than the yield taking into account the shift as advised by the backstress"
            "eg plasticity will occur when sigma_old+dsig > sig_yield + backstress"
            
           # xback=back(a,c,xbackprev,plas,epsp,plasprev,epspprev,xbackprevcy,nu)
        
            "total backstress"
            xsum=np.sum(xback)
            
            "isotropic hardening"
           #rs= iso(Qs,bs,plas,plasprev,rsprev)
            sigy=sigy0+rs
            
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
            
            enew.append(etotal[i])
            snew.append(sig)
        
            de=(estart-etotal[i])/steps
            
            
            for k in range(1,steps+1):
                #develop increments of strain starting from previous total strain
                if (k==steps):
                    el=lcnd[i+1]
                else:
                    el=etotal[i]+de*k
                    
                inc[i] = inc[i]+ 1
                    
                #xback=back(a,c,xbackprev,plas,epsp)
                
        
                
                epsp,plas,xback,rs=newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,xbackprevcy,rs,rsprev)
                
               
                    #Use this plastic strain to calculate the stress. Note: since this is a 
                    #strain-controlled analysis, the increment in strain is the increment of
                    #used in the Newton-Raphson input.
                
                
                xbackprev=xback
                epspprev=epsp
                rsprev=rs
                plasprev=plas
                
                    
                enew.append(el)
                snew.append(xmod*(el-epsp))
                
               
                    
             

            sigend=xmod*(el-epsp)      
                    
                     
                
        plt.plot(enew,snew)
        plt.plot(fstrain,fstress,'ro')
        plt.show()
        
oneD_plasticity.AF_Model()
            