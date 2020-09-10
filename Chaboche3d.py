# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 09:55:55 2020

@author: mi19356
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy 
from sympy import *
from sympy.solvers import solve
import copy
class oneD_plasticity:
    def __init__(self):
        pass
    def af():    
        v=0.3
        Emod=51000.
        
        
        "need to include here how you want the uniaxial loading to be, either "
        "uniaxial stress or strain which is different in terms of components of "
        "stress as explained in https://csmbrannon.net/2012/08/02/distinction-between-uniaxial-stress-and-uniaxial-strain/"
        straincontrol=0.012
        
        
        sigy=200.
        aval=[7500,1000,30]
        cval=[310,300,4]
        
        xback=[[]*3]*6
        sig=np.zeros(6)
        lam=np.zeros(6)
        sig=np.zeros(6)
        sigend=np.zeros(6)
        xbackprev=[[0]*3]*6
        epsp=np.zeros(6)
        plas=0
        epspprev=np.zeros(6)
        plasprev=0
        depsp=0.0000001
        
        stresvm=np.zeros((3,3))
        
        sigtol=1E-10
        toler=1E-10
        steps=10
        
        max_it=10000000
        
        enew=[]
        snew=[]
        
         
        
        
        
        
        
         #read in necessary information
        straincontrol=np.array(np.genfromtxt("loading_condition.txt",unpack=True))
        
        hookes=np.array([[1-v,v,v,0.,0.,0.],[v,1-v,v,0.,0.,0.],[v,v,1-v,0.,0.,0.],
                        [0.,0.,0.,(1.-(2.*v))/2.,0.,0.],[0.,0.,0.,0.,(1.-(2.*v))/2.,0.],
                        [0.,0.,0.,0.,0.,(1.-(2.*v))/2.]])*(Emod/((1+v)*(1-2*v)))
        
        hookesst_st=np.array([[1,-v,-v,0.,0.,0.],[-v,1,-v,0.,0.,0.],[-v,-v,1,0.,0.,0.],
                        [0.,0.,0.,2.+(2.*v),0.,0.],[0.,0.,0.,0.,2.+(2.*v),0.],
                        [0.,0.,0.,0.,0.,2.+(2.*v)]])*(1./Emod)
        
        
       #lcnd=np.array([-v*straincontrol,straincontrol,-v*straincontrol,np.zeros(len(straincontrol)),np.zeros(len(straincontrol)),np.zeros(len(straincontrol))])
        lcnd=np.array([np.zeros(len(straincontrol)),straincontrol,np.zeros(len(straincontrol)),np.zeros(len(straincontrol)),np.zeros(len(straincontrol)),np.zeros(len(straincontrol))])
       
        strain=np.diff(lcnd,axis=1)
        
        
        for i in range(0,np.size(strain,axis=1)):
            
            dsig=np.matmul(hookes,strain[:,i])
            
            "starting strain"
            estart=lcnd[:,i+1]
            
            "direction of strain increment"
            nu=np.sign(strain[:,i])
            
            "update with previous values"
            sig=sigend
           # xback=xbackprev
            
            
            "total backstress"
            
            xsum=np.sum(xback,axis=1)
            
            "solve the strain at the yeild stress"
            sigyvec=np.array([0.,sigy,0.,0.,0.,0.])
            sig=nu*sigyvec+xsum-sig
            
            etotal=np.matmul(np.linalg.inv(hookes),sig)
            
            "check to see if there is plasticity"
           #zeroind=np.where(abs(dsig)>sigtol)[0]
           # lam[zeroind] = (nu[zeroind]*sigy + xsum[zeroind] - sig[zeroind])/dsig[zeroind]
            
            "caclulate the stress and strain at yield point"
                    
          #  etotal = lcnd[:,i] + lam*strain[:,i]
            #sig = sig + lam*np.matmul(hookes,strain[:,i])
                   
            "add the stress at yield"     
            enew.append(etotal)
            snew.append(sig)
            
            "increment of strain"
            de=(estart-etotal)/steps
            
            "distribute increments across number of steps"
            xbackval=[[0]*3]*6
            xbackval[0]=[[0]*3]*6
            
            for k in range(1,steps+1):
                       #develop increments of strain starting from previous total strain
                if (k==steps):
                    el=lcnd[:,i+1]
                else:
                    el=etotal+de*k
                    
                           
                #"call newton-raphson algorithm to solve for plastic increment"       
                epsp,plas,xback=newtraphson(aval,cval,epsp,plas,nu,el,xbackprev,xback,depsp,epspprev,plasprev,max_it,toler,hookes,sigy,xbackval,k)
                
               
                xbackprev=copy.deepcopy(xback)
                epspprev=copy.deepcopy(epsp)
                plasprev=copy.deepcopy(plas)
                        
                            
                enew.append(el)
                snew.append(np.matmul(hookes,el-epsp))
                

                
def back(aval,cval,xbackprev,plas,epsp,plasprev,epspprev,xback):
    for i in range(len(epsp)):
        xback[i]=((2./3.)*np.dot(aval,(epsp[i]-epspprev[i])) - np.dot(cval,(xbackprev[i]))*(plas-plasprev))+xbackprev[i]
    
    return xback
            
"derivative of backstress w.r.t plastic strain "   
def dback(aval,cval,xbackprev,normal,depsp):
    "set up symbols for derivative"
    depsp1,depsp2,depsp3,depsp3,depsp4,depsp5,depsp6,normal1,normal2,normal3,normal4,normal5,normal6=symbols('depsp1 depsp2 depsp3 depsp3 depsp4 depsp5 depsp6 normal1 normal2 normal3 normal4 normal5 normal6',real=True)
    "inialise arrays"
    depsp=np.array([depsp1,depsp2,depsp3,depsp4,depsp5,depsp6])
    normal=np.array([normal1,normal2,normal3,normal4,normal5,normal6])
    epsp=depsp*normal
    "calculate increment in effective plastic strain"
    plas=(np.dot((3./2.)*epsp,epsp))**0.5
    
    "calculate backstress"
    func=[]
    for i in range(len(epsp)):
        func.append(aval*epsp[i]-cval[i]*np.dot(xbackprev[i],plas))
        
    #dxback=a - c*(xbackprev)*nu
    return 

def vonmises1(sig,xsum):
    #"find deviatoric stress using symmetry to convert"
        #"stress vector to tensor"
  
    #"hydrostatic stress"
    stresvm[0,0]=np.trace(sig)/3.0
    stresvm[1,1]=np.trace(sig)/3.0
    stresvm[2,2]=np.trace(sig)/3.0
    # "deviatoric stress"
    stressprime=sig-stresvm
        
    #"convert stress prime to vector using symmetry"
    stressprimevec=np.array([stressprime[0,0],stressprime[1,1],stressprime[2,2],
                                 stressprime[0,1],stressprime[1,0],stressprime[2,0]])
        
           #"yeild function"
    func=np.dot((3./2.)*stressprimevec-xsum,stressprimevec-xsum)**0.5
    return func

def vonmises2(sig_x):
    
    func=((((sig_x[0,0]-sig_x[1,1])**2.)+((sig_x[1,1]-sig_x[2,2])**2.)+
           ((sig_x[2,2]-sig_x[0,0])**2.)+6.*((sig_x[0,1]**2.0)+(sig_x[1,2]**2.0)
                                             +(sig_x[2,0]**2.0)))/2.0)**0.5
    return func

def normal_vonmises(sig,sigy,xsum):
    sig_x1,sig_x2,sig_x3,sig_x12,sig_x23,sig_x31=symbols('sig_x1 sig_x2 sig_x3 sig_x12 sig_x23 sig_x31',real=True)

    sig_total=np.array([sig_x1,sig_x2,sig_x3,sig_x12,sig_x23,sig_x31])
    
    func=(((((sig_x1-sig_x2)**2.)+((sig_x2-sig_x3)**2.)+
           ((sig_x3-sig_x1)**2.)+6.*((sig_x12**2.0)+(sig_x23**2.0)
                                             +(sig_x31**2.0)))/2.0)**0.5)-sigy
    df_dsig=diff(func,sig_total)
    
    "substitute values"
    
    df_dsigreal=df_dsig.subs({sig_x1:sig[0,0]-xsum[0],sig_x2:sig[1,1]-xsum[1],
                              sig_x3:sig[2,2]-xsum[2],sig_x12:sig[0,1]-xsum[3],
                              sig_x23:sig[1,2]-xsum[4],sig_x31:sig[2,0]-xsum[5]})
    return df_dsigreal

def dfunc(hookes,el,sigy,epspr,aval,cval,xbackprevt,depspr,normalr,xback2,epspprev,plasprev):
    depsp,normal1,normal2,normal3,normal4,normal5,normal6=symbols('depsp normal1 normal2 normal3 normal4 normal5 normal6',real=True)
    #epsp1,epsp2,epsp3,epsp4,epsp5,epsp6=symbols('epsp1 epsp2 epsp3 epsp4 epsp5 epsp6')
    epsp2=depsp*np.array([normal1,normal2,normal3,normal4,normal5,normal6])+epspr
    #epsp=np.array([epsp1,epsp2,epsp3,epsp4,epsp5,epsp6])
    

    
    "calculate increment in effective plastic strain"
    plas=(np.dot((2./3.)*epsp2,epsp2))**0.5
    
    "calculate backstress"
   
    for i in range(len(epsp2)):
       #xback2[i]=np.sum(np.dot((2./3.),np.dot(aval,epsp2[i]))-cval*np.dot(xbackprevt[i],plas))
       xback2[i]=np.sum(((2./3.)*np.dot(aval,(epsp2[i]-epspprev[i])) - np.dot(cval,(xbackprevt[i]))*(plas-plasprev))+xbackprevt[i])
    
    
    
    "check to find zeros in real array"
    #zeroind=np.where(epspr==0)[0]
    "don't use these values in derivative"
   # epsp2[zeroind]=0.    
    "check to see if the whole array is zero, if so, only use 3 components"
    #if (len(zeroind)==6):
    #    #epsp=np.array([epsp1,epsp2,epsp3,0,0,0])
        
     #   epsp2=np.array([epsp2[0],epsp2[1],epsp2[2],0,0,0])
        
            
    #depsp=np.array([depsp1,depsp2,depsp3,depsp4,depsp5,depsp6])
    #normal=np.array([normal1,normal2,normal3,normal4,normal5,normal6])
    "since this is uniaxial, the vectors are relatated"
    #eppvec=np.array([(el[1]-epsp2[1])*-0.3,el[1]-epsp2[1],(el[1]-epsp2[1])*-0.3,0,0,0])
    sig=np.matmul(hookes,el-epsp2)
    
    "must change sigy to a function for isotropic hardening"
    "remove the backstress from the stress"
    sig_x=sig-xback2
    
    func=((((sig_x[0]-sig_x[1])**2.)+((sig_x[1]-sig_x[2])**2.)+
           ((sig_x[2]-sig_x[0])**2.)+6.*((sig_x[3]**2.0)+(sig_x[4]**2.0)
                                             +(sig_x[5]**2.0)))/2.0)**0.5-sigy
    "find derivative with respect to increment of strain"
   # if len(zeroind)==6:
        #epsp=np.array([epsp1,epsp2,epsp3])
    #    epsp2=epsp2[:-3]
     #   dfunc_epsp=diff(func,depsp)
        
    #else:  
     #   dfunc_epsp=diff(func,depsp)
        
    dfunc_epsp=diff(func,depsp)
    "substitute epsp"
    #dfunc_epspreal=dfunc_epsp.subs({depsp1:depspr[0],depsp2:depspr[1],depsp3:depspr[2],depsp4:depspr[3],
    #                                depsp5:depspr[4],depsp6:depspr[5],normal1:normalr[0],normal2:normalr[1],
    #                                normal3:normalr[2],normal4:normalr[3],normal5:normalr[4],normal6:normalr[5]})
   # if len(zeroind)==6:
        
   #     dfunc_epspreal=dfunc_epsp.subs({epsp1:epspr[0],epsp2:epspr[1],epsp3:epspr[2]})
   # else:
      #  dfunc_epspreal=dfunc_epsp.subs({epsp1:epspr[0],epsp2:epspr[1],epsp3:epspr[2],epsp4:epspr[3],epsp5:epspr[4],epsp6:epspr[5]})

    #dfunc_epspreal=dfunc_epsp.subs({epsp1:epspr[0],epsp2:epspr[1],epsp3:epspr[2],epsp4:epspr[3],epsp5:epspr[4],epsp6:epspr[5]})
    #xbackprevt=xbackprevt[0].subs({depsp:depspr,normal1:normalr[0],normal2:normalr[1],normal3:normalr[2],normal4:normalr[3],normal5:normalr[4],normal6:normalr[5]})
    dfunc_epspreal=dfunc_epsp.subs({depsp:depspr,normal1:normalr[0],normal2:normalr[1],normal3:normalr[2],normal4:normalr[3],normal5:normalr[4],normal6:normalr[5]})
    
   
    return dfunc_epspreal

    
def newtraphson(aval,cval,epspn,plas,nu,el,xbackprevn,xback,depsp,epspprev,plasprev,max_it,toler,hookes,sigy,xbackval,kval):
    
    
    for n in range(0,max_it):
    
        "test current increment in plastic strain"
        
        sign=np.matmul(hookes,el-epspn)
        xsum=np.sum(xback,axis=1)
        
        "remove backstress"
        sign=sign-xsum
        sigten=np.array([[sign[0],sign[3],sign[5]],[sign[3],sign[1],sign[4]],
                         [sign[5],sign[4],sign[2]]])
        
        funcn=vonmises2(sigten)-sigy
        
        normal=normal_vonmises(sigten,sigy,xsum)
         
        if(abs(funcn)<toler):
            return epspn,plas,xback
        
        df=dfunc(hookes,el,sigy,epspn,aval,cval,xbackprevn,depsp,normal,xback,epspprev,plasprev)

        depsp=-funcn/df
        
       
          
                   
        epspn =epspn+depsp*normal
        plas =(np.dot((2./3.)*epspn,epspn))**0.5
        
        xback=back(aval,cval,xbackprevn,plas,epspn,plasprev,epspprev,xback)
       
        
                    
            
                
            
    return epspn,plas,xback
                
oneD_plasticity.af()