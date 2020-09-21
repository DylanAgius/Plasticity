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
import sympy as sym
from tqdm import tqdm


class oneD_plasticity:
   
    def __init__(self,vkin,backstress,iso,sigy0,xmod,readfile,steps):
        "save the model type"
        self.ktype=vkin[0]
    
        "type of kinematic hardening"
        if (vkin[0]=='MAF'):
        
            "find number of backstress"
            nbs=int(len(backstress)/2.)
            "add backstress parameters to individual arrays"
            aval=np.array(backstress[0:nbs])
            cval=np.array(backstress[nbs:len(backstress)])
            self.aval=aval
            self.cval=cval
        elif (vkin[0]=='MAFM'):
            "save number of backstresses"
            self.MAFMnum=vkin[1]
            "check to see the number of multiplicative backstresses"
            nbs=int((len(backstress)-int(vkin[1]*4))/2.)
            "add backstress parameters to individual arrays"
            aval=np.array(backstress[0:nbs])
            cval=np.array(backstress[nbs+int(vkin[1]*2):(nbs*2)+int(vkin[1]*2)])
            "now create arrays for the multiplicative backstress"
            am=[]
            cm=[]
            totallen=int(nbs+int(vkin[1]*2))
            for i in range(0,int((vkin[1]*2))):
                am.append(backstress[nbs+i]/backstress[totallen+nbs+i])
                cm.append(backstress[totallen+nbs+i])
            self.aval=aval
            self.cval=cval
            self.am=am
            self.cm=cm
                
        
        "determine if isotropic hardening is to be used"
        if (iso[0]=='yes'):
            Qs=iso[1]
            bs=iso[2]
        else:
            Qs=0.0
            bs=0.0
        self.Qs=Qs
        self.bs=bs
        
        "save yield strength"
        self.sigy0=sigy0
        
        "save elastic modulus"
        self.xmod=xmod
        
        "determine if reading turning points from file"
        if (readfile[0]=='no'):
            lcnd=np.tile([readfile[2],readfile[1]],readfile[3])
            lcnd=np.insert(lcnd,0,0)
            lcnd=np.insert(lcnd,1,readfile[1])
        else:
            readtp=open(readfile[1],"r")
            lcnd=readtp.read().splitlines()
            readtp.close()
            lcnd=np.asarray(lcnd).astype(np.float64)
               
    
        
        self.lcnd=lcnd
        
        "number of data points per branch"
        self.steps=steps


    def plotter_totalvals(self):
        totalstrain=self.strain
        totalstress=self.stress
       
        "plott total stress and strain"
        figtotal= plt.figure()
        
        plt.plot(totalstrain,totalstress)
        plt.xlabel('Strain')
        plt.ylabel('Stress (MPa)')
        
        plt.show()
        
    def plotter_meanstress(self):
        steps=self.steps
        
        totalstrain=self.strain
        totalstress=self.stress
        
        "calculate mean stress"
        maxstress=totalstress[steps+1::2*(steps+1)]
        minstress=totalstress[2*(steps+1)::2*(steps+1)]
        
        meanstress=[]
        for i in range(len(maxstress)-1):
            meanstress.append((maxstress[i]+minstress[i])/2.)
        
        cyclecount=range(1,len(maxstress))
        
        "plott total stress and strain"
        figtotal= plt.figure()
        
        plt.plot(cyclecount,meanstress)
        
        plt.xlabel('Cycle')
        plt.ylabel('Mean Stress (MPa)')
        
        plt.show()
        

    def Plast_Model(self):
        
        "initlise backstress parameters"
        ktype=self.ktype
        if ktype=='MAF':
            a=self.aval
            c=self.cval
        elif ktype=='MAFM':
            MAFMnum=self.MAFMnum
            a=self.aval
            c=self.cval
            am=self.am
            cm=self.cm
            astar=am[MAFMnum:MAFMnum*2]
            cstar=cm[MAFMnum:MAFMnum*2]
            am=self.am[0:MAFMnum]
            cm=self.cm[0:MAFMnum]
           
        
        Qs=self.Qs
        bs=self.bs
        sigy0=self.sigy0
        xmod=self.xmod
        lcnd=self.lcnd
        steps=self.steps
      
        
        
        
        #read in necessary information
        #lcnd=np.genfromtxt("loading_condition.txt",unpack=True)
        
        #data calculated using fortran
        #fstrain,fstress=np.genfromtxt("af_stress_strain.txt", unpack=True)
        
        "initialise values"
        "netwon raphson parameters"
        max_it=int(1e5)
        toler = 1e-9
        depsp=0
        epsp=0
        sig=0
        epspprev=0
        plasprev=0
        sigend=0
        plas=0
        rsprev=0
        plasprev=0
        
        "intialise backstresses"
        if ktype=='MAF':
            xbackprev=np.zeros(len(a))
        elif ktype=='MAFM':
            xbackprev=np.zeros(len(a)+len(am))
            xbackmafmprev=np.zeros(len(am))
            xbackstarprev=np.zeros(len(astar))
        "initialise arrays"
        estrain=np.zeros(len(lcnd))
        nu=np.zeros(len(lcnd))
        dsig=np.zeros(len(lcnd))
        inc=[0 for i in range(10000)]
        enew=[0]
        snew=[0]
        etotal=np.zeros(10000)
        el=np.zeros(10000)
        enewtotal=[]
        snewtotal=[]
       
        
        "nonlinear backstress calculation using forward Euler integration"
        
        def back(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu):
            xback=(a*(epsp-epspprev) -c*nu*(xbackprev)*(epsp-epspprev))+xbackprev
            #xback=a*(epsp) -c*(xbackprev)*(epsp)
            return xback
        
        "nonlinear MAFM backstress calculation using forward Euler integration"
        
        def backMAFM(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstarprev,xbackstar):
            "classic basktresses"
            xback[0:len(a)]=(a*(epsp-epspprev) -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)])*(epsp-epspprev))+xbackprev[0:(len(xbackprev)-MAFMnum)]
            xbackstar=cstar*(astar - nu*xbackstarprev)*(epsp-epspprev) + xbackstarprev
            xback[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]=(((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]))*(epsp-epspprev)) + xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]
        
            return xback,xbackstar
        
        def backMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstarprev,xbackstar):
            "classic basktresses"
            #am=np.array(am)
            #cm=np.array(cm)
            #astar=np.array(astar)
            #cstar=np.array(cstar)
           
            k1=a*(epsp-epspprev) -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)])*(epsp-epspprev)
            k2=a*(epsp-epspprev) -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)]+(k1/2.))*((epsp-epspprev)/2)
            k3=a*(epsp-epspprev) -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)]+(k2/2.))*((epsp-epspprev)/2)
            k4=a*(epsp-epspprev) -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)]+k3)*(epsp-epspprev)
            
            xback[0:len(a)]=xbackprev[0:(len(xbackprev)-MAFMnum)] + (1/6)*(k1+(2.*k2)+(2.*k3)+k4)
            
            "multiplicative backstress"
           
            k11=cstar*(astar - nu*xbackstarprev)*(epsp-epspprev)
            k22=cstar*(astar - nu*(xbackstarprev+(k11/2)))*((epsp-epspprev)/2)
            k33=cstar*(astar - nu*(xbackstarprev+(k22/2)))*((epsp-epspprev)/2)
            k44=cstar*(astar - nu*(xbackstarprev+k33))*(epsp-epspprev)
           
            xbackstar=xbackstarprev + (1/6)*(k11+(2.*k22)+(2.*k33)+k44)
            
            "fourth backstress"
            k13=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]))*(epsp-epspprev)
            k23=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*((xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])+(k13/2))))*((epsp-epspprev)/2)
            k33=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*((xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])+(k23/2))))*((epsp-epspprev)/2)
            k43=(((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*((xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])+k33)))*(epsp-epspprev))

            xback[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]=xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum] + (1/6)*(k13+(2.*k23)+(2.*k33)+k43)
            return xback,xbackstar
        
        def dbackMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstarprev,xbackstar):
           
            epspt=symbols('epspt',real=True)
            "classic basktresses"
            am=np.array(am)
            cm=np.array(cm)
            astar=np.array(astar)
            cstar=np.array(cstar)
           
            k1=a*epspt -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)])*epspt
            k2=a*epspt -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)]+(k1/2.))*(epspt/2)
            k3=a*epspt -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)]+(k2/2.))*(epspt/2)
            k4=a*epspt -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)]+k3)*epspt
            
            xback=xbackprev[0:(len(xbackprev)-MAFMnum)] + (1/6)*(k1+(2.*k2)+(2.*k3)+k4)
            
         
            
            "fourth backstress"
            k11=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]))*epspt
            k22=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*((xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])+(k11/2))))*(epspt/2)
            k33=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*((xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])+(k22/2))))*(epspt/2)
            k44=((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*((xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])+k33)))*epspt

            xback=np.append(xback,(xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum] + (1/6)*(k11+(2.*k22)+(2.*k33)+k44)))
            
            dfunc_epsp=diff(xback,epspt)
            
            dxback=dfunc_epsp.subs({epspt:(epsp-epspprev)})
            return dxback
            
        "derivative of backstress w.r.t plastic strain "   
        def dback(a,c,xbackprev,plas,nu):
        
            dxback=a - c*(xbackprev)*nu
            return dxback
        
        "derivative of backstress w.r.t plastic strain for MAFM "   
        def dbackMAFM(a,c,am,cm,astar,cstar,xbackprev,xbackstarprev,plas,nu,back,xbackstar,epsp,epspprev):
        
            dxback=a - c*(xbackprev[0:(len(xbackprev)-MAFMnum)])*nu
            #dbackstar=cstar*(astar - nu*xbackstarprev)
         
            #xbackstar=cstar*(astar - nu*xbackstarprev)*(epsp-epspprev) + xbackstarprev
            #dxbackmafm=((cm+cstar*(astar - nu*xbackstar))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]))
            #dxbackmafm=np.dot(((cm+cstar*(astar - nu*xbackstar))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])),-cstar*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])*dbackstar)

            epspt=symbols('epspt',real=True)
            
            #am=np.array(am)
            #cm=np.array(cm)
            #astar=np.array(astar)
            #cstar=np.array(cstar)
           
            backmafmsolv=(cm+cstar*(astar - nu*xbackstarprev))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum])*epspt + xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]
            dfunc_epsp=diff(backmafmsolv,epspt)
            
            dfunc_epspreal=dfunc_epsp.subs({epspt:(epsp-epspprev)})
            
            dxback=np.append(dxback,dfunc_epspreal)
            return dxback
        "nonlinear isotropic hardening using forward Euler integration"
        def iso(Qs,bs,plas,plasprev,rsprev):
            rs=bs*(Qs-rsprev)*(plas-plasprev)+rsprev
            return rs
        
        "derivative of isotropic hardening w.r.t plastic strain"
        def diso(Qs,bs,rsprev,nu):
            drs=bs*(Qs-rsprev)*nu
            return drs
        
        "nonlinear isotropic hardening using forward Euler integration"
        def isoback(Qs,bs,plas,plasprev,rsprev):
            rs=symbols('rs',real=True)
            
            funcrs=bs*(Qs-rs)*(plas-plasprev)+rsprev -rs
            
            funcsol=solve(funcrs,rs)
            
            return funcsol
        
        "derivative of isotropic hardening w.r.t plastic strain"
        def disoback(Qs,bs,rsprev,plas,plasprev):
            rs,plast=symbols('rs plast',real=True)
            
            funcrs=bs*(Qs-rs)*plast+rsprev -rs
            drsfunc=diff(funcrs,plast)
            
            funcsol=solve(drsfunc,rs)
            
           
            #drs=funcsol.subs({plast:(plas-plasprev)})
            
           
            return funcsol[0]
        
        "runge-kutta"
        def AFRK4(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu):
            
            k1=a*(epsp-epspprev) -c*nu*(xbackprev)*(epsp-epspprev)
            k2=a*(epsp-epspprev) -c*nu*(xbackprev+(k1/2.))*(epsp-epspprev)
            k3=a*(epsp-epspprev) -c*nu*(xbackprev+(k2/2.))*(epsp-epspprev)
            k4=a*(epsp-epspprev) -c*nu*(xbackprev+k3)*(epsp-epspprev)
            
            xback=xbackprev+ (1/6)*(k1+(2.*k2)+(2.*k3)+k4)
            
            return xback
        
        "runge-kutta"
        def dAFRK4(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu):
            epspt=symbols('epspt',real=True)
            
            k1=a*epspt -c*nu*(xbackprev)*epspt
            k2=a*epspt -c*nu*(xbackprev+(k1/2.))*epspt
            k3=a*epspt -c*nu*(xbackprev+(k2/2.))*epspt
            k4=a*epspt-c*nu*(xbackprev+k3)*epspt
          
            
            dxbackfunc= xbackprev+ (1/6)*(k1+ (2.*k2)+(2.*k3)+k4)
            
            dfunc_epsp=diff(dxbackfunc,epspt)
            
            dxback=dfunc_epsp.subs({epspt:(epsp-epspprev)})
            
            
            
            return dxback
                        
        
        "newton-raphson method"
        def newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev):
        
        
            for n in range(0,max_it):
                
                sig=xmod*(el-epsp)
                "von mises stress invariant"
              
              
                strvm=abs(sig-np.sum(xback))
                "calculate the yield stress"
                sigy=sigy0+rs
                "check to see if the point remains in the yield surface"
                func=strvm-sigy
           
                if(abs(func)<toler):
                    return epsp,plas,xback,rs
                else:
                    
                    dxback=dback(a,c,xbackprev,plas,nu)
                    #dxback=dAFRK4(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu)
                    dxsum=np.sum(dxback)
                    drs=diso(Qs,bs,rsprev,nu)
                    dfunc = nu*(-xmod - dxsum)-drs
                    depsp=-func/dfunc
                    
                    #epsp=(epsp)+flow*depsp
                    epsp += depsp
                    plas += nu*depsp
                    
                    "update backstress using new plastic strain increment"
                    xback=back(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu)
                    #xback=AFRK4(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu)
                    "update isotropic hardening using new plastic strain increment"
                    rs= iso(Qs,bs,plas,plasprev,rsprev)
                    
            return epsp,plas,xback,rs
        
        "newton-raphson method MAFM"
        def newtraphsonMAFM(a,c,am,cm,astar,cstar,MAFMnum,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,xbackstarprev,xbackstar):
       
        
            for n in range(0,max_it):
                
                sig=xmod*(el-epsp)
                "von mises stress invariant"
              
              
                strvm=abs(sig-np.sum(xback))
                "calculate the yield stress"
                sigy=sigy0+rs
                "check to see if the point remains in the yield surface"
                func=strvm-sigy
           
                if(abs(func)<toler):
                    return epsp,plas,xback,rs,xbackstar
                else:
                    
                    dxback=dbackMAFM(a,c,am,cm,astar,cstar,xbackprev,xbackstarprev,plas,nu,back,xbackstar,epsp,epspprev)
                    #dxback=dbackMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstar,xbackstarprev)
                    
                    dxsum=np.sum(dxback)
                    drs=diso(Qs,bs,rsprev,nu)
                    
                    #drs=disoback(Qs,bs,rsprev,plas,plasprev)
                    dfunc = nu*(-xmod - dxsum)-drs
                    depsp=-func/dfunc
                    
                   
                    epsp += depsp
                    plas += nu*depsp
                    
                    "update backstress using new plastic strain increment"
                    xback,xbackstar=backMAFM(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstar,xbackstarprev)
                    #xback,xbackstar=backMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstar,xbackstarprev)
                    
                    
                    "update isotropic hardening using new plastic strain increment"
                    rs= iso(Qs,bs,plas,plasprev,rsprev)
                    #rsback=iso(Qs,bs,plas,plasprev,rsprev)
                    
            return epsp,plas,xback,rs,xbackstar
                
        estrain=np.diff(lcnd)
        
        
        "create a loop to loop through the increments of strain"
        "add pogress bar"
        pbar = tqdm(total=len(estrain))
        for i in range(0,len(estrain)):
            pbar.update((len(estrain))/(len(estrain)))
            #pbar.update(int(((i+1)/(len(estrain)-1))*100))
            "starting strain"
            estart=lcnd[i+1]
            
        
            "calculate the current increment in stress from increment of strain"
            dsig=xmod*estrain[i]
        
            "update values at the end of previous branch to the start of the current branch of loading"
            if ktype=='MAF':
                xback=xbackprev
            elif ktype=='MAFM':
                xback=xbackprev
                xbackstar=xbackstarprev
            sig=sigend
            rs=rsprev
        
            "loading direction provided by sign of the strain increment"
            nu=np.sign(estrain[i])
        
      
        
            "now we know need to check with the current increment in stress is greater"
            "than the yield taking into account the shift as advised by the backstress"
            "eg plasticity will occur when sigma_old+dsig > sig_yield + backstress"

        
            "total backstress"
            xsum=np.sum(xback)
            
            "update yield stress"
          
            sigy=sigy0+rs
            
            "check to see if there is plasticity"
            lam = (nu*sigy + xsum - sig)/dsig
        
            "if lam> 1 then the increment is elastic, therefore the total stress can be"
            "calculated directly from this increment"
        
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
                "develop increments of strain starting from previous total strain"
                if (k==steps):
                    el=lcnd[i+1]
                else:
                    el=etotal[i]+de*k
                    
                inc[i] = inc[i]+ 1
                    
                "call on the newton-raphson method to solve for the increment of plastic strain"
                if ktype=='MAF':
                    epsp,plas,xback,rs=newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev)
                
                
                    xbackprev=xback
                
                elif ktype=='MAFM':
                    epsp,plas,xback,rs,xbackstar=newtraphsonMAFM(a,c,am,cm,astar,cstar,MAFMnum,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,xbackstarprev,xbackstar)
                                             
       
                
                    xbackprev=xback
                    xbackstarprev=xbackstar
                epspprev=epsp
                rsprev=rs
                plasprev=plas
                
                "calculate the new stress for the current increment of strain"    
                enew.append(el)
                snew.append(xmod*(el-epsp))
           

            sigend=xmod*(el-epsp)
            
            self.strain=enew
            self.stress=snew
        pbar.close()
"type of kinematic hardening model to use"
"if traditiona chaboche use: MAF; if multiplicative model use: MAFM"
"if using MAFM, please specify the number of multiplicative backstresses being used."
"This example uses 1"
kinv=['MAFM',1]             

"order required C1,C2,C3,Cn..,gamma1,gamma2,gamma3,gamman.. for the MAF model"
"if using the MAFM model, specify the C and gamma values for the multiplicative backstress "
"at the end:C1,C2,C3,C4,Cstar..,gamma1,gamma2,gamma3,gamma4,gammastar..  "
#backstress=[440,1200,5100,100,800,0.1,10000,1500,18000,500]
#backstress=[40000,6900,7000,1090,1100,800,440,450,60,1]
backstress=[4400,120000,5100,117000,800,4,1000,15,1800,5000]


"if you are using isotropic hardening need to answer yes"
"order required Q,b"
iso=['no',-50, 1]

"material yield stress"
sigy0=200.0

"material elastic modulus"
xmod=69000

"if you want to read turning points from file answer 'yes' if not specify max"
"and minimum values interested in and the number of cycles"

"this reads turning points from a text file (example file given)"
readfile=['yes','loading_condition.txt']

"this reads turning points and number of cycles"
#readfile=['no',0.02,0.002,60]

"number of data points per branch"
steps=100
    
             
        
loadData=oneD_plasticity(kinv,backstress,iso,sigy0,xmod,readfile,steps)
loadData.Plast_Model()
loadData.plotter_totalvals()
#loadData.plotter_meanstress()
            