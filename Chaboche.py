# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt


class oneD_plasticity:
   
    def __init__(self,backstress,iso,sigy0,xmod,readfile,steps):
        
        "find number of backstress"
        nbs=int(len(backstress)/2.)
        "add backstress parameters to individual arrays"
        aval=np.array(backstress[0:nbs])
        cval=np.array(backstress[nbs:len(backstress)])
        self.aval=aval
        self.cval=cval
        
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
        
        self.lcnd=lcnd
        
        "number of data points per branch"
        self.steps=steps


    def plotter_totalvals(self):
        totalstrain=self.strain
        totalstress=self.stress
       
        "plott total stress and strain"
        figtotal= plt.figure()
        
        plt.plot(totalstrain,totalstress)
        
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
        
        plt.show()
        

    def AF_Model(self):
        
        a=self.aval
        c=self.cval
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
        max_it=int(1e10)
        toler = 1e-10
        depsp=0
        
       
        epsp=0
        sig=0
        epspprev=0
        plasprev=0
        sigend=0
        xbackprev=np.zeros(len(a))
        plas=0
        xbackprevcy=0
        rsprev=0
        plasprev=0
        
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
        
        def back(a,c,xbackprev,plas,epsp,plasprev,epspprev,xbackprevcy,nu):
            xback=(a*(epsp-epspprev) -c*nu*(xbackprev)*(epsp-epspprev))+xbackprev
            #xback=a*(epsp) -c*(xbackprev)*(epsp)
            return xback
            
        "derivative of backstress w.r.t plastic strain "   
        def dback(a,c,xbackprev,plas,nu,xbackprevcy):
        
            dxback=a - c*(xbackprev)*nu
            return dxback
        
        "nonlinear isotropic hardening using forward Euler integration"
        def iso(Qs,bs,plas,plasprev,rsprev):
            rs=bs*(Qs-rsprev)*(plas-plasprev)+rsprev
            return rs
        
        "derivative of isotropic hardening w.r.t plastic strain"
        def diso(Qs,bs,rsprev):
            drs=bs*(Qs-rsprev)
            return drs
        
        "newton-raphson method"
        def newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,xbackprevcy,rs,rsprev):
        
        
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
                    
                    dxback=dback(a,c,xbackprev,plas,nu,xbackprevcy)
                    dxsum=np.sum(dxback)
                    drs=diso(Qs,bs,rsprev)
                    dfunc = nu*(-xmod - dxsum)-drs
                    depsp=-func/dfunc
                    
                    #epsp=(epsp)+flow*depsp
                    epsp += depsp
                    plas += nu*depsp
                    
                    "update backstress using new plastic strain increment"
                    xback=back(a,c,xbackprev,plas,epsp,plasprev,epspprev,xbackprevcy,nu)
                    
                    "update isotropic hardening using new plastic strain increment"
                    rs= iso(Qs,bs,plas,plasprev,rsprev)
                    
            return epsp,plas,xback,rs
                
        estrain=np.diff(lcnd)
        
        
        "create a loop to loop through the increments of strain"
        for i in range(0,len(estrain)):
            
            "starting strain"
            estart=lcnd[i+1]
            
        
            "calculate the current increment in stress from increment of strain"
            dsig=xmod*estrain[i]
        
            "update values at the end of previous branch to the start of the current branch of loading"
        
            xback=xbackprev
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
                epsp,plas,xback,rs=newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,xbackprevcy,rs,rsprev)
                
                
                xbackprev=xback
                epspprev=epsp
                rsprev=rs
                plasprev=plas
                
                "calculate the new stress for the current increment of strain"    
                enew.append(el)
                snew.append(xmod*(el-epsp))
           

            sigend=xmod*(el-epsp)
            
            self.strain=enew
            self.stress=snew
                    
             

"order required C1,C2,C3,Cn..,gamma1,gamma2,gamma3,gamman.."
backstress=[7500,1000,30,310,300,4]


"if you are using isotropic hardening need to answer yes"
"order required Q,b"
iso=['yes',-90.,1.8]

"material yield stress"
sigy0=200.0

"material elastic modulus"
xmod=69000.0

"if you want to read turning points from file answer 'yes' if not specify max"
"and minimum values interested in and the number of cycles"

readfile=['no',0.052,0.032,100]

"number of data points per branch"
steps=80
    
             
        
loadData=oneD_plasticity(backstress,iso,sigy0,xmod,readfile,steps)
loadData.AF_Model()
loadData.plotter_totalvals()
loadData.plotter_meanstress()
            