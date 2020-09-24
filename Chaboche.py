# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy 
from sympy import diff, symbols
from sympy.solvers import solve
import copy
import sympy as sym
from tqdm import tqdm
import xlwt 


class straincontrol:
   
    def __init__(self,vkin,backstress,isoparam,sigy0,xmod,readfile,steps):
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
        if (isoparam[0]=='yes'):
            Qs=isoparam[1]
            bs=isoparam[2]
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
        
        "control method"
        self.control='strain'
     
        


    def plotter_totalvals(self):
        totalstrain=self.strain
        totalstress=self.stress
       
        "plott total stress and strain"
 
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
        
        "plot total stress and strain"

        plt.plot(cyclecount,meanstress)
        
        plt.xlabel('Cycle')
        plt.ylabel('Mean Stress (MPa)')
        
        plt.show()
        
        self.meanstress=meanstress
        self.cyclecount=cyclecount
        
    def dataextract(self):
        "check to see which data has been asked for"
        #peake=self.peake
        #cyclecount=self.cyclecount
        #meanstress=self.meanstress
        totalstrain=np.asarray(self.strain).astype(np.float)
        totalstress=np.asarray(self.stress).astype(np.float)
        
        workbook = xlwt.Workbook()  
  
        sheet1 = workbook.add_sheet("Total_stress_strain") 
        style = xlwt.easyxf('font: bold 1') 
  
       
        sheet1.write(0, 0, 'Strain', style) 
        sheet1.write(0,1, 'Stress (MPa)', style)
        for i in range(0,len(totalstress)):
            
            sheet1.write(1+i,0, totalstrain[i])
            sheet1.write(1+i,1,totalstress[i])
            workbook.save("simulated_data.xls") 
            
        "check to see if mean stress present"
        if 'meanstress' in dir(self):
            sheet2=workbook.add_sheet("Mean_Stress") 
            cyclecount=np.asarray(self.cyclecount).astype(np.float)
            meanstress=np.asarray(self.meanstress).astype(np.float)
            sheet2.write(0, 0, 'Cycle', style) 
            sheet2.write(0,1, 'Mean Stress (MPa)',style)
            
            for i in range(0,len(meanstress)):
                sheet2.write(1+i,0, cyclecount[i])
                sheet2.write(1+i,1,meanstress[i])
                workbook.save("simulated_data.xls")
                
        "check to see if peak strain present"  
        if 'peake' in dir(self):
            sheet3=workbook.add_sheet("Peak_Strain") 
            cyclecount=np.asarray(self.cyclecount).astype(np.float)
            peake=np.asarray(self.peake).astype(np.float)
            sheet3.write(0, 0, 'Cycle', style) 
            sheet3.write(0,1, 'Peak Strain', style)
            
            for i in range(0,len(peake)):
                sheet3.write(1+i,0, cyclecount[i])
                sheet3.write(1+i,1,peake[i])
                workbook.save("simulated_data.xls")
        

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
        control=self.control
      
  
        
        "initialise values"
        "netwon raphson parameters"
        max_it=int(1e5)
        toler = 1e-6
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
                xback=copy.deepcopy(xbackprev)
            elif ktype=='MAFM':
                xback=copy.deepcopy(xbackprev)
                xbackstar=copy.deepcopy(xbackstarprev)
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
                    epsp,plas,xback,rs=newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,max_it,toler,Qs,bs,control,xmod,sigy0)
                
                
                    xbackprev=copy.deepcopy(xback)
                
                elif ktype=='MAFM':
                    epsp,plas,xback,rs,xbackstar=newtraphsonMAFM(a,c,am,cm,astar,cstar,MAFMnum,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,xbackstarprev,xbackstar,max_it,toler,Qs,bs,control,xmod,sigy0)
                                             
       
                
                    xbackprev=copy.deepcopy(xback)
                    xbackstarprev=copy.deepcopy(xbackstar)
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
        
class stresscontrol:
   
    def __init__(self,vkin,backstress,isoparam,sigy0,xmod,readfile,steps):
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
        if (isoparam[0]=='yes'):
            Qs=isoparam[1]
            bs=isoparam[2]
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
        self.control='stress'
        
       


    def plotter_totalvals(self):
        totalstrain=self.strain
        totalstress=self.stress
       
        "plott total stress and strain"
 
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
        
        "plot total stress and strain"

        plt.plot(cyclecount,meanstress)
        
        plt.xlabel('Cycle')
        plt.ylabel('Mean Stress (MPa)')
        
        plt.show()
        
        self.meanstress=meanstress
        self.cyclecount=cyclecount
        
    def plotter_peakstrain(self):
        steps=self.steps
        
        totalstrain=self.strain
        totalstress=self.stress
        
        "calculate peak strain"
        maxstrain=totalstrain[steps+1::2*(steps+1)]
        
        
        peake=[]
        for i in range(len(maxstrain)-1):
            peake.append(maxstrain[i]*100)
        
        cyclecount=range(1,len(maxstrain))
        
        "plot total stress and strain"

        plt.plot(cyclecount,peake)
        
        plt.xlabel('Cycle')
        plt.ylabel('Peak Strain, %')
        
        plt.show()
        
        self.peake=peake
        self.cyclecount=cyclecount
        
    def dataextract(self):
        "check to see which data has been asked for"
        #peake=self.peake
        #cyclecount=self.cyclecount
        #meanstress=self.meanstress
        totalstrain=np.asarray(self.strain).astype(np.float)
        totalstress=np.asarray(self.stress).astype(np.float)
        
        workbook = xlwt.Workbook()  
  
        sheet1 = workbook.add_sheet("Total_stress_strain") 
        style = xlwt.easyxf('font: bold 1') 
  
       
        sheet1.write(0, 0, 'Strain', style) 
        sheet1.write(0,1, 'Stress (MPa)', style)
        for i in range(0,len(totalstress)):
            
            sheet1.write(1+i,0, totalstrain[i])
            sheet1.write(1+i,1,totalstress[i])
            workbook.save("simulated_data.xls") 
            
        "check to see if mean stress present"
        if 'meanstress' in dir(self):
            sheet2=workbook.add_sheet("Mean_Stress") 
            cyclecount=np.asarray(self.cyclecount).astype(np.float)
            meanstress=np.asarray(self.meanstress).astype(np.float)
            sheet2.write(0, 0, 'Cycle', style) 
            sheet2.write(0,1, 'Mean Stress (MPa)',style)
            
            for i in range(0,len(meanstress)):
                sheet2.write(1+i,0, cyclecount[i])
                sheet2.write(1+i,1,meanstress[i])
                workbook.save("simulated_data.xls")
                
        "check to see if peak strain present"  
        if 'peake' in dir(self):
            sheet3=workbook.add_sheet("Peak_Strain") 
            cyclecount=np.asarray(self.cyclecount).astype(np.float)
            peake=np.asarray(self.peake).astype(np.float)
            sheet3.write(0, 0, 'Cycle', style) 
            sheet3.write(0,1, 'Peak Strain', style)
            
            for i in range(0,len(peake)):
                sheet3.write(1+i,0, cyclecount[i])
                sheet3.write(1+i,1,peake[i])
                workbook.save("simulated_data.xls")

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
        control=self.control
      
  
        
        "initialise values"
        "netwon raphson parameters"
        max_it=int(1e5)
        toler = 1e-6
        depsp=0
        epsp=0
        sig=0
        epspprev=0
        plasprev=0
        sigend=0
        estrainend=0
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
        #etotal=np.zeros(10000)
        etotal=0
        el=np.zeros(10000)
    
        stressinc=np.diff(lcnd)
        
        
        "create a loop to loop through the increments of strain"
        "add pogress bar"
        pbar = tqdm(total=len(stressinc))
        for i in range(0,len(stressinc)):
            pbar.update((len(stressinc))/(len(stressinc)))
            #pbar.update(int(((i+1)/(len(estrain)-1))*100))
            "starting stress"
            stressstart=lcnd[i+1]
            
        
            "increment of stress based on loading"
            dsig=stressinc[i]
            
            
            "update values at the end of previous branch to the start of the current branch of loading"
            if ktype=='MAF':
                xback=copy.deepcopy(xbackprev)
            elif ktype=='MAFM':
                xback=copy.deepcopy(xbackprev)
                xbackstar=copy.deepcopy(xbackstarprev)
            sig=sigend
            etotal= estrainend
            rs=rsprev
        
            "loading direction provided by sign of the strain increment"
            nu=np.sign(stressinc[i])
        
      
        
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
                etotal=etotal + stressinc[i]/xmod
                sig=sig + stressinc[i]
                continue
              
            "caclulate the stress and strain at yield point"
        
            etotal = etotal + (lam/xmod)*stressinc[i]
            sig = sig + lam*stressinc[i]    
            
            enew.append(etotal)
            snew.append(sig)
        
            de=(stressstart-sig)/steps
            
            
            for k in range(1,steps+1):
                "develop increments of stress starting from previous total stress"
                if (k==steps):
                    el=lcnd[i+1]
                else:
                    el=sig+de*k
                    
                inc[i] = inc[i]+ 1
                    
                "call on the newton-raphson method to solve for the increment of plastic strain"
                if ktype=='MAF':
                    epsp,plas,xback,rs=newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,max_it,toler,Qs,bs,control,xmod,sigy0)
                
                
                    xbackprev=copy.deepcopy(xback)
                
                elif ktype=='MAFM':
                    epsp,plas,xback,rs,xbackstar=newtraphsonMAFM(a,c,am,cm,astar,cstar,MAFMnum,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,xbackstarprev,xbackstar,max_it,toler,Qs,bs,control,xmod,sigy0)
                                             
       
                 
                    xbackprev=copy.deepcopy(xback)
                    xbackstarprev=copy.deepcopy(xbackstar)
                epspprev=epsp
                rsprev=rs
                plasprev=plas
                
                "calculate the new stress for the current increment of strain"    
                enew.append(epsp + el/xmod)
                snew.append(el)
           

            sigend=el
            estrainend=epsp + el/xmod
            
            self.strain=enew
            self.stress=snew
        pbar.close()

"newton-raphson method MAFM"
def newtraphsonMAFM(a,c,am,cm,astar,cstar,MAFMnum,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,xbackstarprev,xbackstar,max_it,toler,Qs,bs,control,xmod,sigy0):
   

    for n in range(0,max_it):
        if control=='stress':
            sig=el
        else:
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
            
            dxback=dbackMAFM(a,c,am,cm,astar,cstar,xbackprev,xbackstarprev,plas,nu,back,xbackstar,epsp,epspprev,MAFMnum)
            #dxback=dbackMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstar,xbackstarprev)
            
            dxsum=np.sum(dxback)
            drs=diso(Qs,bs,rsprev,nu)
            
            #drs=disoback(Qs,bs,rsprev,plas,plasprev)
            dfunc = nu*(-xmod - dxsum)-drs
            depsp=-func/dfunc
            
           
            epsp += depsp
            plas += nu*depsp
            
            "update backstress using new plastic strain increment"
            xback,xbackstar=backMAFM(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstar,xbackstarprev,xback)
            #xback,xbackstar=backMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstar,xbackstarprev,xback)
            
            
            "update isotropic hardening using new plastic strain increment"
            rs= iso(Qs,bs,plas,plasprev,rsprev)
            #rsback=iso(Qs,bs,plas,plasprev,rsprev)
            
    return epsp,plas,xback,rs,xbackstar

"newton-raphson method"
def newtraphson(a,c,epsp,plas,nu,el,xbackprev,inc,xback,depsp,epspprev,plasprev,rs,rsprev,max_it,toler,Qs,bs,control,xmod,sigy0):


    for n in range(0,max_it):
        if control=='stress':
            sig=el
        else:
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

"nonlinear backstress calculation using forward Euler integration"

def back(a,c,xbackprev,plas,epsp,plasprev,epspprev,nu):
    xback=(a*(epsp-epspprev) -c*nu*(xbackprev)*(epsp-epspprev))+xbackprev
    #xback=a*(epsp) -c*(xbackprev)*(epsp)
    return xback

"nonlinear MAFM backstress calculation using forward Euler integration"

def backMAFM(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstarprev,xbackstar,xback):
    "classic basktresses"
    xback[0:len(a)]=(a*(epsp-epspprev) -c*nu*(xbackprev[0:(len(xbackprev)-MAFMnum)])*(epsp-epspprev))+xbackprev[0:(len(xbackprev)-MAFMnum)]
    xbackstar=cstar*(astar - nu*xbackstarprev)*(epsp-epspprev) + xbackstarprev
    xback[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]=(((cm + cstar*(astar - nu*xbackstarprev))*(am - nu*xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]))*(epsp-epspprev)) + xbackprev[(len(xbackprev)-MAFMnum):len(xbackprev)+MAFMnum]

    return xback,xbackstar

def backMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstarprev,xbackstar,xback):
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

def dbackMAFMRK4(a,c,am,cm,astar,cstar,xbackprev,plas,epsp,plasprev,epspprev,nu,MAFMnum,xbackstarprev,xbacksta):
   
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
def dbackMAFM(a,c,am,cm,astar,cstar,xbackprev,xbackstarprev,plas,nu,back,xbackstar,epsp,epspprev,MAFMnum):

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
 


            
