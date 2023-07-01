import numpy as np
import math as m
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize as opt
import ROOT

from iminuit import Minuit,cost
import uncertainties 
from uncertainties import ufloat


#Program to fit the experimental data with theoretical formulas to obtain the peak and total efficiency curve of a Ge detector. Ratios between the net counts under the peaks of the energy spectra and 
#the expected number of counts for the Co, Cs and Ba radioactive sources is used. Ratio between the primary and secondary transition yields of the 15O nucleus are also fitted. Total efficiency of the Cs source is 
#used.


#Branching ratios to use 

b662=0.851
sigmab662=0.02

b662float= ufloat(b662,sigmab662) 


bBR1173=0.998720096
sigmabBR1173=0.00000005488406

bBR2505=1.99744815399133e-08

sigmabBR2505=2.99662685163321e-06


bBR1173float= ufloat(bBR1173,sigmabBR1173)



BR1332=1.
sigmaBR1332=0.

BR1332float= ufloat(BR1332,sigmaBR1332)

a1BR356=0.63196
sigmaa1BR356=0.017934


a1BR356float= ufloat(a1BR356,sigmaa1BR356)

a1BR276=0.075152
sigmaa1BR276=0.0023912


a1BR276float= ufloat(a1BR276,sigmaa1BR276)

a1BR53=0.14518
sigmaa1BR53=0.024766


a1BR53float= ufloat(a1BR53,sigmaa1BR53)

a2BR302=0.09686
sigmaa2BR302=0.0012615


a2BR302float= ufloat(a2BR302,sigmaa2BR302)

a2BR383=0.045385
sigmaa2BR383=0.0010005


a2BR383float= ufloat(a2BR383,sigmaa2BR383)


BR53=0.17
sigmaBR53=0.029

BR53float= ufloat(BR53,sigmaBR53)



BR276=0.088
sigmaBR276=0.0028

BR276float= ufloat(BR276,sigmaBR276)



BR356=0.74
sigmaBR356=0.021

BR356float= ufloat(BR356,sigmaBR356)


BR223=0.018
sigmaBR223=0.0010

BR223float= ufloat(BR223,sigmaBR223)

BR302=0.668
sigmaBR302=0.0087

BR302float= ufloat(BR302,sigmaBR302)

BR383=0.313
sigmaBR383=0.0069

BR383float= ufloat(BR383,sigmaBR383)



BR79=0.88
sigmaBR79=0.019

BR79float= ufloat(BR79,sigmaBR79)


BR160=0.11
sigmaBR160=0.018

BR160float= ufloat(BR160,sigmaBR160)




BR80=1.
sigmaBR80=0.

BR80float= ufloat(BR80,sigmaBR80)


BR763=0.231742380678551

BR1380=0.575043128234618

BR2373=0.15813686026452

BR7554=0.0350776308223117


#energies

E_662=0.661657
E_1173=1.173228
E_1332=1.332501
    
E_276=0.2763989
E_302=0.3028508
E_356=0.3560129
E_383=0.3838485
    
E_53=0.0531622
E_223=0.2232368
E_79=0.0796142
E_80=0.0809979
    
E_160=0.160612
    
    
E_5182=5.182
E_6174=6.1749
E_6791=6.7914
E_763=0.7634
E_1380=1.3801
E_2373=2.373
E_7554=7.5545

E_2505=2.505729




def func(a,b,c,k1,k2,k3):

    #Quantities to be fitted to get the efficiency curve: 
    # net_counts/expected_counts Cs, Co,  Ba, ratio yields of primary to secondary transitions 15O, efftot Cs,  net_counts/expected_counts     sum peak 2505keV Co
    
    y=[0.0184935292224484,0.0160318628603875,0.0149047791402305,    0.00163673927231275,0.00409705178262663,0.0136803518467206,0.00193083363308166,   4.84416414487933,3.26210127533426,2.1554788894365,   0.0202542116329361    , 0.000323278497210588  ]   
    y=np.array(y)

    sigmay=[0.000234827870379168,0.000171066043808209,0.000160966323226344,5.1719202605436e-05,0.000124775008543149,0.000412011695708118,5.94355075089929E-05,            0.226700812445768,    0.0930221298414698, 0.269550934929455,   0.000253878146240374   ,    1.6701023052691e-05]
    sigmay=np.array(sigmay) 
    
    sigmatot=np.zeros(len(y),float)
    
    
        
    val=np.zeros(len(y),float)
    
    val[0]=b662*eff_peak(a,b,c,E_662)
    
    #temp0=eff_peak(a,b,c,E_662)*sigmab662
    
    #sigmatot[0]=np.sqrt(  temp0**2. + sigmay[0]**2. )
    
    temp0=b662float*eff_peak(a,b,c,E_662)
    
    sigmatot[0]=np.sqrt(  uncertainties.std_dev(temp0)**2. + sigmay[0]**2. )
    
    
    
    val[1]=bBR1173*eff_peak(a,b,c,E_1173)*(1.-eff_tot(eff_peak(a,b,c,E_1332),k1,k2,k3,E_1332)*BR1332)
    
    #temp1=  eff_peak(a,b,c,E_1173)*(1.-eff_tot(eff_peak(a,b,c,E_1332),k1,k2,k3,E_1332)*BR1332 )*sigmabBR1173 
    
    #sigmatot[1]=np.sqrt(   temp1**2. + sigmay[1]**2.  )
    
    temp1=bBR1173float*eff_peak(a,b,c,E_1173)*(1.-eff_tot(eff_peak(a,b,c,E_1332),k1,k2,k3,E_1332)*BR1332)
    
    sigmatot[1]=np.sqrt(   uncertainties.std_dev(temp1)**2. + sigmay[1]**2.  )
    
    
    
    
    val[2]=bBR1173*BR1332*eff_peak(a,b,c,E_1332)*(1-eff_tot(eff_peak(a,b,c,E_1173),k1,k2,k3,E_1173))
    
    #temp2=eff_peak(a,b,c,E_1332)*eff_tot(eff_peak(a,b,c,E_1173),k1,k2,k3,E_1173)*BR1332*sigmabBR1173
    
    #sigmatot[2]=np.sqrt(  temp2**2. + sigmay[2]**2.  )
    
    temp2=bBR1173float*BR1332*eff_peak(a,b,c,E_1332)*(1-eff_tot(eff_peak(a,b,c,E_1173),k1,k2,k3,E_1173))
    
    sigmatot[2]=np.sqrt(  uncertainties.std_dev(temp2)**2. + sigmay[2]**2.  )
    
    
    
    val[3]=a1BR276*eff_peak(a,b,c,E_276)*(1.-eff_tot(eff_peak(a,b,c,E_160),k1,k2,k3,E_160)*BR160 - eff_tot(eff_peak(a,b,c,E_79),k1,k2,k3,E_79)*BR79- eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80*BR79)
    
    #temp3=eff_peak(a,b,c,E_276)*np.sqrt( ( (1.-eff_tot(eff_peak(a,b,c,E_160),k1,k2,k3,E_160)*BR160 - eff_tot(eff_peak(a,b,c,E_79),k1,k2,k3,E_79)*BR79 - eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR79*BR80)*sigmaa1BR276 )**2.   +   ( a1BR276*eff_tot(eff_peak(a,b,c,E_160),k1,k2,k3,E_160)*sigmaBR160    )**2  +( a1BR276*eff_tot(eff_peak(a,b,c,E_79),k1,k2,k3,E_79)*sigmaBR79      )**2  +   ( a1BR276*eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80*sigmaBR79       )**2                                  )
    
    #sigmatot[3]=np.sqrt(  temp3**2. + sigmay[3]**2.  )
    
    temp3=a1BR276float*eff_peak(a,b,c,E_276)*(1.-eff_tot(eff_peak(a,b,c,E_160),k1,k2,k3,E_160)*BR160float - eff_tot(eff_peak(a,b,c,E_79),k1,k2,k3,E_79)*BR79float- eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80*BR79float)
    
    sigmatot[3]=np.sqrt(  uncertainties.std_dev(temp3)**2. + sigmay[3]**2.  ) 
    
    val[4]=a1BR53*BR302*eff_peak(a,b,c,E_302)*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53)- eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80) + a2BR302*eff_peak(a,b,c,E_302)*(1-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80) 
    
    #temp4=np.sqrt(  ( BR302*eff_peak(a,b,c,E_302)*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53)-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80)*sigmaa1BR53    )**2.  +( a1BR53*eff_peak(a,b,c,E_302)*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53)-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80)*sigmaBR302  )**2.   +( eff_peak(a,b,c,E_302)*(1.-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80)*sigmaa2BR302        )**2.                        )
    
    #sigmatot[4]=np.sqrt(  temp4**2. + sigmay[4]**2.  )  
    
    temp4=a1BR53float*BR302float*eff_peak(a,b,c,E_302)*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53)- eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80) + a2BR302float*eff_peak(a,b,c,E_302)*(1-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80) 
    
    sigmatot[4]=np.sqrt(  uncertainties.std_dev(temp4)**2. + sigmay[4]**2.  ) 
    
    
    val[5]=a1BR356*eff_peak(a,b,c,E_356)*(1.-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80)
    
    #temp5=  eff_peak(a,b,c,E_356)*(1-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80)*sigmaa1BR356
    
    #sigmatot[5]=np.sqrt(  temp5**2. + sigmay[5]**2.  )
    
    temp5=a1BR356float*eff_peak(a,b,c,E_356)*(1.-eff_tot(eff_peak(a,b,c,E_80),k1,k2,k3,E_80)*BR80)
    
    sigmatot[5]=np.sqrt(  uncertainties.std_dev(temp5)**2. + sigmay[5]**2.  ) 
    
    
    
    val[6]=a1BR53*BR383*eff_peak(a,b,c,E_383)*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53)) + a2BR383*eff_peak(a,b,c,E_383)
    
    #temp6=eff_peak(a,b,c,E_383)*np.sqrt( (BR383*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53))*sigmaa1BR53  )**2. +  (a1BR53*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53))*sigmaBR383  )**2.  + ( sigmaa2BR383   )**2                                   )
    
    #sigmatot[6]=np.sqrt(  temp6**2. + sigmay[6]**2.  ) 
    
    temp6=a1BR53float*BR383float*eff_peak(a,b,c,E_383)*(1.-eff_tot(eff_peak(a,b,c,E_53),k1,k2,k3,E_53)) + a2BR383float*eff_peak(a,b,c,E_383)
    
    sigmatot[6]= np.sqrt(  uncertainties.std_dev(temp6)**2. + sigmay[6]**2.  ) 
    
    val[7]=(eff_peak(a,b,c,E_763)*(1.-eff_tot(eff_peak(a,b,c,E_6791),k1,k2,k3,E_6791)) ) / (eff_peak(a,b,c,E_6791)*(1.-eff_tot(eff_peak(a,b,c,E_763),k1,k2,k3,E_763)))
    
    sigmatot[7]=sigmay[7]
    
    val[8]=(eff_peak(a,b,c,E_1380)*(1.-eff_tot(eff_peak(a,b,c,E_6174),k1,k2,k3,E_6174))) / (eff_peak(a,b,c,E_6174)*(1.-eff_tot(eff_peak(a,b,c,E_1380),k1,k2,k3,E_1380)))
    
    sigmatot[8]=sigmay[8]
    
    val[9]=(eff_peak(a,b,c,E_2373)*(1.-eff_tot(eff_peak(a,b,c,E_5182),k1,k2,k3,E_5182))) / (eff_peak(a,b,c,E_5182)*(1.-eff_tot(eff_peak(a,b,c,E_2373),k1,k2,k3,E_2373)))
    
    sigmatot[9]=sigmay[9]
    
    val[10]=eff_tot(eff_peak(a,b,c,E_662),k1,k2,k3,E_662)
    
    sigmatot[10]=sigmay[10]
    
    
    val[11]=bBR1173*BR1332*eff_peak(a,b,c,E_1173)*eff_peak(a,b,c,E_1332)
    
    #temp11= sigmabBR1173*BR1332*eff_peak(a,b,c,E_1173)*eff_peak(a,b,c,E_1332)  
    
    #sigmatot[11]=np.sqrt(   temp11**2. + sigmay[11]**2.  )
    
    temp11=bBR1173float*BR1332*eff_peak(a,b,c,E_1173)*eff_peak(a,b,c,E_1332)
    
    sigmatot[11]=np.sqrt(  uncertainties.std_dev(temp11)**2. + sigmay[11]**2.  )
    
    
    
    
    
    chi2=0.
    
    chi2par=np.zeros(len(y),float)
    
    
    
    
    #chiquadro yield
    for i in range(len(val)):
    
       chi2par[i]= pow((y[i]-val[i])/sigmatot[i],2.);
       chi2+=chi2par[i];
      
         
    print("Val calculated  " , val,'\n')
    print("        \n")
    print("Sigma yield exp  ",sigmay,'\n')
    print("         \n")
    print('Total sigma ',sigmatot,'\n')
    print("        \n")    
    print("Chiquadro  ",chi2,'\n')
    print("        \n")
    return (chi2);
    
    
    

   

#calculates peak efficiency
def eff_peak(a,b,c,Ex):
   d0=2.0
   d1=5.4
   d2=0.1
   d=1.35
   
   D=(1.-np.exp(-( (d+d0)/(d1+d2*np.sqrt(Ex))  )))/(d+d0)**2.
   
   A=0.0379
   eff_peak=A*np.exp(a+b*np.log(Ex)+c*np.log(Ex)**2)
   
   return(eff_peak)
   
#calculates total efficiency   
def eff_tot(eff_p,k1,k2,k3,Ex):
   eff_tot=eff_p*np.exp(-(k1+k2*np.log(Ex)+k3*np.log(Ex)*np.log(Ex)))  
   
   return(eff_tot)







 
#parameters J
#parameters=[0.08,-0.57,-0.103,-1.47,-0.6,-0.1] 


#parameters D
parameters=[0.099169882251258,-0.587111677130061,-0.0944993174452277,-1.37134385962633,-0.554341534254716,-0.0433382652728538]

parameters=np.array(parameters)

#lsq=opt.least_squares(residuals, parameters, args=(x,y), xtol=1e-07, loss='cauchy')


m = Minuit(func,a=parameters[0],b=parameters[1],c=parameters[2],k1=parameters[3],k2=parameters[4],k3=parameters[5])

#m.limits = [(0.030,0.2), (-0.7,-0.5),(-0.340,-0.042),(-1.6,-1.2),(-0.580,-0.300),(-0.190,-0.04)]

#m.fixed["k1"] = True
#m.fixed["k2"] = True
#m.fixed["k3"] = True

m.migrad()  # run optimiser
m.hesse()   # run covariance estimator


print(m.values) 
print(m.errors)  

a=m.values[0]
b=m.values[1]
c=m.values[2]
k1=m.values[3]
k2=m.values[4]
k3=m.values[5]

 
prova=func(a,b,c,k1,k2,k3)

   
   
   

