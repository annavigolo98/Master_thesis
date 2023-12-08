import numpy as np
import math as m
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize as opt
import scipy
import ROOT
#ciao prova1

from iminuit import Minuit,cost
import uncertainties

#needs total stopping power files to work   prova1


graphStopPowerNitrogen = ROOT.TGraph( "e_tot_stopping_N.txt" )
graphStopPowerTantalum = ROOT.TGraph( "e_tot_stopping_Ta.txt" )
graphStopPowerArgon = ROOT.TGraph( "e_tot_stopping_Ar.txt" )
graphStopPower15N = ROOT.TGraph( "e_tot_stopping_15N.txt" )

#yield profile txt file
graphYield=ROOT.TGraph("Yield_TaN_dec_1_simulation.txt")

#layer composition: thickness in at/cm^2, percentage of nitrogen, percentage of tantalum 
#COMP LAYERS dec_1
#layer_dx=[0.547524, 0.140000, 0.100000, 0.100000 ] #10^18at/cm^2
#percN=[55.0132, 41.2155 , 12.0000 ,  3.1721  ]
#percTa=[100.-55.0132,100.-41.2155,100.-12.0000,100.-3.1721]


#COMP LAYERS dec_1 prova 1
#layer_dx=[0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ] #10^18at/cm^2
#percN=[42.6,50.1,51.9,52.2,59.7,48.3,41.4,24.4,9.7,3.7 ]
#percTa=[]


#COMP LAYERS dec_1 prova
layer_dx=[0.887524/30.,0.887524*2./30.,0.887524/10.,0.887524/10.,0.887524/10. ,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ] #10^18at/cm^2
percN=[31.1,47.6,50.1,51.9,52.2,59.7,48.3,41.4,24.4,9.7,3.7 ]
percTa=[]

for i in range(len(percN)):
   percTa.append(100.-percN[i])


  
#COMP LAYERS dec_2
   
#layer_dx=[0.100234,0.134804,0.103019,0.229821,0.115302,0.104322,0.111275] #10^18at/cm^2
#percN=[72.6546,53.7575, 53.6916, 47.4629,34.8149, 22.8344, 11.9323    ]
#percTa=[100.-72.6546,100.-53.7575,100.-53.6916,100.-47.4629,100.-34.8149,100.-22.8344, 100.-11.9323]
  
#COMP LAYERS mar_3
#layer_dx=[0.10001,0.27856,0.102153,0.375279,0.202666 ] #10^18at/cm^2
#percN=[29.2815, 42.9672, 32.3414, 20.0711, 6.8396   ]
#percTa=[100.-29.2815,100.-42.9672,100.-32.3414,100.-20.0711,100.-6.8396]























def stopping_eff(ind,E_x): #ind=index of layer to use, E_x energy for which we want to calculate the stopping power and the effective stopping power.

  
   e_totN14 = graphStopPowerNitrogen.Eval( E_x ) 
   e_totTa  = graphStopPowerTantalum.Eval( E_x ) 
   
   
   e_eff_TaN= e_totN14 + (percTa[ind]/percN[ind])* e_totTa 

   return(e_eff_TaN)




#introduce beam straggling effects:  makes a convolution between the yield profile and a gaussian for each point of the yield profile, 
#given the energy values and the parameter used to scale the standard deviation of the gaussian 
def straggling2( x_, par  ): 
    
    conv=np.zeros(len(x_),float)
    sigma=np.zeros(len(x_),float)
    
    for i in range(len(x_)):
      
      
      sigma1   = par*np.sqrt(x_[i]-278.+0.001)
      
      sigma[i]=sigma1

      if (x_[i]>278. and x_[i]<319.5):   #straggling starts from the resonance energy 
        E_inf=x_[i]-20.
        E_sup=x_[i]+20. 
       
        INT= integrate.quad(lambda x: (1./(np.sqrt(2.*np.pi)*sigma1))*np.exp( -(x-x_[i])**2 / ( sigma1*sigma1*2 ) )*graphYield.Eval(x),E_inf ,E_sup)
        conv[i]=INT[0]    
      else:
        conv[i]=graphYield.Eval(x_[i])
    return conv,sigma


    
    
    
#MAIN FUNCTION      
a=14./15. #conversion factor from laboratory to CM reference frame

Ep=np.arange(275.,300.,0.1) 
Nlayer=len(layer_dx) # number of layers
E_0=275. #minimum energy
    


#fnameYield='Yield_TaN_dec_1_simulation.txt'
#Ep, Yield = np.genfromtxt(fnameYield,dtype='float',comments='#',usecols=(0,1),unpack=True) 

#fnameYield='Yield_TaN_dec_2_simulation.txt'
#Ep, Yield = np.genfromtxt(fnameYield,dtype='float',comments='#',usecols=(0,1),unpack=True) 

    
fnameYield='Yield_TaN_dec_1_simulation.txt'
Ep, Yield = np.genfromtxt(fnameYield,dtype='float',comments='#',usecols=(0,1),unpack=True)    

#standard deviation array for the convolution and convoluted yield

sigma=np.zeros(len(Yield),float) 

Yield2=np.zeros(len(Ep),float)

#experimental data

fnameexpdata='yieldexp2.txt' #data 27 mar
Eexp, Yieldexp, Ysigmaexp = np.genfromtxt(fnameexpdata,dtype='float',comments='#',usecols=(0,1,2),unpack=True) 
s=0.3
Yield2, sigma = straggling2( Ep,s)


res=np.zeros(len(Yieldexp),float)

#function for chi-square minimization
def chi2(a):
   
   chitot=0.
   chi2=np.zeros(len(Yieldexp),float)
   
   for i in range( len(Yieldexp) ):
        idx = (np.abs(Ep - Eexp[i])).argmin()
        res[i]=Yieldexp[i] - Yield2[idx]*a
        chi2[i] = pow( ( Yieldexp[i] - Yield2[idx]*a )/Ysigmaexp[i], 2 )
        chitot+=chi2[i]
           
   print("Chiquadro tot ",chitot,'\n')
   print("        \n")
   return (chitot);






a=5./(1.44*1e11) #costant to normalize simulated data with experimental data (to be fitted)
 
m = Minuit(chi2,a)
m.simplex()
m.migrad()  # run optimiser
m.hesse()   # run covariance estimator


a=m.values[0]

print('a  ',a,'\n')
print(m.values)
print(m.errors)




#plot residuals of the fit
plt.errorbar(Eexp,res,Ysigmaexp,fmt='None',color='black',marker='o')



plt.xlabel('E_p (keV)')
plt.ylabel('Y_{exp}-Y_{sim} (1e-6 atoms)')
plt.show()

#plot values
plt.plot(Ep,Yield2*a)
plt.plot(Ep,Yield*a)
plt.errorbar(Eexp,Yieldexp,Ysigmaexp,fmt='None',color='green',marker='s')



plt.xlabel('E_p (keV)')
plt.ylabel('Yield (1e-6 atoms)')
plt.show()

#plot standard deviation used for the straggling
plt.scatter(Ep,sigma,marker='o',s=1)


plt.xlabel('E_p (keV)')
plt.ylabel('sigma straggling keV')
plt.show()

    

#writes in a txt file
fout='Yield_TaN_dec_1_simulation_stragg.txt'
f=open(fout,'w')
f.write('#Implanted_target_1 (dec_1) \n')
f.write('#Energy keV, Yield(10^-6 atoms) \n')
for i in range(len(Ep)): 
             f.write(str(Ep[i])+'  '+str(Yield2[i])+'\n')
f.close()    
    
    
fout='sigma_conv_TaN_dec_1_simulation_stragg.txt'
f=open(fout,'w')
f.write('#Implanted_target_1 (dec_1) \n')
f.write('#Energy keV, Sigma_convolution \n')
for i in range(len(Ep)): 
             f.write(str(Ep[i])+'  '+str(sigma[i])+'\n')
f.close()    
        
    
    
    
    
    
    
