import numpy as np
import math as m
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import optimize as opt
import scipy.integrate as integrate
import ROOT


#stopping power

graphStopPowerNitrogen14 = ROOT.TGraph( "e_tot_stopping_N14.txt" )
graphStopPowerTantalum = ROOT.TGraph( "e_tot_stopping_Ta.txt" )
graphStopPowerNitrogen15 = ROOT.TGraph( "e_tot_stopping_N15.txt" )
graphStopPowerZirconium = ROOT.TGraph( "e_tot_stopping_Zr.txt" )

#S-factor and detector total efficiency curevs
graphSfactor=ROOT.TGraph( "s_6793_li.txt" )  # in MeV e MeV/barn
graphEfficiency=ROOT.TGraph( "eff_tot55gradi.txt" )



#calculates efficiency at energy E_x

def calc_efficiency(E_x):
    
    efficiency=graphEfficiency.Eval( E_x )

    return(efficiency)
    
    
    
#writes total stopping power of an element

def write_tot_stopping(e_tot_stopping,E_stopping,name):    
    
    fout='n.txt'
    if (name=='N14'):
      fout='e_tot_stopping_N14.txt'
      
    if (name=='N15'):
      fout='e_tot_stopping_N15.txt'
        
    if (name=='Ta'):
      fout='e_tot_stopping_Ta.txt'
    
    if (name=='Zr'):
      fout='e_tot_stopping_Zr.txt'
      
        
    f=open(fout,'w')
    
    
    
    for i in range(len(E_stoppingN14)): 
        f.write(str(E_stopping[i])+'  '+str(e_tot_stopping[i])+str('\n'))
               
    f.close()
    



#functions to calculate the stopping power at a given energy E_x 


#TaN targets

def stopping_true_TaN(E_x): 
   
   e_totN14 = graphStopPowerNitrogen14.Eval( E_x ) #stopping power totale a E_x per N14
   e_totTa  = graphStopPowerTantalum.Eval( E_x ) #stopping power totale a E_x per Ta
   
   e_true_TaN=(e_totN14+e_totTa)/2. 
   
   return(e_true_TaN)

#ZrN targets
   
def stopping_true_ZrN(E_x): 
   
   e_totN14 = graphStopPowerNitrogen14.Eval( E_x ) #stopping power totale a E_x per N14
   e_totZr  = graphStopPowerZirconium.Eval( E_x ) #stopping power totale a E_x per Ta
   
   e_true_ZrN=(e_totN14+e_totZr)/2. 
   
   return(e_true_ZrN)   

# TaN targets

def stopping_eff_TaN(E_x): 
   
   e_totN14 = graphStopPowerNitrogen14.Eval( E_x )
   e_totN15 = graphStopPowerNitrogen15.Eval( E_x ) 
   e_totTa  = graphStopPowerTantalum.Eval( E_x ) 
   
   e_eff_TaN= e_totN14+0.004*e_totN15+e_totTa
   

   return(e_eff_TaN)
   
   
# ZrN targets

def stopping_eff_ZrN(E_x): 
   
   e_totN14 = graphStopPowerNitrogen14.Eval( E_x )
   e_totN15 = graphStopPowerNitrogen15.Eval( E_x ) 
   e_totZr  = graphStopPowerZirconium.Eval( E_x ) 
   
   e_eff_ZrN= e_totN14+0.004*e_totN15+e_totZr
   

   return(e_eff_ZrN)
   
   



#function to calculate the energy loss in a given target

#TaN Targets


def deltae_TaN(E_x): 
  
  dx=1. #10^18 atoms/cm^2 
  
  
  e_true1=stopping_true_TaN(E_x) #stopping power at E_x
  DE=e_true1*dx #DE a E_x
  
  E_def=E_x-DE/2.
  e_true_def=stopping_true_TaN(E_def) #stopping power at E_x-DE/2
  
  DE=e_true_def*dx # DE a E_x-DE/2
  
  
  return(DE)  
  
  
#ZrN Targets


def deltae_ZrN(E_x): 
  
  dx=1. #10^18 atoms/cm^2
  
  
  e_true1=stopping_true_ZrN(E_x) #stopping power at E_x
  DE=e_true1*dx #DE a E_x
  
  E_def=E_x-DE/2.
  e_true_def=stopping_true_ZrN(E_def) #stopping power at E_x-DE/2
  
  DE=e_true_def*dx # DE a E_x-DE/2
  
  
  return(DE)   
   
#function to integrate the cross section divided by the effective stopping power to calculate the yield profile

#TaN

   
def integral_TaN(E_1,E_2):  #extrema of energy in the laboratory frame
   a=14./15.
   
  
   EMax = E_2 #E_2=E
   EMin = E_1 #E_1=E-DE

   integral = 0.
    
  
   nSteps = 100000

  
   step = (EMax - EMin)/nSteps

  
   E_step = EMin + step/2  # trapezoid
   #E_step = EMin           # classic method
   
   
   
   def calculateCrossSection(E_x):
     
     S_factor= graphSfactor.Eval(  E_x*1e-3 )*1e3 #keVbarn
     cross=  (S_factor/E_x)*np.exp(-212.4/(np.sqrt( E_x )))
     return(cross)
 
     
     
   for i in range( nSteps ):
       stopPower     = stopping_eff_TaN( E_step ) #stopping power in E_step
       
       crossSection  = calculateCrossSection( E_step*a )
       integral     += step*crossSection/stopPower
       E_step        += step
 
   
   
   return(integral)
   
   
#ZrN

def integral_ZrN(E_1,E_2): 
   a=14./15. 

   
   EMax = E_2 #E_2=E
   EMin = E_1 #E_1=E-DE

   
   integral = 0.
   
   nSteps = 100000

   step = (EMax - EMin)/nSteps

  
   E_step = EMin + step/2  # trapezoid
   #E_step = EMin           # classic method
   
   def calculateCrossSection(E_x):
     
     S_factor= graphSfactor.Eval(  E_x*1e-3 )*1e3 #keVbarn
     cross=  (S_factor/E_x)*np.exp(-212.4/(np.sqrt( E_x )))
     return(cross)
     
     
   for i in range( nSteps ):
       stopPower     = stopping_eff_ZrN( E_step ) #stopping power in E_step
       
       crossSection  = calculateCrossSection( E_step*a )
       integral     += step*crossSection/stopPower
       E_step        += step

      
   
   
   return(integral)   

      
    
#MAIN FUNCTION    
   

#6.793MeV transition
fname='s_6793_li.txt'
E_cm, S = np.genfromtxt(fname,dtype='float',comments='#',usecols=(0,1),unpack=True)
E_cm=E_cm*1e3 #energy cm in keV
S=S*1e3 #S in keVbarn


    


#plot S factor vs E_cm
plt.xlabel('E_cm MeV')
plt.ylabel('S(E) MeV barn')
plt.plot(E_cm*1e-3,S*1e-3) #plot in MeV MeV/barn
plt.show()



sigma=np.zeros(len(E_cm),float)
sigma[:]=(S[:]/E_cm[:])*np.exp(-212.4/(np.sqrt(E_cm[:])))




a2=15./14. #conversion from E_cm to E_lab

E_lab=np.zeros(len(E_cm),float)
E_lab=a2*E_cm


#plot sigma vs E_cm
plt.xlabel('E_cm keV')
plt.ylabel('sigma(E) barn')
plt.plot(E_cm,sigma)
plt.show()







# stopping power N14, N15, Ta, Zr


#14N
fname14N='H_in_N14.txt'
E_stoppingN14, e_elN14, e_nucN14 = np.genfromtxt(fname14N,dtype='float',comments='#',usecols=(0,2,3),unpack=True)
e_totN14=e_elN14+e_nucN14

#15N

fnameN15='H_in_N15.txt'
E_stoppingN15, e_elN15, e_nucN15 = np.genfromtxt(fnameN15,dtype='float',comments='#',usecols=(0,2,3),unpack=True)
e_totN15=e_elN15+e_nucN15

#Ta

fnameTa='H_in_Ta.txt'
E_stoppingTa, e_elTa, e_nucTa = np.genfromtxt(fnameTa,dtype='float',comments='#',usecols=(0,2,3),unpack=True)
e_totTa=e_elTa+e_nucTa


#Zr

fnameZr='H_in_Zr.txt'
E_stoppingZr, e_elZr, e_nucZr = np.genfromtxt(fnameZr,dtype='float',comments='#',usecols=(0,2,3),unpack=True)
e_totZr=e_elZr+e_nucZr



#writes total stopping power for N14
#name='N14'
#write_tot_stopping(e_totN14,E_stoppingN14,name)

#writes total stopping power for N15
#name='N15'
#write_tot_stopping(e_totN15,E_stoppingN15,name)

#writes total stopping power for Ta
#name='Ta'
#write_tot_stopping(e_totTa,E_stoppingTa,name)

#writes total stopping power for Zr
#name='Zr'
#write_tot_stopping(e_totZr,E_stoppingZr,name)



E_0=E_lab[0]


Y_TaN=[]
Y_ZrN=[]

#TaN

for i in range(len(E_lab)):
    #calcolo dell'integrale per ogni energia del protone
    DE=0. #deltae
    Yield=0.
    
    #estremi integrale all'inizio
    
    ext_sup=E_lab[i] 
    DE=deltae_TaN(ext_sup)
    #print("DE_i ",DE,'\n')
    
    
    ext_inf=E_lab[i]-DE
    
    
    
    if(ext_inf<E_0):
        #integro da E_0 a ext_sup
         
        Yield+=integral_TaN(E_0,ext_sup)
        
        
      
      
    if(ext_inf>=E_0):
        #integro da ext_inf a ext_sup
        Yield+= integral_TaN(ext_inf,ext_sup)
        
       
    Y_TaN.append(Yield)
     

Y_TaN=np.array(Y_TaN)



#ZrN

for i in range(len(E_lab)):
    #calcolo dell'integrale per ogni energia del protone
    DE=0. #deltae
    Yield=0.
    
    #estremi integrale all'inizio
    
    ext_sup=E_lab[i] 
    DE=deltae_ZrN(ext_sup)
    #print("DE_i ",DE,'\n')
    
    
    ext_inf=E_lab[i]-DE
    
    
    
    if(ext_inf<E_0):
        #integro da E_0 a ext_sup
         
        Yield+=integral_ZrN(E_0,ext_sup)
        
        
      
      
    if(ext_inf>=E_0):
        #integro da ext_inf a ext_sup
        Yield+= integral_ZrN(ext_inf,ext_sup)
        
        
    Y_ZrN.append(Yield) 

Y_ZrN=np.array(Y_ZrN)






#plt.plot(E_lab,Y_TaN)
#plt.ylabel('Yield TaN (1e-6 atoms)')
#plt.xlabel('E_lab (keV)')
#plt.yscale('log')
#plt.show()

#plt.plot(E_lab,Y_ZrN)
#plt.ylabel('Yield ZrN (1e-6 atoms)')
#plt.xlabel('E_lab (keV)')
#plt.yscale('log')
#plt.show()

Ie=624150900000000  #I/e con I=100muA

rate_TaN=1e-6*Y_TaN*Ie
rate_ZrN=1e-6*Y_ZrN*Ie

plt.plot(E_lab,rate_TaN)
plt.ylabel('Rate TaN (1/s)')
plt.xlabel('E_lab (keV)')
plt.yscale('log')
plt.show()

plt.plot(E_lab,rate_ZrN)
plt.ylabel('Rate ZrN (1/s)')
plt.xlabel('E_lab (keV)')
plt.yscale('log')
plt.show()

#transizione da RC a 6793 keV state

efficiency=np.zeros(len(E_lab),float)

E_gamma2=6793.0 #keV #fotone secondario
Q=7297.0 #keV
E_gamma=np.zeros(len(E_lab),float)

for i in range(len(E_lab)):
    E_gamma[i]=E_lab[i]*(1./a2)+Q-E_gamma2
    
    
#print(len(E_gamma),'\n')
#print(E_gamma,'\n')
#for i in range(len(E_lab)):
#   efficiency[i]=  calc_efficiency( E_lab[i]   )

fname='eff_tot55gradi.txt'
E_gammadati, eff_gammadati = np.genfromtxt(fname,dtype='float',comments='#',usecols=(0,1),unpack=True) 

#print(E_gammadati,'   ',eff_gammadati,'\n')


for i in range(len(E_gamma)):
   for j in range(len(E_gammadati)-1): # indici j, j+1
       if(E_gamma[i]<E_gammadati[0]):
            
            efficiency[i]=eff_gammadati[0]
            
       if(E_gamma[i]>E_gammadati[j] and E_gamma[i]<=E_gammadati[j+1]):
                
                
            efficiency[i]=eff_gammadati[j]+(E_gamma[i]-E_gammadati[j])*(  (eff_gammadati[j+1]-eff_gammadati[j])/(E_gammadati[j+1]-E_gammadati[j])  )      
            
             
       if(E_gamma[i]>E_gammadati[len(E_gammadati)-1]):
            efficiency[i]=eff_gammadati[len(E_gammadati)-1]
             
plt.scatter(E_gamma,efficiency)
plt.plot(E_gammadati,eff_gammadati)
plt.xlabel('E_gamma (keV)')
plt.ylabel('Efficiency (tot)')
plt.show()
    


#plot efficienza
#plt.scatter(E_gamma,efficiency)
#plt.xlabel('E_gamma (keV)')
#plt.ylabel('Efficiency (tot)')
#plt.show()



expected_rateTaN=efficiency*rate_TaN
expected_rateZrN=efficiency*rate_ZrN

#plt.xlim(0.2,1.4)
plt.title('Expected rate for TaN target (RC-6.793MeV)')
plt.plot(E_lab*1e-3,expected_rateTaN,linewidth=2)
plt.xlabel('E_lab (MeV)')
plt.ylabel('counts/s')
plt.yscale('log')
plt.grid()
plt.show()


#plt.xlim(0.2,1.4)
plt.title('Expected rate for ZrN target (RC-6.793MeV)')
plt.plot(E_lab*1e-3,expected_rateZrN,linewidth=2)
plt.xlabel('E_lab (MeV)')
plt.ylabel('counts/s')
plt.yscale('log')
plt.grid()
plt.show()


fout='rateRC6.79.txt'
f=open(fout,'w')
f.write('#Energy MeV, rate RC-6.793MeV TaN,  rate RC-6.79 MeV ZrN \n')
for i in range(len(E_lab)): 
             f.write(str(E_lab[i]*1e-3)+'  '+str(expected_rateTaN[i])+'   '+str(expected_rateZrN[i])+'\n')
f.close()












