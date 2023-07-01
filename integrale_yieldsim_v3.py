import numpy as np
import math as m
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize as opt
import ROOT


graphStopPowerNitrogen = ROOT.TGraph( "e_tot_stopping_N.txt" )
graphStopPowerTantalum = ROOT.TGraph( "e_tot_stopping_Ta.txt" )


#COMPOSIZIONE LAYERS dec_1
#layer_dx=[0.547524, 0.140000, 0.100000, 0.100000 ] #10^18at/cm^2
#percN=[55.0132, 41.2155 , 12.0000 ,  3.1721  ]
#percTa=[100.-55.0132,100.-41.2155,100.-12.0000,100.-3.1721]




#COMPOSIZIONE LAYERS dec_1 prova 1
#layer_dx=[0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ] #10^18at/cm^2
#percN=[42.6,50.1,51.9,52.2,59.7,48.3,41.4,24.4,9.7,3.7 ]
#percTa=[]

#COMPOSIZIONE LAYERS dec_1 prova2
layer_dx=[0.887524/30.,0.887524*2./30.,0.887524/10.,0.887524/10.,0.887524/10. ,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ] #10^18at/cm^2
percN=[31.1,47.6,50.1,51.9,52.2,59.7,48.3,41.4,24.4,9.7,3.7 ]
percTa=[]


for i in range(len(percN)):
   percTa.append(100.-percN[i])
#COMPOSIZIONE LAYERS dec_1 fit
#layer_dx=[ 0.08930819164619654 , 0.20886529840734203   ] #10^18at/cm^2
#percN=[  28.028700107292423,41.396862419368325   ]
#percTa=[100.-28.028700107292423   ,100.-41.396862419368325  ]
  
#COMPOSIZIONE LAYERS dec_2
   
#layer_dx=[0.100234,0.134804,0.103019,0.229821,0.115302,0.104322,0.111275] #10^18at/cm^2
#percN=[72.6546,53.7575, 53.6916, 47.4629,34.8149, 22.8344, 11.9323    ]
#percTa=[100.-72.6546,100.-53.7575,100.-53.6916,100.-47.4629,100.-34.8149,100.-22.8344, 100.-11.9323]
  
#COMPOSIZIONE LAYERS mar_3
#layer_dx=[0.10001,0.27856,0.102153,0.375279,0.202666 ] #10^18at/cm^2
#percN=[29.2815, 42.9672, 32.3414, 20.0711, 6.8396   ]
#percTa=[100.-29.2815,100.-42.9672,100.-32.3414,100.-20.0711,100.-6.8396]



#funzioni per calcolare lo stopping power a una energia data E_x; interpolazione lineare ROOT 

def stopping_true(ind,E_x): # E_x energia per la quale voglio calcolare stopping power e effective stopping power.

   e_totN14 = graphStopPowerNitrogen.Eval( E_x ) #stopping power totale a E_x per N14
   e_totTa  = graphStopPowerTantalum.Eval( E_x ) #stopping power totale a E_x per Ta
   
   e_true_TaN= 0.01*percN[ind]*e_totN14 + 0.01*percTa[ind]*e_totTa #stopping power true per TaN target
   
   return(e_true_TaN)


def stopping_eff(ind,E_x): #ind=indice di layer da usare, E_x energia per la quale voglio calcolare stopping power e effective stopping power.

  
   e_totN14 = graphStopPowerNitrogen.Eval( E_x ) #stopping power totale a E_x per N14
   e_totTa  = graphStopPowerTantalum.Eval( E_x ) #stopping power totale a E_x per Ta
   
   
   e_eff_TaN= e_totN14 + (percTa[ind]/percN[ind])* e_totTa # effective stopping power in base al layer

   return(e_eff_TaN)

   

#funzione per calcolare deltaE effettivo a partire da un valore di energia nel layer considerato (NB target TaN)


def deltae(ind,E_x): #ind = indice del target che si vuole usare, E_x = energia del protone nel lab
   
  DE=0.
  
  # Calcoliamo gli estremi dell'integrale
  Xmax = layer_dx[ind] 
  Xmin = 0. #E_1=E-DE
   
  # Definiamo il numero di step
  nSteps = 1000

  # Calcoliamo la grandezza dello step
  step = (Xmax - Xmin)/nSteps
   
  E_step=E_x
 
  for i in range( nSteps ):
       stopPower     = stopping_true(ind, E_step) #stopping power in E_step   nel lab
       E_step        -= step*stopPower
  
  #calcolo DE
  DE=E_x-E_step
  
  
  return(DE)


#funzione per calcolare l'integrale tra due estremi (in un layer)


def integral_L(ind,E_1,E_2):  #ind=indice di layer E_1 ed E_2 sono ext_inf e ext_sup dell'integrale nel main


   
   
   a=14./15.
   # Calcoliamo gli estremi dell'integrale nel lab
   EMax = E_2 #E_2=E
   EMin = E_1 #E_1=E-DE

   # Definiamo le variabili ausiliarie per l'integrale
   integral = 0.
    
   # Definiamo il numero di step
   nSteps = 1000

   # Calcoliamo la grandezza dello step
   step = (EMax - EMin)/nSteps

   # Definiamo la variabile energia
   E_step = EMin + step/2  # Metodo Trapezio
   #E_step = EMin           # Metodo Classico
   
   def calculateCrossSection(E):
     
     T=((0.9893)**2.)/4. #gamma of the resonamce in CM frame (da report)
     E_R=259.56 #keV risonanza nel CM (da report)
     
     
     cross=1./((E-E_R)**2.+T)
     return(cross)
     
     
   for i in range( nSteps ):
       stopPower     = stopping_eff(ind, E_step) #stopping power in E_step   nel lab
       crossSection  = calculateCrossSection( E_step*a ) #cross section in CM
       integral     += step*crossSection/stopPower
       E_step        += step

      
   
   
   return(integral)



    
#writes total stopping power of an element

def write_tot_stopping(e_tot_stopping,E_stopping,name):    
    
    if (name=='N'):
      fout='e_tot_stopping_N.txt'
    if (name=='Ta'):
      fout='e_tot_stopping_Ta.txt'
    if (name=='Ar'):
      fout='e_tot_stopping_Ar.txt'
    if (name=='15N'):
      fout='e_tot_stopping_15N.txt'
        
    f=open(fout,'w')
    
    for i in range(len(E_stoppingN14)): 
        f.write(str(E_stopping[i])+'  '+str(e_tot_stopping[i])+str('\n'))
               
    f.close()
    
#writes effective stopping power

def write_eff_stopping(e_totN14,e_totTa,E_stoppingN14):
  
  fout='eff_stoppingTaN_dec_1_simulation.txt'
  f=open(fout,'w')
  
  
  
  for i in range(len(layer_dx)):
  
      e_eff_TaN= e_totN14 + (percTa[i]/percN[i])* e_totTa
      f.write('#Implanted_target_1 (dec_1) \n')
      f.write('#Energy keV, e_effTaN (ev/10^15at/cm^2) \n')
      for i in range(len(E_stoppingN14)): 
             f.write(str(E_stoppingN14[i])+'  '+str(e_eff_TaN[i])+str('\n'))
      f.write("end layer  "+str('\n'))       
  f.close()
     










#MAIN FUNCTION      
# stopping power TaN 


#14N
fname14N='H_in_N14.txt'
E_stoppingN14, e_elN14, e_nucN14 = np.genfromtxt(fname14N,dtype='float',comments='#',usecols=(0,2,3),unpack=True)

e_totN14=e_elN14+e_nucN14

#Ta

fnameTa='H_in_Ta.txt'
E_stoppingTa, e_elTa, e_nucTa = np.genfromtxt(fnameTa,dtype='float',comments='#',usecols=(0,2,3),unpack=True)

e_totTa=e_elTa+e_nucTa

#Ar

fnameAr='H_in_Ar.txt'
E_stoppingAr, e_elAr, e_nucAr = np.genfromtxt(fnameAr,dtype='float',comments='#',usecols=(0,2,3),unpack=True)

e_totAr=e_elAr+e_nucAr


#15N

fnameN15='H_in_N15.txt'
E_stopping15N, e_el15N, e_nuc15N = np.genfromtxt(fnameN15,dtype='float',comments='#',usecols=(0,2,3),unpack=True)

e_tot15N=e_el15N+e_nuc15N


#WRITE STOPPING POWERS
#writes total stopping for N
#name='N'
#write_tot_stopping(e_totN14,E_stoppingN14,name)

#writes total stopping power for Ta
#name='Ta'
#write_tot_stopping(e_totTa,E_stoppingTa,name)


#writes total stopping power for Ar gas
#name='Ar'
#write_tot_stopping(e_totAr,E_stoppingAr,name)



#writes total stopping power for 15N
#name='15N'
#write_tot_stopping(e_tot15N,E_stopping15N,name)


#writes effective  stopping power
#write_eff_stopping(e_totN14,e_totTa,E_stoppingN14)




#
E_p=np.arange(275.,300.,0.1)  # step di 0.5 keV step della convoluzione
Nlayer=len(layer_dx) #NUMERO LAYERS #4;7;5
E_0=275. #energia minima



a=14./15.


Yield1=[]  #Yield value da aggiornare

partial_yield=np.zeros([Nlayer,len(E_p)],float)



for i in range(len(E_p)):
    #calcolo dell'integrale per ogni energia del protone
    DE=0. #deltae
    Yield=0.
    
    #estremi integrale all'inizio
    
    ext_sup=E_p[i] 
    DE=deltae(0,ext_sup)
    #print("DE_0 ",DE,'\n')
    
    
    ext_inf=E_p[i]-DE
    
    
    for j in range(Nlayer):
      
      
      if(ext_inf<0.):
        #integro da 0 a ext_sup
         
        Yield+=integral_L(j,0.,ext_sup)
        partial_yield[j,i]=integral_L(j,0.,ext_sup)
        #print('break1 \n')
        break   #poi esco dal for dei layers
      
      
      if(ext_inf>=0.):
        #integro da ext_inf a ext_sup
        Yield+= integral_L(j,ext_inf,ext_sup)
        partial_yield[j,i]=integral_L(j,ext_inf,ext_sup) 
        
      if(j==Nlayer-1):
        #print('break2 \n')
        break #esce dal while senza calcolare i nuovi estremi perchè sono finiti
        
        
      ext_sup=ext_inf #E-DE   
      DE=deltae(j+1,ext_sup)
      ext_inf=ext_sup-DE
      
    Yield1.append(Yield)
    

Yield1=np.array(Yield1)







#dati sperimentali
fnamedata='yieldexp2.txt'
a=5./(1.44*1e11) #costante per scalare la simulazione ai dati sperimentali 
E_exp, Y_exp, sigmaY_exp = np.genfromtxt(fnamedata,dtype='float',comments='#',usecols=(0,1,2),unpack=True)


for i in range(Nlayer):
  #plot partial yield 
  plt.plot(E_p,partial_yield[i,:]*a,color='black')


#plot valori

plt.errorbar(E_exp,Y_exp,sigmaY_exp,fmt='None',color='green',marker='s')
plt.scatter(E_p,Yield1*a,marker='o',s=1)

plt.xlabel('E_p (keV)')
plt.ylabel('Yield (1e-6 atoms)')
plt.show()


fout='Yield_TaN_dec_1_simulation.txt'
f=open(fout,'w')
f.write('#Implanted_target_1 (dec_1) \n')
f.write('#Energy keV, Yield(10^-6 atoms) \n')
for i in range(len(E_p)): 
             f.write(str(E_p[i])+'  '+str(Yield1[i])+'\n')
f.close()






