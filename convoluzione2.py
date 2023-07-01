import numpy as np
import math as m
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize as opt
import scipy
import ROOT


from iminuit import Minuit,cost
import uncertainties

graphStopPowerNitrogen = ROOT.TGraph( "e_tot_stopping_N.txt" )
graphStopPowerTantalum = ROOT.TGraph( "e_tot_stopping_Ta.txt" )
graphStopPowerArgon = ROOT.TGraph( "e_tot_stopping_Ar.txt" )
graphStopPower15N = ROOT.TGraph( "e_tot_stopping_15N.txt" )


graphYield=ROOT.TGraph("Yield_TaN_dec_1_simulation.txt")

#COMPOSIZIONE LAYERS dec_1
#layer_dx=[0.547524, 0.140000, 0.100000, 0.100000 ] #10^18at/cm^2
#percN=[55.0132, 41.2155 , 12.0000 ,  3.1721  ]
#percTa=[100.-55.0132,100.-41.2155,100.-12.0000,100.-3.1721]


#COMPOSIZIONE LAYERS dec_1 prova 1
#layer_dx=[0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ] #10^18at/cm^2
#percN=[42.6,50.1,51.9,52.2,59.7,48.3,41.4,24.4,9.7,3.7 ]
#percTa=[]


#COMPOSIZIONE LAYERS dec_1 prova
layer_dx=[0.887524/30.,0.887524*2./30.,0.887524/10.,0.887524/10.,0.887524/10. ,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10.,0.887524/10. ] #10^18at/cm^2
percN=[31.1,47.6,50.1,51.9,52.2,59.7,48.3,41.4,24.4,9.7,3.7 ]
percTa=[]

for i in range(len(percN)):
   percTa.append(100.-percN[i])


  
#COMPOSIZIONE LAYERS dec_2
   
#layer_dx=[0.100234,0.134804,0.103019,0.229821,0.115302,0.104322,0.111275] #10^18at/cm^2
#percN=[72.6546,53.7575, 53.6916, 47.4629,34.8149, 22.8344, 11.9323    ]
#percTa=[100.-72.6546,100.-53.7575,100.-53.6916,100.-47.4629,100.-34.8149,100.-22.8344, 100.-11.9323]
  
#COMPOSIZIONE LAYERS mar_3
#layer_dx=[0.10001,0.27856,0.102153,0.375279,0.202666 ] #10^18at/cm^2
#percN=[29.2815, 42.9672, 32.3414, 20.0711, 6.8396   ]
#percTa=[100.-29.2815,100.-42.9672,100.-32.3414,100.-20.0711,100.-6.8396]























def stopping_eff(ind,E_x): #ind=indice di layer da usare, E_x energia per la quale voglio calcolare stopping power e effective stopping power.

  
   e_totN14 = graphStopPowerNitrogen.Eval( E_x ) #stopping power totale a E_x per N14
   e_totTa  = graphStopPowerTantalum.Eval( E_x ) #stopping power totale a E_x per Ta
   
   
   e_eff_TaN= e_totN14 + (percTa[ind]/percN[ind])* e_totTa # effective stopping power in base al layer

   return(e_eff_TaN)










#calcolo dx integrato da E_p-DE a E_p per un singolo elemento es Ta o N
def integral_dx(E_1,E_2,ind,d):  

   # Calcoliamo gli estremi dell'integrale
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
     
     
   for i in range( nSteps ):
       if(d=='N'):
         stopPower     = stopping_N(E_step,ind) #stopping power in E_step
       if(d=='Ta'):
         stopPower     = stopping_Ta(E_step,ind) #stopping power in E_step
         
       integral     += step*1./stopPower
       E_step        += step

   return(integral)


#gaussiana nel punto x
def gaussian( x, s ):
    
    #mu=(258.+320.)/2.
    return np.exp( -(x)**2 / ( s*s*2 ) )


def gauss( s,len_xfitte ): 
    
   
    N=len_xfitte
    step=0.05 
    DE=N*step/2.
    E_min=-DE
    E_max=DE
    N=len_xfitte
    
    
    dx = np.arange(E_min ,E_max,step)  #step uguale o minore della differenza E_p[i+1]-E_p[i]
    
    #print('len gaussian i',len(dx),'\n')
    #print('len profile',N,'\n')
    #print('N_steps ',N,'\n')
    #print('step',(E_max-E_min)/N,'\n')
    
    
    kernel = gaussian( dx, s )
    return kernel / sum( kernel )

#stopping power nel punto E_x di N e Ta
#interpolazione lineare

def stopping_N(E_x,ind):
    
    e_tot = graphStopPowerNitrogen.Eval( E_x )*0.01*percN[ind]
    
    return(e_tot)
    
def stopping_Ta(E_x,ind):
    
    e_tot = graphStopPowerTantalum.Eval( E_x )*0.01*percTa[ind]
    
    return(e_tot)
    
#sigma straggling a E_p


def calc_DE(E_x,ind,d):
   DE=0.
   # Calcoliamo gli estremi dell'integrale
   Xmax= layer_dx[ind]  
   Xmin = 0. #E_1=E-DE
   
   # Definiamo il numero di step
   nSteps = 1000

   # Calcoliamo la grandezza dello step
   step = (Xmax - Xmin)/nSteps
   
   E_step=E_x
 
   for i in range( nSteps ):
   
        if (d=='N'):
           stopPower     = stopping_N(E_step,ind) #stopping power in E_step   nel lab
           
        if (d=='Ta'):
           stopPower     = stopping_Ta(E_step,ind) #stopping power in E_step   nel lab
           
              
        E_step        -= step*stopPower
  
   #calcolo DE
   DE=E_x-E_step
   
   return(DE)
   
   
    
def calculateSigmaStraggling(E_x):  # E_x= proton energy
  
  
  
  
  
  
  
  E_res=259.56*15./14.  #(278.1)
  sigma_res=( 1./(2.*np.sqrt(2.*np.log(2.))) )*0.9893*15./14.
  E_res=E_res
  
  
   
  spessore1=0.
  DE_N=0.
  d='N'
  
  E_sup=E_x
  DE_N=calc_DE(E_x,0,d)  #calcola DE nel layer
  E_inf=E_sup-DE_N
  
  for j in range(len(layer_dx)):
        
     if(E_inf<E_res and E_sup>E_res):
          
          spessore1+=integral_dx(E_res,E_sup,j,d) 
          print('break',E_x,'\n')        
          break;  
          
     if(E_inf>=E_res):
          spessore1+=integral_dx(E_inf,E_sup,j,d)   
     
     
     if(j==Nlayer-1):
        
        break 
            
     E_sup=E_inf    
     DE_N=calc_DE(E_sup,j+1,d)
     E_inf=E_sup-DE_N
       
   
  A_N=np.sqrt(14.*(spessore1))*1.2 #E in keV dE/dx in keV/(10^18at/cm^2) dx in 10^18at/cm^2   
  
  
 
  
  spessore2=0.
  DE_Ta=0.
  d='Ta'
  
  E_sup=E_x
  DE_Ta=calc_DE(E_x,0,d)  #calcola DE nel layer
  E_inf=E_sup-DE_Ta
  
  for j in range(len(layer_dx)):
        
     if(E_inf<E_res and E_sup>E_res ):
          spessore2+=integral_dx(E_res,E_sup,j,d) 
          print('break',E_x,'\n')   
          break;  
          
          
     if(E_inf>=E_res):
          spessore2+=integral_dx(E_inf,E_sup,j,d)   
     
     
     if(j==Nlayer-1):
       
        break 
            
     E_sup=E_inf #E-DE   
     DE_Ta=calc_DE(E_sup,j+1,d)
     E_inf=E_sup-DE_Ta
       
  
  
  
  A_Ta=np.sqrt(73.*(spessore2))*1.2 #E in keV dE/dx in keV/(10^18at/cm^2)  dx in 10^18at/cm^2  #stopping power a E_p-DE_Ta/2
    
  
  sigma_N=A_N/(2.*np.sqrt(2.*np.log(2.)) )
  
  
  sigma_Ta=A_Ta/(2.*np.sqrt(2.*np.log(2.)) )
  
  sigma=np.sqrt(sigma_N**2.+sigma_Ta**2.)  
       
  return(sigma)


# convoluzione; input array x più fitte delle energie di fascio , y energie di fascio della stessa dimensione delle energie di fascio    
def straggling( x, y ):
    
    # Definiamo un array della larghezza del profilo
    
    conv = np.zeros( shape=[len(y)] )
    
    sigma=np.zeros(len(y),float)
    
    
    for i in range( len( y ) ):
        temp    = np.zeros( shape=[len( y )] )
        temp[i] = y[i]
        
        
        sigma1   = calculateSigmaStraggling( x[i] )
        if(sigma1!=0.):
          temp    =  np.convolve(y[i], gauss( sigma1,len(x)*10 ), mode='full')[::10]  #array con la convoluzione di dimensioni del profilo
        
             
        print("sigma_straggling ",sigma1,'\n')
        conv   += temp #sommo temp a conv
        sigma[i]=sigma1
      
    return conv, sigma

# prende in input l'array delle energie del fascio e calcola lo straggling su un range più fitto
def straggling2( x_, par  ): 
    
    conv=np.zeros(len(x_),float)
    sigma=np.zeros(len(x_),float)
    
    for i in range(len(x_)):
      
      #par=0.589 #fit
      sigma1   = par*np.sqrt(x_[i]-278.+0.001)
      #sigma1   = calculateSigmaStraggling( x_[i] )
      sigma[i]=sigma1

      if (x_[i]>278. and x_[i]<319.5):   #straggling dalla risonanza in poi  
        E_inf=x_[i]-20.
        E_sup=x_[i]+20. #xmax a 320
       
        INT= integrate.quad(lambda x: (1./(np.sqrt(2.*np.pi)*sigma1))*np.exp( -(x-x_[i])**2 / ( sigma1*sigma1*2 ) )*graphYield.Eval(x),E_inf ,E_sup)
        conv[i]=INT[0]    
      else:
        conv[i]=graphYield.Eval(x_[i])
    return conv,sigma


    
    
    
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


a=14./15.    
Ep=np.arange(275.,300.,0.1) 
Nlayer=len(layer_dx) #NUMERO LAYERS #4;7;5
E_0=275. #energia minima    
    


#fnameYield='Yield_TaN_dec_1_simulation.txt'
#Ep, Yield = np.genfromtxt(fnameYield,dtype='float',comments='#',usecols=(0,1),unpack=True) 

#fnameYield='Yield_TaN_dec_2_simulation.txt'
#Ep, Yield = np.genfromtxt(fnameYield,dtype='float',comments='#',usecols=(0,1),unpack=True) 

    
fnameYield='Yield_TaN_dec_1_simulation.txt'
Ep, Yield = np.genfromtxt(fnameYield,dtype='float',comments='#',usecols=(0,1),unpack=True)    

#sigma della convoluzione e yield convoluta

sigma=np.zeros(len(Yield),float) 

Yield2=np.zeros(len(Ep),float)

#dati sper

fnameexpdata='yieldexp2.txt' #dati 27 mar
Eexp, Yieldexp, Ysigmaexp = np.genfromtxt(fnameexpdata,dtype='float',comments='#',usecols=(0,1,2),unpack=True) 
s=0.3
Yield2, sigma = straggling2( Ep,s)
#residui e chiq

res=np.zeros(len(Yieldexp),float)


def chi2(a):
   
   chitot=0.
   chi2=np.zeros(len(Yieldexp),float)
   #calcolo chi2
 
 
   #fit anna lmfit
   #for i in range(len(Y_exp)):
   #       chi2[i]= pow((Y_exp[i]-graphYieldstragg.Eval(E_exp[i]))/sigmaY_exp[i],2.);
   #       chitot+=chi2[i]
 
   #fit jakub lmfit
   for i in range( len(Yieldexp) ):
        idx = (np.abs(Ep - Eexp[i])).argmin()
        res[i]=Yieldexp[i] - Yield2[idx]*a
        chi2[i] = pow( ( Yieldexp[i] - Yield2[idx]*a )/Ysigmaexp[i], 2 )
        chitot+=chi2[i]
           
   print("Chiquadro tot ",chitot,'\n')
   print("        \n")
   return (chitot);






a=5./(1.44*1e11) #costante per scalare la simulazione ai dati sperimentali


 
m = Minuit(chi2,a)
m.simplex()
m.migrad()  # run optimiser
m.hesse()   # run covariance estimator


a=m.values[0]

print('a  ',a,'\n')
print(m.values)
print(m.errors)

#<ValueView a=3.4706666468774385e-11>
#<ErrorView a=1.5905760743682655e-13>


#plot residui
plt.errorbar(Eexp,res,Ysigmaexp,fmt='None',color='black',marker='o')



plt.xlabel('E_p (keV)')
plt.ylabel('Y_{exp}-Y_{sim} (1e-6 atoms)')
plt.show()

#plot valori
plt.plot(Ep,Yield2*a)
plt.plot(Ep,Yield*a)
plt.errorbar(Eexp,Yieldexp,Ysigmaexp,fmt='None',color='green',marker='s')



plt.xlabel('E_p (keV)')
plt.ylabel('Yield (1e-6 atoms)')
plt.show()

#plot sigma straggling
plt.scatter(Ep,sigma,marker='o',s=1)


plt.xlabel('E_p (keV)')
plt.ylabel('sigma straggling keV')
plt.show()



M_0=1.00727647
M_1=14.0067

E_rcm=259.56 #keV

E_rlab=((M_0+M_1)/M_1)*E_rcm

lambda2_2= ( ((M_0+M_1)/M_1)**2.)* (   4.125e-18/(M_0*(E_rlab*1e3))   ) #massa in u energia in eV


moltfactor=[]
for j in range(len(layer_dx)):
    effstop=stopping_eff(j,E_rlab)
    moltfactor.append( (  (M_1+M_0)/M_1 )*effstop*2./lambda2_2)
    
moltfactor=np.array(moltfactor)
    
print('Moltipl factor x layer j ',moltfactor,'\n')    
    


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
        
    
    
    
    
    
    
