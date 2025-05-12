# Master's thesis programs

Some scripts I used for the target analysis for the 14N(p,g)15O reaction, studied in my master thesis (related to the LUNA experiment to be performed in 2022-2023 at Gran Sasso National Laboratories).  https://www.lngs.infn.it/en/luna.


1) Simulation to calculate the yield profile (integral of the cross section over the effective stoppping power between two energy values) starting from the target composition: number of layers,
thickness in atoms/cm^2, percentages of nitrogen and tantalum in each layer. The cross section is calculated apart from a normalization constant (to be fitted). NAMEFILE: 'integrale_yieldsim_v3.py'.

2) Program to add beam straggling effects by convolution of the yield values with a gaussian function. NAMEFILE: 'convoluzione2.py'.


![provasim27-1](https://github.com/user-attachments/assets/9a682f38-50b2-454b-a0dd-9afc0e0730e4)


5) Program to fit the experimental data with theoretical formulas to obtain the peak and total efficiency curve of a Ge detector. Ratios between the net counts under the peaks of the energy spectra and the expected number of counts for the Co, Cs and Ba radioactive sources is used. Ratio between the primary and secondary transition yields of the 15O nucleus are also fitted. Total efficiency of the Cs source is used. NAMEFILE: 'iminuit1.py'.
   

4) Program to calculate the yield and the expected counting rate starting from the S-factor (to calculate the cross section) and the stopping powers of the target elements. For every proton energy
   the cross section divided by the effective stopping power is integrated over the energy lost by the beam in the target to obtain the yield. Then the rate is calculated as R=yield x I/e x      efficiency, assuming a beam current I=100 muA (e is the electric charge) and target thickness of 10^18 atoms/cm^2, the efficiency is of the Ge detector. The transition of the 15O nucleus from direct capture to the excited level at 6.79 MeV is considered. NAMEFILE: 'sfactorRC6.79.py'.




   

3) Program to analyze '.root' files containing energy spectra for gamma rays acquired by a Ge detector. Gets the charge and timestamps form related histograms saved in the '.root' file
   and calculates the areas under the peaks of interest in the energy spectra (also subtracting the background). Calculates the yield (area under a peak divided by the charge) for each
   beam energy and for each gamma ray transition studied for the 15O nucleus. NAMEFILE 'open_tree_multifilev4.C'.
   
![yieldenergyfit](https://github.com/user-attachments/assets/aac7b868-1507-44b5-b7ea-00dee3367001)

6) fit of the yield profile of the targets analyzed with a double arctangent function: a rising edge defines the first set of parameters; the second one defines the second set. NAMEFILE: 'fit_atan.C'.
   
7) Energy calibration of a detector starting from energy in channels and nominal values for the energies in keV of the peaks of interest. Gaussian functions  are used to fit the peaks in the energy spectrum and get the centroids in channels. Linear relation between energies in channels and the ones in keV is assumed. NAMEFILE: 'calibration.C'.   
