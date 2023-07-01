#include <iostream>
#include <cmath>
#include <vector>


struct slimport_data_t {
	
	
	UShort_t        energy;
	UShort_t        channel;
	UShort_t        cfd;
	ULong64_t	timeStamp; //time stamp
	
	
};

// analyzes more files in '.root' format containing energy spectra, charge, timestamps. Calculates the areas under the peaks of interest in the energy spectrum  subtracting the background. 
//Gets the timestamps and the highest value (the total time of acquisition). Calculates the yield as the ratio of the area and the charge. 
//Yields are calculated for each beam energy and gamma ray transition of the 15O nucleus.  


void readchannel()
{
 
 
 
vector<const char*> dati_input{"run616_lf.root","run617_lf.root","run618_lf.root","run619_lf.root","run620_lf.root","run621_lf.root","run622_lf.root","run623_lf.root","run624_lf.root","run625_lf.root","run626_lf.root","run627_lf.root","run628_lf.root","run629_lf.root","run630_lf.root","run631_lf.root","run632_lf.root","run633_lf.root","run634_lf.root"};

vector <double> charge={15441,16676,12007,12881,18769,14855,14068,13349,12054,11498,10488,10534,10569,10648,12410,25149,16395,19951,20535};

string output1="timestamp_scan.txt";

 ofstream outputfile1(output1);
   if(!outputfile1)
    {   cout<<"errore apertura file!!!"<<endl;
        return (-1);                                              }


 
string output2="scan_general.txt";

 ofstream outputfile2(output2);
   if(!outputfile2)
    {   cout<<"errore apertura file!!!"<<endl;
        return (-1);                                              }
        
        
        
string output4="aree_run.txt";

 ofstream outputfile4(output4);
   if(!outputfile4)
    {   cout<<"errore apertura file!!!"<<endl;
        return (-1);                                              }        
 
 
string output3="yield_scan.txt";

 ofstream outputfile3(output3);
   if(!outputfile3)
    {   cout<<"errore apertura file!!!"<<endl;
        return (-1);                                              } 
        
        
        
        
string output5="charge_scan.txt";

 ofstream outputfile5(output5);
   if(!outputfile5)
    {   cout<<"errore apertura file!!!"<<endl;
        return (-1);                                              }        
        
        
//outputfile        
       


for(int a=0; a<dati_input.size(); a++){  

 //timestamp
 vector <double> timestamp_ch={};
  
 slimport_data_t indata;
 TFile *infile = new TFile(dati_input[a]);
 TTree *intree = (TTree*)infile->Get(Form("Board 0"));
 
 TH1D *h_spectrum= new TH1D("h_spectrum","decayspectrum",32768,0.,32768); //stesso binning della scheda, deve venire uguale alle aree con histo
 
 TH1D *hcharge= new TH1D("hcharge","hcharge",32768,0.,32768); //stesso binning della scheda, deve venire uguale alle aree con histo
 
 TH1D *timestamp= new TH1D("timestamp","timestamp",5000,0.,80.*1e9);
 
 

 
 TBranch *inbranch1 = intree->GetBranch(Form("channel"));
 inbranch1->SetAddress(&indata.channel);
 
 intree->SetBranchAddress("energy",&indata.energy);
 
 for(int i=0;i<intree->GetEntries();i++)
 { 
   intree->GetEntry(i);
   inbranch1->GetEntry(i);
   
   if(indata.channel==2)   //germanio
   {h_spectrum->Fill(indata.energy);
   
    
    };
   
  };
  

 
 
  intree->SetBranchAddress("energy",&indata.energy);
 
 for(int i=0;i<intree->GetEntries();i++)
 { 
   intree->GetEntry(i);
   inbranch1->GetEntry(i);
   
   if(indata.channel==7)   //charge
   {hcharge->Fill(indata.energy);
    
     };
   
  };
  
  
 
  
 intree->SetBranchAddress("timeStamp",&indata.timeStamp);
 
 for(int i=0;i<intree->GetEntries();i++)
 { 
   intree->GetEntry(i);
   inbranch1->GetEntry(i);
   
   if(indata.channel==2)  //timestamp ch2 germanio
   {timestamp->Fill(indata.timeStamp);
    timestamp_ch.push_back(indata.timeStamp);
     };
   
  };
  

  
cout<<"Run  "<<a<<""<<endl;  
cout<<"last timestamp of ch_opened in file number "<<a<<"  "<<timestamp_ch[timestamp_ch.size()-1]<<endl;  
outputfile1<<timestamp_ch[timestamp_ch.size()-1]*1e-8<<endl;

 
 
 
//exterma of the peak(s)  to integrate in the energy spectra 



vector <double> min_peak={1120.,2030.,3495.,11140.,  10000.,9090.,7630.}; //channels
vector <double> max_peak={1140.,2060,3520.,11180.,  10050.,9130.,7665.};

vector <double> m1_bkg={50,50.,50.,50.,  50.,50.,1.};   //number of histogram bins
vector <double> m2_bkg={50.,50.,50.,50.,  50.,50.,50.};


//area under the peaks of interest and yield (area divided by the charge)
vector <double> Area={0.,0.,0.,0.,  0.,0.,0.};
vector <double> Area_sigma={0.,0.,0.,0.,  0.,0.,0.};

vector <double> yield={0.,0.,0.,0.,  0.,0.,0.};
vector <double> sigmayield={0.,0.,0.,0.,  0.,0.,0.};


//background subtraction
for (int j=0;j<min_peak.size();j++)
{
double L=min_peak[j];
double U=max_peak[j];

double m1=m1_bkg[j]; //bin background sx
double m2=m2_bkg[j]; //bin background dx


double n=h_spectrum->FindFixBin(U)-h_spectrum->FindFixBin(L);


double x_in= h_spectrum->FindFixBin(L);
double x_fin= h_spectrum->FindFixBin(U);
//cout<<"L, U "<<L<<"  "<<U<<endl;
//cout<<"extrema  "<<x_in<<"  "<<x_fin<<endl;

double G= h_spectrum->Integral(x_in,x_fin);

x_in=h_spectrum->FindFixBin(L-m1);
x_fin=h_spectrum->FindFixBin(L-1);
//cout<<"L-m1, L-1 "<<L-m1<<"  "<<L-1<<endl;
//cout<<"extrema  "<<x_in<<"  "<<x_fin<<endl;

double bkg1= h_spectrum->Integral(x_in,x_fin);


x_in=h_spectrum->FindFixBin(U+1);
x_fin=h_spectrum->FindFixBin(U+m2);
//cout<<"U+1, U+m2 "<<U+1<<"  "<<U+m2<<endl;
//cout<<"extrema  "<<x_in<<"  "<<x_fin<<endl;

double bkg2= h_spectrum->Integral(x_in,x_fin);

double B=(bkg1/m1+bkg2/m2)*n/2.;
double sigmaB=sqrt( bkg1/(m1*m1) + bkg2/(m2*m2) )*n/2.;


Area[j]=G-B;
Area_sigma[j]=sqrt(G + sigmaB*sigmaB  );


//cout<<"Net counts run  "<<j<<endl;
//cout<<Area[j]<<"   "<<Area_sigma[j]<<endl;

outputfile2<<a<<"   Area (counts)    "<<j<<"   "<<Area[j]<<"   "<<Area_sigma[j]<<endl;





}; //end for j




double carica=charge[a];
//double carica_sigma=0.03*carica;  // 3% statistical error on charge reading
double carica_sigma=0.0;  //no error on the charge

//histogram charge to be compared with counter one

double carica_histo=hcharge->Integral();
//double carica_histo_sigma=0.03*carica_histo;  
double carica_histo_sigma=0.;  

 




for (int j=0.;j<min_peak.size();j++)
{
 yield[j]=  Area[j]/carica;

 sigmayield[j]=  (1./carica)*sqrt( pow(Area_sigma[j],2.) + pow( (carica_sigma*Area[j]/carica) ,2.)      );
 
};



outputfile2<<a<<"   charge (mu Coulomb)    "<<carica<<"   "<<carica_sigma<<endl;
outputfile2<<a<<"   last timestamp (10ns)    "<<timestamp_ch[timestamp_ch.size()-1]<<endl;
outputfile2<<"                   "<<endl;
outputfile2<<"                   "<<endl;

outputfile5<<carica<<"   "<<carica_sigma<<"    "<<carica_histo<<"  "<<carica_histo_sigma<<endl;


for (int i=0;i<Area.size();i++)
{  
  outputfile4<<Area[i]<<"   "<<Area_sigma[i]<<"    ";
  outputfile3<<yield[i]<<"   "<<sigmayield[i]<<"    ";
  
};

outputfile4<<"   "<<endl;
outputfile3<<"   "<<endl;





}; //end for a


};
