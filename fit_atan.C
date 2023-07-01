#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

void grafico(const char *namefile) {

TCanvas *c1= new TCanvas("c1");
TGraphErrors *d= new TGraphErrors(namefile,"%lg %lg %lg ");
//TGraphErrors *d2= new TGraphErrors(namefile2,"%lg %lg %lg %lg");


//TF1  *f1= new TF1("f1", "[0]+x*[1]", 0.,12000);



//d->GetXaxis()->SetTitle("d(cm)");
//d->GetYaxis()->SetTitle("S/N");
//d->GetXaxis()->CenterTitle();
//d->GetYaxis()->CenterTitle();
//d->SetTitle("Signal to noise ratio 511 keV ");

d->SetMarkerStyle(6);
d->SetMarkerSize(1);
d->SetMarkerColor(1);
d->SetLineColor(1);
//d->SetLineColor(kRed);

//d->Draw("AP");

//d2->GetXaxis()->SetTitle("d(cm)");
//d2->GetYaxis()->SetTitle("S/N");
//d2->GetXaxis()->CenterTitle();
//d2->GetYaxis()->CenterTitle();
//d2->SetTitle("Signal to noise ratio 511 keV ");

//d2->SetMarkerStyle(6);
//d2->SetMarkerSize(1);
//d2->SetMarkerColor(1);


//d2->Draw("same");



c1->SetGrid();
   
   
TMultiGraph  *mg  = new TMultiGraph();
mg->Add(d);
//mg->Add(d2);
mg->SetTitle("Yield 'deposition 130' 29/03/23");
mg->GetXaxis()->SetTitle("Energy (keV)");
mg->GetYaxis()->SetTitle("Yield");
mg->Draw("AP");


   //TLegend *legend = new TLegend();
   //legend->SetHeader("Header","C"); 				// option "C" allows to center the header
   //legend->AddEntry(d,"CH0");
  // legend->AddEntry(d2,"CH1");
   //gStyle->SetLegendFillColor(0);

   //legend->Draw();


//fit con arcotangente

// x_bar 0.16297
double sx= 273.;
double dx= 342.;


//inizializzo i parametri di fit salita

//double a_par=1.00961e-01;

//double b_par= 1.96335e+00;

//double c_par=-2.78311e+02;

//double d_par=1.50797e-01;

//plot parametri fit

double a_par= 1.03865e-01;

double b_par= 1.79729e+00;

double c_par=-2.78302e+02 ;

double d_par=-1.99144e-02;



//inizializzo i parametri di fit discesa


//double e_par=1.16070e-01; //par4

//double f_par=3.10781e-01; //par5

//double g_par=-3.32730e+02 ; //par6


//plot parfit

double e_par=1.13099e-01; //par4

double f_par=3.20696e-01; //par5

double g_par=-3.32815e+02 ; //par6



TF1* arctan_func = new TF1("arctan_func","[0]*TMath::ATan( [1]*(x+[2]) )  + [3] -[4]*TMath::ATan( [5]*(x+[6]) )",sx,dx);
	arctan_func-> SetParameter(0,a_par);	
	arctan_func-> SetParameter(1,b_par);
	arctan_func-> SetParameter(2,c_par);
	arctan_func-> SetParameter(3,d_par);
	arctan_func-> SetParameter(4,e_par);
	arctan_func-> SetParameter(5,f_par);
	arctan_func-> SetParameter(6,g_par);
//d->Fit(arctan_func,"R");


arctan_func->Draw("same");

// calcolo residui

const int h=d->GetN();
double resy[h];
double x[h];
double erry[h];

// f->Eval(x_val[i])

for (int i=0; i<d->GetN();i++)
{
   double yval=d->GetY()[i];
   double xval=d->GetX()[i];
   double yvalerr=d->GetEY()[i];
   //cout<<"X, Y, errY values "<<xval<<"   "<<yval<<"  "<<yvalerr<<endl;
   x[i]=xval;
   resy[i]=yval-arctan_func->Eval(xval);
   erry[i]=yvalerr;
   }
   
   
TCanvas *c2=new TCanvas();
c2->SetGrid();
TGraphErrors* res = new TGraphErrors(h,x,resy,0,erry);
res->GetXaxis()->SetTitle("Energy (keV)");
res->GetYaxis()->SetTitle("Yield");
res->SetTitle("Residuals 'deposition 130' 29/03/23");

res -> SetMarkerStyle(20);
res -> SetMarkerColor(1);
res -> SetLineColor(1);

res->GetXaxis()->SetRangeUser(273.,342.);
//res->GetYaxis()->SetRangeUser(-0.01, 0.01);
res->Draw("APE");
TF1* retta_res = new TF1("retta_res","[0]",273.,342.);
retta_res->Draw("same");   


} 
