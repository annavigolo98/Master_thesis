void macrofit(){
    
    
    //Energy calibration of a detector starting from energy in channels and nominal values for the energies in keV of the peaks of interest. Gaussian functions  are used to fit the peaks in the energy spectrum and get the 
    //centroids in channels. Linear relation between energies in channels and the ones in keV is assumed.

    const int h=12;  //number of data points used
    double x_val[h] = {477.612 ,511.0, 609.321, 1460.820, 2614.511,1173.228, 1332.501, 661.657, 276.3989,302.8508,356.0129,383.8485 };//data x 
    double y_val[h] = {7.04832e+02,7.54264e+02 ,8.99298e+02,2.15188e+03,3.85960e+03,   1.73140e+03,1.96660e+03, 9.75924e+02, 4.07591e+02,4.46619e+02,5.25136e+02,5.66180e+02  };//data y
    
    double x_err[h] = {0.003,0.,0.007,0.005,0.010, 0.003,0.005,0.003, 0.0012,0.0005,0.0007,0.0012 };   
    double y_err[h] = {7.04663e-02,5.01543e-03,1.53571e-01,3.55407e-01,7.16339e-01,   1.61831e-02,1.66704e-02,1.03344e-02, 2.17023e-02,1.12082e-02,5.69437e-03,1.49487e-02 }; //errors  y
    
    
    TCanvas *c1=new TCanvas();
    c1->SetGrid();
    TGraphErrors* g = new TGraphErrors(h, x_val, y_val, x_err, y_err);
    TF1* f = new TF1("f", "[0]+[1]*x",0,3000);
    
    f ->SetLineWidth(1);
    g -> SetMarkerStyle(20);
    g -> SetMarkerColor(2);
    g -> SetLineColor(2);
    f->SetLineColor(38);
   
    
    g->Fit(f);
    //f->SetParameter(0,-2.56610e-01 );
    //f->SetParameter(1,1.47608 );
    
    g->SetTitle("Calibration Ge detector; E (keV) ; E (ch)");
    //g-> GetXaxis()-> SetRangeUser(0,12000);
    //g-> GetYaxis()-> SetRangeUser(-100,1600);
    g->Draw("APE");
    f->Draw("same");
    
    //Info on the fit
    cout << "Chi^2 = " << f->GetChisquare() << endl;
    cout << "NDF = " << f->GetNDF() << endl;
    cout << "R =  " << g->GetCorrelationFactor() << endl;
    double rho= g->GetCorrelationFactor();
    double N=g->GetN();
    double tnc=rho*sqrt( (N-2.)/(1.-rho*rho)  );
    cout<<"tnc  "<<tnc<<endl;
    
    //E_ch=q+m*E_kev
    double m=f->GetParameter(1);
    double q=f->GetParameter(0);
    
    double sigmam=f->GetParError(1);
    double sigmaq=f->GetParError(0);
    
    //E_keV=a+b*E_ch
    double a=-q/m;
    double sigmaa=(1./m)*sqrt( pow(sigmaq,2.) + pow((q*sigmam/m),2.)  );
    double b= 1./m;
    double sigmab= sigmam/pow(m,2.);
    
    cout<<"E_kev=a+bE_ch"<<endl;
    cout<<"a pm sigmaa "<<a<<"  "<<sigmaa<<endl;
    cout<<"b pm sigmab "<<b<<"  "<<sigmab<<endl;
    
  
    
  
    

}

