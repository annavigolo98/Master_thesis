void macrofit(){
    
    
    //Dati
    const int h=12;  //numero dati per le dimensioni degli array
    double x_val[h] = {477.612 ,511.0, 609.321, 1460.820, 2614.511,1173.228, 1332.501, 661.657, 276.3989,302.8508,356.0129,383.8485 };//dati x 
    double y_val[h] = {7.04832e+02,7.54264e+02 ,8.99298e+02,2.15188e+03,3.85960e+03,   1.73140e+03,1.96660e+03, 9.75924e+02, 4.07591e+02,4.46619e+02,5.25136e+02,5.66180e+02  };//dati y
    
    double x_err[h] = {0.003,0.,0.007,0.005,0.010, 0.003,0.005,0.003, 0.0012,0.0005,0.0007,0.0012 };   //0,0,0...se nulli //errori sulle x
    double y_err[h] = {7.04663e-02,5.01543e-03,1.53571e-01,3.55407e-01,7.16339e-01,   1.61831e-02,1.66704e-02,1.03344e-02, 2.17023e-02,1.12082e-02,5.69437e-03,1.49487e-02 }; //errori sulle y
    
    //Canvas 1: grafico con errori + fit
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
    
    //Info aggiuntive al fit
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
    
    //Residui: dati
    //double res_x[h] = {686.75,1570.8,1707.1,1782.8};  //dati x
    double res_y[h];
    
    //Se ci sono sia errori sulle x che sulle y, l'errore sui residui dovrebbe essere quello sulle y una volta che anche l'errore x è stato riportato sulle y
    double res_err[h];
    
    double yth_err[h];
    
    
    for(int i = 0; i < h; i++){
        res_y[i] = y_val[i] - f->Eval(x_val[i]); 
        res_err[i] =sqrt( pow(y_err[i],2.) + pow(sigmam*x_val[i],2.) + pow(m*x_err[i],2.) + pow(sigmaq,2.) );
        
    }
    
    //Canvas 2: grafico dei residui
    TCanvas *c2=new TCanvas();
    TGraphErrors* g_res = new TGraphErrors(h, x_val, res_y, 0, res_err);
    TF1* f_m = new TF1("f_m", "[0]",0.,3000.);
    c2->SetGrid();
    g_res -> SetMarkerStyle(20);
    g_res -> SetMarkerColor(2);
    g_res -> SetLineColor(2);
    f_m->SetLineColor(38);
    
    //g_res -> Fit(f_m);
    
    g_res->SetTitle("Residuals; E (keV) ; residuals (ch)");
    
    //g_res-> GetYaxis()-> SetRangeUser(-0.5,0.5);
    
    g_res->Draw("APE");
    f_m->Draw("same");
}

