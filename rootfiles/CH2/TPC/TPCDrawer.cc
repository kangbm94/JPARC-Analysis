
Double_t HypTPCdEdx(Double_t Z, Double_t *x, Double_t *p){
  //x : poq
  //p[0] : converting constant p[1] : density effect correction p[2] : mass
  Double_t me  = 0.5109989461;
  Double_t rho = TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672); //[g cm-3]
  Double_t K = 0.307075; //[MeV cm2 mol-1]
  Double_t ZoverA = 17.2/37.6; //[mol g-1]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]
  Double_t I2 = 0.9*188.0 + 0.1*41.7; I2 = I2*I2; //Mean excitaion energy [eV]
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t mom = 1000.*x[0]*Z; //MeV
  Double_t beta2 = mom*mom/(mom*mom+p[2]*p[2]);
  Double_t gamma2 = 1./(1.-beta2);
  Double_t Wmax = 2*me*beta2*gamma2/((me+p[2])*(me+p[2])+2*me*p[2]*(TMath::Sqrt(gamma2)-1));
  Double_t dedx = p[0]*constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2)-beta2-p[1]);
  return dedx;
}

Double_t HypTPCBethe(Double_t *x, Double_t *p){ return HypTPCdEdx(1, x, p); }
Int_t HypTPCdEdxPID_temp(Double_t dedx, Double_t poq){
  Double_t bethe_par[2] = {7195.92, -10.5616};
  Double_t limit = 0.6; //GeV/c
  Double_t mpi = 139.57039;
  Double_t mk  = 493.677;
  Double_t mp  = 938.2720813;
  Double_t md  = 1875.612762;
  TF1 *f_pim = new TF1("f_pim", HypTPCBethe, -3., 0., 3);
  TF1 *f_km = new TF1("f_km", HypTPCBethe, -3., 0., 3);
  TF1 *f_pip = new TF1("f_pip", HypTPCBethe, 0., 3., 3);
  TF1 *f_kp = new TF1("f_kp", HypTPCBethe, 0., 3., 3);
  TF1 *f_p = new TF1("f_p", HypTPCBethe, 0., 3., 3);
  TF1 *f_d = new TF1("f_d", HypTPCBethe, 0., 3., 3);

  f_pim -> SetParameters(bethe_par[0], bethe_par[1], mpi);
  f_km -> SetParameters(bethe_par[0], bethe_par[1], mk);
  f_pip -> SetParameters(bethe_par[0], bethe_par[1], mpi);
  f_kp -> SetParameters(bethe_par[0], bethe_par[1], mk);
  f_p -> SetParameters(bethe_par[0], bethe_par[1], mp);
  f_d -> SetParameters(bethe_par[0], bethe_par[1], md);

  Int_t pid[3] = {0};
  if(poq >= limit){
    pid[0]=1; pid[1]=1; pid[2]=1;
  }
  else if(limit > poq && poq >= 0.){
    Double_t dedx_d = f_d -> Eval(poq); Double_t dedx_p = f_p -> Eval(poq);
    Double_t dedx_kp = f_kp -> Eval(poq); Double_t dedx_pip = f_pip -> Eval(poq);
    if(dedx_d > dedx && dedx >= dedx_kp) pid[2]=1;
    if(dedx_p > dedx){
      pid[0]=1; pid[1]=1;
    }
  }
  else if(0.> poq && poq >= -limit){
    pid[0]=1; pid[1]=1;
  }
  else{
    pid[0]=1; pid[1]=1;
  }

  delete f_pim;
  delete f_km;
  delete f_pip;
  delete f_kp;
  delete f_p;
  delete f_d;

  Int_t output = pid[0] + pid[1]*2 + pid[2]*4;
  return output;
}


void TPCDrawer(){
	Double_t bethe_pars[2] = {7195.92, -10.5616};
	Double_t mprt  = 938.2720813;
		TF1* f_bethe = new TF1("f_betaP",HypTPCBethe,-3,0.1,3);
		f_bethe->SetLineWidth(2);
		f_bethe->SetParameters(bethe_pars[0],bethe_pars[1],mprt);
	TFile* file = new TFile("run05641_DstTPCHelixTracking.root");
	TTree* tree = (TTree*)file->Get("tpc");
	tree->Draw("dEdx:mom0*charge>>(1000,-2,2,500,0,500)","ntTpc>0","colz");
	f_bethe->Draw("same");
}
