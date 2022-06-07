#include "ToFManager.hh"
ToFManager t;
double FitRange1=-2.,FitRange2=2.;
double CalculatedToF(double p, double path, double mass){
	return path*sqrt(mass*mass+p*p)/p/LightSpeed;
}

void CallibrateBtofOffset(){
	double t1 = -0.4,t2 = 0.4; int nb = 80;
	vector<int>ID = {2,0,0,1,2};
	vector<double>Param={0,1};
	double oldOffset[8];
	t.LoadOldCableOffset("./param/CableOffset_before.txt",oldOffset);
	t.MakeParameterFile("./param/CableOffset_after.txt");

	for(int i=0;i<8;++i){
		cout<<oldOffset[i]<<endl;
	}
	TCanvas* c1 = new TCanvas("c1","c1",1500,900);
	c1->Divide(4,2);
	TH1* hist[8];
	TChain* chain = (TChain*)t.MakePublicChain();
	for(int i=0;i<8;++i){
		ID[2]=i;
		int seg = i+1;
		c1->cd(i+1);
		TString HistName = Form("hist%d",seg);
		hist[i]=new TH1D(HistName,HistName,nb,t1,t2);
		TString arg = Form("btof>>")+HistName;
		TCut SegCut=Form("Bh2Seg==%d",seg);
		TCut NhitCut=Form("nhBh2==1");
		chain->Draw(arg,SegCut);
		hist[i]->Fit("f_gaus");

		double p1 = -f_gaus->GetParameter(1)+oldOffset[i];
		Param[0]=p1;
		t.WriteParameter(ID,Param);
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
	}
}


void CallibrateParticleFromHist(int particle){
	double t1 = -2.,t2 = 2.; int nb = 100;

	vector<int>ID = {7,0,0,1,2};
	vector<double>Param={0,-1};

	TF1* f_gaus = new TF1("f_gaus","gaus",-10,10);
	TH1* h_pi[24];
	double Peak[24];
	double p1;
	double oldVTOF[tofnseg];
	t.LoadOldVTOF("./param/VTOF_before.txt",oldVTOF);
	TCut BtofCut = Form("btof<1&&btof>-1");
	TLine* PeakLine[24];
	TLine* WidthLine[24];
	TCanvas* c1 = new TCanvas("c1","c1",1500,900);
	c1->Divide(6,4);
	t.MakeParameterFile("./param/VTOF_after.txt");
	for(int i=0;i<24;i++){
		int seg=i+1;
		ID[2]=i;
		c1->cd(i+1);
		TString HistName = Form("hist%d",i+1);
		h_pi[i]=(TH1*)t.GetStofHisto(i+1,particle);
		h_pi[i]->SetAxisRange(t1,t2);
		h_pi[i]->Draw();
		cout<<"Segment: "<< i+1<<endl;
		double PeakHeight = h_pi[i]->GetMaximum();
		double PeakPosition=h_pi[i]->GetBinCenter(h_pi[i]->GetMaximumBin());
		if(seg>6&&seg<11){
		//	PeakPosition=0.6;
		}
		f_gaus->SetRange(PeakPosition-0.2,PeakPosition+0.3);
		f_gaus->SetParLimits(1,PeakPosition-0.2,PeakPosition+0.3);
		f_gaus->SetParLimits(2,0.05,1);
		h_pi[i]->Fit("f_gaus","R");
		Peak[i]= f_gaus->GetParameter(1);
		double sigma = f_gaus->GetParameter(2);
		PeakLine[i] = new TLine(Peak[i],0,Peak[i],PeakHeight);
		PeakLine[i]	->SetLineWidth(2); 
		PeakLine[i]	->SetLineColor(kBlue); 
		PeakLine[i]	->Draw("Same"); 
		WidthLine[i] = new TLine(Peak[i]-sigma/2,PeakHeight/2,Peak[i]+sigma/2,PeakHeight/2);
		WidthLine[i]	->SetLineWidth(2); 
		WidthLine[i]	->SetLineColor(kBlue); 
		WidthLine[i]	->Draw("Same"); 
		p1 = -f_gaus->GetParameter(1)+oldVTOF[i];
		if(seg==1){
			Param[0]=oldVTOF[i];
		}
		else{
			Param[0]=p1;
		}
		t.WriteParameter(ID,Param);
		gSystem->ProcessEvents();
	}
}
void CallibrateParticleFromChain(int particle){
	double t1 = -2.,t2 = 2.; int nb = 200;

	vector<int>ID = {7,0,0,1,2};
	vector<double>Param={0,-1};
	TString Data = Form("cstof[0]-tTofCalc[%d]",particle-1);

	TF1* f_gaus = new TF1("f_gaus","gaus",-10,10);
	TH1* h_pi[24];
	double Peak[24];
	double p1;
	double oldVTOF[tofnseg];
	t.LoadOldVTOF("./param/VTOF_before.txt",oldVTOF);
	TCut BtofCut = Form("abs(btof)<1");
	TCut MomCut = Form("pKurama<1");
//	TCut MomCut = Form("");
//	TCut TrigCut = Form("trigflag[20]>0");
//	TCut TrigCut = Form("m2Combi==%d",particle-1);
	TCut TrigCut = "";
	TCut BH2Cut = Form("nhBH2==1");
	TCut KuramaTrackCut = Form("ntKurama==1");
//	TCut Chi2Cut = Form("chisqrKurama<50");
	TCut Chi2Cut = Form("");
	TCut PosCut = Form("abs(xtgt)<25&&abs(ytgt)<20");
//	TCut BH2Cut = Form("");
	double M2_min=0,M2_max=0;
	if(particle==1){
		M2_min=0;M2_max=0.15;
	}else if(particle==2){
		M2_min=0.15;M2_max=0.4;
	}else if(particle==3){
		M2_min=0.5;M2_max=1.2;
	}
	TCut M2Cut = Form("abs(m2)>%f&&abs(m2)<%f",M2_min,M2_max);
	TCut GlobalCut=BtofCut&&MomCut&&TrigCut&&M2Cut&&Chi2Cut&&KuramaTrackCut;
	TLine* PeakLine[24];
	TLine* WidthLine[24];
	TCanvas* c1 = new TCanvas("c1","c1",1500,900);
	c1->Divide(6,4);
	t.MakeParameterFile("./param/VTOF_after.txt");
	TChain* chain = (TChain*)t.MakePublicChain();
	for(int i=0;i<24;i++){
		ID[2]=i;
		int seg = i+1;
		c1->cd(i+1);
		TString HistName = Form("hist%d",seg);
		h_pi[i]=new TH1D(HistName,HistName,nb,t1,t2);
		TString arg = Data+">>"+HistName;
		TCut SegCut=Form("TofSeg==%d",seg);
		chain->Draw(arg,GlobalCut&&SegCut);
		cout<<h_pi[i]->GetEntries()<<endl;
		h_pi[i]->Draw();
		cout<<"Segment: "<< i+1<<endl;
		double PeakHeight = h_pi[i]->GetMaximum();
		double PeakPosition=h_pi[i]->GetBinCenter(h_pi[i]->GetMaximumBin());
		if(abs(PeakPosition)>0.7){
			PeakPosition = 0;
		}
		f_gaus->SetRange(PeakPosition-0.2,PeakPosition+0.2);
		f_gaus->SetParLimits(1,PeakPosition-0.2,PeakPosition+0.2);
		f_gaus->SetParLimits(2,0.05,1.5);
		h_pi[i]->Fit("f_gaus","R");
		Peak[i]= f_gaus->GetParameter(1);
		double sigma = f_gaus->GetParameter(2);
		PeakLine[i] = new TLine(Peak[i],0,Peak[i],PeakHeight);
		PeakLine[i]	->SetLineWidth(2); 
		PeakLine[i]	->SetLineColor(kBlue); 
		PeakLine[i]	->Draw("Same"); 
		WidthLine[i] = new TLine(Peak[i]-sigma/2,PeakHeight/2,Peak[i]+sigma/2,PeakHeight/2);
		WidthLine[i]	->SetLineWidth(2); 
		WidthLine[i]	->SetLineColor(kBlue); 
		WidthLine[i]	->Draw("Same"); 
		p1 = -f_gaus->GetParameter(1)+oldVTOF[i];
		Param[0]=p1;
		if(seg<10){
	//		Param[0]=oldVTOF[i];
		}
		t.WriteParameter(ID,Param);
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
	}
}

void CallibrateFromScattering(){
	double t1 = 1.,t2 = 10.; int nb = 200;
	int particle = 1;
	vector<int>ID = {7,0,0,1,2};
	vector<double>Param={0,-1};
	const double BeamMomentum = 1820;
	TString Data = Form("cstof[%d]-path[0]*sqrt(%f*%f+%f*%f)/%f/%f",particle-1,PionMass,PionMass,BeamMomentum,BeamMomentum,BeamMomentum,LightSpeed);

	TF1* f_gaus = new TF1("f_gaus","gaus",-10,10);
	TH1* h_pi[24];
	double Peak[24];
	double p1;
	double oldVTOF[tofnseg];
	t.LoadOldVTOF("./param/VTOF_before.txt",oldVTOF);
	TCut BtofCut = Form("btof<2&&btof>-2");
	TCut TrigCut = Form("");
	TCut BH2Cut = Form("");
	double M2_min=0,M2_max=0;
	TCut GlobalCut=BtofCut&&TrigCut;
	TLine* PeakLine[24];
	TLine* WidthLine[24];
	TCanvas* c1 = new TCanvas("c1","c1",1500,900);
	c1->Divide(6,4);
	t.MakeParameterFile("./param/VTOF_after.txt");
	t.WriteComment("#########################################################");
	t.WriteComment("# VTOF  Callibrated with Scattering Run                 #");
	t.WriteComment("#########################################################");
	TChain* chain = (TChain*)t.MakePublicChain();
	for(int i=0;i<24;i++){
		ID[2]=i;
		int seg = i+1;
		c1->cd(i+1);
		TString HistName = Form("hist%d",seg);
		h_pi[i]=new TH1D(HistName,HistName,nb,t1,t2);
		TString arg = Data+">>"+HistName;
		if(seg==12){
			TString Datas = Form("cstof[%d]-path[0]*sqrt(%f*%f+%f*%f)/%f/%f",1,PionMass,PionMass,BeamMomentum,BeamMomentum,BeamMomentum,LightSpeed);
			arg = Datas+">>"+HistName;
		}
		TCut SegCut=Form("TofSeg==%d",seg);
		chain->Draw(arg,GlobalCut&&SegCut);
		cout<<h_pi[i]->GetEntries()<<endl;
		h_pi[i]->Draw();
		cout<<"Segment: "<< i+1<<endl;
		double PeakHeight = h_pi[i]->GetMaximum();
		double PeakPosition=h_pi[i]->GetBinCenter(h_pi[i]->GetMaximumBin());
		f_gaus->SetRange(PeakPosition-0.4,PeakPosition+0.4);
		f_gaus->SetParLimits(1,PeakPosition-0.3,PeakPosition+0.3);
		f_gaus->SetParLimits(2,0.05,1.5);
		if(seg<3){
			PeakPosition = 8;
			f_gaus->SetRange(PeakPosition-3,PeakPosition+3);
			f_gaus->SetParLimits(1,PeakPosition-1.5,PeakPosition+1.5);
			f_gaus->SetParLimits(2,0.05,2);
		}
		if(seg==12){
			f_gaus->SetRange(PeakPosition-0.4,PeakPosition+0.4);
			f_gaus->SetParLimits(1,PeakPosition-0.3,PeakPosition+0.3);
			f_gaus->SetParLimits(2,0.05,1.5);
		}
		h_pi[i]->Fit("f_gaus","R");
		Peak[i]= f_gaus->GetParameter(1);
		double sigma = f_gaus->GetParameter(2);
		PeakLine[i] = new TLine(Peak[i],0,Peak[i],PeakHeight);
		PeakLine[i]	->SetLineWidth(2); 
		PeakLine[i]	->SetLineColor(kBlue); 
		PeakLine[i]	->Draw("Same"); 
		WidthLine[i] = new TLine(Peak[i]-sigma/2,PeakHeight/2,Peak[i]+sigma/2,PeakHeight/2);
		WidthLine[i]	->SetLineWidth(2); 
		WidthLine[i]	->SetLineColor(kBlue); 
		WidthLine[i]	->Draw("Same"); 
		p1 = -f_gaus->GetParameter(1)+oldVTOF[i];
		Param[0]=p1;
		t.WriteParameter(ID,Param);
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
	}
}

void ToFManager::FitTimewalk(int seg, int UD, int particle,bool mod){
	double ecut = 0.2,emax=1.8;
	if(seg==7&&UD==0) emax=1.5;
	if(seg==8&&UD==1) emax=1.7;
	if(seg==9&&UD==0) emax=1.5;
	if(seg==9&&UD==1) emax=1.5;
	TString ct = Form("ToF%d_%d",seg,UD);
	slew_func->SetRange(ecut,emax);
	double p0min = 0,p0max=5,p1min=-5,t0=-3,t1=3;
	slew_func->SetParLimits(0,-5,1);
	slew_func->SetParLimits(1,p1min,ecut);
	slew_func->SetParLimits(2,t0,t1);
	TH2D* bfh = GetPHCHisto(seg,UD,particle); 
	bfh->SetOption("col");
	slew_func->SetRange(ecut,emax);
	slew_func->SetParLimits(0,-5,0);
	slew_func->SetParLimits(1,p1min,0);
	slew_func->SetParLimits(2,t0,t1);
	if(mod){
		bfh->Fit("slew_func","QR");
		bfh->Fit("slew_func","R");
		double p0 = slew_func->GetParameter(0);
		double p1 = slew_func->GetParameter(1);
		double p2 = slew_func->GetParameter(2);
	}
	else{
		double epb=0.04;//EnergyPerBin
		int bps=3;
		int bin1 = ecut/epb;
		int bin2 = emax/epb;
		double xp[100];
		double yp[100];
		int cnt=0;
		for(int i=bin1;i<bin2+1;i++){
			xp[cnt]=epb*i;
			TH1D* hy = bfh->ProjectionY("py",i-bps/2,i+bps/2);
			yp[cnt]=hy->GetBinCenter(hy->GetMaximumBin());
			cnt++;
		}
		TGraph* gr = new TGraph(cnt,xp,yp);
		gr->SetTitle(ct);
		gr->Draw("APsame");
		//gr->SetMinimum(-3);
		//gr->SetMaximum(3);
		gr->GetXaxis()->SetLimits(0,4);
		gr->GetYaxis()->SetRangeUser(-3,3);
		bfh->Draw("samecol");
		gr->SetMarkerStyle(2);
		gr->SetMarkerSize(2);
		slew_func->SetLineWidth(4);
		gr->Fit("slew_func","QR");
		gr->Fit("slew_func","R");
		double p0 = slew_func->GetParameter(0);
		double p1 = slew_func->GetParameter(1);
		double p2 = slew_func->GetParameter(2);
	}
}
