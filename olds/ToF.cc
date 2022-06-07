#include "ToFManager.hh"
ToFManager t;
double pos_offset[24]=
{15.1,14.8,16.7,17.3,17.0,15.3,//
	16.3,15.3,15.1,15.3,15.3,14.7,//
	14.8,15.0,15.5,15.5,16.1,16.8,//
	16.9,16.7,16.9,17.5,16,14.5//
};
void ToF(){
	gStyle->SetOptFit(1111);
//	t.LoadFile("./rootfiles/Dsts/BigKuramaHodoscope.root");
	t.LoadFile("./rootfiles/Dsts/BigKKAna.root");
//	t.LoadKHodo();
	t.LoadKK();
//	t.LoadFile("./rootfiles/Dsts/KH_3.root");
	
//	t.LoadFile("./rootfiles/Dsts/DstKuramaHodoscope_wo_TOFPHC.root");
}
void Monitor(int seg){
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);
	c1->Divide(2,2);
	int ScatHistoNumber = 2200+seg;
	int M2HistoNumber = 2100+seg;
	c1->cd(1);
	t.Get2DHistoFromNumber(ScatHistoNumber);
	c1->cd(2);
	t.Get1DHistoFromNumber(M2HistoNumber);
	c1->cd(3);
	TH1* Hist_Pion = (TH1*)t.GetStofHisto(seg,0);
	Hist_Pion-> Draw();
	c1->cd(4);
	TH1* Hist_Proton = (TH1*)t.GetStofHisto(seg,0);
	Hist_Proton-> Draw();
}
void DrawM2FromKH(){
	TString M2 = Form("sqrt(m2[0])>>(1000,0,2)");
	int KuramaCut = 200;
	TCut CutKurama = Form("chisqrKurama<%d",KuramaCut);	
	int xtgtCut = 25,ytgtCut=20;
	TCut CutTgt = Form("sqrt(xtgtKurama*xtgtKurama)<%d&&sqrt(ytgtKurama*ytgtKurama)<%d",xtgtCut,ytgtCut);
	TCut Cut = CutKurama&&CutTgt;
	t.DrawFromKHodo(M2,Cut);
}
void DrawM2FromKK(int seg){
	vector<double> Origin = {4,0,-67};
	vector<double> Size = {25,20,65};
	double CloseDist = 100;
	
	double KuramaChi = 50;
	double K18Chi = 10;

	double Mom = 1.;

	TString Mass = Form("sqrt(m2[0])>>(1000,0,2)");
	TString M2Scat = Form("pKurama:qKurama*m2[0]>>(400,-1.5,2.5,100,0,2.5)");
	TString MassScat = Form("pKurama:qKurama*sqrt(m2[0])>>(400,-1.5,2.5,100,0,2.5)");

	TCut SegCut;
	if(seg!=0){
	 SegCut = Form("tofsegKurama==%d",seg);
	}
	else{
	 SegCut = Form("tofsegKurama!=%d",seg);
	}
	TCut KuramaCut = Form("chisqrKurama<%f",KuramaChi);	
	TCut K18Cut = Form("chisqrK18<%f",KuramaChi);	
	TCut ChiCut = K18Cut&&KuramaCut;
	TCut VertCut= VertexCut(Origin,Size);
	TCut CloseDistCut = Form("closeDist<%f",CloseDist);
	
	TCut MomCut = Form("pKurama<%f",Mom);

	TCut Cut = SegCut&&VertCut&&CloseDistCut&&ChiCut; 
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);
//	c1->Divide(2,2);
/*
	c1->cd(1);
	t.DrawFromKK(Mass,Cut);
	c1->cd(2);
	t.DrawFromKK(Mass,Cut&&MomCut);
	c1->cd(3);
	*/
	t.DrawFromKK(M2Scat,Cut,"colz");
	c1->SaveAs(Form("./Files/PDF/M2Scat%d.pdf",seg));
	//	t.DrawFromKK(MassScat,Cut,"colz");
}
void DoKK(){
	for(int seg=0;seg<25;++seg){
		DrawM2FromKK( seg);
		gSystem->ProcessEvents();
	}
}
void DoProton(){
	TF1* f_gaus = new TF1("f_gaus","gaus",-10,10);
	TH1* h_pi[24];
	TFile* f = new TFile("ProtonCali.root","recreate");
	double t_range_max = 3;
	double t_pi[24];
	double t_range[24];
	double t_bin[24];
	double p1;
	fstream file;
//	file.open("HodoParam_KBM");
	file.open("STOF_Pion");
	double buf[7];
	double param_offset[24];
	TLine* PeakLine[24];
	for(int i=0;i<24;i++){
		ReadTSV(file,buf);
		param_offset[i]=buf[5];
	}
	TCanvas* c1 = new TCanvas("c1","c1",1500,900);
	c1->Divide(6,4);
	t.MakeParameterFile("STOF_Proton");
	for(int i=0;i<24;i++){
		c1->cd(i+1);
		h_pi[i]=(TH1*)t.GetStofHisto(i+1,2);//0 -> Pion, 1 -> Kaon, 2 -> Proton
		double peak = h_pi[i]->GetBinContent(h_pi[i]->GetMaximumBin());
		double PeakPosition= h_pi[i]->GetBinCenter(h_pi[i]->GetMaximumBin());
		f_gaus->SetRange(PeakPosition-0.3,PeakPosition+0.3);
		f_gaus->SetParLimits(1,PeakPosition-0.2,PeakPosition+0.2);
		f_gaus->SetParLimits(2,0.05,0.3);
		h_pi[i]->Fit("f_gaus","R");
		t_pi[i]= f_gaus->GetParameter(1);
		PeakLine[i] = new TLine(PeakPosition,0,PeakPosition,1.1*peak);
		PeakLine[i]	->SetLineWidth(2); 
		PeakLine[i]	->SetLineColor(kBlue); 
		PeakLine[i]	->Draw("Same"); 
//		h_pi[i]->Draw();
		p1 = f_gaus->GetParameter(1);
		if(i>6) p1 = 0; 
		t.WriteParameter(7,0,i,1,2,param_offset[i]-p1,-1);
		f->cd();
		h_pi[i]->Write();
	PeakLine[i]->Write();
	}
}
void DoPion(){
	TF1* f_gaus = new TF1("f_gaus","gaus",-10,10);
	TH1* h_pi[24];
	TFile* f = new TFile("PionCali.root","recreate");
	double t_range_max = 3;
	double t_pi[24];
	double t_range[24];
	double t_bin[24];
	double p1;
	fstream file;
	file.open("HodoParam_KBM");
	double buf[7];
	double param_offset[24];
	TLine* PeakLine[24];
	for(int i=0;i<24;i++){
		ReadTSV(file,buf);
		param_offset[i]=buf[5];
	}
	TCanvas* c1 = new TCanvas("c1","c1",1500,900);
	c1->Divide(6,4);
	t.MakeParameterFile("STOF");
	for(int i=0;i<24;i++){
		c1->cd(i+1);
		h_pi[i]=(TH1*)t.GetStofHisto(i+1,0);
		double peak = h_pi[i]->GetBinContent(h_pi[i]->GetMaximumBin());
		double PeakPosition= h_pi[i]->GetBinCenter(h_pi[i]->GetMaximumBin());
		f_gaus->SetRange(PeakPosition-0.3,PeakPosition+0.3);
		f_gaus->SetParLimits(1,PeakPosition-0.2,PeakPosition+0.2);
		f_gaus->SetParLimits(2,0.05,0.3);
		h_pi[i]->Fit("f_gaus","R");
		t_pi[i]= f_gaus->GetParameter(1);
		PeakLine[i] = new TLine(PeakPosition,0,PeakPosition,1.1*peak);
		PeakLine[i]	->SetLineWidth(2); 
		PeakLine[i]	->SetLineColor(kBlue); 
		PeakLine[i]	->Draw("Same"); 
//		h_pi[i]->Draw();
		p1 = f_gaus->GetParameter(1);
		t.WriteParameter(7,0,i,1,2,param_offset[i]-p1,-1);
		f->cd();
		h_pi[i]->Write();
	PeakLine[i]->Write();
	}
}
void DoPHC(){
	gStyle->SetOptFit(1111);
	t.MakeParameterFile("ToF_PHC");
	t.SaveHisto("ToFPHC.root");
	int particle = -1;
	for(int seg=1;seg<25;seg++){
		if(seg<7){
			particle = 0;
			t.FitTimewalk(seg,0,particle,0);
		}
		else{
			particle = 2;
			t.FitTimewalk(seg,0,particle,0);
		}
	}
	t.NextParameter();
	for(int seg=1;seg<25;seg++){
		if(seg<7){
			particle = 0;
			t.FitTimewalk(seg,1,particle,0);
		}
		else{
			particle = 2;
			t.FitTimewalk(seg,1,particle,0);
		}
	}
}
void DoPHC(int seg,int UD){
	gStyle->SetOptFit(1111);
	t.SaveHisto(Form("ToFPHC_%d.root",seg));
	t.FitTimewalk(seg,UD,-1,0);
}

void FitADC(){
	gStyle->SetOptFit(1111);
	t.LoadFile("./rootfiles/CH2/Hodoscope05477.root");
	TFile* file = new TFile("ToFADC.root","recreate");
	TH1D* h[2*tofnseg];
	TString ht[2*tofnseg];
	TH1D* h_ped[2*tofnseg];
	TString ht_ped[2*tofnseg];
	double ped[2*tofnseg];double mpv[2*tofnseg];
	TF1* f_landau = new TF1("f_landau","landau",0,2000);
	TF1* f_gaus = new TF1("f_gaus","gaus",0,500);
	t.MakeParameterFile("tof_ADC");
	for(int i=0;i<tofnseg;i++){
		int ds = i+tofnseg;
		cout<<i<<endl;
		file->cd();
		ht_ped[i]=Form("ToFUPED%d",i+1);
		ht_ped[ds]=Form("TOFDPED%d",i+1);
		h_ped[i]=new TH1D(ht_ped[i],ht_ped[i],500,0,500);
		h_ped[ds]=new TH1D(ht_ped[ds],ht_ped[ds],500,0,500);
		t.GetPEDHisto(ht_ped[i],i+1,0);
		t.GetPEDHisto(ht_ped[ds],i+1,1);
		h_ped[i]->Draw();
		h_ped[ds]->Draw();
		double mean = h_ped[i]->GetMean();
		f_gaus->SetRange(mean-50,mean+50);
		t.LoadChain("./rootfiles/CH2/KuramaTracking05453.root");
		h_ped[i]->Fit("f_gaus","QR");
		h_ped[i]->Fit("f_gaus","QR");
		ped[i]=f_gaus->GetParameter(1);
		mean = h_ped[ds]->GetMean();
		f_gaus->SetRange(mean-50,mean+50);
		h_ped[ds]->Fit("f_gaus","QR");
		h_ped[ds]->Fit("f_gaus","QR");
		ped[ds]=f_gaus->GetParameter(1);

		h_ped[i]->Write();
		h_ped[ds]->Write();
	}
	for(int i=0;i<tofnseg;i++){
		int ds = i+tofnseg;
		double sig;
		cout<<i<<endl;
		file->cd();
		ht[i]=Form("ToFUADC%d",i+1);
		ht[ds]=Form("TOFDADC%d",i+1);
		h[i]=new TH1D(ht[i],ht[i],500,0,3000);
		h[ds]=new TH1D(ht[ds],ht[ds],500,0,3000);
		t.GetADCHisto(ht[i],i+1,0);
		t.GetADCHisto(ht[ds],i+1,1);
		h[i]->Draw();
		h[ds]->Draw();
		f_landau->SetRange(500,2500);
		h[i]->Fit("f_landau","QR");
		mpv[i]=f_landau->GetParameter(1);
		sig=f_landau->GetParameter(2);
		f_landau-> SetRange(mpv[i]-sig,mpv[i]+5*sig);
		h[i]->Fit("f_landau","QR");
		mpv[i]=f_landau->GetParameter(1);
		f_landau->SetRange(500,2500);
		h[ds]->Fit("f_landau","QR");
		mpv[ds]=f_landau->GetParameter(1);
		sig=f_landau->GetParameter(2);
		f_landau-> SetRange(mpv[ds]-sig,mpv[ds]+5*sig);
		h[ds]->Fit("f_landau","QR");
		mpv[ds]=f_landau->GetParameter(1);
		h[i]->Write();
		h[ds]->Write();
	}
	for(int i=0;i<tofnseg;i++){
		t.WriteParameter(7,0,i,0,0,ped[i],mpv[i]);
	}
	for(int i=0;i<tofnseg;i++){
		t.WriteParameter(7,0,i,0,1,ped[i+tofnseg],mpv[i+tofnseg]);
	}
	file->Save();
}

void FitTDC(){
	TF1* f_gaus = new TF1("f_gaus","gaus",0,1e7);
	t.LoadFile("./rootfiles/CH2/Hodoscope05395.root");
	TH1* h[2*tofnseg];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	double pp1,pp2;
	double width = 2500;
	t.MakeParameterFile("tof_TDC");
	for(int i=0;i<tofnseg;i++){
		cout<<i<<endl;
		h[i] = (TH1*)t.GetHisto(i+1,0);
		c1->cd(1);
		h[i]->Draw();
		pp1 = h[i]->GetBinCenter(h[i]->GetMaximumBin());
		f_gaus->SetParLimits(1,pp1-width,pp1+width);
		f_gaus->SetRange(pp1-width,pp1+width);
		h[i]->Fit("f_gaus","R");
		h[i]->SetAxisRange(pp1-5*width,pp1+5*width);
		int t0 = f_gaus->GetParameter(1);
		t.WriteParameter(7,0,i,1,0,t0,-0.000939002);
	}
	for(int i=0;i<tofnseg;i++){
		cout<<i<<endl;
		h[tofnseg+i] = (TH1*)t.GetHisto(i+1,1);
		c1->cd(2);
		h[tofnseg+i]->Draw();
		pp1 = h[tofnseg+i]->GetBinCenter(h[tofnseg+i]->GetMaximumBin());
		h[tofnseg+i] = (TH1*)t.GetHisto(i+1,1);
		f_gaus->SetParLimits(1,pp1-width,pp1+width);
		f_gaus->SetRange(pp1-width,pp1+width);
		h[tofnseg+i]->Fit("f_gaus","R");
		h[tofnseg+i]->SetAxisRange(pp1-5*width,pp1+5*width);
		int t0 = f_gaus->GetParameter(1);
		t.WriteParameter(7,0,i,1,1,t0,-0.000939002);
	}
}
