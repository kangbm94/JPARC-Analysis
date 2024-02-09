#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
double Gaus(double* x, double* p){
	double mean = x[0]-p[1];
	double up = mean / p[2];
	up=up*up;
	up/=2;
	return p[0] *exp(-up);
}
double DGaus(double* x,double* p){
	double p1[3];
	double p2[3];
	for(int i=0;i<3;++i){
		p1[i]=p[i];
		p2[i]=p[3+i];
	}
	return Gaus(x,p1)+Gaus(x,p2);
}
TF1* fdgaus = new TF1("fdgaus","DGaus",1.5,1.6,6);
TF1* fxistar = new TF1("xistar","Gaus",1.45,1.65,3);
TF1* fxibg = new TF1("xibg","Gaus",1.5,1.6,3);
void DrawMM(){
	SetStyle();
	TFile* f2 = new TFile("TPCInvM_cd7_WS_firstT.root");	
	TTree* tr2 = (TTree*)f2->Get("tree");
	double mm2;
	tr2->SetBranchAddress("MM",&mm2);
	TH1D* MMHist = new TH1D("Missing Mass","Missing Mass",100,1.2,1.7);
	int ent = tr2->GetEntries();
	for(int i=0;i<ent;++i){
		tr2->GetEntry(i);
		MMHist->Fill(mm2);
	}
	fdgaus->SetParameter(0,100);
	fdgaus->SetParameter(1,1.53);
	fdgaus->SetParLimits(1,1.52,1.555);
	fdgaus->SetParameter(2,0.020);
	fdgaus->SetParLimits(2,0.005,0.020);
	fdgaus->SetParameter(3,150);
	fdgaus->SetParameter(4,1.53);
	fdgaus->SetParLimits(4,1.45,1.6);
	fdgaus->SetParameter(5,0.05);
	fdgaus->SetParLimits(5,0.03,0.1);
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	double tsize = 0.05;
	MMHist->Draw();
	MMHist->Fit("fdgaus","R");
	double pars[3];
	double parb[3];
	for(int i=0;i<3;++i){
		pars[i]=fdgaus->GetParameter(i);
		parb[i]=fdgaus->GetParameter(3+i);
		fxistar->SetParameter(i,fdgaus->GetParameter(i));
		fxibg->SetParameter(i,fdgaus->GetParameter(3+i));
	}
	fxistar->SetLineColor(kBlue);
	fxistar->Draw("same");
	fxibg->SetLineColor(kBlack);
	fxibg->Draw("same");
	MMHist->GetYaxis()->SetTitleSize(tsize);
	MMHist->GetXaxis()->SetTitleSize(tsize);
	MMHist->GetXaxis()->SetTitle("M_{p(K-,K+)X} [GeV/c^{2}]");
	MMHist->GetYaxis()->SetTitle("Entries / 10 MeV/c^{2}");
	c1->SaveAs("MissingMass.png");
	double bin_density = 100. / 0.5;
	double XiCount=pars[0]*pars[2]*sqrt(2*acos(-1))*bin_density;
	double BgCount=parb[0]*parb[2]*sqrt(2*acos(-1))*bin_density;
	cout<<"Entries = "<<XiCount<<"Ratio = "<<XiCount/BgCount<<endl;
	cout<<MMHist->GetEffectiveEntries()<<endl;
}

void DrawM2(){
	SetStyle();
	TFile* f2 = new TFile("dstfiles_ws/run05641_DstKScatAna.root");	
	auto H = (TH2D*)f2->Get("h6004");
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	H->Draw("colz");
	H->GetXaxis()->SetTitle("Mass^{2}/Z [GeV^{2}/c^{4}]");
	H->GetYaxis()->SetTitle("P[GeV / c]");
	c1->SaveAs("KuramaM2.png");
}
