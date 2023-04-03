#include "../../../KKManager.hh"
double PI = acos(-1);
double mm1,mm2,MLd,tcm,MXi,DistLd,MomLd;
bool FlgLd,FlgXi,InTargetLd;
TFile* f1;TFile* f2;
TTree* tr1;
TTree* tr2;
double range1=10,range2=250;

void DrawXi(){
	TF1* prop = new TF1("prop","[0]*exp(-x/[1])",range1,range2); 
	TF1* LambdaWBG = new TF1("LambdaWBG","GausWithBG",1.,1.2,6);
	TF1* XiWBG = new TF1("XiWBG","GausWithBG",1.2,1.4,6);
	prop->SetParNames("const","c_tau");
	LambdaWBG->SetParNames("LdConst","LdMass","LdWidth","LdBgConst","LdBgMass","LdBgWidth");
	XiWBG->SetParNames("XiConst","XiMass","XiWidth","XiBgConst","XiBgMass","XiBgWidth");
	double par[6];

	gStyle->SetOptFit(000);
	gStyle->SetOptStat(10);
	gStyle->SetTitleFont(132,"x");//,“t”);
	gStyle->SetTitleFont(132,"t");//,“t”);
	gStyle->SetTitleFont(132,"y");//,“t”);
	f1 = new TFile("SelectedEvents.root");
//	f2 = new TFile("TPCInvOld.root");
	f2 = new TFile("TPCInv12.root");
	tr1 = (TTree*)f1->Get("tree");
	tr1->SetBranchAddress("XiMM",&mm1);
	tr1->SetBranchAddress("XiThetaCM",&tcm);
	tr2 = (TTree*)f2->Get("tree");
	tr2->SetBranchAddress("MM",&mm2);
	tr2->SetBranchAddress("InvMLd",&MLd);
	tr2->SetBranchAddress("InvMXiCor",&MXi);
//	tr2->SetBranchAddress("InvMXi",&MXi);
	tr2->SetBranchAddress("FlgLd",&FlgLd);
	tr2->SetBranchAddress("FlgXi",&FlgXi);
	tr2->SetBranchAddress("InTargetLd",&InTargetLd);
	tr2->SetBranchAddress("DistLd",&DistLd);
	tr2->SetBranchAddress("MomLd",&MomLd);
	TH1D* LdHist = new TH1D("Lambda","#Lambda Invariant Mass",50,1.07,1.17);
	TH1D* LdXiHist = new TH1D("LambdawXi","#Lambda w #Xi Invariant Mass",50,1.07,1.17);
	TH1D* XiHist = new TH1D("Xi","#Xi Invariant Mass",100,1.24,1.44);
	TH1D* LdDist = new TH1D("LdPropLength","#Lambda FlightLength/#gamma#beta",20,range1,range2);
	TH1D* LdMMHist = new TH1D("LdMissing Mass","#Lambda tagged Missing Mass",100,1,2);
	TH1D* MMHist = new TH1D("Missing Mass","Missing Mass",100,1,2);
	int ent = tr2->GetEntries();
	int n_ximm=0;
	
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
		bool mmfl = true;
		if(abs(mm2-1.315)>0.1){
	//		continue;
			mmfl = false;
		}
		n_ximm++;
		FlgXi=(FlgXi and !InTargetLd and mmfl);
//		FlgLd=(FlgLd and !InTargetLd );
	//	FlgXi=(FlgXi );// !InTargetLd);
	//	if(FlgXi)LdHist->Fill(MLd);
		if(FlgLd)LdHist->Fill(MLd);
		if(FlgXi)XiHist->Fill(MXi);
		if(FlgXi)LdXiHist->Fill(MLd);
		if(FlgXi)LdDist->Fill(DistLd/MomLd*1.115);
	}
	int ldpeak = LambdaWBG->GetMaximum();
	LambdaWBG->SetParLimits(0,0.7*ldpeak,ldpeak);
	LambdaWBG->SetParLimits(1,1.11,1.13);
	LambdaWBG->SetParLimits(2,0.001,0.01);
	LambdaWBG->SetParLimits(3,5,25);
	LambdaWBG->SetParLimits(4,1.10,1.15);
	LambdaWBG->SetParLimits(5,0.01,0.1);
	XiWBG->SetParLimits(0,10,50);
	XiWBG->SetParLimits(1,1.31,1.33);
	XiWBG->SetParLimits(2,0.003,0.03);
	XiWBG->SetParLimits(3,5,20);
	XiWBG->SetParLimits(4,1.28,1.36);
	XiWBG->SetParLimits(5,0.01,0.1);
	cout<<"Xi MM Window: "<<n_ximm<<endl;
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	TCanvas* c4 = new TCanvas("c4","c4",600,600);
	TCanvas* c5 = new TCanvas("c5","c5",600,600);
	TCanvas* c1 = new TCanvas("c1","c1",1800,600);
	c1->Divide(3,1);
	TF1* fgau = new TF1("fgau","Gaussianf",1.1,1.4,3);
	TF1* fgau2 = new TF1("fgau2","Gaussianf",1.1,1.4,3);
	TF1* fgau3 = new TF1("fgau3","Gaussianf",1.1,1.4,3);
	TF1* fgau4 = new TF1("fgau4","Gaussianf",1.1,1.4,3);
	bool drld=true,drxi=true,drprop=true;
	c1->cd(1);
//	LdHist->Draw();
//	LdXiHist->SetLineColor(kRed);
	LdXiHist->Draw("");
	LdXiHist->Fit("LambdaWBG","");
	LdXiHist->GetXaxis()->SetTitle("M_{p #pi^{-}} [GeV/c^{2}]");
	LdXiHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	LdXiHist->GetXaxis()->SetNdivisions(5,false);
	LdXiHist->GetYaxis()->SetNdivisions(10);
	for(int i=0;i<6;++i){
		par[i]= LambdaWBG->GetParameter(i);
	}
	for(int i=0;i<3;++i){
		fgau->SetParameter(i,par[i]);
		fgau2->SetParameter(i,par[i+3]);
	}
	fgau->SetLineColor(kGreen);
	fgau2->SetLineColor(kBlack);
	fgau->Draw("same");
	fgau2->Draw("same");
		double LdMass =LambdaWBG->GetParameter(1);
		double LdWidth =LambdaWBG->GetParameter(2);
//		TText* LdLabel = new TText(1.12,40,Form(" %.0f \n+- %.0f MeV/c2",1000*LdMass,1000*LdWidth));
//		LdLabel->Draw("same");
//	c2->cd();
	c1->cd(2);
	XiHist->Draw();
	XiHist->Fit("XiWBG","");
	XiHist->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}} [GeV/c^{2}]");
	XiHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	XiHist->GetXaxis()->SetNdivisions(5,false);
	XiHist->GetYaxis()->SetNdivisions(10);
	for(int i=0;i<6;++i){
		par[i]= XiWBG->GetParameter(i);
	}
	for(int i=0;i<3;++i){
		fgau3->SetParameter(i,par[i]);
		fgau4->SetParameter(i,par[i+3]);
	}
	fgau3->SetLineColor(kGreen);
	fgau4->SetLineColor(kBlack);
	fgau3->Draw("same");
	fgau4->Draw("same");
	if(drxi){
		double XiMass =XiWBG->GetParameter(1);
		double XiWidth =XiWBG->GetParameter(2);
//		TText* XiLabel = new TText(1.33,20,Form(" %.0f +- %.0f MeV/c2",1000*XiMass,1000*XiWidth));
//		XiLabel->Draw("same");
	}
	c1->cd(3);
	prop->SetParLimits(1,50,90);
	LdDist->Draw();
	LdDist->Fit("prop","");
	if(drprop){
		double ctau = prop->GetParameter(1);
		TText* PropLabel;// = new TText("prh");
//		PropLabel = new TText(150,60,Form("c_tau = %.1f mm",ctau));
//		PropLabel->Draw("same");
	}
	LdDist->GetXaxis()->SetTitle("#Lambda flight length /#gamma#beta [mm]");
	LdDist->GetYaxis()->SetTitle("Entries / 12 mm");
	LdDist->GetXaxis()->SetNdivisions(10);
	LdDist->GetYaxis()->SetNdivisions(10);
	
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
	if(FlgLd and abs(MLd - LdMass)<3*LdWidth)LdMMHist ->Fill(mm2);
	MMHist ->Fill(mm2);
	}
#if 0
	c4->cd();
	MMHist->Draw();
	MMHist->GetXaxis()->SetTitle("MissMass [GeV / c^{2}]");
	MMHist->GetYaxis()->SetTitle("Entries / 10 MeV/c^{2}");
	MMHist->GetXaxis()->SetNdivisions(10);
	MMHist->GetYaxis()->SetNdivisions(10);
	c5->cd();
	LdMMHist->Draw();
	LdMMHist->GetXaxis()->SetTitle("MissMass [GeV / c^{2}]");
	LdMMHist->GetYaxis()->SetTitle("Entries / 10 MeV/c^{2}");
	LdMMHist->GetXaxis()->SetNdivisions(10);
	LdMMHist->GetYaxis()->SetNdivisions(10);
#endif



}



























double LpG(double* x,double*par){
  double gau = par[0]*TMath::Gaus(x[0],par[1],par[2]);
  double lad = par[3]*TMath::Landau(x[0],par[4],par[5]);
  return gau+lad;
}

void DrawDifCr(){
	

	TGraphErrors* g1;
	TGraphErrors* g2;
	TH1D* hist1 = new TH1D("Xi Count vs cos(#theta_{CM})","Xi Count vs cos(#theta_{CM})",10,0.75,1);
	hist1->GetXaxis()->SetTitle("cos(#theta_{CM})");
	hist1->GetYaxis()->SetTitle("Entry [ / 0.025]");
	TH1D* hist2 = new TH1D("Xi* Count vs cos(#theta_{CM})","Xi* Count vs cos(#theta_{CM})",10,0.75,1);
	hist1->GetXaxis()->SetTitle("cos(#theta_{CM})");
	hist1->GetYaxis()->SetTitle("Entry [ / 0.025]");
	int ent = tr1->GetEntries();	
	for(int i=0;i<ent;++i){
		tr1->GetEntry(i);
		if(abs(mm1-1.32657)<0.0145*3) hist1->Fill(cos(tcm/180*PI));
		if(abs(mm1-1.540)<0.0135*3) hist2->Fill(cos(tcm/180*PI));
	}
	double xic[10],xisc[10];
	for(int i=0;i<10;++i){
		xic[i]=hist1->GetBinContent(i+1);
		xisc[i]=hist2->GetBinContent(i+1);
		cout<<Form("%d : %f",i,xic[i])<<endl; 
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	hist1->Draw();
	c1->cd(2);
	hist2->Draw();
}
void FitIM(){
	TH1D* h = new TH1D("h","h",100,1,2);
	int ent = tr2->GetEntries();	
	for(int i=0;i<ent;++i){
		tr2->GetEntry(i);
		h->Fill(MLd);
	}
	  TF1* flg = new TF1("flgaus","LpG",1.07,1.6,6);
	h->Draw();
	flg->SetParLimits(0,100,120);
	flg->SetParLimits(1,1.1,1.12);
	flg->SetParLimits(2,0.01,0.03);
	flg->SetParLimits(4,1,1.2);
	flg->SetParLimits(5,0,1);
	h->Fit("flgaus","R");
	h->GetXaxis()->SetTitle("InvMass [GeV/c2]");
	h->SetTitle("#Lambda Invariant Mass");

}

void DrawXi(int dum){
	TH1D* hist1 = new TH1D("KuramaMM","KuramaMM",200,1,2);
	hist1->GetXaxis()->SetTitle("MissMass [GeV / c2]");
	hist1->GetYaxis()->SetTitle("Entry [ /5 MeV]");
	TH1D* hist2 = new TH1D("KuramaMMCut","KuramaMMCut",200,1,2);
	int ent = tr1->GetEntries();
	for(int i=0;i<ent;++i){
		tr1->GetEntry(i);
		hist1->Fill(mm1);
	}
	ent = tr2->GetEntries();
	for(int i=0;i<ent;++i){
		tr2->GetEntry(i);
		if(abs(MLd-1.12)<0.044)hist2->Fill(mm2);
		//		if(MLd>0)	hist2->Fill(mm2);
	}

	hist2->SetLineColor(kRed);
	TCanvas*c1 = new TCanvas("c1","c1",1200,800);
	c1->cd();
	hist1->Draw();
	hist2->Draw("same");
}
void DrawMMAngular(){
	double mm,theta;
	TFile* f1 = new TFile("SelectedEvents.root");
	TTree* tr1 = (TTree*)f1->Get("tree");
	tr1->SetBranchAddress("XiMM",&mm);
	tr1->SetBranchAddress("XiThetaCM",&theta);
	TH2D* hist1 = new TH2D("MM vs cos(#theta_{CM})","MM vs cos(#theta_{CM})",50,0.5,1,50,1.2,1.7);
	int ent = tr1->GetEntries();
	gStyle->SetOptStat(0000);
	hist1->GetXaxis()->SetTitle("cos(#theta_{CM})");
	hist1->GetYaxis()->SetTitle("MissMass[GeV/c^2]");
	for(int i=0;i<ent;++i){
		tr1->GetEntry(i);
		hist1->Fill(cos(theta/180 * PI),mm);
	}
	hist1->Draw("colz");
}
