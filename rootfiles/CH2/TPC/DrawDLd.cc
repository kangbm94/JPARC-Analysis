#include "../../../KKManager.hh"
double PI = acos(-1);
double mm1,mm2,MLd,tcm,MXi,DistLd,MomLd;
TFile* f1;TFile* f2;
TTree* tr1;
TTree* tr2;
double range1=10,range2=250;
vector<double>* InvMLd = new vector<double>;
vector<bool>* InTargetLd = new vector<bool>;
int nLd=0;
void DrawDLd(){
	TF1* prop = new TF1("prop","[0]*exp(-x/[1])",range1,range2); 
	TF1* LambdaWBG = new TF1("LambdaWBG","GausWithBG",1.,1.2,6);
	TF1* Lambda2WBG = new TF1("Lambda2WBG","GausWithBG",1.,1.2,6);
	prop->SetParNames("const","c_tau");
	LambdaWBG->SetParNames("LdConst","LdMass","LdWidth","LdBgConst","LdBgMass","LdBgWidth");
	Lambda2WBG->SetParNames("Ld2Const","Ld2Mass","Ld2Width","Ld2BgConst","Ld2BgMass","Ld2BgWidth");
	double par[6];

	gStyle->SetOptFit(000);
	gStyle->SetOptStat(10);
	gStyle->SetTitleFont(132,"x");//,“t”);
	gStyle->SetTitleFont(132,"t");//,“t”);
	gStyle->SetTitleFont(132,"y");//,“t”);
	f1 = new TFile("SelectedEvents.root");
//	f2 = new TFile("TPCInvOld.root");
	f2 = new TFile("DoubleLambda.root");
	tr1 = (TTree*)f1->Get("tree");
	tr1->SetBranchAddress("XiMM",&mm1);
	tr1->SetBranchAddress("XiThetaCM",&tcm);
	tr2 = (TTree*)f2->Get("tree");
	tr2->SetBranchAddress("MM",&mm2);
	tr2->SetBranchAddress("InvMLd",&InvMLd);
//	tr2->SetBranchAddress("InvMXi",&MXi);
	tr2->SetBranchAddress("InTargetLd",&InTargetLd);
	tr2->SetBranchAddress("nLd",&nLd);
	TH1D* LdHist = new TH1D("Lambda","#Lambda Invariant Mass",50,1.07,1.17);
	TH1D* LdHist2 = new TH1D("Lambda2","#Lambda Invariant Mass",50,1.07,1.17);
	TH1D* LdDist = new TH1D("LdPropLength","#Lambda FlightLength/#gamma#beta",20,range1,range2);
	TH1D* LdMMHist = new TH1D("LdMissing Mass","2#Lambda tagged Missing Mass",100,1,2);
	TH1D* MMHist = new TH1D("Missing Mass","Missing Mass",100,1,2);
	int ent = tr2->GetEntries();
	int n_ximm=0;
	int ndld = 0;	
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
		double mld=0,mld2=0;
		if(nLd>1){
			mld=InvMLd->at(0);
			LdHist->Fill(InvMLd->at(0));
		}
		if(nLd>1){
			mld2=InvMLd->at(1);
		}
		if(mld == mld2 and mld2!=0){
			cout<<"Warning! double count!"<<endl;
			cout<<mld2<<endl;
			continue;
		}
		if(nLd>1){
			LdHist2->Fill(InvMLd->at(1));
			LdHist2->Fill(InvMLd->at(0));
		}
	}
	int ldpeak = LdHist->GetMaximum();
	LambdaWBG->SetParLimits(0,0.7*ldpeak,ldpeak);
	LambdaWBG->SetParLimits(1,1.11,1.13);
	LambdaWBG->SetParLimits(2,0.001,0.01);
	LambdaWBG->SetParLimits(3,0,5);
	LambdaWBG->SetParLimits(4,1.10,1.15);
	LambdaWBG->SetParLimits(5,0.01,0.1);
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	TF1* fgau = new TF1("fgau","Gaussianf",1.1,1.4,3);
	TF1* fgau2 = new TF1("fgau2","Gaussianf",1.1,1.4,3);
	TF1* fgau3 = new TF1("fgau3","Gaussianf",1.1,1.4,3);
	TF1* fgau4 = new TF1("fgau4","Gaussianf",1.1,1.4,3);
	bool drld=true,drxi=true,drprop=true;
	c1->cd();
	LdHist->Draw();
//	LdXiHist->SetLineColor(kRed);
	LdHist->Fit("LambdaWBG","0");
	LdHist->GetXaxis()->SetTitle("M_{p #pi^{-}} [GeV/c^{2}]");
	LdHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	LdHist->GetXaxis()->SetNdivisions(5,false);
	LdHist->GetYaxis()->SetNdivisions(10);
	for(int i=0;i<6;++i){
		par[i]= LambdaWBG->GetParameter(i);
	}
	for(int i=0;i<3;++i){
		fgau->SetParameter(i,par[i]);
		fgau2->SetParameter(i,par[i+3]);
	}
	fgau->SetLineColor(kRed);
	fgau2->SetLineColor(kBlack);
//	fgau->Draw("same");
//	fgau2->Draw("same");
	double LdMass =LambdaWBG->GetParameter(1);
	double LdWidth =LambdaWBG->GetParameter(2);
//		TText* LdLabel = new TText(1.12,40,Form(" %.0f \n+- %.0f MeV/c2",1000*LdMass,1000*LdWidth));
//		LdLabel->Draw("same");
	c2->cd();
	LdHist2->GetXaxis()->SetTitle("M_{p #pi^{-}} [GeV/c^{2}]");
	LdHist2->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	LdHist2->GetXaxis()->SetNdivisions(5,false);
	LdHist2->GetYaxis()->SetNdivisions(10);
	int ldpeak2 = LdHist2 ->GetMaximum();
	Lambda2WBG->SetParLimits(0,0.5*ldpeak2,ldpeak2);
	Lambda2WBG->SetParLimits(1,1.105,1.13);
	Lambda2WBG->SetParLimits(2,0.001,0.01);
	Lambda2WBG->SetParLimits(3,0,3);
	Lambda2WBG->SetParLimits(4,1.10,1.15);
	Lambda2WBG->SetParLimits(5,0.01,0.1);
	LdHist2->Fit("Lambda2WBG","0");
	for(int i=0;i<6;++i){
		par[i]= Lambda2WBG->GetParameter(i);
	}
	double mldi2=par[1];
	double mldw2=par[2];
	for(int i=0;i<3;++i){
		fgau3->SetParameter(i,par[i]);
		fgau4->SetParameter(i,par[i+3]);
	}
	LdHist2->Draw();
	fgau3->SetLineColor(kRed);
	fgau4->SetLineColor(kBlack);
	fgau3->Draw("same");
//	fgau4->Draw("same");
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
		if(nLd!=2) continue;
		double mld1 = InvMLd->at(0);	
		double mld2 = InvMLd->at(1);	
		if(abs(mld2-mldi2)<3*mldw2 or abs(mld1-mldi2)<3*mldw2)LdMMHist->Fill(mm2);
	}
	c3->cd();
	LdMMHist->Draw();
	LdMMHist->GetXaxis()->SetTitle("MissMass [GeV / c^{2}]");
	LdMMHist->GetYaxis()->SetTitle("Entries / 10 MeV/c^{2}");
	LdMMHist->GetXaxis()->SetNdivisions(10);
	LdMMHist->GetYaxis()->SetNdivisions(10);
	/*
	MMHist->Draw();
	MMHist->GetXaxis()->SetTitle("MissMass [GeV / c^{2}]");
	MMHist->GetYaxis()->SetTitle("Entries / 10 MeV/c^{2}");
	MMHist->GetXaxis()->SetNdivisions(10);
	MMHist->GetYaxis()->SetNdivisions(10);
	LdMMHist->Draw();
*/



}



