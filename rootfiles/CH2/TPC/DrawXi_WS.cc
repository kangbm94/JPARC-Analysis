#include "../../../KKManager.hh"
#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
double PI = acos(-1);
int evnum, runnum;
int evnumX, runnumX;
double mm1,mm2,MLd,MLdCor,tcm,MXi,DistLd,MomLd,MXiCor,mmpi0,MissingM,CheckSumM,CheckSumMCor;
double MissingKpLd,MissingpKpLd;
double CosCM,CosLab;
bool KmCor,KpCor;
bool FlgLd,FlgXi,InTargetLd,isGood;
TFile* f1;TFile* f2;
TTree* tr1;
TTree* tr2;
vector<double>ldlen;
double range1=10,range2=250;
double tsize = 0.05;
double acpt[40];
double acpt_kn[12] = {
	0.0135391, 0.031935, 0.0451718, 0.068755,
	0.0766182,0.0864692,0.0916622,0.101086,
	0.121607,0.159989,0.212732,0.242154};
double angle[40];
void fcn (int& npar, double* grad, double& fval, double* par, int flag){
	double nll = 0;
	for(double len:ldlen){
		double expv = 1./par[1]*exp(-len/par[1]); 
		nll+=-log(expv);
	}
	fval = nll;
};

void DrawXi_WS(){
	for(int i=0;i<40;++i){
		angle[i]= -0.975+ i * 0.05;
		acpt[i]= 0;
		if(i > 27){
			acpt[i] = acpt_kn[i - 28];
		}
	}
	SetStyle();
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.4);
	gStyle->SetStatH(0.5);
	ldlen.clear();
	TF1* prop = new TF1("prop","[0]*exp(-x/[1])",range1,range2); 
	prop->SetParameter(0,70);
	prop->SetParameter(1,77);
	TH2D* MMScat = new TH2D("KuramaMM vs (Kurama-Xi)MM","KuramaMM vs (Kurama-Xi)MM",20,1.2,1.7,20,-0.5,0.5);
	TF1* LambdaWBG = new TF1("LambdaWBG","GausWithBG",1.,1.2,6);
	TF1* XiWBG = new TF1("XiWBG","GausWithBG",1.2,1.4,6);
	prop->SetParNames("const","c_tau");
	LambdaWBG->SetParNames("LdConst","LdMass","LdWidth","LdBgConst","LdBgMass","LdBgWidth");
	XiWBG->SetParNames("XiConst","XiMass","XiWidth","XiBgConst","XiBgMass","XiBgWidth");
	double par[6];

	gStyle->SetOptFit(000);
	gStyle->SetOptStat(10);
	gStyle->SetTitleFontSize(0.1);
	gStyle->SetTitleFont(132,"x");//,“t”);
	gStyle->SetTitleFont(132,"t");//,“t”);
	gStyle->SetTitleFont(132,"y");//,“t”);
																//	f2 = new TFile("TPCInvM_cd15_Cor2ndWPrio_TrackRes.root");	
																//	f2 = new TFile("TPCInvM_cd15_Cor2ndWPrio.root");	
																//	f2 = new TFile("TPCInvM_cd15_NoCorWOPrio.root");	
																//	f2 = new TFile("TPCInvM_cd25_HelMom.root");	
	//f2 = new TFile("TPCInvM_cd7_WS_firstT.root");	
//	f2 = new TFile("TPCInvMKinFit_cd7_WS.root");	
//	f2 = new TFile("TPCInvM_cd7_WS.root");	
	f2 = new TFile("TPCInvMScaled_cd7_WS.root");	
	//	f2 = new TFile("InvM/TPCInvM05641_cd7_WS.root");	
	tr2 = (TTree*)f2->Get("tree");
	tr2->SetBranchAddress("runnum",&runnum);
	tr2->SetBranchAddress("evnum",&evnum);
	tr2->SetBranchAddress("MM",&mm2);
	tr2->SetBranchAddress("mmpi0",&mmpi0);
	tr2->SetBranchAddress("KmCor",&KmCor);
	tr2->SetBranchAddress("KpCor",&KpCor);
	tr2->SetBranchAddress("InvMLd",&MLd);
	tr2->SetBranchAddress("InvMXi",&MXi);
	//	tr2->SetBranchAddress("InvMXi",&MXi);
	tr2->SetBranchAddress("FlgLd",&FlgLd);
	tr2->SetBranchAddress("FlgXi",&FlgXi);
	tr2->SetBranchAddress("InTargetLd",&InTargetLd);
	tr2->SetBranchAddress("DistLd",&DistLd);
	tr2->SetBranchAddress("MomLd",&MomLd);
	tr2->SetBranchAddress("InvMLdCor",&MLdCor);
	tr2->SetBranchAddress("MissingM",&MissingM);
	tr2->SetBranchAddress("MissingKpLd",&MissingKpLd);
	tr2->SetBranchAddress("MissingpKpLd",&MissingpKpLd);
	tr2->SetBranchAddress("CheckSumM",&CheckSumM);
	tr2->SetBranchAddress("CheckSumMCor",&CheckSumMCor);
	tr2->SetBranchAddress("isGood",&isGood);
	tr2->SetBranchAddress("CosCM",&CosCM);
	tr2->SetBranchAddress("CosLab",&CosLab);
	TH1D* LdHist = new TH1D("Lambda","Invariant Mass(p #pi)",50,1.07,1.17);
	TH1D* LdCorHist = new TH1D("Lambda","#Lambda Invariant MassCor",500,1.07,1.17);
	TH1D* LdXiHist = new TH1D("LambdawXi","#Lambda w #Xi Invariant Mass",50,1.07,1.17);
	TH1D* LdXiPropHist = new TH1D("LambdaPropwXi","#Lambda w #Xi Invariant Mass Prop",50,1.07,1.17);
	TH1D* XiHist = new TH1D("Xi","#Xi Invariant Mass",100,1.24,1.44);
	TH1D* XiCorHist = new TH1D("XiCor","Invariant Mass(#Lambda,#pi)",100,1.24,1.44);
	TH1D* XiPropHist = new TH1D("XiProp","#Xi Invariant Mass Prop",100,1.24,1.44);
	TH1D* LdDist = new TH1D("LdPropLength","#Lambda FlightLength/#gamma#beta",20,range1,range2);
	TGraph* LdFlGr;
	TH1D* LdMMHist = new TH1D("LdMissing Mass","#Lambda tagged Missing Mass",100,1,2);
	TH1D* MMHist = new TH1D("Missing Mass","Missing Mass",100,1,2);
	TH2D* MMPi0 = new TH2D("MMPi0","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X}",100,-1.2,0.4,100,1.2,1.7);
	TH2D* MMPi0Cor = new TH2D("MMPi0Cor","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X}",100,-1.2,0.4,100,1.2,1.7);
	TH2D* MMPi0Vtx = new TH2D("MMPi0Vtx","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X}",100,-1.5,1,100,1.1,1.8);
	TH1D* MPi0 = new TH1D("MPi0","Sum_{TPC}-X",50,0,0.3);
	TH1D* MMTPC = new TH1D("MMTPC","Sum_{TPC}-X",200,-2,2);
	TH2D* MMKmKpLd2D = new TH2D("MissingMass(Km,KpLd):KuramaMM","MissingMass(Km,KpLd):KuramaMM",100,1.1,1.8,100,-1,2);
	TH1D* MMKmKpLd = new TH1D("MissingMass(Km,KpLd)","MissingMass(Km,KpLd)",100,-1,2);
	TH2D* MMKmpKpLd2D = new TH2D("MissingMass(Kmp,KpLd):KuramaMM","MissingMass(Kmp,KpLd):KuramaMM",100,1.1,1.8,100,-1,2);
	TH1D* MMKmpKpLd = new TH1D("MissingMass(Kmp,KpLd)","MissingMass(Kmp,KpLd)",100,-1,2);
	int ent = tr2->GetEntries();
	int n_ximm=0;
	int ent2 = tr2->GetEntries();
	TFile* file = new TFile("pkkXiX.root","recreate");
	TTree* tree = new TTree("tree","tree");
	tree->Branch("mmKurama",&mm2);
	tree->Branch("mmTPC",&CheckSumM);
	TH1D* h1 = new TH1D("h1","h1",20,0,1.001);
	double gr_ang[12];
	double gr_ang_er[12];
	double gr_ent[12];
	double gr_err[12];
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
		bool mmfl = true;
		if(isGood and abs(MLd-1.115)<0.015 ){
			MMPi0->Fill(CheckSumM,mm2);
			if(1){
				MMPi0Cor->Fill(CheckSumM,mm2);
			}
			else{
				MMPi0Vtx->Fill(CheckSumMCor,mm2);
			}
			MMTPC->Fill(CheckSumM);
			tree->Fill();
		}
		if(isGood and abs(CheckSumM-0.15)<0.15 and abs(mm2-1.55)<0.1){
			MPi0->Fill(CheckSumM);
		}
		if(abs(mm2-1.321)>0.1){
			//			continue;
			mmfl = false;
		}
		//		if(runnum != 5641) continue;
		n_ximm++;
		bool flFlag = false;
		//		FlgXi=(FlgXi and !InTargetLd and mmfl);
		FlgXi=(FlgXi and mmfl);
		//		FlgXi=(FlgXi and mmfl);
		//		FlgLd=(FlgLd and !InTargetLd );
		//	FlgXi=(FlgXi );// !InTargetLd);
		//	if(FlgXi)LdHist->Fill(MLd);
		if(DistLd > 11.) flFlag = true;
		//		if(FlgLd and  mmfl and !InTargetLd and FlgXi)LdHist->Fill(MLd);
		if(FlgXi )LdHist->Fill(MLd);
		if(FlgLd &&( !FlgXi or 1)){
			if(MissingpKpLd<0 or 1){
				MMKmKpLd->Fill(MissingKpLd);
				MMKmKpLd2D->Fill(mm2,MissingKpLd);
			}
			MMKmpKpLd->Fill(MissingpKpLd);
			MMKmpKpLd2D->Fill(mm2,MissingpKpLd);
		}
		//		if(FlgLd and mmfl)LdHist->Fill(MLd);
		if(FlgXi)XiHist->Fill(MXi);
		if(FlgXi)XiCorHist->Fill(MXi);
		if(FlgXi and !InTargetLd )XiPropHist->Fill(MXi);
		//		if(FlgXi and !InTargetLd)LdXiHist->Fill(MLd);
		//		if(FlgLd)LdXiHist->Fill(MLd);
		if(FlgXi)
		{
			LdCorHist->Fill(MLdCor);
			//			LdXiPropHist->Fill(MLd);
		}
		if(FlgXi and !InTargetLd and flFlag){
			XiPropHist->Fill(MXi);
		}
		if(FlgXi and !InTargetLd and flFlag){
			ldlen.push_back(DistLd/MomLd*1.115);
			LdDist->Fill(DistLd/MomLd*1.115);
		}
		if(abs(MissingM - 1.535) < 0.12){
			h1->Fill(CosCM);
		}
	}
	for(int i=0; i<20;++i){
		double ent = h1->GetBinContent(i+1);
		double cor = acpt[i + 20];
		if(cor==0) continue;
		else{
			gr_ang[i-8] = 0.425 + 0.05*(i-8);
			gr_ang_er[i-8] = 0.05/sqrt(12);
//			gr_ent[i-8] = ent / cor;
			ent *= 0.0886; 
			gr_ent[i-8] = ent *0.6;
			gr_err[i-8] = sqrt(ent*0.6) ;
			cout<<gr_ang[i-8]<<" ," <<acpt[i+20]<<endl;
		}


	}
	int ldpeak = LdHist->GetMaximum();
	LambdaWBG->SetParLimits(0,0.7*ldpeak,ldpeak);
	LambdaWBG->SetRange(1.08,1.12);
	LambdaWBG->SetParLimits(1,1.11,1.13);
	LambdaWBG->SetParLimits(2,0.001,0.01);
	LambdaWBG->SetParLimits(3,0,ldpeak* 0.2);
	LambdaWBG->SetParLimits(4,1.10,1.15);
	LambdaWBG->SetParLimits(5,0.005,0.03);
	int xipeak = XiCorHist->GetMaximum();
	XiWBG->SetParLimits(0,xipeak/2,xipeak*1.3);
	XiWBG->SetParLimits(1,1.31,1.33);
	XiWBG->SetParLimits(2,0.003,0.03);
	XiWBG->SetParLimits(3,1,xipeak/2);
	XiWBG->SetParLimits(4,1.30,1.34);
	XiWBG->SetParLimits(5,0.02,0.07);
	cout<<"Xi MM Window: "<<n_ximm<<endl;
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	TCanvas* c4 = new TCanvas("c4","c4",600,600);
	TCanvas* c5 = new TCanvas("c5","c5",600,600);
	TCanvas* c1 = new TCanvas("c1","c1",1800,600);
	c1->Divide(3,1);
	TF1* fgau = new TF1("fgau","Gaussianf",1.07,1.4,3);
	TF1* fgau2 = new TF1("fgau2","Gaussianf",1.07,1.4,3);
	TF1* fgau3 = new TF1("fgau3","Gaussianf",1.07,1.4,3);
	TF1* fgau4 = new TF1("fgau4","Gaussianf",1.07,1.4,3);
	bool drld=true,drxi=true,drprop=true;
	c1->cd(1);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	LdHist->Draw();
	//	LdXiHist->SetLineColor(kRed);
	//	LdCorHist->Draw(""); //	LdCorHist->SetLineColor(kBlack); //	LdXiHist->Draw("");
	//	LdXiPropHist->Draw("same");
	//	LdXiPropHist->SetLineColor(kRed);
	//	LdHist->Fit("LambdaWBG","R");
	LdHist->GetXaxis()->SetTitle("M_{p #pi^{-}} [GeV/c^{2}]");
	LdHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	LdHist->GetXaxis()->SetNdivisions(5,false);
	LdHist->GetYaxis()->SetNdivisions(10);
	LdHist->GetYaxis()->SetTitleSize(tsize);
	LdHist->GetXaxis()->SetTitleSize(tsize);
	LdHist->SetStats(0);
	cout<<"Gaus"<<endl;
	//	LdHist->Fit("gaus");
	TLatex* LdEnt = new TLatex(1.125,0.9*(LdHist->GetMaximum()),Form("Entries = %g",LdHist->GetEffectiveEntries()));
	//	LdEnt->Draw();
	for(int i=0;i<6;++i){
		par[i]= LambdaWBG->GetParameter(i);
	}
	for(int i=0;i<3;++i){
		fgau->SetParameter(i,par[i]);
		fgau2->SetParameter(i,par[i+3]);
	}
	fgau->SetLineColor(kGreen);
	fgau2->SetLineColor(kBlack);
	//	fgau->Draw("same");
	//	fgau2->Draw("same");
	double LdMass =LambdaWBG->GetParameter(1);
	double LdWidth =LambdaWBG->GetParameter(2);
	//		TText* LdLabel = new TText(1.12,40,Form(" %.0f \n+- %.0f MeV/c2",1000*LdMass,1000*LdWidth));
	//		LdLabel->Draw("same");
	//	c2->cd();
	c1->cd(2);
	XiCorHist->Draw();
	//	XiHist->Draw("same");
	//	XiPropHist->Draw("same");
	//	XiPropHist->Draw("");
	//	XiPropHist->SetLineColor(kRed);
	//	XiCorHist->SetLineColor(kBlack);
	//	XiCorHist->Fit("XiWBG","");
	XiCorHist->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}} [GeV/c^{2}]");
	XiCorHist->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}} [GeV/c^{2}]");
	XiCorHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	XiCorHist->GetXaxis()->SetNdivisions(5,false);
	XiCorHist->GetXaxis()->SetTitleSize(tsize);;
	XiCorHist->GetYaxis()->SetNdivisions(10);
	XiCorHist->GetYaxis()->SetTitleSize(tsize);;
	for(int i=0;i<6;++i){
		par[i]= XiWBG->GetParameter(i);
	}
	for(int i=0;i<3;++i){
		fgau3->SetParameter(i,par[i]);
		fgau4->SetParameter(i,par[i+3]);
	}
	fgau3->SetLineColor(kGreen);
	fgau4->SetLineColor(kBlack);
	//	fgau3->Draw("same");
	//	fgau4->Draw("same");
	if(drxi){
		double XiMass =XiWBG->GetParameter(1);
		double XiWidth =XiWBG->GetParameter(2);
		//		TText* XiLabel = new TText(1.33,20,Form(" %.0f +- %.0f MeV/c2",1000*XiMass,1000*XiWidth));
		//		XiLabel->Draw("same");
	}
	XiCorHist->SetStats(0);
	TLatex* XiEnt = new TLatex(1.33,0.9*(XiCorHist->GetMaximum()),Form("Entries = %g",XiCorHist->GetEffectiveEntries()));
	//	XiEnt->Draw();
	c1->cd(3);
	//	prop->SetParLimits(1,50,90);
	LdDist->Draw("E1");
	//	TLatex* LdDistEnt = new TLatex(125,40,Form("Entries = %g",LdDist->GetEffectiveEntries()));
	//	LdDistEnt->Draw();

	int ierflg=0;
	TMinuit* min = new TMinuit(2);	
	min->SetFCN(fcn);	
	min->mnparm(0,"const",70,0.01,0,0,ierflg);
	min->mnparm(1,"ctau",77,0.001,0,0,ierflg);
	double arglist[2]={1000,0.01};
	min->mnexcm("MIGRAD",arglist,2,ierflg);
	double p0,p0e,p1,p1e;
	min->GetParameter(0,p0,p0e);
	min->GetParameter(1,p1,p1e);
	prop->SetParameter(0,p0);
	prop->SetParameter(1,p1);
	//	prop->Draw("same");
	prop->SetRange(35,200);
	LdDist->Fit("prop","R");
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
	LdDist->GetYaxis()->SetTitleSize(tsize);
	LdDist->GetXaxis()->SetTitleSize(tsize);
	LdDist->SetStats(0);
	cout<<"LdEnt: "<<LdXiHist->GetEffectiveEntries()<<endl;
	cout<<"XiEnt: "<<XiCorHist->GetEffectiveEntries()<<endl;
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
		if(FlgLd and abs(MLd - LdMass)<3*LdWidth)LdMMHist ->Fill(mm2);
		MMHist ->Fill(mm2);
	}
	c1->SaveAs("XiReconHist.png");
	TCanvas* c6 = new TCanvas("c6","c6",900,900);
	gPad->SetMargin(0.15,0.10,0.15,0.10);
	MMPi0Cor->Draw("colz");
	MMPi0Cor->SetStats(0);
	MMPi0Cor->GetXaxis()->SetTitle("MM_{p(K^{-},K^{+}#Xi^{-})X}[GeV/c^{2}]");
	MMPi0Cor->GetYaxis()->SetTitle("MM_{p(K^{-},K^{+})X}[GeV/c^{2}]");
//	MMPi0Cor->GetXaxis()->SetNdivisions(20);
//	MMPi0Cor->GetYaxis()->SetNdivisions(10);
//	MMPi0Cor->GetYaxis()->SetTitleSize(tsize*0.8);
//	MMPi0Cor->GetXaxis()->SetTitleSize(tsize);
	TLine* l1 = new TLine(0,1.45,0,1.65);
	TLine* l2 = new TLine(0,1.65,0.3,1.65);
	TLine* l3 = new TLine(0.3,1.45,0.3,1.65);
	TLine* l4 = new TLine(0,1.45,0.3,1.45);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(kDashed);
	l1->SetLineWidth(4);
	l2->SetLineColor(kRed);
	l2->SetLineStyle(kDashed);
	l2->SetLineWidth(4);
	l3->SetLineColor(kRed);
	l3->SetLineStyle(kDashed);
	l3->SetLineWidth(4);
	l4->SetLineColor(kRed);
	l4->SetLineStyle(kDashed);
	l4->SetLineWidth(4);

	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");

	TLatex* t = new TLatex(0.,1.65,"#pi^{0}");
	t->SetTextColor(kRed);
	t->SetTextSize(0.2);
	//	t->Draw();
	c6->SaveAs("MMPi0.png");

	TCanvas* c7 = new TCanvas("c7","c7",900,900);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	MPi0->SetStats(0);
	TLatex* stats = new TLatex(0.2,7,Form("Entries : %g", MPi0->GetEffectiveEntries()));
	MPi0->Draw();
	//	stats->Draw();
	MPi0->GetXaxis()->SetTitle("M_{#Xi^{-} -X}[GeV/c]");
	TString titleY = Form(" Entries /%g MeV/{c^2}",MPi0->GetXaxis()->GetBinWidth(5)*1000);
	MPi0->SetLineWidth(3);
	MPi0->SetLineColor(kBlack);
	MPi0->GetYaxis()->SetTitle(titleY);
	MPi0->GetYaxis()->SetTitleSize(tsize);
	MPi0->GetYaxis()->SetNdivisions(5);
	MPi0->GetXaxis()->SetNdivisions(10);
	c7->SaveAs("MPi0.png");
	TCanvas* c8 = new TCanvas("c8","c8",900,900);
	MMTPC->SetStats(0);
	TLatex* stats2 = new TLatex(0.2,100,Form("Entries : %g", MMTPC->GetEffectiveEntries()));
	MMTPC->Draw();
	//	stats2->Draw();
	MMTPC->GetXaxis()->SetTitle("M_[#Lambda #pi^{-}]-X[GeV/c]");
	TString titleY2 = Form(" Entries /%g MeV/{c^2}",MMTPC->GetXaxis()->GetBinWidth(5)*1000);
	MMTPC->SetLineWidth(3);
	MMTPC->SetLineColor(kBlack);
	MMTPC->GetYaxis()->SetTitle(titleY2);
	MMTPC->GetYaxis()->SetNdivisions(5);
	MMTPC->GetXaxis()->SetNdivisions(10);
	c8->SaveAs("MMTPC.png");
	file->cd();
	MMPi0->Write();
	MPi0->Write();
	MMTPC->Write();
	file->Write();
	TCanvas* c9 = new TCanvas("c9","c9",1200,600);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	TGraphErrors* DifGr = new TGraphErrors(12,gr_ang,gr_ent,gr_ang_er,gr_err);
	DifGr->Draw("AP");
	DifGr->GetXaxis()->SetTitle("Cos(#theta_{CM})");
	DifGr->GetYaxis()->SetTitle("#Xi(1535)^{-} Count[/0.05] ");
	c9->SaveAs("AcptXiStar.png");
}
