#include "../../../KKManager.hh"
#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
double PI = acos(-1);
int evnum, runnum;
int evnumX, runnumX;
double mm1,mm2,MLd,MLdCor,tcm,MXi,DistLd,MomLd,MXiCor,mmpi0,MissingM,CheckSumM;
double MissingKpLd,MissingpKpLd;
bool FlgLd,FlgXi,InTargetLd,isGood; TFile* f1;TFile* f2;
TTree* tr1;
TTree* tr2;
vector<double>ldlen;
double range1=10,range2=250;
double tsize = 0.05;
void fcn (int& npar, double* grad, double& fval, double* par, int flag){
	double nll = 0;
	for(double len:ldlen){
		double expv = 1./par[1]*exp(-len/par[1]); 
		nll+=-log(expv);
	}
	fval = nll;
};

void DrawXi(){
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
	f1 = new TFile("SelectedEvents.root");
	tr1 = (TTree*)f1->Get("tree");
	tr1->SetBranchAddress("runnum",&runnumX);
	tr1->SetBranchAddress("evnum",&evnumX);
	tr1->SetBranchAddress("XiMM",&mm1);
	tr1->SetBranchAddress("XiThetaCM",&tcm);
//	f2 = new TFile("TPCInvM_cd15_Cor2ndWPrio.root");	
	f2 = new TFile("TPCInvM_cd15.root");	
	tr2 = (TTree*)f2->Get("tree");
	tr2->SetBranchAddress("runnum",&runnum);
	tr2->SetBranchAddress("evnum",&evnum);
	tr2->SetBranchAddress("MM",&mm2);
	tr2->SetBranchAddress("mmpi0",&mmpi0);
	tr2->SetBranchAddress("InvMLd",&MLd);
	tr2->SetBranchAddress("InvMXi",&MXi);
	tr2->SetBranchAddress("InvMXi",&MXiCor);
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
	tr2->SetBranchAddress("isGood",&isGood);
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
	TH2D* MMPi0 = new TH2D("MMPi0","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X.}",100,1.1,1.8,100,-1.5,1);
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
	for(int i = 0;i<ent;++i){
		tr2->GetEntry(i);
		bool mmfl = true;
		/*
			 for(int j=0;j<ent2; ++j){
			 tr1->GetEntry(j);
			 if(evnum == evnumX and runnum == runnumX){
			 mm2 = mm1; 
			 break;
			 }
			 }
			 */
		if(isGood){
			MMPi0->Fill(mm2,CheckSumM);
			MMTPC->Fill(CheckSumM);
			tree->Fill();
		}
		if(isGood and abs(CheckSumM-0.15)<0.15 and abs(mm2-1.55)<0.1){
			MPi0->Fill(CheckSumM);
		}
		//		if(abs(mm2-1.321)>0.1){
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
		if(FlgLd )LdHist->Fill(MLd);
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
		if(FlgXi)XiCorHist->Fill(MXiCor);
		if(FlgXi and !InTargetLd and flFlag)XiPropHist->Fill(MXi);
		if(FlgXi)LdXiHist->Fill(MLd);
		//		if(FlgLd)LdXiHist->Fill(MLd);
		if(FlgXi)
		{
			LdCorHist->Fill(MLdCor);
			LdXiPropHist->Fill(MLd);
		}
		if(FlgXi and !InTargetLd and flFlag){
			XiPropHist->Fill(MXiCor);
		}
		if(FlgXi and !InTargetLd and flFlag){
			ldlen.push_back(DistLd/MomLd*1.115);
			LdDist->Fill(DistLd/MomLd*1.115);
		}
	}
	int ldpeak = LdHist->GetMaximum();
	LambdaWBG->SetParLimits(0,0.7*ldpeak,ldpeak);
	LambdaWBG->SetParLimits(1,1.11,1.13);
	LambdaWBG->SetParLimits(2,0.001,0.01);
	LambdaWBG->SetParLimits(3,5,45);
	LambdaWBG->SetParLimits(4,1.10,1.15);
	LambdaWBG->SetParLimits(5,0.01,0.07);
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
	LdXiPropHist->Draw("same");
	LdXiPropHist->SetLineColor(kRed);
	LdHist->Fit("LambdaWBG","");
	LdHist->GetXaxis()->SetTitle("M_{p #pi^{-}} [GeV/c^{2}]");
	LdHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	LdHist->GetXaxis()->SetNdivisions(5,false);
	LdHist->GetYaxis()->SetNdivisions(10);
	LdHist->GetYaxis()->SetTitleSize(tsize);
	LdHist->GetXaxis()->SetTitleSize(tsize);
	LdHist->SetStats(0);
	TLatex* LdEnt = new TLatex(1.125,180,Form("Entries = %g",LdHist->GetEffectiveEntries()));
	LdEnt->Draw();
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
	XiCorHist->Draw();
	//	XiHist->Draw("same");
//	XiPropHist->Draw("same");
	//	XiPropHist->Draw("");
//	XiPropHist->SetLineColor(kRed);
	//	XiCorHist->SetLineColor(kBlack);
	XiCorHist->Fit("XiWBG","");
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
	fgau3->Draw("same");
	fgau4->Draw("same");
	if(drxi){
		double XiMass =XiWBG->GetParameter(1);
		double XiWidth =XiWBG->GetParameter(2);
		//		TText* XiLabel = new TText(1.33,20,Form(" %.0f +- %.0f MeV/c2",1000*XiMass,1000*XiWidth));
		//		XiLabel->Draw("same");
	}
	XiCorHist->SetStats(0);
	TLatex* XiEnt = new TLatex(1.33,0.9*(XiCorHist->GetMaximum()),Form("Entries = %g",XiCorHist->GetEffectiveEntries()));
	XiEnt->Draw();
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
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	MMPi0->Draw("colz");
	MMPi0->SetStats(0);
	MMPi0->GetYaxis()->SetTitle("MM[p(K^{-},K^{+}#Lambda #pi^{-})X][GeV/c^{2}]");
	MMPi0->GetXaxis()->SetTitle("MM{p(K^{-},K^{+})X}[GeV/c^{2}]");
	MMPi0->GetXaxis()->SetNdivisions(10);
	MMPi0->GetYaxis()->SetNdivisions(10);
	MMPi0->GetYaxis()->SetTitleSize(tsize*0.8);
	MMPi0->GetXaxis()->SetTitleSize(tsize);
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
	/*
	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");
*/
	TLatex* t = new TLatex(0.,1.65,"#pi^{0}");
	t->SetTextColor(kRed);
	t->SetTextSize(0.2);
//	t->Draw();
	c6->SaveAs("MMpKmKpLdPi_12C.png");

	TCanvas* c7 = new TCanvas("c7","c7",900,900);
	MPi0->SetStats(0);
	TLatex* stats = new TLatex(0.2,7,Form("Entries : %g", MPi0->GetEffectiveEntries()));
	MPi0->Draw();
//	stats->Draw();
	MPi0->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}}-X[GeV/c]");
	TString titleY = Form(" Entries /%g MeV/c{^}",MPi0->GetXaxis()->GetBinWidth(5)*1000);
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
	MMTPC->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}}-X[GeV/c]");
	TString titleY2 = Form(" Entries /%g MeV/c{^}",MMTPC->GetXaxis()->GetBinWidth(5)*1000);
	MMTPC->SetLineWidth(3);
	MMTPC->SetLineColor(kBlack);
	MMTPC->GetYaxis()->SetTitle(titleY2);
	MMTPC->GetYaxis()->SetNdivisions(5);
	MMTPC->GetXaxis()->SetNdivisions(10);
	c8->SaveAs("MMTPC.png");
	TCanvas* c9 = new TCanvas("c9","c9",900,900);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	MMKmKpLd->Draw();
	MMKmKpLd->SetStats(0);
	MMKmKpLd->SetLineWidth(3);
	MMKmKpLd->SetLineColor(kBlack);
	MMKmKpLd->GetXaxis()->SetNdivisions(10);
	MMKmKpLd->GetXaxis()->SetTitle("MM[X(K^{-},K^{+})#Lambda]");
	MMKmKpLd->GetYaxis()->SetNdivisions(5);
	MMKmKpLd->GetYaxis()->SetTitle("Entries / 10 MeV");
	MMKmKpLd->GetYaxis()->SetTitleSize(tsize);
	c9->SaveAs("MMKmKpLd_12C.png");
	TCanvas* c10 = new TCanvas("c10","c10",900,900);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	MMKmKpLd2D->Draw("colz");
	MMKmKpLd2D->SetStats(0);
	MMKmKpLd2D->GetYaxis()->SetNdivisions(10);
	MMKmKpLd2D->GetYaxis()->SetTitle("MM[X(K^{-},K^{+})#Lambda]");
	MMKmKpLd2D->GetYaxis()->SetTitleSize(tsize);
	MMKmKpLd2D->GetXaxis()->SetNdivisions(10);
	MMKmKpLd2D->GetXaxis()->SetTitle("MM[p(K^{-},K^{+})X]");
	c10->SaveAs("MMKmKpLd2D_12C.png");
	TCanvas* c11 = new TCanvas("c11","c11",900,900);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	MMKmpKpLd->Draw();
	MMKmpKpLd->SetStats(0);
	MMKmpKpLd->SetLineWidth(3);
	MMKmpKpLd->SetLineColor(kBlack);
	MMKmpKpLd->GetXaxis()->SetNdivisions(10);
	MMKmpKpLd->GetXaxis()->SetTitle("MM[p(K^{-},K^{+}#Lambda) X]");
	MMKmpKpLd->GetYaxis()->SetNdivisions(5);
	MMKmpKpLd->GetYaxis()->SetTitle("Entries / 10 MeV");
	MMKmpKpLd->GetYaxis()->SetTitleSize(tsize);
	c11->SaveAs("MMpKmKpLd_12C.png");
	TCanvas* c12 = new TCanvas("c12","c12",900,900);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	MMKmpKpLd2D->Draw("colz");
	MMKmpKpLd2D->SetStats(0);
	MMKmpKpLd2D->GetYaxis()->SetNdivisions(10);
	MMKmpKpLd2D->GetYaxis()->SetTitle("MM[p(K^{-},K^{+}#Lambda) X]");
	MMKmpKpLd->GetYaxis()->SetNdivisions(5);
	MMKmpKpLd2D->GetYaxis()->SetTitleSize(tsize);
	MMKmpKpLd2D->GetXaxis()->SetNdivisions(10);
	MMKmpKpLd2D->GetXaxis()->SetTitle("MM[p(K^{-},K^{+})X]");
	c12->SaveAs("MMpKmKpLd2D_12C.png");
	file->cd();
	MMPi0->Write();
	MPi0->Write();
	MMTPC->Write();
	file->Write();
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








