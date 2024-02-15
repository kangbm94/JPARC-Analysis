#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
double PI = acos(-1);
TFile* f1;TFile* f2;
TTree* tr1;
TTree* tr2;
vector<double>ldlen;
double range1=10,range2=250;
double tsize = 0.05;
int TAPS = 20,TB=15;
double kpPx,kpPy,kpPz;
double kmPx,kmPy,kmPz;
double MomxXi,MomyXi,MomzXi;
double MomxLd,MomyLd,MomzLd;
double MomxP,MomyP,MomzP;
double dM,dP;
bool FlgXi;
double InvMXi,InvMLd;
double cTh,cPh,cPh1,cPh2,cPh3;
double MissMass,Coplanarity;
double mPi0 = 0.135;
double m2Pi0 = 0.135*0.135;
double m2Pi0Sig = 0.013;
bool TrigAPS,TrigB;
void fcn (int& npar, double* grad, double& fval, double* par, int flag){
	double nll = 0;
	for(double len:ldlen){
		double expv = 1./par[1]*exp(-len/par[1]); 
		nll+=-log(expv);
	}
	fval = nll;
};
double Gaussian(double x, double mean, double sigma, double peak){
	double par=(x-mean)/sigma;
	double val = peak*exp(-par*par/2);
	return val;
}
double Gaussianf(double* x,double* par){
	double amp1 = par[0],mean1=par[1],width1=par[2];
	return Gaussian(x[0],mean1,width1,amp1);
}
double GausWithBG(double* x, double* par){
	double amp1 = par[0],mean1=par[1],width1=par[2];
	double amp2 = par[3],mean2=par[4],width2=par[5];
	double signal = Gaussian(x[0],mean1,width1,amp1);
	double bg = Gaussian(x[0],mean2,width2,amp2);
	return signal + bg;
}
double MMCut = 0.015;
double XiStarMMCut = 0.0365;//sigma = 18.25 MeV
double dM2Cut = 0.012;
double CopCut = 1.10;

void DrawXiGenfit(){
	SetStyle();
//	gStyle->SetOptStat(111);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.4);
	gStyle->SetStatH(0.5);
	ldlen.clear();
	TH2D* MMScat = new TH2D("KuramaMM vs (Kurama-Xi)MM","KuramaMM vs (Kurama-Xi)MM",20,1.2,1.7,20,-0.5,0.5);
	TF1* LambdaWBG = new TF1("LambdaWBG","GausWithBG",1.08,1.15,6);
	TF1* XiWBG = new TF1("XiWBG","GausWithBG",1.28,1.36,6);
	LambdaWBG->SetParNames("LdConst","LdMass","LdWidth","LdBgConst","LdBgMass","LdBgWidth");
	XiWBG->SetParNames("XiConst","XiMass","XiWidth","XiBgConst","XiBgMass","XiBgWidth");
	double par[6];
	bool CH2Target = true;
	gStyle->SetOptFit(000);
	gStyle->SetOptStat(10);
	gStyle->SetTitleFontSize(0.1);
	gStyle->SetTitleFont(132,"x");//,“t”);
	gStyle->SetTitleFont(132,"t");//,“t”);
	gStyle->SetTitleFont(132,"y");//,“t”);
	if(CH2Target){
		f2 = new TFile("Ch2PolaAnal.root");	
	}
	else{
		f2 = new TFile("CarbonPolaAnal.root");	
	}
	TTree* tree = (TTree*)f2->Get("tree");
	tree->SetBranchAddress("TrigAPS",&TrigAPS);
	tree->SetBranchAddress("TrigB",&TrigB);
	tree->SetBranchAddress("MissMass",&MissMass);
	tree->SetBranchAddress("dM",&dM);
	tree->SetBranchAddress("Coplanarity",&Coplanarity);
	tree->SetBranchAddress("cTh",&cTh);
	tree->SetBranchAddress("cPh",&cPh);
	tree->SetBranchAddress("cPh1",&cPh1);
	tree->SetBranchAddress("cPh2",&cPh2);
	tree->SetBranchAddress("cPh3",&cPh3);
	tree->SetBranchAddress("InvMXi",&InvMXi);
	tree->SetBranchAddress("InvMLd",&InvMLd);

	tree->SetBranchAddress("PKm_x",&kmPx);
	tree->SetBranchAddress("PKm_y",&kmPy);
	tree->SetBranchAddress("PKm_z",&kmPz);
	
	tree->SetBranchAddress("PKp_x",&kpPx);
	tree->SetBranchAddress("PKp_y",&kpPy);
	tree->SetBranchAddress("PKp_z",&kpPz);
	
	tree->SetBranchAddress("PXi_x",&MomxXi);
	tree->SetBranchAddress("PXi_y",&MomyXi);
	tree->SetBranchAddress("PXi_z",&MomzXi);
	
	tree->SetBranchAddress("PLd_x",&MomxLd);
	tree->SetBranchAddress("PLd_y",&MomyLd);
	tree->SetBranchAddress("PLd_z",&MomzLd);
	
	tree->SetBranchAddress("PP_x",&MomxP);
	tree->SetBranchAddress("PP_y",&MomyP);
	tree->SetBranchAddress("PP_z",&MomzP);
	TH1D* LdHist = new TH1D("Lambda","Invariant Mass(p #pi)",50,1.07,1.17);
	TH1D* LdHistMMCut = new TH1D("Lambda","Invariant Mass(p #pi)MMCut",50,1.07,1.17);
	TH1D* LdHistdM2Cut = new TH1D("Lambda","Invariant Mass(p #pi)2M2Cut",50,1.07,1.17);
	TH1D* LdXiHist = new TH1D("LambdawXi","#Lambda w #Xi Invariant Mass",50,1.07,1.17);
	TH1D* LdXiPropHist = new TH1D("LambdaPropwXi","#Lambda w #Xi Invariant Mass Prop",50,1.07,1.17);
	TH1D* XiHist = new TH1D("Xi","#Xi Invariant Mass",100,1.25,1.40);
	TH1D* XiHistMMCut = new TH1D("Xi","#Xi Invariant MassMMCut",100,1.25,1.40);
	TH1D* XiHistdM2Cut = new TH1D("Xi","#Xi Invariant MassdM2Cut",100,1.25,1.40);
	TH1D* LdDist = new TH1D("LdPropLength","#Lambda FlightLength/#gamma#beta",20,range1,range2);
	TGraph* LdFlGr;
	TH1D* LdMMHist = new TH1D("LdMissing Mass","#Lambda tagged Missing Mass",100,1,2);
	TH1D* MMHist = new TH1D("Missing Mass","Missing Mass",100,1,2);
	TH2D* MMPi0 = new TH2D("MMPi0","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X.}",1000,-1.5,1,100,1.1,1.8);
	TH2D* MMPi0Cor = new TH2D("MMPi0Cor","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X.}",100,-0.1,0.1,100,1.1,1.8);
	TH2D* MMPi0Vtx = new TH2D("MMPi0Vtx","Sum_{TPC}-X:MM_{p(K^{-},K^{+})X.}",100,-1.5,1,100,1.1,1.8);
	TH1D* MMXiStar = new TH1D("MMXiStar","MMXiStar",30,1.53-1.1*XiStarMMCut,1.53+1.1*XiStarMMCut);
	TH1D* MMXiStarTagged = new TH1D("MMXiStarTagged","MMXiStarTagged",30,1.53-1.1*XiStarMMCut,1.53+1.1*XiStarMMCut);
	TH1D* MPi0 = new TH1D("MPi0","Sum_{TPC}-X",30,m2Pi0 - 2*m2Pi0Sig,m2Pi0+2*m2Pi0Sig);
	TH1D* MMTPC = new TH1D("MMTPC","Sum_{TPC}-X",2000,-0.3,0.3);
	TH2D* MMKmKpLd2D = new TH2D("MissingMass(Km,KpLd):KuramaMM","MissingMass(Km,KpLd):KuramaMM",100,1.1,1.8,100,-1,2);
	TH1D* MMKmKpLd = new TH1D("MissingMass(Km,KpLd)","MissingMass(Km,KpLd)",100,-1,2);
	TH2D* MMKmpKpLd2D = new TH2D("MissingMass(Kmp,KpLd):KuramaMM","MissingMass(Kmp,KpLd):KuramaMM",100,1.1,1.8,100,-1,2);
	TH1D* MMKmpKpLd = new TH1D("MissingMass(Kmp,KpLd)","MissingMass(Kmp,KpLd)",100,-1,2);
	int ent = tree->GetEntries();
	bool acpt;
	int n_ximm=0;
	for(int i = 0;i<ent;++i){
		tree->GetEntry(i);
		bool mmfl = true;
		bool acpt = true;
		bool MMCutFlag = true;
		bool dM2CutFlag = true;
		double dM2 = dM*abs(dM);
		if(1){
			if(abs(InvMLd-1.115)>0.03) acpt = false;
			if(abs(InvMXi-1.321)>0.03) acpt = false;
		}
		if(acpt ){
			MMPi0->Fill(dM2,MissMass);
			if(1){
				MMPi0Cor->Fill(dM2,MissMass);
			}
			else{
				MMPi0Vtx->Fill(dM2,MissMass);
			}
			if(abs(MissMass-1.321)<0.1)MMTPC->Fill(dM2);
		}
		if(abs(dM2-m2Pi0)<2*m2Pi0Sig and dM2>0 and abs(MissMass-1.53)<XiStarMMCut){
			MPi0->Fill(dM2);
		}
		if(abs(MissMass-1.321)>0.1){
			mmfl = false;
		}
		if(abs(MissMass-1.321)>MMCut){
			MMCutFlag = false;
		}
		if(abs(dM2)>dM2Cut){
			dM2CutFlag = false;
		}
		if(abs(MissMass - 1.53) < XiStarMMCut){
			MMXiStar->Fill(MissMass);
			if( dM2 > 0 and abs(dM2-m2Pi0) < 2*m2Pi0Sig){
				MMXiStarTagged->Fill(MissMass);
			}
		}
		//		if(runnum != 5641) continue;
		n_ximm++;
		bool flFlag = false;
		if(mmfl){
			LdHist->Fill(InvMLd);
			XiHist->Fill(InvMXi);
			if(MMCutFlag){
				LdHistMMCut->Fill(InvMLd);
				XiHistMMCut->Fill(InvMXi);
				if(dM2CutFlag){
					LdHistdM2Cut->Fill(InvMLd);
					XiHistdM2Cut->Fill(InvMXi);
				}
			}
		}
		



	}
	int ldpeak = LdHist->GetMaximum();
	LambdaWBG->SetParLimits(0,0.6*ldpeak,ldpeak);
	LambdaWBG->SetRange(1.08,1.14);
	LambdaWBG->SetParLimits(1,1.11,1.13);
	LambdaWBG->SetParLimits(2,0.001,0.01);
	LambdaWBG->SetParLimits(3,0,ldpeak* 0.4);
	LambdaWBG->SetParLimits(4,1.11,1.13);
	LambdaWBG->SetParLimits(5,0.005,0.03);
	int xipeak = XiHist->GetMaximum();
	XiWBG->SetRange(1.28,1.36);
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
	c1->Divide(2,1);
	TF1* fgau = new TF1("fgau","Gaussianf",1.07,1.4,3);
	TF1* fgau2 = new TF1("fgau2","Gaussianf",1.07,1.4,3);
	TF1* fgau3 = new TF1("fgau3","Gaussianf",1.07,1.4,3);
	TF1* fgau4 = new TF1("fgau4","Gaussianf",1.07,1.4,3);
	bool drld=true,drxi=true,drprop=true;
	c1->cd(1);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	LdHist->Draw();
	LdHistMMCut->Draw("same");
	LdHistdM2Cut->Draw("same");
	//	LdXiHist->SetLineColor(kRed);
	//	LdCorHist->Draw(""); //	LdCorHist->SetLineColor(kBlack); //	LdXiHist->Draw("");
//	LdXiPropHist->Draw("same");
//	LdXiPropHist->SetLineColor(kRed);
	LdHist->GetXaxis()->SetTitle("M_{p #pi^{-}} [GeV/c^{2}]");
	LdHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	LdHist->GetXaxis()->SetNdivisions(5,false);
	LdHist->GetYaxis()->SetNdivisions(10);
	LdHist->GetYaxis()->SetTitleSize(tsize);
	LdHist->GetXaxis()->SetTitleSize(tsize);
	LdHist->SetStats(0);
	LdHist->Fit("LambdaWBG","R0");
	cout<<"Gaus"<<endl;
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
	fgau2->SetLineStyle(kDashed);
//	fgau->Draw("same");
//	fgau2->Draw("same");
	double LdMass =LambdaWBG->GetParameter(1);
	double LdWidth =LambdaWBG->GetParameter(2);
	//		TText* LdLabel = new TText(1.12,40,Form(" %.0f \n+- %.0f MeV/c2",1000*LdMass,1000*LdWidth));
	//		LdLabel->Draw("same");
	//	c2->cd();
	c1->cd(2);
	XiHist->Draw();
	XiHistMMCut->Draw("same");
	XiHistdM2Cut->Draw("same");

	LdHist->SetLineColor(kBlack);
	LdHistMMCut->SetLineColor(kBlue);
	LdHistdM2Cut->SetLineColor(kRed);
	XiHist->SetLineColor(kBlack);
	XiHistMMCut->SetLineColor(kBlue);
	XiHistdM2Cut->SetLineColor(kRed);

	TLegend* legXi = new TLegend(0.55,0.7,0.89,0.89);
	legXi->AddEntry(XiHist,Form("Reconstructed Events: %g",XiHist->GetEntries()),"l");
	legXi->AddEntry(XiHistMMCut,Form("p(K^{-},K^{+})X Cut Events: %g",XiHistMMCut->GetEntries()),"l");
	legXi->AddEntry(XiHistdM2Cut,Form("p(K^{-},K^{+}#Xi)X Cut Events: %g",XiHistdM2Cut->GetEntries()),"l");
	legXi->Draw();
	legXi->SetBorderSize(0);
	//	XiHist->Draw("same");
//	XiPropHist->Draw("same");
	//	XiPropHist->Draw("");
//	XiPropHist->SetLineColor(kRed);
	//	XiHist->SetLineColor(kBlack);
	XiHist->Fit("XiWBG","R0");
	XiHist->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}} [GeV/c^{2}]");
	XiHist->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}} [GeV/c^{2}]");
	XiHist->GetYaxis()->SetTitle("Entries / 2 MeV/c^{2}");
	XiHist->GetXaxis()->SetNdivisions(5,false);
	XiHist->GetXaxis()->SetTitleSize(tsize);;
	XiHist->GetYaxis()->SetNdivisions(10);
	XiHist->GetYaxis()->SetTitleSize(tsize);;
	for(int i=0;i<6;++i){
		par[i]= XiWBG->GetParameter(i);
	}
	for(int i=0;i<3;++i){
		fgau3->SetParameter(i,par[i]);
		fgau4->SetParameter(i,par[i+3]);
	}
	fgau3->SetLineColor(kGreen);
	fgau4->SetLineColor(kBlack);
	fgau4->SetLineStyle(kDashed);
//	fgau3->Draw("same");
//	fgau4->Draw("same");
	if(drxi){
		double XiMass =XiWBG->GetParameter(1);
		double XiWidth =XiWBG->GetParameter(2);
		//		TText* XiLabel = new TText(1.33,20,Form(" %.0f +- %.0f MeV/c2",1000*XiMass,1000*XiWidth));
		//		XiLabel->Draw("same");
	}
	XiHist->SetStats(0);
	TLatex* XiEnt = new TLatex(1.33,0.9*(XiHist->GetMaximum()),Form("Entries = %g",XiHist->GetEffectiveEntries()));
//	XiEnt->Draw();
	cout<<"LdEnt: "<<LdHist->GetEffectiveEntries()<<endl;
	cout<<"XiEnt: "<<XiHist->GetEffectiveEntries()<<endl;

	for(int i = 0;i<ent;++i){
		tree->GetEntry(i);
		if((InvMLd - LdMass)<3*LdWidth)LdMMHist ->Fill(MissMass);
		MMHist ->Fill(MissMass);
	}
	c1->SaveAs("XiReconHist.png");
	TCanvas* c6 = new TCanvas("c6","c6",900,900);
	gPad->SetMargin(0.15,0.1,0.1,0.1);
	MMPi0Cor->Draw("colz");
	MMPi0Cor->SetStats(0);
	MMPi0Cor->GetXaxis()->SetTitle("MM^{2}(p(K^{-},K^{+}#Xi)[GeV^{2}/c^{4}]");
	MMPi0Cor->GetYaxis()->SetTitle("MM(p(K^{-},K^{+})X)[GeV/c^{2}]");
	MMPi0Cor->GetXaxis()->SetNdivisions(5);
	MMPi0Cor->GetYaxis()->SetNdivisions(10);
	MMPi0Cor->GetYaxis()->SetTitleSize(tsize*0.8);
	MMPi0Cor->GetXaxis()->SetTitleSize(tsize);
	TLine* l1 = new TLine(0,1.53-XiStarMMCut,0,1.53+XiStarMMCut);
	TLine* l2 = new TLine(0,1.53+XiStarMMCut,m2Pi0+2*m2Pi0Sig,1.53+XiStarMMCut);
	TLine* l3 = new TLine(m2Pi0 + 2*m2Pi0Sig,1.53-XiStarMMCut,m2Pi0+2*m2Pi0Sig,1.53+XiStarMMCut);
	TLine* l4 = new TLine(0,1.53-XiStarMMCut,m2Pi0+2*m2Pi0Sig,1.53-XiStarMMCut);
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

	TLine* l5 = new TLine(0.1,1.321-MMCut,0.1,1.321+MMCut);
	TLine* l6 = new TLine(-0.1,1.321+MMCut,0.1,1.321+MMCut);
	TLine* l7 = new TLine(-0.1,1.321-MMCut,-0.1,1.321+MMCut);
	TLine* l8 = new TLine(-0.1,1.321-MMCut,0.1,1.321-MMCut);
	l5->SetLineColor(kBlack);
	l5->SetLineStyle(kDashed);
	l5->SetLineWidth(4);
	l6->SetLineColor(kBlack);
	l6->SetLineStyle(kDashed);
	l6->SetLineWidth(4);
	l7->SetLineColor(kBlack);
	l7->SetLineStyle(kDashed);
	l7->SetLineWidth(4);
	l8->SetLineColor(kBlack);
	l8->SetLineStyle(kDashed);
	l8->SetLineWidth(4);
//	l5->Draw("same");
	l6->Draw("same");
//	l7->Draw("same");
	l8->Draw("same");

	TLatex* t = new TLatex(0.,1.65,"#pi^{0}");
	t->SetTextColor(kRed);
	t->SetTextSize(0.2);
//	t->Draw();
	if(CH2Target)c6->SaveAs("MMPi0.png");
	else c6->SaveAs("MMPi0Dia.png");
	TCanvas* c7 = new TCanvas("c7","c7",900,900);
//	MPi0->SetStats(0);
//	TLatex* stats = new TLatex(0.2,7,Form("Entries : %g", MPi0->GetEffectiveEntries()));
	MPi0->Draw();
	MPi0->Fit("gaus");
//	stats->Draw();
	MPi0->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}}-X[GeV/c]");
	TString titleY = Form(" Entries /%.2g MeV/c^{2}",MPi0->GetXaxis()->GetBinWidth(5)*1000);
	MPi0->SetLineWidth(3);
	MPi0->SetLineColor(kBlack);
	MPi0->GetYaxis()->SetTitle(titleY);
	MPi0->GetYaxis()->SetTitleSize(tsize);
	MPi0->GetYaxis()->SetNdivisions(5);
	MPi0->GetXaxis()->SetNdivisions(10);
	MPi0->GetXaxis()->SetTitle("Mass^{2}[GeV^{2}/c^{4}]");
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
	c4->cd();
	MMXiStar->Draw();
	MMXiStar->SetStats(0);
	MMXiStar->GetYaxis()->SetNdivisions(5);
	MMXiStar->GetXaxis()->SetNdivisions(10);
	MMXiStarTagged->SetLineColor(kRed);
	MMXiStarTagged->Draw("same");
	TLegend* legXiStar = new TLegend(0.55,0.7,0.89,0.89);
	legXiStar->AddEntry(MMXiStar,Form("#Xi Reconstructed Events: %g",MMXiStar->GetEntries()),"l");
	legXiStar->AddEntry(MMXiStarTagged,Form("#pi^{0} Tagged Events: %g",MMXiStarTagged->GetEntries()),"l");
	legXiStar->Draw();
	legXiStar->SetBorderSize(0);
	c4->SaveAs("XiStar.png");
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








