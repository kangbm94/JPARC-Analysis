#include "PolarizationAnal.hh"
#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
TH1D* HTh;
TH1D* HPh;
TH1D* HPh1;
TH1D* HPh2;
TH1D* HPh3;
TH1D* HCM;
TH1D* HCMWe;
vector<double>*DecayMomx;
vector<double>*DecayMomy;
vector<double>*DecayMomz;
vector<double> *MM;
vector<double>*PKp;
vector<double>*uKp;
vector<double>*vKp;
vector<double>*PKm;
vector<double>*uKm;
vector<double>*vKm;
vector<int>* trigflag;
int TAPS = 20,TB=15;
double kpPx,kpPy,kpPz;
double kmPx,kmPy,kmPz;
double MomxXi,MomyXi,MomzXi;
double MomxLd,MomyLd,MomzLd;
double MomxP,MomyP,MomzP;
bool FlgXi;
double InvMXi,InvMLd;
double cTh,cPh,cPh1,cPh2,cPh3;
double MissMass,Coplanarity,dM;
bool TrigAPS,TrigB;
void RealPola(int i);
TTree* tree2;
vector<double>datap;
void FCN(int &npar, double* gin, double &f, double * par, int iflag){
	f=0;
	double par1 = par[0];
	for(auto dp:datap){
		double delta = 1.0 - dp * par1;
		f+= delta* delta;
	}
//	cout<<Form("ND = %lu,chi2 = %g, par = %g",datap.size(),f,par1)<<endl;
}
void PolaSummary(){
	bool UnbinnedFit = false;
	TFile* Acpt = new TFile("KpAcceptance.root");	
	TH2D* hAcpt = (TH2D*)Acpt->Get("BeamAcceptance");
	int nbin = 4;
	int nbinCM = 50;
	double MMCut = 0.055;
	double dM2Cut = 0.912;
	MMCut *=1;
	dM2Cut*=1;
	double CopCut = 1.10;
	SetStyle();
	gStyle->SetTitleSize(0.08,"XY");
	HTh = new TH1D("HistTh","HistTh",nbin,-1.,1.);
	HPh = new TH1D("HistPh","HistPh",nbin,-1.,1.);
	HPh1 = new TH1D("HistPh1","HistPh1",nbin,-1.,1.);
	HPh2 = new TH1D("HistPh2","HistPh2",nbin,-1.,1.);
	HPh3 = new TH1D("HistPh3","HistPh3",nbin,-1.,1.);
	HCM = new TH1D("HistCMTheta","HistCMTheta",100,0.7,1.);
	HCMWe = new TH1D("HistCMThetaWeigh","HistCMThetaWeigh",100,0.7,1.);
	TF1* fLin = new TF1("fLin","[0]+[1]*x",-1,1);
	TF1* fPLin = new TF1("fPLin","0.5+[0]*x",-1,1);
	
	TH1D* HMM = new TH1D("HistMM","HistMM",100,1.2,1.7);
	TH1D* HMMCut = new TH1D("HistMMCut","HistMMCut",100,1.2,1.7);
	int rn_s = 5641;
	int rn_e = 5666;
//	rn_s = 5667;
//	rn_e = 5697;
	
///	TFile* file = new TFile(Form("Ch2PolaAnal.root"));
	TFile* file = new TFile(Form("CarbonPolaAnal.root"));
//	TFile* file = new TFile(Form("DstPolaAnal_%d_%d.root",rn_s,rn_e));
	TTree* tree = (TTree*)file->Get("tree");
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

	vector<double>VecTh;
	vector<double>VecPh;
	vector<double>VecPh1;
	vector<double>VecPh2;
	vector<double>VecPh3;
	int ent = tree->GetEntries();
	double WSum=0;
	int NDat = 0;
	double PolAngle[100]={0};
	double PolAngleE[100]={0};
	double PolTh[100]={0};
	double PolPh1[100]={0};
	double PolPh2[100]={0};
	double PolPh3[100]={0};
	double PolThE[100]={0};
	double PolPh1E[100]={0};
	double PolPh2E[100]={0};
	double PolPh3E[100]={0};

	double CMAngle[1000]={0};
	double CMAngleE[1000]={0};
	double CMTh[1000]={0};
	double CMThE[1000]={0};
	double lbCM=0.8,hbCM=1;
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		TVector3 TVKP(kpPx,kpPy,kpPz);
		TVector3 TVKM(kmPx,kmPy,kmPz);
		double kpP = TVKP.Mag(); 
		if(kpP == 0)continue;
		double kmP = TVKM.Mag(); 
		auto CosKP = kpPz/kpP; 
		double KpTheta = acos(CosKP)*180./acos(-1);
		int BinAcpt = hAcpt-> FindBin(KpTheta,kpP);
		double content = hAcpt->GetBinContent(BinAcpt);
		double w = 1./content;
		if(!isnan(cTh) and abs(Coplanarity)< CopCut	and abs(MissMass - 1.321)<MMCut and dM*dM<dM2Cut){
//			cout<<Form("PK = %g, th = %g, w = %g",kpP,KpTheta,w)<<endl;
			HMMCut->Fill(MissMass);
			VecTh.push_back(cTh);
			VecPh.push_back(cPh);
			VecPh1.push_back(cPh1);
			VecPh2.push_back(cPh2);
			VecPh3.push_back(cPh3);
			double mk = 0.493;
			double mp = 0.938;
			TLorentzVector LVKP(TVKP,hypot(kpP,mk));
			TLorentzVector LVKM(TVKM,hypot(kmP,mk));
			TLorentzVector LVTarget(0,0,0,mp);
			auto LVCM = LVKM +LVTarget;
			auto V = LVCM.BoostVector();
			auto KMLab = LVKM.Vect();
			auto KPLab = LVKP.Vect();
			double KMPLab = KMLab.Mag();
			double KPPLab = KPLab.Mag();
			double LabTh = (KMLab*KPLab)/KMPLab/KPPLab;	
			LVKM.Boost(-V);
			LVKP.Boost(-V);
			auto KMCM = LVKM.Vect();
			auto KPCM = LVKP.Vect();
			double KMPCM = KMCM.Mag();
			double KPPCM = KPCM.Mag();
			double CMThe = (KMCM*KPCM)/KMPCM/KPPCM;	
//			cout<<"Lab "<<LabTh<<"CM "<<CMThe<<endl;
			WeightedFill(HCM,CMThe,1);	
			if(w >20)continue;
				//continue;
			int BinCM = GetBin(nbinCM,0.7,1,CMThe);
			CMTh[BinCM]+= w;
			CMThE[BinCM]+= w*w;
			WSum+=w;
			NDat++;
			WeightedFill(HTh,cTh,w);	
			WeightedFill(HPh,cPh,w);	
			WeightedFill(HPh1,cPh1,w);	
			WeightedFill(HPh2,cPh2,w);	
			WeightedFill(HPh3,cPh3,w);	
			int BinTh = GetBin(nbin,-1,1,cTh);
			PolTh[BinTh]+= w;
			PolThE[BinTh]+= w*w;
			int BinPh1 = GetBin(nbin,-1,1,cPh1);
			PolPh1[BinPh1]+= w;
			PolPh1E[BinPh1]+= w*w;
			int BinPh2 = GetBin(nbin,-1,1,cPh2);
			PolPh2[BinPh2]+= w;
			PolPh2E[BinPh2]+= w*w;
			int BinPh3 = GetBin(nbin,-1,1,cPh3);
			PolPh3[BinPh3]+= w;
			PolPh3E[BinPh3]+= w*w;
			WeightedFill(HMM,MissMass,w);
		}
	}
	for(int ib=0;ib<nbin;++ib){
		PolAngle[ib]=GetBinCenter(nbin,-1,1,ib);	
		PolAngleE[ib]= 2./nbin/sqrt(12);
		PolThE[ib]= sqrt(PolThE[ib]);
		PolPh1E[ib]= sqrt(PolPh1E[ib]);
		PolPh2E[ib]= sqrt(PolPh2E[ib]);
		PolPh3E[ib]= sqrt(PolPh3E[ib]);
	}
	for(int ib=0;ib<nbinCM;++ib){
		CMAngle[ib] = GetBinCenter(nbinCM,lbCM,hbCM,ib);
		CMAngleE[ib]= (hbCM-lbCM)/nbinCM/sqrt(12);
		CMThE[ib] = sqrt(CMThE[ib]);	
	}

	TString YTitle = Form("Counts / %.2g ",2./nbin);
	TString YTitleCM = Form("Counts / %.4g ",(hbCM-lbCM)/nbinCM);
	TGraphErrors* gTh = new TGraphErrors(nbin,PolAngle,PolTh,PolAngleE,PolThE);

	TGraphErrors* gPh1 = new TGraphErrors(nbin,PolAngle,PolPh1,PolAngleE,PolPh1E);
	TGraphErrors* gPh2 = new TGraphErrors(nbin,PolAngle,PolPh2,PolAngleE,PolPh2E);
	TGraphErrors* gPh3 = new TGraphErrors(nbin,PolAngle,PolPh3,PolAngleE,PolPh3E);
	TGraphErrors* gCM = new TGraphErrors(nbinCM,CMAngle,CMTh,CMAngleE,CMThE);
	TMinuit* min = new TMinuit(1);
	min->SetPrintLevel(-1);
	int ierflg = 0;
	double arglist[2]={500,0.01};
	double par1,par0,bins,par1E,par0E;
//	double par1,par1E;
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(2,2);
	c1->cd(1);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	gTh->Draw("AP");
	gTh->GetXaxis()->SetTitle("cos#theta");
	gTh->GetYaxis()->SetTitle(YTitle);
	bins = WSum / nbin;
	cout<<NDat<<endl;
	cout<<WSum<<endl;
	fLin->FixParameter(0,bins);
	if(UnbinnedFit){
		datap = VecTh;
		min->SetFCN(FCN);	
		min->mnparm(0,"p1",0.,0.01,-1,1,ierflg);
		min->mnexcm("MIGRAD",arglist,2,ierflg);
		min->GetParameter(0,par1,par1E);
		par0=1;
		fLin->SetParameter(0,bins);
		fLin->SetParameter(1,par1*bins);
		auto fLin1 = (TF1*)fLin->Clone("fLin1");
		fLin1->Draw("same");
	}
	else{
		gTh->Fit("fLin");
		par0 = fLin->GetParameter(0)/bins;
		par1 = fLin->GetParameter(1)/bins;
	}	
	cout<<"N = "<<par0<<endl;
	cout<<"P = "<<par1<<endl;
	c1->cd(2);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	gPh1->Draw("AP");
	gPh1->GetXaxis()->SetTitle("cos#phi_{#beta}");
	gPh1->GetYaxis()->SetTitle(YTitle);
	if(UnbinnedFit){
		datap = VecPh1;
		min->mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
		min->mnexcm("MIGRAD",arglist,2,ierflg);
		min->GetParameter(0,par1,par1E);
		fLin->SetParameter(0,bins);
		fLin->SetParameter(1,par1*bins);
		cout<<"P = "<<par1<<endl;
		auto fLin2 = (TF1*)fLin->Clone("fLin2");
		fLin2->Draw("same");
	}
	else{
		gPh1->Fit("fLin");
		par0 = fLin->GetParameter(0)/bins;
		par1 = fLin->GetParameter(1)/bins;
		cout<<"N = "<<par0<<endl;
		cout<<"P = "<<par1<<endl;
	}

	c1->cd(3);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	gPh2->Draw("AP");
	gPh2->GetXaxis()->SetTitle("cos#phi_{#gamma}");
	gPh2->GetYaxis()->SetTitle(YTitle);
	if(UnbinnedFit){
		datap = VecPh2;
		min->mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
		min->mnexcm("MIGRAD",arglist,2,ierflg);
		min->GetParameter(0,par1,par1E);
		fLin->SetParameter(0,bins);
		fLin->SetParameter(1,par1*bins);
		cout<<"P = "<<par1<<endl;
		auto fLin3 = (TF1*)fLin->Clone("fLin3");
		fLin3->Draw("same");
	}
	else{
		gPh2->Fit("fLin");
		par0 = fLin->GetParameter(0)/bins;
		par1 = fLin->GetParameter(1)/bins;
		cout<<"N = "<<par0<<endl;
		cout<<"P = "<<par1<<endl;
	}
	c1->cd(4);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	gPh3->Draw("AP");
	gPh3->GetXaxis()->SetTitle("cos#phi_{#alpha}");
	gPh3->GetYaxis()->SetTitle(YTitle);
	if(UnbinnedFit){
		datap = VecPh3;
		min->mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
		min->mnexcm("MIGRAD",arglist,2,ierflg);
		min->GetParameter(0,par1,par1E);
		fLin->SetParameter(0,bins);
		fLin->SetParameter(1,par1*bins);
		cout<<"P = "<<par1<<endl;
		auto fLin4 = (TF1*)fLin->Clone("fLin4");
		fLin4->Draw("same");
	}
	else{
		gPh3->Fit("fLin");
		par0 = fLin->GetParameter(0)/bins;
		par1 = fLin->GetParameter(1)/bins;
		double par1Er = fLin->GetParError(1)/bins;
		cout<<"N = "<<par0<<endl;
		cout<<"P = "<<par1<<" +- "<<par1Er<<endl;
	}
	c1->SaveAs("Polarization.png");
	/*
	datap = VecPh3;
	min->SetFCN(FCN);	
	min->mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
	min->mnexcm("MIGRAD",arglist,2,ierflg);
	min->GetParameter(0,par1,par1E);
	fLin->SetParameter(0,HThDen);
	fLin->SetParameter(1,par1*HThDen);
	cout<<"P = "<<par1<<endl;
	auto fLin4 = (TF1*)fLin->Clone("fLin4");
	fLin4->Draw("same");
*/
	
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	gPad->SetMargin(0.15,0.10,0.15,0.05);
	gPh3->Draw("AP");
	gPh3->GetXaxis()->SetTitleSize(0.07);
	gPh3->GetXaxis()->SetLabelSize(0.05);
	gPh3->GetYaxis()->SetTitleSize(0.05);
	gPh3->GetYaxis()->SetLabelSize(0.05);
	gPh3->GetYaxis()->SetNdivisions(5);
//	gPh3->GetYaxis()->SetRangeUser(0,150);
//	gPh3->SetTitleOffset(0.7);
	TLatex* lt = new TLatex(-0.7,40,Form("1 + #alpha_{#Xi}#alpha_{#Lambda}cos#phi_{#alpha}"));
	lt->SetTextSize(0.1);
	lt->Draw();
	c2->SaveAs("LambdaPolarization.png");
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	gPad->SetMargin(0.15,0.05,0.15,0.05);
	gCM->Draw("AP");	
	gCM->GetXaxis()->SetTitle("cos#theta_{K^{-}K^{+},CM}");
	gCM->GetXaxis()->SetTitleSize(0.07);
	gCM->GetXaxis()->SetLabelSize(0.05);
	gCM->GetYaxis()->SetTitleSize(0.05);
	gCM->GetYaxis()->SetLabelSize(0.05);
	gCM->GetYaxis()->SetNdivisions(5);
	gCM->GetYaxis()->SetTitle(YTitleCM);
	c3->SaveAs("CosCM.png");
//	gCM->SetLineColor(kBlack);
	cout<<Form("CopRat: %d/%d",(int)HCM->GetEntries(),(int)HMM->GetEntries())<<endl;
//	TCanvas* c4 = new TCanvas("c4","c4",600,600);
/*
	gPh3->Draw("ap");
	gPh3->Fit("fLin");
	par0 = fLin->GetParameter(0)/bins;
	par1 = fLin->GetParameter(1)/bins;
	double par1Er = fLin->GetParError(1)/bins;
	cout<<"N = "<<par0<<endl;
	cout<<"P = "<<par1<<" +- "<<par1Er<<endl;
*/
}
