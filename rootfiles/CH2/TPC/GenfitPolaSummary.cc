#include "PolarizationAnal.hh"
#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
TH1D* HTh;
TH1D* HPh;
TH1D* HPh1;
TH1D* HPh2;
TH1D* HPh3;
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
double MissMass,Coplanarity;
bool TrigAPS,TrigB;
void RealPola(int i);
TTree* tree2;
vector<double>datap;
void FCN(int &npar, double* gin, double &f, double * par, int iflag){
	f=0;
	double par1 = par[0];
	for(auto dp:datap){
		double delta = 0.5 - dp * par1;
		f+= delta* delta;
	}
//	cout<<Form("ND = %lu,chi2 = %g, par = %g",datap.size(),f,par1)<<endl;
}
void GenfitPolaSummary(){
	int nbin = 20;
	SetStyle();
	gStyle->SetTitleSize(0.08,"XY");
	HTh = new TH1D("HistTh","HistTh",nbin,-1.,1.);
	HPh = new TH1D("HistPh","HistPh",nbin,-1.,1.);
	HPh1 = new TH1D("HistPh1","HistPh1",nbin,-1.,1.);
	HPh2 = new TH1D("HistPh2","HistPh2",nbin,-1.,1.);
	HPh3 = new TH1D("HistPh3","HistPh3",nbin,-1.,1.);
	TF1* fLin = new TF1("fLin","[0]+[1]*x",-1,1);
	TF1* fPLin = new TF1("fPLin","0.5+[0]*x",-1,1);
	
	TH1D* HMM = new TH1D("HistMM","HistMM",100,1.2,1.7);
	TH1D* HMMCut = new TH1D("HistMMCut","HistMMCut",100,1.2,1.7);
	int rn_s = 5641;
	int rn_e = 5666;
//	int rn_s = 5667;
//	int rn_e = 5697;
	TFile* file = new TFile(Form("GenfitPolaAnal_%d_%d.root",rn_s,rn_e));
	TTree* tree = (TTree*)file->Get("tree");
	tree->SetBranchAddress("TrigAPS",&TrigAPS);
	tree->SetBranchAddress("TrigB",&TrigB);
	tree->SetBranchAddress("MissMass",&MissMass);
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
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		if(abs(InvMLd-1.115) > 0.03 or abs(InvMXi-1.321)>0.03) continue;
		if(!isnan(cTh) and abs(Coplanarity)< 0.05	and abs(MissMass - 1.321)<0.1){
			HTh->Fill(cTh);
			HPh->Fill(cPh);
			HPh1->Fill(cPh1);
			HPh2->Fill(cPh2);
			HPh3->Fill(cPh3);
			HMMCut->Fill(MissMass);
			VecTh.push_back(cTh);
			VecPh.push_back(cPh);
			VecPh1.push_back(cPh1);
			VecPh2.push_back(cPh2);
			VecPh3.push_back(cPh3);
		}
		if((MissMass - 1.321)<0.1){
			HMM->Fill(MissMass);
		}
	}
	HTh->GetXaxis()->SetRangeUser(-1.1,1.1);
	HPh->GetXaxis()->SetRangeUser(-1.1,1.1);
	HPh1->GetXaxis()->SetRangeUser(-1.1,1.1);
	HPh2->GetXaxis()->SetRangeUser(-1.1,1.1);
	HPh3->GetXaxis()->SetRangeUser(-1.1,1.1);
	TMinuit* min = new TMinuit(1);
	min->SetPrintLevel(-1);
	int ierflg = 0;
	double arglist[2]={500,0.01};
	double par1,par1E;
	datap = VecTh;
	min->SetFCN(FCN);	
	min->mnparm(0,"p1",0.,0.01,-1,1,ierflg);
	min->mnexcm("MIGRAD",arglist,2,ierflg);
	min->GetParameter(0,par1,par1E);
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	c1->cd(1);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	HTh->Draw();
	HTh->GetXaxis()->SetTitle("cos#theta");
	HTh->GetYaxis()->SetTitle("Counts / 0.05 ");
	//	HTh->Fit("fLin");
	double HThDen = (HTh->GetEntries() / nbin);
	fLin->SetParameter(0,HThDen);
	fLin->SetParameter(1,par1*HThDen);
	cout<<"P = "<<par1<<endl;
	auto fLin1 = (TF1*)fLin->Clone("fLin1");
	fLin1->Draw("same");
	c1->cd(2);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	HPh1->Draw();
	HPh1->GetXaxis()->SetTitle("cos#phi_{#beta}");
	HPh1->GetYaxis()->SetTitle("Counts / 0.05 ");
	datap = VecPh1;
	min->mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
	min->mnexcm("MIGRAD",arglist,2,ierflg);
	min->GetParameter(0,par1,par1E);
	fLin->SetParameter(0,HThDen);
	fLin->SetParameter(1,par1*HThDen);
	cout<<"P = "<<par1<<endl;
	auto fLin2 = (TF1*)fLin->Clone("fLin2");
	fLin2->Draw("same");
	c1->cd(3);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	HPh2->Draw();
	HPh2->GetXaxis()->SetTitle("cos#phi_{#gamma}");
	HPh2->GetYaxis()->SetTitle("Counts / 0.05 ");
	datap = VecPh2;
	min->mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
	min->mnexcm("MIGRAD",arglist,2,ierflg);
	min->GetParameter(0,par1,par1E);
	fLin->SetParameter(0,HThDen);
	fLin->SetParameter(1,par1*HThDen);
	cout<<"P = "<<par1<<endl;
	auto fLin3 = (TF1*)fLin->Clone("fLin3");
	fLin3->Draw("same");

	c1->cd(4);
	gPad->SetMargin(0.15,0.1,0.2,0.1);
	HPh3->Draw();
	HPh3->GetXaxis()->SetTitle("cos#phi_{#alpha}");
	HPh3->GetYaxis()->SetTitle("Counts / 0.05 ");
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

	
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	HMM->Draw();
	HMMCut->SetLineColor(kRed);
	HMMCut->Draw("same");
	cout<<Form("CopRat: %d/%d",(int)HMMCut->GetEntries(),(int)HMM->GetEntries())<<endl;
}
