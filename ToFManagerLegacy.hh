#include "Utils.hh"
#include "stof.hh"
static const int tofnseg = 24;
static const double ToFPos0 = -900;
static const double ToFWidth = 80;
static const double ToFPitch = 75;

TCut XCut(int seg){
	double pos0 = ToFPos0 +ToFPitch*(seg-1);
	return Form("xtofKurama>%f&&xtofKurama<%f",pos0,pos0+ToFWidth);
}
TCut YCut(int seg,double ycut=100){
	return Form("ytofKurama>%f&&ytofKurama<%f",-ycut,ycut);
}
TCut YTCut(int seg,double slope, double offset){
	return Form("utTofSeg[%d]-%f*ytofKurama<%f",seg-1,slope,offset);
}
TCut MTCut(int seg, double offset){
	return Form("utTofSeg[%d]/2+dtTofSeg[%d]/2<%f",seg-1,seg-1,offset);
}
TString BF2D(int seg, int UD, TString ht){
	TString val;
	val=Form("%stTofSeg[%d]",UorD(UD).Data(),seg-1) + (TString)":"+Form("%sdeTofSeg[%d]",UorD(UD).Data(),seg-1)+(TString)">>"+ht;
	cout<< "BF2D : "<<val<<endl;
	return val;
};
TString AF2D(int seg, int UD, TString ht,double p0,double p1,double p2){
	TString val;
	val=Form("%stTofSeg[%d]",UorD(UD).Data(),seg-1) + (TString)Form("-(%f/sqrt(%sdeTofSeg[%d]-%f))+%f",p0,UorD(UD).Data(),seg-1,p1,p2)+":"+Form("%sdeTofSeg[%d]",UorD(UD).Data(),seg-1)+">>"+ht;
	return val;
};
class ToFManager{
	private:
		fstream pf;
		TFile* file;
		TFile* hist_file;
		TChain* chain = new TChain("kurama");
		double ucut[tofnseg] = {1,		1,	1,	1,	1,	1,  //1 	2		3 	4 	5 	6 
														-0.06,-0.04,	-0.01,	1,	1,	1,	//7		8		9		10	11	12
															1,	1,	1,	1,	1,	1,	//13	14	15	16	17	18
															1,	1,	1,	1,	1,	1};	//19	20	21	22	23	24
	public:
		ToFManager(){};
		void LoadFile(TString filename);
		void SaveHisto(TString filename);
		void LoadChain(TString filename);
		void MakeParameterFile(string title);
		void WriteParameter(int Cid, int Plid, int Segid, int AorT, int UorD, double p0, double p1);
		void WriteParameter(int Cid, int Plid, int Segid, int Type, int UorD, int nParam, double p0, double p1,double p2);
		void NextParameter();
		TH1* GetHisto(int seg, int UD);
		void GetADCHisto(TString ht, int seg, int UD);
		void GetPEDHisto(TString ht, int seg, int UD);
		void YTPositionFit(int seg,int UD); 
		TH1* GetStofHisto(int seg, int particle);


		void FitTimewalk(int seg, int UD,double offset);
};
void ToFManager::LoadFile(TString filename){
	file = new TFile(filename,"READ");
}
void ToFManager::SaveHisto(TString filename){
	hist_file = new TFile(filename,"recreate");
}
void ToFManager::LoadChain(TString filename){
	chain->Add(filename);
}
void ToFManager::MakeParameterFile(string title){
	pf.open(title,fstream::out);
}
void ToFManager::NextParameter(){
	pf<<"###############"<<endl;
}
void ToFManager::WriteParameter(int Cid, int Plid, int Segid, int AorT, int UorD, double p0, double p1){
	pf<<Cid<<"\t"	;
	pf<<Plid<<"\t"	;
	pf<<Segid<<"\t"	;
	pf<<AorT<<"\t"	;
	pf<<UorD<<"\t"	;
	pf<<p0<<"\t"	;
	pf<<p1<<endl	;
}
void ToFManager::WriteParameter(int Cid, int Plid, int Segid, int Type, int UorD ,int nParam, double p0, double p1,double p2){
	pf<<Cid<<"\t"	;
	pf<<Plid<<"\t"	;
	pf<<Segid<<"\t"	;
	pf<<Type<<"\t"	;
	pf<<UorD<<"\t"	;
	pf<<nParam<<"\t"	;
	pf<<p0<<"\t"	;
	pf<<p1<<"\t"	;
	pf<<p2<<endl	;
}
TH1* ToFManager::GetHisto(int seg, int UD){// 0 -> Up, 1 -> Down
	int hn = 6*10000+100*seg + 3+UD;
	TH1* h = (TH1*)file->Get(Form("h%d",hn));
	int peak = h->GetBinContent(h->GetMaximumBin());
	int peak_ref = 20;
	int nbin = 200/peak + 1;
	if(nbin> 1)	h->Rebin(nbin);
	return h;
}
void ToFManager::GetADCHisto(TString ht,int seg, int UD){// 0 -> Up, 1 -> Down
	TString hn ;
	if(UD==0){
		hn = Form("tofua[%d]",seg-1);
	}
	else{
		hn = Form("tofda[%d]",seg-1);
	}
	cout<<"Drawing : " << hn <<endl;
	TCut cut = Form("tofhitpat==%d",seg);
	TTree* tr = (TTree*)file->Get("tree");
	tr->Draw(hn+">>"+ht,cut);
}
void ToFManager::GetPEDHisto(TString ht,int seg, int UD){// 0 -> Up, 1 -> Down
	TString hn ;
	if(UD==0){
		hn = Form("tofua[%d]",seg-1);
	}
	else{
		hn = Form("tofda[%d]",seg-1);
	}
	cout<<"Drawing : " << hn <<endl;
	TCut cut = Form("tofhitpat!=%d",seg);
	TTree* tr = (TTree*)file->Get("tree");
	tr->Draw(hn+">>"+ht,cut);
}
void ToFManager::YTPositionFit(int seg, int UD){
	TF1* f_lin = new TF1("f_lin","pol1",-500,500);
	TString dr = Form("%stTofSeg[%d]:ytofKurama>>(100,-1000,1000,100,10,25)",UorD(UD).Data(),seg-1);
	TCut dcut = Form("ntSdcOut==1&&tofsegKurama==%d",seg);
}

void ToFManager::FitTimewalk(int seg, int UD,double offset){
	double mbf,maf;
	double t0=10,t1=18;
	double ecut=0.6,p0max=15,p0min=0.3,p1min=-5,p2min= 30;
	TF1* slew_func= new TF1("slew_func",SlewFunc,ecut,3,3);
	slew_func->SetRange(ecut,5);
	slew_func->SetParLimits(0,p0min,p0max);
	slew_func->SetParLimits(1,p1min,ecut);
	slew_func->SetParLimits(2,t0-20,t1+10);
	TString bhn = Form("TOF%s%d_before",UorD(UD).Data(),seg);
	TString ahn = Form("TOF%s%d_after",UorD(UD).Data(),seg);
	TH2D* bfh = new TH2D(bhn,bhn,100,0,5,1000,t0,t1);
	TH2D* afh = new TH2D(ahn,ahn,100,0,5,1000,(t0-t1)/2,(t1-t0)/2);
	TString ct = Form("c%d",canv);
	TCanvas* c1 = new TCanvas(ct,ct,1200,600);
	c1->Divide(2,1);
//	TCut cut = Form("TofSeg==%d",seg);
	TCut cut = Form("tofsegKurama==%d",seg);
//	cut= cut&&"ntSdcOut==1";
//	cut = cut&&Form("utofKurama<%f",ucut[seg-1]);
//	cut = cut&&Form("pKurama<%f",1.);
	cut = cut&&"nhTof==1";
	cut=cut&&XCut(seg);//&&YCut(seg);
//	if(seg>6)	cut=cut&&YCut(seg,100);
//	cut=cut&&YCut(seg,50);
//	cut=cut&&YTCut(seg,0.0064,offset);
	cut=cut&&MTCut(seg,offset);
	c1->cd(1);
	hist_file->cd();
	chain->Draw(BF2D(seg,UD,bhn),cut,"colz");
	mbf = bfh->GetMean(2);
	cout<< "bhn : "<<bhn<<endl;
	slew_func->SetRange(ecut,5);
	slew_func->SetParLimits(0,p0min,p0max);
	slew_func->SetParLimits(1,p1min,ecut);
	slew_func->SetParLimits(2,t0-p2min,t1);
	bfh->Fit("slew_func","QR");
	bfh->Fit("slew_func","R");
	double p0 = slew_func->GetParameter(0);
	double p1 = slew_func->GetParameter(1);
	double p2 = slew_func->GetParameter(2);
	c1->cd(2);
	chain->Draw(AF2D(seg,UD,ahn,p0,p1,p2),cut,"colz");
	maf = afh->GetMean(2);
	p2+=(mbf-maf);
	bfh->Write();
	afh->Write();
	WriteParameter(7,0,seg-1,UD,1,3,p0,p1,p2);
}

TH1* ToFManager::GetStofHisto(int seg, int particle){
	int hn = 10000+100*seg+particle+1;
	TH1* h = (TH1*)file->Get(Form("h%d",hn));
	h->SetAxisRange(0,10);
	int peak = h->GetBinContent(h->GetMaximumBin());
	int peak_ref = 10;
	int nbin = 10/peak + 1;
	if(nbin> 1)	h->Rebin(nbin);
	h->SetAxisRange(0,10);
	return h;
}
