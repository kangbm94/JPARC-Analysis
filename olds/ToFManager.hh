#include "Utils.hh"
#include "stof.hh"
static const int tofnseg = 24;
class ToFManager{
	private:
		fstream pf;
		TFile* file;
		TFile* hist_file;
		TChain* chain = new TChain("kurama");
		TChain* khodo;
		TChain* kk;
	public:
		ToFManager(){};
		void LoadFile(TString filename);
		void SaveHisto(TString filename);
		void LoadChain(TString filename);
		void LoadKHodo();
		void LoadKK();
		void MakeParameterFile(string title);
		void WriteParameter(int Cid, int Plid, int Segid, int AorT, int UorD, double p0, double p1);
		void WriteParameter(int Cid, int Plid, int Segid, int Type, int UorD, int nParam, double p0, double p1,double p2);
		void NextParameter();
		TH1* GetHisto(int seg, int UD);
		void GetADCHisto(TString ht, int seg, int UD);
		void GetPEDHisto(TString ht, int seg, int UD);
		TH1* GetStofHisto(int seg, int particle);
		TH2D* GetPHCHisto(int seg, int UD,int particle);
		void FitTimewalk(int seg, int UD,int particle,bool mod);
		void Get1DHistoFromNumber(int num){
			TH1* htemp = (TH1*)file->Get(Form("h%d",num));
			htemp->Draw("colz");
		}
		void Get2DHistoFromNumber(int num){
			TH2* htemp = (TH2*)file->Get(Form("h%d",num));
			htemp->Draw("colz");
		}
		void DrawFromKHodo(TString Data,TCut Cut, TString Options = ""){
			khodo->Draw(Data,Cut, Options);
		}
		void DrawFromKK(TString Data,TCut Cut, TString Options = ""){
			kk->Draw(Data,Cut, Options);
		}
};
TCut VertexCut(vector <double> Origin, vector <double> Size){
	TCut CutX = Form("abs(vtx-%f)<%f",Origin[0],Size[0]);
	TCut CutY = Form("abs(vty-%f)<%f",Origin[1],Size[1]);
	TCut CutZ = Form("abs(vtx-%f)<%f",Origin[2],Size[2]);
	return CutX&&CutY&&CutZ;
}
void ToFManager::LoadFile(TString filename){
	file = new TFile(filename,"READ");
}
void ToFManager::SaveHisto(TString filename){
	hist_file = new TFile(filename,"recreate");
}
void ToFManager::LoadChain(TString filename){
	chain->Add(filename);
}
void ToFManager::LoadKHodo(){
	khodo = (TChain*)file->Get("khodo");
}
void ToFManager::LoadKK(){
	kk = (TChain*)file->Get("kk");
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
TH2D* ToFManager::GetPHCHisto(int seg, int UD,int particle){
	int hn = 4*10000+100*seg +10*UD+ particle+1;
	TH2D* h = (TH2D*)file->Get(Form("h%d",hn));
	return h;
};
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
void ToFManager::FitTimewalk(int seg, int UD, int particle,bool mod){
	double ecut = 0.2,emax=1.5;
	if(seg==7&&UD==0) emax = 1.3;
	if(seg==7&&UD==1) ecut = 0.6;
	if(seg==8&&UD==1) ecut = 0.8;
	if(seg==9&&UD==0) ecut = 0.5;
	if(seg<5) ecut=0.6;
	if(seg==24&&UD==0) emax=1.8;
	if(seg==24&&UD==1) emax=1.8;
	TString ct = Form("ToF%d_%d",seg,UD);
	TCanvas* c1 = new TCanvas(ct,ct,1200,600);
	TF1* slew_func= new TF1("slew_func",SlewFunc,ecut,5,3);
	slew_func->SetRange(ecut,emax);
	double p0min = 0,p0max=5,p1min=-3,t0=-3,t1=3;
	slew_func->SetParLimits(0,-5,0);
	slew_func->SetParLimits(1,p1min,ecut);
	slew_func->SetParLimits(2,t0,t1);
	TH2D* bfh = GetPHCHisto(seg,UD,particle); 
	slew_func->SetRange(ecut,4);
	slew_func->SetParLimits(0,-5,0);
	slew_func->SetParLimits(1,p1min,ecut);
	slew_func->SetParLimits(2,t0,t1);
	if(mod){
		bfh->Fit("slew_func","QR");
		bfh->Fit("slew_func","R");
		double p0 = slew_func->GetParameter(0);
		double p1 = slew_func->GetParameter(1);
		double p2 = slew_func->GetParameter(2);
		hist_file->cd();
		bfh->Write();
		WriteParameter(7,0,seg-1,UD,1,3,-p0,p1,-p2);
	}
	else{
		double epb=0.04;
		int bps=3;
		int bin1 = ecut/epb;
		int bin2 = emax/epb;
		double xp[100];
		double yp[100];
		int cnt=0;
		for(int i=bin1;i<bin2+1;i++){
			xp[cnt]=epb*i;
			TH1D* hy = bfh->ProjectionY("py",i-bps/2,i+bps/2);
			yp[cnt]=hy->GetBinCenter(hy->GetMaximumBin());
			cnt++;
		}
		TGraph* gr = new TGraph(cnt,xp,yp);
		gr->SetTitle(ct);
		gr->Draw("APsame");
		//gr->SetMinimum(-3);
		//gr->SetMaximum(3);
		gr->GetXaxis()->SetLimits(0,4);
		gr->GetYaxis()->SetRangeUser(-3,3);
		bfh->Draw("same");
		gr->SetMarkerStyle(2);
		gr->SetMarkerSize(2);
		slew_func->SetLineWidth(4);
		gr->Fit("slew_func","QR");
		gr->Fit("slew_func","R");
		double p0 = slew_func->GetParameter(0);
		double p1 = slew_func->GetParameter(1);
		double p2 = slew_func->GetParameter(2);
		hist_file->cd();
		c1->Write();
		WriteParameter(7,0,seg-1,UD,1,3,-p0,p1,-p2);
	}
}

TH1* ToFManager::GetStofHisto(int seg, int particle){
	int hn = 10000+100*seg+particle+1;
	TH1* h = (TH1*)file->Get(Form("h%d",hn));
	h->SetAxisRange(-2,2);
	return h;
}
