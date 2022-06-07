#include "Utils.hh"
#include "ParLimits.hh"
int bh2seg=1,bh1seg=1,htofseg=20,tofseg=1;
bool DrawCanvas=false;
class HodoscopeManager{
	private: 
		static const int bh1nseg=11; 
		static const int bh2nseg=8; 
		static const int htofnseg=34;
		static const int mtdepth= 16;
		static const int corrections=4;
		int line = 0;
		bool baccut=false;
		double lsb;
		double range;
		TFile* monfile;
		TTree* runtr;
		TFile* slewfile;
		int rn;double gn[16];double gne[16];double sg[16];double sge[16];double ch2[16];int ndfs[16];double pd[16];int hit[16];double cbtof;
	fstream f;
	fstream fo;
	fstream fs;
	public:
		HodoscopeManager(){
		}
		void MakePHCFile(string title,int conf, int runnum);
		void MakePHCFile(string title);
		void LoadSPHCFile(string title){
			fs.open(title,fstream::in);
	if(!fs){
		cout<<"File not found"<<endl;
		exit(0);
	}
		}
		void LoadPHCFile(string title){
			fo.open(title,fstream::in);
	if(!fo){
		cout<<"File not found"<<endl;
		exit(0);
	}
		}
		void BACCut(bool cut);
		void ScanSegment(TFile* file,int bh2segment, int bacconf, int* segs);
		void ScanSegments(TFile* file,int bacconf,int* bh1segs,int* htofsegs);
		void ScanSegment(TChain* chain,int bh2segment, int bacconf, int* segs);
		void ScanSegments(TChain* chain,int bacconf,int* bh1segs,int* htofsegs);
		void FitTimewalk(TChain* chain, int bac,int det, int seg,int conf);
		void MakeMonitorFile(TString name);
		void MonitorBH2(int run, int bac,double* gain,double* sig,double* chi2,int* ndf,double* ped); 
		void WriteMonitorFile();
		void BH1_U2D();
		void BH1_2_BH2();
		void BH2_U2D();
		void OtherLine();
		void AddSTOF();
};
void HodoscopeManager::MakeMonitorFile(TString name){
	monfile = new TFile(name,"recreate");
	runtr= new TTree("run_monitor","run_monitor");
	runtr->Branch("runnum",&rn,"runnum/I");
	runtr->Branch("MPV",gn,"MPV[16]/D");
	runtr->Branch("MPVErr",gne,"MPVErr[16]/D");
	runtr->Branch("Sig",sg,"Sig[16]/D");
	runtr->Branch("SigErr",sge,"SigErr[16]/D");
	runtr->Branch("Chi2",ch2,"Chi2[16]/D");
	runtr->Branch("NDF",ndfs,"NDF[16]/I");
	runtr->Branch("ped",pd,"ped[16]/D");
	runtr->Branch("cbtof",&cbtof,"cbtof/D");
}
void HodoscopeManager::WriteMonitorFile(){
	monfile->Write();
}

void HodoscopeManager::MonitorBH2(int run, int bac,double* gain, double* sig, double* chi2,int* ndf,double* ped){
	gStyle->SetOptFit(1111);
  	if(gSystem->AccessPathName(Filename(run))){
		cout<<"File "+Filename(run)+" do not exist!"<<endl;
		return;
	}
	TFile* file = new TFile(Filename(run),"READ");
	if(file->IsZombie()){
		gSystem->Exec(Form("Run %d >> ZombieRunLog.txt",run));
		return;
	}
	cout<<Form("Run%d",run)<<endl;
	if(file->GetSize()<1e6){
		cout<<"File size too small"<<endl;
		return;
	}
	TDirectory *dir = monfile->mkdir(Form("Hist_Run%d",run));
	dir->cd();
	cout<<"load tree"<<endl;
	TTree* chain = (TTree*)file->Get("tree");
	int ents = chain->GetEntries();
	if(ents<10000){
		cout<<"Entries too small"<<endl;
		return;
	}
	cout<<ents<<endl;
	TString lt = Form("bh2_landau_run_%d",run);
	TCanvas* c1 = new TCanvas("cl",lt,1800,1200);
	c1->Divide(4,4);
	bac=1-bac;
//	TCut cut = Form("bacnhits!=%d",bac);
	TCut cut="";
	TH1D* h[16];
	TCut entcut[8];
	double sig_cut=2;	
	double sig_range=10;	
	int rb;
	double mpv_tmp,sig_tmp;
	TF1* land = new TF1("land","landau",0,4000);
	TF1* pedf = new TF1("pedf","gaus",0,4000);
	TF1* gausf = new TF1("gausf","gaus",-3,3);
	TH1D* cbtof_hist = new TH1D(Form("cbtof_%d",run),Form("cbtof_%d",run),1000,-5,5);
	cout<<"process data"<<endl;
	TCanvas* c3 = new TCanvas("c3","cbtof",1800,1200);
	for(int i=0;i<8;i++){
		entcut[i] = Form("bh2hitpat[0]==%d",i+1);
		TString uht = Form("BH2UADC[%d]_run%d",i+1,run);
		h[i]=new TH1D(uht,uht,3600,0,3600);
		c1->cd(i+1);
		chain->Draw(Form("bh2ua[%d]>>",i)+uht,cut&&entcut[i]);
		rb=RebinHist(h[i]);
		h[i]->Rebin(rb);
		h[i]->Draw();
		land->SetRange(0,4000);
		h[i]->Fit("land","QRM");
		mpv_tmp=land->GetParameter(1);
		sig_tmp=land->GetParameter(2);
		land->SetRange(mpv_tmp-sig_cut*sig_tmp,mpv_tmp+sig_range*sig_tmp);
		gain[i]=land->GetParameter(1);
		gne[i]=land->GetParError(1);
		sig[i]=land->GetParameter(2);
		sge[i]=land->GetParError(2);
		chi2[i]=land->GetChisquare();
		ndf[i]=land->GetNDF();
		h[i]->Fit("land","QRM");
		
		TString dht = Form("BH2DADC[%d]_run%d",i+1,run);
		h[8+i]=new TH1D(dht,dht,3600,0,3600);
		h[8+i]->Rebin(rb);
		c1->cd(i+9);
		chain->Draw(Form("bh2da[%d]>>",i)+dht,cut&&entcut[i]);
		land->SetRange(0,4000);
		h[8+i]->Fit("land","QRM");
		mpv_tmp=land->GetParameter(1);
		sig_tmp=land->GetParameter(2);
		land->SetRange(mpv_tmp-sig_cut*sig_tmp,mpv_tmp+sig_range*sig_tmp);
		h[8+i]->Fit("land","QRM");
		gain[8+i]=land->GetParameter(1);
		gne[8+i]=land->GetParError(1);
		sig[8+i]=land->GetParameter(2);
		sge[8+i]=land->GetParError(2);
		chi2[8+i]=land->GetChisquare();
		ndf[8+i]=land->GetNDF();
	}
	c1->Modified();
	c1->Update();
	c1->SetBatch(kTRUE);
	gSystem->ProcessEvents();
	c3->cd();
	chain->Draw(Form("cbtof>>cbtof_%d",run));
	cbtof_hist->Fit("gausf","QRM");
	cbtof=gausf->GetParameter(2);
	
	TString pt = Form("bh2_ped_run_%d",run);
	TCanvas* c2 = new TCanvas("cp",pt,1800,1200);
	c2->Divide(4,4);
	TH1D* ph[16];
	double pedrange=10;
	for(int i=0;i<8;i++){
//		entcut[i] = Form("bh2hitpat[0]!=%d&&bh2hitpat[0]!=%d&&bh2hitpat[0]!=%d",i+1,i,i+2);
		entcut[i] = Form("bh2hitpat[0]!=%d",i+1);
		TString uht = Form("BH2UPED[%d]_run%d",i+1,run);
		ph[i]=new TH1D(uht,uht,720,0,720);
		c2->cd(i+1);
		chain->Draw(Form("bh2ua[%d]>>",i)+uht,cut&&entcut[i]);
		rb=RebinHist(ph[i]);
		ph[i]->Rebin(rb);
		ph[i]->Draw();
		pedf->SetRange(upm[i]-pedrange,upm[i]+pedrange);
		ph[i]->Fit("pedf","QRML");
		ped[i]=pedf->GetParameter(1);	
		TString dht = Form("BH2DPED[%d]_run%d",i+1,run);
		ph[8+i]=new TH1D(dht,dht,720,0,720);
		ph[8+i]->Rebin(rb);
		c2->cd(i+9);
		chain->Draw(Form("bh2da[%d]>>",i)+dht,cut&&entcut[i]);
		pedf->SetRange(dpm[i]-pedrange,dpm[i]+pedrange);
		ph[8+i]->Fit("pedf","QRML");
		ped[8+i]=pedf->GetParameter(1);	
	}
	c2->Modified();
	c2->Update();
	c2->SetBatch(kTRUE);
	gSystem->ProcessEvents();
	rn=run;
	for(int i=0;i<16;i++){
		gn[i]=gain[i];
		sg[i]=sig[i];
		ch2[i]=chi2[i];
		ndfs[i]=ndf[i];
		pd[i]=ped[i];
	}
	runtr->Fill();
	delete chain;
}
void HodoscopeManager::FitTimewalk(TChain* chain, int bac,int det, int seg,int conf){
	double mbf,maf;
	TString bhn = Detname[det]+UorD(conf)+(TString)to_string(seg)+"_before";
	TString ahn = Detname[det]+UorD(conf)+(TString)to_string(seg)+"_after";
	TH2D* bfh = new TH2D(bhn,bhn,100,0,5,1000,-5,5);
	TH2D* afh = new TH2D(ahn,ahn,100,0,5,1000,-5,5);
	TString ct = Form("c%d",canv);
	TCanvas* c1 = new TCanvas(ct,ct,1200,600);
	c1->Divide(2,1);
	TF1* slew_func= new TF1("slew_func",SlewFunc,ecut,5,3);
	Segment[det]=seg;
	bac=1-bac;
	TCut baccut = Form("bacnhits!=%d",bac);
	slew_func->SetRange(ecut,5);
	slew_func->SetParLimits(0,p0min,p0max);
	slew_func->SetParLimits(1,-p1min,ecut);
	slew_func->SetParLimits(2,-5,5);
	c1->cd(1);
	cout<<BF2D(det,bhn,conf)<<endl;
	chain->Draw(BF2D(det,bhn,conf),TimeCut(det,cutl,cuth)&&baccut,"colz");
	bfh->Fit("slew_func","R");
	mbf=bfh->GetMean(2);
	double p1;
	double p2;
	double p3;
	p1=slew_func->GetParameter(0);	
	p2=slew_func->GetParameter(1);	
	p3=slew_func->GetParameter(2);	
	c1->cd(2);
	//	cout<<AF2D(det,ahn,p1,p2,p3,conf)<<endl;
	chain->Draw(AF2D(det,ahn,p1,p2,p3,conf),TimeCut(det,cutl,cuth)&&baccut,"colz");
	maf=afh->GetMean(2);
	p3=p3+mbf-maf;
	chain->Draw(AF2D(det,ahn,p1,p2,p3,conf),TimeCut(det,cutl,cuth)&&baccut,"colz");
	int ids[6]={det,0,seg-1,conf,1,3};
	PHCWritter(f,ids,p1,p2,p3);
	line++;
	c1->Update();
	gSystem->ProcessEvents();
	canv++;
}

void HodoscopeManager::BACCut(bool cut){
	baccut = cut;
}

void HodoscopeManager::ScanSegment(TFile* file, int bh2segment,int bacconf,int* segs){
	TTree* tree = (TTree*)file->Get("tree");
	int ncnt[bh1nseg][htofnseg]={0};
	int maxbin=0;
	TString ht= Form("2DMtx_Hodoscope_%d",bh2segment);
	TH2I* Hist=new TH2I(ht,ht,bh1nseg,1,bh1nseg+1,htofnseg,1,htofnseg+1);
	TCut cut=Form("bh1nhits>0&&htofnhits>0&&bh2hitpat[0]==%d",bh2segment);
	if(bacconf==0){
		cut=cut&&"bacnhits==0";
	}
	else if(bacconf==1){
		cut=cut&&"bacnhits!=0";
	}
	tree->Draw("htofhitpat[0]:bh1hitpat[0]>>"+ht,cut,"colz");
	for(int j=0;j<bh1nseg;j++){
		for(int k=0;k<htofnseg;k++){
			int x = j+1,y=k+1;
			ncnt[j][k] = (int)Hist->GetBinContent(x,y);
			maxbin=Max(maxbin,ncnt[j][k]);
		}
	}
	for(int j=0;j<bh1nseg;j++){
		for(int k=0;k<htofnseg;k++){
			if(ncnt[j][k]==maxbin){
				cout<<"Hodoscope segment : "<<bh2segment<<" MaxBins: "<<maxbin<<endl;
				segs[0]=j+1;segs[1]=k+1;
			}
		}
	}
}
void HodoscopeManager::ScanSegments(TFile* file,int bacconf, int* bh1,int* htof){
	TString ct = Form("c%d",canv);
	TCanvas* c = new TCanvas(ct,ct,1200,600);
	c->Divide(4,2);
	for(int i=1;i<bh2nseg+1;i++){
		c->cd(i);
		int segs[2];
		ScanSegment(file,i,bacconf,segs);
		bh1[i-1]=segs[0];
		htof[i-1]=segs[1];
		cout<<"Hodoscope Seg: "<<i<<" (BH1,HTOF) : ("<<bh1[i-1]<<" , "<<htof[i-1]<<" )"<<endl;
	}
	canv++;
}


void HodoscopeManager::ScanSegment(TChain* chain, int bh2segment,int bacconf,int* segs){
	int ncnt[bh1nseg][htofnseg]={0};
	int maxbin=0;
	TString ht= Form("2DMtx_Hodoscope_%d",bh2segment);
	TH2I* Hist=new TH2I(ht,ht,bh1nseg,1,bh1nseg+1,htofnseg,1,htofnseg+1);
	TCut cut=Form("bh1nhits>0&&htofnhits>0&&bh2hitpat[0]==%d",bh2segment);
	if(bacconf==0){
		cut=cut&&"bacnhits==0";
	}
	else if(bacconf==1){
		cut=cut&&"bacnhits!=0";
	}
	chain->Draw("htofhitpat[0]:bh1hitpat[0]>>"+ht,cut,"colz");
	for(int j=0;j<bh1nseg;j++){
		for(int k=0;k<htofnseg;k++){
			int x = j+1,y=k+1;
			ncnt[j][k] = (int)Hist->GetBinContent(x,y);
			maxbin=Max(maxbin,ncnt[j][k]);
		}
	}
	for(int j=0;j<bh1nseg;j++){
		for(int k=0;k<htofnseg;k++){
			if(ncnt[j][k]==maxbin){
				cout<<"Hodoscope segment : "<<bh2segment<<" MaxBins: "<<maxbin<<endl;
				segs[0]=j+1;segs[1]=k+1;
			}
		}
	}
}
void HodoscopeManager::ScanSegments(TChain* chain,int bacconf, int* bh1,int* htof){
	TString ct = Form("c%d",canv);
	TCanvas* c = new TCanvas(ct,ct,1200,600);
	c->Divide(4,2);
	for(int i=1;i<bh2nseg+1;i++){
		c->cd(i);
		int segs[2];
		ScanSegment(chain,i,bacconf,segs);
		bh1[i-1]=segs[0];
		htof[i-1]=segs[1];
		cout<<"Hodoscope Seg: "<<i<<" (BH1,HTOF) : ("<<bh1[i-1]<<" , "<<htof[i-1]<<" )"<<endl;
	}
	canv++;
}
void HodoscopeManager::MakePHCFile(string title){
	f.open(title,fstream::out);
}
void HodoscopeManager::MakePHCFile(string title, int conf,int runnum){
	TString dirs= "~/PHCHists/";
	TString fname = Form("PHCHists_run%d",runnum);
	slewfile = new TFile(dirs+fname+".root","recreate");
	cout << title << endl;
	f.open(title,fstream::out);
	if(!f){
		cout<<"File not created"<<endl;
		exit(0);
	}
	if(conf==1){
		f<<
			"#######################################################################"
			<<endl;line++;
		f<<
			"# CId	PlId	Seg	UorD	Type	nParam	p0	p1	..."
			<<endl;line++;
		f<<
			"#"
			<<endl;line++;
		f<<
"#	Type=0 no correction will be performed"
			<<endl;line++;
		f<<
"#	Type=1 ct=t-p[0]/sqrt(a-p[1])+p[2]"
			<<endl;line++;
		f<<
"#	Type=2 ct=p[0]*a*a + p[1]*a + p[2] for Fiber"
			<<endl;line++;
		f<<
"#"
			<<endl;line++;
		f<<
"#######################################################################"
			<<endl;line++;
		f<<
"# BH1 Up                                                              #"
			<<endl;line++;
		f<<
"#######################################################################"
			<<endl;line++;
	}
}
void HodoscopeManager::BH1_U2D(){
		f<<
"#######################################################################"
			<<endl;line++;
		f<<
"# BH1 Down                                                            #"
			<<endl;line++;
		f<<
"#######################################################################"
			<<endl;line++;
}
void HodoscopeManager::BH1_2_BH2(){
		f<<
"#######################################################################"
			<<endl;line++;
		f<<
"# BH2 Up                                                              #"
			<<endl;line++;
		f<<
"#######################################################################"
			<<endl;line++;
}
void HodoscopeManager::BH2_U2D(){
		f<<
"#######################################################################"
			<<endl;line++;
		f<<
"# BH2 Down                                                            #"
			<<endl;line++;
		f<<
"#######################################################################"
			<<endl;line++;
}
void HodoscopeManager::OtherLine(){
	string lines;
	for(int i=0;i<line;i++){
		getline(fo,lines);
	}
	while(getline(fo,lines)){
		f<<lines<<endl;line++;
	}
	slewfile->Write();
}
void HodoscopeManager::AddSTOF(){
	int phclines=57;
	string line_phc;
	string line_stof;
	for(int i=0;i<phclines;i++){
		getline(fs,line_stof);//fs = HodoPHCParam_StofCorrection
		getline(fo,line_phc);//fo=HodoPHCParamRBR
		f<<line_phc<<endl;
	}
	while(getline(fs,line_stof)){
		f<<line_stof<<endl;
	}
}
