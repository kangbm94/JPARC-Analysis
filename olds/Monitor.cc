#include "HodoscopeManager.hh"
void Monitor(){
	cout<<"SortRuns(int run_start,int run_end)"<<endl;
	cout<<"MonitorRuns()"<<endl;
}
void MonitorRuns(){
	TChain* chain = new TChain("run_monitor");
	TString f1= "RunMonitor5500_5599.root";
	TString f2= "RunMonitor5600_5700.root";
	chain->Add(f1);
//	chain->Add(f2);
	const int ent = chain->GetEntries();
	TGraphErrors* gu;
//	TGraphErrors* gd;
	int run;double gain[16];double gainerr[16];
	chain->SetBranchAddress("runnum",&run);
	chain->SetBranchAddress("MPV",gain);
	chain->SetBranchAddress("MPVErr",gainerr);
	double runs[ent];double mpv[ent];double mpve[ent];double re[ent];
	for(int i=0;i<ent;i++){
		chain->GetEntry(i);
		re[i]=0.5;
		runs[i]=run;
		mpv[i]=gain[0];
		mpve[i]=gainerr[0];
	}
	gu=new TGraphErrors(ent,runs,mpv,re,mpve);
	gu->GetYaxis()->SetRangeUser(0,2000);
	gu->Draw();
}
void SortRuns(int run_start, int run_end){
	HodoscopeManager h;
	double gain[16];double sig[16];double chi2[16];int ndf[16];double ped[16];
	h.MakeMonitorFile(Form("RunMonitor%d_%d.root",run_start,run_end));
	for(int i=run_start;i<run_end+1;i++){
		h.MonitorBH2(i,0,gain,sig,chi2,ndf,ped);
	}
	h.WriteMonitorFile();
}
void SortRun(int run){
  	gStyle->SetOptFit(1111);
	HodoscopeManager h;
	double gain[16];double sig[16];double chi2[16];int ndf[16];double ped[16];
	TFile* file = new TFile(Filename(run),"READ");
	if(file->IsZombie()){
		gSystem->Exec(Form("Run %d >> ZombieRunLog.txt",run));
		return;
	}
	if(file->GetSize()<1e6){
		cout<<"File size too small"<<endl;
		return;
	}
	cout<<"Creating file "<<Form("../rootfiles/RunMonitor%d.root",run)<<endl;
	h.MakeMonitorFile(Form("../rootfiles/RunMonitor%d.root",run));
	h.MonitorBH2(run,0,gain,sig,chi2,ndf,ped);
	h.WriteMonitorFile();
}
void Monitor(int run){
	SortRun(run);
}
