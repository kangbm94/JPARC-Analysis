double cbtof[11][8];
int bacnhits;
int trigflag[32];
int trigpat[32];
TChain* gchain;
TTree* gtree;
TFile* gfile;
TFile* gfileOut;
TString gcut;
int KCount,PiCount,BGCount,rn;
double KTime,PiTime;
double KWidth,PiWidth,BGWidth;
TH1D* Hist_tof[2000];
double dgaus(double* x, double* p){
	double val=0;
	val=Gaussian(x[0],p[0],p[1],p[2])+Gaussian(x[0],p[3],p[4],p[5]);//Mean,sig, amp
	val+=Gaussian(x[0],p[0],p[6],p[7]);
	return val;
}
double sgaus(double* x, double* p){
	double val=0;
	val=Gaussian(x[0],p[0],p[1],p[2]);//Mean,sig, amp
	return val;
}
TF1* func_dgaus=new TF1("func_dgaus","dgaus",-5,5,8);
TF1* func_kgaus=new TF1("func_kgaus","sgaus",-5,5,3);
TF1* func_pigaus=new TF1("func_pigaus","sgaus",-5,5,3);
TF1* func_bgaus=new TF1("func_bgaus","sgaus",-5,5,3);
TChain* AssignChain(TString filename){
	TChain* chain=new TChain("tree");
	chain->Add(filename);
	chain->SetBranchAddress("trigflag",trigflag);
	chain->SetBranchAddress("trigpat",trigpat);
	chain->SetBranchAddress("bacnhits",&bacnhits);
	chain->SetBranchAddress("cbtof",cbtof);
	return chain;
}
void CountKPi(int run){
	TCanvas* c1 = new TCanvas("c1","c1",1800,1200);
	TString name_base = "Hodoscope0";
	fstream f;
	TString dir = "./params/";
	f.open((string)dir+"Run0"+to_string(run)+"Counts",fstream::out);
	TString filename = name_base+Form("%d.root",run);
	TChain* chain=new TChain("tree");
	chain=AssignChain((TString)filename);
	double tofm=-5,tofM=5;
	TString ht = Form("H_cbtof%d",run);
	TCut cut="trigpat==11";
	int nbin=1000;
	int scale = nbin/(tofM-tofm);
	int hn=run-5000;
	gfileOut->cd("Do Not Touch");
	Hist_tof[hn] = new TH1D(ht,ht,nbin,tofm,tofM);
	//	TH1D* Hist_bac = new TH1D(bt,bt,2,0,2);
	chain->Draw("cbtof>>"+ht,cut);
	func_dgaus->SetParLimits(0,-2,-1);
	func_dgaus->SetParLimits(1,0.05,0.2);
	//	func_dgaus->SetParLimits(2,em,eM);
	func_dgaus->SetParLimits(3,-0.5,1);
	func_dgaus->SetParLimits(4,0.05,0.2);
	//	func_dgaus->SetParLimits(5,kpr*em,kpr*eM);
	func_dgaus->SetParLimits(6,1,10);
	Hist_tof[hn]->Fit("func_dgaus","R");
	c1->SaveAs(dir+ht+".png");
	double par[8];
	for(int i=0;i<8;i++){
		par[i]=func_dgaus->GetParameter(i);
		cout<<par[i]<<endl;
	}
	func_kgaus->SetParameter(0,par[0]);
	func_kgaus->SetParameter(1,par[1]);
	func_kgaus->SetParameter(2,par[2]);
	
	func_pigaus->SetParameter(0,par[3]);
	func_pigaus->SetParameter(1,par[4]);
	func_pigaus->SetParameter(2,par[5]);

	func_bgaus->SetParameter(0,par[0]);
	func_bgaus->SetParameter(1,par[6]);
	func_bgaus->SetParameter(2,par[7]);

	func_kgaus->SetLineColor(4);
	func_pigaus->SetLineColor(3);
	func_bgaus->SetLineColor(1);
	func_bgaus->SetLineStyle(10);

	func_kgaus->Draw("Same");
	func_pigaus->Draw("Same");
	func_bgaus->Draw("Same");

	rn=run;
	KTime=par[0];
	KWidth=par[1];
	KCount=int(scale*par[2]);
	PiTime=par[3];
	PiWidth=par[4];
	PiCount=int(scale*par[5]);
	BGWidth=par[6];
	BGCount=int(scale*par[7]);




	gtree->Fill();
	f<<"TimeUnit: ps, "<<"Cut:"<< cut<<endl;
	f<<"KaonMeantime	"<<1000*par[0]<<endl;
	f<<"KaonSigma	"<<1000*par[1]<<endl;
	f<<"PionMeantime	"<<1000*par[3]<<endl;
	f<<"PionSigma	"<<1000*par[4]<<endl;
	f<<"BackgroundSigma	"<<1000*par[6]<<endl;
	f<<"KaonEntries	"<<int(scale*par[2])<<endl;
	f<<"PionEntries	"<<int(scale*par[5])<<endl;
	f<<"BackgroundEntries	"<<int(scale*par[7])<<endl;
	cout<<Form("Kaon: %d, Pion: %d",int(scale*par[2]),int(scale*par[5]))<<endl;
}
void Beam(){

	func_dgaus->SetParName(0,"KMeanTime[ns]");
	func_dgaus->SetParName(1,"KSpread(Sigma)[ns]");
	func_dgaus->SetParName(2,"KCount/100");
	func_dgaus->SetParName(3,"PiMeanTime");
	func_dgaus->SetParName(4,"PiSpread(Sigma)");
	func_dgaus->SetParName(5,"PiCount/100");
	func_dgaus->SetParName(6,"BGWidth");
	func_dgaus->SetParName(7,"BGCount/100");
	gStyle->SetOptFit(1111);
	cout<<"CountKPi(int run)"<<endl;
	//	for(int i=0;i<300;i++){
	gfileOut=new TFile("KPiCount.root","recreate");	
	gtree = new TTree("tree","tree");
	gtree->Branch("run",&rn,"run/I");
	gtree->Branch("KTime",&KTime,"KTime/D");
	gtree->Branch("KWidth",&KWidth,"KWidth/D");
	gtree->Branch("KCount",&KCount,"KCount/I");
	gtree->Branch("PiTime",&PiTime,"PiTime/D");
	gtree->Branch("PiWidth",&PiWidth,"PiWidth/D");
	gtree->Branch("PiCount",&PiCount,"PiCount/I");
	gtree->Branch("BGWidth",&BGWidth,"BGWidth/D");
	gtree->Branch("BGCount",&BGCount,"BGCount/I");
	gfileOut->mkdir("Do Not Touch");
	CountKPi(5764);
	gfileOut->Write();
	//	}
}
