#include "KKMethod.hh"
#include "KKLegacy.hh"
#include "DetectorID.hh"
void KKAna(){
//		KM.LoadFile("./rootfiles/CH2/AllKKAna.root");
//		KM.LoadFile("./rootfiles/Production/AllKKAna.root");
//	KM.LoadFile("./rootfiles/CH2/DstKKAna05641.root");
	//		KM.LoadFile("../Other/E07_data/DstKKAna_CH2_Phase2.root");
	//	KM.LoadFile("~/WS_data/ch2target/run05666_KKAnaTest.root");
//	KM.LoadKK();
}

	
double dt = 0.025;


void CountKPlus(){
	gStyle->SetOptFit(111);
	int nb=100;
	double m1= 0.35, m2 = 0.6;
	TChain* chain = KM.GetPublicChain();
	double KaonM2 = (KaonMass*KaonMass)/1000/1000;
	TH1D* hist = new TH1D("Kaons","Kaons",nb,m1,m2);
	chain->Draw("sqrt(m2)>>Kaons","pKurama<1.4&&inside==1&&qKurama>0");
	TF1* func = new TF1("GausWithPol","GausWithPol",0.3,0.7,5);
	TF1* func2 = new TF1("GausWithPol2","GausWithPol",0.3,0.7,5);
	double peak = hist->GetMaximum();
	double base = (hist->GetBinContent(1)+hist->GetBinContent(nb-1))/2;
	func->SetParLimits(0,peak/2,peak*2);
	func->SetParLimits(1,0.45,0.55);
	func->SetParLimits(2,0.01,0.05);
	func->SetRange(0.35,0.6);
	func->SetParLimits(3,base/2,base*2);
	hist->Fit("GausWithPol","R");
	double kpeak = func->GetParameter(0);
	double kmass = func->GetParameter(1);
	double ksig = func->GetParameter(2);
	double p1 = func->GetParameter(3);
	double p2 = func->GetParameter(4);
	double kparset[5]={kpeak,kmass,ksig,0,0};
	double bparset[5]={0,kmass,ksig,p1,p2};
	double nKaon = kpeak*sqrt(2*Pi())*ksig*nb/(m2-m1);
	func->SetParameters(kparset);
	func->SetLineColor(kBlue);
	func->Draw("same");

	func2->SetParameters(bparset);
	func2->SetLineColor(kBlack);
	func2->Draw("same");
	cout<<"Effective Entries: "<<hist->GetEffectiveEntries()<<endl;
	cout<<"Kaon Counts: "<<nKaon<<endl;
	cout<<"Background Counts: "<<(p1+p2*(m1/2+m2/2-kmass))*nb<<endl;

}


double CountXi(double* par_,bool drawing,int i){
	double par[6]={0};
	TFile* file = new TFile("SelectedEvents.root","read");
//	TFile* file = new TFile("TPCInv.root","read");
	TCut cut = Form("cos(XiThetaCM/180*3.141592)<%f &&cos(XiThetaCM/180*3.141592)>%f",1-dt*i,1-dt*(i+1));
	TTree* tree = (TTree*)file->Get("tree");
	TH1D* h_xi = new TH1D("MissMass","MissMass",200,1,2);
	tree->Draw("XiM2>>MissMass",cut);
//	tree->Draw("MM>>MissMass","");
//	tree->Draw("MM>>MissMass","abs(InvMLd-1.12)<0.044");
	TH1D* h = KM.XiMinusFit(par,h_xi);
//	int nbin = h->GetXaxis()->GetNbinsX();
	int nbin = h->GetNbinsX();
	double bin_density = nbin/1;	
	TF1* XiSpectra = new TF1(Form("XiSpectra%d",i),"Gaussianf",1.2,1.4,3);
	TF1* XiBackground = new TF1(Form("XiBackground%d",i),"Gaussianf",1.2,1.4,3);
	TF1* GausWithBGf = new TF1(Form("GausWithBGf_%d",i),"GausWithBG",1.2,1.4,6);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetLineColor(kGreen);
	XiBackground->SetParameters(par[3],par[4],par[5]);
	XiBackground->SetLineColor(kBlue);
	GausWithBGf->SetParameters(par);
	GausWithBGf->SetLineColor(kRed);
	double XiCount = par[0]*par[2]*sqrt(2*Pi())*bin_density;
	int b1 = h_xi->FindBin(par[1]-3*par[2]); 
	int b2 = h_xi->FindBin(par[1]+3*par[2]);
	int Window = h_xi->Integral(b1,b2);
//	cout<<Form("Number of Xis: %f/%d",XiCount*bin_density,Window)<<endl;
	cout<<"Xi Mass: "<<1000*par[1]<<" MeV/c2"<<endl;
	cout<<"Xi Width: "<<1000*par[2]<<" MeV/c2"<<endl;
	for(int i=0;i<6;i++){
		par_[i]=par[i];
	}
	if(drawing){
		h->Draw();
		XiSpectra->Draw("same");
		XiBackground->Draw("same");
		GausWithBGf->SetRange(1.2,1.4);
		GausWithBGf->Draw("same");
	}
	return XiCount;
}


void CountXi(){
	double par[6];bool draw = true;
//	CountXi(par,draw);
}
double CountXiStar(double* par_,bool drawing,int i){
	double par[6];
	TFile* file = new TFile("SelectedEvents.root","read");
//	TFile* file = new TFile("TPCInv.root","read");
	TCut cut = Form("cos(XiThetaCM/180*3.141592)<%f &&cos(XiThetaCM/180*3.141592)>%f",1-dt*i,1-dt*(i+1));
	TTree* tree = (TTree*)file->Get("tree");
	TH1D* h_xi = new TH1D("MissMass","MissMass",200,1,2);
	
	tree->Draw("XiM2>>MissMass",cut);
//	tree->Draw("MM>>MissMass","abs(InvMLd-1.12)<0.044");
	TH1D* h = KM.XiStarFit(par,h_xi);
//	int nbin = h->GetXaxis()->GetNbinsX();
	int nbin = h->GetNbinsX();
	double bin_density = nbin/1;	
	TF1* XiSpectra = new TF1("XiSpectra","Gaussianf",1.4,1.6,3);
	TF1* XiBackground = new TF1("XiBackground","Gaussianf",1.4,1.6,3);
	TF1* GausWithBGf2 = new TF1(Form("GausWithBGf2_%d",i),"GausWithBG",1.2,1.4,6);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetLineColor(kGreen);
	XiBackground->SetParameters(par[3],par[4],par[5]);
	XiBackground->SetLineColor(kBlue);
	GausWithBGf2->SetParameters(par);
	GausWithBGf2->SetLineColor(kRed);
	double XiCount = par[0]*par[2]*sqrt(2*Pi())*bin_density;
//	cout<<"Number of XiStars: "<<XiCount*bin_density<<endl;
	cout<<"Xi StarMass: "<<1000*par[1]<<" MeV/c2"<<endl;
	cout<<"Xi StarWidth: "<<1000*par[2]<<" MeV/c2"<<endl;
	for(int i=0;i<6;i++){
		par_[i]=par[i];
	}
	if(drawing){
		h->Draw();
		XiSpectra->Draw("same");
		XiBackground->Draw("same");
		GausWithBGf2->SetRange(1.4,1.6);
		GausWithBGf2->Draw("same");
	}
	return XiCount;
}
void CountXiStar(){
	double par[6];bool draw = true;
//	CountXiStar(par,draw);
}

void DrawPred(int runnum){
	TString datafile,predfile;
	datafile = Form("./rootfiles/CH2/DstKKAna0%d.root",runnum);
	predfile = Form("PredictedDataReal0%d.root",runnum);
	TFile* prf = new TFile(predfile,"read");
	TTree* tree = (TTree*)prf->Get("tree");
	int pred,evnum;
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("pred",&pred);
	KM.LoadFile(datafile);
	KM.LoadKK();
	TChain* chain = KM.GetPublicChain();
	KKEvent Event(chain);
	int ent = chain->GetEntries();
	TH2D* HistMQ = new TH2D("HistMQ","HistMQ",300,-1,2,100,0,2);
	TH1D* HistMM = new TH1D("HistMM","HistMM",400,0,2);
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);	
	c1->Divide(2,1);
	for(int i=0;i<ent;++i){	
		Event.LoadEvent(i);
		tree->GetEntry(i);	
		if(i!=evnum){
			Event.Clear();
			continue;
		}
		if(pred==0){
			Event.Clear();
			continue;
		}
		Event.CutChiSqr();	
		Event.CutVertex();	
		Event.CutMomentum(1.1);	
		Event.CutCharge(-1);	
//		Event.CutM2();	
//		Event.CutM2(0.0,0.5);	
//		Event.CutM2(0,0.12);	
		for(int ikk=0;ikk<Event.Getnkk();++ikk){
			double p = Event.GetMomentum(ikk);	
			double q = Event.GetCharge(ikk);	
			double mm = Event.GetMissMass(ikk);	
			double m2 = Event.GetM2(ikk);	
			HistMQ->Fill(sqrt(m2)/q,p);	
			HistMM->Fill(mm);	
		}
		Event.Clear();
	}
	c1->cd(1);
	HistMQ->Draw("colz");
	c1->cd(2);
	HistMM->Draw();
}


void GetXiSpectra(int runnum){
	TString infile,outfile;
	if(runnum==0){
		infile= "./rootfiles/CH2/AllKKAna.root";
		outfile= "./SelectedEvents.root";
	}
	else{
		infile= "./rootfiles/CH2/run0"+to_string(runnum)+"_DstHSKKAna.root";
		outfile= "./KPTagged0"+to_string(runnum)+".root";
	}
	KM.LoadFile(infile);
	KM.LoadKK();
	TChain* chain = KM.GetPublicChain();
	KKEvent Event(chain);
	
	int ent = chain->GetEntries();
	int cnt =0;	
	TH1D* h = new TH1D("MissMass","MissMass",200,1,2);
	int XiEv,XiRun;
	double XiM2,XiP,XiU,XiV,XiTheta,XiThetaCM;

	TFile* file = new TFile(outfile,"RECREATE");
	TTree* tree = new TTree("tree","tree");
	tree->Branch("runnum",&XiRun,"runnum/I");	
	tree->Branch("evnum",&XiEv,"evnum/I");	
	tree->Branch("XiM2",&XiM2,"XiM2/D");	
	tree->Branch("XiMM",&XiM2,"XiMM/D");	
	tree->Branch("XiP",&XiP,"XiP/D");	
	tree->Branch("XiU",&XiU,"XiU/D");	
	tree->Branch("XiV",&XiV,"XiV/D");	
	tree->Branch("XiTheta",&XiTheta,"XiTheta/D");	
	tree->Branch("XiThetaCM",&XiThetaCM,"XiThetaCM/D");	

	for(int i=0;i<ent;++i){
		Indicator(i,ent);
		Event.LoadEvent(i);
//		Event.CutTrig(22);//23 = DPS;
		Event.CutChiSqr();	
		Event.CutVertex();	
		Event.CutMomentum();	
		Event.CutCharge();	
		Event.CutM2();
		for(int ikk=0;ikk<Event.Getnkk();++ikk){
			XiRun=Event.GetRunNum();
			XiEv=Event.GetEvNum();
			XiM2=Event.GetMissMass(ikk);
			XiP=Event.GetMomentum(ikk);
			XiU=Event.GetUKP(ikk);
			XiV=Event.GetVKP(ikk);
			XiTheta=Event.GetTheta(ikk);
			XiThetaCM=Event.GetThetaCM(ikk);
			cnt++;
			tree->Fill();
			h->Fill(XiM2);
			if(ikk){
				cout<<"Xi Multiplicity!"<<endl;
			}
		}
		Event.Clear();
	}
	file->Write();
	file->Close();
	cout<<cnt<<endl;
	h->Draw();
}
/*
void TriggerAnalyze(int runnum){
	TString infile,outfile;
	if(runnum==0){
		infile= "./rootfiles/CH2/AllKKAna.root";
		outfile= "./AllTriggerTagged.root";
	}
	else{
		infile= "./rootfiles/CH2/DstKKAna0"+to_string(runnum)+".root";
		outfile= "./TriggerTagged0"+to_string(runnum)+".root";
	}
	KM.LoadFile(infile);
	KM.LoadKK();
	TChain* chain = KM.GetPublicChain();
	KKEvent Event(chain);
	TVector3 Pos(9,0,-43);
	TVector3 RealSize(30.3,50,20.2);

	for(int i=0;i<ent;++i){
		Indicator(i,ent);
		Event.LoadEvent(i);
		Event.CutChiSqr();	
		Event.CutVertex();	
		Event.CutMomentum();	
		Event.CutCharge();	
		Event.CutM2();
		for(int ikk=0;ikk<Event.Getnkk();++ikk){
			XiEv=i;
			XiM2=Event.GetMissMass(ikk);
			XiP=Event.GetMomentum(ikk);
			XiU=Event.GetUKP(ikk);
			XiV=Event.GetVKP(ikk);
			XiTheta=Event.GetTheta(ikk);
			cnt++;
			tree->Fill();
			h->Fill(Event.GetMissMass(ikk));
			if(ikk){
				cout<<"Xi Multiplicity!"<<endl;
			}
		}
		Event.Clear();
	}




}
*/
void GetDiffCr(){
	double pxi[20][6],pxis[10][6];
	double nxi[20],nxis[20],nx[20];
	double nxie[20],nxise[20];
	double cth[20],cte[20];
	TCanvas* c1 = new TCanvas("c1","c1",1200,900);
	c1->Divide(5,4);
	int np = 20;
	int lim = 12;
	for(int i=0;i<np;++i){
		c1->cd(i+1);
		nxi[i]=CountXi(pxi[i],true,i);
		if(i<10) nxie[i]=1/sqrt(nxi[i]);
		else nxie[i]=0;
		cout<<nxie[i]<<endl;
		cth[i]=1-dt*(i+0.5);
		cte[i]=0;
		nx[i]=1;
	}
	TCanvas* c2 = new TCanvas("c2","c2",1200,900);
	c2->Divide(5,4);
	for(int i=0;i<np;++i){
		c2->cd(i+1);
		nxis[i]=CountXiStar(pxis[i],true,i);	
		if(i<lim)nxise[i]=1/sqrt(nxis[i]);
		else nxise[i]=0;
	}
	TCanvas* c3 = new TCanvas("c3","c3",1000,600);
	c3->Divide(2,1);
//	TGraphErrors* g1= new TGraphErrors(np,cth,nxi,cte,nxie);
//	TGraphErrors* g2= new TGraphErrors(np,cth,nxis,cte,nxise);
	TGraphErrors* g1= new TGraphErrors(10,cth,nx,cte,nxie);
	TGraphErrors* g2= new TGraphErrors(12,cth,nx,cte,nxise);
	c3->cd(1);
	g1->Draw("AP");
	g1->SetTitle("#Xi vs cos(#theta)");
	g1->GetXaxis()->SetTitle("cos(#theta)");
	g1->GetYaxis()->SetTitle("Yield Error Bar");
	g1->SetMarkerStyle(31);
	g1->SetMarkerSize(1);
	g1->GetYaxis()->SetRangeUser(0,2);
	g1->GetXaxis()->SetLimits(0.6,1);
	c3->cd(2);
	g2->Draw("AP");
	g2->SetTitle("#Xi* vs cos(#theta)");
	g2->GetXaxis()->SetTitle("cos(#theta)");
	g2->GetYaxis()->SetTitle("Yield Error Bar");
	g2->SetMarkerStyle(31);
	g2->SetMarkerSize(1);
	g2->GetYaxis()->SetRangeUser(0,2);
	g2->GetXaxis()->SetLimits(0.6,1);
}
