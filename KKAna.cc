#include "KKMethod.hh"
#include "KKLegacy.hh"
void KKAna(){
//		KM.LoadFile("./rootfiles/CH2/AllKKAna.root");
		KM.LoadFile("./rootfiles/Production/AllKKAna.root");
//	KM.LoadFile("./rootfiles/CH2/DstKKAna05641.root");
	//		KM.LoadFile("../Other/E07_data/DstKKAna_CH2_Phase2.root");
	//	KM.LoadFile("~/WS_data/ch2target/run05666_KKAnaTest.root");
	KM.LoadKK();
}



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


void CountXi(double* par_,bool drawing){
	double par[6];
	TH1D* h = KM.XiMinusFit(par);
//	int nbin = h->GetXaxis()->GetNbinsX();
	int nbin = h->GetNbinsX();
	cout<<nbin<<endl;
	int bin_density = nbin/1;	
	TF1* XiSpectra = new TF1("XiSpectra","Gaussianf",1.2,1.4,3);
	TF1* XiBackground = new TF1("XiBackground","Gaussianf",1.2,1.4,3);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetLineColor(kGreen);
	XiBackground->SetParameters(par[3],par[4],par[5]);
	XiBackground->SetLineColor(kBlue);
	GausWithBGf->SetParameters(par);
	GausWithBGf->SetLineColor(kRed);
	double XiCount = par[0]*par[2]*sqrt(2*Pi());
	cout<<"Number of Xis: "<<XiCount*bin_density<<endl;
	cout<<"Xi Mass: "<<1000*par[1]<<" MeV/c2"<<endl;
	cout<<"Xi Width: "<<1000*par[2]<<" MeV/c2"<<endl;
	for(int i=0;i<6;i++){
		par_[i]=par[i];
	}
	if(drawing){
		TCanvas* c1;
		c1= new TCanvas("XiCanv","XiCanv",1200,800);
		c1->cd();
		h->Draw();
		XiSpectra->Draw("same");
		XiBackground->Draw("same");
		GausWithBGf->Draw("same");
	}
}
void CountXi(){
	double par[6];bool draw = true;
	CountXi(par,draw);
}
void CountXiStar(double* par_,bool drawing){
	double par[6];
	TH1D* h = KM.XiStarFit(par);
//	int nbin = h->GetXaxis()->GetNbinsX();
	int nbin = h->GetNbinsX();
	cout<<nbin<<endl;
	int bin_density = nbin/1;	
	TF1* XiSpectra = new TF1("XiSpectra","Gaussianf",1.4,1.6,3);
	TF1* XiBackground = new TF1("XiBackground","Gaussianf",1.4,1.6,3);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetLineColor(kGreen);
	XiBackground->SetParameters(par[3],par[4],par[5]);
	XiBackground->SetLineColor(kBlue);
	GausWithBGf->SetParameters(par);
	GausWithBGf->SetLineColor(kRed);
	double XiCount = par[0]*par[2]*sqrt(2*Pi());
	cout<<"Number of XiStars: "<<XiCount*bin_density<<endl;
	cout<<"Xi StarMass: "<<1000*par[1]<<" MeV/c2"<<endl;
	cout<<"Xi StarWidth: "<<1000*par[2]<<" MeV/c2"<<endl;
	for(int i=0;i<6;i++){
		par_[i]=par[i];
	}
	if(drawing){
		TCanvas* c1;
		c1= new TCanvas("XiCanv","XiCanv",1200,800);
		c1->cd();
		h->Draw();
		XiSpectra->Draw("same");
		XiBackground->Draw("same");
		GausWithBGf->Draw("same");
	}
}
void CountXiStar(){
	double par[6];bool draw = true;
	CountXiStar(par,draw);
}


void GetXiSpectra(){
	TChain* chain = KM.GetPublicChain();
	KKEvent Event(chain);
	
	int ent = chain->GetEntries();
//	cout<<ent<<endl;
	int cnt =0;	
	TH1D* h = new TH1D("MissMass","MissMass",200,1,2);
	int XiEv;
	double XiM2,XiP,XiU,XiV,XiTheta;
	
	TFile* file = new TFile("Xi_CH2.root","RECREATE");
	TTree* tree = new TTree("tree","tree");
	tree->Branch("XiEv",&XiEv,"XiEv/I");	
	tree->Branch("XiM2",&XiM2,"XiM2/D");	
	tree->Branch("XiP",&XiP,"XiP/D");	
	tree->Branch("XiU",&XiU,"XiU/D");	
	tree->Branch("XiV",&XiV,"XiV/D");	
	tree->Branch("XiTheta",&XiTheta,"XiTheta/D");	

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
	file->Write();
	file->Close();
	cout<<cnt<<endl;
	h->Draw();
}
