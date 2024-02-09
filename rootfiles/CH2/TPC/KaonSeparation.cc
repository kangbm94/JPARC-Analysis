#include "Math.hh"
double mpi = 139.570 / 1000; 
double mk = 493.677 / 1000; 
double mp = 938.272/1000;
double m2[5];
double pKurama[5];
double qKurama[5];
int inside[5];
double PI = acos(-1);
void SeparateKaon(int i);
void KaonSeparation(){
	SeparateKaon(5641);
}
double Gaus(double x, double amp,double mean, double sig){
	double PI = 3.141592;
	double X = (x - mean)/ sig;
	double N = sig * sqrt(2 * PI);
	return amp* 1./N * exp( -0.5*X*X);
//	return amp*exp( -0.5*X*X);
}
double Gauss(double* var, double* p){
	double x = var[0];
	double Amp = p[0];
	double Mean = p[1];
	double Sig = p[2];
	return Gaus(x,Amp,Mean,Sig);
}
double Spectra(double* x, double* p){
	double m = x[0];
	double PiAmp = p[0];
	double PiMass = p[1];
	double PiWidth = p[2];
	double KAmp = p[3];
	double KMass = p[4];
	double KWidth = p[5];
	double PAmp = p[6];
	double PMass = p[7];
	double PWidth = p[8];
//	double PiPar[3] ={ PiAmp,PiMass,PiWidth};
//	double KPar[3] ={ KAmp,KMass,KWidth};
//	double PPar[3] ={ PAmp,PMass,PWidth};

	double val = Gaus(m,PiAmp,PiMass,PiWidth)
	+ Gaus(m,KAmp,KMass,KWidth)
	+ Gaus(m,PAmp,PMass,PWidth);
	return  val;
}
TF1* fSpectra = new TF1("fSpectra","Spectra",0,3,9);
TF1* fGauss = new TF1("fGauss","Gauss",0,3,9);
void SeparateKaon(int i){
	TString filename = Form("dstfiles_ws/CH2/run0%d_DstKScatAna.root",i);
	TFile* file = new TFile(filename);
	TTree* tree = (TTree*)file->Get("kk");
	int ent = tree->GetEntries();
	TH1D* hist[100];
	TH2D* hist2D[100];
	TH1D* histPi[100];
	TH1D* histK[100];
	TH1D* histP[100];
	int pbin = 1;
	int nbin = 3000;
	for(int ih=0;ih<pbin;++ih){
		TString title = Form("Hist%d",ih);
		TString title2D = Form("Hist2D%d",ih);
		hist[ih] = new TH1D(title,title,nbin,-1,2);	
		hist2D[ih] = new TH2D(title2D,title2D,300,-1,2,250,0,2.5);	
		TString titlePi = Form("HistPion%d",ih);
		TString titleK = Form("HistKaon%d",ih);
		TString titleP = Form("HistProton%d",ih);
	}
	tree->SetBranchAddress("m2Org",m2);
	tree->SetBranchAddress("pKurama",pKurama);
	tree->SetBranchAddress("qKurama",qKurama);
	tree->SetBranchAddress("inside",inside);
	for(int iev=0;iev<ent;++iev){
		tree->GetEntry(iev);
		if(!inside[0]) continue;
		for(int ih=0;ih<pbin;++ih){
			hist[ih]->Fill(m2[0]*qKurama[0]);
			hist2D[ih]->Fill(m2[0]*qKurama[0],pKurama[0]);
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	hist2D[0]->Draw("colz");
	

	TCanvas* c2 = new TCanvas("c2","c2",1200,800);
	hist[0]->Draw("colz");
	
	int nh = hist[0]->GetEntries();
	nh *= 3./nbin;
	
	fGauss->SetParameter(0,nh /3);
	fGauss->SetParameter(1,mpi*mpi);
	fGauss->SetParameter(2,0.01);
	cout<<"NH = " << nh<<endl;	

	fSpectra->SetParameter(0,nh /10);
	fSpectra->SetParameter(1,mpi*mpi);
	fSpectra->SetParameter(2,0.01);
	fSpectra->SetParLimits(2,0.005,0.03);
	fSpectra->SetParameter(3,nh / 30);
	fSpectra->SetParameter(4,mk*mk);
	fSpectra->SetParameter(5,0.04);
	fSpectra->SetParLimits(5,0.03,0.10);
	fSpectra->SetParameter(6,nh );
	fSpectra->SetParameter(7,mp*mp);
	fSpectra->SetParameter(8,0.1);
	fSpectra->SetParLimits(8,0.07,0.12);

	hist[0]->Fit("fSpectra","R");
//	fSpectra->Draw("same");
//	fGauss->Draw("same");
}
