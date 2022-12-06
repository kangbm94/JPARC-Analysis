#ifndef KKManager_h 
#define KKManager_h 1
#include "ToFManager.hh"
#include "KKTrack.hh"
double Beampx,Beampy,Beampz,Kmpx,Kmpy,Kmpz,Kppx,Kppy,Kppz,Phipx,Phipy,Phipz,MissMassPhi,MissMassLambda;
int trigpat[32];
double Gaussianf(double* x,double* par){
	double amp1 = par[0],mean1=par[1],width1=par[2];
	return Gaussian(x[0],mean1,width1,amp1);
}
double GausWithBG(double* x, double* par){
	double amp1 = par[0],mean1=par[1],width1=par[2];
	double amp2 = par[3],mean2=par[4],width2=par[5];
	double signal = Gaussian(x[0],mean1,width1,amp1);
	double bg = Gaussian(x[0],mean2,width2,amp2);
	return signal + bg;
}
TF1* GausWithBGf = new TF1("GausWithBGf","GausWithBG",1,1.5,6);
TF1* GausWithBGf2 = new TF1("GausWithBGf2","GausWithBG",1,1.5,6);


double GausWithPol(double* x, double* par){
	double peak = par[0], mean = par[1], sigma = par[2];
	double gaus = Gaussian(x[0],mean,sigma,peak);
	int Order = 1;
	double poly=0;
	for(int i=0;i<Order+1;i++){
		poly+=pow(x[0]-mean,i)*par[3+i];
	}
	return gaus+poly;
}
double Mass[3] = {PionMass/1000,KaonMass/1000,ProtonMass/1000};
double MassWindow= 0.05;
void GausMassFit(TH1* hist, int particle,int charge,double* par){
	double mean = charge*Mass[particle];
	double r1 = mean-MassWindow, r2 = mean+MassWindow;
	double peak = hist->GetBinContent(mean); f_gaus->SetRange(r1,r2);
//	f_gaus->SetParLimits(0,0.5*peak,1.5*peak);
	f_gaus->SetParLimits(1,r1,r2);
	f_gaus->SetParLimits(2,0.01,0.1);
	hist->Fit("f_gaus","QR0");
	for(int i=0;i<3;i++){
		par[i]=f_gaus->GetParameter(i);
	}
	gPad->cd();
};


class KKManager: public ToFManager{
	public: 
		KKManager(){}
		void LoadKK();
		TH1D* XiMinusFit(double* par,TH1D* hist = NULL);			
		TH1D* XiStarFit(double* par,TH1D* hist = NULL);			
};

TH1D* KKManager::XiMinusFit(double* par,TH1D* hist = NULL){
	TH1D* h;
	if(!hist)h= (TH1D*)GetHistogram(6315);
	else h=hist;
	int peak= -1;
	int b1 = h->FindBin(XiMinusMass-0.1);
	int b2 = h->FindBin(XiMinusMass+0.1);
	for(int i=b1;i<b2;++i){ 
	int b = (int) h->GetBinContent(i);
	if(b>peak) peak=b;
	}
	GausWithBGf->SetRange(XiMinusMass-0.1,XiMinusMass+0.1);
	GausWithBGf->SetParLimits(0,peak/2,peak);
	GausWithBGf->SetParLimits(1,XiMinusMass-0.01,XiMinusMass+0.01);
	GausWithBGf->SetParLimits(2,1e-2,1e-1);
	GausWithBGf->SetParLimits(3,peak/10,peak/2);
	GausWithBGf->SetParLimits(4,XiMinusMass-0.1,XiMinusMass+0.1);
	GausWithBGf->SetParLimits(5,0.01,0.3);
	h->Fit("GausWithBGf","QR0");
	for(int i=0;i<6;++i){
		par[i]=GausWithBGf->GetParameter(i);
	}
	
	return (TH1D*)h;	
}


TH1D* KKManager::XiStarFit(double* par,TH1D* hist = NULL){
	TH1D* h;
	if(!hist)h= (TH1D*)GetHistogram(6315);
	else h=hist;
	int peak= -1;
	int b1 = h->FindBin(XiStarMass-0.1);
	int b2 = h->FindBin(XiStarMass+0.1);
	for(int i=b1;i<b2;++i){ 
	int b = (int) h->GetBinContent(i);
	if(b>peak) peak=b;
	}
	GausWithBGf2->SetRange(XiStarMass-0.1,XiStarMass+0.1);
	GausWithBGf2->SetParLimits(0,peak/3,peak);
	GausWithBGf2->SetParLimits(1,XiStarMass-0.01,XiStarMass+0.01);
	GausWithBGf2->SetParLimits(2,1e-3,2e-2);
	GausWithBGf2->SetParLimits(3,peak/10,peak/2);
	GausWithBGf2->SetParLimits(4,XiStarMass-0.1,XiStarMass+0.1);
	GausWithBGf2->SetParLimits(5,1e-2,0.3);
	h->Fit("GausWithBGf2","QR0");
	for(int i=0;i<6;++i){
		par[i]=GausWithBGf2->GetParameter(i);
	}
	return (TH1D*)h;	
}







void KKManager::LoadKK(){
	LoadChain("kk");
	DataChain->SetBranchAddress("ntK18",&ntK18);
	DataChain->SetBranchAddress("pK18",pK18);
	DataChain->SetBranchAddress("trigpat",trigpat);
	DataChain->SetBranchAddress("chisqrK18",chisqrK18);
	DataChain->SetBranchAddress("ntKurama",&ntKurama);
	DataChain->SetBranchAddress("chisqrKurama",chisqrKurama);
	DataChain->SetBranchAddress("nKm",&nKm);
	DataChain->SetBranchAddress("nKp",&nKp);
	DataChain->SetBranchAddress("nKK",&nKK);
	DataChain->SetBranchAddress("m2",m2);
	DataChain->SetBranchAddress("pKurama",pKurama);
	DataChain->SetBranchAddress("qKurama",qKurama);
	DataChain->SetBranchAddress("pOrg",pOrg);
	DataChain->SetBranchAddress("ukp",ukp);
	DataChain->SetBranchAddress("vkp",vkp);
	DataChain->SetBranchAddress("ukm",ukm);
	DataChain->SetBranchAddress("vkm",vkm);
	DataChain->SetBranchAddress("inside",inside);

}
KKManager KM;
#endif
