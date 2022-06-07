TChain* chain;
void KuramaStofCorrection(){
	cout<<"CorrectSegment(int seg, double cor)"<<endl;
	cout<<"LoadFile(TString filename)"<<endl;
	chain = new TChain("kurama");
	TString filename = "KuramaTracking_stof05477.root";
	chain->Add(filename);
}

void LoadFile(TString filename){
	chain->Add(filename);
}
double pol2sqr(double* x, double* p){
	return p[0]*pow((x[0]-p[1]),2);
}
TF1* pol2f = new TF1("pol2f","pol2sqr",0,2,2);
void CorrectSegment(int seg, double cor){
//	TString m2 = Form("pKurama:qKurama*pKurama*pKurama*(pow(299.792458*(stof+%f)/path,2)-1)",cor);
//	TString m2 = Form("pKurama:pKurama*pKurama*(pow(299.792458*(stof+%f)/path,2)-1)",cor);
	TString m2 = Form("pKurama:qKurama*pKurama*pKurama*(pow(299.792458*(stof+%f)/path,2)-1)",cor);
	int nbin1 = 1000,nbin2=100;
	double m2min = -1,m2max=2,pmin=0,pmax=2.5;
	TH2D* hist = new TH2D("ht","ht",nbin1,m2min,m2max,nbin2,pmin,pmax);
//	TH2D* hist = new TH2D("ht","ht",nbin2,pmin,pmax,nbin1,m2min,m2max);
	TCut cut = Form("tofsegKurama==%d",seg);
	chain->Draw(m2+">>ht",cut,"colz");	
//	pol2f->SetParLimits(1,0,0.4);
//	hist->Fit("pol2f");
}
