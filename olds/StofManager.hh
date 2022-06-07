double m_pi = CPionMass()/1000;//MeV to GeV
double m_pi2= m_pi*m_pi;
double m_k = CKaonMass()/1000;//MeV to GeV
double m_k2= m_k*m_k;
double m_p = ProtonMass()/1000;//MeV to GeV
double m_p2= m_p*m_p;
double m2res = 8.82222e-02;
const int ntkurama=20;

double stof_to_m2(double stof, double p, double path){
	Double_t beta = path/stof/LightSpeed();//LightSpeed: mm / ns
	return p*p*(1.-beta*beta)/beta/beta;//stof: ns,  p: GeV/c, path: mm
}
double m2_to_stof(double m2, double p, double path){
	return path*sqrt(m2+p*p)/p/LightSpeed();
}


int ntKurama;
double pKurama[ntkurama];double qKurama[ntkurama];double m2[ntkurama];double tofsegKurama[ntkurama];double path[ntkurama];double stofs[ntkurama];double chisqrKurama[ntkurama];
int segn;
double stof_dt;
double chisqr;

class StofManager{
	private:
		TChain* chain;
		fstream f;
		TFile* m2file;
		TTree* tree;
		TH2D* h[1000];
		TH2D* h2[1000];
		TH1D* h1[100];
		int hcnt=0;
		int h1cnt=0;
		int h2cnt=0;
		double c2cut = 200;
		TString outfilename;
		TDirectory* dir[24];
	public:
		StofManager(){
			chain= new TChain("kurama"); 
			chain->SetBranchAddress("ntKurama",&ntKurama);
			chain->SetBranchAddress("path",path);
			chain->SetBranchAddress("pKurama",pKurama);
			chain->SetBranchAddress("qKurama",qKurama);
			chain->SetBranchAddress("tofsegKurama",tofsegKurama);
			chain->SetBranchAddress("chisqrKurama",chisqrKurama);
			chain->SetBranchAddress("stof",stofs);
			//			chain->SetBranchAddress("m2",m2);
		};
		void LoadKurama(TString filename){
			chain->Add(filename);
			cout<<"ent: "<<chain->GetEntries()<<endl;
			cout<<c2cut<<endl;
		}
		void MakeFile(TString filename);
		double DrawScatterPlot(int seg,double stof_offset,double m2_min,double m2_max);
		double DrawM2(int seg,double stof_offset,double m2_min,double m2_max);
};
double StofManager::DrawScatterPlot(int seg, double stof_offset, double m2_min, double m2_max){
	if(!gSystem->AccessPathName(outfilename)){
		gDirectory->cd(Form("Hist_m2%d",seg));
	}
	TString ht = Form("momentum_tofseg%d:protonm2",seg);		
//	h[hcnt]= new TH2D(ht,ht,1000,0,3,100,-5,5);
	h[hcnt]= new TH2D(ht,ht,1000,-1.5,1.5,1000,0,2);
	TString ht2 = Form("momentum_tofseg%d:charge*m2",seg);		
	double chi2 = 0;
	double cnt=0;
	double m2_mean=0;
	TF1* flin = new TF1("flin","flin",0,3,2);
	if(stof_offset==0){
//		TCut cut = Form("chisqrKurama<%f&&tofsegKurama==%d",c2cut,seg);
		TCut cut = Form("chisqrKurama<%f",c2cut);
//		chain->Draw("qKurama*m2:pKurama>>"+ht,cut,"colz");
		chain->Draw("pKurama:qKurama*m2>>"+ht,cut,"colz");
		chi2= 1E9;
	}
	else{
		int ent = chain->GetEntries();
		for(int i=0;i<ent;i++){
			Indicator(i,ent);
			chain->GetEntry(i);
			for(int j=0;j<ntKurama;j++){
//				if(chisqrKurama[j]<c2cut&&(int)tofsegKurama[j]==seg){
				if(chisqrKurama[j]<c2cut){
					double m2_calc = stof_to_m2(stofs[j]+stof_offset,pKurama[j],path[j]);
					if(m2_calc<m2_max&&m2_calc>m2_min&&pKurama[j]>0){
//						h[hcnt]->Fill(pKurama[j],qKurama[j]*m2_calc);
						h[hcnt]->Fill(qKurama[j]*m2_calc,pKurama[j]);
						chi2+=square(m_p-sqrt(m2_calc))/m2res;
						m2_mean+=m2_calc;
						cnt++;
						//							cout<<i<<" th evt m2 accept: "<<m2_calc<<endl;
					}
				}
			}
		}
		cout<<"Done!"<<endl;
		chi2=chi2/cnt;
		m2_mean=m2_mean/cnt;
		flin->SetParameter(0,m2_mean);
		flin->SetParLimits(0,m2_mean-0.01,m2_mean+0.01);
		flin->SetParLimits(1,-0.1,0.1);
		h[hcnt]->Fit("flin","RQ");
		cout<<"Slope : "<<flin->GetParameter(1)<<endl;
		cout<<"Protonmass_mean = "<<sqrt(m2_mean)<<endl;
		h[hcnt]->Draw("colz");
	}
	hcnt++;
	return chi2;
};
double StofManager::DrawM2(int seg, double stof_offset, double m2_min, double m2_max){
	TString ht = Form("charge*m2_seg%d",seg);		
	h1[h1cnt]= new TH1D(ht,ht,1000,m2_min,m2_max);
	double chi2 = 0;
	if(stof_offset==0){
		TCut cut = Form("chisqrKurama<%f&&tofsegKurama==%d",c2cut,seg);
		chain->Draw("qKurama*m2:pKurama>>"+ht,cut,"colz");
		chi2= 1E9;
	}
	else{
		int ent = chain->GetEntries();
		for(int i=0;i<ent;i++){
			Indicator(i,ent);
			chain->GetEntry(i);
			for(int j=0;j<ntKurama;j++){
				if(chisqrKurama[j]<c2cut&&(int)tofsegKurama[j]==seg){
					double m2_calc = stof_to_m2(stofs[j]+stof_offset,pKurama[j],path[j]);
					if(m2_calc<m2_max&&m2_calc>m2_min){
						//							cout<<i<<" th evt m2 accept: "<<m2_calc<<endl;
						h1[h1cnt]->Fill(qKurama[j]*m2_calc);
					}
				}
			}
		}
		cout<<"Done!"<<endl;
		h1[h1cnt]->Draw();
	}
	h1cnt++;
	return chi2;
};
void StofManager::MakeFile(TString filename){
	m2file = new TFile(filename);
	outfilename=filename;
	for(int i=0;i<24;i++){
		dir[i] = m2file->mkdir(Form("Hist_m2%d",i+1));
	}
};
