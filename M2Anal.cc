#include "StofManager.hh"
TString dir = "rootfiles/CH2/";
int trigpat[32];
int trigflag[32];
TF1* fgaus = new TF1("fgaus","gaus",-50,50);
double prev_offset[24];




double xtofKurama[ntkurama];
double ytofKurama[ntkurama];
double utofKurama[ntkurama];
double vtofKurama[ntkurama];
void M2Anal(){
	cout<<"FitSTOF( int run)"<<endl;
	cout<<"STOFCorrection( int run)"<<endl;
	cout<<"StofOffset(int run)"<<endl;
	cout<<"Scat run: 5395"<<endl;
}






void M2Correction(int run){
	int seg =3;
	double c2cut = 200;
	double mom= 1.82;
	double m2min=-5,m2max=5,pmin=0,pmax=3;
	TString filename = dir+Form("KuramaTracking_stofcorrection0%d.root",run);
	TFile* file = new TFile(filename,"READ");
	TTree* tree = (TTree*)file->Get("kurama");
	tree->SetBranchAddress("ntKurama",&ntKurama);
	tree->SetBranchAddress("path",path);
	tree->SetBranchAddress("pKurama",pKurama);
	tree->SetBranchAddress("qKurama",qKurama);
	tree->SetBranchAddress("tofsegKurama",tofsegKurama);
	tree->SetBranchAddress("chisqrKurama",chisqrKurama);
	tree->SetBranchAddress("stof",stofs);
	//	TString filename = dir+Form("SdcOutTracking0%d.root",run);
}





void StofOffset(int run){
	double c2cut = 200;
	double m2min=-2,m2max=2,pmin=0,pmax=2.5;
	double mom=1.82;
	double path_cut = 3451;
	TString filename = dir+Form("KuramaTracking_stofcorrection0%d.root",run);
	TFile* file = new TFile(filename,"READ");
	TTree* tree = (TTree*)file->Get("kurama");
	tree->SetBranchAddress("ntKurama",&ntKurama);
	tree->SetBranchAddress("path",path);
	tree->SetBranchAddress("pKurama",pKurama);
	tree->SetBranchAddress("qKurama",qKurama);
	tree->SetBranchAddress("tofsegKurama",tofsegKurama);
	tree->SetBranchAddress("chisqrKurama",chisqrKurama);
	tree->SetBranchAddress("stof",stofs);
	TH1D* stof_hist = new TH1D("h","h",1000,-50,50);	
	TH1D* stof_calchist = new TH1D("h2","h2",1000,-50,50);	
	TH1D* res_hist = new TH1D("h3","h3",1000,-50,50);	
	TCanvas* c= new TCanvas("c1","c1",1200,600);
	c->Divide(3,1);
	int ent = tree->GetEntries();
	for(int i=0;i<ent;i++){
		Indicator(i,ent);
		tree->GetEntry(i);
		for(int j=0;j<ntKurama;j++){
			if(chisqrKurama[j]<c2cut&&path[j]>path_cut){
				double time = m2_to_stof(m_pi2,mom,path[j]); 	
				//		cout<<time<<endl;
				//		cout<<(int)tofsegKurama[0]<<endl;
				stof_calchist->Fill(time);	
				stof_hist->Fill(stofs[j]);	
				res_hist->Fill(time-stofs[j]);
			}
		}
	}
	double stof_offset = -2.0378;
	double mean,std,cmean,cstd,rmean,rstd;
	c->cd(1);
	mean=stof_hist->GetMean();
	std=stof_hist->GetStdDev();
	fgaus->SetRange(mean-std,mean+std);
	stof_hist->Fit("fgaus","R");
	stof_hist->Draw();
	mean=fgaus->GetParameter(1);
	cout<<"stof: "<<mean<<endl;
	c->cd(2);
	cmean=stof_calchist->GetMean();
	cstd=stof_calchist->GetStdDev();
	fgaus->SetRange(cmean-cstd,cmean+cstd);
	stof_calchist->Fit("fgaus","R");
	stof_calchist->Draw();
	cmean=fgaus->GetParameter(1);
	cout<<"ctof: "<<cmean<<endl;
	c->cd(3);
	rmean=res_hist->GetMean();
	rstd=res_hist->GetStdDev();
	fgaus->SetRange(rmean-rstd,rmean+rstd);
	res_hist->Fit("fgaus","R");
	res_hist->Draw();
	rmean=fgaus->GetParameter(1);
	cout<<"res_correction: "<<rmean+stof_offset<<endl;
}











void STOFCorrection(int run){
	double c2cut = 200;
	double m2min=-2,m2max=2,pmin=0,pmax=2.5;
	//	TString filename = dir+Form("SdcOutTracking0%d.root",run);
	//TString filename = dir+Form("KuramaTracking_stof0%d.root",run);
	TString filename = dir+Form("KuramaTracking_stofrerecal0%d.root",run);
	cout<<filename<<endl;
	TFile* file = new TFile(filename,"READ");
	TTree* tree = (TTree*)file->Get("kurama");
	tree->SetBranchAddress("ntKurama",&ntKurama);
	tree->SetBranchAddress("path",path);
	tree->SetBranchAddress("pKurama",pKurama);
	tree->SetBranchAddress("qKurama",qKurama);
	tree->SetBranchAddress("tofsegKurama",tofsegKurama);
	tree->SetBranchAddress("chisqrKurama",chisqrKurama);
	tree->SetBranchAddress("stof",stofs);
	/*
		 tree->SetBranchAddress("m2",m2);
	*/
	fstream f;
	f.open("stof_offset",fstream::out);
	int ent = tree->GetEntries();
	cout<<ent<<endl;
	TH1D* stof_hist[24];
	TH1D* stof_calchist[24];
	TH1D* res_hist[24];
	TCanvas* c= new TCanvas("c1","c1",1200,600);
	TCanvas* c2= new TCanvas("c2","c2",1200,600);
	TCanvas* c3= new TCanvas("c3","c3",1200,600);
	c->Divide(6,4);
	c2->Divide(6,4);
	c3->Divide(6,4);
	double mean[24];double std[24];
	double cmean[24];double cstd[24];
	double rmean[24];
	double mom = 1.82;
	for(int i=0;i<24;i++){
		TString title = Form("stof_%d",i+1);
		stof_hist[i]= new TH1D(title,title,1000,-50,50);
		TString title2 = Form("stofcalc_%d",i+1);
		stof_calchist[i]= new TH1D(title2,title2,1000,-50,50);
		TString title3 = Form("res_%d",i+1);
		res_hist[i]= new TH1D(title3,title3,1000,-50,50);
		/*
			 c->cd(i+1);
			 tree->Draw("stof>>"+title,Form("TofSeg==%d",i+1));
			 mean[i]=sth[i]->GetMean();
			 std[i]=sth[i]->GetStdDev();
			 fgaus->SetRange(mean[i]-sth[i],mean[i]+std[i]);
			 sth[i]->Fit("fgaus","QR");
			 mean[i]=fgaus->GetParameter(1);
			 cout<<title+" : "<<mean[i]<<endl;
			 f<<7<<"\t";
			 f<<i<<"\t";
			 f<<1<<"\t";
			 f<<2<<"\t";
			 f<<p[i]<<"\t";
			 f<<-1<<endl;
			 */
	}
	cout<<"Histos initialized"<<endl;
	for(int i=0;i<ent;i++){
		Indicator(i,ent);
		tree->GetEntry(i);
		for(int j=0;j<1;j++){
			if(chisqrKurama[j]<c2cut){
				double time = m2_to_stof(m_pi2,mom,path[j]); 	
				//		cout<<time<<endl;
				//		cout<<(int)tofsegKurama[0]<<endl;
				stof_calchist[(int)tofsegKurama[j]-1]->Fill(time);	
				stof_hist[(int)tofsegKurama[j]-1]->Fill(stofs[j]-prev_offset[(int)tofsegKurama[j]-1]);	
				res_hist[(int)tofsegKurama[j]-1]->Fill(time-stofs[j]+prev_offset[(int)tofsegKurama[j]-1]);
			}
		}
	}
	double ct;
	for(int i=0;i<24;i++){
		c->cd(i+1);
		mean[i]=stof_hist[i]->GetMean();
		std[i]=stof_hist[i]->GetStdDev();
		fgaus->SetRange(mean[i]-std[i],mean[i]+std[i]);
		stof_hist[i]->Draw();
		stof_hist[i]->Fit("fgaus","QR");
		mean[i]=stof_hist[i]->GetBinCenter(stof_hist[i]->GetMaximumBin());
		//	mean[i]=fgaus->GetParameter(1);
		cout<<Form("stof_%d : ",i+1)<<mean[i]<<endl;

		c2->cd(i+1);
		cmean[i]=stof_calchist[i]->GetMean();
		cstd[i]=stof_calchist[i]->GetStdDev();
		fgaus->SetRange(cmean[i]-cstd[i],cmean[i]+cstd[i]);
		stof_calchist[i]->Draw();
		stof_calchist[i]->Fit("fgaus","QR");
		//		cmean[i]=fgaus->GetParameter(1);
		cmean[i]=stof_calchist[i]->GetBinCenter(stof_calchist[i]->GetMaximumBin());
		cout<<Form("stofcalc_%d : ",i+1)<<cmean[i]<<endl;
		c3->cd(i+1);
		res_hist[i]->Draw();
		rmean[i]=res_hist[i]->GetBinCenter(res_hist[i]->GetMaximumBin());
		cout<<Form("res_%d : ",i+1)<<cmean[i]-mean[i]<<endl;
		ct=rmean[i];
		f<<7<<"\t";
		f<<i<<"\t";
		f<<1<<"\t";
		f<<2<<"\t";
		f<<ct<<"\t";
		f<<-1<<endl;
	}
}
/*
	 void FitSTOF(int run){
	 double c2cut = 200;
	 double m2min=-2,m2max=2,pmin=0,pmax=2.5;
//	TString filename = dir+Form("SdcOutTracking0%d.root",run);
TString filename = dir+Form("KuramaTracking_stofcorrection0%d.root",run);
cout<<filename<<endl;
TFile* file = new TFile(filename,"READ");
TTree* tree = (TTree*)file->Get("kurama");
tree->SetBranchAddress("path",path);
tree->SetBranchAddress("pKurama",pKurama);
tree->SetBranchAddress("qKurama",qKurama);
tree->SetBranchAddress("tofsegKurama",tofsegKurama);
//		tree->SetBranchAddress("stof",stofs);
//	double prev_offset=-2.0378;
fstream f;
f.open("stof_offset",fstream::out);
int ent = tree->GetEntries();
cout<<ent<<endl;
TString tt = "ht";
double mom = 1.82;
TH1D* sth[24];
TH1D* pth[24];
TCanvas* c= new TCanvas("c1","c1",1200,600);
c->Divide(6,4);
double mean[24];double std[24];
for(int i=0;i<24;i++){
TString title = Form("stof_%d",i+1);
sth[i]= new TH1D(title,title,1000,-50,50);
TString title2 = Form("stofcalc_%d",i+1);
pth[i]= new TH1D(title2,title2,1000,-50,50);
c->cd(i+1);
tree->Draw("stof>>"+title,Form("TofSeg==%d",i+1),Form("chisqrKurama<%f"),c2cut);
mean[i]=sth[i]->GetMean();
std[i]=sth[i]->GetStdDev();
fgaus->SetRange(mean[i]-std[i],mean[i]+std[i]);
sth[i]->Fit("fgaus","QR");
mean[i]=fgaus->GetParameter(1)-prev_par[i];
cout<<title+" : "<<mean[i]<<endl;
}
for(int i=0;i<ent;i++){
tree->GetEntry(i);

}
//	for(int 
}
void FitSTOF_allseg(int run){
run = 5395;
double c2cut = 200;
double m2min=-2,m2max=2,pmin=0,pmax=2.5;
//	TString filename = dir+Form("SdcOutTracking0%d.root",run);
TString filename = dir+Form("KuramaTracking_stof0%d.root",run);
cout<<filename<<endl;
TFile* file = new TFile(filename,"READ");
TTree* tree = (TTree*)file->Get("kurama");
tree->SetBranchAddress("path",path);
tree->SetBranchAddress("pKurama",pKurama);
tree->SetBranchAddress("qKurama",qKurama);
tree->SetBranchAddress("tofsegKurama",tofsegKurama);
double prev_offset=-2.0378;
fstream f;
f.open("stof_offset",fstream::out);
int ent = tree->GetEntries();
cout<<ent<<endl;
TString tt = "ht";
TH2D* h=new TH2D(tt,tt,100,m2min,m2max,100,pmin,pmax);
TH1D* sth[24];
TCanvas* c= new TCanvas("c1","c1",1200,600);
c->Divide(6,4);
double mean[24];double std[24];
for(int i=0;i<24;i++){
TString title = Form("stof_%d",i+1);
sth[i]= new TH1D(title,title,1000,-50,50);
c->cd(i+1);
tree->Draw("stof>>"+title,Form("TofSeg==%d",i+1));
tree->Draw("stof>>"+title,Form("TofSeg==%d",i+1));
mean[i]=sth[i]->GetMean();
std[i]=sth[i]->GetStdDev();
fgaus->SetRange(mean[i]-std[i],mean[i]+std[i]);
sth[i]->Fit("fgaus","QR");
mean[i]=fgaus->GetParameter(1);
cout<<title+" : "<<mean[i]<<endl;
}
}
*/
