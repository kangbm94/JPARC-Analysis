#include "StofManager.hh"
TString dir = "rootfiles/CH2/";
void StofTest(){
	cout<<"StofTest(double offset)"<<endl;
	cout<<"units are in ns"<<endl;
}
void StofTest(double offset){
	StofManager S;
	int run = 5477;
	TString filename = dir+Form("KuramaTracking_stofrecal0%d.root",run);
	TCanvas* c= new TCanvas("c1","c1",1200,600);
	S.LoadKurama(filename);
//	S.DrawScatterPlot(10,offset,-10,10);
	double c2=S.DrawScatterPlot(10,offset,-0.5,2);
	cout<<"chi2 = "<<c2<<endl;
//	S.DrawScatterPlot(15,0.000000001);
}
void StofAutomatic(int run){
	gStyle->SetOptFit(1111);
	TString filename = dir+Form("KuramaTracking_stof0%d.root",run);
	StofManager S;
	S.LoadKurama(filename);
	int seg = 24;
	for(int i=0;i<seg;i++){
		
	}
}
