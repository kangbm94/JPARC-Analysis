#include "ChamberMethod.hh"
ChamberManager c;
double sig = 1;//sig = 1.2 for Sdc4;
void SdcOut(){
	c.LoadFile("rootfiles/Callibration/AllSdcOutTracking.root");
	c.LoadSdcOut();
}
void Sdc3Drift(){
	TFile* fl = new TFile("Sdc3Drift.root","recreate");
//	Sdc1dt2=120;
	double t1=-10,t2=Sdc2dt2;
	bool track = true;
	track = true;
	TGraphErrors* g[4];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/Sdc3Drift.txt");
	double t2_list[4] = {100,100,110,100};
	double t1_list[4] = {-5,-3,-6,-5};
	double t0_list[4] = {-5,-3,-5,-3};
	for(int i=1;i<5;++i){
		cout<<i<<" Start!"<<endl;
		c1->cd(i);
		t1=t1_list[i-1];
		t2=t2_list[i-1];
		double t0 = t0_list[i-1];
		f_pol5->SetRange(t1+10,t2-10);
		f_pol5->SetParLimits(0,t0,0);
//		f_pol5->SetParameter(0,t0/2);
		f_pol5->SetParameter(1,0.00);
		f_pol5->SetParameter(2,0.001);
		f_pol5->SetParameter(3,0.0);
		f_pol5->SetParameter(4,0.0);
		f_pol5->SetParameter(5,0.0);
		f_pol5->SetLineWidth(5);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,Sdc3DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		g[i-1]->Fit("f_pol5","RM");
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			if(j==0){
				par.push_back(-f_pol5->GetParameter(j));
			}
			else{
				par.push_back(f_pol5->GetParameter(j));
			}
		}
		c.WriteDriftParameter(30+i,0,6,6,par);
		cout<<6+i<<" End!"<<endl;
	}
	fl->Save();
}
void Sdc4Drift(){
	TFile* fl = new TFile("Sdc4Drift.root","recreate");
//	Sdc1dt2=120;
	double t1=-10,t2=Sdc2dt2;
	bool track = true;
	track = true;
	TGraphErrors* g[4];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-2e-2,4e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/Sdc4Drift.txt");
	double t2_list[4] = {250,250,250,250};
	double t1_list[4] = {-10,-10,-10,-10};
	double t0_list[4] = {-5,-5,-10,-5};
	for(int i=1;i<5;++i){
		cout<<i<<" Start!"<<endl;
		c1->cd(i);
		t1=t1_list[i-1];
		t2=t2_list[i-1];
		double t0 = t0_list[i-1];
		f_pol5->SetRange(t1+15,t2-20);
		f_pol5->SetParLimits(0,t0,0);
//		f_pol5->SetParameter(0,t0/2);
		f_pol5->SetParameter(1,0.0);
		f_pol5->SetParameter(2,0.002);
		f_pol5->SetParameter(3,0);
		f_pol5->SetParameter(4,0.0);
		f_pol5->SetParameter(5,0.0);
		f_pol5->SetLineWidth(5);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(4+i,Sdc4DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		g[i-1]->Fit("f_pol5","RM");
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			if(j==0){
				par.push_back(-f_pol5->GetParameter(j));
			}
			else{
				par.push_back(f_pol5->GetParameter(j));
			}
		}
		c.WriteDriftParameter(34+i,0,6,6,par);
		cout<<6+i<<" End!"<<endl;
	}
	fl->Save();
}
