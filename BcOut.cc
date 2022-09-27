#include "ChamberManager.hh"
ChamberManager c;
double sig = 1;
void BcOut(){
	c.LoadFile("rootfiles/Callibration/AllBcOutTracking.root");
}
double DriftCurve(double* x, double* p){
	double dt = x[0]-p[0];
	double val=0;
	for(int i=1;i<6;i++){
		val+=pow(dt,i)*p[i];
	}
	return val;
}
void BcDrift(){
	TFile* fl = new TFile("BcDrift.root","recreate");
	Sdc1dt2=80;
	double t1=-10,t2=Sdc1dt2;
	bool track = true;
	TGraphErrors* g[12];
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(4,3);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/BcDrift.txt");
	double t2_list[12] = {40,40,40,40,40,40,
		40,40,40,40,40,40};
	double t1_list[12] = {-5,-5,-5,-5,-5,
		-5,-5,-5,-5,-5,-5,-5};
	double t0_list[12] = {-2,-4,-5,-5,-3,
		-5,-5,-5,-5,-5,-5,-5};
	for(int i=1;i<13;++i){
		if(i==7|i==8) continue;
		cout<<i<<" Start!"<<endl;
		c1->cd(i);
		t1=t1_list[i-1];
		t2=t2_list[i-1];
		double t0 = t0_list[i-1];
		f_pol5->SetParLimits(0,t0,0);
		f_pol5->SetRange(t1,t2);
		f_pol5->SetParameter(0,t0/2);
		f_pol5->SetParameter(1,0.0);
		f_pol5->SetParameter(2,0.001);
		f_pol5->SetParameter(3,0.0);
		f_pol5->SetParameter(4,0.0);
		f_pol5->SetParameter(5,0.0);
		f_pol5->SetLineWidth(5);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,BcDL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		g[i-1]->Fit("f_pol5","R");
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
		c.WriteDriftParameter(112+i,0,6,6,par);
		cout<<i<<" End!"<<endl;
	}
	fl->Save();
}
void Compare(){
	TString base_dir = "./rootfiles/Defocus/old/";
	int runnum = 5754;
	ChamberManager my;
	ChamberManager ws;
	my.LoadFile(base_dir+Form("run0%d_BcOutTracking.root",runnum));
	my.LoadBcOut();
	ws.LoadFile(base_dir+Form("run0%d_BcOutTrackingWs.root",runnum));
	ws.LoadBcOut();
	int ent = my.GetEntries();
	TH2D* histX = new TH2D("x_dif","x_dif",40,-250,250,100,-5,5);
	TH2D* histY = new TH2D("y_dif","y_dif",40,-250,250,100,-5,5);
	for(int i=0;i<ent;++i){
		my.GetEvent(i);
		ws.GetEvent(i);
		auto myTrack = my.GetTrack(0);
		auto wsTrack = ws.GetTrack(0);
		if(i%10000==0){
			cout<<i<<" th"<<endl;
		}
		for(int j=0;j<40;++j){
			double z = -250+j*(500./40.);
			TVector2 mypos = myTrack.GetPosition(z);
			TVector2 wspos = wsTrack.GetPosition(z);
			TVector2 dif = mypos-wspos;
			double x = dif.X(),y=dif.Y();
			histX->Fill(z,x);
			histY->Fill(z,y);
		}
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	histX->Draw("colz");
	c1->cd(2);
	histY->Draw("colz");



}
