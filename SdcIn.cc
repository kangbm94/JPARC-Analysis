#include "ChamberMethod.hh"
ChamberManager c;
double sig = 1;
void SdcIn(){
	c.LoadFile("rootfiles/Callibration/AllSdcInTracking.root");
}
void Sdc1Drift(){
	TFile* fl = new TFile("Sdc1Drift.root","recreate");
	Sdc1dt2=80;
	double t1=-10,t2=Sdc1dt2;
	bool track = true;
	track = false;
	TGraphErrors* g[6];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/Sdc1Drift.txt");
	double t2_list[6] = {80,80,80,80,90,80};
	double t1_list[6] = {-10,-10,-10,-10,-10,-10};
	double t0_list[6] = {-5,-7,-8,-5,-7,-8};
	for(int i=1;i<7;++i){
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
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,Sdc1DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		g[i-1]->Fit("f_pol5","R");
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			par.push_back(f_pol5->GetParameter(j));
		}
		c.WriteDriftParameter(i,0,6,6,par);
		cout<<i<<" End!"<<endl;
	}
	fl->Save();
}
void Sdc2Drift(){
	TFile* fl = new TFile("Sdc2Drift.root","recreate");
	Sdc1dt2=120;
	double t1=-10,t2=Sdc2dt2;
	bool track = true;
	track = false;
	TGraphErrors* g[4];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/Sdc2Drift.txt");
	double t2_list[4] = {120,120,120,120};
	double t1_list[4] = {-5,-5,-5,-10};
	double t0_list[4] = {-5,-5,-5,-3};
	for(int i=1;i<5;++i){
		cout<<i<<" Start!"<<endl;
		c1->cd(i);
		t1=t1_list[i-1];
		t2=t2_list[i-1];
		double t0 = t0_list[i-1];
		f_pol5->SetRange(t1,t2);
		f_pol5->SetParLimits(0,t0,0);
//		f_pol5->SetParameter(0,t0/2);
		f_pol5->SetParameter(1,0.00);
		f_pol5->SetParameter(2,0.001);
		f_pol5->SetParameter(3,0.0);
		f_pol5->SetParameter(4,0.0);
		f_pol5->SetParameter(5,0.0);
		f_pol5->SetLineWidth(5);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(6+i,Sdc2DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		g[i-1]->Fit("f_pol5","R");
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			par.push_back(f_pol5->GetParameter(j));
		}
		c.WriteDriftParameter(6+i,0,6,6,par);
		cout<<6+i<<" End!"<<endl;
	}
	fl->Save();
}

/*
void Sdc2T0(){
	int rm = 930,rM = 980;

	TF1* f_lin = new TF1("f_lin","[0]*(x-[1])",rm,rM);	
	TFile* fl = new TFile("Sdc2T0.root","recreate");
	TF1* f_err = new TF1("f_err","Err_function",rm,rM,4);	
	TH1* h[800];
	double lsb=0.833;
	f_lin->SetParLimits(1,rm,rM);
	c.MakeParameterFile("Sdc2T0");
	TCanvas* can[4];
	TCanvas* cantot[4];
	double p0,p1,p2,p3;
	for(int j=0;j<4;j++){
		c.Sdc1Layer(6+j);
		cantot[j]=new TCanvas(Form("Cantot%d",j),Form("Cantot%d",j),1200,600);
		cout<<Form("Getting Drift Histo %d",j)<<endl;
		h[200*(j+1)-1]=c.GetHisto(j+7);
		int peak = h[200*(j+1)-1]->GetMaximum();
		//		f_err->SetRange(850,900);
		f_err->SetParLimits(0,peak/2,peak);
		f_err->SetParLimits(1,930,970);
		f_err->SetParLimits(2,0.1,100);
		h[200*(j+1)-1]->SetAxisRange(800,970);
		h[200*(j+1)-1]->Draw();
		fl->cd();
		h[200*(j+1)-1]->Write();
		h[200*(j+1)-1]->Fit("f_err","RM");
		p2 = f_err->GetParameter(2);
		f_err->FixParameter(2,p2);
		can[j]=new TCanvas(Form("Can%d",j),Form("Can%d",j),1200,600);
		can[j]->Divide(10,7);
		if(j<2){
			for(int i=0;i<Sdc2XNW;++i){
				can[j]->cd(i+1);
				cout<<Form("L%dWire %d",7+j,i+1)<<endl;
				h[200*j+i]=(TH1*)c.GetHisto(j+7,i+1);
				h[200*j+i]->SetAxisRange(800,970);
				int peakw = h[200*j+i]->GetMaximum();
				f_err->SetParLimits(0,peakw/2,peakw);
				//			h[200*j+i]->Fit("f_lin","RM");
				h[200*j+i]->Fit("f_err","RM");
				fl->cd();
				h[200*j+i]->Draw();
				fl->cd();
				h[200*j+i]->Write();
				double par = f_err->GetParameter(1)+sig*p2;
				c.WriteT0Parameter(7+j,i+1,lsb,-par);
			}
		}
		else {
			for(int i=0;i<Sdc2YNW;++i){
				can[j]->cd(i+1);
				cout<<Form("L%dWire %d",7+j,i+1)<<endl;
				h[200*j+i]=(TH1*)c.GetHisto(j+7,i+1);
				h[200*j+i]=(TH1*)c.GetHisto(j+7,i+1);
				h[200*j+i]->SetAxisRange(800,970);
				int peakw = h[200*j+i]->GetMaximum();
				f_err->SetParLimits(0,peakw/2,peakw);
				//			h[200*j+i]->Fit("f_lin","RM");
				h[200*j+i]->Fit("f_err","RM");
				fl->cd();
				h[200*j+i]->Draw();
				fl->cd();
				h[200*j+i]->Write();
				double par = f_err->GetParameter(1)+sig*p2;
				c.WriteT0Parameter(7+j,i+1,lsb,-par);
			}
		}
		cout<<7+j<<" Complete"<<endl;
	}
	fl->Save();
}
*/
