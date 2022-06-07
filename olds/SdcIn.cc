#include "ChamberManager.hh"
ChamberManager c;
double Err_function(double* x, double* p){
	double dx =(p[1] -x[0])/p[2];
	return p[0]*((erf(dx)+1)/2)+p[3];
}
double sig = 1;
void SdcIn(){
//	c.LoadFile("rootfiles/CH2/SdcInTracking05453.root");
	c.LoadFile("rootfiles/CH2/BigSdcIn.root");
}
double DriftCurve(double* x, double* p){
	double dt = x[0]-p[0];
	double val=0;
	for(int i=1;i<6;i++){
		val+=pow(dt,i)*p[i];
	}
	return val;
}

void Sdc1T0(){
	int rm = 410,rM = 470;

	TF1* f_lin = new TF1("f_lin","[0]*(x-[1])",rm,rM);	
	TFile* fl = new TFile("Sdc1T0.root","recreate");
	TF1* f_err = new TF1("f_err","Err_function",rm,rM,4);	
	TH1* h[1200];
	double lsb=0.833;
	f_lin->SetParLimits(1,rm,rM);
	c.MakeParameterFile("Sdc1T0");
	TCanvas* can[6];
	TCanvas* cantot[6];
	double p0,p1,p2,p3;
	for(int j=0;j<6;j++){
		c.Sdc1Layer(j);
		cantot[j]=new TCanvas(Form("Cantot%d",j),Form("Cantot%d",j),1200,600);
		h[200*(j+1)-1]=c.GetHisto(j+1);
//		h[200*(j+1)-1]=c.GetDTTrackHisto(j+1);

		int peak = h[200*(j+1)-1]->GetMaximum();
		//		f_err->SetRange(850,900);
		f_err->SetParLimits(0,peak/2,peak);
		f_err->SetParLimits(1,400,470);
		f_err->SetParLimits(2,0.1,100);
		h[200*(j+1)-1]->SetAxisRange(200,470);
		h[200*(j+1)-1]->Draw();
		fl->cd();
		h[200*(j+1)-1]->Write();
		h[200*(j+1)-1]->Fit("f_err","RM");
		p2 = f_err->GetParameter(2);
		f_err->FixParameter(2,p2);
		can[j]=new TCanvas(Form("Can%d",j),Form("Can%d",j),1200,600);
		can[j]->Divide(8,8);
		for(int i=0;i<Sdc1NW;++i){
			can[j]->cd(i+1);
			h[200*j+i]=(TH1*)c.GetHisto(j+1,i+1);
			h[200*j+i]=(TH1*)c.GetHisto(j+1,i+1);
			h[200*j+i]->SetAxisRange(200,470);
			int peakw = h[200*j+i]->GetMaximum();
			f_err->SetParLimits(1,400,470);
			f_err->SetParLimits(0,peakw/2,peakw);
			//			h[200*j+i]->Fit("f_lin","RM");
			h[200*j+i]->Fit("f_err","QRM");
			fl->cd();
			h[200*j+i]->Draw();
			fl->cd();
			h[200*j+i]->Write();
			double par = f_err->GetParameter(1)+sig*p2;// 1.645 = 90 %
			if(j==0) par+=4;
			c.WriteT0Parameter(1+j,i+1,lsb,-par);
		}
	}
	fl->Save();
}

void Sdc1Drift(){
	TFile* fl = new TFile("Sdc1Drift.root","recreate");
	Sdc1dt2=90;
	double t1=0,t2=Sdc1dt2;
	bool track = true;
	track = false;
	TGraphErrors* g[6];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->FixParameter(0,0);
//	f_pol5->SetParLimits(0,0,3);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("Sdc1Drift");
	c.WriteTParameter(t2);
	double t2_list[6] = {80,70,80,80,80,80};
	for(int i=1;i<7;++i){
		cout<<i<<" Start!"<<endl;
		c1->cd(i);
		t2=t2_list[i-1];
		f_pol5->SetRange(t1,t2);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,Sdc1DL,t1,t2, track);
		g[i-1]->Draw("AP");
		f_pol5->FixParameter(0,0);
		g[i-1]->Fit("f_pol5","R");
		f_pol5->SetParLimits(0,-5,5);
		g[i-1]->Fit("f_pol5","R");
		double par[6];
		for(int j=0;j<6;j++){
			par[j] = f_pol5->GetParameter(j);
		}
		c.WriteDriftParameter(i,0,6,6,par);
		cout<<i<<" End!"<<endl;
	}
	fl->Save();
}

void Sdc1Drift(int dum){
	TFile* fl = new TFile("Sdc1Drift.root","recreate");
	Sdc1dt2=120;
	double t1=0,t2=Sdc1dt2;
	bool track = true;
	TH2D* h[6];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	TF1* fgaus = new TF1("fgaus","gaus",0,7);
	f_pol5->FixParameter(0,0);
//	f_pol5->SetParLimits(0,0,3);
	f_pol5->SetParLimits(1,-3e-1,3e-1);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("Sdc1Drift_dtdl");
	c.WriteTParameter(t2);
	double t1_list[6] = {0,0,0,0,0,0};
	double t2_list[6] = {70,70,70,70,70,70};
	int bw=3;
	TH1* h1[1000];
	double dtm[1000];
	double dlm[1000];
	TGraph* g;
	for(int i=1;i<7;++i){
		c1->cd(i);
		t1=t1_list[i-1];
		t2=t2_list[i-1];
		f_pol5->SetRange(t1,t2);
		h[i-1] =	(TH2D*)c.GetDTDLHisto(i);
		int nbinx=h[i-1] -> GetNbinsX();
		TH1* hdum = h[i-1]->ProjectionX();
		for(int j=0;j<nbinx-bw;j++){
			TString ht = Form("proj_%d",j);
			h1[j]= h[i-1]->ProjectionY(ht,j-1,j-1+bw);
			dlm[j] = h1[j]->GetBinCenter(h1[j]->GetMaximumBin());
			dtm[j] = hdum->GetBinCenter(j+1+bw/2);
		}
		g = new TGraph(nbinx-bw-1,dtm,dlm);
	//		h[i-1] ->SetAxisRange(0.,5.,"Y");
		h[i-1]->Draw("col");
		f_pol5->FixParameter(0,0);
		g->Draw("SAME");
//		h[i-1]->Fit("f_pol5","QR");
		g->Fit("f_pol5","QR");
		cout<<Form("Layer%d",i)<<endl;
		f_pol5->FixParameter(0,0);
//		h[i-1]->Fit("f_pol5","R");
		g->Fit("f_pol5","R");
		double par[6];
	
		for(int j=0;j<6;j++){
			par[j]=0;
			par[j] = f_pol5->GetParameter(j);
		}
		c.WriteDriftParameter(i,0,6,6,par);
	}
	fl->Save();

}

void Sdc2Drift(){
	TFile* fl = new TFile("Sdc2Drift.root","recreate");
	double t1=0,t2=Sdc2dt2;
	bool track = true;
	TGraphErrors* g[4];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->FixParameter(0,0);
	//	f_pol5->SetParLimits(0,-3,3);
	f_pol5->SetParLimits(1,-1e-1,1e-1);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-9,3e-9);
	c.MakeParameterFile("Sdc2Drift");
	c.WriteTParameter(t2);
	double t2_list[4] = {120,120,120,120};
	for(int i=1;i<5;++i){
		c1->cd(i);
		t2=t2_list[i-1];
		f_pol5->SetRange(t1,t2);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(6+i,Sdc2DL,t1,t2, track);
		g[i-1]->Draw("AP");
		f_pol5->FixParameter(0,0);
		g[i-1]->Fit("f_pol5","RM");
//		f_pol5->SetParLimits(0,-15,15);
		g[i-1]->Fit("f_pol5","R");
		double par[6];
		for(int j=0;j<6;j++){
			if(i!=0){
			par[j] = f_pol5->GetParameter(j);
			}
			else{
				par[j] = f_pol5->GetParameter(j)*pow(-1,j);
			}
		}
		c.WriteDriftParameter(6+i,0,6,6,par);
	}
	fl->Save();
}

void Sdc2Drift(int dum){
	TFile* fl = new TFile("Sdc2Drift.root","recreate");
	double t1=0,t2=Sdc2dt2;
	bool track = true;
	TH2D* h[4];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->FixParameter(0,0);
	//	f_pol5->SetParLimits(0,-3,3);
	f_pol5->SetParLimits(1,-1e-1,1e-1);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-9,3e-9);
	c.MakeParameterFile("Sdc2Drift_dtdl");
	c.WriteTParameter(t2);
	double t1_list[4] = {5,5,5,5};
	double t2_list[4] = {120,120,120,120};
	for(int i=1;i<5;++i){
		c1->cd(i);
		t2=t2_list[i-1];
		f_pol5->SetRange(t1,t2);
		h[i-1] =	(TH2D*)c.GetDTDLHisto(6+i);
		h[i-1]->Draw("col");
		f_pol5->FixParameter(0,0);
		h[i-1]->Fit("f_pol5","RM");
//		f_pol5->SetParLimits(0,-15,15);
		h[i-1]->Fit("f_pol5","R");
		double par[6];
		for(int j=0;j<6;j++){
			if(i!=0){
			par[j] = f_pol5->GetParameter(j);
			}
			else{
				par[j] = f_pol5->GetParameter(j)*pow(-1,j);
			}
		}
		c.WriteDriftParameter(6+i,0,6,6,par);
	}
	fl->Save();
}
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

