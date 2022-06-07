#include "ChamberManager.hh"
ChamberManager c;
double Err_function(double* x, double* p){
	double dx =(p[1] -x[0])/p[2];
	return p[0]*((erf(dx)+1)/2)+p[3];
}
double sig = 1;
void BcOut(){
	c.LoadFile("rootfiles/CH2/BcOutTracking05453.root");
}
double DriftCurve(double* x, double* p){
	double dt = x[0]-p[0];
	double val=0;
	for(int i=1;i<6;i++){
		val+=pow(dt,i)*p[i];
	}
	return val;
}

void BcT0(){
	int rm = 440,rM = 465;

	TF1* f_lin = new TF1("f_lin","[0]*(x-[1])",rm,rM);	
	TFile* fl = new TFile("BcT0.root","recreate");
	TF1* f_err = new TF1("f_err","Err_function",rm,rM,4);	
	TH1* h[2400];
	double lsb=0.833;
	f_lin->SetParLimits(1,rm,rM);
	c.MakeParameterFile("BcT0");
	TCanvas* can[12];
	TCanvas* cantot[12];
	double p0,p1,p2,p3;
	for(int j=0;j<12;j++){
		c.BcLayer(j);
		cantot[j]=new TCanvas(Form("Cantot%d",j),Form("Cantot%d",j),1200,600);
		h[200*(j+1)-1]=c.GetHisto(j+1);
//		h[200*(j+1)-1]=c.GetDTTrackHisto(j+1);

		int peak = h[200*(j+1)-1]->GetMaximum();
		//		f_err->SetRange(850,900);
		f_err->SetParLimits(0,peak/2,peak);
		f_err->SetParLimits(1,445,460);
		f_err->SetParLimits(2,0.1,100);
		h[200*(j+1)-1]->SetAxisRange(300,470);
		h[200*(j+1)-1]->Draw();
		fl->cd();
		h[200*(j+1)-1]->Write();
		h[200*(j+1)-1]->Fit("f_err","RM");
		p2 = f_err->GetParameter(2);
		f_err->FixParameter(2,p2);
		can[j]=new TCanvas(Form("Can%d",j),Form("Can%d",j),1200,600);
		can[j]->Divide(8,8);
		for(int i=0;i<BcNW;++i){
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
//			if(j==0) par+=4;
			c.WriteT0Parameter(Bcid+j,i+1,lsb,-par);
		}
	}
	fl->Save();
}

void BcDrift(){
	TFile* fl = new TFile("BcDrift.root","recreate");
	Bcdt2=90;
	double t1=0,t2=Bcdt2;
	bool track = true;
	track = false;
	TGraphErrors* g[12];
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(4,3);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->FixParameter(0,0);
//	f_pol5->SetParLimits(0,0,3);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("BcDrift");
	c.WriteTParameter(t2);
	double t1_list[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double t2_list[12] = {60,60,60,60,60,60,60,60,60,60,60,60};
	for(int i=1;i<13;++i){
		cout<<i<<" Start!"<<endl;
		c1->cd(i);
		t1=t1_list[i-1];
		t2=t2_list[i-1];
		f_pol5->SetRange(t1,t2);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,BcDL,t1,t2, track);
		g[i-1]->Draw("AP");
		f_pol5->FixParameter(0,0);
		g[i-1]->Fit("f_pol5","R");
		f_pol5->SetParLimits(0,-5,5);
		g[i-1]->Fit("f_pol5","R");
		double par[6];
		for(int j=0;j<6;j++){
			par[j] = f_pol5->GetParameter(j);
		}
		c.WriteDriftParameter(Bcid+i-1,0,6,6,par);
		cout<<i<<" End!"<<endl;
	}
	fl->Save();
}

void BcDrift(int dum){
	TFile* fl = new TFile("BcDrift.root","recreate");
	Bcdt2=90;
	bool track = true;
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(4,3);
	double t1=10,t2=Bcdt2;
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	TF1* fgaus = new TF1("fgaus","gaus",0,7);
	f_pol5->FixParameter(0,0);
//	f_pol5->SetParLimits(0,0,3);
	f_pol5->SetParLimits(1,-3e-1,3e-1);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("BcDrift_dtdl");
	c.WriteTParameter(t2);
	double t1_list[12] = {5,5,5,5,5,5,5,5,5,5,5,5};
	double t2_list[12] = {35,35,35,35,35,35,35,35,35,35,35,35};
	int bw=3;
	TH1* h1[1200];
	TH2D* h[12];
	double dtm[1200];
	double dlm[1200];
	TGraph* g;
	for(int i=1;i<13;++i){
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
			dtm[j] = hdum->GetBinCenter(j+bw/2);
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
		c.WriteDriftParameter(Bcid-1+i,0,6,6,par);
	}
	fl->Save();

}

