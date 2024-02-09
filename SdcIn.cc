#include "ChamberMethod.hh"
ChamberManager c;
double sig = 1;
void SdcIn(){
	c.LoadFile("rootfiles/CH2/run05000_SdcInTracking.root");
}
void Restore(){
	c.MakeParameterFile("./param/Sdc1T0.txt");
	ifstream pfile;
	pfile.open("DCTdcCalib_KBM");
	double dat[20];
	TString dat_line;
	int read_flag = ReadConfLine(pfile,dat,dat_line);
	vector<double>params;
	while(read_flag){
		if(read_flag==1){
			if(dat[0] -1 > -1 and dat[0] - 6 < 1){
				params.push_back(dat[2]);
			}
		}
		read_flag = ReadConfLine(pfile,dat,dat_line);
	}
	int p = 0;
	TString Layer[6] = {
	"SDC1-V1","SDC1-V2","SDC1-X1","SDC1-X2","SDC1-U1","SDC1-U2"
	};
	for(int i = 0; i<6;++i){
		c.WriteComment("#"+Layer[i]);
		for(int j=0;j<64;++j){
			vector<int>ID = {i+1,j+1};
			vector<double> pars;
			pars.push_back(0.833);
			pars.push_back(params.at(p));
			c.WriteParameter(ID,pars);
			p++;
		}
	}
}
void Sdc1T0(){
	c.MakeParameterFile("./param/Sdc1T0.txt");
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	vector<vector<double>>TotalParams;
	TString Layer[6] = {
	"SDC1-V1","SDC1-V2","SDC1-X1","SDC1-X2","SDC1-U1","SDC1-U2"
	};
	for(int i=0;i<6;++i){
		c1->cd(i+1);
		vector<double>params;
		int layer = i+1;
		auto* h = c.GetTDCHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(420,460);
		f_errf->SetParameter(0,peakH);
		f_errf->SetParameter(2,-peakH/20);
		f_errf->SetParLimits(0,peakH/2,peakH);
		f_errf->SetParLimits(2,-peakH/10,0);
		f_errf->SetParLimits(3,0,peakH/10.);
		h->Fit("f_errf","R");
		for(int j=0;j<4;++j){
			params.push_back(f_errf->GetParameter(j));
		}
		TotalParams.push_back(params);
	}
	
	vector<vector<vector<double>>>AllParams;
	TCanvas* canv[6];
	for(int i=0;i<6;++i){
		TString ct = Form("canv%d",i);
		canv[i] = new TCanvas(ct,ct,1400,800);
		canv[i]->Divide(16,4);
		int layer = i+1;
		vector<vector<double>>WireParams;
		f_errf->SetRange(420,460);
		for(int w = 0; w<Sdc1NW;w++){
			vector<double>params;
			canv[i]->cd(w+1);
			int wire = w+1;
			auto* h = c.GetTDCHisto(layer,wire);
			int peakH = h->GetMaximum();
			auto peak = h->GetBinCenter(h->GetMaximumBin());
			f_errf->SetParameter(0,peakH);
			f_errf->SetParameter(2,-peakH/20);
			f_errf->SetParLimits(0,peakH/2,peakH);
			f_errf->SetParLimits(1,435,455);
			f_errf->SetParLimits(2,-peakH/10,0);
			f_errf->SetParLimits(3,0,peakH/10.);
			h->Fit("f_errf","R");
			for(int j=0;j<4;++j){
				params.push_back(f_errf->GetParameter(j));
			}
			WireParams.push_back(params);	
		}
		AllParams.push_back(WireParams);
	}
	TCanvas* c3 = new TCanvas("c3","c3",100,100,800,600);
	c3->Divide(3,2);
	for(int i=0;i<6;++i){
		c3->cd(i+1);
		TGraph* gr = new TGraph();
		auto WP = AllParams.at(i);
		double gt = TotalParams.at(i)[2];
		cout<<"sl = "<<gt<<endl;
		gt=0;
		c.WriteComment("#"+Layer[i]);
		for(int w=0;w<WP.size();++w){
			auto par = WP.at(w);
			gr->AddPoint(w+1,par.at(1)-gt);
			vector<int>ID = {i+1,w+1};
			vector<double> pars;
			pars.push_back(0.833);
			pars.push_back(-par.at(1)+gt);
			c.WriteParameter(ID,pars);
		}
		gr->Draw("AP");
		gr->SetMarkerStyle(23);
		gr->SetMarkerSize(2);
	}

}
void Sdc1Drift(){
	TFile* fl = new TFile("Sdc1Drift.root","recreate");
	Sdc1dt2=80;
	double t1=-10,t2=Sdc1dt2;
	bool track = true;
//	track = false;
	TGraphErrors* g[6];
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/Sdc1Drift.txt");
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(3,2);
	double t0_list[6] ;
	double slope[6];
	for(int i=0;i<6;++i){
		int layer = i+1;
		c2->cd(i+1);
		auto h = c. GetDTTrackHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(peak-35,peak-5);
		f_errf->SetParameter(0,peakH);
//		f_errf->SetParameter(2,5);
		f_errf->SetParLimits(2,0.1,10);
		f_errf->SetParLimits(3,0,peakH/10.);
		h->Fit("f_errf","QR");
		auto t0 = f_errf->GetParameter(1);
		f_errf->SetRange(t0-7,t0+7);
		h->Fit("f_errf","QR");
		t0_list[i]=f_errf->GetParameter(1);
		slope[i]=f_errf->GetParameter(2);
	}
	ifstream pfile;
	pfile.open("DCDriftParam_KBM");
	double dat[20];
	TString dat_line;
	vector<vector<double>> params;
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	int read_flag = ReadConfLine(pfile,dat,dat_line);
	while(read_flag){
		if(read_flag==1){
			if(dat[0] -1 > -1 and dat[0] - 6 < 1){
				vector<double>param;
				for(int i=0;i<dat[3];++i){
					param.push_back(dat[4+i]);
				}
				params.push_back(param);
			}
		}
		read_flag = ReadConfLine(pfile,dat,dat_line);
	}
	for(int i=1;i<7;++i){
		c1->cd(i);
		double t0 = t0_list[i-1];
		double sl = slope[i-1];
		t1 = 0;
		t2 = t1 + 120;
		f_pol5->SetParLimits(0,t0,t0+sl);
		auto param = params.at(i-1);
		double dt0 = 0;
		dt0 = param.at(0);
		for(int j=0;j<param.size();++j){
//			f_pol5->SetParameter(j,param.at(j));
		}
		f_pol5->SetParLimits(1,2e-2,4e-2);
		f_pol5->SetParLimits(2,0,3e-3);
		f_pol5->SetParLimits(3,-3e-4,0);
		f_pol5->SetParLimits(4,0,3e-6);
		f_pol5->SetParLimits(5,-3e-10,0);
		f_pol5->FixParameter(0,0);
		f_pol5->SetParameter(1,5e-2);
		f_pol5->SetRange(t1,t2-30);
		f_pol5->SetLineWidth(4);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,Sdc1DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		cout<<"Fitting Graph"<<endl;
		f_pol5->SetLineColor(kRed);
		f_pol5->SetLineStyle(0);
		g[i-1]->Fit("f_pol5","R");
		f_pol5->SetRange(t1-10,t2+10);
		f_pol5->SetLineWidth(3);
		f_pol5->SetLineStyle(kDashed);
		f_pol5->SetLineColor(kBlue);
		f_pol5->Draw("same");
		c1->Update();
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			if(j==0)par.push_back(f_pol5->GetParameter(j));
			//	+ params.at(i-1)[0]);
			else par.push_back(f_pol5->GetParameter(j));
		}
		if(i==1)c.WriteComment(Form("##dt_length = %f ns",t2-t1));
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
//	track = false;
	TGraphErrors* g[4];
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	f_pol5->SetParLimits(1,-3e-2,3e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
	c.MakeParameterFile("./param/Sdc2Drift.txt");
	double t0_list[4];
	double slope[4];
	double t2_list[4] = {120,120,120,120};
	double t1_list[4] = {-5,-5,-5,-10};
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,2);
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	for(int i=0;i<4;++i){
		int layer = i+7;
		c2->cd(i+1);
		auto h = c. GetDTTrackHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(peak-35,peak-5);
		f_errf->SetParameter(0,peakH);
		f_errf->SetParameter(2,5);
		f_errf->SetParLimits(3,0,peakH/10.);
		h->Fit("f_errf","QR");
		auto t0 = f_errf->GetParameter(1);
		f_errf->SetRange(t0-15,t0+15);
		h->Fit("f_errf","QR");
		t0_list[i]=f_errf->GetParameter(1);
		slope[i]=f_errf->GetParameter(2);
	}



	ifstream pfile;
	pfile.open("DCDriftParam_KBM");
	double dat[20];
	TString dat_line;
	vector<vector<double>> params;
	int read_flag = ReadConfLine(pfile,dat,dat_line);
	while(read_flag){
		if(read_flag==1){
			if(dat[0] -7 > -1 and dat[0] - 10 < 1){
				vector<double>param;
				for(int i=0;i<dat[3];++i){
					param.push_back(dat[4+i]);
				}
				params.push_back(param);
			}
		}
		read_flag = ReadConfLine(pfile,dat,dat_line);
	}
	for(int i=1;i<5;++i){
		c1->cd(i);
		double t0 = t0_list[i-1];
		double sl = slope[i-1];
		cout<<sl<<endl;
		t1 = 0 ;
		t2 = t1 + 180;
		f_pol5->SetRange(t1+20,t2);
		f_pol5->SetParLimits(1,3.5e-2,5e-2);
		f_pol5->SetParLimits(2,0,3e-3);
		f_pol5->SetParLimits(3,-1e-5,0);
		f_pol5->SetParLimits(4,0,1e-7);
		f_pol5->SetParLimits(5,-3e-10,-3e-11);
		f_pol5->FixParameter(0,0);
		f_pol5->SetParameter(2,1e-3);
		f_pol5->SetParameter(3,-1e-6);
		f_pol5->SetParameter(4,1e-7);
		f_pol5->SetParameter(5,-1e-10);
		f_pol5->SetRange(t1,t2-30);
		cout<<i<<" Start!"<<Form("t0 = %f, t1 = %f, t2 = %f",t0,t1,t2)<<endl;
//		f_pol5->SetParLimits(0,t0,0);
		auto param = params.at(i-1);

//		for(int j=0;j<param.size();++j){
//			f_pol5->SetParameter(j,param.at(j));
//		}
//		f_pol5->SetParameter(0,t0/2);
	/*
		f_pol5->SetParameter(1,0.00);
		f_pol5->SetParameter(2,0.001);
		f_pol5->SetParameter(3,0.0);
		f_pol5->SetParameter(4,0.0);
		f_pol5->SetParameter(5,0.0);
		f_pol5->SetLineWidth(5);
*/
		f_pol5->SetLineWidth(5);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(6+i,Sdc2DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		g[i-1]->Fit("f_pol5","R");
//		f_pol5->Draw("same");
//		f_pol5->Draw("");
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
//			if(j==0)par.push_back(f_pol5->GetParameter(j) + params.at(i-1)[0]);
//			else par.push_back(f_pol5->GetParameter(j));
			par.push_back(f_pol5->GetParameter(j));
		}
		if(i==1)c.WriteComment(Form("##dt_length = %f ns",t2-t1));
		c.WriteDriftParameter(6+i,0,6,6,par);
		cout<<6+i<<" End!"<<endl;
	}
	fl->Save();
}
void Sdc2T0(){
	c.MakeParameterFile("./param/Sdc2T0.txt");
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	vector<vector<double>>TotalParams;
	TString Layer[4] = {
	"SDC2-X1","SDC2-X2","SDC2-Y1","SDC2-Y2"
	};
	for(int i=0;i<4;++i){
		c1->cd(i+1);
		vector<double>params;
		int layer = i+7;
		auto* h = c.GetTDCHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(920,980);
		f_errf->SetParameter(0,peakH);
		f_errf->SetParameter(2,-peakH/40);
		f_errf->SetParLimits(0,peakH/2,peakH);
		f_errf->SetParLimits(2,-peakH/20,0);
		f_errf->SetParLimits(3,0,peakH/10.);
		h->Fit("f_errf","R");
		for(int j=0;j<4;++j){
			params.push_back(f_errf->GetParameter(j));
		}
		TotalParams.push_back(params);
	}
	
	vector<vector<vector<double>>>AllParams;
	TCanvas* canv[4];
	for(int i=0;i<4;++i){
		TString ct = Form("canv%d",i);
		canv[i] = new TCanvas(ct,ct,1400,800);
		canv[i]->Divide(10,7);
		int layer = i+7;
		vector<vector<double>>WireParams;
		f_errf->SetRange(920,980);
		int Sdc2NW;
		if(i<2) Sdc2NW =  Sdc2XNW;
		else Sdc2NW =  Sdc2YNW;
		for(int w = 0; w<Sdc2NW;w++){
			vector<double>params;
			canv[i]->cd(w+1);
			int wire = w+1;
			auto* h = c.GetTDCHisto(layer,wire);
			int peakH = h->GetMaximum();
			auto peak = h->GetBinCenter(h->GetMaximumBin());
			f_errf->SetParameter(0,peakH);
			f_errf->SetParameter(2,-peakH/40);
			f_errf->SetParLimits(0,peakH/2,peakH);
			f_errf->SetParLimits(1,920,970);
			f_errf->SetParLimits(2,-peakH/3,0);
			f_errf->SetParLimits(3,0,peakH/10.);
			h->Fit("f_errf","QR");
			for(int j=0;j<4;++j){
				params.push_back(f_errf->GetParameter(j));
			}
			WireParams.push_back(params);	
		}
		AllParams.push_back(WireParams);
	}
	TCanvas* c3 = new TCanvas("c3","c3",100,100,800,600);
	c3->Divide(2,2);
	for(int i=0;i<4;++i){
		c3->cd(i+1);
		TGraph* gr = new TGraph();
		auto WP = AllParams.at(i);
		double gt = TotalParams.at(i)[2];
		cout<<"sl = "<<gt<<endl;
		gt=0;
		c.WriteComment("#"+Layer[i]);
		for(int w=0;w<WP.size();++w){
			auto par = WP.at(w);
			gr->AddPoint(w+1,par.at(1)-gt);
			vector<int>ID = {i+7,w+1};
			vector<double> pars;
			pars.push_back(0.833);
			pars.push_back(-par.at(1)+gt);
			c.WriteParameter(ID,pars);
		}
		gr->Draw("AP");
		gr->SetMarkerStyle(23);
		gr->SetMarkerSize(2);
	}

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
