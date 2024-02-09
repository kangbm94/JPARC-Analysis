#include "ChamberMethod.hh"
ChamberManager c;
double sig = 1;//sig = 1.2 for Sdc4;
void SdcOut(){
	c.LoadFile("rootfiles/SdcOutTracking.root");
	c.LoadSdcOut();
}
void Sdc3Drift(){
	TFile* fl = new TFile("Sdc3Drift.root","recreate");
//	Sdc1dt2=120;
	double t1=-10,t2=Sdc2dt2;
	bool track = true;
	track = true;
	TGraphErrors* g[4];
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,2);
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	c.MakeParameterFile("./param/Sdc3Drift.txt");
	double t2_list[4] = {110,110,110,110};
	double t1_list[4] = {-10,-10,-10,-10};
	double t0_list[4] = {-5,-3,-5,-3};
	double slope[4];
	for(int i=0;i<4;++i){
		int layer = i+1;
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
			if(dat[0] -31 > -1 and dat[0] - 34 < 1){
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
		t1 = sl;
		t2 = t1 + 150;
		cout<<i<<" Start!"<<Form("t0 = %f, t1 = %f, t2 = %f",t0,t1,t2)<<endl;
		f_pol5->SetRange(t1+20,t2-40);
		cout<<sl<<endl;
//		f_pol5->FixParameter(0,-sl);
		auto param = params.at(i-1);
		for(int j=0;j<param.size();++j){
//			f_pol5->SetParameter(j,param.at(j));
		}
		f_pol5->SetParLimits(0,-sl-2,-sl+2);
		f_pol5->SetParLimits(1,1e-3,4e-2);	
		f_pol5->SetParLimits(2,0,3e-3);	
		f_pol5->SetParLimits(3,-3e-4,-0e-5);	
		f_pol5->SetParLimits(4,3e-8,3e-7);	
		f_pol5->SetParLimits(5,-3e-11,0);
		f_pol5->SetParameter(2,1e-3);
		f_pol5->SetParameter(3,-1e-5);
		f_pol5->SetParameter(4,5e-8);
		f_pol5->SetParameter(5,-2e-11);
		

		f_pol5->SetLineWidth(5);
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(i,Sdc3DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		g[i-1]->Draw("AP");
		f_pol5->SetRange(t1+10,t2-30);
		g[i-1]->Fit("f_pol5","R0M");
		f_pol5->SetRange(t1-10,t2+10);
		f_pol5->Draw("same");
		c1->Update();
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			if(j==0)par.push_back(f_pol5->GetParameter(j));
				//+ param.at(0));
			else par.push_back(f_pol5->GetParameter(j));
		}
		if(i==1)c.WriteComment(Form("##dt_length = %f ns",t2-t1));
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
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,2);
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TF1* f_pol5 = new TF1("f_pol5",DriftCurve,t1,t2,6);
	/*
	f_pol5->SetParLimits(1,-2e-2,4e-2);
	f_pol5->SetParLimits(2,-3e-3,3e-3);
	f_pol5->SetParLimits(3,-3e-4,3e-4);
	f_pol5->SetParLimits(4,-3e-6,3e-6);
	f_pol5->SetParLimits(5,-3e-8,3e-8);
*/
	c.MakeParameterFile("./param/Sdc4Drift.txt");
	double t2_list[4] = {280,280,280,280};
	double t1_list[4] = {-10,-10,-10,-10};
	double t0_list[4] = {-5,-5,-10,-5};
	double slope[4];
	for(int i=0;i<4;++i){
		int layer = i+5;
		c2->cd(i+1);
		auto h = c. GetDTTrackHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(peak-65,peak-10);
		f_errf->SetParameter(0,peakH);
		f_errf->SetParameter(2,5);
		f_errf->SetParLimits(3,0,peakH/10.);
		h->Fit("f_errf","QR");
		auto t0 = f_errf->GetParameter(1);
		f_errf->SetRange(t0-25,t0+25);
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
			if(dat[0] -35 > -1 and dat[0] - 38 < 1){
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
		double sl = t0 - slope[i-1];
		t1 = 0;
		t2 = t1 + 200;
		cout<<i<<" Start!"<<Form("t0 = %f, t1 = %f, t2 = %f",t0,t1,t2)<<endl;
		g[i-1] =	(TGraphErrors*)c.DLDTGraph(4+i,Sdc4DL,t1,t2, track);
		cout<<"Graph Loaded"<<endl;
		f_pol5->SetRange(t1,t2-20);
		auto param = params.at(i-1);
		for(int j=0;j<param.size();++j){
//			f_pol5->SetParameter(j,param.at(j));
		}
		f_pol5->FixParameter(0,0);
		f_pol5->SetParLimits(1,0.03,0.07);
		f_pol5->SetParLimits(2,-3e-3,3e-3);
		f_pol5->SetParLimits(3,-3e-5,3e-5);
		
		f_pol5->SetParameter(2,1e-4);
		f_pol5->SetParameter(3,-1e-6);
		f_pol5->SetParameter(4,1e-9);
		f_pol5->SetParameter(5,-1e-12);
		
		
//		f_pol5->SetParLimits(4,-3e8,3e8);
//		f_pol5->SetParLimits(5,-3e10,3e10);
/*
		f_pol5->SetParLimits(3,-3e6,3e6);
	*/	
		f_pol5->SetLineWidth(3);
		g[i-1]->Draw("AP");
		f_pol5->SetLineColor(kRed);
		g[i-1]->Fit("f_pol5","RL");
		f_pol5->SetRange(t1-10,t2+10);
//		f_pol5->SetLineColor(kBlue);
//		f_pol5->Draw("same");
		cout<<"Fit Graph"<<endl;
		vector<double> par;
		cout<<"Getting Parameters"<<endl;
		for(int j=0;j<6;j++){
			if(j==0)par.push_back(f_pol5->GetParameter(j));
			//+ param.at(0));
			else par.push_back(f_pol5->GetParameter(j));
		}
		if(i==1)c.WriteComment(Form("##dt_length = %f ns",t2-t1));
		c.WriteDriftParameter(34+i,0,6,6,par);
		cout<<6+i<<" End!"<<endl;
	}
	fl->Save();
}
void Sdc3T0(){
	c.MakeParameterFile("./param/Sdc3T0.txt");
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	vector<vector<double>>TotalParams;
	TString Layer[4] = {
	"SDC3-X1","SDC3-X2","SDC3-Y1","SDC3-Y2"
	};
	for(int i=0;i<4;++i){
		c1->cd(i+1);
		vector<double>params;
		int layer = i+1;
		auto* h = c.GetTDCHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(800,900);
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
		canv[i]->Divide(16,6);
		int layer = i+1;
		vector<vector<double>>WireParams;
		f_errf->SetRange(800,900);
		for(int w = 0; w<Sdc3NW;w++){
			vector<double>params;
			canv[i]->cd(w+1);
			int wire = w+1;
			auto* h = c.GetTDCHisto(layer,wire);
			int peakH = h->GetMaximum();
			auto peak = h->GetBinCenter(h->GetMaximumBin());
			f_errf->SetParameter(0,peakH);
			f_errf->SetParameter(2,-peakH/40);
			f_errf->SetParLimits(0,peakH/2,peakH);
			f_errf->SetParLimits(1,830,920);
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
//		gt=0;
		c.WriteComment("#"+Layer[i]);
		for(int w=0;w<WP.size();++w){
			auto par = WP.at(w);
			gr->AddPoint(w+1,par.at(1)-gt);
			vector<int>ID = {i+31,w+1};
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
void Sdc4T0(){
	c.MakeParameterFile("./param/Sdc4T0.txt");
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	vector<vector<double>>TotalParams;
	TString Layer[4] = {
	"SDC4-Y1","SDC4-Y2","SDC4-X1","SDC4-X2"
	};
	for(int i=0;i<4;++i){
		c1->cd(i+1);
		vector<double>params;
		int layer = i+5;
		auto* h = c.GetTDCHisto(layer);
		int peakH = h->GetMaximum();
		auto peak = h->GetBinCenter(h->GetMaximumBin());
		f_errf->SetRange(800,900);
		f_errf->SetParameter(0,peakH);
		f_errf->SetParameter(1,850);
		f_errf->SetParameter(2,-peakH/50);
//		f_errf->SetParLimits(0,peakH/2,peakH);
//		f_errf->SetParLimits(2,-peakH/20,0);
//		f_errf->SetParLimits(3,0,peakH/10.);
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
		if(i<2)canv[i]->Divide(10,7);
		else canv[i]->Divide(12,8);
		int layer = i+5;
		vector<vector<double>>WireParams;
		f_errf->SetRange(800,900);
		int Sdc4NW;
		if(i<2)Sdc4NW=Sdc4YNW;
		else Sdc4NW=Sdc4XNW;
		for(int w = 0; w<Sdc4NW;w++){
			vector<double>params;
			canv[i]->cd(w+1);
			int wire = w+1;
//			cout<<layer<<" , "<<wire<<endl;
			auto* h = c.GetTDCHisto(layer,wire);
			int peakH = h->GetMaximum();
			auto peak = h->GetBinCenter(h->GetMaximumBin());
			f_errf->SetParameter(0,peakH);
			f_errf->SetParameter(2,-peakH/40);
			f_errf->SetParLimits(0,peakH/2,peakH);
			f_errf->SetParLimits(1,800,900);
			f_errf->SetParLimits(2,-peakH/3,0);
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
	c3->Divide(2,2);
	for(int i=0;i<4;++i){
		c3->cd(i+1);
		TGraph* gr = new TGraph();
		auto WP = AllParams.at(i);
		double gt = TotalParams.at(i)[2];
		cout<<"sl = "<<gt<<endl;
		gt=0;
		c.WriteComment("###"+Layer[i]);
		c.WriteComment("#");
		for(int w=0;w<WP.size();++w){
			auto par = WP.at(w);
			gr->AddPoint(w+1,par.at(1)-1.28*gt);
			vector<int>ID = {i+35,w+1};
			vector<double> pars;
			pars.push_back(0.833);
			pars.push_back(-par.at(1)+1.28*gt);
			c.WriteParameter(ID,pars);
		}
		gr->Draw("AP");
		gr->SetMarkerStyle(23);
		gr->SetMarkerSize(2);
	}

}
