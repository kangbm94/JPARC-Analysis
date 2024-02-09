TF1* fgaus = new TF1("fgaus","gaus",0,2);
TString nocor = "NoCor/";
//TString nocor = "Cor2nd/Resolution_Fixed/";
TString cor2 = "Cor2nd/";
double nw = 10;
double GetMomentum(TString run){
	int n = run.Length();
	double mom = atof(&run[n-3]);
	return mom;
}
void MomResol(){
	TString RunList[12] = {"pm300","pm400","pm500","pm800","pi400","pi500","pi600","pi800","p400","p500","p600","p800" };
//	TString RunList[12] = {"pm300","pm400","pm500","pm800","p400","p500","p600","p800","pi400","pi500","pi600","pi800" };
	TH1D* HistMom[24];
	TH1D* HistHSMom[24];
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(4,3);
	TCanvas* c2 = new TCanvas("c2","c2",50,50,800,800);
	c2->Divide(4,3);
	double dum[4] = {0};
	double mom_pmnocor[4];
	double off_pmnocor[4];
	double res_pmnocor[4];
	double mom_pmcor[4];
	double off_pmcor[4];
	double res_pmcor[4];
	double mom_HSpmnocor[4];
	double off_HSpmnocor[4];
	double res_HSpmnocor[4];
	double mom_HSpmcor[4];
	double off_HSpmcor[4];
	double res_HSpmcor[4];
	double mom_pinocor[4];
	double off_pinocor[4];
	double res_pinocor[4];
	double mom_picor[4];
	double off_picor[4];
	double res_picor[4];
	double mom_HSpinocor[4];
	double off_HSpinocor[4];
	double res_HSpinocor[4];
	double mom_HSpicor[4];
	double off_HSpicor[4];
	double res_HSpicor[4];
	for(int i=0;i<12;++i){
		TFile* file = new TFile(nocor+RunList[i]+".root");
		TTree* tree = (TTree*)file->Get("tpc");
		double mom = 0.001*GetMomentum(RunList[i]);
		TString title = RunList[i]+"hist";
		HistMom[i] = new TH1D(title,title,100,0.7 * mom, 1.3 * mom);
		TString cut = "helix_flag == 100 && ntTpc==1&&nhtrack>20&&abs(hypot(helix_cx,helix_cy)-helix_r)<10";
		c1->cd(i+1);
		cout<<RunList[i]<<endl;
		tree->Draw("mom0>>"+title,cut);
		double w = HistMom[i]->GetBinWidth(1);
		double peak = HistMom[i]->GetBinCenter(HistMom[i]->GetMaximumBin());
		fgaus->SetRange(peak - nw*w, peak + nw *w);
		HistMom[i]->Fit("fgaus","R");
		cout<<RunList[i]<<endl;
		double mom_fit = fgaus->GetParameter(1);
		double mom_res = fgaus->GetParameter(2);
		if(i<4){
			mom_pmnocor[i] = mom;
			off_pmnocor[i] = mom_fit - mom;
			res_pmnocor[i] = mom_res/mom;
		}
		else if( i < 8){
			mom_pinocor[i-4] = mom;
			off_pinocor[i-4] = mom_fit - mom;
			res_pinocor[i-4] = mom_res/mom;
		}
	}
	for(int i=0;i<12;++i){
		TFile* file = new TFile(cor2+RunList[i]+".root");
		TTree* tree = (TTree*)file->Get("tpc");
		double mom = 0.001*GetMomentum(RunList[i]);
		TString title = RunList[i]+"hist";
		HistMom[12+i] = new TH1D(title,title,100,0.7 * mom, 1.3 * mom);
		TString cut = "helix_flag == 100 && ntTpc==1&&nhtrack>20&&abs(hypot(helix_cx,helix_cy)-helix_r)<10";
		c2->cd(i+1);
		cout<<RunList[i]<<endl;
		tree->Draw("mom0>>"+title,cut);
		double w = HistMom[12+i]->GetBinWidth(1);
		double peak = HistMom[12+i]->GetBinCenter(HistMom[12+i]->GetMaximumBin());
		fgaus->SetRange(peak - nw*w, peak + nw *w);
		HistMom[12+i]->Fit("fgaus","R");
		cout<<RunList[i]<<endl;
		double mom_fit = fgaus->GetParameter(1);
		double mom_res = fgaus->GetParameter(2);
		if(i<4){
			mom_pmcor[i] = mom;
			off_pmcor[i] = mom_fit - mom;
			res_pmcor[i] = mom_res/mom;
		}
		else if( i < 8){
			mom_picor[i-4] = mom;
			off_picor[i-4] = mom_fit - mom;
			res_picor[i-4] = mom_res/mom;
		}
	}
		fgaus->SetLineColor(kBlue);
	for(int i=0;i<12;++i){
		TFile* file = new TFile(nocor+"HS"+RunList[i]+".root");
		TTree* tree = (TTree*)file->Get("tpc");
		double mom = 0.001*GetMomentum(RunList[i]);
		TString title = "HS"+RunList[i]+"hist";
		HistHSMom[i] = new TH1D(title,title,100,0.7 * mom, 1.3 * mom);
		TString cut = "helix_flag == 200 && ntTpc==1&&nhtrack>20&&abs(hypot(helix_cx,helix_cy)-helix_r)<10";
		c1->cd(i+1);
		cout<<RunList[i]<<endl;
		tree->Draw("mom0>>"+title,cut,"SAME");
		HistHSMom[i]->SetLineColor(kRed);
		double w = HistHSMom[i]->GetBinWidth(1);
		double peak = HistHSMom[i]->GetBinCenter(HistHSMom[i]->GetMaximumBin());
		fgaus->SetRange(peak - nw*w, peak + nw *w);
		HistHSMom[i]->Fit("fgaus","R");
		cout<<RunList[i]<<endl;
		double mom_fit = fgaus->GetParameter(1);
		double mom_res = fgaus->GetParameter(2);
		if(i<4){
			mom_HSpmnocor[i] = mom;
			off_HSpmnocor[i] = mom_fit - mom;
			res_HSpmnocor[i] = mom_res/mom;
		}
		else if( i < 8){
			mom_HSpinocor[i-4] = mom;
			off_HSpinocor[i-4] = mom_fit - mom;
			res_HSpinocor[i-4] = mom_res/mom;
		}
	}
	for(int i=0;i<12;++i){
		TFile* file = new TFile(cor2+"HS"+RunList[i]+".root");
		TTree* tree = (TTree*)file->Get("tpc");
		double mom = 0.001*GetMomentum(RunList[i]);
		TString title = RunList[i]+"hist";
		HistHSMom[12+i] = new TH1D(title,title,100,0.7 * mom, 1.3 * mom);
		TString cut = "helix_flag == 200 && ntTpc==1&&nhtrack>20&&abs(hypot(helix_cx,helix_cy)-helix_r)<10";
		c2->cd(i+1);
		cout<<RunList[i]<<endl;
		tree->Draw("mom0>>"+title,cut,"SAME");
		HistHSMom[12+i]->SetLineColor(kRed);
		double w = HistHSMom[12+i]->GetBinWidth(1);
		double peak = HistHSMom[12+i]->GetBinCenter(HistHSMom[12+i]->GetMaximumBin());
		fgaus->SetRange(peak - nw*w, peak + nw *w);
		HistHSMom[12+i]->Fit("fgaus","R");
		cout<<RunList[i]<<endl;
		double mom_fit = fgaus->GetParameter(1);
		double mom_res = fgaus->GetParameter(2);
		if(i<4){
			mom_HSpmcor[i] = mom;
			off_HSpmcor[i] = mom_fit - mom;
			res_HSpmcor[i] = mom_res/mom;
		}
		else if( i < 8){
			mom_HSpicor[i-4] = mom;
			off_HSpicor[i-4] = mom_fit - mom;
			res_HSpicor[i-4] = mom_res/mom;
		}
	}
	TGraphErrors* picor = new TGraphErrors(4,mom_picor,res_picor,dum,dum);
	TGraphErrors* pinocor = new TGraphErrors(4,mom_pinocor,res_pinocor,dum,dum);
	TGraphErrors* HSpicor = new TGraphErrors(4,mom_HSpicor,res_HSpicor,dum,dum);
	TGraphErrors* HSpinocor = new TGraphErrors(4,mom_HSpinocor,res_HSpinocor,dum,dum);
	picor->SetMarkerColor(kRed);
	picor->SetMarkerStyle(23);
	pinocor->SetMarkerStyle(23);
	HSpicor->SetMarkerColor(kBlue);
	HSpinocor->SetMarkerColor(kMagenta);
	HSpicor->SetMarkerStyle(23);
	HSpinocor->SetMarkerStyle(23);
	picor->SetMarkerSize(3);
	pinocor->SetMarkerSize(3);
	HSpicor->SetMarkerSize(3);
	HSpinocor->SetMarkerSize(3);
	TMultiGraph* g = new TMultiGraph();
	TLegend* leg = new TLegend(0.1,0.65,0.7,0.9);
	leg->AddEntry(pinocor,"pi WO Correction");
	leg->AddEntry(HSpinocor,"HSpi WO Correction");
	leg->AddEntry(picor,"pi W Correction");
	leg->AddEntry(HSpicor,"HSpi W Correction");
	g->GetXaxis()->SetTitle("Momentum [GeV/c]");
	g->GetYaxis()->SetTitle("dP / P [%]");
	g->GetYaxis()->SetRangeUser(0,0.12);
	g->Add(picor);
	g->Add(pinocor);
	g->Add(HSpicor);
	g->Add(HSpinocor);
	TCanvas* c3 = new TCanvas("c3","c3",1200,800);
	g->Draw("AP");
	leg->Draw("same");
	TCanvas* c4 = new TCanvas("c4","c4",1200,800);
	TGraphErrors* pmcor = new TGraphErrors(4,mom_pmcor,res_pmcor,dum,dum);
	TGraphErrors* pmnocor = new TGraphErrors(4,mom_pmnocor,res_pmnocor,dum,dum);
	TGraphErrors* HSpmcor = new TGraphErrors(4,mom_HSpmcor,res_HSpmcor,dum,dum);
	TGraphErrors* HSpmnocor = new TGraphErrors(4,mom_HSpmnocor,res_HSpmnocor,dum,dum);
	pmcor->SetMarkerColor(kRed);
	pmcor->SetMarkerStyle(23);
	pmnocor->SetMarkerStyle(23);
	HSpmcor->SetMarkerColor(kBlue);
	HSpmnocor->SetMarkerColor(kMagenta);
	HSpmcor->SetMarkerStyle(23);
	HSpmnocor->SetMarkerStyle(23);
	pmcor->SetMarkerSize(3);
	pmnocor->SetMarkerSize(3);
	HSpmcor->SetMarkerSize(3);
	HSpmnocor->SetMarkerSize(3);
	TMultiGraph* g2 = new TMultiGraph();
	TLegend* leg2 = new TLegend(0.1,0.65,0.7,0.9);
	leg2->AddEntry(pmnocor,"pm WO Correction");
	leg2->AddEntry(HSpmnocor,"HSpm WO Correction");
	leg2->AddEntry(pmcor,"pm W Correction");
	leg2->AddEntry(HSpmcor,"HSpm W Correction");
	g2->GetXaxis()->SetTitle("Momentum [GeV/c]");
	g2->GetYaxis()->SetTitle("dP / P [%]");
	g2->GetYaxis()->SetRangeUser(0,0.12);
	g2->Add(pmcor);
	g2->Add(pmnocor);
	g2->Add(HSpmcor);
	g2->Add(HSpmnocor);
	g2->Draw("AP");
	leg2->Draw("same");
}
