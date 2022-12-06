void Fitter(){
	TFile* file1 = new TFile("run05754_DstTPCBcSim.root","READ");
	TFile* file2 = new TFile("run05754_DstTPCBcOut.root","READ");
	double par1[32];
	double par2[32];
	TCanvas* c1 = new TCanvas("c1","c1",1600,700);
	c1->Divide(3,2);
	TCanvas* c2 = new TCanvas("c2","c2",1600,700);
	c2->Divide(3,2);
	TF1* fGaus = new TF1("fGaus","gaus",-15,15);
	int j=1;
	double k[32];
	for(int i=0;i<32;++i){
		TString ht = Form("h%d",200002+1000*(i+1));
		auto h = (TH2D*)file1->Get(ht);
		k[i]=i;	
		if(i%6==0) {c1->cd(j);j++; h->Draw("colz");}
		h->Fit("fGaus","0");
		par1[i]=fGaus->GetParameter(2);
	}
	j=1;
	for(int i=0;i<32;++i){
		TString ht = Form("h%d",200002+1000*(i+1));
		auto h = (TH2D*)file2->Get(ht);
		
		if(i%6==0) {c2->cd(j);j++; h->Draw("colz");}
		h->Fit("fGaus","0");
		par2[i]=fGaus->GetParameter(2);
	}
	TGraph* gr1 = new TGraph(32,k,par1);
	TGraph* gr2 = new TGraph(32,k,par2);
	TCanvas* c3 = new TCanvas("c3","c3",1200,600);
	c3->cd();
	gr1->SetTitle("Sim");
	gr1->Draw("ALP");
	gr1->SetLineColor(kRed);
	gr1->SetLineWidth(3);
	gr1->GetYaxis()->SetRangeUser(1,6);
	gr2->SetTitle("Real");
	gr2->Draw("LPSame");
	gr2->SetLineColor(kBlue);
	gr2->SetLineWidth(3);
	
	auto leg = new TLegend(0.1,0.7,0.48,0.9);
	leg->AddEntry(gr1,"Sim","l");
	leg->AddEntry(gr2,"Real","l");
	leg->Draw("same");
}
