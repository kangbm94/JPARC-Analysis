double PI = acos(-1);
double mm1,mm2,MLd,tcm;
TFile* f1;TFile* f2;
TTree* tr1;
TTree* tr2;
void DrawXi(){
	f1 = new TFile("SelectedEvents.root");
	f2 = new TFile("TPCInv.root");
	tr1 = (TTree*)f1->Get("tree");
	tr1->SetBranchAddress("XiM2",&mm1);
	tr1->SetBranchAddress("XiThetaCM",&tcm);
	tr2 = (TTree*)f2->Get("tree");
	tr2->SetBranchAddress("MM",&mm2);
	tr2->SetBranchAddress("InvMLd",&MLd);
}

void DrawDifCr(){
	

	TGraphErrors* g1;
	TGraphErrors* g2;
	TH1D* hist1 = new TH1D("Xi Count vs cos(#theta_{CM})","Xi Count vs cos(#theta_{CM})",10,0.75,1);
	hist1->GetXaxis()->SetTitle("cos(#theta_{CM})");
	hist1->GetYaxis()->SetTitle("Entry [ / 0.025]");
	TH1D* hist2 = new TH1D("Xi* Count vs cos(#theta_{CM})","Xi* Count vs cos(#theta_{CM})",10,0.75,1);
	hist1->GetXaxis()->SetTitle("cos(#theta_{CM})");
	hist1->GetYaxis()->SetTitle("Entry [ / 0.025]");
	int ent = tr1->GetEntries();	
	for(int i=0;i<ent;++i){
		tr1->GetEntry(i);
		if(abs(mm1-1.32657)<0.0145*3) hist1->Fill(cos(tcm/180*PI));
		if(abs(mm1-1.540)<0.0135*3) hist2->Fill(cos(tcm/180*PI));
	}
	double xic[10],xisc[10];
	for(int i=0;i<10;++i){
		xic[i]=hist1->GetBinContent(i+1);
		xisc[i]=hist2->GetBinContent(i+1);
		cout<<Form("%d : %f",i,xic[i])<<endl; 
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	hist1->Draw();
	c1->cd(2);
	hist2->Draw();
}
void DrawXi(int dum){
	TH1D* hist1 = new TH1D("KuramaMM","KuramaMM",200,1,2);
	hist1->GetXaxis()->SetTitle("MissMass [GeV / c2]");
	hist1->GetYaxis()->SetTitle("Entry [ /5 MeV]");
	TH1D* hist2 = new TH1D("KuramaMMCut","KuramaMMCut",200,1,2);
	int ent = tr1->GetEntries();
	for(int i=0;i<ent;++i){
		tr1->GetEntry(i);
		hist1->Fill(mm1);
	}
	ent = tr2->GetEntries();
	for(int i=0;i<ent;++i){
		tr2->GetEntry(i);
		if(abs(MLd-1.12)<0.044)hist2->Fill(mm2);
		//		if(MLd>0)	hist2->Fill(mm2);
	}

	hist2->SetLineColor(kRed);
	TCanvas*c1 = new TCanvas("c1","c1",1200,800);
	c1->cd();
	hist1->Draw();
	hist2->Draw("same");
}
void DrawMMAngular(){
	double mm,theta;
	TFile* f1 = new TFile("SelectedEvents.root");
	TTree* tr1 = (TTree*)f1->Get("tree");
	tr1->SetBranchAddress("XiMM",&mm);
	tr1->SetBranchAddress("XiThetaCM",&theta);
	TH2D* hist1 = new TH2D("MM vs cos(#theta_{CM})","MM vs cos(#theta_{CM})",50,0.5,1,50,1.2,1.7);
	int ent = tr1->GetEntries();
	gStyle->SetOptStat(0000);
	hist1->GetXaxis()->SetTitle("cos(#theta_{CM})");
	hist1->GetYaxis()->SetTitle("MissMass[GeV/c^2]");
	for(int i=0;i<ent;++i){
		tr1->GetEntry(i);
		hist1->Fill(cos(theta/180 * PI),mm);
	}
	hist1->Draw("colz");
}
