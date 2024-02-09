int GetBinContent(TH1D* hist, int lb, int hb){
	int ent = 0;
	for(int i=lb;i<hb+1;++i){
		ent +=hist -> GetBinContent(i);
	}
	return ent;
}
void Compare(){
	cout<<"ComparePolarity()"<<endl;
	cout<<"CompareCoplanarity()"<<endl;
}
void ComparePolarity(){
	int rn_s = 5641;
	int rn_e = 5666;
	TFile* fileCh2 = new TFile(Form("GenfitPolaAnal_%d_%d.root",rn_s,rn_e));
	TTree* treeCh2 = (TTree*)fileCh2->Get("tree");
	rn_s = 5667;
	rn_e = 5697;
	TFile* fileCarbon = new TFile(Form("GenfitPolaAnal_%d_%d.root",rn_s,rn_e));
	TTree* treeCarbon = (TTree*)fileCarbon->Get("tree");
}
void CompareCoplanarity(){
	int rn_s = 5641;
	int rn_e = 5666;
	TFile* fileCh2 = new TFile(Form("GenfitPolaAnal_%d_%d.root",rn_s,rn_e));
	TTree* treeCh2 = (TTree*)fileCh2->Get("tree");
	rn_s = 5667;
	rn_e = 5697;
	TFile* fileCarbon = new TFile(Form("GenfitPolaAnal_%d_%d.root",rn_s,rn_e));
	TTree* treeCarbon = (TTree*)fileCarbon->Get("tree");

	double CopCh2,CopCarbon;
	double MissMassCh2,MissMassCarbon;
	treeCh2->SetBranchAddress("Coplanarity",&CopCh2);
	treeCh2->SetBranchAddress("MissMass",&MissMassCh2);
	
	treeCarbon->SetBranchAddress("Coplanarity",&CopCarbon);
	treeCarbon->SetBranchAddress("MissMass",&MissMassCarbon);
	int entCh2 = treeCh2 ->GetEntries();
	int entCarbon = treeCarbon ->GetEntries();
	int nbin = 200;
	TH1D* HistCopCh2 = new TH1D("HCh2","HCh2",nbin,-1,1);
	TH1D* HistCopCarbon = new TH1D("HCarbon","HCarbon",nbin,-1,1);
	for(int i=0;i<entCh2;++i){
		treeCh2->GetEntry(i);
		if(abs(MissMassCh2 - 1.321)<0.1){
			HistCopCh2->Fill(CopCh2);
		}
	}
	for(int i=0;i<entCarbon;++i){
		treeCarbon->GetEntry(i);
		if(abs(MissMassCarbon - 1.321)<0.1){
			HistCopCarbon->Fill(CopCarbon);
		}
	}
	int eCh2 = HistCopCh2->GetEntries();
	int eCarbon = HistCopCarbon->GetEntries();
	double rat = eCh2 * 1./ eCarbon;
	cout<<eCh2<<" , "<<eCarbon<<endl;
	cout<<"CH2/Carbon ratio : "<<rat<<endl;
	HistCopCarbon->Scale(rat);
	HistCopCh2->SetLineColor(kRed);
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	HistCopCh2->Draw("");
	gStyle->SetErrorX(0);
	HistCopCarbon->Draw("Samehist");
	double bw = HistCopCh2->GetBinWidth(1);
	double CopW[39];
	double AcptW[39];
	for(int i=1;i<40;++i){
		int lb = nbin/2 + 1 - i;
		int hb = nbin/2 + 1 + i;
		int HitCh2 = GetBinContent(HistCopCh2,lb,hb);
		int HitCarbon = GetBinContent(HistCopCarbon,lb,hb);
		CopW[i-1]= 2*i * bw;	
		AcptW[i-1]= {(HitCh2 -HitCarbon)* 1./HitCh2 };	
	}
	TCanvas* c2 = new TCanvas("c2","c2",1200,800);
	TGraph* g1 = new TGraph(29,CopW,AcptW);
	g1->SetTitle("P Separation Power");
	g1->GetXaxis()->SetTitle("Coplanarity Window");
	g1->GetYaxis()->SetTitle("#frac{N_{CH2}-N_{Dia}}{N_{CH2}}");
	g1->SetMarkerStyle(25);
	g1->SetMarkerSize(1);
	g1->Draw("AP");
}
