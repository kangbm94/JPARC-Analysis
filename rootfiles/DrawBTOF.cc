void DrawBTOF(){
	cout<<"DrawBTOF"<<endl;
}
TString hn = "HodoscopeNoMT05477.root";
void DrawBTOF(int run){
	gStyle->SetOptFit(1111);
//	TFile* file = new TFile(hn,"READ");
//	TTree* tree = (TTree*)file->Get("tree");
	TChain* tree = new TChain("tree");
	tree->Add("Hodoscope05475.root");
	tree->Add("Hodoscope05476.root");
	tree->Add("Hodoscope05477.root");
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(4,2);
	TH1D* ht[8];
	for(int i=0;i<8;i++){
		c1->cd(i+1);
		ht[i]=new TH1D(Form("ht%d",i+1),Form("ht%d",i+1),120,-2,-0.8);
		cout<<"BH2Seg = "<< i+1<<endl;
		tree->Draw(Form("cbtof>>ht%d",i+1),Form("bacnhits==0&&bh2nhits==1&&bh2hitpat[0]==%d",i+1));
		gPad->SetGrid();
		ht[i]->Fit("gaus");
	}
}
