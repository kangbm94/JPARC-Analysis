#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
void DrawK18Mom(){
	SetStyle();
	TFile* file = new TFile("run05641_DstHSKKAna.root");
	TTree* tree= (TTree*)file->Get("kk");
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	TH1D* h = new TH1D("pk18", "K^{-} momentum",100,1.6,2.0);
	h->GetYaxis()->SetNdivisions(5);
	h->GetXaxis()->SetTitle("BeamMomntum [GeV / c]");
	tree->Draw("pK18>>pk18");
	c1->SaveAs("pk18.png");
}
