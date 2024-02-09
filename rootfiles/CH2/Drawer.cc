	TChain* chain = new TChain("kurama");
	TChain* chain2 = new TChain("kk");
void Drawer(){
	chain2->Add("run05641_DstHSKKAna.root");
}
void DrawHS(){
	TH1D* h = new TH1D("hist","P_{K^{-} Beam}",100,1.6,2);
	chain2->Draw("pK18>>hist");
	h->Fit("gaus");
	
}
void Draw(int seg,double slope,double offset){
	TCanvas* c1 = new TCanvas("c1","c1",600,400);
	TString dr = Form("utTofSeg[%d]:ytofKurama>>(100,-1000,1000,100,10,25)",seg-1);
	TString dr2 = Form("dtTofSeg[%d]:ytofKurama>>(100,-1000,1000,100,10,25)",seg-1);
	TString dr3 = Form("utTofSeg[%d]:dtTofSeg[%d]>>(100,10,25,100,10,25)",seg-1,seg-1);
	TString dr4 = Form("utTofSeg[%d]/2+dtTofSeg[%d]/2:ytofKurama>>(100,-1000,1000,100,10,25)",seg-1,seg-1);
	TCut dcut = Form("ntSdcOut==1&&tofsegKurama==%d",seg);
	TCut poscut = Form("utTofSeg[%d]-%f*ytofKurama>%f",seg-1,slope,offset);
	TCut poscut2 = Form("dtTofSeg[%d]-%f*ytofKurama>%f",seg-1,-slope,offset);
//	dcut = dcut&&poscut&&poscut2;
	
	c1->cd();
	chain->Draw(dr4,dcut,"colz");
	TF1* f = new TF1("func",Form("%f*x+%f",-slope,offset),-500,500);
	f->Draw("same");
//	chain->Draw(dr2,dcut,"colz");
//	chain->Draw(dr3,dcut,"colz");
}
