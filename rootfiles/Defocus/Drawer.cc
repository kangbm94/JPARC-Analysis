void Drawer(){
}
void Drawer(TString name){
	//gStyle->SetOptStat(0);
	TFile* file = new TFile(name,"READ");
	TTree* tree = (TTree*)file->Get("tpc");
//	TString title = "xBcOut-cluster_x:cluster_z";
	TString title2 = "Before";
	TString title = "After";
	TH2D* h = new TH2D(title,title,50,-250,250,100,-5,5);
	TH2D* h2 = new TH2D(title2,title2,50,-250,250,100,-5,5);
	double cx=50,cy=5;
	TString cut = Form("abs(cluster_x-%f)<25&&abs(cluster_y-%f)<15&&ntBcOut==1&&nhTpc>8",cx,cy);
//	cut = Form("ntBcOut==1&&nhTpc>8");
	TCanvas* c1 = new TCanvas("c1","c1",1600,800);
	c1->Divide(2,1);
	c1->cd(1);
	tree->Draw(Form("xCorVec:cluster_z>>%s",title2.Data()),cut,"col");
	h2->GetXaxis()->SetTitle("z [mm]");
	h2->GetYaxis()->SetTitle("x residual [mm]");
	c1->cd(2);
	tree->Draw(Form("xCorRes:cluster_z>>%s",title.Data()),cut,"col");
	h->GetXaxis()->SetTitle("z [mm]");
	h->GetYaxis()->SetTitle("x residual [mm]");
}
