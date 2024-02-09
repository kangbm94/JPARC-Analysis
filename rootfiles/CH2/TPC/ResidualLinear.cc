TF1* ResDe = new TF1("ResDe","[0]+[1]/x",90,3000);
TF1* fGaus = new TF1("fGaus","gaus",-10,10);
TF1* ResDl = new TF1("ResDl","[0]+[1]*abs(x)",-300,300);
TF1* ResDl2 = new TF1("ResDl2","[0]+[1]*x*x",-300,300);
//ToDoList : CheckLinearTracking in 3D display
int runnum = 5764;
void MakeResidualHists(){

	TFile* file = new TFile("run05764_DstTPCTracking.root");
	TTree* tree = (TTree*)file->Get("tpc");
//	TFile* file2 = new TFile(Form("run0%d_hists_nocut.root",runnum),"recreate");
	TFile* file2 = new TFile(Form("run0%d_hists.root",runnum),"recreate");
	TH2D* ResXDe = new TH2D("ResX_De","ResX_De",1000,0,3000,1000,-10,10);
	TH2D* ResYDe = new TH2D("ResY_De","ResY_De",1000,0,3000,1000,-10,10);
	TH2D* ResZDe = new TH2D("ResZ_De","ResZ_De",1000,0,3000,1000,-3,3);
	TH2D* ResXDeCl = new TH2D("ResX_DeCl","ResX_DeCl",1000,0,3000,1000,-10,10);
	TH2D* ResYDeCl = new TH2D("ResY_DeCl","ResY_DeCl",1000,0,3000,1000,-10,10);
	TH2D* ResZDeCl = new TH2D("ResZ_DeCl","ResZ_DeCl",1000,0,3000,1000,-3,3);
	TH2D* ResXDl = new TH2D("ResX_Dl","ResX_Dl",1000,-300,300,1000,-10,10);
	TH2D* ResYDl = new TH2D("ResY_Dl","ResY_Dl",1000,-300,300,1000,-10,10);
	TH2D* ResZDl = new TH2D("ResZ_Dl","ResZ_Dl",1000,-300,300,1000,-3,3);
	TH2D* ResRDl = new TH2D("ResR_Dl","ResR_Dl",1000,-300,300,1000,-10,10);
	TH2D* ResXDlCl = new TH2D("ResX_DlCl","ResX_DlCl",1000,-300,300,1000,-10,10);
	TH2D* ResYDlCl = new TH2D("ResY_DlCl","ResY_DlCl",1000,-300,300,1000,-10,10);
	TH2D* ResZDlCl = new TH2D("ResZ_DlCl","ResZ_DlCl",1000,-300,300,1000,-3,3);
	TH2D* ResRDlCl = new TH2D("ResR_DlCl","ResR_DlCl",1000,-300,300,1000,-10,10);

	TH2D* ResRdzCl = new TH2D("ResR_dzCl","ResR_dzCl",1000,-2,2,1000,-10,10);
	TH2D* ResRdz = new TH2D("ResR_dz","ResR_dz",1000,-2,2,1000,-10,10);

	TH2D* ResXdzCl = new TH2D("ResX_dzCl","ResX_dzCl",1000,-2,2,1000,-10,10);
	TH2D* ResXdz = new TH2D("ResX_dz","ResX_dz",1000,-2,2,1000,-10,10);
	TH2D* ResYdzCl = new TH2D("ResY_dzCl","ResY_dzCl",1000,-2,2,1000,-10,10);
	TH2D* ResYdz = new TH2D("ResY_dz","ResY_dz",1000,-2,2,1000,-10,10);
	TH2D* ResZdzCl = new TH2D("ResZ_dzCl","ResZ_dzCl",1000,-2,2,1000,-10,10);
	TH2D* ResZdz = new TH2D("ResZ_dz","ResZ_dz",1000,-2,2,1000,-10,10);
	
	TString XDe = Form("residual_x:track_cluster_de>>ResX_De");
	TString YDe = Form("residual_y:track_cluster_de>>ResY_De");
	TString ZDe = Form("residual_z:track_cluster_de>>ResZ_De");
	TString XDeCl = Form("residual_x:track_cluster_de>>ResX_DeCl");
	TString YDeCl = Form("residual_y:track_cluster_de>>ResY_DeCl");
	TString ZDeCl = Form("residual_z:track_cluster_de>>ResZ_DeCl");
	TString XDl = Form("residual_x:calpos_y>>ResX_Dl");
	TString YDl = Form("residual_y:calpos_y>>ResY_Dl");
	TString ZDl = Form("residual_z:calpos_y>>ResZ_Dl");
	TString RDl = Form("sqrt(residual_x*residual_x+residual_z*residual_z):calpos_y>>ResR_Dl");
	TString XDlCl = Form("residual_x:calpos_y>>ResX_DlCl");
	TString YDlCl = Form("residual_y:calpos_y>>ResY_DlCl");
	TString ZDlCl = Form("residual_z:calpos_y>>ResZ_DlCl");
	TString RDlCl = Form("sqrt(residual_x*residual_x+residual_z*residual_z):calpos_y>>ResR_DlCl");
	
	
	TString Xdz = Form("residual_x:helix_dz>>ResX_dz");
	TString Zdz = Form("residual_z:helix_dz>>ResZ_dz");
	TString Rdz = Form("sqrt(residual_x*residual_x+residual_z*residual_z):helix_dz>>ResR_dz");
	TString Ydz = Form("residual_y:helix_dz>>ResY_dz");
	TString XdzCl = Form("residual_x:helix_dz>>ResX_dzCl");
	TString ZdzCl = Form("residual_z:helix_dz>>ResZ_dzCl");
	TString RdzCl = Form("sqrt(residual_x*residual_x+residual_z*residual_z):helix_dz>>ResR_dzCl");
	TString YdzCl = Form("residual_y:helix_dz>>ResY_dzCl");

	TCut ClSize	="track_cluster_size<2";
	TCut ClSize1	="track_cluster_size>1";
	TCut ChiCut = "chisqr<20";
//	TCut ChiCut = "chisqrTpc<20";
	TCut HitCut = "nhtrack > 3";
	TCut ntCut = "ntTpc<5";
//	TCut VertCut = "abs(x0Tpc)<25&&abs(y0Tpc)<25&&(abs(v0Tpc)>0.1||abs(u0Tpc)>0.1)";
	TCut VertCut = "hough_flag<400&&mom0<0.6&&abs(helix_dz)>0.00&&path>200&&(abs(hitpos_x)>25||abs(hitpos_y)>25)";
	ClSize = ClSize&&ChiCut&&HitCut&&VertCut&&ntCut;
	ClSize1 = ClSize1&&ChiCut&&HitCut&&VertCut&&ntCut;
	auto tag1 = new TNamed("Cut1",ClSize);
	auto tag2 = new TNamed("Cut2",ClSize1);
	tag1->Write();
	tag2->Write();
	/*
	tree->Draw(XDeCl,ClSize,"colz");
	tree->Draw(YDeCl,ClSize,"colz");
	tree->Draw(ZDeCl,ClSize,"colz");
	tree->Draw(XDe,ClSize1,"colz");
	tree->Draw(YDe,ClSize1,"colz");
	tree->Draw(ZDe,ClSize1,"colz");
*/
//	tree->Draw(XDlCl,ClSize,"colz");
	tree->Draw(YDlCl,ClSize,"colz");
//	tree->Draw(ZDlCl,ClSize,"colz");
	tree->Draw(RDlCl,ClSize,"colz");
//	tree->Draw(XDl,ClSize1,"colz");
	tree->Draw(YDl,ClSize1,"colz");
//	tree->Draw(ZDl,ClSize1,"colz");
	tree->Draw(RDl,ClSize1,"colz");
	
	tree->Draw(Xdz,ClSize1,"colz");
	tree->Draw(Ydz,ClSize1,"colz");
	tree->Draw(Zdz,ClSize1,"colz");
	tree->Draw(Rdz,ClSize1,"colz");
	
	tree->Draw(XdzCl,ClSize,"colz");
	tree->Draw(YdzCl,ClSize,"colz");
	tree->Draw(ZdzCl,ClSize,"colz");
	tree->Draw(RdzCl,ClSize,"colz");
	file2->Write();
}
void ResidualLinear(){
}
void Residual(int dum){
//	TFile* file = new TFile("run05721_hists.root");
	TFile* file = new TFile(Form("run0%d_hists.root",runnum));
//	TFile* file = new TFile(Form("run0%d_hists_nocut.root",runnum));
	TH2D* ResXDe = (TH2D*)file->Get("ResX_De");
	TH2D* ResYDe = (TH2D*)file->Get("ResY_De");
	TH2D* ResZDe = (TH2D*)file->Get("ResZ_De");
	TH2D* ResXDeCl = (TH2D*)file->Get("ResX_DeCl");
	TH2D* ResYDeCl = (TH2D*)file->Get("ResY_DeCl");
	TH2D* ResZDeCl = (TH2D*)file->Get("ResZ_DeCl");

	TH2D* ResXDl = (TH2D*)file->Get("ResX_Dl");
	TH2D* ResYDl = (TH2D*)file->Get("ResY_Dl");
	TH2D* ResZDl = (TH2D*)file->Get("ResZ_Dl");
	TH2D* ResXDlCl = (TH2D*)file->Get("ResX_DlCl");
	TH2D* ResYDlCl = (TH2D*)file->Get("ResY_DlCl");
	TH2D* ResZDlCl = (TH2D*)file->Get("ResZ_DlCl");
	TH2D* ResYDz = (TH2D*)file->Get("ResY_dz");
	TH2D* ResYDzCl = (TH2D*)file->Get("ResY_DlCl");
	double rangede[13]= {150,200,250,300,350,400,450,500,550,600,650,700,800};
	double rangedl[13]= {-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250};
	double rangedz[13]= {-1,-0.75,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.4,0.5,0.75,1};
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(6,2);
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(6,2);
	TH1D* ResXde[12];
	TH1D* ResXdl[12];
	TH1D* ResXCldl[12];
	TH1D* ResYdl[12];
	TH1D* ResYCldl[12];
	TH1D* ResYdz[12];
	TH1D* ResYCldz[12];
	TH1D* ResZdl[12];
	TH1D* ResZCldl[12];
	TGraph* RDe = new TGraph();
	TGraph* RXDl = new TGraph();
	TGraph* RYDl = new TGraph();
	TGraph* RZDl = new TGraph();
	TGraph* RXDlCl = new TGraph();
	TGraph* RYDlCl = new TGraph();
	TGraph* RZDlCl = new TGraph();

	TGraph* RYDzCl = new TGraph();
	TGraph* RYDz = new TGraph();
	c1->cd(1);
	ResXDl->Draw("colz");
	c1->cd(2);
	ResYDl->Draw("colz");
	c1->cd(3);
	ResZDl->Draw("colz");
	c1->cd(4);
	ResXDlCl->Draw("colz");
	c1->cd(5);
	ResYDlCl->Draw("colz");
	c1->cd(6);
	ResZDlCl->Draw("colz");
	for(int i=0;i<12;++i){
		c1->cd(i+1);
		TString title = Form("ResX_%d",i);
		ResXde[i] = new TH1D(title,title,1000,-3,3);
		int b1 = ResXDe->GetXaxis()->FindBin(rangede[i]);
		int b2 = ResXDe->GetXaxis()->FindBin(rangede[i+1]);
		cout<<"Bin  "<<b1<<" , "<<b2<<endl;
//		ResXDe->ProjectionY(title,b1,b2);
		double mean = ResXde[i]->GetMean();
		double sig = ResXde[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
//		ResXde[i]->Fit("fGaus","R");
		double de = (rangede[i]+rangede[i+1])/2;
		double res = fGaus->GetParameter(2);
//		RDe->AddPoint(de,res*res);
		RDe->AddPoint(de,res);
	
	}
/*
	for(int i=0;i<12;++i){
		TString titlex = Form("ResX_%d",i);
		TString titlexcl = Form("ResXCl_%d",i);
		TString titley = Form("ResY_%d",i);
		TString titleycl = Form("ResYCl_%d",i);
		TString titlez = Form("ResZ_%d",i);
		TString titlezcl = Form("ResZCl_%d",i);
		ResXdl[i] = new TH1D(titlex,titlex,100,-3,3);
		ResXCldl[i] = new TH1D(titlex,titlex,100,-3,3);
		ResYdl[i] = new TH1D(titley,titley,100,-3,3);
		ResYCldl[i] = new TH1D(titleycl,titleycl,100,-3,3);
		ResZdl[i] = new TH1D(titlez,titlez,100,-3,3);
		ResZCldl[i] = new TH1D(titlezcl,titlezcl,100,-3,3);
		int b1 = ResYDl->GetXaxis()->FindBin(rangedl[i]);
		int b2 = ResYDl->GetXaxis()->FindBin(rangedl[i+1]);
		ResXDl->ProjectionY(titlex,b1,b2);
		ResYDl->ProjectionY(titley,b1,b2);
		ResZDl->ProjectionY(titlez,b1,b2);
		ResXDlCl->ProjectionY(titlexcl,b1,b2);
		ResYDlCl->ProjectionY(titleycl,b1,b2);
		ResZDlCl->ProjectionY(titlezcl,b1,b2);
		double mean = ResXdl[i]->GetMean();
		double sig = ResXdl[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		ResXdl[i]->Fit("fGaus","R");
		double dl = (rangedl[i]+rangedl[i+1])/2;
		double res = fGaus->GetParameter(2);
		RXDl->AddPoint(dl,res);
		mean = ResYdl[i]->GetMean();
		sig = ResYdl[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		cout<<"ResYdl"<<i<<endl;
		ResYdl[i]->Fit("fGaus","R0");
		res = fGaus->GetParameter(2);
		RYDl->AddPoint(dl,res);
		mean = ResZdl[i]->GetMean();
		sig = ResZdl[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		cout<<"ResZdl"<<i<<endl;
		ResZdl[i]->Fit("fGaus","R0");
		res = fGaus->GetParameter(2);
		RZDl->AddPoint(dl,res);
		c2->cd(i+1);
		ResYdl[i]->Draw("");
		ResYdl[i]->SetLineColor(kBlue);
		ResZdl[i]->Draw("same");
		ResZdl[i]->SetLineColor(kRed);
		mean = ResXCldl[i]->GetMean();
		sig = ResXCldl[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		ResXCldl[i]->Fit("fGaus","R");
		res = fGaus->GetParameter(2);
		RXDlCl->AddPoint(dl,res);
		mean = ResYCldl[i]->GetMean();
		sig = ResYCldl[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		ResYCldl[i]->Fit("fGaus","R");
		res = fGaus->GetParameter(2);
		RYDlCl->AddPoint(dl,res);
		mean = ResZCldl[i]->GetMean();
		sig = ResZCldl[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		ResZCldl[i]->Fit("fGaus","R");
		res = fGaus->GetParameter(2);
		RZDlCl->AddPoint(dl,res);
	}
	*/
	TCanvas* c3 = new TCanvas("c3","c3",100,100,600,600);
	c3->Divide(3,2);
	/*
	c3->cd(1);
	RXDl->Draw("AP");
	RXDl->SetTitle("Dl:ResX[clsize<2]");
	RXDl->SetMarkerStyle(25);
	RXDl->SetMarkerSize(1);
	c3->cd(2);
	RYDl->Draw("AP");
	RYDl->SetTitle("Dl:ResY[clsize<2]");
	RYDl->SetMarkerStyle(25);
	RYDl->SetMarkerSize(1);
	c3->cd(3);
	RZDl->Draw("AP");
	RZDl->SetTitle("Dl:ResZ[clsize<2]");
	RZDl->SetMarkerStyle(25);
	RZDl->SetMarkerSize(1);
	c3->cd(4);
	RXDlCl->Draw("AP");
	RXDlCl->SetTitle("Dl:ResX[clsize>1]");
	RXDlCl->SetMarkerStyle(25);
	RXDlCl->SetMarkerSize(1);
	c3->cd(5);
	RYDlCl->Draw("AP");
	RYDlCl->SetTitle("Dl:ResY[clsize>1]");
	RYDlCl->SetMarkerStyle(25);
	RYDlCl->SetMarkerSize(1);
	c3->cd(6);
	RZDlCl->Draw("AP");
	RZDlCl->SetTitle("Dl:ResZ[clsize>1]");
	RZDlCl->SetMarkerStyle(25);
	RZDlCl->SetMarkerSize(1);
*/	
	for(int i=0;i<12;++i){
		TString titley = Form("ResYdz_%d",i);
		TString titleycl = Form("ResYCldz_%d",i);
		ResYdz[i] = new TH1D(titley,titley,100,-3,3);
		ResYCldz[i] = new TH1D(titleycl,titleycl,100,-3,3);
		int b1 = ResYDz->GetXaxis()->FindBin(rangedz[i]);
		int b2 = ResYDzCl->GetXaxis()->FindBin(rangedz[i+1]);
		c2->cd(i+1);
		ResYDz->ProjectionY(titley,b1,b2);
		ResYDzCl->ProjectionY(titleycl,b1,b2);
		double mean = ResYdz[i]->GetMean();
		double sig = ResYdz[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		ResYdz[i]->Fit("fGaus","R");
		double dz = (rangedz[i]+rangedz[i+1])/2;
		double res = fGaus->GetParameter(2);
		RYDz->AddPoint(dz,res);
		
		mean = ResYCldz[i]->GetMean();
		sig = ResYCldz[i]->GetStdDev();
		fGaus->SetRange(mean-sig,mean+sig);
		cout<<"ResYdz"<<i<<endl;
		ResYCldz[i]->Fit("fGaus","R0");
		res = fGaus->GetParameter(2);
		RYDzCl->AddPoint(dz,res);
		
		
	}
	c3->cd(1);
	RYDz->Draw("AP");
	c3->cd(2);
	RYDzCl->Draw("AP");
}
