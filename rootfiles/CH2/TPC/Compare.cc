#include "/Users/MIN/ROOTSharedLibs/MyStyle.hh"
double mk = 493.677/1000;
double mXi = 1321.71/1000;
double mp = 938.272/1000;
double CopCh2,CopCarbon;
double MissMassCh2,MissMassCarbon;
double dMCh2,dMCarbon;
double PKm_xCh2,PKm_yCh2,PKm_zCh2;
double PKp_xCh2,PKp_yCh2,PKp_zCh2;
double PXi_xCh2,PXi_yCh2,PXi_zCh2;
double PP_xCh2,PP_yCh2,PP_zCh2;
double PKm_xCarbon,PKm_yCarbon,PKm_zCarbon;
double PKp_xCarbon,PKp_yCarbon,PKp_zCarbon;
double PXi_xCarbon,PXi_yCarbon,PXi_zCarbon;
double PP_xCarbon,PP_yCarbon,PP_zCarbon;
TTree* treeCh2;
TTree* treeCarbon;
TFile* fileCh2; 
TFile* fileCarbon; 


double MMCut = 0.015;
int GetBinContent(TH1D* hist, int lb, int hb){
	int ent = 0;
	for(int i=lb;i<hb+1;++i){
		ent +=hist -> GetBinContent(i);
	}
	return ent;
}
double ToSigma(double p){
	double range = 0.5+0.5*p;
	double sigma = TMath::NormQuantile(range);
	return sigma;
}
void ComparePolarity(){
}
void CompareCoplanarity(){

	
	int entCh2 = treeCh2 ->GetEntries();
	int entCarbon = treeCarbon ->GetEntries();
	int nbin = 200;
	TH1D* HistCopCh2 = new TH1D("HCh2","HCh2",nbin,-1,1);
	TH1D* HistCopCarbon = new TH1D("HCarbon","HCarbon",nbin,-1,1);
	double dM2Cut = 0.011;
	for(int i=0;i<entCh2;++i){
		treeCh2->GetEntry(i);
		if(abs(MissMassCh2 - 1.321)<MMCut and abs(dMCh2*dMCh2)<dM2Cut){
			HistCopCh2->Fill(CopCh2);
		}
	}
	for(int i=0;i<entCarbon;++i){
		treeCarbon->GetEntry(i);
		if(abs(MissMassCarbon - 1.321)<MMCut and (dMCarbon*dMCarbon)<dM2Cut){
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
	g1->GetYaxis()->SetTitle("Separation[#sigma]");
	g1->SetMarkerStyle(25);
	g1->SetMarkerSize(1);
	g1->Draw("AP");
}
void CompareMissMass(){
	
	double dM2Cut = 0.012;
	int entCh2 = treeCh2 ->GetEntries();
	int entCarbon = treeCarbon ->GetEntries();
	int nbin = 200;
	int nbin2d = 100;
	TH1D* HistCopCh2 = new TH1D("HCopCh2","HCopCh2",nbin,-1,1);
	TH1D* HistDMMCh2 = new TH1D("HMMCh2","HMMCh2",nbin,-0.1,0.1);
	TH2D* Hist2DmCh2 = new TH2D("HCop:MMCh2","HCop:MMCh2",nbin2d,-1,1,nbin2d,-0.1,0.1);

	TH1D* HistCopCarbon = new TH1D("HCopCarbon","HCopCarbon",nbin,-1,1);
	TH1D* HistDMMCarbon = new TH1D("HMMCarbon","HMMCarbon",nbin,-0.1,0.1);
	TH2D* Hist2DmCarbon = new TH2D("HCop:MMCarbon","HCop:MMCarbon",nbin2d,-1,1,nbin2d,-0.1,0.1);
	for(int i=0;i<entCh2;++i){
		treeCh2->GetEntry(i);
		if(abs(MissMassCh2 - 1.321)<MMCut){
			double dM = dMCh2;
			HistDMMCh2->Fill(dM*abs(dM));
			HistCopCh2->Fill(CopCh2);
			Hist2DmCh2->Fill(CopCh2,dM*abs(dM));
		}
	}
	for(int i=0;i<entCarbon;++i){
		treeCarbon->GetEntry(i);
		if(abs(MissMassCarbon - 1.321)<MMCut){
			double dM = dMCarbon;
			HistDMMCarbon->Fill(dM*abs(dM));
			HistCopCarbon->Fill(CopCarbon);
			Hist2DmCarbon->Fill(CopCarbon,dM*abs(dM));
		}
	}
	int eCh2 = HistDMMCh2->GetEntries();
	int eCarbon = HistDMMCarbon->GetEntries();
	double rat = eCh2 * 1./ eCarbon;
	cout<<eCh2<<" , "<<eCarbon<<endl;
	cout<<"CH2/Carbon ratio : "<<rat<<endl;
	HistDMMCarbon->Scale(rat);
	HistDMMCh2->SetLineColor(kRed);
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	gPad->SetMargin(0.15,0.05,0.15,0.05);
	HistDMMCh2->Draw("");
	HistDMMCh2->GetXaxis()->SetTitle("MM^{2}p(K^{-},K^{+}#Xi)^{2} [ GeV^{2}/c^{4}]");
	HistDMMCh2->GetXaxis()->SetLabelSize(0.05);
	HistDMMCh2->GetYaxis()->SetLabelSize(0.05);
	HistDMMCh2->GetYaxis()->SetTitleSize(0.05);
	HistDMMCh2->GetYaxis()->SetTitle("Counts/0.001  [1/GeV^{2}/c^{4}]");
	HistDMMCh2->GetYaxis()->SetNdivisions(5);
	HistDMMCarbon->Draw("Samehist");
//	gPad->SetLogy();
	TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
	leg->AddEntry(HistDMMCh2,"CH_{2} Target","l");
	leg->AddEntry(HistDMMCarbon,"Diamond Target","l");
	leg->Draw();
	leg->SetBorderSize(0);
	TArrow* Left = new TArrow(-dM2Cut,20,-dM2Cut,0);
	TArrow* Right = new TArrow(dM2Cut,20,dM2Cut,0);
	Left->Draw();
	Right->Draw();
	c1->SaveAs("dMComparison.png");
	gStyle->SetErrorX(0);
	double bw = HistDMMCh2->GetBinWidth(1);
	double DMMW[39];
	double AcptW[39];
	for(int i=1;i<40;++i){
//		int lb = 0 ;
//		int hb = 1 + i;
		int lb = nbin/2 + 1 - i;
		int hb = nbin/2 + 1 + i;
		int HitCh2 = GetBinContent(HistDMMCh2,lb,hb);
		int HitCarbon = GetBinContent(HistDMMCarbon,lb,hb);
		DMMW[i-1]= i * bw;	
		double p = {(HitCh2 -HitCarbon)* 1./HitCh2 };	
		AcptW[i-1]= p;
//			ToSigma(p);
	}
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	gPad->SetMargin(0.15,0.05,0.15,0.05);
	TGraph* g1 = new TGraph(29,DMMW,AcptW);
	g1->SetTitle("Proton Purity");
	g1->GetXaxis()->SetTitle("MM^{2} Window [ GeV^2/c^{4}]");
	g1->GetXaxis()->SetLabelSize(0.05);
	g1->GetYaxis()->SetLabelSize(0.05);
	g1->GetYaxis()->SetTitleSize(0.05);
//	g1->GetYaxis()->SetTitle("#frac{N_{CH2}-N_{Dia}}{N_{CH2}}");
	g1->GetYaxis()->SetTitle("Separation Power[#sigma]");
	g1->GetYaxis()->SetNdivisions(5);
	g1->SetMarkerStyle(25);
	g1->SetMarkerSize(1);
	g1->Draw("AP");
	g1->GetYaxis()->SetRangeUser(0,2);
	TArrow* Window = new TArrow(dM2Cut,1,dM2Cut,0);
	Window->Draw();
	c2->SaveAs("Separation.png");
	/*
	TCanvas* c3 = new TCanvas("c3","c3",1200,800);
	c3->Divide(2,1);
	c3->cd(1);
	Hist2DmCh2->Draw("colz");
	c3->cd(2);
	Hist2DmCarbon->Draw("colz");
*/
}

/*
void CompareMom(){
	double MMCut = 0.05;
	int entCh2 = treeCh2 ->GetEntries();
	int entCarbon = treeCarbon ->GetEntries();
	int nbin = 500;
	int nbin2d = 100;
	TH1D* HistCopCh2 = new TH1D("HCopCh2","HCopCh2",nbin,-1,1);
	TH1D* HistDMMCh2 = new TH1D("HMMCh2","HMMCh2",nbin,-0.1,0.1);
	TH2D* Hist2DmCh2 = new TH2D("HCop:MMCh2","HCop:MMCh2",nbin2d,-1,1,nbin2d,-0.1,0.1);

	TH1D* HistCopCarbon = new TH1D("HCopCarbon","HCopCarbon",nbin,-1,1);
	TH1D* HistDMMCarbon = new TH1D("HMMCarbon","HMMCarbon",nbin,-0.1,0.1);
	TH2D* Hist2DmCarbon = new TH2D("HCop:MMCarbon","HCop:MMCarbon",nbin2d,-1,1,nbin2d,-0.1,0.1);
	for(int i=0;i<entCh2;++i){
		treeCh2->GetEntry(i);
		if(abs(MissMassCh2 - 1.321)<MMCut){
			double dP = dPCh2;
			//			double dP = (LVXi-LVMM).Mag();
			HistDMMCh2->Fill(dP*dP);
			HistCopCh2->Fill(CopCh2);
			Hist2DmCh2->Fill(CopCh2,dP*dP);
		}
	}
	for(int i=0;i<entCarbon;++i){
		treeCarbon->GetEntry(i);
		if(abs(MissMassCarbon - 1.321)<MMCut){
			double dP = dPCarbon;
			
//			double dP = (LVXi-LVMM).Mag();
			HistDMMCarbon->Fill(dP*dP);
			HistCopCarbon->Fill(CopCarbon);
			Hist2DmCarbon->Fill(CopCarbon,dP*dP);
		}
	}
	int eCh2 = HistDMMCh2->GetEffectiveEntries();
	int eCarbon = HistDMMCarbon->GetEffectiveEntries();
	double rat = eCh2 * 1./ eCarbon;
	cout<<eCh2<<" , "<<eCarbon<<endl;
	cout<<"CH2/Carbon ratio : "<<rat<<endl;
	HistDMMCarbon->Scale(rat);
	HistDMMCh2->SetLineColor(kRed);
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	HistDMMCh2->Draw("");
	gStyle->SetErrorX(0);
	HistDMMCarbon->Draw("Samehist");
	double bw = HistDMMCh2->GetBinWidth(1);
	double DMMW[39];
	double AcptW[39];
	for(int i=1;i<40;++i){
		int lb = nbin/2 + 1 - i;
		int hb = nbin/2 + 1 + i;
		int HitCh2 = GetBinContent(HistDMMCh2,lb,hb);
		int HitCarbon = GetBinContent(HistDMMCarbon,lb,hb);
		DMMW[i-1]= 2*i * bw;	
		AcptW[i-1]= {(HitCh2 -HitCarbon)* 1./HitCh2 };	
	}
	TCanvas* c2 = new TCanvas("c2","c2",1200,800);
	TGraph* g1 = new TGraph(29,DMMW,AcptW);
	g1->SetTitle("P Separation Power");
	g1->GetXaxis()->SetTitle("dP");
//	g1->GetYaxis()->SetTitle("#frac{N_{CH2}-N_{Dia}}{N_{CH2}}");
	g1->GetYaxis()->SetTitle("P Separation Power");
	g1->SetMarkerStyle(25);
	g1->SetMarkerSize(1);
	g1->Draw("AP");
	TCanvas* c3 = new TCanvas("c3","c3",1200,800);
	c3->Divide(2,1);
	c3->cd(1);
	Hist2DmCh2->Draw("colz");
	c3->cd(2);
	Hist2DmCarbon->Draw("colz");
}
*/
void Compare(){
	SetStyle();
	cout<<"ComparePolarity()"<<endl;
	cout<<"CompareCoplanarity()"<<endl;
	cout<<"CompareMissMass()"<<endl;
	cout<<"CompareMom()"<<endl;
	fileCh2 = new TFile(Form("Ch2PolaAnal.root"));
	treeCh2 = (TTree*)fileCh2->Get("tree");
	treeCh2->SetBranchAddress("Coplanarity",&CopCh2);
	treeCh2->SetBranchAddress("MissMass",&MissMassCh2);
	treeCh2->SetBranchAddress("dM",&dMCh2);
	treeCh2->SetBranchAddress("PKm_x",&PKm_xCh2);
	treeCh2->SetBranchAddress("PKm_y",&PKm_yCh2);
	treeCh2->SetBranchAddress("PKm_z",&PKm_zCh2);
	treeCh2->SetBranchAddress("PKp_x",&PKp_xCh2);
	treeCh2->SetBranchAddress("PKp_y",&PKp_yCh2);
	treeCh2->SetBranchAddress("PKp_z",&PKp_zCh2);
	treeCh2->SetBranchAddress("PXi_x",&PXi_xCh2);
	treeCh2->SetBranchAddress("PXi_y",&PXi_yCh2);
	treeCh2->SetBranchAddress("PXi_z",&PXi_zCh2);
	fileCarbon = new TFile(Form("CarbonPolaAnal.root"));
 	treeCarbon = (TTree*)fileCarbon->Get("tree");

	treeCarbon->SetBranchAddress("Coplanarity",&CopCarbon);
	treeCarbon->SetBranchAddress("MissMass",&MissMassCarbon);
	treeCarbon->SetBranchAddress("dM",&dMCarbon);
	treeCarbon->SetBranchAddress("PKm_x",&PKm_xCarbon);
	treeCarbon->SetBranchAddress("PKm_y",&PKm_yCarbon);
	treeCarbon->SetBranchAddress("PKm_z",&PKm_zCarbon);
	treeCarbon->SetBranchAddress("PKp_x",&PKp_xCarbon);
	treeCarbon->SetBranchAddress("PKp_y",&PKp_yCarbon);
	treeCarbon->SetBranchAddress("PKp_z",&PKp_zCarbon);
	treeCarbon->SetBranchAddress("PXi_x",&PXi_xCarbon);
	treeCarbon->SetBranchAddress("PXi_y",&PXi_yCarbon);
	treeCarbon->SetBranchAddress("PXi_z",&PXi_zCarbon);
}
