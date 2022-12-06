#include "Math.hh"
int ntTpc,runnum,evnum;
vector<int>* pid = new vector<int>;
vector<int>* charge = new vector<int>;
vector<int>* isBeam = new vector<int>;
vector<double>* mom0 = new vector<double>;
vector<double>* mom_vtx = new vector<double>;
vector<double>* mom_vty = new vector<double>;
vector<double>* mom_vtz = new vector<double>;
vector<double>* helix_cx = new vector<double>;
vector<double>* helix_cy = new vector<double>;
vector<double>* helix_z0 = new vector<double>;
vector<double>* helix_r = new vector<double>;
vector<double>* helix_dz = new vector<double>;
vector<double>* closeDist = new vector<double>;
vector<double>* chisqr = new vector<double>;
vector<double>* vtx = new vector<double>;
vector<double>* vty = new vector<double>;
vector<double>* vtz = new vector<double>;
vector<int>* combi_id = new vector<int>;
TFile* file = new TFile("SelectedHelix.root");
TTree* tree = (TTree*)file->Get("tree");


void Clear(){
	ntTpc=0;
	isBeam->clear();
	pid->clear();
	charge->clear();
	chisqr->clear();
	mom0->clear();
	mom_vtx->clear();
	mom_vty->clear();
	mom_vtz->clear();
	closeDist->clear();
	vtx->clear();
	vty->clear();
	vtz->clear();
	helix_cx->clear();
	helix_cy->clear();
	helix_z0->clear();
	helix_r->clear();
	helix_dz->clear();
}
/*
	 */
TLorentzVector particle(double m,TVector3 mom){
	double E = sqrt(m*m+mom.Mag2());
	return TLorentzVector(mom,E);
}
void Test2(){

	gStyle->SetOptStat(0);

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("ntTpc",&ntTpc);
	tree->SetBranchAddress("isBeam",&isBeam);
	tree->SetBranchAddress("chisqr",&chisqr);
	tree->SetBranchAddress("helix_cx",&helix_cx);
	tree->SetBranchAddress("helix_cy",&helix_cy);
	tree->SetBranchAddress("helix_z0",&helix_z0);
	tree->SetBranchAddress("helix_r",&helix_r);
	tree->SetBranchAddress("helix_dz",&helix_dz);
	tree->SetBranchAddress("mom_vtx",&mom_vtx);
	tree->SetBranchAddress("mom_vty",&mom_vty);
	tree->SetBranchAddress("mom_vtz",&mom_vtz);
	tree->SetBranchAddress("vtx",&vtx);
	tree->SetBranchAddress("vty",&vty);
	tree->SetBranchAddress("vtz",&vtz);
	tree->SetBranchAddress("closeDist",&closeDist);
	tree->SetBranchAddress("mom0",&mom0);
	tree->SetBranchAddress("pid",&pid);
	tree->SetBranchAddress("charge",&charge);
	int ent = tree->GetEntries();
	TH1D* hist = new TH1D("nt","nt",20,0,20);

	
	int Xirunnum,Xievnum;
	double Ximm;
	TFile* file2 = new TFile("SelectedEvents.root");
	TTree* tree2 = (TTree*)file2->Get("tree");
	tree2->SetBranchAddress("runnum",&Xirunnum);
	tree2->SetBranchAddress("evnum",&Xievnum);
	tree2->SetBranchAddress("XiM2",&Ximm);
	int xient = tree2->GetEntries();
	int nbin=200;
	TH1D* hist2 = new TH1D("KuramaMM","KuramaMM( )",nbin,1,2);
	//	TH1D* hist3 = new TH1D("VertexY","VertexY",100,-300,300);
	//	TH1D* hist3 = new TH1D("pid1","pid1",30,-1,2);
	TH1D* hist3 = new TH1D("LdIM","LdIM",nbin,1,2);
	TH1D* hist4 = new TH1D("XiIM","XiIM",nbin,1,2);
	TF1* fgaus = new TF1("fgaus","gaus",mL-0.05,mL+0.05);
	double cd_cut = 8;
	int cd_Count=0;
	double chi_cut = 50;
	TFile* Out = new TFile("TPCInv.root","recreate");
	TTree* outtr = new TTree("tree","tree");
	double inv = NAN;
	double invXi = 0;
	double lp = 0;
	double pmom,pimom,ldmom;
	outtr->Branch("runnum",&runnum);
	outtr->Branch("evnum",&evnum);
	outtr->Branch("MM",&Ximm);
	outtr->Branch("InvMLd",&inv);
	outtr->Branch("Pmom",&pmom);
	outtr->Branch("Pimom",&pimom);
	outtr->Branch("Ldmom",&ldmom);
	outtr->Branch("InvMXi",&invXi);
	outtr->Branch("LambdaP",&lp);
	cout<<"Processing..."<<endl;
	for(int i=0;i<ent;++i){
		Clear();
		tree->GetEntry(i);
		if(ntTpc<1) continue;
		for(int j=0;j<xient;j++){
			tree2->GetEntry(j);
			if(Xirunnum==runnum and Xievnum == evnum) break;
		}
		vector<Vertex> verts;
		vector<Particle> parts;
		for(int nt1 = 0; nt1<ntTpc;++nt1){
			if(chisqr->at(nt1)>chi_cut) continue; 
			if(isBeam->at(nt1)) continue; 
			int nh = helix_cx->size();
			double hcx = helix_cx->at(nt1);
			double hcy = helix_cy->at(nt1);
			double hz0 = helix_z0->at(nt1);
			double hr = helix_r->at(nt1);
			double hdz = helix_dz->at(nt1);
			double par1[5] = {hcx,hcy,hz0,hr,hdz};
			int id1 = pid->at(nt1);
			double q1 = charge->at(nt1);
			parts.push_back(Particle(id1,q1,par1,nt1));
		}
		int np = parts.size();
		if(np<1) continue;
		hist->Fill(np);
		for(int nt1=0;nt1<np;++nt1){
			Vertex f(parts[nt1],nt1);
			for(int nt2=nt1+1;nt2<np;++nt2){
				f.AddParticle(parts[nt2],nt2);	
			}
			if(f.NParticle()>1) verts.push_back(f);
		}
		int nvt = verts.size();
		vector<Recon>LdCand;
//		vector<double>XiCand;
		for(auto vt: verts){
			vt.SearchCombination();
			auto Ldc = vt.GetLd();
//			if(Ldc.Exist()) LdCand.push_back(Ldc);
			if(Ldc.Proper()) LdCand.push_back(Ldc);
		}
		int nld= LdCand.size();
		Recon Ld;
		double comp = 9999;
		for(auto m:LdCand){
			if( abs(mL-m.Mass())<comp) {comp=abs(mL-m.Mass());Ld=m;}
		}
		comp = 9999;
		inv=Ld.Mass();
		hist3->Fill(inv);
		outtr->Fill();
	}//evt
	Out->Write();
	TCanvas* cv1 = new TCanvas("c1","c1",1200,800);
	cv1->Divide(2,2);
	cv1->cd(1);
	hist->Draw();
	cv1->cd(2);
	hist2->Draw();
	cv1->cd(3);
	hist3->Draw();
	hist3->Fit("fgaus","R");
	cv1->cd(4);
	hist4->Draw();
	//	cout<<hist2->GetEffectiveEntries()<<endl;
	cout<<cd_Count<<endl;
	double p1 = fgaus->GetParameter(0);
	double sig = fgaus->GetParameter(2);
	int cnt = p1*sig*sqrt(2*PI)*nbin;
	cout<<"NL = "<<cnt<<endl;
}

