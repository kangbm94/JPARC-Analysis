#include "Dir.hh"
int ntTpc,runnum,evnum,nKm,nKp,nKK;
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
vector<int>* helix_flag = new vector<int>;
vector<double>* closeDist = new vector<double>;
vector<double>* chisqr = new vector<double>;
double vtx[5],vty[5],vtz[5];
vector<int>* combi_id = new vector<int>;
int ntKurama=0,ntK18=0,ntAcc=0,ntK18_=0,ntKurama_=0,ent_=0;
double MissPx[5],MissPy[5],MissPz[5];
//TFile* file = new TFile("SelectedHelix12.root");
TFile* file;
TString invname;
//TFile* file = new TFile("SelectedHelixOld.root");
TTree* tree; 
//TTree* tree = (TTree*)file->Get("tpc");
bool isAccidental(int flag){
	if(flag > 399 and flag < 500) return true;
	else return false;
}
bool isKurama(int flag){
	if(flag > 299 and flag < 400) return true;
	else return false;
}
bool isK18(int flag){
	if(flag > 199 and flag < 300) return true;
	else return false;
}
void VertexAnalyze(int i);
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
	helix_cx->clear();
	helix_cy->clear();
	helix_z0->clear();
	helix_r->clear();
	helix_dz->clear();
	helix_flag->clear();
}
/*
	 */
TLorentzVector particle(double m,TVector3 mom){
	double E = sqrt(m*m+mom.Mag2());
	return TLorentzVector(mom,E);
}
void Test2(){
	cout<<"VertexAnalyze(int runnum)"<<endl;
//	for(int i=5641;i< 5642; ++i){
	for(int i=5641;i< 5667; ++i){
		cout<<"Run0"<<i<<endl;
		VertexAnalyze(i);
	}
}
void VertexAnalyze(int runnum){
//	file = new TFile(Form("~/k18-analyzer/rootfiles/run0%d_DstTPCHelixTracking.root",runnum));
	file = new TFile(dir + Form("./run0%d_DstTPCHSKuramaSelectedHelixTracking.root",runnum));
	if(file->GetSize()< 1e5)	return;
	tree = (TTree*)file->Get("tpc");
	bool SearchPi0 = true;
	double cd_cut = 30;
	int cd_ = cd_cut;
	bool KinematicFit = false; 
	invname = Form("TPCInv0%d_cd%2d.root",runnum,cd_);
	if(KinematicFit){
		invname = Form("TPCInvKinFit0%d_cd%2d.root",runnum,cd_);
	}
	cout<<invname<<endl;
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
	tree->SetBranchAddress("helix_flag",&helix_flag);
	tree->SetBranchAddress("closeDist",&closeDist);
	tree->SetBranchAddress("mom0",&mom0);
	tree->SetBranchAddress("pid",&pid);
	tree->SetBranchAddress("charge",&charge);
	tree->SetBranchAddress("nKm",&nKm);
	tree->SetBranchAddress("nKp",&nKp);
	tree->SetBranchAddress("nKK",&nKK);
	tree->SetBranchAddress("vtx",vtx);//"vtx[nKK]/D");
	tree->SetBranchAddress("vty",vty);//"vty[nKK]/D");
	tree->SetBranchAddress("vtz",vtz);//"vtz[nKK]/D");
	tree->SetBranchAddress("MissMomx",MissPx);//"MissPx[nKK]/D");
	tree->SetBranchAddress("MissMomy",MissPy);//"MissPy[nKK]/D");
	tree->SetBranchAddress("MissMomz",MissPz);//"MissPz[nKK]/D");
	int ent = tree->GetEntries();
	TH1D* hist = new TH1D("nt","nt",20,0,20);

	
	int Xirunnum,Xievnum;
	double Ximm;
	TFile* file2 = new TFile("SelectedEvents.root");
	TTree* tree2 = (TTree*)file2->Get("tree");
//	TTree* tree2 = (TTree*)file2->Get("tpc");
	tree2->SetBranchAddress("runnum",&Xirunnum);
	tree2->SetBranchAddress("evnum",&Xievnum);
	tree2->SetBranchAddress("XiM2",&Ximm);
	int xient = tree2->GetEntries();
	int nbin=200;
	TH1D* hist2 = new TH1D("LdVertZ","LdVertZ",nbin,-250,0);
	//	TH1D* hist3 = new TH1D("VertexY","VertexY",100,-300,300);
	//	TH1D* hist3 = new TH1D("pid1","pid1",30,-1,2);
	TH1D* hist3 = new TH1D("LdIM","LdIM",nbin,1,2);
	TH1D* hist4 = new TH1D("XiIM","XiIM",nbin,1,2);
	TH1D* hist5 = new TH1D("XiIMCor","XiIMCor",nbin,1,2);
	hist5->SetLineColor(kRed);
	TF1* fgaus = new TF1("fgaus","gaus",mL-0.05,mL+0.05);
	int cd_Count=0;
	double chi_cut = 50;
	TFile* Out = new TFile(invname,"recreate");
//	TFile* Out = new TFile("TPCInvOld.root","recreate");
	TTree* outtr = new TTree("tree","tree");
	double mmpi0 = NAN;
	double inv = NAN;
	double xiinv = NAN;
	double xi0inv = NAN;
	double xiCorinv = NAN;
	double ldCorinv = NAN;
	double lp = 0;
	double lagmulti = 0;
	double pmom,pimom,ldmom,ldvtx,ldvty,ldvtz,ldpx,ldpy,ldpz,lddist,ldp;
	double ximom,xivtx,xivty,xivtz,xipx,xipy,xipz,xip;
	double xiCormom,xiCorvtx,xiCorvty,xiCorvtz,xiCorpx,xiCorpy,xiCorpz,xiCorp;
	double ldCormom,ldCorvtx,ldCorvty,ldCorvtz,ldCorpx,ldCorpy,ldCorpz,ldCorp;
	double cdLd,cdXi;
	bool Inside,InsideXi,InsideXiCor,ldflg,xiflg,xi0flg,pi0flg,ldflgInvCut,xiflgInvCut;
	
	outtr->Branch("runnum",&runnum);
	outtr->Branch("evnum",&evnum);
	outtr->Branch("ntKurama",&ntKurama);
	outtr->Branch("ntK18",&ntK18);
	outtr->Branch("ntAcc",&ntAcc);
	outtr->Branch("MM",&Ximm);
	outtr->Branch("InvMLd",&inv);
	outtr->Branch("Pmom",&pmom);
	outtr->Branch("Pimom",&pimom);
	outtr->Branch("FlgLd",&ldflg);
	outtr->Branch("FlgLdInvCut",&ldflgInvCut);
	outtr->Branch("VtxLd",&ldvtx);
	outtr->Branch("VtyLd",&ldvty);
	outtr->Branch("VtzLd",&ldvtz);
	outtr->Branch("CdLd",&cdLd);
	outtr->Branch("MomxLd",&ldpx);
	outtr->Branch("MomyLd",&ldpy);
	outtr->Branch("MomzLd",&ldpz);
	outtr->Branch("MomLd",&ldp);
	outtr->Branch("DistLd",&lddist);
	outtr->Branch("InTargetLd",&Inside);
	outtr->Branch("FlgXi",&xiflg);
	outtr->Branch("FlgXiInvCut",&xiflgInvCut);
	outtr->Branch("InvMXi",&xiinv);
	outtr->Branch("VtxXi",&xivtx);
	outtr->Branch("VtyXi",&xivty);
	outtr->Branch("VtzXi",&xivtz);
	outtr->Branch("CdXi",&cdXi);
	outtr->Branch("MomxXi",&xipx);
	outtr->Branch("MomyXi",&xipy);
	outtr->Branch("MomzXi",&xipz);
	outtr->Branch("MomXi",&xip);
	outtr->Branch("InTargetXi",&InsideXi);
	outtr->Branch("InvMXiCor",&xiCorinv);
	outtr->Branch("VtxXiCor",&xiCorvtx);
	outtr->Branch("VtyXiCor",&xiCorvty);
	outtr->Branch("VtzXiCor",&xiCorvtz);
	outtr->Branch("MomxXiCor",&xiCorpx);
	outtr->Branch("MomyXiCor",&xiCorpy);
	outtr->Branch("MomzXiCor",&xiCorpz);
	outtr->Branch("MomXiCor",&xiCorp);
	outtr->Branch("InTargetXiCor",&InsideXiCor);
	outtr->Branch("LagMulti",&lagmulti);
	outtr->Branch("InvMLdCor",&ldCorinv);
	outtr->Branch("VtxLdCor",&ldCorvtx);
	outtr->Branch("VtyLdCor",&ldCorvty);
	outtr->Branch("VtzLdCor",&ldCorvtz);
	outtr->Branch("MomxLdCor",&ldCorpx);
	outtr->Branch("MomyLdCor",&ldCorpy);
	outtr->Branch("MomzLdCor",&ldCorpz);
	outtr->Branch("MomLdCor",&ldCorp);
	outtr->Branch("mmpi0",&mmpi0);
	outtr->Branch("Flgpi0",&pi0flg);
	outtr->Branch("InvMXi0",&xi0inv);
	outtr->Branch("FlgXi0",&xi0flg);
	
	cout<<"Processing..."<<endl;
	vector<int>XiRuns;
	vector<int>XiEvs;
	vector<double>XiMM;
	for(int j=0;j<xient;j++){
		tree2->GetEntry(j);
		if(Xirunnum != runnum) continue;
		XiEvs.push_back(Xievnum);
		XiMM.push_back(Ximm);
	}
	ent = XiEvs.size();
	ntKurama_=0;
	ntK18_=0;
	ent_=0;
	for(int i=0;i<ent;++i){
		ntKurama=0;
		ntK18=0;
		ntAcc=0;
		ldflgInvCut=false;
		xiflgInvCut=false;
		if(i%10==0)cout<<i<<" th event"<<endl;
		Clear();
		int event_num = XiEvs.at(i);
		tree->GetEntry(event_num);
		Ximm = XiMM.at(i);
//		cout<<Ximm<<endl;
		bool go = true;
		if(ntTpc<1) continue;
		/*
		for(int j=0;j<xient;j++){
			Xirunnum = XiRuns.at(j);
			Xievnum = XiEvs.at(j);
			if(Xirunnum==runnum and Xievnum == evnum){
				go = true;
				break;
			}
		}
		*/
		if(!go) continue;
		vector<Vertex> verts;
		vector<Track> parts;
		vector<Track> KuramaTracks;
		for(int nt1 = 0; nt1<ntTpc;++nt1){
//			if(isBeam->at(nt1)) continue; 
			if(isAccidental(helix_flag->at(nt1))){
				ntAcc++;
				continue;
			}
			if(isK18(helix_flag->at(nt1))){
				ntK18++;
				ntK18_++;
				continue;
			}
			if(chisqr->at(nt1)>chi_cut) continue; 
			int nh = helix_cx->size();
			double hcx = helix_cx->at(nt1);
			double hcy = helix_cy->at(nt1);
			double hz0 = helix_z0->at(nt1);
			double hr = helix_r->at(nt1);
			double hdz = helix_dz->at(nt1);
			double par1[5] = {hcx,hcy,hz0,hr,hdz};
			int id1 = pid->at(nt1);
			double q1 = charge->at(nt1);
			if(isKurama(helix_flag->at(nt1))){ 
				ntKurama++;
				KuramaTracks.push_back(Track(id1,q1,par1,nt1));
				ntKurama_++;
				continue;
			}
			parts.push_back(Track(id1,q1,par1,nt1));
		}
		if(ntKurama>1){
			for(auto t : KuramaTracks){
				parts.push_back(t);
			}
		}
		int np = parts.size();
		if(np<1) continue;
		hist->Fill(np);
		for(int nt1=0;nt1<np;++nt1){
			Vertex f(parts[nt1]);
			f.SetCdCut(cd_cut);
			for(int nt2=nt1+1;nt2<np;++nt2){
				f.AddTrack(parts[nt2]);	
			}
			//if(f.NTrack()>1) verts.push_back(f);
			verts.push_back(f);
		}
		int nvt = verts.size();
		vector<Recon>LdCand;
//		vector<double>XiCand;
		for(auto vt: verts){
			vt.TrustChargeInfo(true);
			vt.SearchLdCombination();
			auto Ldc = vt.GetLd();
			LdCand.push_back(Ldc);
		}
		int nld= LdCand.size();
		Recon Ld;
		double comp = 9999;
		for(auto m:LdCand){
			if( abs(mL-m.Mass())<comp) {comp=abs(mL-m.Mass());Ld=m;}
		}
		comp = 9999;
		VertexLH V(Ld,KinematicFit);
		V.SetCdCut(cd_cut);	
		V.TrustChargeInfo(true);
		for(auto p : parts){
			V.AddTrack(p);
		}
		ldflg=Ld.Exist();

		if(ldflg)V.SearchXiCombination();	
		auto Xi = V.GetXi();
//		auto XiCor = V.GetXiCor();
		lagmulti = V.GetLagMulti();
		xiflg=Xi.Exist();
		lddist = 0;
		if(ldflg and xiflg) lddist = (Ld.Vertex()-Xi.Vertex()).Mag();
		inv=Ld.Mass();
		if(abs(inv-1.12)<0.05 and ldflg) ldflgInvCut = true;
		ldvtx=Ld.Vertex().X();
		ldvty=Ld.Vertex().Y();
		ldvtz=Ld.Vertex().Z();
		ldpx=Ld.Momentum().X();
		ldpy=Ld.Momentum().Y();
		ldpz=Ld.Momentum().Z();
		ldp=Ld.Momentum().Mag();
		Inside = InTarget(Ld.Vertex());
		cdLd = Ld.GetCD();

		xiinv=Xi.Mass();
		if(abs(xiinv-1.3)<0.1 and xiflg) xiflgInvCut = true;
		xivtx=Xi.Vertex().X();
		xivty=Xi.Vertex().Y();
		xivtz=Xi.Vertex().Z();
		xipx=Xi.Momentum().X();
		xipy=Xi.Momentum().Y();
		xipz=Xi.Momentum().Z();
		xip=Xi.Momentum().Mag();
		InsideXi = InTarget(Xi.Vertex());
		cdXi = Xi.GetCD();
		Recon LdCor;
		if(xiflg) LdCor = V.GetLd();
		ldCorinv = LdCor.Mass();
		ldCorpx=LdCor.Momentum().X();
		ldCorpy=LdCor.Momentum().Y();
		ldCorpz=LdCor.Momentum().Z();
		ldCorp=LdCor.Momentum().Mag();
		Recon XiStar;
		Recon Pi0;
		pi0flg = false;
		xi0flg = false;

		if(xiflg and SearchPi0 ){
			TVector3 XiP(MissPx[0],MissPy[0],MissPz[0]);
			double EXiS = sqrt(Ximm*Ximm+XiP.Mag2());
			TVector3 XiV(vtx[0],vty[0],-143);
			TLorentzVector XiS(XiP,EXiS);
			auto XiLV = Xi.GetLV();
			XiStar.SetLV(XiS);
			XiStar.SetVertex(XiV);
			Pi0 = Recon(XiStar,Xi,mXiStar,mXi);
			mmpi0= Pi0.Mass();
			pi0flg=true;
		}
		if(ldflg and !xiflg and SearchPi0 ){
			TVector3 XiP(MissPx[0],MissPy[0],MissPz[0]);
			double EXiS = sqrt(Ximm*Ximm+XiP.Mag2());
			TVector3 XiV(vtx[0],vty[0],-143);
			TLorentzVector XiS(XiP,EXiS);
			auto XiLV = Xi.GetLV();
			XiStar.SetLV(XiS);
			XiStar.SetVertex(XiV);
			VertexXiPi VXiPi(XiStar);
			VXiPi.SetCdCut(cd_cut);
			for(auto p : parts){
				VXiPi.AddTrack(p);
			}
			VXiPi.SearchXi0Combination();
			auto Xi0 = VXiPi.GetXi0();
			xi0flg = Xi0.Exist();
			if(xi0flg){
				xi0inv = Xi0.Mass();
				Pi0 = Recon(XiStar,Ld,mXiStar,mL);
				mmpi0= Pi0.Mass();
				pi0flg=true;
			}
//			if(abs(ximm - mXiStar)>0.1 ) xi0flg = false;
		}
		
		/*
		xiCorinv=XiCor.Mass();
		xiCorvtx=XiCor.Vertex().X();
		xiCorvty=XiCor.Vertex().Y();
		xiCorvtz=XiCor.Vertex().Z();
		xiCorpx=XiCor.Momentum().X();
		xiCorpy=XiCor.Momentum().Y();
		xiCorpz=XiCor.Momentum().Z();
		xiCorp=XiCor.Momentum().Mag();
		InsideXiCor = InTarget(XiCor.Vertex());
		*/
		hist2->Fill(ldvtz);
		if(!Inside and ldflg)hist3->Fill(inv);
		if(!InsideXi and xiflg){hist4->Fill(xiinv);}//hist5->Fill(xiCorinv);}
		outtr->Fill();
		ent_++;
	}//evt
	Out->Write();
	cout<<Form("ntK18,ntKurama, nevt (%d,%d,%d)",ntK18_,ntKurama_,ent_)<<endl; 
/*
	TCanvas* cv1 = new TCanvas("c1","c1",1200,800);
//	cv1->Divide(2,2);
	cv1->cd(1);
//	hist->Draw();
	cv1->cd(2);
//	hist2->Draw();
	cv1->cd(3);
//	hist3->Draw();
	cv1->cd(4);
//	hist4->Draw();
//	hist5->Draw("same");
*/
	hist3->Fit("fgaus","R0");
	cout<<hist2->GetEffectiveEntries()<<endl;
	cout<<cd_Count<<endl;
	double p1 = fgaus->GetParameter(0);
	double sig = fgaus->GetParameter(2);
	int cnt = p1*sig*sqrt(2*PI)*nbin;
	cout<<"NL = "<<cnt<<endl;
}

