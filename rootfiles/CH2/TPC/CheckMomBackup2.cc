#include "Dir.hh"
#include "/Users/MIN/Desktop/JPARC/E42/TPC/src/TPCManager.cc"
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
int ntKurama=0,ntK18=0,ntAcc=0;
double MissPx[5],MissPy[5],MissPz[5],MissMass[5];
double KMPX[5],KMPY[5],KMPZ[5];
double KPPX[5],KPPY[5],KPPZ[5];
int inside[5];
double mpx,mpy,mpz,mm_;
double kmx,kmy,kmz,kpx,kpy,kpz; 
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
bool SortTrack(Track A,Track B){
	return A.GetCD()<B.GetCD();
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
void CheckMom(){
	cout<<"VertexAnalyze(int runnum)"<<endl;
//	for(int i=5641;i< 5642; ++i){
	for(int i=5641;i< 5667; ++i){
		cout<<"Run0"<<i<<endl;
		VertexAnalyze(i);
	}
}
void VertexAnalyze(int runnum){
//	file = new TFile(Form("~/k18-analyzer/rootfiles/run0%d_DstTPCHelixTracking.root",runnum));
	TString filename = dir + Form("./dstfiles/run0%d_DstTPCHSKuramaSelectedHelixTracking.root",runnum);
	gTPCManager.LoadClusterFile(filename);
	file = new TFile(dir + Form("./dstfiles/run0%d_DstTPCHSKuramaSelectedHelixTracking.root",runnum));
//	file = new TFile(dir + Form("./dstfiles/before_hough_mod/run0%d_DstTPCHSKuramaSelectedHelixTracking.root",runnum));
//	file = new TFile(dir + Form("./dstfiles/HelixFitNoItr/run0%d_DstTPCHSKuramaSelectedHelixTracking.root",runnum));
//	if(file->GetSize()< 1e5)	return;
	bool SearchPi0 = true;
	double cd_cut = 15;
	int cd_ = cd_cut;
	bool KinematicFit = false; 
	int ent_=0,ntK18_=0,ntKurama_=0,nXi_=0;
	int nP = 0,nPi=0,nt=0;
	invname = Form("TPCInvM0%d_cd%2d.root",runnum,cd_);
	if(KinematicFit){
		invname = Form("TPCInvMKinFit0%d_cd%2d.root",runnum,cd_);
	}
	cout<<invname<<endl;
	gStyle->SetOptStat(0);
	tree->SetBranchAddress("pid",&pid);
	tree->SetBranchAddress("charge",&charge);
	tree->SetBranchAddress("nKm",&nKm);
	tree->SetBranchAddress("nKp",&nKp);
	tree->SetBranchAddress("nKK",&nKK);
	tree->SetBranchAddress("vtx",vtx);//"vtx[nKK]/D");
	tree->SetBranchAddress("vty",vty);//"vty[nKK]/D");
	tree->SetBranchAddress("vtz",vtz);//"vtz[nKK]/D");
	tree->SetBranchAddress("inside",inside);
	tree->SetBranchAddress("KMPX",KMPX);//"MissPx[nKK]/D");
	tree->SetBranchAddress("KMPY",KMPY);//"MissPx[nKK]/D");
	tree->SetBranchAddress("KMPZ",KMPZ);//"MissPx[nKK]/D");
	tree->SetBranchAddress("KPPX",KPPX);//"MissPx[nKK]/D");
	tree->SetBranchAddress("KPPY",KPPY);//"MissPx[nKK]/D");
	tree->SetBranchAddress("KPPZ",KPPZ);//"MissPx[nKK]/D");
	tree->SetBranchAddress("MissMass",MissMass);//"MissPx[nKK]/D");
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
	bool Inside,InsideXi,InsideXiCor,ldflg,xiflg,xi0flg,pi0flg,ldflgInvCut,xiflgInvCut,isGood;
	
	double TPCSumM,TPCSumP,TPCSumPx,TPCSumPy,TPCSumPz;
	double MissingM,MissingP,MissingPx,MissingPy,MissingPz;
	double KmP,KmPx,KmPy,KmPz;
	double KpP,KpPx,KpPy,KpPz;
	bool KmCor,KpCor;
	double CheckSumM,CheckSumP,CheckSumPx,CheckSumPy,CheckSumPz;
	double CosXT;
	double MissingMCor,MissingPCor,MissingPxCor,MissingPyCor,MissingPzCor;
	double CheckSumMCor,CheckSumPCor,CheckSumPxCor,CheckSumPyCor,CheckSumPzCor;
	double Coplanarity = 0;
	outtr->Branch("runnum",&runnum);
	outtr->Branch("evnum",&evnum);
	outtr->Branch("nP",&nP);
	outtr->Branch("nPi",&nPi);
	outtr->Branch("nt",&nt);
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
	outtr->Branch("isGood",&isGood);
	outtr->Branch("MissMassCor",&MissMass[0]);//"MissPx[nKK]/D");
	outtr->Branch("MissMomxCor",&MissPx[0]);//"MissPx[nKK]/D");
	outtr->Branch("MissMomyCor",&MissPy[0]);//"MissPy[nKK]/D");
	outtr->Branch("MissMomzCor",&MissPz[0]);//"MissPz[nKK]/D");
	outtr->Branch("TPCSumM",&TPCSumM);
	outtr->Branch("TPCSumP",&TPCSumP);
	outtr->Branch("TPCSumPx",&TPCSumPx);
	outtr->Branch("TPCSumPy",&TPCSumPy);
	outtr->Branch("TPCSumPz",&TPCSumPz);
	outtr->Branch("MissingM",&MissingM);
	outtr->Branch("MissingP",&MissingP);
	outtr->Branch("MissingPx",&MissingPx);
	outtr->Branch("MissingPy",&MissingPy);
	outtr->Branch("MissingPz",&MissingPz);
	outtr->Branch("KmCor",&KmCor);
	outtr->Branch("KmP",&KmP);
	outtr->Branch("KmPx",&KmPx);
	outtr->Branch("KmPy",&KmPy);
	outtr->Branch("KmPz",&KmPz);
	outtr->Branch("KpCor",&KpCor);
	outtr->Branch("KpP",&KpP);
	outtr->Branch("KpPx",&KpPx);
	outtr->Branch("KpPy",&KpPy);
	outtr->Branch("KpPz",&KpPz);
	outtr->Branch("CheckSumM",&CheckSumM);
	outtr->Branch("CheckSumP",&CheckSumP);
	outtr->Branch("CheckSumPx",&CheckSumPx);
	outtr->Branch("CheckSumPy",&CheckSumPy);
	outtr->Branch("CheckSumPz",&CheckSumPz);
	outtr->Branch("CosXT",&CosXT);
	outtr->Branch("MissingMCor",&MissingMCor);
	outtr->Branch("MissingPCor",&MissingPCor);
	outtr->Branch("MissingPxCor",&MissingPxCor);
	outtr->Branch("MissingPyCor",&MissingPyCor);
	outtr->Branch("MissingPzCor",&MissingPzCor);
	outtr->Branch("CheckSumMCor",&CheckSumMCor);
	outtr->Branch("CheckSumPCor",&CheckSumPCor);
	outtr->Branch("CheckSumPxCor",&CheckSumPxCor);
	outtr->Branch("CheckSumPyCor",&CheckSumPyCor);
	outtr->Branch("CheckSumPzCor",&CheckSumPzCor);
	outtr->Branch("Coplanarity",&Coplanarity);

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
		if(i%100==0)cout<<i<<" th event"<<endl;
		Clear();
		int event_num = XiEvs.at(i);
//		tree->GetEntry(event_num);
		gTPCManager.SetEvent(event_num);
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
		for(int ik = 0; ik < nKK;++ik){
			if(inside[ik]){
				kmx = KMPX[ik];
				kmy = KMPY[ik];
				kmz = KMPZ[ik];
				kpx = KPPX[ik];
				kpy = KPPY[ik];
				kpz = KPPZ[ik];
				mpx = MissPx[ik];
				mpy = MissPy[ik];
				mpz = MissPz[ik];
				mm_ = MissMass[ik];
			}
		}
		vector<Vertex> verts;
		vector<Track> parts;
		vector<Track> K18Tracks;
		vector<Track> KuramaTracks;
		for(int nt1 = 0; nt1<ntTpc;++nt1){
//			if(isBeam->at(nt1)) continue; 
			if(isAccidental(helix_flag->at(nt1))){
				ntAcc++;
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
			if(isK18(helix_flag->at(nt1))){
				ntK18++;
				K18Tracks.push_back(Track(id1,q1,par1,nt1));
				ntK18_++;
				continue;
			}
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
		nt = parts.size();
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
		nP = 0;
		nPi= 0;
		for(auto p:parts){
			if(p.IsP())nP++;
			if(p.IsPi())nPi++;
		}
		int nvt = verts.size();
		vector<Recon>LdCand;
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
		if(xiflg){
			LdCor = V.GetLd();
			nXi_++;
		}
		ldCorinv = LdCor.Mass();
		ldCorpx=LdCor.Momentum().X();
		ldCorpy=LdCor.Momentum().Y();
		ldCorpz=LdCor.Momentum().Z();
		ldCorp=LdCor.Momentum().Mag();
		Recon XiStar;
		Recon Pi0;
		pi0flg = false;
		xi0flg = false;
		TVector3 MissingV(vtx[0],vty[0],-143);
		TVector3 MissingTV(mpx,mpy,mpz);
		double MissingE = sqrt(mm_*mm_+MissingTV.Mag2());
		TLorentzVector MissingLV(MissingTV,MissingE);
		TVector3 K18P(kmx,kmy,kmz);
		TVector3 K18Dir;
		TVector3 KuramaDir;
		KmCor = false;
		double dif = 999;
		for(auto t:K18Tracks){
			auto P = CalcCircleMom(t.GetPar(),MissingV);
			if(P.z()<0) P=-P;
			double PM = P.Mag();
			if(  abs(PM - K18P.Mag())< dif ){
				dif = abs(PM-K18P.Mag());
				K18Dir = P*(1./PM); 
				KmCor = true;
			}
		}
		if(KmCor){
				if(K18P*K18Dir/K18P.Mag() < 0.9)cout<<"K18Correction Fail! cos = "<<K18P*K18Dir/K18P.Mag()<<endl;
				else K18P = K18Dir*K18P.Mag();
		}
		
		dif = 999;
		KpCor = false;
		TVector3 KuramaP(kpx,kpy,kpz);
		for(auto t:KuramaTracks){
			auto P = CalcCircleMom(t.GetPar(),MissingV);
			if(P.z()<0) P=-P;
			double PM = P.Mag();
			if(  abs(PM - KuramaP.Mag())< dif and t.GetQ()==1){
				dif = abs(PM-KuramaP.Mag());
				KuramaDir = P*(1./PM); 
				KpCor = true;
			}
		}
		if(KpCor){
			if(KuramaP*KuramaDir/KuramaP.Mag() < 0.9)cout<<"KuramaCorrection Fail! cos = "<<KuramaP*KuramaDir/KuramaP.Mag()<<endl;
			else KuramaP=KuramaDir*KuramaP.Mag();
		}
		double MissMP = MissingTV.Mag();
		TLorentzVector targetLV = TLorentzVector(0,0,0,mp);
		TLorentzVector K18LV(K18P,sqrt(K18P.Mag2()+mk*mk));
		TLorentzVector KuramaLV(KuramaP,sqrt(KuramaP.Mag2()+mk*mk));
		MissingLV = targetLV + K18LV - KuramaLV;
		if(xiflg and SearchPi0 ){
			auto XiLV = Xi.GetLV();
			XiStar.SetLV(MissingLV);
			XiStar.SetVertex(MissingV);
			Pi0 = Recon(XiStar,Xi,mXiStar,mXi);
			mmpi0= Pi0.Mass();
			pi0flg=true;
		}
		if(ldflg and !xiflg and SearchPi0 ){

			auto XiLV = Xi.GetLV();
			XiStar.SetLV(MissingLV);
			XiStar.SetVertex(MissingV);
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
		isGood = false;
		TLorentzVector TPCSum;
		if(ldflg and xiflg){
			TPCSum = Xi.GetLV();
			isGood = true;
		}
		if(ldflg and !xiflg){
			int pid1 = Ld.GetID1();
			int pid2 = Ld.GetID2();
			TPCSum = Ld.GetLV();	
			vector<double>cd;
			vector<Track>Tracks;
			double targetcd = 9999;
			for(auto p:parts){
				int id = p.GetID();
				if(id == pid1 or id == pid2) continue;
				double cd_;
				cd_ = MinHelixDistance(p.GetPar(),MissingV);
				p.SetCD(cd_);
				Tracks.push_back(p);
			}
			sort(Tracks.begin(),Tracks.end(),SortTrack);
			int ncnt = 0;
			for(auto t : Tracks){
				if(ncnt>0) continue;
				auto P = CalcHelixMom(t.GetPar(),MissingV.y());
				if(t.IsPi() and t.GetQ()==-1){
					TPCSum += TLorentzVector(P,sqrt(mpi*mpi+P.Mag2()));
					targetcd = t.GetCD();
				}
				else ncnt--;
				ncnt ++;
			}
//			if(ncnt==1 and targetcd < 30) isGood = true;
		}
		if (!ldflg){
			vector<double>cd;
			vector<Track>Tracks;
			double targetcd[3] = {9999,9999,9999};
			for(auto p:parts){
				int id = p.GetID();
				double cd_;
				cd_ = MinHelixDistance(p.GetPar(),MissingV);
				p.SetCD(cd_);
				Tracks.push_back(p);
			}
			sort(Tracks.begin(),Tracks.end(),SortTrack);
			int ncnt = 0;
			int npi=0;
			int np=0;
			for(auto t : Tracks){
				if(ncnt>3 or (npi == 2 and np == 1)) continue;
				auto P = CalcHelixMom(t.GetPar(),MissingV.y());
				if(t.IsPi() and t.GetQ()==-1 and npi < 2){
					targetcd[npi] = t.GetCD();
					TPCSum += TLorentzVector(P,sqrt(mpi*mpi+P.Mag2()));	
					npi++;
				}
				else if(t.IsP() and t.GetQ()==1 and np < 1){
					TPCSum += TLorentzVector(P,sqrt(mp*mp+P.Mag2()));	
					targetcd[np+2] = t.GetCD();
					np++;
				}
				else{
					ncnt--;
				}
				ncnt ++;
			}
			bool vcut = (targetcd[2]< 30 and targetcd[1]<30);
//			if(npi == 2 and np == 1 and vcut)isGood = true;
		}
		TPCSumM = TPCSum.Mag();
		TPCSumP = TPCSum.Vect().Mag();
		TPCSumPx = TPCSum.Vect().x();
		TPCSumPy = TPCSum.Vect().y();
		TPCSumPz = TPCSum.Vect().z();
		MissingM = MissingLV.Mag();		
		MissingP = MissingLV.Vect().Mag();		
		MissingPx = MissingLV.Vect().x();		
		MissingPy = MissingLV.Vect().y();		
		MissingPz = MissingLV.Vect().z();		
		auto CheckSum = MissingLV - TPCSum;
		CheckSumM = CheckSum.Mag();	
		CheckSumP = CheckSum.Vect().Mag();
		CheckSumPx = CheckSum.Vect().x();
		CheckSumPy = CheckSum.Vect().y();
		CheckSumPz = CheckSum.Vect().z();
		CosXT = (TPCSum.Vect()*MissingLV.Vect())/TPCSumP/MissingP;
		auto MissingTVCor =TPCSum.Vect()*(MissingP/TPCSumP);	
		TLorentzVector MissingLVCor(MissingTVCor,sqrt(MissingM*MissingM+MissingTVCor.Mag2()));
		MissingMCor = MissingLVCor.Mag();
		MissingPCor = MissingLVCor.Vect().Mag();
		MissingPxCor = MissingLVCor.Vect().x();
		MissingPyCor = MissingLVCor.Vect().y();
		MissingPzCor = MissingLVCor.Vect().z();
		auto CheckSumCor = MissingLVCor - TPCSum;
		CheckSumMCor = CheckSumCor.Mag();
		CheckSumPCor = CheckSumCor.Vect().Mag();
		CheckSumPxCor = CheckSumCor.Vect().x();
		CheckSumPyCor = CheckSumCor.Vect().y();
		CheckSumPzCor = CheckSumCor.Vect().z();
		TVector3 plane = K18P.Cross(KuramaP);
		plane = plane * (1./plane.Mag());
		Coplanarity = plane * TPCSum.Vect()*(1./TPCSumP);
		Coplanarity = acos(Coplanarity);
		//		Coplanarity = 

		//		hist2->Fill(ldvtz);
		if(!Inside and ldflg)hist3->Fill(inv);
//		if(!InsideXi and xiflg){hist4->Fill(xiinv);}//hist5->Fill(xiCorinv);}
		outtr->Fill();
		ent_++;
	}//evt
	Out->Write();
	cout<<Form("ntK18,ntKurama,nXi, nevt (%d,%d,%d,%d)",ntK18_,ntKurama_,nXi_,ent_)<<endl; 
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

