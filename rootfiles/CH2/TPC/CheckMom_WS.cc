#include "Dir.hh"
#include "Branch.hh"
#include "Functions.hh"
#include "Utils.hh"
//TFile* file = new TFile("SelectedHelix12.root");
TFile* file;
TString invname;
//TFile* file = new TFile("SelectedHelixOld.root");
TTree* tree; 
//TTree* tree = (TTree*)file->Get("tpc");
TH1D* CloseDistsLd = new TH1D("closeDistLd","closeDistLd",1000,0,500);
TH1D* CloseDistsXi = new TH1D("closeDistXi","closeDistXi",1000,0,500);
bool isAccidental(int flag){
	if(flag > 399 and flag < 500) return true;
	else return false;
}
bool SortTrack(Track A,Track B){
	return A.GetCD()<B.GetCD();
}
void VertexAnalyze(int i);
bool DoKinematicFit = false; 
bool DoKinematicFitXi = false; 
bool Carbon = true;
/*
*/
TLorentzVector particle(double m,TVector3 mom){
	double E = sqrt(m*m+mom.Mag2());
	return TLorentzVector(mom,E);
}
double Data[100];
vector<double>RPPt;
vector<double>RPTh;
vector<double>RPPh;
vector<double>RPi1Pt;
vector<double>RPi1Th;
vector<double>RPi1Ph;
vector<double>RLdTh;
vector<double>RLdPh;
TMatrixD Cor(8,8);
double ResScale[8] = {
	//Mom	Th	Ph
	1.00,1.00,//Ld
	2.40,2.40,2.40,//P
	2.40,2.40,2.40//Pi1
};
int MinHit = 1e5;
int MaxHit = -1;
TH1D* HistPullLd[10];
TH1D* HistChi2Ld;
TH1D* HistCLLd;
TH1D* HistPullXi[10];
TH1D* HistChi2Xi;
TH1D* HistCLXi;

void CheckMom_WS(){
	cout<<"VertexAnalyze(int runnum)"<<endl;
	fstream resfile;
	resfile.open("ResParam",fstream::in);
	int status = 0;
	double MaxLd = -1;
	ReadConfLine(resfile,Data);
	while(ReadConfLine(resfile,Data)){
		if(Data[0] == -1111){
			status = 1;//dL:len
			continue;
		} if(Data[0] == -2222){
			status = 2;//dP:nHit
			continue;
		}
		if(Data[0] == -3333){ 
			status = 3;//dTh,Ph:dPP
			continue;
		}
		if(Data[0] == -4444){
			status = 4;//dTh,Ph:dPPi1
			continue;
		}
		if(Data[0] == -9999){
			status = 0;
			break;
		}
		if(status==1){
			if(Data[0]+Data[1]>MaxLd){
				MaxLd = Data[0]+Data[1];
			}
			RLdTh.push_back(Data[2]);
			RLdPh.push_back(Data[3]);
		}
		if(status==2){
			if(Data[0]<MinHit) MinHit = Data[0];
			if(Data[0]>MaxHit) MaxHit = Data[0];
			RPPt.push_back(Data[1]);
			RPi1Pt.push_back(Data[2]);
		}
		if(status==3){
			RPTh.push_back(Data[2]);
			RPPh.push_back(Data[3]);
		}
		if(status==4){
			RPi1Th.push_back(Data[2]);
			RPi1Ph.push_back(Data[3]);
		}
	}
	fstream covfile;
	covfile.open("CovMat",fstream::in);
	int Ccol = 0;
	cout<<"Loading Corelations"<<endl;
	while(ReadConfLine(covfile,Data)){
		if(Data[0] > -1){
			for(int Crow = 0; Crow < 8;++Crow){
				if(abs(Data[Crow+1])>0.0)Cor(Ccol,Crow) = Data[Crow+1];
			}
			Ccol++;
		}
	}
	TString VarLd[6] = {
		"PP","PTh","PPh",
		"PiP","PiTh","PiPh"
	};
	for(int i=0;i<6;++i){
		TString title = VarLd[i];
		HistPullLd[i] = new TH1D(title,title,100,-5,5);
	}
	HistChi2Ld = new TH1D("Chi2Ld","Chi2Ld",100,0,1);
	HistCLLd = new TH1D("ConfLLd","Chi2Ld",100,0,1);
	TString VarXi[6] = {
		"LdP","LdTh","LdPh",
		"Pi2P","Pi2Th","Pi2Ph"
	};
	for(int i=0;i<6;++i){
		TString title = VarXi[i];
		HistPullXi[i] = new TH1D(title,title,100,-5,5);
	}
	HistChi2Xi = new TH1D("Chi2Xi","Chi2Xi",100,0,1);
	HistCLXi = new TH1D("ConfLXi","Chi2Xi",100,0,1);


	//		for(int i=5641;i< 5642; ++i){
//	for(int i=5641;i< 5667; ++i){
	for(int i=5667;i< 5698; ++i){
		//	for(int i=5641;i< 5641; ++i){
		//	for(int i=5667;i< 5682; ++i){
		cout<<"Run0"<<i<<endl;
		VertexAnalyze(i);
	}
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);
	c1->Divide(3,2);
	for(int i=0;i<6;++i){
		c1->cd(i+1);
		HistPullLd[i]->Draw();
		HistPullLd[i]->Fit("gaus");
	}
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,1);
	c2->cd(1);
	HistChi2Ld->Draw();
	c2->cd(2);
	HistCLLd->Draw();
	TCanvas* c3 = new TCanvas("c3","c3",1500,700);
	c3->Divide(3,2);
	for(int i=0;i<6;++i){
		c3->cd(i+1);
		HistPullXi[i]->Draw();
		HistPullXi[i]->Fit("gaus");
	}
	TCanvas* c4 = new TCanvas("c4","c4",1200,600);
	c4->Divide(2,1);
	c4->cd(1);
	HistChi2Xi->Draw();
	c4->cd(2);
	HistCLXi->Draw();
	}
	void VertexAnalyze(int runnum){
		//	file = new TFile(Form("~/k18-analyzer/rootfiles/run0%d_DstTPCHelixTracking.root",runnum));
		//file = new TFile(dir + Form("./dstfiles_ws/run0%d_DstTPCKuramaK18Tracking.root",runnum));
		if(!Carbon)file = new TFile(dir + Form("./dstfiles_ws/run0%d_DstE42.root",runnum));
		file = new TFile(dir + Form("./dstfiles_ws/prod/run0%d_DstE42.root",runnum));
		if(file->GetSize()< 1e5)	return;
		tree = (TTree*)file->Get("tpc");
		LoadInputBranches(tree);
		bool SearchPi0 = true;
		double cd_cut = 7;
		int cd_ = cd_cut;
		//	cd_cut = 30;
		TString dir = "InvM/";
		if(Carbon) dir = "InvMProd/";
		invname = dir + Form("TPCInvM0%d_cd%d_WS.root",runnum,cd_);
		if(DoKinematicFit){
			if(DoKinematicFitXi){
				invname = dir + Form("TPCInvMKinFit0%d_cd%d_WS.root",runnum,cd_);
			}
			else{
				invname = dir + Form("TPCInvMLdKinFit0%d_cd%d_WS.root",runnum,cd_);
			}
		}
		cout<<invname<<endl;
		//	gStyle->SetOptStat(0);
		gStyle->SetOptFit(111);
		int ent = tree->GetEntries();
		TH1D* hist = new TH1D("nt","nt",20,0,20);


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
		SetOutputBranches(outtr);
		cout<<"Processing..."<<endl;
		ntKurama_=0;
		ntK18_=0;
		ent_=0;
		TMatrixD Corr = Cor.GetSub(2,7,2,7);
		for(int i=0;i<ent;++i){
			ntKurama=0;
			ntK18=0;
			ntAcc=0;
			ldflgInvCut=false;
			xiflgInvCut=false;
			if(i%10000==0)cout<<i<<" th event"<<endl;
			Clear();
			tree->GetEntry(i);
			if(MissMass->size()>0)MM = MissMass->at(0);
			else MM = -9999;
			bool go = true;
			if(ntTpc<1) continue;
			TLorentzVector MissingLV;
			TVector3 TVKPLab,TVKMLab;
			for(int ik = 0; ik < 1;++ik){
				double kmp = PKM->at(0);
				double ukm = ubTPC->at(0);
				double vkm = vbTPC->at(0);
				double nkm = 1./sqrt(ukm*ukm+vkm*vkm+1); 
				kmx = kmp*ukm;
				kmy = kmp*vkm;
				kmz = kmp*nkm;
				double kpp = PKP->at(0);
				double ukp = usTPC->at(0);
				double vkp = vsTPC->at(0);
				double nkp = 1./sqrt(ukp*ukp+vkp*vkp+1); 
				kpx = kpp*ukp;
				kpy = kpp*vkp;
				kpz = kpp*nkp;

				TLorentzVector LVKM(kmx,kmy,kmz,hypot(kmp,mk));
				TLorentzVector LVP(0,0,0,mp);
				TLorentzVector LVKP(kpx,kpy,kpz,hypot(kpp,mk));
				TLorentzVector LVMM = LVKM+LVP-LVKP;	
				MissingLV = LVMM;
				auto VMM = LVMM.Vect();
				mpx = VMM.x();
				mpy = VMM.y();
				mpz = VMM.z();
				KmP = sqrt(kmx*kmx+kmy*kmy+kmz*kmz);
				KpP = sqrt(kpx*kpx+kpy*kpy+kpz*kpz);
				KmPHelix = 0;
				KpPHelix = 0;
				KmPChisqr = 0;
				KpPChisqr = 0;
				TLorentzVector LVLab = LVKM+LVP;	
				auto FrameCM = LVLab.BoostVector();
				TLorentzVector LVCM = LVLab;
				LVLab.Boost(-FrameCM);
				TLorentzVector KPCM = LVKP,KMCM = LVKM;
				KPCM.Boost(-FrameCM);
				KMCM.Boost(-FrameCM);
				auto TVKPCM = KPCM.Vect();
				auto TVKMCM = KMCM.Vect();
				CosCM = TVKPCM*TVKMCM*(1./(TVKPCM.Mag()*TVKMCM.Mag()));		

				TVKPLab = LVKP.Vect();
				TVKMLab = LVKM.Vect();
				CosLab = TVKPLab*TVKMLab*(1./(TVKPLab.Mag()*TVKMLab.Mag()));		
			}
			vector<Vertex> verts;
			vector<Track> parts;
			vector<Track> K18Tracks;
			vector<Track> KuramaTracks;
			for(int nt1 = 0; nt1<ntTpc;++nt1){
				double p0 = mom0->at(nt1);
				double hcx = helix_cx->at(nt1);
				double hcy = helix_cy->at(nt1);
				double hz0 = helix_z0->at(nt1);
				double hr = helix_r->at(nt1);
				int nh = helix_t->at(nt1).size();
				double ht = helix_t->at(nt1).at(0);
				double hdz = helix_dz->at(nt1);
				double par1[5] = {hcx,hcy,hz0,hr,hdz};
				int id1 = pid->at(nt1);
				double q1 = charge->at(nt1);
				if(isK18->at(nt1)){
					ntK18++;
					K18Tracks.push_back(Track(id1,q1,par1,nt1,-1,nh));
					if(ntK18 == 1){
						KmPHelix = mom0->at(nt1);
						KmPChisqr = chisqr->at(nt1);
					}
					ntK18_++;
					continue;
				}
				if(isKurama->at(nt1)){ 
					ntKurama++;
					KuramaTracks.push_back(Track(id1,q1,par1,nt1,-1,nh));
					if(ntKurama == 1){
						KpPHelix = mom0->at(nt1);
						KpPChisqr = chisqr->at(nt1);
					}
					ntKurama_++;
					continue;
				}
				//			parts.push_back(Track(id1,q1,par1,nt1));
				//			parts.push_back(Track(id1,q1,par1,nt1,p0));
				parts.push_back(Track(id1,q1,par1,nt1,p0,nh,ht));
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
				vt.SearchLdCombination();
				auto Ldc = vt.GetBestCandidate();
				LdCand.push_back(Ldc);
			}
			int nld= LdCand.size();
			Recon Ld;
			double comp = 9999;
			for(auto m:LdCand){
				if( abs(mL-m.Mass())<comp) {comp=abs(mL-m.Mass());Ld=m;}
			}
			ldflg=Ld.Exist();
			Chi2Ld = -1;
			Chi2Xi = -1;
			ProbLd=-1;
			ProbXi=-1;
			double var[8];
			TMatrixD cov;
			if(ldflg and DoKinematicFit){
				int nhP = Ld.GetDaughterNhit(0);
				TVector3 TVP = Ld.GetDaughter(0).Vect();
				double inputP[3] = {TVP.x(),TVP.y(),TVP.z()};
				int inputIP[3] = {nhP,MinHit,MaxHit};
				double ResProton[3];
				vector<vector<double>>RPArr = {RPPt,RPTh,RPPh};
				GetPResolution(inputP,inputIP,RPArr,ResScale,ResProton);
				double ResPP = ResProton[0];
				double ResThP = ResProton[1];
				double ResPhP = ResProton[2];

				int nhPi1 = Ld.GetDaughterNhit(1);
				TVector3 TVPi = Ld.GetDaughter(1).Vect();
				double ResPion1[3];
				double inputPi1[3] = {TVPi.x(),TVPi.y(),TVPi.z()};
				int inputIPi1[3] = {nhPi1,MinHit,MaxHit};
				vector<vector<double>>RPi1Arr = {RPi1Pt,RPi1Th,RPi1Ph};
				GetPiResolution(inputPi1,inputIPi1,RPi1Arr,ResScale,ResPion1);
				double ResPPi1 = ResPion1[0];	
				double ResThPi1 = ResPion1[1];	
				double ResPhPi1 = ResPion1[2];	
				double var[6] = {ResPP,ResThP,ResPhP,ResPPi1,ResThPi1,ResPhPi1 };
				TMatrixD ResMat(6,6);
				for(int jj=0;jj<6;++jj){
					ResMat(jj,jj) = var[jj];
					var[jj]=var[jj]*var[jj];
				}
				auto Cov = ResMat*Corr*ResMat;


				Chi2Ld = Ld.DoKinematicFit(mL,false,var,Cov);
				ProbLd = ROOT::Math::chisquared_cdf(Chi2Ld,Ld.GetKFNDF());
				HistChi2Ld->Fill(Chi2Ld);
				HistCLLd->Fill(ProbLd);
				auto PullLd = Ld.GetKFPull();
				auto PullPP = PullLd.at(0);
				auto PullThP = PullLd.at(1);
				auto PullPhP = PullLd.at(2);
				auto PullPPi1 = PullLd.at(3);
				auto PullThPi1 = PullLd.at(4);
				auto PullPhPi1 = PullLd.at(5);
				for(int iv=0;iv<6;++iv){
					HistPullLd[iv]->Fill(PullLd.at(iv));
				}
			}
			comp = 9999;
			VertexLH V(Ld);
			V.SetCdCut(cd_cut);	
			V.TrustChargeInfo(true);
			for(auto p : parts){
				V.AddTrack(p);
			}

			if(ldflg)V.SearchXiCombination();	

			auto Xi = V.GetXi();
			//		auto XiCor = V.GetXiCor();
			xiflg=Xi.Exist();
			if(xiflg and DoKinematicFit and DoKinematicFitXi){
				int nhPi2 = Xi.GetDaughterNhit(1);
				TVector3 TVPi = Xi.GetDaughter(1).Vect();
				double ResPion2[3];
				double inputPi2[3] = {TVPi.x(),TVPi.y(),TVPi.z()};
				int inputIPi2[3] = {nhPi2,MinHit,MaxHit};
				vector<vector<double>>RPi1Arr = {RPi1Pt,RPi1Th,RPi1Ph};
				GetPiResolution(inputPi2,inputIPi2,RPi1Arr,ResScale,ResPion2);
				double ResPPi2 = ResPion2[0];	
				double ResThPi2 = ResPion2[1];	
				double ResPhPi2 = ResPion2[2];	
				auto UV = Ld.GetVariance();
				double ResPLd,ResThLd,ResPhLd;
				ResPLd = sqrt(UV(0,0));
				ResThLd = sqrt(UV(1,1));
				ResPhLd = sqrt(UV(2,2));
				double varianceXi[8] = {ResPLd,ResThLd,ResPhLd,ResPPi2,ResThPi2,ResPhPi2 };
				for(int j=0;j<6;++j){
					varianceXi[j]=varianceXi[j]*varianceXi[j];
				}
				auto CorPi2 = Cor.GetSub(5,7,5,7);
				TMatrixD ResMatPi2(3,3);
				for(int id=0;id<3;++id){
					ResMatPi2(id,id) = sqrt(varianceXi[id+3]);
				}
				auto CovMatPi2 = ResMatPi2*CorPi2*ResMatPi2;
				auto CovMatLd = MergedOffDiagonalMatrix(UV,CovMatPi2);
				Chi2Xi = Xi.DoKinematicFit(mXi,false,varianceXi,CovMatLd);
				ProbXi = ROOT::Math::chisquared_cdf(Chi2Xi,Xi.GetKFNDF());
				HistChi2Xi->Fill(Chi2Xi);
				HistCLXi->Fill(ProbXi);
				auto PullXi = Xi.GetKFPull();
				for(int iv=0;iv<6;++iv){
					HistPullXi[iv]->Fill(PullXi.at(iv));
				}
			}
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
			if(ldflg){
			auto Proton = Ld.GetDaughter(0);
			ppx = Proton.X();
			ppy = Proton.Y();
			ppz = Proton.Z();
			pp = Proton.Vect().Mag();
			}
			Inside = InTarget(Ld.Vertex());
			cdLd = Ld.GetCD();

			xiinv=Xi.Mass();
			if(abs(xiinv-1.321)<0.1 and xiflg) xiflgInvCut = true;
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
			Vtx = vtx->at(0);
			Vty = vty->at(0);
			Vtz = vtz->at(0);
			TVector3 MissingV(Vtx,Vty,Vtz);
			TVector3 MissingTV(mpx,mpy,mpz);
			double MissingE = sqrt(mm_*mm_+MissingTV.Mag2());
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
					//				KmCor = true;
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
				auto pars = t.GetPar();
				double cx = pars[0];
				double cy = pars[1];
				double cr = pars[3];
				if(  abs(PM - KuramaP.Mag())< dif and t.GetQ()==1 and abs(hypot(cx,cy)-cr) < 30){
					dif = abs(PM-KuramaP.Mag());
					KuramaDir = P*(1./PM); 
					//				KpCor = true;
				}
			}
			if(KpCor){
				if(KuramaP*KuramaDir/KuramaP.Mag() < 0.9)cout<<"KuramaCorrection Fail! cos = "<<KuramaP*KuramaDir/KuramaP.Mag()<<endl;
				else KuramaP=KuramaDir*KuramaP.Mag();
			}
			cdX = -9999;
			if(KmCor and KpCor){
				auto tk18 = K18Tracks.at(0);		
				auto tkurama = KuramaTracks.at(0);
				double t1,t2;
				auto VertX = VertexPointHelix(tk18.GetPar(),tkurama.GetPar(),cdX,t1,t2);
				VtxCor = VertX.x();
				VtyCor = VertX.y();
				VtzCor = VertX.z();

				K18Dir = CalcCircleMom(tk18.GetPar(),VertX);
				KuramaDir = CalcCircleMom(tkurama.GetPar(),VertX);
				K18Dir = K18Dir * (1./K18Dir.Mag());
				if(K18Dir.Z()<0) K18Dir = -K18Dir;
				KuramaDir = KuramaDir * (1./KuramaDir.Mag());
				if(KuramaDir.Z()<0) KuramaDir = -KuramaDir;
				K18P = K18P.Mag()*K18Dir;
				KuramaP = KuramaP.Mag()*KuramaDir;

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
			MMpKmKp = MissingLV.Mag();

			TLorentzVector pKpLdMM;
			MissingpKmKpLd = -9999;
			MissingpKmKpLdpm = -9999;
			MissingpKmKpXi = -9999;
			MissingpKmKppm = -9999;
			MissingpKmKpKm = -9999;
			TLorentzVector KpLdMM;
			auto LdLV = Ld.GetLV();

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
			TLorentzVector TPCSumCor;

			if(ldflg ){
				auto KpLdLV = LdLV + KuramaLV;
				auto InLV =  K18LV + targetLV;
				auto pKmKpLdLV = InLV - KpLdLV;
				MissingpKmKpLd = pKmKpLdLV.Mag();	

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
				TLorentzVector pmLV;
				for(auto t : Tracks){
					if(ncnt>0) continue;
					auto P = CalcHelixMom(t.GetPar(),MissingV.y());
					if(t.IsPi() and t.GetQ()==-1){
						pmLV = TLorentzVector(P,sqrt(mpi*mpi+P.Mag2()));
						TPCSum += pmLV; 
						targetcd = t.GetCD();
					}
					else ncnt--;
					ncnt ++;
				}
				//			if(ncnt==1 and targetcd < 30) {
				//				isGood = true;
				auto pKmKpLdpmLV = pKmKpLdLV - pmLV;	
				if(pmLV.Mag()!= 0) MissingpKmKpLdpm = pKmKpLdpmLV.Mag(); 
				//			}
			}
			if(ldflg and xiflg){
				TPCSum = Xi.GetLV();
				double IMXi = Xi.Mass();
				TVector3 XiPrVt(VtxCor,VtyCor,VtzCor);
				auto pXiCor = CalcCircleMom(Xi.GetPar(),XiPrVt);	
				//			if(pXiCor.z()<0){
				pXiCor = -pXiCor;
				//			}
				isGood = true;
				auto XiLvCor = TLorentzVector(pXiCor,IMXi);
				TPCSumCor = XiLvCor;
				auto MissingpKmKpXiLV = MissingLV - Xi.GetLV();	
				MissingpKmKpXi = MissingpKmKpXiLV.Mag();
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
			{
				double cd_min = 9999;
				int tr_id = -1;
				TLorentzVector pmLV;
				for(int id = 0; id < parts.size();++id){
					auto p = parts.at(id);
					if( !( p.IsPi() and p.GetQ()==-1)) continue;
					double cd_;
					cd_ = MinHelixDistance(p.GetPar(),MissingV);
					if(cd_<cd_min){
						cd_min = cd_;
						tr_id = id;
						auto pmTV = CalcHelixMom(p.GetPar(),MissingV.y());
						pmLV = TLorentzVector(pmTV,hypot(PionMass/1000,pmTV.Mag()));
					}
				}
				auto MissingpKmKppmLV = MissingLV - pmLV;
				//			cout<<"pm mass = "<<pmLV.Mag()<<endl;
				cd_pKmKppm = cd_min;
				if(cd_min < 9999){
					MissingpKmKppm=MissingpKmKppmLV.Mag();
				}
			}
			TPCSumM = TPCSum.Mag();
			TPCSumP = TPCSum.Vect().Mag();
			TPCSumPx = TPCSum.Vect().x();
			TPCSumPy = TPCSum.Vect().y();
			TPCSumPz = TPCSum.Vect().z();
			if(abs(MissingLV.Mag()-MM)>0.1){
				MissingLV = TLorentzVector(0,0,0,-9999);
			}
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
			auto CheckSumCor = MissingLV - TPCSumCor;
			CheckSumMCor = CheckSumCor.Mag();
			CheckSumPCor = CheckSumCor.Vect().Mag();
			CheckSumPxCor = CheckSumCor.Vect().x();
			CheckSumPyCor = CheckSumCor.Vect().y();
			CheckSumPzCor = CheckSumCor.Vect().z();
			TVector3 plane = K18P.Cross(KuramaP);
			plane = plane * (1./plane.Mag());
			Coplanarity = plane * TPCSum.Vect()*(1./TPCSumP);
			//		Coplanarity = acos(Coplanarity);
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

