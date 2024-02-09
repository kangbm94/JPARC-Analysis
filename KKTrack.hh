#include "PhysicalConstants.hh"
#include "KKBranch.hh"
#ifndef KKTrack_hh
#define KKTrack_hh
double PiUpperThreshold=0.12;
double KLowerThreshold= 0.15;
double KUpperThreshold= 0.4;
double PLowerThreshold= 0.5;
double PUpperThreshold= 1.8;

static double K18Cut=20;
static double KuramaCut=200;

class KKBeam{
	protected:
		double U_,V_,P_,ChiSqr_,Norm_;
	public:
		KKBeam(){};
		KKBeam(double u,double v,double p,double chisqr){
			U_=u;V_=v;P_=p;ChiSqr_=chisqr;Norm_=sqrt(1+u*u+v*v);
		}
		double P(){
			return P_;
		}
		double Px(){
			return P_*U_/Norm_;
		}
		double Py(){
			return P_*V_/Norm_;
		}
		double Pz(){
			return P_/Norm_;
		}
		bool CutChiSqr(double chi2){
			return ChiSqr_<chi2;
		}
		double GetU(){
			return U_;
		}
		double GetV(){
			return U_;
		}
		double GetMomentum(){
			return P_;
		}
};
class KKScat: public KKBeam{
	protected:
		double M2_;
		double Q_;
		//		double Vtxx_,Vtxy_,Vtxz_;
	public:
		KKScat(){};
		KKScat(double u,double v,double p,double chisqr,double m2,double q):KKBeam(u,v,p,chisqr){
			M2_=m2;Q_=q;
		};
		int ParticleID();
		double En(){
			return sqrt(P_*P_+M2_);
		}
		TLorentzVector FourVector(){
			return TLorentzVector(Px(),Py(),Pz(),En());
		};
		bool CutM2(double m2_min,double m2_max){
			return m2_min<M2_&&M2_<m2_max;
		}
		double GetM2(){
			return M2_;
		}
		double GetCharge(){
			return Q_;
		}
		bool CutMomentum(double P_cut){
//			return P_<P_cut and 1.1<P_;
			return P_<P_cut;
		}
		bool CutCharge(double Q){
			return Q+Q_;//qKurama always should be 1 or -1;
		}
};
int KKScat::ParticleID(){
	int ID=0;
	if(M2_<PiUpperThreshold){
		ID=PionID;	
	}
	else if(KLowerThreshold<M2_&&M2_<KUpperThreshold){
		ID=KaonID;
	}else if(PLowerThreshold<M2_&&M2_<PUpperThreshold){
		ID=ProtonID;
	}
	else{
		ID=0;
	}
	return Q_*ID;
}
class KKTrack{
	private:
		KKBeam KM_;
		KKScat KP_;
		bool Inside_;
		double MissMass_;
		double Theta_;
		double ThetaCM_;
		double MissMassCalc_;
	public:
		KKTrack(KKBeam km,KKScat kp,bool inside,double missmass,double theta,double thetaCM){
			KM_=km;KP_=kp;Inside_=inside;MissMass_=missmass;Theta_=theta;ThetaCM_=thetaCM;
			auto KMLV = TLorentzVector(KM_.Px(),KM_.Py(),KM_.Pz(),sqrt(KaonMass*KaonMass/1000./1000.+KM_.P()*KM_.P()));
			auto PLV = TLorentzVector(0,0,0,ProtonMass/1000);
			auto KPLV = TLorentzVector(KP_.Px(),KP_.Py(),KP_.Pz(),sqrt(KaonMass*KaonMass/1000./1000.+KP_.P()*KP_.P()));
			auto XLV = KMLV +PLV - KPLV;
			MissMassCalc_ = XLV.Mag();
		}
		bool CutChiSqr(double k18cut,double kuramacut){
			return KM_.CutChiSqr(k18cut)&&KP_.CutChiSqr(kuramacut); };
		bool CutVertex(){
			return Inside_;
		};
		bool CutVertex(TVector3 Pos,TVector3 Size){
			return Inside_;
		};
		bool CutMomentum(double P_cut){
			return KP_.CutMomentum(P_cut);
		};
		bool CutCharge(double Q){
			return KP_.CutCharge(Q);
		};
		bool CutM2(double m2_min,double m2_max){
			return KP_.CutM2(m2_min,m2_max);
		}
		double GetMissMass(){
			return MissMass_;
		}
		double GetMissMassCalc(){
			return MissMassCalc_;
		}
		double GetTheta(){
			return Theta_;
		}
		double GetThetaCM(){
			return ThetaCM_;
		}
		double GetBeamMomentum(){
			return KM_.GetMomentum();
		}
		double GetMomentum(){
			return KP_.GetMomentum();
		}
		double GetCharge(){
			return KP_.GetCharge();
		}
		double GetM2(){
			return KP_.GetM2();
		}
		double GetUKP(){
			return KP_.GetU();
		}
		double GetVKP(){
			return KP_.GetV();
		}
};












class KKEvent{
	private:	
		TChain* kkChain;
		int nKm_,nKp_,nKK_,ntKurama_,ntK18_,runnum,evnum;
		int inside_[25];
		int trigpat[32];
		double pKurama_[25];double qKurama_[25];double xkp_[25];double ykp_[25];double ukp_[25];double vkp_[25];double pOrg_[25];double m2_[25];double pK18_[25];double utgtK18_[25];double vtgtK18_[25];double ukm_[25];double vkm_[25];double vtz_[25];double vtx_[25];double vty_[25];double MissMass_[25];double theta_[25];double thetaCM_[25];double closeDist_[25];
		double chisqrK18_[25];
		double chisqrKurama_[25];
		int nkk=0;

	private:
		int evt_num;
		vector<KKTrack> Tracks;
	public:
		KKEvent(TChain* chain);
		void LoadEvent(int event);
		void Clear(){
			Tracks.clear();nkk=0;
		};
		void CutChiSqr(double k18cut=K18Cut,double kuramacut=KuramaCut);
		void CutTrig(int trig);	
		void CutVertex();	
		void CutMomentum(double P_cut=1.4);	
		void CutCharge(double Q=1);	
		void CutM2(double m2_min = 0.15,double m2_max=0.4);	
		int GetRunNum(){return runnum;}
		int GetEvNum(){return evnum;}
		int Getnkk(){
			return nkk;
		}
		double GetMissMass(int ikk){
			return Tracks[ikk].GetMissMass();
		}
		double GetMissMassCalc(int ikk){
			return Tracks[ikk].GetMissMassCalc();
		}
		
		double GetTheta(int ikk){
			return Tracks[ikk].GetTheta();
		}
		double GetThetaCM(int ikk){
			return Tracks[ikk].GetThetaCM();
		}
		double GetUKP(int ikk){
			return Tracks[ikk].GetUKP();
		}
		double GetVKP(int ikk){
			return Tracks[ikk].GetVKP();
		}
		double GetMomentum(int ikk){
			return Tracks[ikk].GetMomentum();
		}
		double GetM2(int ikk){
			return Tracks[ikk].GetM2();
		}
		double GetCharge(int ikk){
			return Tracks[ikk].GetCharge();
		}
};

void KKEvent::CutChiSqr(double k18cut=K18Cut,double kuramacut=KuramaCut){
	vector<KKTrack> TrackCut;
	int nkk_cut=0;
	for(int ikk=0;ikk<nkk;++ikk){
		KKTrack Track = Tracks[ikk];
		if(Track.CutChiSqr(k18cut,kuramacut)){
				TrackCut.push_back(Track);
				nkk_cut++;
		}
	}
	Tracks.clear();
	for(int i=0;i<nkk_cut;++i){
		Tracks.push_back(TrackCut[i]);
	}
	nkk=nkk_cut;
}
	
void KKEvent::CutTrig(int trig){
	vector<KKTrack> TrackCut;
	int nkk_cut=0;
	for(int ikk=0;ikk<nkk;++ikk){
		if(trigpat[trig]>0){	
			KKTrack Track = Tracks[ikk];
			TrackCut.push_back(Track);
			nkk_cut++;
		}
	}
	Tracks.clear();
	for(int i=0;i<nkk_cut;++i){
		Tracks.push_back(TrackCut[i]);
	}
	nkk=nkk_cut;
}

void KKEvent::CutVertex(){	vector<KKTrack> TrackCut;
	int nkk_cut=0;
	for(int ikk=0;ikk<nkk;++ikk){
		KKTrack Track = Tracks[ikk];
		if(Track.CutVertex()){
				TrackCut.push_back(Track);
				nkk_cut++;
		}
	}
	Tracks.clear();
	for(int i=0;i<nkk_cut;++i){
		Tracks.push_back(TrackCut[i]);
	}
	nkk=nkk_cut;
}

void KKEvent::CutMomentum(double P_cut=1.4){
	vector<KKTrack> TrackCut;
	int nkk_cut=0;
	for(int ikk=0;ikk<nkk;++ikk){
		KKTrack Track = Tracks[ikk];
		if(Track.CutMomentum(P_cut)){
				TrackCut.push_back(Track);
				nkk_cut++;
		}
	}
	Tracks.clear();
	for(int i=0;i<nkk_cut;++i){
		Tracks.push_back(TrackCut[i]);
	}
	nkk=nkk_cut;
}
void KKEvent::CutCharge(double Q=1){
	vector<KKTrack> TrackCut;
	int nkk_cut=0;
	for(int ikk=0;ikk<nkk;++ikk){
		KKTrack Track = Tracks[ikk];
		if(Track.CutCharge(Q)){
				TrackCut.push_back(Track);
				nkk_cut++;
		}
	}
	Tracks.clear();
	for(int i=0;i<nkk_cut;++i){
		Tracks.push_back(TrackCut[i]);
	}
	nkk=nkk_cut;
}

void KKEvent::CutM2(double m2_min = 0.15,double m2_max=0.4){
	vector<KKTrack> TrackCut;
	int nkk_cut=0;
	for(int ikk=0;ikk<nkk;++ikk){
		KKTrack Track = Tracks[ikk];
		if(Track.CutM2(m2_min,m2_max)){
				TrackCut.push_back(Track);
				nkk_cut++;
		}
	}
	Tracks.clear();
	for(int i=0;i<nkk_cut;++i){
		Tracks.push_back(TrackCut[i]);
	}
	nkk=nkk_cut;
}














KKEvent::KKEvent(TChain* chain){
	kkChain = chain ;
	cout<<"Loading Chain..."<<endl;
	kkChain->SetBranchAddress("evnum",&evnum);
	kkChain->SetBranchAddress("runnum",&runnum);
	kkChain->SetBranchAddress("trigpat",trigpat);
	kkChain->SetBranchAddress("ntK18",&ntK18_);
	kkChain->SetBranchAddress("pK18",pK18_);
	kkChain->SetBranchAddress("chisqrK18",chisqrK18_);
	kkChain->SetBranchAddress("ntKurama",&ntKurama_);
	kkChain->SetBranchAddress("chisqrKurama",chisqrKurama_);
	kkChain->SetBranchAddress("nKm",&nKm_);
	kkChain->SetBranchAddress("nKp",&nKp_);
	kkChain->SetBranchAddress("nKK",&nKK_);
	kkChain->SetBranchAddress("m2",m2_);
	kkChain->SetBranchAddress("pKurama",pKurama_);
	kkChain->SetBranchAddress("qKurama",qKurama_);
	kkChain->SetBranchAddress("vtx",vtx_);
	kkChain->SetBranchAddress("vty",vty_);
	kkChain->SetBranchAddress("vtz",vtz_);
	kkChain->SetBranchAddress("closeDist",closeDist_);
	kkChain->SetBranchAddress("pOrg",pOrg_);
	kkChain->SetBranchAddress("ukp",ukp_);
	kkChain->SetBranchAddress("vkp",vkp_);
	kkChain->SetBranchAddress("ukm",ukm_);
	kkChain->SetBranchAddress("vkm",vkm_);
	kkChain->SetBranchAddress("MissMassCorr",MissMass_);
	kkChain->SetBranchAddress("theta",theta_);
	kkChain->SetBranchAddress("thetaCM",thetaCM_);
	kkChain->SetBranchAddress("inside",inside_);
}


void KKEvent::LoadEvent(int event){
	kkChain->GetEntry(event);
	double ukm_offset = 0.1;
	for(int ikp=0;ikp<nKp_;++ikp){
		for(int ikm=0;ikm<nKm_;++ikm){
			KKBeam Beam(ukm_[nkk]+ukm_offset,vkm_[nkk],pK18_[nkk],chisqrK18_[nkk]);
			KKScat Scat(ukp_[ikp],vkp_[ikp],pKurama_[ikp],chisqrKurama_[ikp],m2_[ikp],qKurama_[ikp]);
			KKTrack Track(Beam,Scat,inside_[nkk],MissMass_[nkk],theta_[nkk],thetaCM_[nkk]);
			Tracks.push_back(Track);
			nkk++;
		}
	}
}












class Track{
	private:

		double BeamU;
		double BeamV;
		double BeamP;
		double BeamNorm;

		double M2;
		double P;
		double Q;
		double U;
		double V;
		double Norm;
		int Inside;
		double ChisqrK18;
		double ChisqrKurama;
	public:
		Track(){};
		Track(double m2, double p, double q,double u,double v,int inside,double chisqrK18,double chisqrKurama){
			M2=m2;
			P=p;
			Q=q;
			U=u;
			V=v;
			Norm=sqrt(1+u*u+v*v);
			Inside = inside;
			ChisqrK18=chisqrK18;
			ChisqrKurama=chisqrKurama;
		};
		bool IsInside(){
			return Inside;
		}
		TLorentzVector GetFourVector(){
			return TLorentzVector(U*P/Norm,V*P/Norm,P/Norm,P*P+M2);
		};
		int ParticleID();
		double GetPx(){
			return P*U/Norm;
		}
		double GetPy(){
			return P*V/Norm;
		}
		double GetPz(){
			return P/Norm;
		}
		double GetMomentum(){
			return P;
		}
		double GetKaonE(){
			return sqrt(P*P+KaonMass/1000*KaonMass/1000);
		}
		double GetE(){
			return sqrt(P*P+M2);
		}
		double GetM2(){
			return M2;
		}
		double GetQM2(){
			return Q*M2;
		}
		int GetCharge(){
			return Q;
		}
		int CutChiSqr(double ChiCutK18,double ChiCutKurama){
			if(ChisqrK18<ChiCutK18&&ChisqrKurama<ChiCutKurama){
				return 1;
			}
			else{
				return 0;
			}
		}
		void SetBeam(double u,double v, double p){
			BeamU=u;BeamV=v;BeamP=p;BeamNorm=sqrt(1+u*u+v*v);
		}

		double GetBeamPx(){
			return BeamP*BeamU/BeamNorm;
		};
		double GetBeamPy(){
			return BeamP*BeamV/BeamNorm;
		};
		double GetBeamPz(){
			return BeamP/BeamNorm;
		};
		double GetBeamMomentum(){
			return BeamP;
		};
};
int Track::ParticleID(){
	int ID=0;
	if(M2<PiUpperThreshold){
		ID=PionID;	
	}
	else if(KLowerThreshold<M2&&M2<KUpperThreshold){
		ID=KaonID;
	}else if(PLowerThreshold<M2&&M2<PUpperThreshold){
		ID=ProtonID;
	}
	else{
		ID=0;
	}
	return Q*ID;
}
#endif
