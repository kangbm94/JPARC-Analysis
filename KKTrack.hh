#include "PhysicalConstants.hh"
#include "KKBranch.hh"
double PiUpperThreshold=0.12;
double KLowerThreshold= 0.15;
double KUpperThreshold= 0.4;
double PLowerThreshold= 0.5;
double PUpperThreshold= 1.8;
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
};
class KKScat: public KKBeam{
	protected:
		double M2_;
		double Q_;
		int Inside_;
		//		double Vtxx_,Vtxy_,Vtxz_;
	public:
		KKScat(){};
		KKScat(double u,double v,double p,double chisqr,double m2,double q,int inside):KKBeam(u,v,p,chisqr){
			M2_=m2;Q_=q;Inside_=inside;
		};
		int ParticleID();
		double En(){
			return sqrt(P_*P_+M2_);
		}
		double GetKaonE(){
			return sqrt(P_*P_+KaonMass/1000*KaonMass/1000);
		}
		TLorentzVector FourVector(){
			return TLorentzVector(Px(),Py(),Pz(),En());
		};
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
	public:
		KKTrack(KKBeam km,KKScat kp){
			KM_=km;KP_=kp;
		}
		bool CutChiSqr(double k18cut,double kuramacut){
			return KM_.CutChiSqr(k18cut)&&KP_.CutChiSqr(kuramacut);
		};
};


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








int nKaon;
double BeamPx[25];
double BeamPy[25];
double BeamPz[25];
double BeamP[25];

double KPx[25];
double KPy[25];
double KPz[25];
double KP[25];
double KM2[25];
int KQ[25];



class KKEvent{
	private:	
		TChain* kkChain;
		int nKm_,nKp_,nKK_,ntKurama_,ntK18_;
		int inside_[25];
		double pKurama_[25];double qKurama_[25];double xkp_[25];double ykp_[25];double ukp_[25];double vkp_[25];double pOrg_[25];double m2_[25];double pK18_[25];double utgtK18_[25];double vtgtK18_[25];double ukm_[25];double vkm_[25];double vtz_[25];double vtx_[25];double vty_[25];
		double chisqrK18_[25];
		double chisqrKurama_[25];
		int nkk=0;

	private:
		double K18Cut=20;
		double KuramaCut=200;
		int evt_num;
		vector<KKTrack> Tracks;
	public:
		KKEvent(TChain* chain,int evt_num);
		
};
















KKEvent::KKEvent(TChain* chain,int evt_num){
	kkChain = chain ;
	kkChain->GetEntry(evt_num);
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
	kkChain->SetBranchAddress("pOrg",pOrg_);
	kkChain->SetBranchAddress("ukp",ukp_);
	kkChain->SetBranchAddress("vkp",vkp_);
	kkChain->SetBranchAddress("ukm",ukm_);
	kkChain->SetBranchAddress("vkm",vkm_);
	kkChain->SetBranchAddress("inside",inside_);
	for(int ikp=0;ikp<nKp_;++ikp){
		for(int ikm=0;ikm<nKm_;++ikm){
			KKBeam Beam(ukm_[nkk],vkm_[nkk],pK18_[nkk],chisqrK18_[nkk]);
			KKScat Scat(ukp_[ikp],vkp_[ikp],pKurama_[ikp],chisqrKurama_[ikp],m2_[ikp],qKurama_[ikp],inside_[ikp]);
			KKTrack Track(Beam,Scat);
			Tracks.push_back(Track);
			nkk++;
		}
	}
}




