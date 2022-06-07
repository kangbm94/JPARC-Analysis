#include "PhysicalConstants.hh"
double PiUpperThreshold=0.12;
double KLowerThreshold= 0.15;
double KUpperThreshold= 0.4;
double PLowerThreshold= 0.5;
double PUpperThreshold= 1.8;

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







