#include "Math.hh"
class PolaAnal{
	private:
		double mXi = 1.32171;
		double mp = 938.272/1000;
		double mL = 1115.683/1000;
		TVector3 Km;
		TVector3 Kp;
		TVector3 Xi;
		TVector3 Ld;
		TVector3 P;
		TVector3 PolXi;
		TVector3 PolLd;
		double cTh,Th,Ph,Ph1,Ph2,Ph3;
		TVector3 XLd;
		TVector3 YLd;
		TVector3 ZLd;

	public:
		PolaAnal(TVector3 PKm,TVector3 PKp,TVector3 PXi,TVector3 PLd,TVector3 PP){
			Km=PKm,Kp=PKp,Xi=PXi,Ld=PLd,P=PP;
			PolXi= NormalCross(Km,Kp);
			TLorentzVector LVXi(Xi,hypot(Xi.Mag(),mXi));
			TLorentzVector LVLd(Ld,hypot(Ld.Mag(),mL));
			TLorentzVector LVP(P,hypot(P.Mag(),mP));
			auto BoostXi = LVXi.BoostVector();
			auto LVLdCM = LVLd;
			LVLdCM.Boost(-BoostXi);
			auto LdCM = LVLdCM.Vect();
			cTh = NormalDot(LdCM,PolXi);
			Th=acos(cTh);
			ZLd = LdCM*(1./LdCM.Mag());
			XLd = NormalCross(PolXi,ZLd);
			YLd = NormalCross(ZLd,XLd);
			double PLx = BXi*sin(Th)/(1+AXi*cTh);
			double PLy = CXi*sin(Th)/(1+AXi*cTh);
			double PLz = (AXi+cTh)/(1+AXi*cTh);
	
			PolLd = PLx*XLd+PLy*YLd+PLz*ZLd;

			auto BoostLd = LVLd.BoostVector();
			auto LVPCM = LVP;
			LVPCM.Boost(-BoostLd);
			auto PCM = LVPCM.Vect();
			Ph = acos(NormalDot(PolLd,PCM));
			Ph1 = acos(NormalDot(XLd,PCM));
			Ph2 = acos(NormalDot(YLd,PCM));
			Ph3 = acos(NormalDot(ZLd,PCM));
		
		}
		double GetTheta(){
			return Th;
		}
		double GetPhi(){
			return Ph;
		}
		double GetPhi1(){
			return Ph1;
		}
		double GetPhi2(){
			return Ph2;
		}
		double GetPhi3(){
			return Ph3;
		}

};



