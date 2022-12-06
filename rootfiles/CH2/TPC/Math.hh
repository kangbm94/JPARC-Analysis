double ZTarget=-143;
const double& HS_field_0 = 0.9860;
const double& HS_field_Hall_calc = 0.90607;
const double& HS_field_Hall=  0.873800000;

double PI = acos(-1);
double mpi = 139.570/1000;//GeV
double mk = 493.677/1000;
double mp = 938.272/1000;
double mL = 1115.68/1000;
double mXi = 1321.71/1000;
double mXiStar = 1535/1000;


std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
TF1 fint("fint",s_tmp.c_str(),-4.,4.);


TVector3 VertexPointHelix(const Double_t par1[5], const Double_t par2[5],
		Double_t& dist, Double_t& t1, Double_t& t2)
{
	//helix function 1
	//x = [0] + [3]*cos(t);
	//y = [1] + [3]*sin(t);
	//z = [2] + [3]*[4]*t;

	//helix function 2
	//x = [5] + [8]*cos(t);
	//y = [6] + [8]*sin(t);
	//z = [7] + [8]*[9]*t;

	TF2 fvert_helix("fvert_helix",
			"pow(([0]+[3]*cos(x))-([5]+[8]*cos(y)),2)+pow(([1]+[3]*sin(x))-([6]+[8]*sin(y)),2)+pow(([2]+[3]*[4]*x)-([7]+[8]*[9]*y),2)",
			-5.,5.,-5.,5.);

	fvert_helix.SetParameter(0, par1[0]);
	fvert_helix.SetParameter(1, par1[1]);
	fvert_helix.SetParameter(2, par1[2]);
	fvert_helix.SetParameter(3, par1[3]);
	fvert_helix.SetParameter(4, par1[4]);
	fvert_helix.SetParameter(5, par2[0]);
	fvert_helix.SetParameter(6, par2[1]);
	fvert_helix.SetParameter(7, par2[2]);
	fvert_helix.SetParameter(8, par2[3]);
	fvert_helix.SetParameter(9, par2[4]);

	Double_t close_zin, close_zout;
	fvert_helix.GetMinimumXY(close_zin, close_zout);
	t1 = close_zin;
	t2 = close_zout;
	dist = TMath::Sqrt(fvert_helix.GetMinimum());

	Double_t xin = par1[0]+par1[3]*cos(close_zin);
	Double_t xout = par2[0]+par2[3]*cos(close_zout);
	Double_t yin =  par1[1]+par1[3]*sin(close_zin);
	Double_t yout = par2[1]+par2[3]*sin(close_zout);
	Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
	Double_t zout =  par2[2]+par2[3]*par2[4]*close_zout;

	// Double_t vx = (par1[0]+par1[3]*cos(close_zin) + par2[0]+par2[3]*cos(close_zout))/2.;
	// Double_t vy = (par1[1]+par1[3]*sin(close_zin) + par2[1]+par2[3]*sin(close_zout))/2.;
	// Double_t vz = (par1[2]+par1[3]*par1[4]*close_zin + par2[2]+par2[3]*par2[4]*close_zout)/2.;
	Double_t vx = (xin+xout)/2.;
	Double_t vy = (yin+yout)/2.;
	Double_t vz = (zin+zout)/2.;

	Double_t dist2 = sqrt(pow(xin-xout,2)
			+pow(yin-yout,2)
			+pow(zin-zout,2));
	dist = dist2;
	Double_t vertx = -1.*vx;
	Double_t verty = vz;
	Double_t vertz = vy + ZTarget;
	return TVector3(vertx, verty, vertz);
}


TVector3 CalcHelixMom(double par[5], double y)
{

	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

	double t = (y-par[2])/(par[3]*par[4]);
	double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c

	//From here!!!!
	double tmp_px = pt*(-1.*sin(t));
	double tmp_py = pt*(cos(t));
	double tmp_pz = pt*(par[4]);
	double px = -tmp_px*0.001;
	double py = tmp_pz*0.001;
	double pz = tmp_py*0.001;
	return TVector3(px,py,pz);
}



double GetTcal(double* par,TVector3 pos){
	TVector3 pos_(-pos.X(),
			pos.Z()-ZTarget,
			pos.Y());
	double fpar[8];
	for(int ip=0; ip<5; ++ip){
		fpar[ip] = par[ip];
	}
	fpar[5] = pos_.X();
	fpar[6] = pos_.Y();
	fpar[7] = pos_.Z();
	fint.SetParameters(fpar);
	double min_t = fint.GetMinimumX();
	return min_t;
}

TVector3 HelixPos(double* par, double t){
	double cx = par[0],cy=par[1],z0 = par[2],r = par[3],dz = par[4];
	double x = cx+r*cos(t);
	double y = cy+r*sin(t);
	double z = z0+r*dz*t;
	return TVector3(-x,z,y+ZTarget);
}

double MinHelixDistance(double* par, TVector3 pos){
	auto t = GetTcal(par,pos);
	auto v = HelixPos(par,t);
	return (pos-v).Mag();
}
class Particle :public TLorentzVector{
	private:
		int Q_=0;
		int pid_=0;
		double cd_ = 1e9;
		double hpar[5]={0};
		int tid_ = -1;
	public: 
		Particle(int pid, int Q, double* par,int tid){
			Q_  = Q;
			pid_= pid;
			tid_= tid;
			for(int i=0;i<5;++i) hpar[i] = par[i];
		}

		/*Particle(double m, int Q, TVector3 V,double* par){
			this->SetXYZM(V.X(),V.Y(),V.Z(),m);
			Q_=Q;
			for(int i=0;i<5;++i) hpar[i] = par[i];
			}
			*/
		bool IsP(){
			if(pid_%8/4) return true;
			else return false;
		}
		bool IsK(){
			if(pid_%4/2) return true;
			else return false;
		}
		bool IsPi(){
			if(pid_%2/1) return true;
			else return false;
		}
		bool GetQ(){
			return Q_;
		}
		bool GetID(){
			return tid_;
		}
		void SetCD(double cd){
			cd_ = cd;
		}
		double GetCD(){
			return cd_;
		}
		double* Gethpar(){
			return hpar;
		}
};
bool InTarget(TVector3 Vect){

}

class Recon{
	private:
		vector<TLorentzVector>Daughters;
		TLorentzVector LV;
		TVector3 Vert;
		bool exist = false;
		bool proper = false;
	public:
		Recon(vector<TLorentzVector> D){
			LV.SetXYZM(0,0,0,0);Daughters = D;
			for(auto lv:D) LV+=lv;	
			exist = true;
		}
		Recon(){}
		TVector3 GetVertex(){
			return Vert;
		}
		bool Exist(){
			return exist;
		}
		TLorentzVector GetLV(){
			return LV;
		}
		double Mass(){
			return LV.Mag();
		}
		TVector3 GetMomentum(){
			return LV.Vect();
		}
		bool Proper(){
			if(vert);	
		}
};


class Vertex{
	private:
		vector<Particle> Particles;
		TVector3 vert;
		vector<TVector3> verts;
		double cdcut = 10;
		void SetVert(){
			int n = verts.size();
			vert=TVector3(0,0,0);
			for(int i=0;i<n;++i)vert+=verts[i] * (1./n);
		}
		int Vert_id=0;
		vector<Recon>LdCand;
		vector<TVector3>Pmom;
		vector<TVector3>Pimom;
		vector<TVector3>Ldmom;
		vector<double>XiCand;
	public:
		Vertex(Particle p,int nt){
			Particles.push_back(p);Vert_id+=pow(2,nt);
//			cout<<"Vertex"<<endl;
		}
		int NParticle(){
			return Particles.size();}
		void AddParticle(Particle p,int nt){
			int np = NParticle();
			if(np==1){
				auto par1 = Particles[0].Gethpar();
				auto par2 = p.Gethpar();
				double cd,t1,t2;
				auto pos = VertexPointHelix(par1,par2,cd,t1,t2);
				if(cd<cdcut){
					Particles.push_back(p);
					verts.push_back(pos);Vert_id+=pow(2,nt);SetVert();}}
			else{
				double t = GetTcal(p.Gethpar(),vert);
				auto point = HelixPos(p.Gethpar(),t);
				double cd = (vert-point).Mag();
				double ct,t1,t2;
				if(cd<cdcut){
					vector<TVector3> cand;
					for(auto pt:Particles){
						cand.push_back(VertexPointHelix(pt.Gethpar(),p.Gethpar(),ct,t1,t2));
					}
					TVector3 vec(0,0,0);
					for(auto v:cand){vec+=v*(1./cand.size());	
						Particles.push_back(p);Vert_id+=pow(2,nt);verts.push_back(vec);
						SetVert();
					}
				}
			}
		}
		void SearchCombination(){
			int np = NParticle();
			vector<Particle>PCand;PCand.clear();
			vector<Particle>PiCand;PiCand.clear();
			for(auto p:Particles){
				if(p.IsP())PCand.push_back(p);
				if(p.IsPi())PiCand.push_back(p);
			}
			for(auto p:PCand){
				for(auto pi:PiCand){
					double cd_,t1_,t2_;
					if(p.GetID() == pi.GetID()) continue;
					auto ppi = VertexPointHelix(p.Gethpar(),pi.Gethpar(),cd_,t1_,t2_); 
					auto p1 = CalcHelixMom(p.Gethpar(),ppi.y());
					auto p2 = CalcHelixMom(pi.Gethpar(),ppi.y());
					auto pLV = TLorentzVector(p1,sqrt(mp*mp+p1.Mag2()));
					auto piLV = TLorentzVector(p2,sqrt(mpi*mpi+p2.Mag2()));
					vector<TLorentzVector> lv1 = {pLV,piLV};
					LdCand.push_back(Recon(lv1));
					auto piLVInv = TLorentzVector(-p2,sqrt(mpi*mpi+p2.Mag2()));
					vector<TLorentzVector> lv2 = {pLV,piLVInv};
					LdCand.push_back(Recon(lv2));
				}
			}
		}
		/*
		double GetLdIM(double& pm,double& pim,double& ldm){
			double comp = 9999;
			double val = -1;
			if(LdCand.size()==0) val= -1;
			for(int i=0 ; LdCand.size();++i){
				auto m = LdCand[i];
				if( abs(mL-m)<comp) {
					comp=abs(mL-m);
					val=m;
//					pm=Pmom[i];
//					pim=Pmom[i];
//					ldm=Pmom[i];
				}
			}
			return val;
		}
		*/
		Recon GetLd(){
			double comp = 9999;
			Recon val ;
			int num = 0;
			for(auto ldc : LdCand){ if( abs(mL-ldc.Mass())<comp) {comp=abs(mL-ldc.Mass());val=ldc;}}
			return val;
		}

		double GetXiIM(){
			double comp = 9999;
			double val = -1;
			if(XiCand.size()==0) val= -1;
			for(auto m : XiCand) if( abs(mXi-m)<comp) {comp=abs(mXi-m);val=m;}
			return val;
		}
};
