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
double PI = acos(-1);
double mpi = 139.570/1000;//GeV
double mk = 493.677/1000;
double mp = 938.272/1000;
double mL = 1115.68/1000;
double mXi = 1321.71/1000;
double mXiStar = 1535/1000;
//	TFile* file = new TFile("run05000_DstTPCHelixTracking.root");
//	TTree* tree = (TTree*)file->Get("tpc");
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
	 class Particle : public TLorentzVector{
	 private:
	 int q_ = 0;
	 public:
	 Particle(double m,int q, TVector3 mom){
	 q_=q;
	 this.SetXYZM(mom.X(),mom.Y(),mom.Z(),m); 
	 }
	 int GetCharge(){
	 return q_;
	 }
	 };
	 */
TLorentzVector particle(double m,TVector3 mom){
	double E = sqrt(m*m+mom.Mag2());
	return TLorentzVector(mom,E);
}
bool isP(TLorentzVector v){
	if(abs(v.Mag()-mp)<0.0100) return true;
	else return false;
}
bool isPi(TLorentzVector v){
	if(abs(v.Mag()-mpi)<0.0100) return true;
	else return false;
}
bool isK(TLorentzVector v){
	if(abs(v.Mag()-mk)<0.0100) return true;
	else return false;
}
void Test2(){
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
	TH1D* hist = new TH1D("cd","cd",1000,0,100);


	int Xirunnum,Xievnum;
	double Ximm;
	TFile* file2 = new TFile("SelectedEvents.root");
	TTree* tree2 = (TTree*)file2->Get("tree");
	tree2->SetBranchAddress("runnum",&Xirunnum);
	tree2->SetBranchAddress("evnum",&Xievnum);
	tree2->SetBranchAddress("XiM2",&Ximm);
	int xient = tree2->GetEntries();
	int nbin=200;
	TH1D* hist2 = new TH1D("LInvMass","LInvMass",nbin,1,2);
	TH1D* hist3 = new TH1D("KuramaMM","KuramaMM( |IM-MLd| < 0.042 )",nbin,1,2);
	//	TH1D* hist3 = new TH1D("VertexY","VertexY",100,-300,300);
	//	TH1D* hist3 = new TH1D("pid1","pid1",30,-1,2);
	TH1D* hist4 = new TH1D("LMom","LMom",100,0,2);
	TF1* fgaus = new TF1("fgaus","gaus",mL-0.05,mL+0.05);
	double cd_cut = 8;
	int cd_Count=0;
	double chi_cut = 50;
	TFile* Out = new TFile("TPCInv.root","recreate");
	TTree* outtr = new TTree("tree","tree");
	double inv = 0;
	double lp = 0;
	outtr->Branch("runnum",&runnum);
	outtr->Branch("evnum",&evnum);
	outtr->Branch("MM",&Ximm);
	outtr->Branch("InvM",&inv);
	outtr->Branch("LambdaP",&lp);
	for(int i=0;i<ent;++i){
		Clear();
		tree->GetEntry(i);
		if(ntTpc<1) continue;
		for(int j=0;j<xient;j++){
			tree2->GetEntry(j);
			if(Xirunnum==runnum and Xievnum == evnum) break;
		}

		vector<TLorentzVector>LCand;
		LCand.clear();
		vector<double> LIM;
		LIM.clear();
		vector<double> LP;
		LP.clear();

		for(int nt1 = 0; nt1<ntTpc;++nt1){
			if(chisqr->at(nt1)>chi_cut) continue; 
			int nh = helix_cx->size();
			double hcx = helix_cx->at(nt1);
			double hcy = helix_cy->at(nt1);
			double hz0 = helix_z0->at(nt1);
			double hr = helix_r->at(nt1);
			double hdz = helix_dz->at(nt1);
			double par1[5] = {hcx,hcy,hz0,hr,hdz};
			for(int nt2 = 0; nt2<ntTpc;++nt2){
				if( not (nt2<nt1)) continue;
				if(chisqr->at(nt2)>chi_cut) continue; 
				if(isBeam->at(nt1) or isBeam->at(nt2)) continue; 
				double hcx2 = helix_cx->at(nt2);
				double hcy2 = helix_cy->at(nt2);
				double hz02 = helix_z0->at(nt2);
				double hr2 = helix_r->at(nt2);
				double hdz2 = helix_dz->at(nt2);
				double par2[5] = {hcx2,hcy2,hz02,hr2,hdz2};
				double cd,t1,t2;
				auto vert = VertexPointHelix(par1,par2,cd,t1,t2);
				hist->Fill(cd);
				int id1 = pid->at(nt1);
				int id2 = pid->at(nt2);
				vector<TLorentzVector> cand1;
				vector<TLorentzVector> cand2;
				cand1.clear();
				cand2.clear();
				if(cd<cd_cut){
					cd_Count++;
					auto p1 = CalcHelixMom(par1,vert.y());
					auto p2 = CalcHelixMom(par2,vert.y());
					double q1 = charge->at(nt1);
					double q2 = charge->at(nt2);
					if((id1%8)/4==1)cand1.push_back(particle(mp,p1));	
					if((id1%4)/2==1)cand1.push_back(particle(mk,p1));	
					if((id1%2)==1)cand1.push_back(particle(mpi,p1));	
					if((id2%8)/4==1)cand2.push_back(particle(mp,p2));	
					if((id2%4)/2==1)cand2.push_back(particle(mk,p2));	
					if((id2%2)==1)cand2.push_back(particle(mpi,p2));
					if(q1 > 0 and q2 < 0){
						for(int ic1=0;ic1<cand1.size();++ic1){
							for(int ic2=0;ic2<cand2.size();++ic2){
								auto c1 = cand1[ic1];
								auto c2 = cand2[ic2];
								if(isP(c1) and isPi(c2)){
									LCand.push_back(c1+c2);
								}
							}	
						}
					}
					if(q1 < 0 and q2 > 0){
						for(int ic1=0;ic1<cand1.size();++ic1){
							for(int ic2=0;ic2<cand2.size();++ic2){
								auto c1 = cand1[ic1];
								auto c2 = cand2[ic2];
								if(isPi(c1) and isP(c2)){
									LCand.push_back(c1+c2);
								}
							}	
						}
					}
					for(auto h : LCand){
						LIM.push_back(h.Mag());
						LP.push_back(h.Vect().Mag());
						hist2->Fill(h.Mag());
					}
				}//cd Cut
			}//nt2
		}//nt1
		bool multi = false;
		for(auto m : LIM){
			if(multi) cout<<"MultiLambda?"<<endl;
			if(abs(m-mL)<0.042 and not multi){
				hist3->Fill(Ximm);
				inv = m;
				lp = LP[0];
				hist4->Fill(lp);
				outtr->Fill();
				multi = true;
			};
		}




	}//evt
	Out->Write();
	TCanvas* cv1 = new TCanvas("c1","c1",1200,800);
	cv1->Divide(2,2);
	cv1->cd(1);
	hist->Draw();
	cv1->cd(2);
	hist2->Draw();
	hist2->Fit("fgaus","R");
	cv1->cd(3);
	hist3->Draw();
	cv1->cd(4);
	hist4->Draw();
	//	cout<<hist2->GetEffectiveEntries()<<endl;
	cout<<cd_Count<<endl;
	double p1 = fgaus->GetParameter(0);
	double sig = fgaus->GetParameter(2);
	int cnt = p1*sig*sqrt(2*PI)*nbin;
	cout<<"NL = "<<cnt<<endl;
}



void Summary(){
	//	delete file;
	//	delete tree;
	file = new TFile("run05000_DstTPCHelixTracking.root");
	tree = (TTree*)file->Get("tpc");
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
	cout<<ent<<endl;
	TFile* out = new TFile("SelectedHelixx.root","recreate");
	TTree* trout = new TTree("tree","tree");
	trout->Branch("runnum",&runnum);
	trout->Branch("evnum",&evnum);
	trout->Branch("ntTpc",&ntTpc);
	trout->Branch("isBeam",&isBeam);
	trout->Branch("chisqr",&chisqr);
	trout->Branch("helix_cx",&helix_cx);
	trout->Branch("helix_cy",&helix_cy);
	trout->Branch("helix_z0",&helix_z0);
	trout->Branch("helix_r",&helix_r);
	trout->Branch("helix_dz",&helix_dz);
	trout->Branch("mom_vtx",&mom_vtx);
	trout->Branch("mom_vty",&mom_vty);
	trout->Branch("mom_vtz",&mom_vtz);
	trout->Branch("vtx",&vtx);
	trout->Branch("vty",&vty);
	trout->Branch("vtz",&vtz);
	trout->Branch("closeDist",&closeDist);
	trout->Branch("mom0",&mom0);
	trout->Branch("pid",&pid);
	trout->Branch("charge",&charge);
	for(int i=0;i<ent;++i){
		Clear();
		tree->GetEntry(i);
		if(i%10000==0) cout<<i<<endl;
		if(ntTpc<1) continue;
		trout->Fill();
		/*
			 for(int nt1 = 0; nt1<ntTpc;++nt1){
			 TVector3 vert(vtx->at(nt1),vty->at(nt1),vtz->at(nt1));
			 for(int nt2 = 0; nt2<ntTpc;++nt2){
			 if(nt1 == nt2) continue;
			 TVector3 vert2(vtx->at(nt2),vty->at(nt2),vtz->at(nt2));
			 }
			 }
			 */

	}
	out->Write();
}
