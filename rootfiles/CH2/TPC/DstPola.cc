#include "PolarizationAnal.hh"
double mk = 493.677/1000;
double mXi = 1321.71/1000;
double mp = 938.272/1000;
TH1D* HTh;
TH1D* HPh;
TH1D* HPh1;
TH1D* HPh2;
TH1D* HPh3;
vector<double>*DecayMomx;
vector<double>*DecayMomy;
vector<double>*DecayMomz;
vector<double> *MM;
vector<double>*PKp;
vector<double>*uKp;
vector<double>*vKp;
vector<double>*PKm;
vector<double>*uKm;
vector<double>*vKm;
vector<int>* trigflag;
int TAPS = 20,TB=15;
double kpPx,kpPy,kpPz;
double kmPx,kmPy,kmPz;
double MomxXi,MomyXi,MomzXi;
double MomxLd,MomyLd,MomzLd;
double MomxP,MomyP,MomzP;
double dM,dP;
double dM_x,dM_y,dM_z,dM_e,MM_Fermi;
bool FlgXi;
double InvMXi,InvMLd;
double cTh,cPh,cPh1,cPh2,cPh3;
double MissMass,Coplanarity;
bool TrigAPS,TrigB;
void RealPola(int i);
TTree* tree2;
void DstPola(){
	HTh = new TH1D("HistTh","HistTh",22,-1.1,1.1);
	HPh = new TH1D("HistPh","HistPh",22,-1.1,1.1);
	HPh1 = new TH1D("HistPh1","HistPh1",22,-1.1,1.1);
	HPh2 = new TH1D("HistPh2","HistPh2",22,-1.1,1.1);
	HPh3 = new TH1D("HistPh3","HistPh3",22,-1.1,1.1);
	TF1* fLin = new TF1("fLin","[0]+[1]*x",-1,1);
	
	int rn_s = 5641;
	int rn_e = 5666;
	rn_s = 5567;
//	rn_e = 5639;
//	rn_s = 5667;
	rn_e = 5697;
	TFile* file2 = new TFile(Form("DstPolaAnal_%d_%d.root",rn_s,rn_e),"recreate");
	tree2 = new TTree("tree","tree");
	tree2->Branch("TrigAPS",&TrigAPS);
	tree2->Branch("TrigB",&TrigB);
	tree2->Branch("MissMass",&MissMass);
	tree2->Branch("Coplanarity",&Coplanarity);
	tree2->Branch("cTh",&cTh);
	tree2->Branch("cPh",&cPh);
	tree2->Branch("cPh1",&cPh1);
	tree2->Branch("cPh2",&cPh2);
	tree2->Branch("cPh3",&cPh3);
	tree2->Branch("dM",&dM);	
	tree2->Branch("dM_x",&dM_x);	
	tree2->Branch("dM_y",&dM_y);	
	tree2->Branch("dM_z",&dM_z);	
	tree2->Branch("dM_e",&dM_e);	
	tree2->Branch("MM_Fermi",&MM_Fermi);	
	
	tree2->Branch("InvMXi",&InvMXi);
	tree2->Branch("InvMLd",&InvMLd);
	
	tree2->Branch("PKm_x",&kmPx);
	tree2->Branch("PKm_y",&kmPy);
	tree2->Branch("PKm_z",&kmPz);
	
	tree2->Branch("PKp_x",&kpPx);
	tree2->Branch("PKp_y",&kpPy);
	tree2->Branch("PKp_z",&kpPz);
	
	tree2->Branch("PXi_x",&MomxXi);
	tree2->Branch("PXi_y",&MomyXi);
	tree2->Branch("PXi_z",&MomzXi);
	
	tree2->Branch("PLd_x",&MomxLd);
	tree2->Branch("PLd_y",&MomyLd);
	tree2->Branch("PLd_z",&MomzLd);
	
	tree2->Branch("PP_x",&MomxP);
	tree2->Branch("PP_y",&MomyP);
	tree2->Branch("PP_z",&MomzP);



	for(int i=rn_s;i<rn_e+1;++i){

		cout<<"Run "<<i<<endl;
		if(5641 <= i and i <= 5666) continue;
		RealPola(i);
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	c1->cd(1);
	HTh->Draw();
	HTh->Fit("fLin");
	double p1 = fLin->GetParameter(1);
	double nb = (HTh->GetEntries())*1./20;
	cout<<p1/nb<<endl;
	c1->cd(2);
	HPh1->Draw();
	HPh1->Fit("fLin");
	p1 = fLin->GetParameter(1);
	nb = (HPh1->GetEntries())/20;
	cout<<p1/nb<<endl;
	c1->cd(3);
	HPh2->Draw();
	HPh2->Fit("fLin");
	p1 = fLin->GetParameter(1);
	nb = (HPh2->GetEntries())/20;
	cout<<p1/nb<<endl;
	c1->cd(4);
	HPh3->Draw();
	HPh3->Fit("fLin");
	p1 = fLin->GetParameter(1);
	nb = (HPh3->GetEntries())/20;
	cout<<p1/nb<<endl;
	file2->Write();
}
void RealPola(int rn){
	TFile* file =nullptr;
	file = TFile::Open(Form("genfitfiles/Prod/run0%d_GenfitXiSearch.root",rn));
//	TFile* file = new TFile(Form("genfitfiles/Ch2/run0%d_GenfitXiSearch.root",rn));
	if(!file or file->IsZombie()) return;
	TTree* tree = (TTree*)file->Get("tpc");

	int ent = tree->GetEntries();
	tree->SetBranchAddress("trigflag",&trigflag);
//	tree->SetBranchAddress("MissMassCorrDETPC",&MM);
	tree->SetBranchAddress("MissMassCorrDE",&MM);
	tree->SetBranchAddress("pK18",&PKm);
	tree->SetBranchAddress("utgtK18",&uKm);
	tree->SetBranchAddress("vtgtK18",&vKm);
	
	tree->SetBranchAddress("pCorrDETPC",&PKp);
	tree->SetBranchAddress("utgtTPCKurama",&uKp);
	tree->SetBranchAddress("vtgtTPCKurama",&vKp);

	tree->SetBranchAddress("XiMom_x",&MomxXi);
	tree->SetBranchAddress("XiMom_y",&MomyXi);
	tree->SetBranchAddress("XiMom_z",&MomzXi);
	tree->SetBranchAddress("LambdaMom_x",&MomxLd);
	tree->SetBranchAddress("LambdaMom_y",&MomyLd);
	tree->SetBranchAddress("LambdaMom_z",&MomzLd);
	tree->SetBranchAddress("DecaysMom_x",&DecayMomx);
	tree->SetBranchAddress("DecaysMom_y",&DecayMomy);
	tree->SetBranchAddress("DecaysMom_z",&DecayMomz);
	tree->SetBranchAddress("Xiflag",&FlgXi);
	tree->SetBranchAddress("XiMass",&InvMXi);
	tree->SetBranchAddress("LambdaMass",&InvMLd);


	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		if(!FlgXi	)continue;
		int TrgAPS = trigflag->at(TAPS);
		int TrgB = trigflag->at(TB);
		TrigAPS = TrgAPS;
		TrigB = TrgB;
		//		cout<<MM->size()<<endl;
		MissMass = MM->at(0);
//		if(abs(InvMLd-1.115)>0.02) continue;
//		if(abs(InvMXi-1.321)>0.04) continue;
		double MomKm = PKm->at(0);
		double UKm = uKm->at(0);
		double VKm = vKm->at(0);
		double NKm = 1./hypot(1,hypot(UKm,VKm));
		kmPx = UKm*MomKm;
		kmPy = VKm*MomKm;
		kmPz = NKm*MomKm;
		TVector3 Km(kmPx,kmPy,kmPz);
		
		double MomKp = PKp->at(0);
		double UKp = uKp->at(0);
		double VKp = vKp->at(0);
		double NKp = 1./hypot(1,hypot(UKp,VKp));
		kpPx = UKp*MomKp;
		kpPy = VKp*MomKp;
		kpPz = NKp*MomKp;
		TVector3 Kp(kpPx,kpPy,kpPz);
		TVector3 Xi(MomxXi,MomyXi,MomzXi);
		TVector3 Ld(MomxLd,MomyLd,MomzLd);

		TLorentzVector LVKp(Kp,hypot(Kp.Mag(),mk));
		TLorentzVector LVKm(Km,hypot(Km.Mag(),mk));
		TLorentzVector LVXi(Xi,hypot(Xi.Mag(),mXi));
		TLorentzVector LVP(0,0,0,mp);
		TLorentzVector LVMM = LVKm+LVP-LVKp;
		
		dM = (LVMM-LVXi).Mag();
		dM_x = (LVMM-LVXi).X();
		dM_y = (LVMM-LVXi).Y();
		dM_z = (LVMM-LVXi).Z();
		dM_e = (LVMM-LVXi).E();
		
		auto LV_Fermi = LVP+LVXi - LVMM;
		MM_Fermi = LV_Fermi.Mag();

		dP = LVMM.Vect().Mag() - LVXi.Vect().Mag();
		dP = LVMM.Vect().Mag() - LVXi.Vect().Mag();
		auto Plane = Km.Cross(Kp);
		Plane = Plane *(1./Plane.Mag());
		auto XiDir = Xi*(1./Xi.Mag());
		Coplanarity = XiDir*Plane;

		MomxP = DecayMomx->at(0);
		MomyP = DecayMomy->at(0);
		MomzP = DecayMomz->at(0);
		TVector3 P(MomxP,MomyP,MomzP);
		

		PolaAnal Pol(Km,Kp,Xi,Ld,P);
		
		double Th = Pol.GetTheta();
		double Ph = Pol.GetPhi();
		double Ph1 = Pol.GetPhi1();
		double Ph2 = Pol.GetPhi2();
		double Ph3 = Pol.GetPhi3();
		cTh = cos(Th);
		cPh = cos(Ph);
		cPh1 = cos(Ph1);
		cPh2 = cos(Ph2);
		cPh3 = cos(Ph3);
		tree2->Fill();
		if(abs(MissMass-1.321)>0.1)continue;
		if(TrigAPS){
		}
		if(TrigB){
			HTh->Fill(cTh);
			HPh->Fill(cPh);
			HPh1->Fill(cPh1);
			HPh2->Fill(cPh2);
			HPh3->Fill(cPh3);
		}
	}
}
