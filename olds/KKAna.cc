#include "KKManager.hh"
KKManager KM;
void KKAna(){
	//	KM.LoadFile("./rootfiles/CH2/AllKKAna.root");
	//	KM.LoadFile("./rootfiles/CH2/BigKKAna.root");
	KM.LoadFile("~/WS_data/ch2target/run05666_KKAnaTest.root");
	KM.LoadKK();
}
void GetPhi(){
	TChain* chain = KM.GetPublicChain();
	int ent = chain->GetEntries();
	cout<<ent<<endl;
	chain->SetBranchAddress("ntK18",&ntK18);
	chain->SetBranchAddress("pK18",pK18);
	chain->SetBranchAddress("utgtK18",utgtK18);
	chain->SetBranchAddress("vtgtK18",vtgtK18);
	chain->SetBranchAddress("ntKurama",&ntKurama);
	chain->SetBranchAddress("nKm",&nKm);
	chain->SetBranchAddress("nKp",&nKp);
	chain->SetBranchAddress("m2",m2);
	chain->SetBranchAddress("pKurama",pKurama);
	chain->SetBranchAddress("qKurama",qKurama);
	chain->SetBranchAddress("ukp",ukp);
	chain->SetBranchAddress("vkp",vkp);
	chain->SetBranchAddress("inside",inside);
	bool pass;
	int cnt=0;
	double MaxP=700/1000;
	double MinP=200/1000;
	double KaonWindow=0.15;
	double p_cut=0.75;
	double KaonM2 = (KaonMass*KaonMass)/1000/1000;
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TH1D* h = new TH1D("Phi","Phi",700,0.8,1.5);
	TH1D* hp = new TH1D("mom","mom",700,0.1,3);
	TH1D* hm = new TH1D("mass","mass",100,0.05,0.45);
	//	ent=ent/20;
	TLorentzVector Kaons[25];
	double MC[25];
	int idx[25]={0};
	TFile* file = new TFile("SortedPhi.root","recreate");
	TTree* tree = new TTree("tree","tree");
	tree->Branch("Beampx",&Beampx,"Beampx/D");
	tree->Branch("Beampy",&Beampy,"Beampy/D");
	tree->Branch("Beampz",&Beampz,"Beampz/D");
	tree->Branch("Kmpx",&Kmpx,"Kmpx/D");
	tree->Branch("Kmpy",&Kmpy,"Kmpy/D");
	tree->Branch("Kmpz",&Kmpz,"Kmpz/D");
	tree->Branch("Kppx",&Kppx,"Kppx/D");
	tree->Branch("Kppy",&Kppy,"Kppy/D");
	tree->Branch("Kopz",&Kppz,"Kppz/D");
	tree->Branch("Phipx",&Phipx,"Phipx/D");
	tree->Branch("Phipy",&Phipy,"Phipy/D");
	tree->Branch("Phipz",&Phipz,"Phipz/D");
	tree->Branch("PhiMass",&MissMassPhi,"PhiMass/D");
	tree->Branch("MissMassLambda",&MissMassLambda,"MissMassLambda/D");


	ent=ent/10;
	for(int i=0;i<ent;++i){
		Indicator(i,ent);
		chain->GetEntry(i);
		pass = true;
		double BeamNorm=sqrt(1+utgtK18[0]*utgtK18[0]+vtgtK18[0]*vtgtK18[0]);
		Beampx=pK18[0]*utgtK18[0]/BeamNorm;
		Beampy=pK18[0]*vtgtK18[0]/BeamNorm;
		Beampz=pK18[0]/BeamNorm;

		TLorentzVector KBeam(Beampx,Beampy,Beampz,sqrt(KaonM2+pK18[0]*pK18[0]));
		TLorentzVector PTarget(0,0,0,ProtonMass/1000);
		if(ntK18==nKm&&ntKurama==nKp&&nKm==1&&nKp==2){
			pass=false;//nKK
		}
		for(int j=0;j<nKp;++j){
			MC[j]=0;
			MC[j]=abs(m2[j]-KaonM2);
		}
		int nKaon=0;
		TMath::Sort(nKp,MC,idx,0);
		for(int j=0;j<nKp;++j){
			if(MC[idx[j]]<KaonWindow&&pKurama[idx[j]]<p_cut){
				nKaon++;
				//				cout<<MC[idx[j]]<<endl;
			}
		}
		if(nKaon<2){
			pass+=true;
		}
		double charge[25]={0};
		int counter = 0;
		double Pcont[25]={0};
		for(int j=0;j<nKaon;++j){
			int k = idx[j];
			if(pKurama[idx[j]]>p_cut){
				//				cout<<"CutByMomentum!"<<endl;
				//				cout<<pKurama[idx[j]]<<endl;
				//				cout<<MC[idx[j]]<<endl;
				continue;
			}
			double norm = sqrt(1+ukp[k]*ukp[k]+vkp[k]*vkp[k]);
			double px = pKurama[k]*ukp[k]/norm;
			double py = pKurama[k]*vkp[k]/norm;
			double pz = pKurama[k]/norm;
			double pt = sqrt(pKurama[k]*pKurama[k]+KaonM2);
			Kaons[counter]=TLorentzVector(px,py,pz,pt);
			charge[counter]=qKurama[k];
			Pcont[counter]=pKurama[k];
			counter++;
		}
		if(pass){
			continue;
		}
		nKaon=counter;
		//		cout<<Form("Kaons: %d",nKaon)<<endl;
		TLorentzVector Phi[25];
		double pairs[25][2];
		int nPhi=0;
		for(int j=0;j<nKaon;j++){
			for(int k=j+1;k<nKaon;k++){
				//				cout<<Form("(%d,%d,%d)",nKaon,j,k)<<endl;
				if(charge[j]*charge[k]<0){
					//				cout<<charge[j]<<" "<<charge[k]<<endl;
					Phi[nPhi]=Kaons[j]+Kaons[k];
					nPhi++;
					if(charge[j]<0){
						pairs[nPhi][0]=j;
						pairs[nPhi][1]=k;
					}
					else{
						pairs[nPhi][0]=k;
						pairs[nPhi][1]=j;
					}
				}
			}
		}
		//		cout<<Form("Phis: %d",nPhi)<<endl;
		double MassRes[25]={0};
		double MissMass[25]={0};
		for(int j=0;j<nPhi;++j){
			MassRes[j]=abs(Phi[j].Mag()-PhiMass/1000);
			MissMass[j]=Phi[j].Mag();
		}
		int phidx[25]={0};
		TMath::Sort(nPhi,MassRes,phidx,0);
		int Kmidx=pairs[phidx[0]][0];
		int Kpidx=pairs[phidx[0]][1];
		cnt++;
		//		cout<<MissMass[phidx[0]]<<endl;
		MissMassPhi=MissMass[phidx[0]];

		TLorentzVector MissLambda = (KBeam+PTarget) -Phi[phidx[0]];

		h->Fill(MissMassPhi);
		hp->Fill(Kaons[Kmidx].P());
		hm->Fill(m2[0]);
		Kmpx=Kaons[Kmidx].Px();
		Kmpy=Kaons[Kmidx].Py();
		Kmpz=Kaons[Kmidx].Pz();
		Kppx=Kaons[Kpidx].Px();
		Kppy=Kaons[Kpidx].Py();
		Kppz=Kaons[Kpidx].Pz();
		Phipx=Phi[phidx[0]].Px();
		Phipy=Phi[phidx[0]].Py();
		Phipz=Phi[phidx[0]].Pz();
		MissMassLambda = MissLambda.Mag();
		tree->Fill();
	}
	cout<<cnt<<endl;
	c1->cd(1);
	h->Draw();
	c1->cd(2);
	hp->Draw();
	c1->cd(3);
	hm->Draw();
	file->Write();
}

